#!/usr/bin/env python3

"""
Script to filter narrowPeak files for optimal motif analysis with MEME
Usage: python filter_peaks_for_motif.py input.narrowPeak output.narrowPeak [options]
"""

import argparse
import pandas as pd
import numpy as np
import sys

def filter_peaks_for_motif_analysis(input_file, output_file, 
                                   top_n=500, 
                                   min_score=50, 
                                   min_signal=5.0,
                                   min_pvalue=5.0,
                                   max_length=100,
                                   min_length=50,
                                   summit_window=200,
                                   remove_blacklist=True):
    """
    Filter peaks for optimal motif discovery
    
    Parameters:
    - top_n: Number of top peaks to keep (default: 500)
    - min_score: Minimum peak score (column 5)
    - min_signal: Minimum signal value (column 7)
    - min_pvalue: Minimum -log10(pvalue) (column 8)
    - max_length: Maximum peak length
    - min_length: Minimum peak length
    - summit_window: Window around summit to extract (if using summit)
    - remove_blacklist: Remove common blacklisted regions
    """
    
    # Read narrowPeak file
    columns = ['chrom', 'start', 'end', 'name', 'score', 'strand', 
               'signalValue', 'pValue', 'qValue', 'summit']
    
    try:
        df = pd.read_csv(input_file, sep='\t', names=columns, header=None)
        print(f"Loaded {len(df)} peaks from {input_file}")
    except Exception as e:
        print(f"Error reading file: {e}")
        return
    
    initial_count = len(df)
    
    # 1. Basic quality filters
    df = df[df['score'] >= min_score]
    print(f"After score filter (>={min_score}): {len(df)} peaks")
    
    df = df[df['signalValue'] >= min_signal]
    print(f"After signal filter (>={min_signal}): {len(df)} peaks")
    
    df = df[df['pValue'] >= min_pvalue]
    print(f"After p-value filter (>={min_pvalue}): {len(df)} peaks")
    
    # 2. Length filters
    df['length'] = df['end'] - df['start']
    df = df[(df['length'] >= min_length) & (df['length'] <= max_length)]
    print(f"After length filter ({min_length}-{max_length}bp): {len(df)} peaks")
    
    # 3. Remove common blacklisted regions (basic filter)
    if remove_blacklist:
        # Remove centromeres, heterochromatin, and repetitive regions
        blacklist_patterns = ['chrY', 'chrM', 'chr.*_random', 'chr.*_alt']
        for pattern in blacklist_patterns:
            df = df[~df['chrom'].str.contains(pattern, na=False)]
        print(f"After blacklist filter: {len(df)} peaks")
    
    # 4. Remove peaks in highly repetitive regions (simple heuristic)
    # Keep only standard chromosomes
    standard_chroms = [f'chr{i}' for i in range(1, 23)] + ['chrX']
    df = df[df['chrom'].isin(standard_chroms)]
    print(f"After standard chromosome filter: {len(df)} peaks")
    
    # 5. Sort by multiple criteria and take top peaks
    # Primary: p-value, Secondary: signal value, Tertiary: score
    df = df.sort_values(['pValue', 'signalValue', 'score'], 
                       ascending=[False, False, False])
    
    # Take top N peaks
    if len(df) > top_n:
        df = df.head(top_n)
        print(f"Selected top {top_n} peaks")
    
    # 6. Optional: Create summit-centered regions for better motif finding
    if summit_window > 0:
        df_summit = df.copy()
        # Summit is relative to start position
        df_summit['summit_abs'] = df_summit['start'] + df_summit['summit']
        df_summit['start'] = df_summit['summit_abs'] - summit_window // 2
        df_summit['end'] = df_summit['summit_abs'] + summit_window // 2
        
        # Ensure coordinates are positive
        df_summit['start'] = np.maximum(df_summit['start'], 0)
        
        # Recalculate summit relative position for new coordinates
        df_summit['summit'] = summit_window // 2
        
        # Drop the temporary column
        df_summit = df_summit.drop('summit_abs', axis=1)
        
        # Save summit-centered version
        summit_output = output_file.replace('.narrowPeak', '_summit.narrowPeak')
        df_summit[columns].to_csv(summit_output, sep='\t', header=False, index=False)
        print(f"Summit-centered peaks ({summit_window}bp) saved to: {summit_output}")
    
    # 7. Save filtered peaks
    df[columns].to_csv(output_file, sep='\t', header=False, index=False)
    
    print(f"\nFiltering Summary:")
    print(f"Input peaks: {initial_count}")
    print(f"Output peaks: {len(df)}")
    print(f"Retention rate: {len(df)/initial_count*100:.1f}%")
    print(f"Filtered peaks saved to: {output_file}")
    
    # Print some statistics
    print(f"\nPeak Statistics:")
    print(f"Score range: {df['score'].min():.1f} - {df['score'].max():.1f}")
    print(f"Signal range: {df['signalValue'].min():.2f} - {df['signalValue'].max():.2f}")
    print(f"P-value range: {df['pValue'].min():.2f} - {df['pValue'].max():.2f}")
    print(f"Length range: {df['length'].min()} - {df['length'].max()} bp")
    print(f"Median length: {df['length'].median():.0f} bp")

def main():
    parser = argparse.ArgumentParser(
        description="Filter narrowPeak files for optimal motif analysis",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Filter strategies for motif analysis:

1. CONSERVATIVE (high confidence):
   python filter_peaks_for_motif.py input.narrowPeak output.narrowPeak \\
          --top-n 200 --min-score 100 --min-signal 10 --min-pvalue 10

2. BALANCED (recommended for MYC):
   python filter_peaks_for_motif.py input.narrowPeak output.narrowPeak \\
          --top-n 500 --min-score 50 --min-signal 5 --summit-window 200

3. LIBERAL (more peaks, may include noise):
   python filter_peaks_for_motif.py input.narrowPeak output.narrowPeak \\
          --top-n 1000 --min-score 20 --min-signal 2

Examples:
   # Basic filtering for MYC peaks
   python filter_peaks_for_motif.py H82_c_MYC_peaks.narrowPeak H82_c_MYC_filtered.narrowPeak

   # High-confidence peaks only
   python filter_peaks_for_motif.py H82_c_MYC_peaks.narrowPeak H82_c_MYC_high_conf.narrowPeak \\
          --top-n 200 --min-score 100 --summit-window 150
        """
    )
    
    parser.add_argument('input', help='Input narrowPeak file')
    parser.add_argument('output', help='Output filtered narrowPeak file')
    parser.add_argument('--top-n', type=int, default=500,
                       help='Number of top peaks to keep (default: 500)')
    parser.add_argument('--min-score', type=float, default=50,
                       help='Minimum peak score (default: 50)')
    parser.add_argument('--min-signal', type=float, default=5.0,
                       help='Minimum signal value (default: 5.0)')
    parser.add_argument('--min-pvalue', type=float, default=5.0,
                       help='Minimum -log10(pvalue) (default: 5.0)')
    parser.add_argument('--max-length', type=int, default=1000,
                       help='Maximum peak length (default: 1000)')
    parser.add_argument('--min-length', type=int, default=50,
                       help='Minimum peak length (default: 50)')
    parser.add_argument('--summit-window', type=int, default=200,
                       help='Window around summit (0 to disable, default: 200)')
    parser.add_argument('--no-blacklist', action='store_true',
                       help='Skip blacklist filtering')
    
    args = parser.parse_args()
    
    filter_peaks_for_motif_analysis(
        args.input, 
        args.output,
        top_n=args.top_n,
        min_score=args.min_score,
        min_signal=args.min_signal,
        min_pvalue=args.min_pvalue,
        max_length=args.max_length,
        min_length=args.min_length,
        summit_window=args.summit_window,
        remove_blacklist=not args.no_blacklist
    )

if __name__ == "__main__":
    main()
