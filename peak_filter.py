
# Step 1: Read the Hs_EPDnew.bed file to get TSS positions
bed_df = pd.read_csv("Hs_EPDnew.bed", sep=" ", header=None, 
                     names=["chromosome", "region_start", "region_end", "name", "score", "strand", "thick_start", "thick_end"])

# Step 2: Calculate TSS position for each transcript based on strand information in the BED file
tss_positions = {}
for idx, row in bed_df.iterrows():
    transcript = row["name"]
    chrom = row["chromosome"]
    
    # For + strand, TSS is at thick_start
    # For - strand, TSS is at thick_end
    if row["strand"] == "+":
        tss = row["thick_start"]
    else:  # "-" strand
        tss = row["thick_end"]
    
    tss_positions[transcript] = {
        "chromosome": chrom,
        "tss": tss,
        "strand": row["strand"]
    }

# Step 3: Process the original dataframe (df) to identify merged peaks
transcript_peaks = {}
for transcript, group in df.groupby('f'):
    # Skip transcripts that don't have TSS information
    if transcript not in tss_positions:
        continue
        
    # Sort by chromosome and position
    group = group.sort_values(['a', 'b'])
    
    # Process each chromosome separately
    peaks = []
    for chrom, chrom_group in group.groupby('a'):
        # Skip if chromosome doesn't match the one in TSS data
        if chrom != tss_positions[transcript]["chromosome"]:
            continue
            
        chrom_group = chrom_group.sort_values('b')
        
        # Initialize with first region
        if len(chrom_group) > 0:
            current_peak = {
                "chromosome": chrom,
                "start": chrom_group.iloc[0]['b'],
                "end": chrom_group.iloc[0]['c'],
                "height": chrom_group.iloc[0]['e']
            }
            
            # Merge adjacent regions
            for idx, row in chrom_group.iloc[1:].iterrows():
                if row['b'] == current_peak["end"]:
                    # Merge this region with the current peak
                    current_peak["end"] = row['c']
                    current_peak["height"] = max(current_peak["height"], row['e'])
                else:
                    # Store the previous peak and start a new one
                    peaks.append(current_peak)
                    current_peak = {
                        "chromosome": chrom,
                        "start": row['b'],
                        "end": row['c'],
                        "height": row['e']
                    }
            
            # Add the last peak
            peaks.append(current_peak)
    
    transcript_peaks[transcript] = peaks

# Step 4: Identify transcripts with peaks ONLY in the +/-500 range (not in 500-3000 range)
transcripts_to_retain = []
for transcript, peaks in transcript_peaks.items():
    tss = tss_positions[transcript]["tss"]
    
    # Check for peaks within +/-500 of TSS
    has_peaks_near_tss = False
    has_peaks_in_forbidden_range = False
    
    for peak in peaks:
        # Check if peak is within +/-500 of TSS
        if peak["start"] <= tss + 500 and peak["end"] >= tss - 500:
            has_peaks_near_tss = True
        
        # Check if peak is in the forbidden range: 500-3000 from TSS (either direction)
        # Two cases to check:
        # 1. Peak starts in the +500 to +3000 range
        if peak["start"] > tss + 500 and peak["start"] <= tss + 3000:
            has_peaks_in_forbidden_range = True
        # 2. Peak ends in the -3000 to -500 range
        if peak["end"] < tss - 500 and peak["end"] >= tss - 3000:
            has_peaks_in_forbidden_range = True
    
    # Retain transcript if it has peaks near TSS AND doesn't have peaks in the forbidden range
    if has_peaks_near_tss and not has_peaks_in_forbidden_range:
        transcripts_to_retain.append(transcript)

# Step 5: Filter the original dataframe
filtered_df = df[df['f'].isin(transcripts_to_retain)]

# Step 6: Display results
print(f"Found {len(transcripts_to_retain)} transcripts with:")
print(f"  - At least one peak within +/-500 of TSS")
print(f"  - NO peaks in the 500-3000 range of TSS")
print(f"Original dataframe had {len(df['f'].unique())} unique transcripts")
print(f"Filtered dataframe has {len(filtered_df['f'].unique())} unique transcripts")
