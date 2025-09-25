#PBS -l walltime=20:00:00
#PBS -l nodes=1:ppn=4
#PBS -A PAS1348
echo $PWD
sh ~/.bash_profile
source activate macs3
cd /fs/scratch/PAS1348/tmp/mapping/bam/rmdups/

#source ~/.bash_profile

#module load miniconda3

#plotFingerprint -b *.bam --plotFileFormat svg --plotFile test.svg --ignoreDuplicates --minMappingQuality 30 --skipZeros --smartLabels -p 20
#plotFingerprint -b A1.rmdups.bam A2.rmdups.bam B1.rmdups.bam C1.rmdups.bam Y1.rmdups.bam Z1.rmdups.bam -T "H82 input" --plotFileFormat svg --plotFile H82_input.svg --ignoreDuplicates --minMappingQuality 30 --smartLabels -p 20
#plotFingerprint -b S1.rmdups.bam T1.rmdups.bam U1.rmdups.bam E2.rmdups.bam F2.rmdups.bam G2.rmdups.bam -T "H446 input" --plotFileFormat svg --plotFile H446_input.svg --ignoreDuplicates --minMappingQuality 30 --smartLabels -p 20
#plotFingerprint -b D1.rmdups.bam E1.rmdups.bam F1.rmdups.bam B2.rmdups.bam C2.rmdups.bam D2.rmdups.bam G1.rmdups.bam H1.rmdups.bam I1.rmdups.bam J1.rmdups.bam K1.rmdups.bam L1.rmdups.bam --labels IgG#4 IgG#5 IgG#6 H3K4me3#28 H3K4me3#29 H3K4me3#30 E2F1#7 E2F1#8  E2F1#9 c_MYC#10 c_MYC#11  c_MYC#12 -T "H82 E2F1 and c-MYC TF's" --plotFileFormat svg --plotFile H82_E2F1_c-MYC.svg --ignoreDuplicates --minMappingQuality 30 --smartLabels -p 20

#macs3 callpeak -t D1.rmdups.bam E1.rmdups.bam F1.rmdups.bam -c A1.rmdups.bam A2.rmdups.bam B1.rmdups.bam C1.rmdups.bam Y1.rmdups.bam Z1.rmdups.bam -f BAMPE -g 2.9e9 -n H82_IgG --outdir macs3 -B
#macs3 callpeak -t B2.rmdups.bam C2.rmdups.bam D2.rmdups.bam -c A1.rmdups.bam A2.rmdups.bam B1.rmdups.bam C1.rmdups.bam Y1.rmdups.bam Z1.rmdups.bam -f BAMPE -g 2.9e9 -n H82_H3K4me3 --outdir macs3 -B
macs3 callpeak -t G1.rmdups.bam H1.rmdups.bam -c A1.rmdups.bam A2.rmdups.bam B1.rmdups.bam C1.rmdups.bam Y1.rmdups.bam Z1.rmdups.bam -f BAMPE -g 2.9e9 -n H82_E2F1_082724 --outdir macs3 -B
#macs3 callpeak -t J1.rmdups.bam K1.rmdups.bam L1.rmdups.bam -c A1.rmdups.bam A2.rmdups.bam B1.rmdups.bam C1.rmdups.bam Y1.rmdups.bam Z1.rmdups.bam -f BAMPE -g 2.9e9 -n H82_c_MYC --outdir macs3 -B


#macs3 callpeak -t V1.rmdups.bam W1.rmdups.bam X1.rmdups.bam -c S1.rmdups.bam T1.rmdups.bam U1.rmdups.bam E2.rmdups.bam F2.rmdups.bam G2.rmdups.bam -f BAMPE -g 2.9e9 -n H446_IgG --outdir macs3 -B
#macs3 callpeak -t H2.rmdups.bam I2.rmdups.bam J2.rmdups.bam -c S1.rmdups.bam T1.rmdups.bam U1.rmdups.bam E2.rmdups.bam F2.rmdups.bam G2.rmdups.bam -f BAMPE -g 2.9e9 -n H446_H3K4me3 --outdir macs3 -B
#macs3 callpeak -t M1.rmdups.bam N1.rmdups.bam O1.rmdups.bam -c S1.rmdups.bam T1.rmdups.bam U1.rmdups.bam E2.rmdups.bam F2.rmdups.bam G2.rmdups.bam -f BAMPE -g 2.9e9 -n H446_E2F1 --outdir macs3 -B
#macs3 callpeak -t P1.rmdups.bam Q1.rmdups.bam R1.rmdups.bam -c S1.rmdups.bam T1.rmdups.bam U1.rmdups.bam E2.rmdups.bam F2.rmdups.bam G2.rmdups.bam -f BAMPE -g 2.9e9 -n H446_c_MYC --outdir macs3 -B
#plotFingerprint -b V1.rmdups.bam W1.rmdups.bam X1.rmdups.bam H2.rmdups.bam I2.rmdups.bam J2.rmdups.bam M1.rmdups.bam N1.rmdups.bam O1.rmdups.bam P1.rmdups.bam Q1.rmdups.bam R1.rmdups.bam --labels IgG#22 IgG#23 IgG#24 H3K4me3#34 H3K4me3#35 H3K4me3#36 E2F1#13 E2F1#14  E2F1#15 c_MYC#16 c_MYC#17  c_MYC#18 -T "H446 E2F1 and c-MYC TF's" --plotFileFormat svg --plotFile H446_E2F1_c-MYC.svg --ignoreDuplicates --minMappingQuality 30 --smartLabels -p 20

source deactivate
