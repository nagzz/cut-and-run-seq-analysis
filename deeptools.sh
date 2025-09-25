#PBS -l walltime=10:00:00
#PBS -l nodes=1:ppn=2
#PBS -A PAS1348
#PBS -t 1-6
#PBS -M srik05@osumc.edu
#PBS -V
echo $PWD
sh ~/.bash_profile
source activate deeptools
cd /fs/ess/PAS1348/nagesh/cutandrun_hannahs/output/mapping/bam/rmdups/macs3

#INFILE1=`awk "NR==$PBS_ARRAY_INDEX" treatmentbwFiles.txt`
INFILE1=`awk "NR==$PBS_ARRAY_INDEX" bwFiles.txt`
#INFILE2=`awk "NR==$PBS_ARRAY_INDEX" controlbwFiles.txt`
#source ~/.bash_profile
current_time=$(date "+%Y_%m_%d_%H_%M_%S")
#module load miniconda3
echo $current_time
#ope="subtract"
#echo "$ope operation"
#computeMatrix reference-point --referencePoint center -S H446_c_MYC.bw H446_E2F1.bw H446_IgG.bw -R test.bed --binSize 1 -b 3000 -a 3000 --skipZeros -o test_H446.mat.gz -p 5
#plotProfile -m H446.mat.gz -out test_H446.svg --perGroup --plotTitle "Cut&Run seq signal plot" --refPointLabel promoter
#computeMatrix reference-point --referencePoint center -S H82_c_MYC.bw H82_E2F1.bw H82_IgG.bw -R test.bed --binSize 1 -b 3000 -a 3000 --skipZeros -o test_H82.mat.gz -p 5
#plotProfile -m H82.mat.gz -out test_H82.svg --perGroup --plotTitle "Cut&Run seq signal plot" --refPointLabel promoter
#bigwigCompare -b1 H446_E2F1.bw -b2 H446_IgG.bw --operation first -bl black_list_regions.bed --binSize 1 -p 10 -o H446_E2F1_IgG.bw -of bigwig
#bigwigCompare -b1 H446_IgG.bw -b2 H446_E2F1.bw --operation subtract --binSize 1 -p 10 -o H446_IgG_E2F1.bw -of bigwig

#echo "Processing files ${INFILE1} ${INFILE2}"
echo "Processing files ${INFILE1}"

#echo "bigwigCompare initiated"
#bigwigCompare -b1 ${INFILE1} -b2 ${INFILE2} -bl black_list_regions.bed --operation $ope -p 25 -o ${INFILE1%%.bw}_bigwigcompare_${current_time}.bw -of bigwig
#bigwigCompare -b1 H82_E2F1_082724.bw -b2 H82_IgG.bw -bl black_list_regions.bed --operation $ope -p 25 -o H82_E2F1_082724_bigwigcompare_${current_time}.bw -of bigwig

echo "converting bigwig to wig"
#bigWigToWig ${INFILE1%%.bw}_bigwigcompare_${current_time}.bw ${INFILE1%%.bw}_bigwigcompare_${current_time}.wig
bigWigToWig ${INFILE1} ${INFILE1%%bigwigcompare2*.bw}bigwigcompare2_${current_time}.wig

#bigWigToWig H82_E2F1_082724_bigwigcompare_${current_time}.bw H82_E2F1_082724_bigwigcompare_${current_time}.wig

echo "converting wig to bed format"
#convert2bed --input=wig -x <${INFILE1%%.bw}_bigwigcompare_${current_time}.wig> ${INFILE1%%.bw}_bigwigcompare_${current_time}.bed
convert2bed --input=wig -x <${INFILE1%%bigwigcompare2*.bw}bigwigcompare2_${current_time}.wig> ${INFILE1%%bigwigcompare2*.bw}bigwigcompare2_${current_time}.bed

#convert2bed --input=wig -x <H82_E2F1_082724_bigwigcompare_${current_time}.wig> H82_E2F1_082724_bigwigcompare_${current_time}.bed

#awk 'BEGIN {FS="\t"} {if ($5<=0) print $1"\t"$2"\t"$3}' ${INFILE1%%.bw}_bigwigcompare_${current_time}.bed >> ${INFILE1%%.bw}_black_list_regions_${current_time}.bed

#awk 'BEGIN {FS="\t"} {if ($5<=0) print $1"\t"$2"\t"$3}' H82_E2F1_082724_bigwigcompare_${current_time}.bed >> H82_E2F1_082724_black_list_regions_${current_time}.bed

#echo "bigwigCompare initiated"
#bigwigCompare -b1 ${INFILE1} -b2 ${INFILE2} -bl ${INFILE1%%.bw}_black_list_regions_${current_time}.bed --operation $ope -p 25 -o ${INFILE1%%.bw}_bigwigcompare2_${current_time}.bw -of bigwig
#bigwigCompare -b1 H82_E2F1_082724.bw -b2 H82_IgG.bw -bl H82_E2F1_082724_black_list_regions_${current_time}.bed --operation $ope -p 25 -o H82_E2F1_082724_bigwigcompare2_${current_time}.bw -of bigwig

#echo "computeMatrix initiated"
#computeMatrix reference-point --referencePoint center -S ${INFILE1}_bigwigcompare2_${current_time}.bw -bl ${INFILE1}_black_list_regions_{current_time}.txt -R promoterFiltered.bed --binSize 1 -b 3000 -a 3000 --skipZeros -o ${INFILE1}_${current_time}.mat.gz -p 25

#computeMatrix reference-point --referencePoint center -S H82_c_MYC_bigwigcompare2_2024_08_27_12_35_48.bw H82_E2F1_082724_bigwigcompare2_${current_time}.bw -R promoterFiltered.bed --binSize 1 -b 3000 -a 3000 --skipZeros -o H82_${current_time}.mat.gz -p 25

#echo "plotProfile initiated"
#plotProfile -m ${INFILE1}_${current_time}.mat.gz -out ${INFILE1}_${current_time}.svg --perGroup --plotTitle "Cut&Run seq signal plot" --refPointLabel promoter
#plotProfile -m H82_${current_time}.mat.gz -out H82_${current_time}.svg --perGroup --plotTitle "Cut&Run seq signal plot" --refPointLabel promoter

echo "done"

source deactivate
