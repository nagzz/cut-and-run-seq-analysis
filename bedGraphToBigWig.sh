#PBS -l walltime=10:00:00
#PBS -l nodes=1:ppn=1
#PBS -A PAS1348
#PBS -t 1 
#PBS -V
echo $PWD
sh ~/.bash_profile
source activate bdgtobw
cd /fs/ess/PAS1348/nagesh/cutandrun_hannahs/output/mapping/bam/rmdups/macs3

#source ~/.bash_profile
INFILE=( *_treat_pileup.bdg )
i="${INFILE}"
echo "$i"
#module load minico
#bedGraphToBigWig "$i" hg38.chrom.sizes "${i%_treat_pileup.bdg}.bw"
bedGraphToBigWig H82_E2F1_082724_treat_pileup.bdg hg38.chrom.sizes H82_E2F1_082724.bw

source deactivate
