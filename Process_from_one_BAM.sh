#!/bin/sh

#SBATCH -p short
#SBATCH -t 0-6:00
#SBATCH --mem=32G
#SBATCH -c 8

module load gcc/6.2.0
module load samtools/1.15.1
module load java/jdk-11.0.11
module load picard

# Process_from_one_BAM.sh

# sbatch Process_from_2_BAM.sh /n/scratch3/users/c/cht269/Marwan_spon_data/10JB_5_11_2017 10JB_5_11_2017 /n/scratch3/users/c/cht269/Marwan_prev_CLL_data/10JB_Analysis old_10JB_WES_Tumor 10JB /n/scratch3/users/c/cht269/Marwan_spon_data/10JB_Analysis
# sbatch Process_from_2_BAM.sh /n/scratch3/users/c/cht269/Marwan_spon_data/OJ802_4_21_2022 OJ802_4_21_2022 /n/scratch3/users/c/cht269/Marwan_prev_CLL_data/OJ802_Analysis old_OJ802_WES_Tumor OJ802 /n/scratch3/users/c/cht269/Marwan_spon_data/OJ802_Analysis


echo "bam1 directory = " $1 #/n/scratch3/users/c/cht269/Marwan_spon_data/10JB_5_11_2017
echo "bam1 name = " $2 #13JH-11-5-16
echo "sample name = " $3 #10JB
echo "output directory = " $4 #/n/scratch3/users/c/cht269/Marwan_spon_data/10JB_Analysis

cd $4

# Extract unmapped reads
samtools view -b -f 4 -f 0x1 -f 0x2 $1/$2.bam > $3_WES_Tumor_unmapped_reads.bam

# Extract reads mapping to chromosome 2, 14 and 22
samtools view -b -h -P -f 0x1 -f 0x2 -L /n/scratch3/users/c/cht269/gencode_v19_whole_chr2_14_22.bed $1/$2.bam > $3_WES_Tumor_Ig_chr_reads.bam

# Merge BAMs
samtools merge -@ 8 -o $3_WES_Tumor_subset.bam $3_WES_Tumor_Ig_chr_reads.bam $3_WES_Tumor_unmapped_reads.bam

# Validate SAM file
java -jar $PICARD/picard.jar ValidateSamFile I=$3_WES_Tumor_subset.bam MODE=SUMMARY

# Index
samtools index $3_WES_Tumor_subset.bam

# Run MixCR
mixcr analyze exome-full-length \
    --species hsa \
    --threads 8 \
        $3_WES_Tumor_subset.bam $3_WES_Tumor_subset_MixCR_Results

mkdir $3_WES_Tumor_subset_MixCR
mv $3_WES_Tumor_subset_MixCR* $3_WES_Tumor_subset_MixCR

# Run IgCaller
source /n/scratch3/users/c/cht269/py3_IgCaller/bin/activate
python3 /n/scratch3/users/c/cht269/IgCaller/IgCaller -I /n/scratch3/users/c/cht269/IgCaller/IgCaller_reference_files/ -V hg19 -C ensembl -T $3_WES_Tumor_subset.bam -R /n/scratch3/users/c/cht269/Homo_sapiens_assembly19.fasta -seq wes -@ 8 -d 4 -ad 4 -o . -bq 20 -mq 10





