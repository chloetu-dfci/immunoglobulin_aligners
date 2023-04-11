#!/bin/sh

#SBATCH -p short
#SBATCH -t 0-6:00
#SBATCH --mem=32G
#SBATCH -c 8

module load gcc/6.2.0
module load samtools/1.15.1
module load java/jdk-11.0.11
module load picard

# Process_from_TWO_BAM.sh

# sbatch Process_from_2_BAM.sh /n/scratch3/users/c/cht269/Marwan_spon_data/10JB_5_11_2017 10JB_5_11_2017 /n/scratch3/users/c/cht269/Marwan_prev_CLL_data/10JB_Analysis old_10JB_WES_Tumor 10JB /n/scratch3/users/c/cht269/Marwan_spon_data/10JB_Analysis
# sbatch Process_from_2_BAM.sh /n/scratch3/users/c/cht269/Marwan_spon_data/OJ802_4_21_2022 OJ802_4_21_2022 /n/scratch3/users/c/cht269/Marwan_prev_CLL_data/OJ802_Analysis old_OJ802_WES_Tumor OJ802 /n/scratch3/users/c/cht269/Marwan_spon_data/OJ802_Analysis


echo "bam1 directory = " $1 #/n/scratch3/users/c/cht269/Marwan_spon_data/10JB_5_11_2017
echo "bam1 name = " $2 #13JH-11-5-16
echo "bam2 directory = " $3 #/n/scratch3/users/c/cht269/Marwan_prev_CLL_data/10JB_Analysis
echo "bam2 name = " $4 #old_10JB_WES_Tumor
echo "sample name = " $5 #10JB
echo "output directory = " $6 #/n/scratch3/users/c/cht269/Marwan_spon_data/10JB_Analysis

cd $6

# Merge the two sorted BAMs together
samtools merge -@ 8 -o new_$5_WES_Tumor.bam $1/$2.bam $3/$4.bam

# Fixing “missing read group error”
java -jar $PICARD/picard.jar AddOrReplaceReadGroups \
      I=new_$5_WES_Tumor.bam  \
      O=fixed_$5_WES_Tumor.bam  \
      RGID=$5 \
      RGLB=lib1 \
      RGPL=ILLUMINA \
      RGPU= unit1 \
      RGSM= $5

# Fixmate
 java -jar $PICARD/picard.jar FixMateInformation \
       I=$1/fixed_$5_WES_Tumor.bam \
       O=$1/new_$5_WES_Tumor.bam \
       ADD_MATE_CIGAR=true

# Index
samtools index new_$5_WES_Tumor.bam

# Extract unmapped reads
samtools view -b -f 4 -f 0x1 -f 0x2 new_$5_WES_Tumor.bam > new_$5_WES_Tumor_unmapped_reads.bam

# Extract reads mapping to chromosome 2, 14 and 22
samtools view -b -h -P -f 0x1 -f 0x2 -L /n/scratch3/users/c/cht269/gencode_v19_whole_chr2_14_22.bed new_$5_WES_Tumor.bam > new_$5_WES_Tumor_Ig_chr_reads.bam

# Merge BAMs
samtools merge -@ 8 -o new_$5_WES_Tumor_subset.bam new_$5_WES_Tumor_Ig_chr_reads.bam new_$5_WES_Tumor_unmapped_reads.bam

# Validate SAM file
java -jar $PICARD/picard.jar ValidateSamFile I=new_$5_WES_Tumor_subset.bam MODE=SUMMARY

# Index
samtools index new_$5_WES_Tumor_subset.bam

# Run MixCR
mixcr analyze exome-full-length \
    --species hsa \
    --threads 8 \
        new_$5_WES_Tumor_subset.bam new_$5_WES_Tumor_subset_MixCR_Results

mkdir new_$5_WES_Tumor_subset_MixCR
mv new_$5_WES_Tumor_subset_MixCR* new_$5_WES_Tumor_subset_MixCR

# Run IgCaller
source /n/scratch3/users/c/cht269/py3_IgCaller/bin/activate
python3 /n/scratch3/users/c/cht269/IgCaller/IgCaller -I /n/scratch3/users/c/cht269/IgCaller/IgCaller_reference_files/ -V hg19 -C ensembl -T new_$5_WES_Tumor_subset.bam -R /n/scratch3/users/c/cht269/Homo_sapiens_assembly19.fasta -seq wes -@ 8 -d 4 -ad 4 -o . -bq 20 -mq 10


