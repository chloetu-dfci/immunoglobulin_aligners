#!/bin/sh

#SBATCH -p short
#SBATCH -t 0-11:00
#SBATCH --mem=32G
#SBATCH -c 8

# Process_from_BAM_and_FASTQ.sh

# sbatch Process_from_BAM_and_FASTQ.sh /n/scratch3/users/c/cht269/Marwan_spon_data/13JH_11_5_16 13JH-11-5-16 /n/scratch3/users/c/cht269/Marwan_spon_data/13JH_4_7_2022 13JH_4_7_2022 13JH /n/scratch3/users/c/cht269/Marwan_spon_data/13JH_Analysis
# sbatch Process_from_BAM_and_FASTQ.sh /n/scratch3/users/c/cht269/Marwan_spon_data/110BJ_3_1_13 110BJ-3-1-13 /n/scratch3/users/c/cht269/Marwan_spon_data/110BJ_4_6_2022 110BJ_4_6_2022 110BJ /n/scratch3/users/c/cht269/Marwan_spon_data/110BJ_Analysis
# sbatch Process_from_BAM_and_FASTQ.sh /n/scratch3/users/c/cht269/Marwan_spon_data/PC863_25_4_16 PC863-25-4-16 /n/scratch3/users/c/cht269/Marwan_spon_data/PC863_4_21_2022 PC863_4_21_2022 PC863 /n/scratch3/users/c/cht269/Marwan_spon_data/PC863_Analysis
# sbatch Process_from_BAM_and_FASTQ.sh /n/scratch3/users/c/cht269/Marwan_spon_data/CDs106_25_4_16 CDs106-25-4-16 /n/scratch3/users/c/cht269/Marwan_spon_data/CDs106_4_5_2022 Cds106_4_5_2022 CDs106 /n/scratch3/users/c/cht269/Marwan_spon_data/CDS106_Analysis


module load gcc/6.2.0
module load bwa/0.7.17
module load samtools/1.15.1
module load java/jdk-11.0.11
module load picard


echo "fastq directory = " $1 #/n/scratch3/users/c/cht269/Marwan_spon_data/13JH_11_5_16
echo "fastq name = " $2 #13JH-11-5-16
echo "bam2 directory = " $3 #/n/scratch3/users/c/cht269/Marwan_spon_data/13JH_4_7_2022
echo "bam2 name = " $4 #13JH_4_7_2022
echo "sample name = " $5 #13JH
echo "output directory = " $6 #/n/scratch3/users/c/cht269/Marwan_spon_data/13JH_Analysis

cd $6

# Concatenate 
cat $1/*R1_001.fastq.gz > $1/$2_R1.fastq.gz
cat $1/*R2_001.fastq.gz > $1/$2_R2.fastq.gz

echo $1/$2_R1.fastq.gz

# Sort and index the old fastq files
bwa mem -t 8 /n/scratch3/users/c/cht269/Homo_sapiens_assembly19.fasta $1/$2_R1.fastq.gz $1/$2_R2.fastq.gz | samtools sort -@ 8 -o $1/$2_WES_Tumor.bam 

samtools index $1/$2_WES_Tumor.bam 

# Merge the two sorted BAMs together
samtools merge -@ 8 -o  $1/new_$5_WES_Tumor.bam $1/$2.bam $3/$4.bam

# Fixing “missing read group error”
java -jar $PICARD/picard.jar AddOrReplaceReadGroups \
      I=$1/new_$5_WES_Tumor.bam  \
      O=$1/fixed_$5_WES_Tumor.bam  \
      RGID=$5 \
      RGLB=lib1 \
      RGPL=ILLUMINA \
      RGPU=unit1 \
      RGSM=$5

# Fixmate
 java -jar $PICARD/picard.jar FixMateInformation \
       I=$1/fixed_$5_WES_Tumor.bam \
       O=$1/new_$5_WES_Tumor.bam \
       ADD_MATE_CIGAR=true

# Index
samtools index  $1/new_$5_WES_Tumor.bam

# Extract unmapped reads
samtools view -@ 8 -b -f 4 -f 0x1 -f 0x2 $1/new_$5_WES_Tumor.bam > $1/new_$5_WES_Tumor_unmapped_reads.bam

# Extract reads mapping to chromosome 2, 14 and 22
samtools view -@ 8 -b -h -P -f 0x1 -f 0x2 -L /n/scratch3/users/c/cht269/gencode_v19_whole_chr2_14_22.bed $1/new_$5_WES_Tumor.bam > $1/new_$5_WES_Tumor_Ig_chr_reads.bam

# Merge BAMs
samtools merge -@ 8 -o $1/new_$5_WES_Tumor_subset.bam new_$5_WES_Tumor_Ig_chr_reads.bam $1/new_$5_WES_Tumor_unmapped_reads.bam

# Validate SAM file
java -jar $PICARD/picard.jar ValidateSamFile I=$1/new_$5_WES_Tumor_subset.bam MODE=SUMMARY

# Index
samtools index $1/new_$5_WES_Tumor_subset.bam

# Run MixCR
mixcr analyze exome-full-length \
    --species hsa \
    --threads 8 \
        $1/new_$5_WES_Tumor_subset.bam new_$5_WES_Tumor_subset_MixCR_Results

mkdir new_$5_WES_Tumor_subset_MixCR
mv new_$5_WES_Tumor_subset_MixCR* new_$5_WES_Tumor_subset_MixCR

# Run IgCaller
source /n/scratch3/users/c/cht269/py3_IgCaller/bin/activate
python3 /n/scratch3/users/c/cht269/IgCaller/IgCaller -I /n/scratch3/users/c/cht269/IgCaller/IgCaller_reference_files/ -V hg19 -C ensembl -T $1/new_$5_WES_Tumor_subset.bam -R /n/scratch3/users/c/cht269/Homo_sapiens_assembly19.fasta -seq wes -@ 8 -d 4 -ad 4 -o . -bq 20 -mq 10

