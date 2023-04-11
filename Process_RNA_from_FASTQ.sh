#!/bin/sh

#SBATCH -p short
#SBATCH -t 0-5:00
#SBATCH --mem=28G
#SBATCH -c 8

# Process_RNA_from_FASTQ.sh
# Align, Index and run MixCR on RNA BAMs

module load conda3/latest
source activate conda_env
module load gcc/6.2.0
module load bwa/0.7.17
module load samtools/1.15.1
module load java/jdk-11.0.11
module load picard

echo "bam directory = " $1 #/n/scratch3/users/c/cht269/Marwan_spon_data/110BJ_RNA
echo "sample name = " $2 #110BJ
echo "output directory = " $3 #/n/scratch3/users/c/cht269/Marwan_spon_data/110BJ_Analysis

cd $1

cat *R1_001.fastq.gz > $2_R1.fastq.gz
cat *R2_001.fastq.gz > $2_R2.fastq.gz

bwa mem -t 8 /n/scratch3/users/c/cht269/Homo_sapiens_assembly19.fasta $2_R1.fastq.gz $2_R2.fastq.gz | samtools sort -@ 8 -o $2_RNA_Tumor.bam

# Fixing “missing read group error”
java -jar $PICARD/picard.jar AddOrReplaceReadGroups \
      I=$2_RNA_Tumor.bam  \
      O=fixed_$2_RNA_Tumor.bam  \
      RGID=$2 \
      RGLB=lib1 \
      RGPL=ILLUMINA \
      RGPU= unit1 \
      RGSM= $2

rm $2_RNA_Tumor.bam $2_RNA_Tumor.bam.bai

mv fixed_$2_RNA_Tumor.bam $2_RNA_Tumor.bam

java -jar $PICARD/picard.jar ValidateSamFile I=$2_RNA_Tumor.bam MODE=SUMMARY

samtools index $2_RNA_Tumor.bam

mixcr analyze rnaseq-full-length \
    --species hsa \
    --threads 8 \
        $2_RNA_Tumor.bam $2_RNA_Tumor_MixCR_Results

mkdir $3/$2_RNA_Tumor_MixCR
mv $2_RNA_Tumor_MixCR_Results* $3/$2_RNA_Tumor_MixCR