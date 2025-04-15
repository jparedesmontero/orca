#!/bin/bash
#SBATCH --job-name=orca_genome_assembly
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --array=0-58
#SBATCH --time=1-00:00:00
#SBATCH --output=logs/bowtie_snpcall.out  # Changed to match your SNP calling output
#SBATCH --mail-user=akmedler@svsu.edu
#SBATCH --mail-type=ALL

# Load necessary modules
module load anaconda3
eval "$(conda shell.bash hook)"
conda activate orca-genomics

module load htslib
module load bowtie2
module load samtools
module load bcftools

# Define sample names based on _1.fastq file names
SAMPLES=($(ls fastq_files/*_1.fastq | sed 's|.*/||' | sed 's/_1.fastq//' | sort))

# Get the sample name corresponding to the SLURM array task ID
SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID]}

echo "Processing $SAMPLE..."

# Align reads to the reference
bowtie2 --very-fast-local -p 16 -x orca_index \
  -1 fastq_files/${SAMPLE}_1.fastq \
  -2 fastq_files/${SAMPLE}_2.fastq \
  -S ${SAMPLE}.sam

# Convert SAM to BAM and sort
echo "Converting and sorting BAM in one step..."
samtools view -@ 16 -bS ${SAMPLE}.sam | samtools sort -@ 16 -m 16G -T /scratch/tmp_sort -o ${SAMPLE}.bam
rm ${SAMPLE}.sam

# Index BAM
echo "Indexing BAM file..."
samtools index ${SAMPLE}.bam

# SNP Calling
echo "Calling SNPs with BCFtools..."
bcftools mpileup -Ou -f reference.fasta ${SAMPLE}.bam | bcftools call -mv -Ob -o ${SAMPLE}.bcf

# Convert BCF to VCF
echo "Converting BCF to VCF..."
bcftools view -Ov -o ${SAMPLE}.vcf ${SAMPLE}.bcf

# Compress and index VCF
echo "Compressing and indexing VCF..."
bgzip -c ${SAMPLE}.vcf > ${SAMPLE}.vcf.gz
bcftools index ${SAMPLE}.vcf.gz

echo "Pipeline completed successfully!"

