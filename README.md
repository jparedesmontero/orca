# Orca_Final_Project
---
# DONE IN SHARED FOLDER
```
cd /ocean/projects/agr250001p/shared/orcas
```

```
module load anaconda3
conda create -n sra-tools -c bioconda -c conda-forge sra-tools
conda activate sra-tools
```
#Add download.sh script here
#Downlaod orca reference from NCBI
```
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/937/001/465/GCF_937001465.1_mOrcOrc1.1/GCF_937001465.1_mOrcOrc1.1_genomic.fna.gz
unzip  GCF_937001465.1_mOrcOrc1.1_genomic.fna.gz
mv GCF_937001465.1_mOrcOrc1.1_genomic.fna reference.fasta
```
---
# THESE STEPS ARE RUN IN YOUR OCEAN FOLDER
#Data is in the shared folder
#Create symlink to use data from your ocean folder
```
myocean
ln -s /ocean/projects/agr250001p/shared/orcas/fastq_files .
```
#### Confirm your symlink was created
```
ls
```
#use Bowtie2 to assemble genomes
```
vim assembly.slurm
```
- Type `I`

```
#!/bin/bash
#SBATCH --job-name=orca_genome_assembly
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --array=0-58
#SBATCH --time=1-00:00:00
#SBATCH --output=logs/bowtie.out
#SBATCH --mail-user=akmedler@svsu.edu
#SBATCH --mail-type=ALL

#Define samples names
SAMPLES=($(ls fastq_files/*_1.fastq | sed 's|.*/||' | sed 's/_1.fastq//' | sort))
# Get sample name for this task ID
SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID]}

echo "Processing $SAMPLE..."

module load bowtie2
bowtie2 --very-fast-local -p 16 -x orca_index -1 fastq_files/${SAMPLE}_1.fastq -2 fastq_files/${SAMPLE}_2.fastq -S ${SAMPLE}.sam

module load samtools

echo "Converting and Sorting BAM in one step..."
samtools view -@ 16 -bS ${SAMPLE}.sam | samtools sort -@ 16 -m 16G -T /scratch/tmp_sort -o ${SAMPLE}.bam
rm ${SAMPLE}.sam

echo "Indexing BAM file..."
samtools index ${SAMPLE}.bam

echo "Calling SNPs with BCFtools..."
bcftools mpileup -Ou -f reference.fasta ${SAMPLE}.bam | bcftools call -mv -Ob -o ${SAMPLE}.bcf

echo "Converting BCF to VCF..."
bcftools view -Ov -o ${SAMPLE}.vcf ${SAMPLE}.bcf

echo "Compressing and indexing VCF..."
bgzip -c ${SAMPLE}.vcf > ${SAMPLE}.vcf.gz
bcftools index ${SAMPLE}.vcf.gz

echo "Pipeline completed successfully!"
```
- Execute

```
sbatch assembly.slurm
```

---
# Look for selection across genes
- Download gene annotations
```
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/937/001/465/GCF_937001465.1_mOrcOrc1.1/GCF_937001465.1_mOrcOrc1.1_genomic.gff.gz
unzip GCF_937001465.1_mOrcOrc1.1_genomic.gff.gz
mv GCF_937001465.1_mOrcOrc1.1_genomic.gff annotations.gff
```
---



### This project will analyze which genese are under positive selection in orcas whales, and what roles they play in unique traits, such as social behaviors and hunting strategies.

## Hypothesis: Genes under positive selection in orcas are associated with specific social behaviors and hunting strategies, and that ecotypes with similiar behaviors will share common genetic traits. 

## Approach: 
#  1) Sequence alignment
#          - Map sequencing reads to reference orca genome
#  2) SNP calling 
#          - Identify single nucleotide polymorphism (SNPs), insertions, & deletion (INDELs)
#  3) VCF tools/files
#          - Find genetic differences that may contribute to ecotype-specific traits
#  4) Phylogenetic Analysis and PCA 
#          - Visualize genomic variation among ecotypes

# DATA:
#      reference genome -> 
#       samples -> PRJNA531206




##### 1 creating space for fastq files
mkdir fastq
##### 2 
cd fastq
#### 3 moving srr accessions into fastq directory
mv srr_accessions_orca.txt fastq


 #####
 for acc in $(cat srr_accessions_orca.txt); do echo "Downloading $acc..."; fastq-dump --split-files --gzip $acc; done


module load anaconda3

conda create -n sra-tools -c bioconda -c conda-forge sra-tools

conda activate sra-tools

Now you can run fasterq-dump in different CPUs, but first ask the HPC for 16 CPUs for starters:

salloc --cpus-per-task=16 --mem=16G --time=2:00:00

Now run fasterq-dump:

fasterq-dump fastq_files/SRR* -e 8 -O fastq_dir/

Your fastq.gz files will be stored in this folder: fastq_dir
