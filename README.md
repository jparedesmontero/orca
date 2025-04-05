# Orca_Final_Project

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
#      reference genome -> PRJNA531206
#      NCBI samples -> 




##### 1 creating space for fastq files
mkdir fastq
##### 2 
cd fastq
#### 3 moving srr accessions into fastq directory
mv srr_accessions_orca.txt fastq

module load anaconda3

conda create -n sra-tools -c bioconda -c conda-forge sra-tools
conda activate sra-tools

Now you can run fasterq-dump in different CPUs, but first ask the HPC for 16 CPUs for starters:

salloc --cpus-per-task=16 --mem=16G --time=2:00:00

Now run fasterq-dump:

fasterq-dump SRR* -e 8 -O fastq_dir/

Your fastq.gz files will be stored in this folder: fastq_dir
