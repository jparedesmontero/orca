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




##### 1
mkdir fastq
##### 2 
cd fastq
#### 3
mv srr_accessions_orca.txt fastq/
#### 4
for acc in $(cat srr_accessions_orca.txt); do
    fasterq-dump $acc --split-files -O .
done
#### 5 
for acc in $(cat srr_accessions_orca.txt); do     echo "Downloading $acc...";     fastq-dump --split-files --gzip $acc; done
