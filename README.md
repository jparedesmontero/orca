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
gunzip  GCF_937001465.1_mOrcOrc1.1_genomic.fna.gz
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
# Use Bowtie2 to assemble genomes
- Build bowtie index
```
bowtie2-build reference.fasta orca_index
```
- Create script to run bowtie
```
vim assembly.slurm
```
- Type `I`

```
#!/bin/bash
#SBATCH --job-name=orca_genome_assembly
#SBATCH --cpus-per-task=16
#SBATCH --mem=40G
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
module load bcftools
module load htslib

echo "Converting and Sorting BAM in one step..."
samtools view -@ 8 -bS ${SAMPLE}.sam | samtools sort -@ 8 -m 4G -T /scratch/tmp_sort -o ${SAMPLE}.bam
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
gunzip GCF_937001465.1_mOrcOrc1.1_genomic.gff.gz
mv GCF_937001465.1_mOrcOrc1.1_genomic.gff annotations.gff
```
---
# Index VCF files 
module load bcftools
module load htslib
tabix -f -p vcf *.vcf.gz

# Make sure you have the metadata.tsv file in your current directory (file can be downloaded from this repository
# Make sure you have working vcf files (not empty)
# Make sure you have the annotations.gff file

# Merge vcf files
```
vi merge_vcfs_from_csv.sh
```
- Type i and paste:
```
#!/bin/bash

module load bcftools

META="metadata.tsv"  # Updated filename if changed

# Extract unique ecotypes from the cleaned metadata
ecotypes=$(awk -F',' 'NR > 1 { print $4 }' "$META" | sort | uniq)

for ecotype in $ecotypes; do
    echo "Merging VCFs for ecotype: $ecotype"

    # All ecotype labels are now already filename-safe
    label="$ecotype"

    # Get all sample IDs for this ecotype
    sample_ids=$(awk -F',' -v e="$ecotype" 'NR > 1 && $4 == e { print $1 }' "$META")

    vcf_list=""
    for id in $sample_ids; do
        vcf_file="${id}.vcf.gz"
        if [[ -f "$vcf_file" ]]; then
            vcf_list="$vcf_list $vcf_file"
        else
            echo "Skipping: $vcf_file not found"
        fi
    done

    if [[ -z "$vcf_list" ]]; then
        echo "No VCFs found for ecotype: $ecotype"
        continue
    fi

    echo "Merging into: ${label}.merged.vcf.gz"
    bcftools merge -m all -Oz -o "${label}.merged.vcf.gz" $vcf_list
    bcftools index "${label}.merged.vcf.gz"
done

echo "All ecotype-based merges completed."

```
# Annotate vcfs
```
module load bedtools/2.30.0

for vcf in *.merged.vcf.gz; do
  prefix=$(basename "$vcf" .merged.vcf.gz)
  bedtools intersect -header -wa -wb \
    -a "$vcf" \
    -b annotations.gff > ${prefix}.annotated.tsv
done
```

# Obtain Tajima's D per gene per merged file
- Creating bed file
```
awk '$3=="gene" { print $1"\t"$4-1"\t"$5"\t"$9 }' annotations.gff > genes.bed
```
- Calculating pi diversity per 1000 bp windows
```
vi pi.sh
```
- Type I and paste
```
mkdir -p pi_window_output
module load vcftools/0.1.16
module load htslib
for vcf in *.merged.vcf.gz; do
    base=$(basename "$vcf" .merged.vcf.gz)
    echo "Calculating nucleotide diversity for $base in 1000bp windows..."
    vcftools --gzvcf "$vcf" \
             --window-pi 1000 \
             --window-pi-step 1000 \
             --out pi_window_output/"$base"
done
```
- Calculating Tajima's D
```
vi tajd.sh
```
- Type I and paste
```
#!/bin/bash

mkdir -p tajimaD_output
module load vcftools/0.1.16
module load htslib
for vcf in *.merged.vcf.gz; do
    base=$(basename "$vcf" .merged.vcf.gz)
    echo "Calculating Tajima's D for $base in 1000 bp windows..."

    vcftools --gzvcf "$vcf" \
             --TajimaD 1000 \
             --out tajimaD_output/"$base"
done
```



# make BED of gene region
# extract only gene lines, convert to 0‑based BED
grep -P "\tgene\t" annotations.gff \
  | awk 'BEGIN{FS=OFS="\t"}{
      # parse ID from the 9th column
      split($9,a,";"); 
      id=a[1]; gsub("ID=","",id);
      print $1, $4-1, $5, id
    }' > genes.bed
# List sample populations
# one sample name per line, matching the VCF header
pop1.txt  
pop2.txt

# Install dependencies 
```
pip install scikit-allel numpy pandas
```

# The script 
```
#!/usr/bin/env python3
import allel
import numpy as np
import pandas as pd
import sys
```

# --- User parameters ---
```
vcf_path   = 'input.vcf.gz'
genes_bed  = 'genes.bed'
pop1_file  = 'pop1.txt'
pop2_file  = 'pop2.txt'
out_csv    = 'gene_fst.csv'
```
# -----------------------

# Load sample lists
pop1 = [s.strip() for s in open(pop1_file) if s.strip()]
pop2 = [s.strip() for s in open(pop2_file) if s.strip()]

# Read gene coordinates
```
genes = pd.read_csv(genes_bed, sep='\t',
                    names=['chrom','start','end','gene'],
                    dtype={'chrom':str, 'start':int, 'end':int, 'gene':str})
```
# Read VCF (only once!)
```
print("Loading VCF…")
callset = allel.read_vcf(vcf_path,
                          fields=['samples','calldata/GT','variants/CHROM','variants/POS'])
samples   = callset['samples']
genos     = allel.GenotypeArray(callset['calldata/GT'])
chroms    = callset['variants/CHROM']
positions = callset['variants/POS']
```

# Map sample names → indices
pop1_idx = [i for i,s in enumerate(samples) if s in pop1]
pop2_idx = [i for i,s in enumerate(samples) if s in pop2]
if not pop1_idx or not pop2_idx:
    sys.exit("Error: could not match any samples in pop1 or pop2.")

# Prepare output
results = []

print("Computing FST per gene…")
for _, row in genes.iterrows():
    chrom, start, end, gene_id = row
    # mask SNPs in this gene
    m = (chroms == chrom) & (positions >= start) & (positions <= end)
    if m.sum() == 0:
        fst_val = np.nan
    else:
        subgeno = genos[m]
        # allele counts per subpop
        ac1 = subgeno.count_alleles(subpop=pop1_idx)
        ac2 = subgeno.count_alleles(subpop=pop2_idx)
        # compute per-site numerator & denominator
        num, den = allel.weir_cockerham_fst(ac1, ac2)
        # weighted FST = sum(num) / sum(den)
        fst_val = np.sum(num) / np.sum(den) if den.sum() > 0 else np.nan
    results.append((gene_id, fst_val))

# Write out
df = pd.DataFrame(results, columns=['gene','fst'])
df.to_csv(out_csv, index=False)
print(f"Done → {out_csv}")
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






