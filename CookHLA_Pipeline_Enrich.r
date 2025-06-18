#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -V
#$ -r yes
#$ -l mem_requested=16G
#$ -N LBIO_SNP_HLA
#$ -q short.q
#$ -e /directflow/SCCGGroupShare/projects/jenli3/EnRICH_SNP_Data/error.stderr
#$ -o /directflow/SCCGGroupShare/projects/jenli3/EnRICH_SNP_Data/output.stdout
#$ -m ae
#$ -M j.lim@garvan.org.au

# -----------------------------------------
# CookHLA Pipeline - LBIO EnRICH Cohort
# -----------------------------------------

# Move to CookHLA working directory
cd /directflow/SCCGGroupShare/projects/jenli3/software/CookHLA

# Activate CookHLA conda environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate CookHLA

# -----------------------------------------
# Step 1: Convert VCF to PLINK format
# -----------------------------------------
plink --vcf /directflow/SCCGGroupShare/projects/SNP_Genotype_Processing/Imputation/13_14_15_Powell_011024/results/vcf_all_merged/imputed_hg38.vcf.gz \
      --make-bed \
      --out /directflow/SCCGGroupShare/projects/jenli3/software/CookHLA/imputed_hg38

# -----------------------------------------
# Step 2: Create Genetic Map
# -----------------------------------------
python -m MakeGeneticMap \
  -i /directflow/SCCGGroupShare/projects/jenli3/software/CookHLA/Enrich_SNP \
  -hg 38 \
  -ref /directflow/SCCGGroupShare/projects/jenli3/software/CookHLA/1000G_REF/1000G_REF.AMR.chr6.hg18.29mb-34mb.inT1DGC \
  -o /directflow/SCCGGroupShare/projects/jenli3/software/CookHLA/output/genetic_map/1000G

# -----------------------------------------
# Step 3: Impute HLA Haplotypes
# -----------------------------------------
python CookHLA.py \
  -i /directflow/SCCGGroupShare/projects/jenli3/software/CookHLA/Enrich_SNP \
  -hg 38 \
  -o /directflow/SCCGGroupShare/projects/jenli3/software/CookHLA/output/HLA_Enrich \
  -ref /directflow/SCCGGroupShare/projects/jenli3/software/CookHLA/1000G_REF/1000G_REF.AMR.chr6.hg18.29mb-34mb.inT1DGC \
  -gm /directflow/SCCGGroupShare/projects/jenli3/software/CookHLA/output/genetic_map/1000G.mach_step.avg.clpsB \
  -ae /directflow/SCCGGroupShare/projects/jenli3/software/CookHLA/output/genetic_map/1000G.aver.erate
