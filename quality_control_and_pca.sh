## Base data quality control
wc -l /gpfs/home4/hdavis/ADHD_meta_Jan2022_iPSYCH1_iPSYCH2_deCODE_PGC.meta #Checking number of SNPs before QC

# Filter based on INFO
awk 'NR==1 || ($1>=1 && $1<=22 && $8>=0.8)' /gpfs/home4/hdavis/ADHD_meta_Jan2022_iPSYCH1_iPSYCH2_deCODE_PGC.meta > ADHD.QC1.meta

# Filter based on MAF
awk 'NR==1 || ( ($7>=0.01 && $7<=0.99) )' ADHD.QC1.meta > ADHD.QC2.meta

# Remove Duplicate SNPS
awk 'NR==1 || !seen[$2]++' ADHD.QC2.meta > ADHD.QC3.meta

# Remove ambiguous SNPs
awk 'NR==1 || !( ($4=="A" && $5=="T") || ($4=="T" && $5=="A") || ($4=="C" && $5=="G") || ($4=="G" && $5=="C") )' ADHD.QC3.meta > ADHD.QC4.meta

wc -l ADHD.QC4.meta #Checking number of SNPs after QC

## Target data quality control
# Checking number of participants
wc -l /projects/einf2700/NeuroIMAGE_genetics/RICOPILI_merged/ricopili_imp_SELECTIVE_UNTAR_4tim/cobg_dir_genome_wide/adh_amd1_eur_es-postpc2.hg19.ch.fl.bg.fam

# Standard GWAS QC
plink2 --bfile /projects/einf2700/NeuroIMAGE_genetics/RICOPILI_merged/ricopili_imp_SELECTIVE_UNTAR_4tim/cobg_dir_genome_wide/adh_amd1_eur_es-postpc2.hg19.ch.fl.bg \
  --maf 0.01 \
  --hwe 1e-6 \
  --geno 0.01 \
  --mind 0.01 \
  --write-snplist \
  --make-just-fam \
  --out genome.QC
  
awk 'NR==1 || !seen[$2]++' genome.QC.fam > genome.QC1.fam

##########################################################
## PCA of ancestry

# Pruning
plink2 \
  --bfile /projects/einf2700/NeuroIMAGE_genetics/RICOPILI_merged/ricopili_imp_SELECTIVE_UNTAR_4tim/cobg_dir_genome_wide/adh_amd1_eur_es-postpc2.hg19.ch.fl.bg \
  --keep genome.QC.fam \
  --extract genome.QC.snplist \
  --indep-pairwise 200 50 0.25 \
  --out genome.QC

# PCA
plink2 \
  --bfile /gpfs/home4/hdavis/genome.QC \
  --extract /gpfs/home4/hdavis/genome.QC.prune.in \
  --pca 10 \
  --out /scratch-shared/hdavis/PCA


