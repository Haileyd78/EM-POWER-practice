### Shell
chmod u+x /gpfs/home4/hdavis/PRSice_linux/PRSice_linux #Giving linux permission to use PRSice

# Performing polygenic scoring using PRScise
/gpfs/home4/hdavis/PRSice_linux/PRSice_linux \ #PRScise application location
  --base /gpfs/home4/hdavis/ADHD.Transformed.meta \ #Summary statistics file
  --snp SNP --chr CHR --bp BP --A1 A1 --A2 A2 --stat BETA --pvalue P \ #Mapping summary statistic columns
  --target /projects/einf2700/NeuroIMAGE_genetics/RICOPILI_merged/ricopili_imp_SELECTIVE_UNTAR_4tim/cobg_dir_genome_wide/adh_amd1_eur_es-postpc2.hg19.ch.fl.bg \ #loading target data
  --binary-target T \ #Specify if outcome phenotype is binary or not (T or F)
  --no-regress \ #No regression 
  --fastscore \ 
  --all-score \
  --clump-r2 0.2 \ #Clumping at 0.2
  --extract /gpfs/home4/hdavis/PRS.valid \ #Extracting only valid SNPs
  --out /gpfs/home4/hdavis/PRS
