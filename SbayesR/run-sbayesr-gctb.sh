#!/bin/bash
#SBATCH --job-name=SBayesR
#SBATCH --output=SBayesR.log
#SBATCH --error=SBayesR.err
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=hailey.davis@radboudumc.nl
#SBATCH --time=12:00:00         
#SBATCH --mem=200G               
#SBATCH --cpus-per-task=2        

# Bayesian updating
# If using mldm make sure you’re running .sh in the same directory as the files
/gpfs/home4/hdavis/gctb_2.04.3_Linux/gctb --sbayes R \
  --mldm /gpfs/home4/hdavis/band_ukb_10k_hm3/band_ukb_10k_hm3/ukb10k.mldm \
  --gwas-summary /gpfs/home4/hdavis/ADHD_sumstats_flipped.ma \
  --pi 0.95,0.02,0.02,0.01 \ # Mixture component starting values
  --gamma 0.0,0.01,0.1,1 \ # Gamma values for variance in genetic effects
  --chain-length 5000 \ # Number of MCMC iterations
  --burn-in 1000 \
  --out-freq 10 \
  --out SBayesR_results

# Scoring 
plink2 \
  --bfile /gpfs/home4/hdavis/genome.QC \
  --score SBayesR_results.snpRes 2 5 8 header \ # Cols for SNP ID, A1, and effect size
  --out sbayesr_scores
