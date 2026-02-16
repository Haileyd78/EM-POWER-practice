.libPaths(c("~/Rlibs", .libPaths()))

library(bigsnpr)
library(bigreadr)
library(dplyr)
library(ggplot2)
library(R.utils)
library(RhpcBLASctl)
library(Cairo)
library(pscl)
library(bigparallelr)

NCORES <- 1   # Number of cores

#Creating directories:
base_dir <- "/gpfs/home4/hdavis/abcd_pgs"
tmp_dir <- Sys.getenv("TMPDIR")


geno_tmp <- file.path(tmp_dir, "bigsnpr")
out_dir      <- file.path(base_dir, "output")
plots_dir    <- file.path(out_dir, "plots")

# Creating bigsnp object of genotypes
bigsnp <- snp_readBed("/gpfs/home4/hdavis/genome.QC.bed",
                      backingfile = file.path(geno_tmp, "genotypes"))

# Loading bigsnp object
obj.bigsnp <- snp_attach(file.path(geno_tmp, "genotypes.rds"))

# Putting genetic data into separate object
G <- obj.bigsnp$genotypes 
G <- snp_fastImputeSimple(G) #- Maybe do this?

# Loading in summary statistics
df_beta <- read.table("/gpfs/home4/hdavis/ADHD.Transformed.meta",header=TRUE)

# Renaming columns for lassosum2
df_beta <- data.frame(
  chr = df_beta$CHR,
  pos = df_beta$BP,
  rsid = df_beta$SNP,
  a1 = df_beta$A1,
  a0 = df_beta$A2,
  beta = df_beta$BETA,
  beta_se = df_beta$SE,
  p = df_beta$P
)

# GWAS effective sample size for binary traits (4 / (1 / n_case + 1 / n_control))
# For quantitative traits, just use the total sample size for `n_eff`.
df_beta$n_eff <- 4 / (1 / 38691  + 1 / 186843 ) 


# Create empty objects for LD loop
ld <- NULL
corr <- NULL

for (chr in 1:22) {
  print(chr)
  # Load saved correlation
  corr0 <- readRDS(paste0("/gpfs/home4/hdavis/abcd_pgs/output/corr_chr/chr", chr, ".rds"))
  
  # Compute LD score
  if (is.null(ld)) {
    ld <- Matrix::colSums(corr0^2)
    # Convert to SFBM
    corr <- as_SFBM(corr0, "/scratch-shared/hdavis/corr_sfbm", compact = TRUE)
  } else {
    ld <- c(ld, Matrix::colSums(corr0^2))
    # Append columns to SFBM
    corr$add_columns(corr0, nrow(corr))
  }
}

# Polygenic score calculation
beta_lassosum2 <- snp_lassosum2(
  corr,
  df_beta,
  ncores = NCORES,
  delta = c(10, 100, 1000, 10000),
  nlambda = 10,
  lambda.min.ratio = 1e-4,
  maxiter = 1000,
  tol = 1e-2
)  

pred_grid2 <- big_prodMat(G, beta_lassosum2, ind.col = df_beta[["_NUM_ID_"]],
                          ncores = NCORES)

params2 <- attr(beta_lassosum2, "grid_param")

# Loading phenotypic data (ADHD)
phenos <- read.table(
  "/projects/einf2700/NeuroIMAGE_genetics/RICOPILI_merged/ricopili_imp_SELECTIVE_UNTAR_4tim/cobg_dir_genome_wide/traits/Mphen_binary.csv",
  header = TRUE,
  sep = ",",      
  stringsAsFactors = FALSE
)

# Loading principal components
pcs <- read.table("/scratch-shared/hdavis/PCA.eigenvec", header=F)

colnames(pcs) <- c("famID", "Subject", paste0("PC",1:10))

pheno<-merge(phenos, pcs, by="famID")


# Models
results <- numeric(ncol(pred_grid2))

for (i in 1:ncol(pred_grid2)) {
  scores <- pred_grid2[, i]
  if (all(is.na(scores))) {
    results[i] <- NA
  } else {
    model <- glm(ADHD_EVER ~ scores + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 
                 + PC9 +PC10, data=pheno, family = "binomial")
    null.model <- glm(ADHD_EVER ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 
                      + PC9 + PC10, data=pheno, family = "binomial")
    results[i] <- 1 - (as.numeric(logLik(model)) / as.numeric(logLik(null.model)))
  }
}    

model_fit <- which.max(results)

CairoPNG(filename = file.path(plots_dir, "R2_hist2.png"),
         width = 800, height = 600)

hist(results, main = "McFadden R² distribution", xlab = "McFadden R²")

dev.off()

saveRDS(list(
  beta = beta_lassosum2,
  pred_grid = pred_grid2,
  params = params2,
  r2_values = results,
  best_model = model_fit
), file.path(out_dir, "lassosum2_results2.rds"))
