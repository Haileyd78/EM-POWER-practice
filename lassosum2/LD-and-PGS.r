.libPaths(c("~/Rlibs", .libPaths()))

library(bigsnpr)
library(bigreadr)
library(dplyr)
library(ggplot2)
library(R.utils)
library(Cairo)
library(pscl)
library(bigparallelr)

NCORES <- 1  

# Creating directories
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

# Creating SNP map data frame
map <- dplyr::transmute(obj.bigsnp$map,
                        chr = chromosome, pos = physical.pos,
                        a0 = allele2, a1 = allele1)

# Loading in gwas data
sumstats <- read.table("/gpfs/home4/hdavis/ADHD.Transformed.meta",header=TRUE)

Nca <- 38691 #N for cases
Nco <- 186843 #N for controls

# organizing sumstats and renaming columns
sumstats <- data.frame(
  chr = sumstats$CHR,
  pos = sumstats$BP,
  rsid = sumstats$SNP,
  a1 = sumstats$A1,
  a0 = sumstats$A2,
  beta = sumstats$BETA,
  beta_se = sumstats$SE,
  freq = (sumstats$FRQ_A_38691 * Nca + sumstats$FRQ_U_186843 * Nco) / (Nca + Nco),
  p = sumstats$P
)

# Calculating effective sample size
sumstats$n_eff <- 4 / (1 / 38691  + 1 / 186843 ) 

df_beta <- snp_match(sumstats, map, return_flip_and_rev = TRUE) %>% 
  mutate(freq = ifelse(`_REV_`, 1 - freq, freq), 
         `_REV_` = NULL, `_FLIP_`= NULL) 

head(df_beta)

##Calculating LD
# Loop to calculate LD correlation matrix by chromosome
for (chr in 1:22) {  
  
  print(chr)
  
  corr0 <- runonce::save_run({
    
    ## indices in 'sumstats'
    ind.chr <- which(df_beta$chr == chr)
    ## indices in 'G'
    ind.chr2 <- df_beta$`_NUM_ID_`[ind.chr]
    
    #Genetic positions (in cM)
    POS2 <- snp_asGeneticPos(
      map$chr[ind.chr2],
      map$pos[ind.chr2],
      dir = "/scratch-shared/hdavis/genetic_maps"
    )
    
    #Computing the banded correlation matrix in sparse matrix format
    snp_cor(G, ind.col = ind.chr2, size = 3 / 1000, infos.pos = POS2, 
            ncores = NCORES)
    
  }, file = paste0("/scratch-shared/hdavis/cor", chr, ".rds"))
  # transform to SFBM (on-disk format) on the fly
  if (chr == 1) {
    ld <- Matrix::colSums(corr0^2)
    corr <- as_SFBM(corr0, "/scratch-shared/hdavis/corr", compact = TRUE)
  } else {
    ld <- c(ld, Matrix::colSums(corr0^2))
    corr$add_columns(corr0, nrow(corr))
  }
}

##Calculating polygenic scores
# Running lassosum2
beta_lassosum2 <- snp_lassosum2(
  corr,
  df_beta,
  ncores = NCORES,
  delta = c(1, 10, 100, 1000, 10000),
  nlambda = 10,
  lambda.min.ratio = 1e-4,
  maxiter = 1000,
  tol = 1e-2
)  

# Calculating polygenic scores based on lassosum2
pred_grid2 <- big_prodMat(G, beta_lassosum2, ind.col = df_beta[["_NUM_ID_"]],
                          ncores = NCORES)

params2 <- attr(beta_lassosum2, "grid_param")

##Testing predictiveness of models
# Loading phenotype data (ADHD)
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

results <- numeric(ncol(pred_grid2))

# Regression of scores on phenotype
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

# Saving regression results and model parameters
saveRDS(list(
  beta = beta_lassosum2,
  pred_grid = pred_grid2,
  params = params2,
  r2_values = results,
  best_model = model_fit
), file.path(out_dir, "lassosum2_results2.rds"))

ind_match <- match(obj.bigsnp$fam$family.ID, phenos$famID)  

# Reorder phenos
phenos_ordered <- phenos[ind_match, ]

# Confirm alignment
all(obj.bigsnp$fam$family.ID == phenos_ordered$famID)  # should be TRUE

# Extract best PGS
best_pgs <- pred_grid2[, model_fit]

# Combine with IDs 
best_pgs_df <- data.frame(
  FamID = phenos_ordered$famID, 
  PGS   = best_pgs
)

# Saving best PGS model
saveRDS(best_pgs_df, file = file.path(out_dir, "best_pgs_aligned.rds"))
