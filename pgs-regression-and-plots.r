library(dplyr)
library(ggplot2)

#Reading in pcs and phenotype
pcs <- read.table("/scratch-shared/hdavis/PCA.eigenvec", header=F) #PCAs of ancestry
phenotype <- read.csv("/projects/einf2700/NeuroIMAGE_genetics/RICOPILI_merged/ricopili_imp_SELECTIVE_UNTAR_4tim/cobg_dir_genome_wide/traits/Mphen_binary.csv") #phenotype data

# Reading in PGSs
lasso<-readRDS("/gpfs/home4/hdavis/abcd_pgs/output/best_pgs_aligned.rds")
bayes<- read.table("/gpfs/home4/hdavis/band_ukb_10k_hm3/band_ukb_10k_hm3/sbayesr_scores.sscore")
prspca<-readRDS("/gpfs/home4/hdavis/ADHD_PRS_PCA.rds")

# Renaming columns
colnames(pcs) <- c("famID", "Subject", paste0("PC",1:10))
colnames(lasso)[colnames(lasso) == "FamID"] <- "famID"
colnames(bayes)[colnames(bayes) == "V1"] <- "famID"
colnames(prspca)[colnames(prspca) == "FID"] <- "famID"

colnames(lasso)[colnames(lasso) == "PGS"] <- "lasso_pgs"
colnames(bayes)[colnames(bayes) == "V6"] <- "bayes_pgs"
colnames(prspca)[colnames(prspca) == "ADHD.prs.pc"] <- "prspca_pgs"
colnames(prspca)[colnames(prspca) == "Pt_0.2"] <- "pt_pgs"

#Merging
pheno <- phenotype[, c("famID", "ADHD_EVER")]
lasso_pgs <- lasso[, c("famID", "lasso_pgs")]
bayes_pgs <- bayes[, c("famID", "bayes_pgs")]
prspca_pgs <- prspca[, c("famID", "prspca_pgs", "pt_pgs")] 

merged_pgs <- pcs %>%
  left_join(pheno, by = "famID") %>%
  left_join(lasso_pgs, by = "famID") %>%
  left_join(bayes_pgs, by = "famID") %>%
  left_join(prspca_pgs, by = "famID")

#Calculating R2
  prs_glm <- function(data, outcome, prs.col, covars = paste0("PC", 1:10)) {
    # data: data frame with phenotype, PRS, and covariates
    # outcome: string, name of binary phenotype
    # prs.col: string, name of PRS column
    # covars: vector of covariates (default PC1-PC10)
    
    # Build formulas
    null.formula <- as.formula(
      paste(outcome, "~", paste(covars, collapse = " + "))
    )
    prs.formula <- as.formula(
      paste(outcome, "~", prs.col, "+", paste(covars, collapse = " + "))
    )
    
    # Fit models
    null.model <- glm(null.formula, data = data, family = binomial)
    model <- glm(prs.formula, data = data, family = binomial)
    
    # Pseudo R2 (McFadden's)
    r2 <- 1 - (as.numeric(logLik(model)) / as.numeric(logLik(null.model)))
    
    # Extract PRS coefficient
    coef <- summary(model)$coeff[prs.col, ]
    beta <- as.numeric(coef[1])
    se <- as.numeric(coef[2])
    p <- as.numeric(coef[4])
    
    # Return as data frame
    return(data.frame(PRS = prs.col, R2 = r2, BETA = beta, SE = se, P = p))
  }
  
  pgs_cols <- c("lasso_pgs", "bayes_pgs", "prspca_pgs", "pt_pgs")
  
  prs_results <- do.call(rbind, lapply(pgs_cols, function(x) {
    prs_glm(data = merged_pgs, outcome = "ADHD_EVER", prs.col = x)
  }))
  
  print(prs_results)
  

# Model calculation 
null.model <- glm(ADHD_EVER ~  PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = merged_pgs, family = binomial)
model <- glm(ADHD_EVER ~ pt_pgs + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = merged_pgs, family = binomial)

# Extracting PRS-PCA association results
prs.result <- NULL
prs.r2<- 1 - (as.numeric(logLik(model)) / as.numeric(logLik(null.model)))
prs.coef <- summary(model)$coeff["pt_pgs",]
prs.beta <- as.numeric(prs.coef[1])
prs.se <- as.numeric(prs.coef[2])
prs.p <- as.numeric(prs.coef[4])


# Storing the results
prs.result <- rbind(prs.result, data.frame(R2=prs.r2, P=prs.p, BETA=prs.beta,SE=prs.se))

print(prs.result) 

#Corr test
# Select only the PGS columns
pgs_matrix <- merged_pgs[, c("lasso_pgs", "bayes_pgs", "prspca_pgs", "pt_pgs")]

# Compute correlation matrix (Pearson by default)
pgs_corr <- cor(pgs_matrix, use = "pairwise.complete.obs")

# Print correlation matrix
print(round(pgs_corr, 3))  # rounded for easier viewing

#########################################################
#P+T graphs

prs.result<-readRDS("C:/Users/haidav/Downloads/Rstudio/prs_result.rds")

df_plot <- data.frame(
  p_thr = prs.result$Threshold,        
  R2    = prs.result$R2,        
  pval  = prs.result$P
)

df_plot <- df_plot %>%
  mutate(
    p = gsub("^Pt_", "", p_thr),  # remove "PT" at the start
    sig = case_when(
      pval < 1e-8 ~ "< 1e-8",
      pval < 1e-6 ~ "< 1e-6",
      pval < 1e-4 ~ "< 1e-4",
      pval < 0.01 ~ "< 0.01",
      pval < 0.05 ~ "< 0.05",
      TRUE        ~ "NS"
    )
  )

sig_colors <- c(
  "< 1e-8" = "#08306b",  
  "< 1e-6" = "#08519c",
  "< 1e-4" = "#2171b5",
  "< 0.01" = "#4292c6",
  "< 0.05" = "#6baed6",
  "NS"     = "#c6dbef"   
)

df_plot

# Bar plot
ggplot(df_plot, aes(x = factor(p), y = R2, fill = sig)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = sig_colors, name = "Model p-value") + # change palette
  theme_minimal(base_size = 14) +
  labs(x = "P-value threshold", y = "McFadden R²", 
       title = "PGS model performance (PRSise)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

################################################
#lassosum2 plot

lasso<-readRDS("lassosum2_results2.rds")

params<-lasso$params
r2<-lasso$r2_values
lasso_dat<-cbind(params,r2)

head(lasso_dat)


ggplot(lasso_dat, aes(x = lambda, y = r2, color = factor(delta))) +
  geom_line() +
  geom_point() +
  scale_x_log10() +  # log scale for lambda
  labs(x = "Lambda", y = expression(R^2), color = "Delta") +
  theme_minimal()
