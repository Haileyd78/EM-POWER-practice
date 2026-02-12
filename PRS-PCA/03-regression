# Regression of ADHD on PRS-PCA
# Reading in data
pcs <- read.table("/scratch-shared/hdavis/PCA.eigenvec", header=F) #PCAs of ancestry
phenotype <- read.csv("/projects/einf2700/NeuroIMAGE_genetics/RICOPILI_merged/ricopili_imp_SELECTIVE_UNTAR_4tim/cobg_dir_genome_wide/traits/Mphen_binary.csv") #phenotype data
prs <- readRDS("ADHD_PRS_PCA.rds")

# Renaming columns
colnames(pcs) <- c("famID", "Subject", paste0("PC",1:10))
colnames(prs)[colnames(prs) == "FID"] <- "famID"
colnames(prs)[colnames(prs) == "IID"] <- "Subject"

# Merging files
pheno<-merge(phenotype, pcs, by="famID")
pheno.prs <- merge(pheno, prs[,c("famID", "ADHD.prs.pc")], by="famID")
pheno.prs <- merge(pheno, prs[,c("famID", "ADHD.prs.pc")], by="famID")

# Model calculation 
null.model <- glm(ADHD_EVER ~  PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = pheno.prs, family = binomial)
model <- glm(ADHD_EVER ~ ADHD.prs.pc + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = pheno.prs, family = binomial)

# Extracting PRS-PCA association results
prs.result <- NULL
prs.r2<- 1 - (as.numeric(logLik(model)) / as.numeric(logLik(null.model)))
prs.coef <- summary(model)$coeff["ADHD.prs.pc",]
prs.beta <- as.numeric(prs.coef[1])
prs.se <- as.numeric(prs.coef[2])
prs.p <- as.numeric(prs.coef[4])


# Storing the results
prs.result <- rbind(prs.result, data.frame(R2=prs.r2, P=prs.p, BETA=prs.beta,SE=prs.se))

print(prs.result) 

