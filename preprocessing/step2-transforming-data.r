##Transforming BETA to OR in sumstats
# Read in summary statistics
dat <- read.table("ADHD.QC4.meta", header=TRUE, stringsAsFactors=FALSE)

# Convert odds ratio to beta
dat$BETA <- log(dat$OR)
dat$SE <- dat$SE / dat$OR


# Saving transformed summary statistics 
write.table(dat, "ADHD.Transformed.meta", quote=F, row.names=F, sep="\t")

##Renaming columns in PC file
# Reading PC file
pcs <- read.table("/scratch-shared/hdavis/PCA.eigenvec", header=F) #PCAs of ancestry

# Renaming columns
colnames(pcs) <- c("famID", "Subject", paste0("PC",1:10))

# Saving updated file
write.table(pcs, 
            file = "/gpfs/home4/hdavis/PCA.txt", 
            row.names = FALSE,    
            quote = FALSE,         
            sep = "\t")           

##Cleaning phenotype file
phenos <- read.csv("/projects/einf2700/NeuroIMAGE_genetics/RICOPILI_merged/ricopili_imp_SELECTIVE_UNTAR_4tim/cobg_dir_genome_wide/traits/Mphen_binary.csv")

phenos_clean <- data.frame(
  famID = phenos$famID,
  Subject = phenos$Subject,
  ADHD_EVER = phenos$ADHD_EVER
)

# Save file
write.table(phenos_clean,
            "/gpfs/home4/hdavis/phen_adhd.txt",
            sep = "\t",
            row.names = FALSE,
            quote = FALSE)
