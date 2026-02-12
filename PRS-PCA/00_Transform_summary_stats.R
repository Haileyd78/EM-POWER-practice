# Read in summary statistics
dat <- read.table("/gpfs/home4/hdavis/ADHD_meta_Jan2022_iPSYCH1_iPSYCH2_deCODE_PGC.meta", header=TRUE, stringsAsFactors=FALSE)

# Convert odds ratio to beta
dat$BETA <- log(dat$OR)

# Saving transformed summary statistics 
write.table(dat, "ADHD.Transformed.meta", quote=F, row.names=F, sep="\t")
