# Read in summary statistics
dat <- read.table("ADHD.QC4.meta", header=TRUE, stringsAsFactors=FALSE)

# Convert odds ratio to beta
dat$BETA <- log(dat$OR)

# Saving transformed summary statistics 
write.table(dat, "ADHD.Transformed.meta", quote=F, row.names=F, sep="\t")
