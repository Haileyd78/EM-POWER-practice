# Reading in sumstats
sumstats <- read.table("/gpfs/home4/hdavis/ADHD.Transformed.meta",header=TRUE)

# Adding shared frequency (since it's a case-control GWAS)
sumstats$FRQ_shared <- with(sumstats,
(FRQ_A_38691 * Nca + FRQ_U_186843 * Nco) / (Nca + Nco))

# Creating .ma format 
# Flipping alleles, beta, FRQ to align with the reference panel (optional)
ma <- data.frame(
SNP = sumstats$SNP, 
A1 = sumstats$A2,
A2 = sumstats$A1,
freq = (1-sumstats$FRQ_shared), 
b = (-1*sumstats$BETA), 
se = sumstats$SE,
p = sumstats$P, # adjust column name
N = (4 / (1 / 38691  + 1 / 186843 ))
)

# Saving .ma file
write.table(ma, "/gpfs/home4/hdavis/ADHD_sumstats_flipped.ma", quote=FALSE, row.names=FALSE, col.names=FALSE, sep=" ")



