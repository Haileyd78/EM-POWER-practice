# Filter based on INFO
awk 'NR==1 || ($1>=1 && $1<=22 && $8>=0.8)' /gpfs/home4/hdavis/ADHD_meta_Jan2022_iPSYCH1_iPSYCH2_deCODE_PGC.meta > ADHD.QC1.meta

# Filter based on MAF
awk 'NR==1 || ( ($7>=0.01 && $7<=0.99) )' ADHD.QC1.meta > ADHD.QC2.meta

# Remove Duplicate SNPS
awk 'NR==1 || !seen[$2]++' ADHD.QC2.meta > ADHD.QC3.meta

# Remove ambiguous SNPs
awk 'NR==1 || !( ($4=="A" && $5=="T") || ($4=="T" && $5=="A") || ($4=="C" && $5=="G") || ($4=="G" && $5=="C") )' ADHD.QC3.meta > ADHD.QC4.meta


