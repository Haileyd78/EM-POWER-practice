# Code for PRS-PCA taken from Coombes et al. 2020
# "A principal component approach to improve association testing with polygenic risk scores"
# https://doi.org/10.1002/gepi.22339

###### FUNCTION TO PERFORM PRS-PCA ##########################################
# INPUTS:
dat <- read.table("/gpfs/home4/hdavis/PRS.all_score", header=TRUE) 
x = "ADHD"
# OUTPUTS:
# list of
#  - data = dataframe with cols (ID , PRS-PCA1 , PRS-PCA2)
#  - r2 = variance explained by each PC of the PRS matrix
#  - loadings = the PC-loadings used to create PRS-PCA1
prs.pc <- function(dat,x){
  xo <- scale(as.matrix(dat[, -c(1, 2)]))  ## scale cols of matrix of only PRSs (remove ID)
  g <- prcomp(xo)   ## perform principal components
  pca.r2 <- g$sdev^2/sum(g$sdev^2)    ## calculate variance explained by each PC
  pc1.loadings <- g$rotation[,1];     ## loadings for PC1
  pc2.loadings <- g$rotation[,2]  	## loadings for PC2
  ## flip direction of PCs to keep direction of association
  ## (sign of loadings for PC1 is arbitrary so we want to keep same direction)
  if (mean(pc1.loadings>0)==0){ 	
      pc1.loadings <- pc1.loadings*(-1)
      pc2.loadings <- pc2.loadings*(-1)
  }
  ## calculate PRS-PCA (outputs PC1 and PC2 even though PC1 sufficient)
  pc1 <- xo %*% pc1.loadings
  pc2 <- xo %*% pc2.loadings
  dat[,paste0(x,".prs.pc")] <- scale(pc1)   ## rescales PRS-PCA1
  dat[,paste0(x,".prs.pc2")] <- scale(pc2)  ## rescales PRS-PCA2
  return(list(data=dat,r2=pca.r2,loadings=pc1.loadings))
}

result <- prs.pc(dat, x)

saveRDS(result$data, "ADHD_PRS_PCA.rds")
