
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ConsensusClusterPlus")

options("repos"=c(CRAN="https://mirrors.nju.edu.cn/CRAN"))
options(BioC_mirror="https://mirrors.nju.edu.cn/bioconductor")


setwd("G:\\cervical_cancer_multiomics\\results\\TMT_protein\\top_variance")

library(ConsensusClusterPlus)
dc = read.csv("top_variance_tumor_z_df_dropna_3000_1674.csv",sep=",",row.names = 1,check.names=FALSE)
dc = as.matrix(dc)
rcc = ConsensusClusterPlus(dc,maxK=5,reps=1000,pItem=0.8,pFeature=1,title="euclidean_km",
                           distance="euclidean",clusterAlg="km",plot='png',seed=1262118322)

cluster <- rcc[[3]]$consensusClass

write.csv(cluster,file="G:\\cervical_cancer_multiomics\\results\\TMT_protein\\top_variance\\euclidean_km\\cluster_3.csv")


