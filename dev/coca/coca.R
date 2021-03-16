library(coca)
library(reshape2)

cancers <- c("AML", "BIC", "COAD", "GBM", "KIRC", "LIHC", "LUSC", "OV", "SARC", "SKCM")

for (cancer in cancers) {
  dataOM_folder = paste('dataOnlyMulti/', cancer, sep="")
  filesOM = list.files(dataOM_folder)
  
  ########################  MULTI TO MULTI
  files_multi = filesOM[which(grepl("multiomics", filesOM))]
  out_multi = paste("dataOnlyMulti/coca/", cancer, "_multiomics_COCA.clst", sep="")

  # Build Matrix Of Clusters
  data <- data.frame()
  for (file in paste(dataOM_folder, files_multi, sep="/")) {
    clust = read.table(file, header=T, sep="\t")
    clust$Cluster = paste(gsub("\\.clst", "", gsub(".*_", "", file)), clust$Cluster, sep="_")
    data = rbind(data, clust)
  }

  moc = dcast(data, Patient ~ Cluster)
  rownames(moc) = moc$Patient
  moc = moc[,-1]
  moc = ifelse(is.na(moc), 0, 1)

  coca <- coca::coca(moc)
  coca_clust = coca$clusterLabels
  names(coca_clust) = rownames(moc)
  coca_clust = as.data.frame(coca_clust)
  coca_clust$Patient = rownames(coca_clust)
  coca_clust = coca_clust[,c(2,1)]
  names(coca_clust) = c("Patient", "Cluster")

  write.table(coca_clust, out_multi, col.names = T, row.names = F, quote=F, sep="\t")
  
  ######################## SINGLE TO MULTI OnlyMulti
  files_singleOM = filesOM[which(grepl("expression", filesOM)|grepl("mirna", filesOM)|grepl("methylation", filesOM))]
  out_singleOM = paste("dataOnlyMulti/coca/", cancer, "_singleToMulti_COCA.clst", sep="")
  
  # Build Matrix Of Clusters 
  data <- data.frame()
  for (file in paste(dataOM_folder, files_singleOM, sep="/")) {
    clust = read.table(file, header=T, sep="\t") 
    clust$Cluster = paste(gsub("\\.clst", "", sub(".*?_", "", file)), clust$Cluster, sep="_")
    data = rbind(data, clust)
  }
  
  moc = dcast(data, Patient ~ Cluster)
  rownames(moc) = moc$Patient
  moc = moc[,-1]
  moc = ifelse(is.na(moc), 0, 1)
  
  coca <- coca::coca(moc)
  coca_clust = coca$clusterLabels
  names(coca_clust) = rownames(moc)
  coca_clust = as.data.frame(coca_clust)
  coca_clust$Patient = rownames(coca_clust)
  coca_clust = coca_clust[,c(2,1)]
  names(coca_clust) = c("Patient", "Cluster")
  
  write.table(coca_clust, out_singleOM, col.names = T, row.names = F, quote=F, sep="\t")
  
  ########################  SINGLE TO MULTI All
  dataAll_folder = paste('dataAll/', cancer, sep="")
  filesAll = list.files(dataAll_folder)
  files_singleAll = filesAll[which(grepl("expression", filesAll)|grepl("mirna", filesAll)|grepl("methylation", filesAll))]
  out_singleAll = paste("dataAll/coca/", cancer, "_singleToMulti_COCA.clst", sep="")
  
  # Build Matrix Of Clusters 
  data <- data.frame()
  for (file in paste(dataAll_folder, files_singleAll, sep="/")) {
    clust = read.table(file, header=T, sep="\t") 
    clust$Cluster = paste(gsub("\\.clst", "", sub(".*?_", "", file)), clust$Cluster, sep="_")
    data = rbind(data, clust)
  }
  
  moc = dcast(data, Patient ~ Cluster)
  rownames(moc) = moc$Patient
  moc = moc[,-1]
  moc = ifelse(is.na(moc), 0, 1)
  
  coca <- coca::coca(moc)
  coca_clust = coca$clusterLabels
  names(coca_clust) = rownames(moc)
  coca_clust = as.data.frame(coca_clust)
  coca_clust$Patient = rownames(coca_clust)
  coca_clust = coca_clust[,c(2,1)]
  names(coca_clust) = c("Patient", "Cluster")
  
  write.table(coca_clust, out_singleAll, col.names = T, row.names = F, quote=F, sep="\t")
}