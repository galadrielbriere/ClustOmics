cancers <- c("AML", "BIC", "COAD", "GBM", "KIRC", "LIHC", "LUSC", "OV", "SARC", "SKCM")
datatypes = c("expression", "mirna", "methylation","multiomics")

for (cancer in cancers) {
  for (datatype in datatypes) {
    mklFolder = paste("dev/get_raw_clusterings/rMKL_similarity_matrices", cancer, datatype, sep="/")
    mklFiles = list.files(mklFolder)
    inFile = mklFiles[which(grepl("cls\\-", mklFiles))]
    
    df = read.table(paste(mklFolder, inFile, sep="/"), header=T, sep = ",")
    names(df) = c("Patient", "Cluster")
    
    outFile = paste("data", cancer, paste(cancer, datatype, "rMKL.clst", sep="_"), sep="/")
    write.table(df, outFile, sep="\t", row.names = F, col.names = T, quote = F)
  }
}
