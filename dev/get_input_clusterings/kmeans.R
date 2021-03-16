library(cluster)

cancers <- c("AML", "BIC", "COAD", "GBM", "KIRC", "LIHC", "LUSC", "OV", "SARC", "SKCM")

# Set to FALSE to include patients that were not measure for all omics
only_multiomic_patients = FALSE
# Only primary tumor tissue
only_primary = TRUE

for (cancer in cancers) {
  data_path <- paste("./raw_data", cancer, sep="/")
  exp_data <- paste(data_path, "exp", sep="/")
  mirna_data <- paste(data_path, "mirna", sep="/")
  met_data <- paste(data_path, "methy", sep="/")
  
  out_path <- paste("./data", cancer, sep="/")

  exp <- read.table(exp_data, header = T)
  mirna <- read.table(mirna_data, header = T)
  met <- read.table(met_data, header = T)
  
  if (!only_primary) {
    keep_tissues = c('01', '02', '03', '05', '06', '07')
  } else {
    if (cancer == 'AML') { keep_tissues = c('03') }
    else if (cancer == 'SKCM') { keep_tissues = c('06') }
    else { keep_tissues = c('01') }
  }
  
  exp <- exp[,substring(colnames(exp), 14, 15) %in% keep_tissues]
  mirna <- mirna[,substring(colnames(mirna), 14, 15) %in% keep_tissues]
  met <- met[,substring(colnames(met), 14, 15) %in% keep_tissues]
  
  # Filter features with 0 variance
  exp <- exp[apply(exp, 1, var) > 0,]
  mirna <- mirna[apply(mirna, 1, var) > 0,]
  met <- met[apply(met, 1, var) > 0,]
  
  # Log Transform sequence data (not methylation data)
  exp = log(1+exp)
  mirna = log(1+mirna)
  
  if (only_multiomic_patients) {
    print("only_multi")
    # Keep multi-omic patients only
    com_pat = Reduce(intersect, list(colnames(exp), colnames(mirna), colnames(met)))
    exp = exp[,com_pat]
    mirna = mirna[,com_pat]
    met = met[,com_pat]
  }
  
  for (dataset in c("exp", "mirna", "met")) {
    if (dataset == "exp") {
      dataset = exp
      out <- paste(out_path, paste(cancer, "expression", "kmeans.clst", sep="_"), sep="/")
    } else if (dataset == "mirna") {
      dataset = mirna
      out <- paste(out_path, paste(cancer, "mirna", "kmeans.clst", sep="_"), sep="/")
    } else {
     dataset = met 
     out <-  paste(out_path, paste(cancer, "methylation", "kmeans.clst", sep="_"), sep="/")
    }
  
    # Choose K using silhouette
    sil <- rep(0, 20)
    clusterings = list()
    #repeat k-means for k 1:20 and extract silhouette:
    for(i in 2:20){
      k1to20 <- kmeans(t(dataset), centers = i)#, nstart = 25, iter.max = 20)
      clusterings[[i]] = k1to20
      ss <- silhouette(k1to20$cluster, dist(t(dataset)))
      sil[i] <- mean(ss[, 3])
    }
    
    kopt = which.max(sil)
    
    # Plot the  average silhouette width
    plot(1:20, sil, type = "b", pch = 19, xlab = "Number of clusters k", ylab="Average silhouette width")
    abline(v = kopt, lty = 2)
    
    clustopt = clusterings[[kopt]]  
    
    df = data.frame('Patient' = names(clustopt$cluster), 'Cluster' = clustopt$cluster)
    write.table(df, out, col.names = T, row.names=F, quote=F, sep="\t")
    
  }
}
