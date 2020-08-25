library("PMA")
library(matrixStats)
library(cluster)

normalize_data <- function(data) {
  data <- data - rowMeans(data)
  keep <- (apply(data, 1, sd) != 0)
  data <- (data/apply(data, 1, sd))[keep,]
  return(data)
}

get.elbow <- function(values, is.max) {
  second.derivatives = c()
  for (i in 2:(length(values) - 1)) {
    second.derivative = values[i + 1] + values[i - 1] - 2 * values[i]
    second.derivatives = c(second.derivatives, second.derivative)
  }
  print(second.derivatives)
  if (is.max) {
    return(which.max(second.derivatives) + 1)
  } else {
    return(which.min(second.derivatives) + 1)
  }
}

cancers <- list(c("AML", "BIC", "COAD", "GBM", "KIRC", "LIHC", "LUSC", "OV", "SARC", "SKCM"))
only_primary <- 1


for (i in 1:length(cancers)) {
  cancer <- cancers[i]

  data_path <- paste("./raw_data", cancer, sep="/")
  exp_data <- paste(data_path, "exp", sep="/")
  mirna_data <- paste(data_path, "mirna", sep="/")
  met_data <- paste(data_path, "methy", sep="/")

  out_path <- paste("./data", cancer, sep="/")
  out_multi_cl <-  paste(out_path, paste(cancer, "multiomics", "MCCA.clst", sep="_"), sep="/")

  exp <- read.table(exp_data, header = T)
  mirna <- read.table(mirna_data, header = T)
  met <- read.table(met_data, header = T)

  if (only_primary == 0) {
    keep_tissues = c('01', '02', '03', '05', '06', '07')
  } else {
    if (cancer == 'AML') { keep_tissues = c('03') }
    else if (cancer == 'SKCM') { keep_tissues = c('06') }
    else { keep_tissues = c('01') }
  }

  exp <- exp[,substring(colnames(exp), 14, 15) %in% keep_tissues]
  mirna <- mirna[,substring(colnames(mirna), 14, 15) %in% keep_tissues]
  met <- met[,substring(colnames(met), 14, 15) %in% keep_tissues]

  com_pat = Reduce(intersect, list(row.names(t(exp)), row.names(t(mirna)), row.names(t(met))))
  exp = exp[,com_pat]
  mirna = mirna[,com_pat]
  met = met[,com_pat]

  # Filter features with 0 variance
  exp <- exp[apply(exp, 1, var) > 0,]
  mirna <- mirna[apply(mirna, 1, var) > 0,]
  met <- met[apply(met, 1, var) > 0,]

  # Log Transform sequence data (not methylation data)
  exp = log(1+exp)
  mirna = log(1+mirna)

  # Keep 2000 features with highest variance for exp and met omics
  nfeatures <- 2000
  if (nrow(exp) > nfeatures) {
    exp_feat_vars <- apply(exp, 1, var)
    exp_threshold <- exp_feat_vars[order(exp_feat_vars, decreasing = T)][nfeatures]
    exp <- exp[exp_feat_vars >= exp_threshold,]
  }
  if (nrow(met) > nfeatures) {
    met_feat_vars <- apply(met, 1, var)
    met_threshold <- met_feat_vars[order(met_feat_vars, decreasing = T)][nfeatures]
    met <- met[met_feat_vars >= met_threshold,]
  }

  # Normalize data to have zero mean and standard deviation 1
  exp <- normalize_data(exp)
  mirna <- normalize_data(mirna)
  met <- normalize_data(met)

  # Transpose matrices
  exp <- t(exp)
  mirna <- t(mirna)
  met <- t(met)

  # Run MultiCCA
  om_list <- list(exp, met, mirna)
  names(om_list) <- c("exp", "met", "mirna")
  max_num_clusters <- 15
  #cca_perm <- MultiCCA.permute(om_list)
  cca_ret = MultiCCA(om_list, ncomponents = max_num_clusters)
  # cca_ret$ws = sparse canonical variates for each omic and each ncomponent
  sample_rep = exp %*% cca_ret$ws[[1]] # reduced
  # sample_rep_mi = mirna %*% cca_ret$ws[[3]]
  # sample_rep_met = met %*% cca_ret$ws[[2]]

  explained_vars = sapply(1:max_num_clusters,
                          function(i) sum(unlist(apply(sample_rep[1:i,,drop=F], 2, var))))
  # sum of variances for i canonical variates = v(i) with 1 <= i <= 15 : used to find optimal elbow (to choose nb of clusters)

  # explained_vars_mi = sapply(1:max_num_clusters,
  #                         function(i) sum(unlist(apply(sample_rep_mi[1:i,,drop=F], 2, var))))
  #
  # explained_vars_met = sapply(1:max_num_clusters,
  #                            function(i) sum(unlist(apply(sample_rep_met[1:i,,drop=F], 2, var))))

  dimension <- get.elbow(explained_vars, is.max=F) # Use the second derivative as aproximation to automatically find the elbow
  # dimension_mi <- get.elbow(explained_vars_mi, is.max=F)
  # dimension_met <- get.elbow(explained_vars_met, is.max=F)


  sample_rep = sample_rep[,1:dimension]
  # sample_rep_mi = sample_rep_mi[,1:dimension_mi]
  # sample_rep_met = sample_rep_met[,1:dimension_met]

  sils = c()
  clustering.per.num.clusters = list()

  for (num_clusters in 2:max_num_clusters) {
    cur_clustering = kmeans(sample_rep, num_clusters, iter.max=100, nstart=30)$cluster
    distmat <- as.matrix(dist(sample_rep))
    distmat <- distmat^2 # Use the same matrix as Kmeans > squared euclidean distance
    sil <- silhouette(cur_clustering, dmatrix=distmat)[,3] # Get silhouette WIDTH
    sil <- mean(sil) # ASW
    sils = c(sils, sil)
    clustering.per.num.clusters[[num_clusters - 1]] = cur_clustering
  }

  #cca.clustering = clustering.per.num.clusters[[which.min(sils)]]
  cca.clustering = clustering.per.num.clusters[[which.max(sils)]]

  df = data.frame('Patient' = names(cca.clustering), 'Cluster' = cca.clustering)
  write.table(df, out_multi_cl, col.names = T, row.names=F, quote=F, sep="\t")

}
