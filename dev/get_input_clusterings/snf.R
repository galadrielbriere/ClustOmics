library(SNFtool)
library(matrixStats)
library(stringr)

cancers <- c("AML", "BIC", "COAD", "GBM", "KIRC", "LIHC", "LUSC", "OV", "SARC", "SKCM")

# Set to FALSE to include patients that were not measure for all omics
only_multiomic_patients = TRUE
# Only primary tumor tissue
only_primary = TRUE

for (i in 1:length(cancers)) {
  cancer <- cancers[i]

  data_path <- paste("./raw_data", cancer, sep="/")
  exp_data <- paste(data_path, "exp", sep="/")
  mirna_data <- paste(data_path, "mirna", sep="/")
  met_data <- paste(data_path, "methy", sep="/")

  out_path <- paste("./data", cancer, sep="/")
  out_exp_cl <- paste(out_path, paste(cancer, "expression", "SNF.clst", sep="_"), sep="/")
  out_mirna_cl <- paste(out_path, paste(cancer, "mirna", "SNF.clst", sep="_"), sep="/")
  out_met_cl <-  paste(out_path, paste(cancer, "methylation", "SNF.clst", sep="_"), sep="/")
  out_multi_cl <-  paste(out_path, paste(cancer, "multiomics", "SNF.clst", sep="_"), sep="/")

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

  # Data normalisation
  exp <- t(standardNormalization(t(exp)))
  mirna <- t(standardNormalization(t(mirna)))
  met <- t(standardNormalization(t(met)))
  
  if (only_multiomic_patients) {
    print("only_multi")
    # Keep multi-omic patients only
    com_pat = Reduce(intersect, list(colnames(exp), colnames(mirna), colnames(met)))
    exp = exp[,com_pat]
    mirna = mirna[,com_pat]
    met = met[,com_pat]
  }

  # Calculate the pair-wise distance
  sq_dist_exp <- dist2(as.matrix(t(exp)), as.matrix(t(exp)))
  sq_dist_mirna <- dist2(as.matrix(t(mirna)), as.matrix(t(mirna)))
  sq_dist_met <- dist2(as.matrix(t(met)), as.matrix(t(met)))

  # Construct similarity graphs
  sigma <- 0.5
  k_exp <- round(1/10 * ncol(sq_dist_exp))  # Number of nearest neighbors = 1/10 of the nb of samples
  k_met <- round(1/10 * ncol(sq_dist_met))
  k_mirna <- round(1/10 * ncol(sq_dist_mirna))

  aff_exp <- affinityMatrix(sq_dist_exp, k_exp, sigma)
  aff_met <- affinityMatrix(sq_dist_met, k_met, sigma)
  aff_mirna <- affinityMatrix(sq_dist_mirna, k_mirna, sigma)

  # Optimal number of clusters with rotation method
  num_clusters_exp <- estimateNumberOfClustersGivenGraph(aff_exp, 2:15)[[3]]  # between 2 and 15 clusters, rotation cost method
  num_clusters_met <- estimateNumberOfClustersGivenGraph(aff_met, 2:15)[[3]]
  num_clusters_mirna <- estimateNumberOfClustersGivenGraph(aff_mirna, 2:15)[[3]]

  # Spectral clustering
  cl_exp <- spectralClustering(aff_exp, num_clusters_exp)
  names(cl_exp) <- row.names(aff_exp)

  cl_met <- spectralClustering(aff_met, num_clusters_met)
  names(cl_met) <- row.names(aff_met)

  cl_mirna <- spectralClustering(aff_mirna, num_clusters_mirna)
  names(cl_mirna) <- row.names(aff_mirna)

  df_exp = data.frame('Patient' = names(cl_exp), 'Cluster' = cl_exp)
  df_met = data.frame('Patient' = names(cl_met), 'Cluster' = cl_met)
  df_mirna = data.frame('Patient' = names(cl_mirna), 'Cluster' = cl_mirna)

  write.table(df_exp, out_exp_cl, col.names = T, row.names=F, quote=F, sep="\t")
  write.table(df_met, out_met_cl, col.names = T, row.names=F, quote=F, sep="\t")
  write.table(df_mirna, out_mirna_cl, col.names = T, row.names=F, quote=F, sep="\t")

  # Multi-omic
  com_pat = Reduce(intersect, list(row.names(t(exp)), row.names(t(mirna)), row.names(t(met))))
  exp = exp[,com_pat]
  mirna = mirna[,com_pat]
  met = met[,com_pat]

  # Calculate the pair-wise distance
  sq_dist_exp <- dist2(as.matrix(t(exp)), as.matrix(t(exp)))
  sq_dist_mirna <- dist2(as.matrix(t(mirna)), as.matrix(t(mirna)))
  sq_dist_met <- dist2(as.matrix(t(met)), as.matrix(t(met)))

  # Construct similarity graphs
  sigma <- 0.5
  k_exp <- round(1/10 * ncol(sq_dist_exp))  # Number of nearest neighbors = 1/10 of the nb of samples
  k_met <- round(1/10 * ncol(sq_dist_met))
  k_mirna <- round(1/10 * ncol(sq_dist_mirna))

  aff_exp <- affinityMatrix(sq_dist_exp, k_exp, sigma)
  aff_met <- affinityMatrix(sq_dist_met, k_met, sigma)
  aff_mirna <- affinityMatrix(sq_dist_mirna, k_mirna, sigma)

  # Similarity Networks Fusion
  iterations <- 30 # number of iterations to create the fused network

  overall_mat = SNF(list(aff_exp, aff_met, aff_mirna), k_exp, iterations)

  # Optimal number of clusters
  num_clusters_multi <- estimateNumberOfClustersGivenGraph(overall_mat, 2:15)[[3]]

  # Clustering
  cl_multi <- spectralClustering(overall_mat, num_clusters_multi)
  names(cl_multi) <- row.names(overall_mat)

  df_multi = data.frame('Patient' = names(cl_multi), 'Cluster' = cl_multi)

  write.table(df_multi, out_multi_cl, col.names = T, row.names=F, quote=F, sep="\t")
}
