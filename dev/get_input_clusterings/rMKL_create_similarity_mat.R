library(matrixStats)
library(stringr)

cancers <- c("AML", "BIC", "COAD", "GBM", "KIRC", "LIHC", "LUSC", "OV", "SARC", "SKCM")

# Set to FALSE to include patients that were not measure for all omics
only_multiomic_patients = TRUE
# Only primary tumor tissue
only_primary = TRUE

normalize_data <- function(data) {
  data <- data - rowMeans(data)
  keep <- (apply(data, 1, sd) != 0)
  data <- (data/apply(data, 1, sd))[keep,]
  return(data)
}

make_similarity_matrices <- function(data, out_dir, datatype, cancer) {
  npatients = ncol(data)
  gammas = 10 ** seq(-6, 6, by=3)

  for (k in 1:length(gammas)) {
    gamma = gammas[[k]] / (2*nrow(data)**2)
    output.mat = matrix(0, ncol=npatients, nrow=npatients)
    for (i in 1:npatients) {
      for (j in 1:npatients) {
        output.mat[i, j] = exp(-norm(as.matrix(data[,i] - data[,j]), type = 'F')**2 * gamma)
      }
    }
    D = apply(output.mat, 2, sum) / npatients
    E = sum(D) / npatients
    J = matrix(1, nrow=npatients, ncol=1) %*% D
    ret = output.mat - J - t(J) + E * matrix(1, ncol=npatients, nrow=npatients)
    ret = diag(1/sqrt(diag(ret))) %*% ret %*% diag(1/sqrt(diag(ret)))
    write.table(ret, paste(out_dir, paste(cancer, datatype, gammas[[k]], sep="_"), sep='/'), col.names = F, row.names=F, quote=F, sep="\t")
  }
}

for (cancer in cancers) {
  
  data_path <- paste("./raw_data", cancer, sep="/")
  exp_data <- paste(data_path, "exp", sep="/")
  mirna_data <- paste(data_path, "mirna", sep="/")
  met_data <- paste(data_path, "methy", sep="/")

  out_path <- paste("./dev/get_raw_clusterings/rMKL_similarity_matrices", cancer, sep="/")
  out_exp <- paste(out_path, "expression", sep="/")
  out_mirna <- paste(out_path, "mirna", sep="/")
  out_met <-  paste(out_path, "methylation", sep="/")
  out_multi <-  paste(out_path, "multiomics", sep="/")
  dir.create(out_path)
  dir.create(out_exp)
  dir.create(out_met)
  dir.create(out_mirna)
  dir.create(out_multi)

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

  # Normalize data
  exp <- normalize_data(exp)
  mirna <- normalize_data(mirna)
  met <- normalize_data(met)
  
  if (only_multiomic_patients) {
    print("only_multi")
    # Keep multi-omic patients only
    com_pat = Reduce(intersect, list(colnames(exp), colnames(mirna), colnames(met)))
    exp = exp[,com_pat]
    mirna = mirna[,com_pat]
    met = met[,com_pat]
  }

  # Single-omic similarity matrices
  make_similarity_matrices(exp, out_exp, 'expression', cancer)
  make_similarity_matrices(mirna, out_mirna, 'miRNA', cancer)
  make_similarity_matrices(met, out_met, 'methylation', cancer)

  # Write ids files
  write.table(names(exp), paste(out_path, "ids_expression.txt", sep="/"), row.names=F, col.names=F, quote=F)
  write.table(names(mirna), paste(out_path, "ids_mirna.txt", sep="/"), row.names=F, col.names=F, quote=F)
  write.table(names(met), paste(out_path, "ids_methylation.txt", sep="/"), row.names=F, col.names=F, quote=F)

  # Multi-omic similarity matrice
  com_pat = Reduce(intersect, list(names(exp), names(mirna), names(met)))
  exp = exp[,com_pat]
  met = met[,com_pat]
  mirna = mirna[,com_pat]

  make_similarity_matrices(exp, out_multi, 'expression', cancer)
  make_similarity_matrices(mirna, out_multi, 'miRNA', cancer)
  make_similarity_matrices(met, out_multi, 'methylation', cancer)

  write.table(com_pat, paste(out_path, "ids_multiomics.txt", sep="/"), row.names=F, col.names=F, quote=F)
}
