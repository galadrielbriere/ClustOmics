library(SNFtool)
library(matrixStats)
library(stringr)
source('NEMO.R')

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
  out_exp_cl <- paste(out_path, paste(cancer, "expression", "NEMO.clst", sep="_"), sep="/")
  out_mirna_cl <- paste(out_path, paste(cancer, "mirna", "NEMO.clst", sep="_"), sep="/")
  out_met_cl <-  paste(out_path, paste(cancer, "methylation", "NEMO.clst", sep="_"), sep="/")
  out_multi_cl <-  paste(out_path, paste(cancer, "multiomics", "NEMO.clst", sep="_"), sep="/")

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
  
  # Run NEMO
  cl_exp = nemo.clustering(list(exp))
  cl_mirna = nemo.clustering(list(mirna))
  cl_met = nemo.clustering(list(met))
  cl_multi = nemo.clustering(list(exp, mirna, met))

  df_exp =  data.frame('Patient' = names(cl_exp), 'Cluster' = cl_exp)
  df_mirna =  data.frame('Patient' = names(cl_mirna), 'Cluster' = cl_mirna)
  df_met =  data.frame('Patient' = names(cl_met), 'Cluster' = cl_met)
  df_multi =  data.frame('Patient' = names(cl_multi), 'Cluster' = cl_multi)


  write.table(df_exp, out_exp_cl, col.names = T, row.names = F, quote=F, sep="\t")
  write.table(df_met, out_met_cl, col.names = T, row.names = F, quote=F, sep="\t")
  write.table(df_mirna, out_mirna_cl, col.names = T, row.names = F, quote=F, sep="\t")
  write.table(df_multi, out_multi_cl, col.names = T, row.names = F, quote=F, sep="\t")
}
