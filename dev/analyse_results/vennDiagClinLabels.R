library(ggpubr)
library(ggvenn)

pval = 0.01
#tot_nb_clinical_labels = 79
clin_labels_per_cancer = c('AML'=4, 'BIC'=11, 'COAD'=10, 'GBM'=5, 'KIRC'=9, 'LIHC'=14, 'LUSC'=10, 'SARC'=4, 'SKCM'=9, 'OV'=3)
cancers <- c("AML", "BIC", "COAD", "GBM", "KIRC", "LIHC", "LUSC", "SARC", "SKCM", "OV")
rel = c('_MULTI_MCCA_NEMO_PINS_SNF_rMKL', "_EXP_MIRNA_MET_NEMO_PINS_SNF_rMKL_kmeans_all")

pins = list()
mcca = list()
snf = list()
mkl = list()
nemo = list()
multi = list()
coca_multi = list()
single = list()
coca_single = list()

pins_exp = list()
snf_exp = list()
mkl_exp = list()
nemo_exp = list()
kmeans_exp = list()
pins_mirna = list()
snf_mirna = list()
mkl_mirna = list()
nemo_mirna = list()
kmeans_mirna = list()
pins_met = list()
snf_met = list()
mkl_met = list()
nemo_met = list()
kmeans_met = list()

get_clin = function(clin_file) {
  clin = read.table(clin_file, header=F, sep="\t")
  clin = clin[which(clin$V2 <= pval),1]
  return(clin)
}

for (cancer in cancers){
  pins_clin_file = paste("outOnlyMulti/clinical", cancer, paste(cancer, "multiomics", "PINS.pval", sep="_"), sep="/")
  pins_clin = get_clin(pins_clin_file)
  pins_clin = as.character(pins_clin)
  if (length(pins_clin) != 0) {
    pins = append(pins, paste(cancer,pins_clin, sep=""))
  }
  
  
  snf_clin_file = paste("outOnlyMulti/clinical", cancer, paste(cancer, "multiomics", "SNF.pval", sep="_"), sep="/")
  snf_clin = get_clin(snf_clin_file)
  snf_clin = as.character(snf_clin)
  if (length(snf_clin) != 0) {
    snf = append(snf, paste(cancer,snf_clin, sep=""))
  }
  
  mcca_clin_file = paste("outOnlyMulti/clinical", cancer, paste(cancer, "multiomics", "MCCA.pval", sep="_"), sep="/")
  mcca_clin = get_clin(mcca_clin_file)
  mcca_clin = as.character(mcca_clin)
  if (length(mcca_clin) != 0) {
    mcca = append(mcca, paste(cancer,mcca_clin, sep=""))
  }
  
  nemo_clin_file = paste("outOnlyMulti/clinical", cancer, paste(cancer, "multiomics", "NEMO.pval", sep="_"), sep="/")
  nemo_clin = get_clin(nemo_clin_file)
  nemo_clin = as.character(nemo_clin)
  if (length(nemo_clin) != 0) {
    nemo = append(nemo, paste(cancer,nemo_clin, sep=""))
  }
  
  mkl_clin_file = paste("outOnlyMulti/clinical", cancer, paste(cancer, "multiomics", "rMKL.pval", sep="_"), sep="/")
  mkl_clin = get_clin(mkl_clin_file)
  mkl_clin = as.character(mkl_clin)
  if (length(mkl_clin) != 0) {
    mkl = append(mkl, paste(cancer,mkl_clin, sep=""))
  }
  
  coca_clin_file = paste("outOnlyMulti/clinical/coca",  paste(cancer, "multiomics", "COCA.pval", sep="_"), sep="/")
  coca_clin = get_clin(coca_clin_file)
  coca_clin = as.character(coca_clin)
  if (length(coca_clin) != 0) {
    coca_multi = append(coca_multi, paste(cancer,coca_clin, sep=""))
  }
  
  coca_clin_file = paste("outAll/clinical/coca",  paste(cancer, "singleToMulti", "COCA.pval", sep="_"), sep="/")
  coca_clin = get_clin(coca_clin_file)
  coca_clin = as.character(coca_clin)
  if (length(coca_clin) != 0) {
    coca_single = append(coca_single, paste(cancer,coca_clin, sep=""))
  }
}

for (cancer in cancers){
  rel_name_multi = paste(cancer, rel[1], sep="")
  rel_name_single = paste(cancer, rel[2], sep="")
  
  multi_clin_file = paste("outOnlyMulti/clinical", cancer, paste(cancer, rel_name_multi, "ClustOmicsClustering.pval", sep="."), sep="/")
  multi_clin = get_clin(multi_clin_file)
  multi_clin = as.character(multi_clin)
  if (length(multi_clin) != 0) {
    multi = append(multi, paste(cancer,multi_clin, sep=""))
  }
  
  single_clin_file = paste("outAll/clinical", cancer, paste(cancer, rel_name_single, "ClustOmicsClustering.pval", sep="."), sep="/")
  single_clin = get_clin(single_clin_file)
  single_clin = as.character(single_clin)
  if (length(single_clin) != 0) {
    single = append(single, paste(cancer,single_clin, sep=""))
  }
}

for (cancer in cancers){
  pins_clin_file_exp = paste("outAll/clinical", cancer, paste(cancer, "expression", "PINS.pval", sep="_"), sep="/")
  pins_clin_exp = get_clin(pins_clin_file_exp)
  pins_clin_exp = as.character(pins_clin_exp)
  if (length(pins_clin_exp) != 0) {
    pins_exp = append(pins_exp, paste(cancer,pins_clin_exp, sep=""))
  }

  pins_clin_file_met = paste("outAll/clinical", cancer, paste(cancer, "methylation", "PINS.pval", sep="_"), sep="/")
  pins_clin_met = get_clin(pins_clin_file_met)
  pins_clin_met = as.character(pins_clin_met)
  if (length(pins_clin_met) != 0) {
    pins_met = append(pins_met, paste(cancer,pins_clin_met, sep=""))
  }
  
  pins_clin_file_mirna = paste("outAll/clinical", cancer, paste(cancer, "mirna", "PINS.pval", sep="_"), sep="/")
  pins_clin_mirna = get_clin(pins_clin_file_mirna)
  pins_clin_mirna = as.character(pins_clin_mirna)
  if (length(pins_clin_mirna) != 0) {
    pins_mirna = append(pins_mirna, paste(cancer,pins_clin_mirna, sep=""))
  }
  
  nemo_clin_file_exp = paste("outAll/clinical", cancer, paste(cancer, "expression", "NEMO.pval", sep="_"), sep="/")
  nemo_clin_exp = get_clin(nemo_clin_file_exp)
  nemo_clin_exp = as.character(nemo_clin_exp)
  if (length(nemo_clin_exp) != 0) {
    nemo_exp = append(nemo_exp, paste(cancer,nemo_clin_exp, sep=""))
  }
  
  nemo_clin_file_met = paste("outAll/clinical", cancer, paste(cancer, "methylation", "NEMO.pval", sep="_"), sep="/")
  nemo_clin_met = get_clin(nemo_clin_file_met)
  nemo_clin_met = as.character(nemo_clin_met)
  if (length(nemo_clin_met) != 0) {
    nemo_met = append(nemo_met, paste(cancer,nemo_clin_met, sep=""))
  }
  
  nemo_clin_file_mirna = paste("outAll/clinical", cancer, paste(cancer, "mirna", "NEMO.pval", sep="_"), sep="/")
  nemo_clin_mirna = get_clin(nemo_clin_file_mirna)
  nemo_clin_mirna = as.character(nemo_clin_mirna)
  if (length(nemo_clin_mirna) != 0) {
    nemo_mirna = append(nemo_mirna, paste(cancer,nemo_clin_mirna, sep=""))
  }
  
  mkl_clin_file_mirna = paste("outAll/clinical", cancer, paste(cancer, "mirna", "rMKL.pval", sep="_"), sep="/")
  mkl_clin_mirna = get_clin(mkl_clin_file_mirna)
  mkl_clin_mirna = as.character(mkl_clin_mirna)
  if (length(mkl_clin_mirna) != 0) {
    mkl_mirna = append(mkl_mirna, paste(cancer,mkl_clin_mirna, sep=""))
  }
  
  mkl_clin_file_exp = paste("outAll/clinical", cancer, paste(cancer, "expression", "rMKL.pval", sep="_"), sep="/")
  mkl_clin_exp = get_clin(mkl_clin_file_exp)
  mkl_clin_exp = as.character(mkl_clin_exp)
  if (length(mkl_clin_exp) != 0) {
    mkl_exp = append(mkl_exp, paste(cancer,mkl_clin_exp, sep=""))
  }
  
  mkl_clin_file_met = paste("outAll/clinical", cancer, paste(cancer, "methylation", "rMKL.pval", sep="_"), sep="/")
  mkl_clin_met = get_clin(mkl_clin_file_met)
  mkl_clin_met = as.character(mkl_clin_met)
  if (length(mkl_clin_met) != 0) {
    mkl_met = append(mkl_met, paste(cancer,mkl_clin_met, sep=""))
  }
  
  snf_clin_file_exp = paste("outAll/clinical", cancer, paste(cancer, "expression", "SNF.pval", sep="_"), sep="/")
  snf_clin_exp = get_clin(snf_clin_file_exp)
  snf_clin_exp = as.character(snf_clin_exp)
  if (length(snf_clin_exp) != 0) {
    snf_exp = append(snf_exp, paste(cancer,snf_clin_exp, sep=""))
  }  
  
  snf_clin_file_mirna = paste("outAll/clinical", cancer, paste(cancer, "mirna", "SNF.pval", sep="_"), sep="/")
  snf_clin_mirna = get_clin(snf_clin_file_mirna)
  snf_clin_mirna = as.character(snf_clin_mirna)
  if (length(snf_clin_mirna) != 0) {
    snf_mirna = append(snf_mirna, paste(cancer,snf_clin_mirna, sep=""))
  }
  
  snf_clin_file_met = paste("outAll/clinical", cancer, paste(cancer, "methylation", "SNF.pval", sep="_"), sep="/")
  snf_clin_met = get_clin(snf_clin_file_met)
  snf_clin_met = as.character(snf_clin_met)
  if (length(snf_clin_met) != 0) {
    snf_met = append(snf_met, paste(cancer,snf_clin_met, sep=""))
  }  
  
  kmeans_clin_file_exp = paste("outAll/clinical", cancer, paste(cancer, "expression", "kmeans.pval", sep="_"), sep="/")
  kmeans_clin_exp = get_clin(kmeans_clin_file_exp)
  kmeans_clin_exp = as.character(kmeans_clin_exp)
  if (length(kmeans_clin_exp) != 0) {
    kmeans_exp = append(kmeans_exp, paste(cancer,kmeans_clin_exp, sep=""))
  }  
  
  kmeans_clin_file_mirna = paste("outAll/clinical", cancer, paste(cancer, "mirna", "kmeans.pval", sep="_"), sep="/")
  kmeans_clin_mirna = get_clin(kmeans_clin_file_mirna)
  kmeans_clin_mirna = as.character(kmeans_clin_mirna)
  if (length(kmeans_clin_mirna) != 0) {
    kmeans_mirna = append(kmeans_mirna, paste(cancer,kmeans_clin_mirna, sep=""))
  }
  
  kmeans_clin_file_met = paste("outAll/clinical", cancer, paste(cancer, "methylation", "kmeans.pval", sep="_"), sep="/")
  kmeans_clin_met = get_clin(kmeans_clin_file_met)
  kmeans_clin_met = as.character(kmeans_clin_met)
  if (length(kmeans_clin_met) != 0) {
    kmeans_met = append(kmeans_met, paste(cancer,kmeans_clin_met, sep=""))
  }  
}


pins = unlist(pins)
mkl = unlist(mkl)
nemo = unlist(nemo)
mcca = unlist(mcca)
snf = unlist(snf)
multi = unlist(multi)
single = unlist(single)
coca_multi = unlist(coca_multi)
coca_single = unlist(coca_single)
pins_exp = unlist(pins_exp)
snf_exp = unlist(snf_exp)
mkl_exp = unlist(mkl_exp)
nemo_exp = unlist(nemo_exp)
kmeans_exp = unlist(kmeans_exp)
pins_mirna = unlist(pins_mirna)
snf_mirna = unlist(snf_mirna)
mkl_mirna = unlist(mkl_mirna)
nemo_mirna = unlist(nemo_mirna)
kmeans_mirna = unlist(kmeans_mirna)
pins_met = unlist(pins_met)
snf_met = unlist(snf_met)
mkl_met = unlist(mkl_met)
nemo_met = unlist(nemo_met)
kmeans_met = unlist(kmeans_met)

inputs_mtom = Reduce(union, list(pins, mkl, nemo, snf, mcca, mkl))
inputs_stom = Reduce(union, list(pins_exp, pins_mirna, pins_met,
                                 mkl_exp, mkl_mirna, mkl_met,
                                 nemo_exp, nemo_mirna, nemo_met,
                                 snf_exp, snf_mirna, snf_met,
                                 kmeans_exp, kmeans_mirna, kmeans_met,
                                 mkl_exp, mkl_mirna, mkl_met))

mtom = list('ClustOmics\nMtoM' = multi, 'COCA\nMtoM' = coca_multi, "Inputs\nMtoM" = inputs_mtom)
stom = list('ClustOmics\nStoM All' = single, 'COCA\nStoM All' = coca_single, "Inputs\nStoM All" = inputs_stom)

a=ggvenn(
  mtom,
  stroke_size = 0.5, set_name_size = 3.5, fill_alpha=0.5,
  fill_color = c("blue", "yellow", "white")
)

b=ggvenn(
  stom,
  stroke_size = 0.5, set_name_size = 3.5, fill_alpha=0.5,
  fill_color = c("purple", "red", "white")
)

x = list('   ClustOmics\nMtoM' = multi, 'COCA\nMtoM' = coca_multi,
      'ClustOmics\nStoM All' = single, 'COCA\nStoM All' = coca_single)
c=ggvenn(
  x,
  stroke_size = 0.5, set_name_size = 3.5, fill_alpha=0.5,
  fill_color = c("blue", "yellow", "purple", "red")
)


figure = ggarrange(
  c,             
  ggarrange(a, b, ncol = 2, labels = c("B", "C")), 
  nrow = 2, 
  labels = "A"       
) 
figure


# pdf(file = "vennDiags.pdf", width=10, height=10)
# figure
# dev.off()