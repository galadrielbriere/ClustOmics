library(ggplot2)
library(ggrepel)
library(reshape2)
library(RColorBrewer)


get_surv = function(surv_file) {
  surv = as.numeric(read.table(surv_file, header=F))
  if (surv == 0) surv=surv+10^-6
  surv = -log10(surv)
  return(surv)
}

pval = 0.01
#tot_nb_clinical_labels = 79
clin_labels_per_cancer = c('AML'=4, 'BIC'=11, 'COAD'=10, 'GBM'=5, 'KIRC'=9, 'LIHC'=14, 'LUSC'=10, 'SARC'=4, 'SKCM'=9, 'OV'=3)
cancers <- c("AML", "BIC", "COAD", "GBM", "KIRC", "LIHC", "LUSC", "SARC", "SKCM", "OV")
rel = c('_MULTI_MCCA_NEMO_PINS_SNF_rMKL',
        "_EXP_MIRNA_MET_NEMO_PINS_SNF_rMKL_kmeans",
        "_EXP_MIRNA_MET_NEMO_PINS_SNF_rMKL_kmeans_all")

pins = list(surv = c())
mcca = list(surv = c())
snf = list(surv = c())
mkl = list(surv = c())
nemo = list(surv = c())
multi = list(surv = c())
coca_multi = list(surv = c())

# RESULTS MULTI TO MULTI
for (cancer in cancers){
  pins_surv_file = paste("outOnlyMulti/survival", cancer, paste(cancer, "multiomics", "PINS.pval", sep="_"), sep="/")
  pins_surv = get_surv(pins_surv_file)
  pins$surv = append(pins$surv, pins_surv)

  snf_surv_file = paste("outOnlyMulti/survival", cancer, paste(cancer, "multiomics", "SNF.pval", sep="_"), sep="/")
  snf_surv = get_surv(snf_surv_file)
  snf$surv = append(snf$surv, snf_surv)

  mcca_surv_file = paste("outOnlyMulti/survival", cancer, paste(cancer, "multiomics", "MCCA.pval", sep="_"), sep="/")
  mcca_surv = get_surv(mcca_surv_file)
  mcca$surv = append(mcca$surv, mcca_surv)

  nemo_surv_file = paste("outOnlyMulti/survival", cancer, paste(cancer, "multiomics", "NEMO.pval", sep="_"), sep="/")
  nemo_surv = get_surv(nemo_surv_file)
  nemo$surv = append(nemo$surv, nemo_surv)

  mkl_surv_file = paste("outOnlyMulti/survival", cancer, paste(cancer, "multiomics", "rMKL.pval", sep="_"), sep="/")
  mkl_surv = get_surv(mkl_surv_file)
  mkl$surv = append(mkl$surv, mkl_surv)

  coca_surv_file = paste("outOnlyMulti/survival/coca", paste(cancer, "multiomics", "COCA.pval", sep="_"), sep="/")
  coca_surv = get_surv(coca_surv_file)
  coca_multi$surv = append(coca_multi$surv, coca_surv)
}

for (cancer in cancers){
  rel_name_multi = paste(cancer, rel[1], sep="")

  multi_surv_file = paste("outOnlyMulti/survival", cancer, paste(cancer, rel_name_multi, "ClustOmicsClustering.pval", sep="."), sep="/")
  multi_surv = get_surv(multi_surv_file)
  multi$surv = append(multi$surv, multi_surv)
}

all_surv = c(pins = pins["surv"], snf = snf["surv"], mkl = mkl["surv"], mcca = mcca["surv"], nemo = nemo["surv"], coca = coca_multi["surv"],
             multi = multi["surv"])

# SURV MULTI
df_surv_multi = as.data.frame(cancers)
df_surv_multi$PINS = unlist(pins["surv"])
df_surv_multi$SNF = unlist(snf["surv"])
df_surv_multi$rMKL = unlist(mkl["surv"])
df_surv_multi$MCCA = unlist(mcca["surv"])
df_surv_multi$NEMO = unlist(nemo["surv"])
df_surv_multi_melt = melt(df_surv_multi)
df_surv_multi$`ClustOmics MtoM` = unlist(multi["surv"])
df_surv_multi$`COCA MtoM` = unlist(coca_multi["surv"])
df_surv_multi_melt_all = melt(df_surv_multi)

vals=c("#F8766D","#00BA38","#619CFF", "#E6AB02", "#A9A9A9","black", "#774936")
names(vals) = c("PINS", "SNF", "rMKL", "MCCA", "NEMO", "ClustOmics MtoM", "COCA MtoM")

p <- ggplot(df_surv_multi_melt_all, aes(x=cancers, y=value)) +
  geom_boxplot(data=df_surv_multi_melt, fill="white", show.legend=F) +
  geom_point(data=df_surv_multi_melt_all, aes(shape=variable, fill=variable), alpha=0.8, size=3, position=position_dodge(width=0.5)) +
  scale_fill_manual(name = "Method",labels = names(vals),values = vals) +
  scale_shape_manual(name = "Method",labels = names(vals),values = c(PINS=21, SNF=21, NEMO=21, MCCA=21, rMKL=21, 'ClustOmics MtoM'=22, 'COCA MtoM'=23)) +
  labs(x ="Cancer", y = "-log10(survival p-value)", color = "Method", shape = "Type") + #guides(fill=FALSE) +
  geom_hline(yintercept=-log10(pval), linetype="dashed") + theme_bw()
p

svg(file = paste("outOnlyMulti/plots/Multi_Surv_Overall_Qual_Clust_", pval, ".svg", sep=""), width=10, height=7)
print(p)
dev.off()

# RESULTS SINGLE ONLY MULTI
pins_exp = list(surv = c())
snf_exp = list(surv = c())
mkl_exp = list(surv = c())
nemo_exp = list(surv = c())
kmeans_exp = list(surv = c())
pins_mirna = list(surv = c())
snf_mirna = list(surv = c())
mkl_mirna = list(surv = c())
nemo_mirna = list(surv = c())
kmeans_mirna = list(surv = c())
pins_met = list(surv = c())
snf_met = list(surv = c())
mkl_met = list(surv = c())
nemo_met = list(surv = c())
kmeans_met = list(surv = c())
single = list(surv = c())
coca_single = list(surv = c())

for (cancer in cancers){
  pins_surv_file_exp = paste("outOnlyMulti/survival", cancer, paste(cancer, "expression", "PINS.pval", sep="_"), sep="/")
  pins_surv_exp = get_surv(pins_surv_file_exp)
  pins_exp$surv = append(pins_exp$surv, pins_surv_exp)

  pins_surv_file_met = paste("outOnlyMulti/survival", cancer, paste(cancer, "methylation", "PINS.pval", sep="_"), sep="/")
  pins_surv_met = get_surv(pins_surv_file_met)
  pins_met$surv = append(pins_met$surv, pins_surv_met)
  
  pins_surv_file_mirna = paste("outOnlyMulti/survival", cancer, paste(cancer, "mirna", "PINS.pval", sep="_"), sep="/")
  pins_surv_mirna = get_surv(pins_surv_file_mirna)
  pins_mirna$surv = append(pins_mirna$surv, pins_surv_mirna)

  nemo_surv_file_exp = paste("outOnlyMulti/survival", cancer, paste(cancer, "expression", "NEMO.pval", sep="_"), sep="/")
  nemo_surv_exp = get_surv(nemo_surv_file_exp)
  nemo_exp$surv = append(nemo_exp$surv, nemo_surv_exp)
  
  nemo_surv_file_met = paste("outOnlyMulti/survival", cancer, paste(cancer, "methylation", "NEMO.pval", sep="_"), sep="/")
  nemo_surv_met = get_surv(nemo_surv_file_met)
  nemo_met$surv = append(nemo_met$surv, nemo_surv_met)
  
  nemo_surv_file_mirna = paste("outOnlyMulti/survival", cancer, paste(cancer, "mirna", "NEMO.pval", sep="_"), sep="/")
  nemo_surv_mirna = get_surv(nemo_surv_file_mirna)
  nemo_mirna$surv = append(nemo_mirna$surv, nemo_surv_mirna)
  
  mkl_surv_file_exp = paste("outOnlyMulti/survival", cancer, paste(cancer, "expression", "rMKL.pval", sep="_"), sep="/")
  mkl_surv_exp = get_surv(mkl_surv_file_exp)
  mkl_exp$surv = append(mkl_exp$surv, mkl_surv_exp)
  
  mkl_surv_file_met = paste("outOnlyMulti/survival", cancer, paste(cancer, "methylation", "rMKL.pval", sep="_"), sep="/")
  mkl_surv_met = get_surv(mkl_surv_file_met)
  mkl_met$surv = append(mkl_met$surv, mkl_surv_met)
  
  mkl_surv_file_mirna = paste("outOnlyMulti/survival", cancer, paste(cancer, "mirna", "rMKL.pval", sep="_"), sep="/")
  mkl_surv_mirna = get_surv(mkl_surv_file_mirna)
  mkl_mirna$surv = append(mkl_mirna$surv, mkl_surv_mirna)
  
  snf_surv_file_exp = paste("outOnlyMulti/survival", cancer, paste(cancer, "expression", "SNF.pval", sep="_"), sep="/")
  snf_surv_exp = get_surv(snf_surv_file_exp)
  snf_exp$surv = append(snf_exp$surv, snf_surv_exp)
 
  snf_surv_file_met = paste("outOnlyMulti/survival", cancer, paste(cancer, "methylation", "SNF.pval", sep="_"), sep="/")
  snf_surv_met = get_surv(snf_surv_file_met)
  snf_met$surv = append(snf_met$surv, snf_surv_met)
  
  snf_surv_file_mirna = paste("outOnlyMulti/survival", cancer, paste(cancer, "mirna", "SNF.pval", sep="_"), sep="/")
  snf_surv_mirna = get_surv(snf_surv_file_mirna)
  snf_mirna$surv = append(snf_mirna$surv, snf_surv_mirna)
  
  kmeans_surv_file_exp = paste("outOnlyMulti/survival", cancer, paste(cancer, "expression", "kmeans.pval", sep="_"), sep="/")
  kmeans_surv_exp = get_surv(kmeans_surv_file_exp)
  kmeans_exp$surv = append(kmeans_exp$surv, kmeans_surv_exp)
  
  kmeans_surv_file_met = paste("outOnlyMulti/survival", cancer, paste(cancer, "methylation", "kmeans.pval", sep="_"), sep="/")
  kmeans_surv_met = get_surv(kmeans_surv_file_met)
  kmeans_met$surv = append(kmeans_met$surv, kmeans_surv_met)
  
  kmeans_surv_file_mirna = paste("outOnlyMulti/survival", cancer, paste(cancer, "mirna", "kmeans.pval", sep="_"), sep="/")
  kmeans_surv_mirna = get_surv(kmeans_surv_file_mirna)
  kmeans_mirna$surv = append(kmeans_mirna$surv, kmeans_surv_mirna)
  
  coca_surv_file = paste("outOnlyMulti/survival/coca", paste(cancer, "singleToMulti", "COCA.pval", sep="_"), sep="/")
  coca_surv = get_surv(coca_surv_file)
  coca_single$surv = append(coca_single$surv, coca_surv)
}

for (cancer in cancers){
  rel_name_single = paste(cancer, rel[2], sep="")
  
  single_surv_file = paste("outOnlyMulti/survival", cancer, paste(cancer, rel_name_single, "ClustOmicsClustering.pval", sep="."), sep="/")
  single_surv = get_surv(single_surv_file)
  single$surv = append(single$surv, single_surv)
}


# SINGLE ONLY MULTI
all_surv_single = c(pins_exp = pins_exp["surv"], snf_exp = snf_exp["surv"], mkl_exp = mkl_exp["surv"], nemo_exp = nemo_exp["surv"], kmeans_exp = kmeans_exp["surv"],
             pins_met = pins_met["surv"], snf_met = snf_met["surv"], mkl_met = mkl_met["surv"], nemo_met = nemo_met["surv"],  kmeans_met = kmeans_met["surv"],
             pins_mirna = pins_mirna["surv"], snf_mirna = snf_mirna["surv"], mkl_mirna = mkl_mirna["surv"], nemo_mirna = nemo_mirna["surv"],  kmeans_mirna = kmeans_mirna["surv"],
             single = single["surv"], coca_single = coca_single["surv"])

df_surv_single = as.data.frame(cancers)
df_surv_single$PINS_exp = unlist(pins_exp["surv"])
df_surv_single$SNF_exp = unlist(snf_exp["surv"])
df_surv_single$rMKL_exp = unlist(mkl_exp["surv"])
df_surv_single$NEMO_exp = unlist(nemo_exp["surv"])
df_surv_single$kmeans_exp = unlist(kmeans_exp["surv"])
df_surv_single$PINS_mirna = unlist(pins_mirna["surv"])
df_surv_single$SNF_mirna = unlist(snf_mirna["surv"])
df_surv_single$rMKL_mirna = unlist(mkl_mirna["surv"])
df_surv_single$NEMO_mirna = unlist(nemo_mirna["surv"])
df_surv_single$kmeans_mirna = unlist(kmeans_mirna["surv"])
df_surv_single$PINS_met = unlist(pins_met["surv"])
df_surv_single$SNF_met = unlist(snf_met["surv"])
df_surv_single$rMKL_met = unlist(mkl_met["surv"])
df_surv_single$NEMO_met = unlist(nemo_met["surv"])
df_surv_single$kmeans_met = unlist(kmeans_met["surv"])
df_surv_single_melt = melt(df_surv_single)
df_surv_single$`ClustOmics StoM` = unlist(single["surv"])
df_surv_single$`COCA StoM` = unlist(coca_single["surv"])
df_surv_single_melt_all = melt(df_surv_single)

df_surv_single_melt$omic = gsub(".*_", "", df_surv_single_melt$variable)
df_surv_single_melt$omic = gsub("exp", "Expression", df_surv_single_melt$omic)
df_surv_single_melt$omic = gsub("met", "Methylation", df_surv_single_melt$omic)
df_surv_single_melt$omic = gsub("mirna", "miRNA", df_surv_single_melt$omic)

df_surv_single_melt_all$omic = gsub(".*_", "", df_surv_single_melt_all$variable)
df_surv_single_melt_all$omic = gsub("exp", "Expression", df_surv_single_melt_all$omic)
df_surv_single_melt_all$omic = gsub("met", "Methylation", df_surv_single_melt_all$omic)
df_surv_single_melt_all$omic = gsub("mirna", "miRNA", df_surv_single_melt_all$omic)

color_vals=c("#F8766D","#00BA38","#619CFF", "#774936", "black")
names(color_vals) = c("Expression", "Methylation", "miRNA", "COCA StoM", "ClustOmics StoM")
symb_vals = c(21, 21, 21, 23, 22)
names(symb_vals) = names(color_vals)

p <- ggplot(df_surv_single_melt_all, aes(x=cancers, y=value)) +
  geom_boxplot(data=df_surv_single_melt, fill="white", show.legend=F) +
  geom_point(data=df_surv_single_melt_all, aes(shape=omic, fill=omic), size=3, alpha=0.8, position=position_dodge(width=0.5)) +
  scale_fill_manual(name = "Omic", values = color_vals) +
  scale_shape_manual(name = "Omic", values = symb_vals) +
  labs(x ="Cancer", y = "-log10(survival p-value)", color = "Omic", shape = "Omic", fill="Omic") + #guides(fill=FALSE) +
  geom_hline(yintercept=-log10(pval), linetype="dashed") + theme_bw() 
p

svg(file = paste("outOnlyMulti/plots/Single_Surv_Overall_Qual_Clust_", pval, ".svg", sep=""), width=10, height=7)
print(p)
dev.off()



# RESULTS SINGLE ALL
pins_exp = list(surv = c())
snf_exp = list(surv = c())
mkl_exp = list(surv = c())
nemo_exp = list(surv = c())
kmeans_exp = list(surv = c())
pins_mirna = list(surv = c())
snf_mirna = list(surv = c())
mkl_mirna = list(surv = c())
nemo_mirna = list(surv = c())
kmeans_mirna = list(surv = c())
pins_met = list(surv = c())
snf_met = list(surv = c())
mkl_met = list(surv = c())
nemo_met = list(surv = c())
kmeans_met = list(surv = c())
single = list(surv = c())
coca_single = list(surv = c())

for (cancer in cancers){
  pins_surv_file_exp = paste("outAll/survival", cancer, paste(cancer, "expression", "PINS.pval", sep="_"), sep="/")
  pins_surv_exp = get_surv(pins_surv_file_exp)
  pins_exp$surv = append(pins_exp$surv, pins_surv_exp)
  
  pins_surv_file_met = paste("outAll/survival", cancer, paste(cancer, "methylation", "PINS.pval", sep="_"), sep="/")
  pins_surv_met = get_surv(pins_surv_file_met)
  pins_met$surv = append(pins_met$surv, pins_surv_met)
  
  pins_surv_file_mirna = paste("outAll/survival", cancer, paste(cancer, "mirna", "PINS.pval", sep="_"), sep="/")
  pins_surv_mirna = get_surv(pins_surv_file_mirna)
  pins_mirna$surv = append(pins_mirna$surv, pins_surv_mirna)
  
  nemo_surv_file_exp = paste("outAll/survival", cancer, paste(cancer, "expression", "NEMO.pval", sep="_"), sep="/")
  nemo_surv_exp = get_surv(nemo_surv_file_exp)
  nemo_exp$surv = append(nemo_exp$surv, nemo_surv_exp)
  
  nemo_surv_file_met = paste("outAll/survival", cancer, paste(cancer, "methylation", "NEMO.pval", sep="_"), sep="/")
  nemo_surv_met = get_surv(nemo_surv_file_met)
  nemo_met$surv = append(nemo_met$surv, nemo_surv_met)
  
  nemo_surv_file_mirna = paste("outAll/survival", cancer, paste(cancer, "mirna", "NEMO.pval", sep="_"), sep="/")
  nemo_surv_mirna = get_surv(nemo_surv_file_mirna)
  nemo_mirna$surv = append(nemo_mirna$surv, nemo_surv_mirna)
  
  mkl_surv_file_exp = paste("outAll/survival", cancer, paste(cancer, "expression", "rMKL.pval", sep="_"), sep="/")
  mkl_surv_exp = get_surv(mkl_surv_file_exp)
  mkl_exp$surv = append(mkl_exp$surv, mkl_surv_exp)
  
  mkl_surv_file_met = paste("outAll/survival", cancer, paste(cancer, "methylation", "rMKL.pval", sep="_"), sep="/")
  mkl_surv_met = get_surv(mkl_surv_file_met)
  mkl_met$surv = append(mkl_met$surv, mkl_surv_met)
  
  mkl_surv_file_mirna = paste("outAll/survival", cancer, paste(cancer, "mirna", "rMKL.pval", sep="_"), sep="/")
  mkl_surv_mirna = get_surv(mkl_surv_file_mirna)
  mkl_mirna$surv = append(mkl_mirna$surv, mkl_surv_mirna)
  
  snf_surv_file_exp = paste("outAll/survival", cancer, paste(cancer, "expression", "SNF.pval", sep="_"), sep="/")
  snf_surv_exp = get_surv(snf_surv_file_exp)
  snf_exp$surv = append(snf_exp$surv, snf_surv_exp)
  
  snf_surv_file_met = paste("outAll/survival", cancer, paste(cancer, "methylation", "SNF.pval", sep="_"), sep="/")
  snf_surv_met = get_surv(snf_surv_file_met)
  snf_met$surv = append(snf_met$surv, snf_surv_met)
  
  snf_surv_file_mirna = paste("outAll/survival", cancer, paste(cancer, "mirna", "SNF.pval", sep="_"), sep="/")
  snf_surv_mirna = get_surv(snf_surv_file_mirna)
  snf_mirna$surv = append(snf_mirna$surv, snf_surv_mirna)
  
  kmeans_surv_file_exp = paste("outAll/survival", cancer, paste(cancer, "expression", "kmeans.pval", sep="_"), sep="/")
  kmeans_surv_exp = get_surv(kmeans_surv_file_exp)
  kmeans_exp$surv = append(kmeans_exp$surv, kmeans_surv_exp)
  
  kmeans_surv_file_met = paste("outAll/survival", cancer, paste(cancer, "methylation", "kmeans.pval", sep="_"), sep="/")
  kmeans_surv_met = get_surv(kmeans_surv_file_met)
  kmeans_met$surv = append(kmeans_met$surv, kmeans_surv_met)
  
  kmeans_surv_file_mirna = paste("outAll/survival", cancer, paste(cancer, "mirna", "kmeans.pval", sep="_"), sep="/")
  kmeans_surv_mirna = get_surv(kmeans_surv_file_mirna)
  kmeans_mirna$surv = append(kmeans_mirna$surv, kmeans_surv_mirna)
  
  coca_surv_file = paste("outAll/survival/coca", paste(cancer, "singleToMulti", "COCA.pval", sep="_"), sep="/")
  coca_surv = get_surv(coca_surv_file)
  coca_single$surv = append(coca_single$surv, coca_surv)
}

for (cancer in cancers){
  rel_name_single = paste(cancer, rel[3], sep="")
  
  single_surv_file = paste("outAll/survival", cancer, paste(cancer, rel_name_single, "ClustOmicsClustering.pval", sep="."), sep="/")
  single_surv = get_surv(single_surv_file)
  single$surv = append(single$surv, single_surv)
}


# SINGLE ALL
all_surv_single = c(pins_exp = pins_exp["surv"], snf_exp = snf_exp["surv"], mkl_exp = mkl_exp["surv"], nemo_exp = nemo_exp["surv"], kmeans_exp = kmeans_exp["surv"],
                    pins_met = pins_met["surv"], snf_met = snf_met["surv"], mkl_met = mkl_met["surv"], nemo_met = nemo_met["surv"],  kmeans_met = kmeans_met["surv"],
                    pins_mirna = pins_mirna["surv"], snf_mirna = snf_mirna["surv"], mkl_mirna = mkl_mirna["surv"], nemo_mirna = nemo_mirna["surv"],  kmeans_mirna = kmeans_mirna["surv"],
                    single = single["surv"], coca_single = coca_single["surv"])

df_surv_single = as.data.frame(cancers)
df_surv_single$PINS_exp = unlist(pins_exp["surv"])
df_surv_single$SNF_exp = unlist(snf_exp["surv"])
df_surv_single$rMKL_exp = unlist(mkl_exp["surv"])
df_surv_single$NEMO_exp = unlist(nemo_exp["surv"])
df_surv_single$kmeans_exp = unlist(kmeans_exp["surv"])
df_surv_single$PINS_mirna = unlist(pins_mirna["surv"])
df_surv_single$SNF_mirna = unlist(snf_mirna["surv"])
df_surv_single$rMKL_mirna = unlist(mkl_mirna["surv"])
df_surv_single$NEMO_mirna = unlist(nemo_mirna["surv"])
df_surv_single$kmeans_mirna = unlist(kmeans_mirna["surv"])
df_surv_single$PINS_met = unlist(pins_met["surv"])
df_surv_single$SNF_met = unlist(snf_met["surv"])
df_surv_single$rMKL_met = unlist(mkl_met["surv"])
df_surv_single$NEMO_met = unlist(nemo_met["surv"])
df_surv_single$kmeans_met = unlist(kmeans_met["surv"])
df_surv_single_melt = melt(df_surv_single)
df_surv_single$`ClustOmics StoM` = unlist(single["surv"])
df_surv_single$`COCA StoM` = unlist(coca_single["surv"])
df_surv_single_melt_all = melt(df_surv_single)

df_surv_single_melt$omic = gsub(".*_", "", df_surv_single_melt$variable)
df_surv_single_melt$omic = gsub("exp", "Expression", df_surv_single_melt$omic)
df_surv_single_melt$omic = gsub("met", "Methylation", df_surv_single_melt$omic)
df_surv_single_melt$omic = gsub("mirna", "miRNA", df_surv_single_melt$omic)

df_surv_single_melt_all$omic = gsub(".*_", "", df_surv_single_melt_all$variable)
df_surv_single_melt_all$omic = gsub("exp", "Expression", df_surv_single_melt_all$omic)
df_surv_single_melt_all$omic = gsub("met", "Methylation", df_surv_single_melt_all$omic)
df_surv_single_melt_all$omic = gsub("mirna", "miRNA", df_surv_single_melt_all$omic)

color_vals=c("#F8766D","#00BA38","#619CFF", "#774936", "black")
names(color_vals) = c("Expression", "Methylation", "miRNA", "COCA StoM", "ClustOmics StoM")
symb_vals = c(21, 21, 21, 23, 22)
names(symb_vals) = names(color_vals)

p <- ggplot(df_surv_single_melt_all, aes(x=cancers, y=value)) +
  geom_boxplot(data=df_surv_single_melt, fill="white", show.legend=F) +
  geom_point(data=df_surv_single_melt_all, aes(shape=omic, fill=omic), size=3, alpha=0.8, position=position_dodge(width=0.5)) +
  scale_fill_manual(name = "Omic", values = color_vals) +
  scale_shape_manual(name = "Omic", values = symb_vals) +
  labs(x ="Cancer", y = "-log10(survival p-value)", color = "Omic", shape = "Omic", fill="Omic") + #guides(fill=FALSE) +
  geom_hline(yintercept=-log10(pval), linetype="dashed") + theme_bw() 
p
svg(file = paste("outAll/plots/Single_Surv_Overall_Qual_Clust_", pval, ".svg", sep=""), width=10, height=7)
print(p)
dev.off()

# By method
color_vals=c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#774936", "black")
names(color_vals) = c("SNF", "NEMO", "PINS", "rMKL", "k-means", "COCA StoM", "ClustOmics StoM")
symb_vals = c(21, 21, 21, 21, 21, 23, 22)
names(symb_vals) = names(color_vals)

df_surv_single_melt_all$method = gsub("_.*", "", df_surv_single_melt_all$variable)
df_surv_single_melt_all$method = gsub("kmeans", "k-means", df_surv_single_melt_all$method)

p <- ggplot(df_surv_single_melt_all, aes(x=cancers, y=value)) +
  geom_boxplot(data=df_surv_single_melt, fill="white", show.legend=F) +
  geom_point(data=df_surv_single_melt_all, aes(shape=method, fill=method), size=3, alpha=0.8, position=position_dodge(width=0.5)) +
  scale_fill_manual(name = "Omic", values = color_vals) +
  scale_shape_manual(name = "Omic", values = symb_vals) +
  labs(x ="Cancer", y = "-log10(survival p-value)", color = "Omic", shape = "Omic", fill="Omic") + #guides(fill=FALSE) +
  geom_hline(yintercept=-log10(pval), linetype="dashed") + theme_bw() 
p
