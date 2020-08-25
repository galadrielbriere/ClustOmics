library(ggplot2)
library(ggrepel)
library(reshape2)

pval = 0.01
clin_labels_per_cancer = c('AML'=4, 'BIC'=11, 'COAD'=10, 'GBM'=5, 'KIRC'=9, 'LIHC'=14, 'LUSC'=10, 'SARC'=4, 'SKCM'=9, 'OV'=3)
cancers <- c("AML", "BIC", "COAD", "GBM", "KIRC", "LIHC", "LUSC", "SARC", "SKCM", "OV")
rel = c('_MULTI_MCCA_NEMO_PINS_SNF_rMKL', "_EXP_MIRNA_MET_NEMO_PINS_SNF_rMKL")

pins = list(surv = c(), clin = c())
mcca = list(surv = c(), clin = c())
snf = list(surv = c(), clin = c())
mkl = list(surv = c(), clin = c())
nemo = list(surv = c(), clin = c())
multi = list(surv = c(), clin = c())
single = list(surv = c(), clin = c())

pins_exp = list(surv = c(), clin = c())
snf_exp = list(surv = c(), clin = c())
mkl_exp = list(surv = c(), clin = c())
nemo_exp = list(surv = c(), clin = c())
pins_mirna = list(surv = c(), clin = c())
snf_mirna = list(surv = c(), clin = c())
mkl_mirna = list(surv = c(), clin = c())
nemo_mirna = list(surv = c(), clin = c())
pins_met = list(surv = c(), clin = c())
snf_met = list(surv = c(), clin = c())
mkl_met = list(surv = c(), clin = c())
nemo_met = list(surv = c(), clin = c())


get_surv = function(surv_file) {
  surv = as.numeric(read.table(surv_file, header=F))
  surv = -log10(surv)
  return(surv)
}

get_clin = function(clin_file) {
  clin = read.table(clin_file, header=F, sep="\t")
  clin = length(which((clin$V2 <= pval)))
  return(clin)
}

# RESULTS MULTI
for (cancer in cancers){
  pins_surv_file = paste("survival", cancer, paste(cancer, "multiomics", "PINS.pval", sep="_"), sep="/")
  pins_surv = get_surv(pins_surv_file)
  pins$surv = append(pins$surv, pins_surv)
  
  pins_clin_file = paste("clinical", cancer, paste(cancer, "multiomics", "PINS.pval", sep="_"), sep="/")
  pins_clin = get_clin(pins_clin_file)
  pins$clin = append(pins$clin, pins_clin)
  
  snf_surv_file = paste("survival", cancer, paste(cancer, "multiomics", "SNF.pval", sep="_"), sep="/")
  snf_surv = get_surv(snf_surv_file)
  snf$surv = append(snf$surv, snf_surv)
  
  snf_clin_file = paste("clinical", cancer, paste(cancer, "multiomics", "SNF.pval", sep="_"), sep="/")
  snf_clin = get_clin(snf_clin_file)
  snf$clin = append(snf$clin, snf_clin)
  
  mcca_surv_file = paste("survival", cancer, paste(cancer, "multiomics", "MCCA.pval", sep="_"), sep="/")
  mcca_surv = get_surv(mcca_surv_file)
  mcca$surv = append(mcca$surv, mcca_surv)
  
  mcca_clin_file = paste("clinical", cancer, paste(cancer, "multiomics", "MCCA.pval", sep="_"), sep="/")
  mcca_clin = get_clin(mcca_clin_file)
  mcca$clin = append(mcca$clin, mcca_clin)
  
  nemo_surv_file = paste("survival", cancer, paste(cancer, "multiomics", "NEMO.pval", sep="_"), sep="/")
  nemo_surv = get_surv(nemo_surv_file)
  nemo$surv = append(nemo$surv, nemo_surv)
  
  nemo_clin_file = paste("clinical", cancer, paste(cancer, "multiomics", "NEMO.pval", sep="_"), sep="/")
  nemo_clin = get_clin(nemo_clin_file)
  nemo$clin = append(nemo$clin, nemo_clin)
  
  mkl_surv_file = paste("survival", cancer, paste(cancer, "multiomics", "rMKL.pval", sep="_"), sep="/")
  mkl_surv = get_surv(mkl_surv_file)
  mkl$surv = append(mkl$surv, mkl_surv)
  
  mkl_clin_file = paste("clinical", cancer, paste(cancer, "multiomics", "rMKL.pval", sep="_"), sep="/")
  mkl_clin = get_clin(mkl_clin_file)
  mkl$clin = append(mkl$clin, mkl_clin)
}

# RESULTS SINGLE
for (cancer in cancers){
  pins_surv_file_exp = paste("survival", cancer, paste(cancer, "expression", "PINS.pval", sep="_"), sep="/")
  pins_surv_exp = get_surv(pins_surv_file_exp)
  pins_exp$surv = append(pins_exp$surv, pins_surv_exp)
  
  pins_clin_file_exp = paste("clinical", cancer, paste(cancer, "expression", "PINS.pval", sep="_"), sep="/")
  pins_clin_exp = get_clin(pins_clin_file_exp)
  pins_exp$clin = append(pins_exp$clin, pins_clin_exp)
  
  pins_surv_file_met = paste("survival", cancer, paste(cancer, "methylation", "PINS.pval", sep="_"), sep="/")
  pins_surv_met = get_surv(pins_surv_file_met)
  pins_met$surv = append(pins_met$surv, pins_surv_met)
  
  pins_clin_file_met = paste("clinical", cancer, paste(cancer, "methylation", "PINS.pval", sep="_"), sep="/")
  pins_clin_met = get_clin(pins_clin_file_met)
  pins_met$clin = append(pins_met$clin, pins_clin_met)
  
  pins_surv_file_mirna = paste("survival", cancer, paste(cancer, "mirna", "PINS.pval", sep="_"), sep="/")
  pins_surv_mirna = get_surv(pins_surv_file_mirna)
  pins_mirna$surv = append(pins_mirna$surv, pins_surv_mirna)
  
  pins_clin_file_mirna = paste("clinical", cancer, paste(cancer, "mirna", "PINS.pval", sep="_"), sep="/")
  pins_clin_mirna = get_clin(pins_clin_file_mirna)
  pins_mirna$clin = append(pins_mirna$clin, pins_clin_mirna)
  
  nemo_surv_file_exp = paste("survival", cancer, paste(cancer, "expression", "NEMO.pval", sep="_"), sep="/")
  nemo_surv_exp = get_surv(nemo_surv_file_exp)
  nemo_exp$surv = append(nemo_exp$surv, nemo_surv_exp)
  
  nemo_clin_file_exp = paste("clinical", cancer, paste(cancer, "expression", "NEMO.pval", sep="_"), sep="/")
  nemo_clin_exp = get_clin(nemo_clin_file_exp)
  nemo_exp$clin = append(nemo_exp$clin, nemo_clin_exp)
  
  nemo_surv_file_met = paste("survival", cancer, paste(cancer, "methylation", "NEMO.pval", sep="_"), sep="/")
  nemo_surv_met = get_surv(nemo_surv_file_met)
  nemo_met$surv = append(nemo_met$surv, nemo_surv_met)
  
  nemo_clin_file_met = paste("clinical", cancer, paste(cancer, "methylation", "NEMO.pval", sep="_"), sep="/")
  nemo_clin_met = get_clin(nemo_clin_file_met)
  nemo_met$clin = append(nemo_met$clin, nemo_clin_met)
  
  nemo_surv_file_mirna = paste("survival", cancer, paste(cancer, "mirna", "NEMO.pval", sep="_"), sep="/")
  nemo_surv_mirna = get_surv(nemo_surv_file_mirna)
  nemo_mirna$surv = append(nemo_mirna$surv, nemo_surv_mirna)
  
  nemo_clin_file_mirna = paste("clinical", cancer, paste(cancer, "mirna", "NEMO.pval", sep="_"), sep="/")
  nemo_clin_mirna = get_clin(nemo_clin_file_mirna)
  nemo_mirna$clin = append(nemo_mirna$clin, nemo_clin_mirna)
  
  mkl_surv_file_exp = paste("survival", cancer, paste(cancer, "expression", "rMKL.pval", sep="_"), sep="/")
  mkl_surv_exp = get_surv(mkl_surv_file_exp)
  mkl_exp$surv = append(mkl_exp$surv, mkl_surv_exp)
  
  mkl_clin_file_exp = paste("clinical", cancer, paste(cancer, "expression", "rMKL.pval", sep="_"), sep="/")
  mkl_clin_exp = get_clin(mkl_clin_file_exp)
  mkl_exp$clin = append(mkl_exp$clin, mkl_clin_exp)
  
  mkl_surv_file_met = paste("survival", cancer, paste(cancer, "methylation", "rMKL.pval", sep="_"), sep="/")
  mkl_surv_met = get_surv(mkl_surv_file_met)
  mkl_met$surv = append(mkl_met$surv, mkl_surv_met)
  
  mkl_clin_file_met = paste("clinical", cancer, paste(cancer, "methylation", "rMKL.pval", sep="_"), sep="/")
  mkl_clin_met = get_clin(mkl_clin_file_met)
  mkl_met$clin = append(mkl_met$clin, mkl_clin_met)
  
  mkl_surv_file_mirna = paste("survival", cancer, paste(cancer, "mirna", "rMKL.pval", sep="_"), sep="/")
  mkl_surv_mirna = get_surv(mkl_surv_file_mirna)
  mkl_mirna$surv = append(mkl_mirna$surv, mkl_surv_mirna)
  
  mkl_clin_file_mirna = paste("clinical", cancer, paste(cancer, "mirna", "rMKL.pval", sep="_"), sep="/")
  mkl_clin_mirna = get_clin(mkl_clin_file_mirna)
  mkl_mirna$clin = append(mkl_mirna$clin, mkl_clin_mirna)
  
  snf_surv_file_exp = paste("survival", cancer, paste(cancer, "expression", "SNF.pval", sep="_"), sep="/")
  snf_surv_exp = get_surv(snf_surv_file_exp)
  snf_exp$surv = append(snf_exp$surv, snf_surv_exp)
  
  snf_clin_file_exp = paste("clinical", cancer, paste(cancer, "expression", "SNF.pval", sep="_"), sep="/")
  snf_clin_exp = get_clin(snf_clin_file_exp)
  snf_exp$clin = append(snf_exp$clin, snf_clin_exp)
  
  snf_surv_file_met = paste("survival", cancer, paste(cancer, "methylation", "SNF.pval", sep="_"), sep="/")
  snf_surv_met = get_surv(snf_surv_file_met)
  snf_met$surv = append(snf_met$surv, snf_surv_met)
  
  snf_clin_file_met = paste("clinical", cancer, paste(cancer, "methylation", "SNF.pval", sep="_"), sep="/")
  snf_clin_met = get_clin(snf_clin_file_met)
  snf_met$clin = append(snf_met$clin, snf_clin_met)
  
  snf_surv_file_mirna = paste("survival", cancer, paste(cancer, "mirna", "SNF.pval", sep="_"), sep="/")
  snf_surv_mirna = get_surv(snf_surv_file_mirna)
  snf_mirna$surv = append(snf_mirna$surv, snf_surv_mirna)
  
  snf_clin_file_mirna = paste("clinical", cancer, paste(cancer, "mirna", "SNF.pval", sep="_"), sep="/")
  snf_clin_mirna = get_clin(snf_clin_file_mirna)
  snf_mirna$clin = append(snf_mirna$clin, snf_clin_mirna)
}

for (cancer in cancers){
  rel_name_multi = paste(cancer, rel[1], sep="")
  rel_name_single = paste(cancer, rel[2], sep="")
  
  multi_surv_file = paste("survival", cancer, paste(cancer, rel_name_multi, "ClustOmicsClustering.pval", sep="."), sep="/")
  multi_surv = get_surv(multi_surv_file)
  multi$surv = append(multi$surv, multi_surv)
  
  single_surv_file = paste("survival", cancer, paste(cancer, rel_name_single, "ClustOmicsClustering.pval", sep="."), sep="/")
  single_surv = get_surv(single_surv_file)
  single$surv = append(single$surv, single_surv)
  
  multi_clin_file = paste("clinical", cancer, paste(cancer, rel_name_multi, "ClustOmicsClustering.pval", sep="."), sep="/")
  multi_clin = get_clin(multi_clin_file)
  multi$clin = append(multi$clin, multi_clin)
  
  single_clin_file = paste("clinical", cancer, paste(cancer, rel_name_single, "ClustOmicsClustering.pval", sep="."), sep="/")
  single_clin = get_clin(single_clin_file)
  single$clin = append(single$clin, single_clin)
}

all_surv = c(pins = pins["surv"], snf = snf["surv"], mkl = mkl["surv"], mcca = mcca["surv"], nemo = nemo["surv"],
             multi = multi["surv"])
all_enriched = c(pins = pins["clin"], snf = snf["clin"], mkl = mkl["clin"], mcca = mcca["clin"], nemo = nemo["clin"],
                 multi = multi["clin"])

# SURV MULTI
df_surv_multi = as.data.frame(cancers)
df_surv_multi$PINS = unlist(pins["surv"])
df_surv_multi$SNF = unlist(snf["surv"])
df_surv_multi$rMKL = unlist(mkl["surv"])
df_surv_multi$MCCA = unlist(mcca["surv"])
df_surv_multi$NEMO = unlist(nemo["surv"])
df_surv_multi$ClustOmicsMtoM = unlist(multi["surv"])
df_surv_multi = df_surv_multi[,c(1, 7, 2:6)]
df_surv_multi = melt(df_surv_multi)

p <- ggplot(df_surv_multi, aes(x=cancers, y=value, fill=variable)) + geom_boxplot(fill="white") + 
  geom_dotplot(binaxis='y', stackdir='center', dotsize = 0.5, position=position_dodge(width=0.2)) + labs(x ="Cancer", y = "-log10(survival p-value)", fill = "Method") + 
  geom_hline(yintercept=-log10(pval), linetype="dashed") + scale_fill_manual(values=c(PINS = "#F8766D", SNF="#00BA38", rMKL="#619CFF", MCCA="#B79F00", NEMO="#00BFC4", ClustOmicsMtoM="#C77CFF"))

svg(file = paste("plots/Multi_Surv_Overall_Qual_Clust_", pval, ".svg", sep=""), width=10, height=7)
print(p)
dev.off()

# CLIN MULTI
df_clin_multi = as.data.frame(cancers)
df_clin_multi$PINS = unlist(pins["clin"])
df_clin_multi$SNF = unlist(snf["clin"])
df_clin_multi$rMKL = unlist(mkl["clin"])
df_clin_multi$MCCA = unlist(mcca["clin"])
df_clin_multi$NEMO = unlist(nemo["clin"])
df_clin_multi$ClustOmicsMtoM = unlist(multi["clin"])

df_clin_multi = melt(df_clin_multi)

p <- ggplot(df_clin_multi, aes(x=cancers, y=value, fill=variable)) + geom_boxplot(fill="white") + 
  geom_dotplot(binaxis='y', stackdir='center', dotsize = 0.5, position=position_dodge(width=0.2)) + 
  labs(x ="Cancer", y = "# enriched clinical labels", fill = "Method") +
  scale_y_continuous(breaks=c(0,2,4,6,8))
  #+ scale_fill_brewer(palette = "Paired")

svg(file = paste("plots/Multi_Clin_Overall_Qual_Clust_", pval, ".svg", sep=""), width=10, height=7)
print(p)
dev.off()


# SINGLE
all_surv_single = c(pins_exp = pins_exp["surv"], snf_exp = snf_exp["surv"], mkl_exp = mkl_exp["surv"], nemo_exp = nemo_exp["surv"],
             pins_met = pins_met["surv"], snf_met = snf_met["surv"], mkl_met = mkl_met["surv"], nemo_met = nemo_met["surv"],
             pins_mirna = pins_mirna["surv"], snf_mirna = snf_mirna["surv"], mkl_mirna = mkl_mirna["surv"], nemo_mirna = nemo_mirna["surv"],
             single = single["surv"])
all_enriched_single = c(pins_exp = pins_exp["clin"], snf_exp = snf_exp["clin"], mkl_exp = mkl_exp["clin"], nemo_exp = nemo_exp["clin"],
                 pins_met = pins_met["clin"], snf_met = snf_met["clin"], mkl_met = mkl_met["clin"], nemo_met = nemo_met["clin"],
                 pins_mirna = pins_mirna["clin"], snf_mirna = snf_mirna["clin"], mkl_mirna = mkl_mirna["clin"], nemo_mirna = nemo_mirna["clin"],
                 single = single["clin"])
# SURV SINGLE
df_surv_single = as.data.frame(cancers)
df_surv_single$PINS_exp = unlist(pins_exp["surv"])
df_surv_single$SNF_exp = unlist(snf_exp["surv"])
df_surv_single$rMKL_exp = unlist(mkl_exp["surv"])
df_surv_single$NEMO_exp = unlist(nemo_exp["surv"])
df_surv_single$PINS_mirna = unlist(pins_mirna["surv"])
df_surv_single$SNF_mirna = unlist(snf_mirna["surv"])
df_surv_single$rMKL_mirna = unlist(mkl_mirna["surv"])
df_surv_single$NEMO_mirna = unlist(nemo_mirna["surv"])
df_surv_single$PINS_met = unlist(pins_met["surv"])
df_surv_single$SNF_met = unlist(snf_met["surv"])
df_surv_single$rMKL_met = unlist(mkl_met["surv"])
df_surv_single$NEMO_met = unlist(nemo_met["surv"])
df_surv_single$ClustOmicsStoM = unlist(single["surv"])

df_surv_single = melt(df_surv_single)
df_surv_single$method = gsub("_.*", "", df_surv_single$variable)
df_surv_single$omic = gsub(".*_", "", df_surv_single$variable)
df_surv_single$omic = gsub("exp", "Expression", df_surv_single$omic)
df_surv_single$omic = gsub("met", "Methylation", df_surv_single$omic)
df_surv_single$omic = gsub("mirna", "miRNA", df_surv_single$omic)

p <- ggplot(df_surv_single, aes(x=cancers, y=value, fill=omic)) + geom_boxplot(fill="white") + 
  geom_dotplot(binaxis='y', stackdir='center', dotsize = 0.5, position=position_dodge(width=0.2)) + labs(x ="Cancer", y = "-log10(survival)", fill = "Omic") + 
  geom_hline(yintercept=-log10(pval), linetype="dashed") + scale_fill_manual(values=c(Expression = "#F8766D", Methylation="#00BA38", miRNA="#619CFF", ClustOmicsStoM="#C77CFF"))
#+ scale_fill_brewer(palette = "Paired")
#c(Expression = "#F8766D", Methylation="#00BA38", miRNA="#619CFF", Multiomics="#C77CFF")

svg(file = paste("plots/Single_Surv_Overall_Qual_Clust_", pval, ".svg", sep=""), width=10, height=7)
print(p)
dev.off()

# CLIN SINGLE
df_clin_single = as.data.frame(cancers)
df_clin_single$PINS_exp = unlist(pins_exp["clin"])
df_clin_single$SNF_exp = unlist(snf_exp["clin"])
df_clin_single$rMKL_exp = unlist(mkl_exp["clin"])
df_clin_single$NEMO_exp = unlist(nemo_exp["clin"])
df_clin_single$PINS_mirna = unlist(pins_mirna["clin"])
df_clin_single$SNF_mirna = unlist(snf_mirna["clin"])
df_clin_single$rMKL_mirna = unlist(mkl_mirna["clin"])
df_clin_single$NEMO_mirna = unlist(nemo_mirna["clin"])
df_clin_single$PINS_met = unlist(pins_met["clin"])
df_clin_single$SNF_met = unlist(snf_met["clin"])
df_clin_single$rMKL_met = unlist(mkl_met["clin"])
df_clin_single$NEMO_met = unlist(nemo_met["clin"])
df_clin_single$ClustOmicsStoM = unlist(single["clin"])

df_clin_single = melt(df_clin_single)
df_clin_single$method = gsub("_.*", "", df_clin_single$variable)
df_clin_single$omic = gsub(".*_", "", df_clin_single$variable)
df_clin_single$omic = gsub("exp", "Expression", df_clin_single$omic)
df_clin_single$omic = gsub("met", "Methylation", df_clin_single$omic)
df_clin_single$omic = gsub("mirna", "miRNA", df_clin_single$omic)

p <- ggplot(df_clin_single, aes(x=cancers, y=value, fill=omic)) + geom_boxplot(fill="white") + 
  geom_dotplot(binaxis='y', stackdir='center', dotsize = 0.3, position=position_dodge(width=0.6)) + 
  labs(x ="Cancer", y = "# enriched clinical labels", fill = "Omic") +
  scale_y_continuous(breaks=c(0,2,4,6,8))

svg(file = paste("plots/Single_Clin_Overall_Qual_Clust_", pval, ".svg", sep=""), width=10, height=7)
print(p)
dev.off()
