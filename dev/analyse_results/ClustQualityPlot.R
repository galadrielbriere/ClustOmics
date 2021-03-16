library(ggplot2)
library(ggrepel)

pval = 0.01
clin_labels_per_cancer = c('AML'=4, 'BIC'=11, 'COAD'=10, 'GBM'=5, 'KIRC'=9, 'LIHC'=14, 'LUSC'=10, 'SARC'=4, 'SKCM'=9, 'OV'=3)
cancers <- c("AML", "BIC", "COAD", "GBM", "KIRC", "LIHC", "LUSC", "SARC", "SKCM", "OV")
rel = c('_MULTI_MCCA_NEMO_PINS_SNF_rMKL', "_EXP_MIRNA_MET_NEMO_PINS_SNF_rMKL_kmeans", "_EXP_MIRNA_MET_NEMO_PINS_SNF_rMKL_kmeans_all")

pins = list(surv = c(), clin = c())
mcca = list(surv = c(), clin = c())
snf = list(surv = c(), clin = c())
mkl = list(surv = c(), clin = c())
nemo = list(surv = c(), clin = c())
multi = list(surv = c(), clin = c())
single = list(surv = c(), clin = c())
single_all = list(surv = c(), clin = c())
coca_multi = list(surv = c(), clin = c())
coca_single = list(surv = c(), clin = c())
coca_single_all = list(surv = c(), clin = c())

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

for (cancer in cancers){
  pins_surv_file = paste("outOnlyMulti/survival", cancer, paste(cancer, "multiomics", "PINS.pval", sep="_"), sep="/")
  pins_surv = get_surv(pins_surv_file)
  pins$surv = append(pins$surv, pins_surv)
  
  pins_clin_file = paste("outOnlyMulti/clinical", cancer, paste(cancer, "multiomics", "PINS.pval", sep="_"), sep="/")
  pins_clin = get_clin(pins_clin_file)
  pins$clin = append(pins$clin, pins_clin)

  snf_surv_file = paste("outOnlyMulti/survival", cancer, paste(cancer, "multiomics", "SNF.pval", sep="_"), sep="/")
  snf_surv = get_surv(snf_surv_file)
  snf$surv = append(snf$surv, snf_surv)
  
  snf_clin_file = paste("outOnlyMulti/clinical", cancer, paste(cancer, "multiomics", "SNF.pval", sep="_"), sep="/")
  snf_clin = get_clin(snf_clin_file)
  snf$clin = append(snf$clin, snf_clin)
  
  mcca_surv_file = paste("outOnlyMulti/survival", cancer, paste(cancer, "multiomics", "MCCA.pval", sep="_"), sep="/")
  mcca_surv = get_surv(mcca_surv_file)
  mcca$surv = append(mcca$surv, mcca_surv)
  
  mcca_clin_file = paste("outOnlyMulti/clinical", cancer, paste(cancer, "multiomics", "MCCA.pval", sep="_"), sep="/")
  mcca_clin = get_clin(mcca_clin_file)
  mcca$clin = append(mcca$clin, mcca_clin)
  
  nemo_surv_file = paste("outOnlyMulti/survival", cancer, paste(cancer, "multiomics", "NEMO.pval", sep="_"), sep="/")
  nemo_surv = get_surv(nemo_surv_file)
  nemo$surv = append(nemo$surv, nemo_surv)
  
  nemo_clin_file = paste("outOnlyMulti/clinical", cancer, paste(cancer, "multiomics", "NEMO.pval", sep="_"), sep="/")
  nemo_clin = get_clin(nemo_clin_file)
  nemo$clin = append(nemo$clin, nemo_clin)
  
  mkl_surv_file = paste("outOnlyMulti/survival", cancer, paste(cancer, "multiomics", "rMKL.pval", sep="_"), sep="/")
  mkl_surv = get_surv(mkl_surv_file)
  mkl$surv = append(mkl$surv, mkl_surv)
  
  mkl_clin_file = paste("outOnlyMulti/clinical", cancer, paste(cancer, "multiomics", "rMKL.pval", sep="_"), sep="/")
  mkl_clin = get_clin(mkl_clin_file)
  mkl$clin = append(mkl$clin, mkl_clin)
}

for (cancer in cancers) {
  coca_multi_surv_file = paste("outOnlyMulti/survival/coca", paste(cancer, "multiomics", "COCA.pval", sep="_"), sep="/")
  coca_multi_surv = get_surv(coca_multi_surv_file)
  coca_multi$surv = append(coca_multi$surv, coca_multi_surv)
  
  coca_multi_clin_file = paste("outOnlyMulti/clinical/coca", paste(cancer, "multiomics", "COCA.pval", sep="_"), sep="/")
  coca_multi_clin = get_clin(coca_multi_clin_file)
  coca_multi$clin = append(coca_multi$clin, coca_multi_clin)
  
  coca_single_surv_file = paste("outOnlyMulti/survival/coca", paste(cancer, "singleToMulti", "COCA.pval", sep="_"), sep="/")
  coca_single_surv = get_surv(coca_single_surv_file)
  coca_single$surv = append(coca_single$surv, coca_single_surv)
  
  coca_single_clin_file = paste("outOnlyMulti/clinical/coca", paste(cancer, "singleToMulti", "COCA.pval", sep="_"), sep="/")
  coca_single_clin = get_clin(coca_single_clin_file)
  coca_single$clin = append(coca_single$clin, coca_single_clin)
  
  coca_single_all_surv_file = paste("outAll/survival/coca", paste(cancer, "singleToMulti", "COCA.pval", sep="_"), sep="/")
  coca_single_all_surv = get_surv(coca_single_all_surv_file)
  coca_single_all$surv = append(coca_single_all$surv, coca_single_all_surv)
  
  coca_single_all_clin_file = paste("outAll/clinical/coca", paste(cancer, "singleToMulti", "COCA.pval", sep="_"), sep="/")
  coca_single_all_clin = get_clin(coca_single_all_clin_file)
  coca_single_all$clin = append(coca_single_all$clin, coca_single_all_clin)
}

for (cancer in cancers){
  rel_name_multi = paste(cancer, rel[1], sep="")
  rel_name_single = paste(cancer, rel[2], sep="")
  rel_name_single_all = paste(cancer, rel[3], sep="")
  
  multi_surv_file = paste("outOnlyMulti/survival", cancer, paste(cancer, rel_name_multi, "ClustOmicsClustering.pval", sep="."), sep="/")
  multi_surv = get_surv(multi_surv_file)
  multi$surv = append(multi$surv, multi_surv)
  
  single_surv_file = paste("outOnlyMulti/survival", cancer, paste(cancer, rel_name_single, "ClustOmicsClustering.pval", sep="."), sep="/")
  single_surv = get_surv(single_surv_file)
  single$surv = append(single$surv, single_surv)
  
  multi_clin_file = paste("outOnlyMulti/clinical", cancer, paste(cancer, rel_name_multi, "ClustOmicsClustering.pval", sep="."), sep="/")
  multi_clin = get_clin(multi_clin_file)
  multi$clin = append(multi$clin, multi_clin)
  
  single_clin_file = paste("outOnlyMulti/clinical", cancer, paste(cancer, rel_name_single, "ClustOmicsClustering.pval", sep="."), sep="/")
  single_clin = get_clin(single_clin_file)
  single$clin = append(single$clin, single_clin)
  
  single_all_surv_file = paste("outAll/survival", cancer, paste(cancer, rel_name_single_all, "ClustOmicsClustering.pval", sep="."), sep="/")
  single_all_surv = get_surv(single_all_surv_file)
  single_all$surv = append(single_all$surv, single_all_surv)
  
  single_all_clin_file = paste("outAll/clinical", cancer, paste(cancer, rel_name_single_all, "ClustOmicsClustering.pval", sep="."), sep="/")
  single_all_clin = get_clin(single_all_clin_file)
  single_all$clin = append(single_all$clin, single_all_clin)
}

nb_surv = list()
mean_enriched = list()

all_surv = c(pins = pins["surv"], snf = snf["surv"], mkl = mkl["surv"], mcca = mcca["surv"], nemo = nemo["surv"],
             multi = multi["surv"], single = single["surv"], single_all = single_all["surv"],
             coca_multi = coca_multi['surv'], coca_single = coca_single['surv'], coca_single_all = coca_single_all['surv'])
all_enriched = c(pins = pins["clin"], snf = snf["clin"], mkl = mkl["clin"], mcca = mcca["clin"], nemo = nemo["clin"],
                 multi = multi["clin"], single = single["clin"], single_all = single_all["clin"],
                 coca_multi = coca_multi['clin'], coca_single = coca_single['clin'], coca_single_all = coca_single_all['clin'])

for (res in all_surv) {
  nb_surv = append(nb_surv, length(which(res >= -log10(pval))))
}
names(nb_surv) = names(all_surv)

for (res in all_enriched) {
  mean_enriched = append(mean_enriched, (sum(unlist(res))))
}
names(mean_enriched) = names(all_enriched)

names(nb_surv) = gsub("\\..*", "", names(nb_surv))
names(mean_enriched) = gsub("\\..*", "", names(mean_enriched))
df <- rbind(as.data.frame(nb_surv), as.data.frame(mean_enriched))
row.names(df) = c("survival", "enriched")
df = as.data.frame(t(df))
df$labels = c("PINS", "SNF", "rMKL", "MultiCCA", "NEMO", "ClustOmics MtoM", "ClustOmics StoM OnlyMulti",
           "ClustOmics StoM All", "COCA MtoM", "COCA StoM OnlyMulti", "COCA StoM All")
df$color = as.factor(as.character(c("Inputs MtoM","Inputs MtoM","Inputs MtoM","Inputs MtoM","Inputs MtoM",
                                    "MtoM", "StoM OnlyMulti","StoM All", "MtoM",
                                    "StoM OnlyMulti", "StoM All")))
df$shape = as.factor(as.character(c("Inputs MtoM","Inputs MtoM","Inputs MtoM","Inputs MtoM","Inputs MtoM",
                                    "ClustOmics", "ClustOmics","ClustOmics", "COCA",
                                    "COCA", "COCA")))
colors = c("#e5e5e5", "#006d77", "#e29578", "#f6bd60")
names(colors) = c('Inputs MtoM', "MtoM", "StoM OnlyMulti", "StoM All")
shapes = c(21, 22, 23)
names(shapes) = c('Inputs MtoM', "ClustOmics", "COCA")

p = ggplot(data=df, aes(x=survival, y=enriched)) 
p = p + geom_point(size=4, aes(shape=shape, fill=color)) 
p = p + geom_label_repel(aes(label=labels, fill=color), box.padding = 0.35, point.padding = 0.7, segment.color = 'grey50')
p = p + guides(fill = guide_legend(override.aes = aes(label = "", shape=""))) 
p = p +  scale_fill_manual(name = "Integration strategy",values=alpha(colors,0.5))
p = p + scale_shape_manual(name = "Tool", values=shapes)
p = p + labs(title = "Overall cancer clusterings biological quality (10 cancer types, 79 clinical labels)")
p = p + xlab(paste("Number of significant survival logrank pvalues (pval<", pval, ")", sep="")) 
p = p + ylab(paste("Number of enriched (pval<", pval, ") clinical labels overall cancers", sep=""))
p = p +  theme_bw() + scale_y_continuous(breaks=seq(16, 28, 1)) 
p

pdf(file = paste("outOnlyMulti/plots/Overall_Qual_Clust_", pval, ".pdf", sep=""), width=10, height=7)
print(p)
dev.off()

