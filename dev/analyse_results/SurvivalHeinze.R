#!/usr/bin/env Rscript

library(optparse)
library(logrankHeinze)
library(survival)
library(survminer)

option_list = list(
  make_option(c("-c", "--clust"), type="character", default=NULL, 
              help="Clustering file name", metavar="character"),
  make_option(c("-s", "--surv_file"), type="character", default=NULL, 
              help="Survival data", metavar="character"),
  make_option(c("-o", "--out_survival"), type="character", default=NULL, 
              help="output file name (survival pvalue)", metavar="character"),
  make_option(c("-f", "--out_fig"), type="character", default=NULL, 
              help="output survival plot", metavar="character")
) 

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

clust = read.table(opt$clust, header=T, sep="\t")
clustering <- gsub(".*_", "", clust$Cluster)
clustering = gsub("unclassified", 0, clustering)
clustering = as.numeric(clustering)
names(clustering) <- clust$Patient

survival.data = read.table(opt$surv_file, header = TRUE)

patient.names = names(clustering)
patient.names.in.file = as.character(survival.data[, 1])
patient.names.in.file = toupper(gsub('-', '\\.', substring(patient.names.in.file, 1, 15)))

if (!all(patient.names %in% patient.names.in.file)) {
  patient.names <- toupper(gsub('-', '\\.', substring(patient.names, 1, 12)))
  stopifnot(all(patient.names %in% patient.names.in.file))
}

indices = match(patient.names, patient.names.in.file)
ordered.survival.data = survival.data[indices,]
ordered.survival.data["cluster"] <- clustering
rownames(ordered.survival.data) = ordered.survival.data$PatientID
ordered.survival.data = ordered.survival.data[,-1]

# Change cluster names to fit HeinzeSurvival code
raw_clusters = unique(ordered.survival.data$cluster)
new_clusters = c(1:length(raw_clusters))
corres = new_clusters
names(corres) = raw_clusters
ordered.survival.data$new_clust = corres[as.character(ordered.survival.data$cluster)]

ordered.survival.data = ordered.survival.data[which(!is.na(ordered.survival.data$Survival)&!is.na(ordered.survival.data$Death)),]

heinze_results = logrankHeinze::logrank.heinze(Survival + Death ~ new_clust, data=ordered.survival.data, max.num.perms = 100000)
pval <- heinze_results$pvalue

write.table(pval, opt$out_survival, col.names=F, row.names=F, sep="\t", quote=F)

asymp_results = survfit(Surv(Survival, Death) ~ cluster, data=ordered.survival.data)
svg(file = opt$out_fig)
plot = ggsurvplot(asymp_results, data = ordered.survival.data, risk.table = F, pval = TRUE)
print(plot)
dev.off()

