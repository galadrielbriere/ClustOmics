library(parallel)
library(pheatmap)
library(RColorBrewer)
library(edgeR)
library(viridis)
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)

setwd("describeBIC")

save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

cancer = "BIC"

########## EXPRESSION
clustering = "../outAll/results/BIC/BIC.BIC_EXP_MIRNA_MET_NEMO_PINS_SNF_rMKL_kmeans_all.ClustOmicsClustering"
clustering = read.table(clustering, header=T)
clustering$Patient = as.character(clustering$Patient)
order = order(clustering$Cluster)
clustering = clustering[order,]

exp = "../raw_data/BIC/exp"
exp = read.table(exp)
exp = t(exp)
exp = exp[which(rownames(exp) %in% clustering$Patient),]

clustering = clustering[which(clustering$Patient %in% rownames(exp)),]

exp = exp[clustering$Patient,]
exp = log2(1+exp)

# Only Symbols
exp = exp[,which(!grepl("\\?",colnames(exp)))]

scaled = scale(exp)
boxplot(scaled[,1:50],col = rainbow(50))

scaled = data.frame(scaled)

bad <- sapply(scaled, function(x) all(is.nan(x)))
bad = bad[which(bad)]

scaled = scaled[, which(!(colnames(scaled) %in% names(bad)))]

unique(rownames(scaled) == clustering$Patient)
scaled$cluster = as.factor(unlist(clustering$Cluster))

# ANOVA
baseformula <- " ~ cluster"
pvals = as.numeric(mclapply(1:(ncol(scaled)-1), function(i) {
  formula <- paste(colnames(scaled)[i], baseformula, sep="")
  p <- kruskal.test(as.formula(formula), data=scaled)$p.value
  return(p)
}, mc.cores=10))

names(pvals) = colnames(scaled)[-ncol(scaled)]

# Adjust pvals
adjp = p.adjust (pvals, method="fdr") 

# Save results
adjp = data.frame(adjp)
adjp$Gene = names(pvals)
adjp = adjp[,c(2,1)]
names(adjp) = c("Gene", "FDR")
write.table(adjp, "FDRpvalsClustOmicsExpression.csv", row.names = F, col.names = T, quote = F, sep = "\t")

# Keep top pvals
quantile(adjp$FDR)
fdr = adjp$FDR
names(fdr) = adjp$Gene
order = order(fdr)
ofdr = fdr[order]

top = ofdr[1:1000]

## HEATMAP
data = scaled
# For top pvals
data = data[,which(colnames(data) %in% names(top))]
clusters = scaled$cluster
clusters = gsub("1", "A", clusters)
clusters = gsub("2", "B", clusters)
clusters = gsub("3", "C", clusters)
clusters = gsub("4", "D", clusters)
clusters = gsub("6", "E", clusters)
clusters = gsub("8", "F", clusters)
clusters = as.factor(clusters)
# Clusters
annotation_row = data.frame(
  'Consensus Cluster' = clusters
)
rownames(annotation_row) = rownames(data)
# Cluster colors
annot_col = brewer.pal(length(unique(scaled$cluster)), "Set1")
annot_col = c("#808000","#f58231","#ff0000","#FFFF00","#00ffff","#0000ff")
names(annot_col) = unique(clusters)
ann_colors = list(
  'Consensus Cluster' = annot_col
)

# Heatmap colors
quantile_breaks <- function(xs, n = 10) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
  breaks[!duplicated(breaks)]
}
mat_breaks <- quantile_breaks(as.matrix(data), n = 10)

# Heatmap
p = pheatmap(data, fontsize=5, show_rownames=F, show_colnames=F, cluster_rows=F,
             cutree_cols = 6,
             annotation_row=annotation_row, annotation_colors=ann_colors,
             annotation_names_row=T,
             breaks = mat_breaks, color = inferno(length(mat_breaks) - 1) ) 

clustgenes = cutree(p$tree_col, k=6)
oclustgenes = clustgenes[order(clustgenes)]
write.table(oclustgenes, "geneClust.csv", row.names = T, col.names = F, quote=F, sep="\t")

unique(names(data) == names(clustgenes))
annotationGenes =  data.frame(paste("X", unlist(clustgenes), sep=""))
rownames(annotationGenes) = names(clustgenes)
names(annotationGenes) = c("Gene Cluster")

annotation_row = data.frame(clusters)
names(annotation_row) = c("Consensus Cluster")
rownames(annotation_row) = rownames(data)

# PAM 50
classif = read.table("tcga_vs_myPAM50.txt", header=T)
pam = classif$TCGAprediction
names(pam) = rownames(classif)
pamNoNa = lapply(pam, function(x) {
  if (is.na(x)) {
    return(classif[names(x),2])
  } else { return(x) }
})
pamNoNa = unlist(pamNoNa)
names(pamNoNa) = names(pam)
pamNoNa = pamNoNa[rownames(annotation_row)]

annotation_row$PAM50 = unlist(pamNoNa)
annotation_row$PAM50 = gsub("Lum", "Luminal ", annotation_row$PAM50)
annotation_row$PAM50 = gsub("Normal", "Normal-like ", annotation_row$PAM50)

annot_col_genes = brewer.pal(length(unique(annotationGenes$`Gene Cluster`)), "Set2")
names(annot_col_genes) = unique(annotationGenes$`Gene Cluster`)
annot_col = brewer.pal(6, "Dark2")
pam_col = c("red", "#808000", "#ffe119", "#4363d8", "grey")
names(pam_col) = unique(annotation_row$PAM50)
names(annot_col) = unique(clusters)

ann_colors = list(
  'Gene Cluster' = annot_col_genes,
  'Consensus Cluster' = annot_col,
  'PAM50' = pam_col
)

p = pheatmap(t(data), fontsize=5, show_rownames=F, show_colnames=F, cluster_rows=T, cluster_cols = F,
             cutree_rows = 6,
             annotation_col=annotation_row, annotation_colors=ann_colors,
             annotation_row=annotationGenes, 
             annotation_names_row=T, annotation_names_col = T,
             breaks = mat_breaks, color = inferno(length(mat_breaks) - 1) ) #, cutree_rows = 6)

save_pheatmap_pdf(p, "expressionHeatmap.pdf")

# Gene clusters analysis
clust1 = clustgenes[which(clustgenes == 1)]
clust2 = clustgenes[which(clustgenes == 2)]
clust3 = clustgenes[which(clustgenes == 3)]
clust4 = clustgenes[which(clustgenes == 4)]
clust5 = clustgenes[which(clustgenes == 5)]
clust6 = clustgenes[which(clustgenes == 6)]

geneClusters = list(X1=gsub(".*\\.", "", names(clust1)),
                    X2=gsub(".*\\.", "", names(clust2)),
                    X3=gsub(".*\\.", "", names(clust3)),
                    X4=gsub(".*\\.", "", names(clust4)),
                    X5=gsub(".*\\.", "", names(clust5)),
                    X6=gsub(".*\\.", "", names(clust6)))
xx <- compareCluster(geneClusters, fun="enrichGO", OrgDb = org.Hs.eg.db,
                     keyType       = 'ENTREZID',
                     ont           = "BP",
                     pAdjustMethod = "fdr",
                     pvalueCutoff  = 0.01,
                     qvalueCutoff  = 0.05)

pdf(file="GOenriched.pdf", width=10, height=10)
dotplot(xx, showCategory=10, title="GO - Biological Process")
dev.off()

########## HEATMAP miRNA
clustering = "../outAll/results/BIC/BIC.BIC_EXP_MIRNA_MET_NEMO_PINS_SNF_rMKL_kmeans_all.ClustOmicsClustering"
clustering = read.table(clustering, header=T)
clustering$Patient = as.character(clustering$Patient)
order = order(clustering$Cluster)
clustering = clustering[order,]

mirna = "../raw_data/BIC/mirna"
mirna = read.table(mirna)
mirna = t(mirna)
mirna = mirna[which(rownames(mirna) %in% clustering$Patient),]

clustering = clustering[which(clustering$Patient %in% rownames(mirna)),]

mirna = mirna[clustering$Patient,]
mirna = log2(1+mirna)

scaled = scale(mirna)
boxplot(scaled[,1:50],col = rainbow(50))

scaled = data.frame(scaled)

bad <- sapply(scaled, function(x) all(is.nan(x)))
bad = bad[which(bad)]

scaled = scaled[, which(!(colnames(scaled) %in% names(bad)))]

unique(rownames(scaled) == clustering$Patient)
scaled$cluster = as.factor(unlist(clustering$Cluster))

# ANOVA
baseformula <- " ~ cluster"
pvalsMI = as.numeric(mclapply(1:(ncol(scaled)-1), function(i) {
  formula <- paste(colnames(scaled)[i], baseformula, sep="")
  p <- kruskal.test(as.formula(formula), data=scaled)$p.value
  return(p)
}, mc.cores=10))

names(pvalsMI) = colnames(scaled)[-ncol(scaled)]

# Adjust pvals
adjpMI = p.adjust (pvalsMI, method="fdr") 

# Save results
adjpMI = data.frame(adjpMI)
adjpMI$Gene = names(pvalsMI)
adjpMI = adjpMI[,c(2,1)]
names(adjpMI) = c("Gene", "FDR")
write.table(adjpMI, "FDRpvalsClustOmicsmiRNA.csv", row.names = F, col.names = T, quote = F, sep = "\t")

# Keep top pvals
quantile(adjpMI$FDR)
top = adjpMI[which(adjpMI$FDR <= 0.01),1]

# Plot heatmap
data = scaled
# For top pvals
data = data[,which(colnames(data) %in% top)]

clusters = scaled$cluster
clusters = gsub("1", "A", clusters)
clusters = gsub("2", "B", clusters)
clusters = gsub("3", "C", clusters)
clusters = gsub("4", "D", clusters)
clusters = gsub("6", "E", clusters)
clusters = gsub("8", "F", clusters)
clusters = as.factor(clusters)

# Clusters
annotation_row = data.frame(clusters)
names(annotation_row) = c("Consensus Cluster")
rownames(annotation_row) = rownames(data)

# PAM 50
classif = read.table("tcga_vs_myPAM50.txt", header=T)
pam = classif$TCGAprediction
names(pam) = rownames(classif)
pamNoNa = lapply(pam, function(x) {
  if (is.na(x)) {
    return(classif[names(x),2])
  } else { return(x) }
})
pamNoNa = unlist(pamNoNa)
names(pamNoNa) = names(pam)
pamNoNa = pamNoNa[rownames(annotation_row)]

annotation_row$PAM50 = unlist(pamNoNa)
annotation_row$PAM50 = gsub("Lum", "Luminal ", annotation_row$PAM50)
annotation_row$PAM50 = gsub("Normal", "Normal-like ", annotation_row$PAM50)

annot_col = brewer.pal(6, "Dark2")
pam_col = c("red", "#808000", "#ffe119", "#4363d8", "black", "grey")
names(pam_col) = unique(annotation_row$PAM50)
names(annot_col) = unique(clusters)

ann_colors = list(
  'Consensus Cluster' = annot_col,
  'PAM50' = pam_col
)

mat_breaks <- quantile_breaks(as.matrix(data), n = 10)

p=pheatmap(t(data), fontsize=5, show_rownames=F, show_colnames=F, cluster_rows=T, cluster_cols = F,
         annotation_col=annotation_row, annotation_colors=ann_colors,
         annotation_names_row=T, annotation_names_col = T,
         breaks = mat_breaks, color = inferno(length(mat_breaks) - 1) ) 

save_pheatmap_pdf(p, "miRNAHeatmap.pdf")

########## METHYLATION
clustering = "../outAll/results/BIC/BIC.BIC_EXP_MIRNA_MET_NEMO_PINS_SNF_rMKL_kmeans_all.ClustOmicsClustering"
clustering = read.table(clustering, header=T)
clustering$Patient = as.character(clustering$Patient)
order = order(clustering$Cluster)
clustering = clustering[order,]

methy = "../raw_data/BIC/methy"
methy = read.table(methy)
methy = t(methy)
methy = methy[which(rownames(methy) %in% clustering$Patient),]

clustering = clustering[which(clustering$Patient %in% rownames(methy)),]

methy = methy[clustering$Patient,]

scaled = scale(methy)
boxplot(scaled[,1:50],col = rainbow(50))

scaled = data.frame(scaled)

bad <- sapply(scaled, function(x) all(is.nan(x)))
bad = bad[which(bad)]

scaled = scaled[, which(!(colnames(scaled) %in% names(bad)))]

unique(rownames(scaled) == clustering$Patient)
scaled$cluster = as.factor(unlist(clustering$Cluster))

# ANOVA
baseformula <- " ~ cluster"
pvalsMET = as.numeric(mclapply(1:(ncol(scaled)-1), function(i) {
  formula <- paste(colnames(scaled)[i], baseformula, sep="")
  p <- kruskal.test(as.formula(formula), data=scaled)$p.value
  return(p)
}, mc.cores=10))

names(pvalsMET) = colnames(scaled)[-ncol(scaled)]

# Adjust pvals
adjpMET = p.adjust (pvalsMET, method="fdr") 

# Save results
adjpMET = data.frame(adjpMET)
adjpMET$Gene = names(pvalsMET)
adjpMET = adjpMET[,c(2,1)]
names(adjpMET) = c("Gene", "FDR")
write.table(adjpMET, "FDRpvalsClustOmicsMethy.csv", row.names = F, col.names = T, quote = F, sep = "\t")

# # Keep top pvals
quantile(adjpMET$FDR)
top = adjpMET[which(adjpMET$FDR <= 0.01),1]

fdr = adjpMET$FDR
names(fdr) = adjpMET$Gene
order = order(fdr)
ofdr = fdr[order]

top = ofdr[1:1000]

# Plot heatmap
data = scaled
# For top pvals
data = data[,which(colnames(data) %in% names(top))]

clusters = scaled$cluster
clusters = gsub("1", "A", clusters)
clusters = gsub("2", "B", clusters)
clusters = gsub("3", "C", clusters)
clusters = gsub("4", "D", clusters)
clusters = gsub("6", "E", clusters)
clusters = gsub("8", "F", clusters)
clusters = as.factor(clusters)

# Clusters
annotation_row = data.frame(clusters)
names(annotation_row) = c("Consensus Cluster")
rownames(annotation_row) = rownames(data)

# PAM 50
classif = read.table("tcga_vs_myPAM50.txt", header=T)
pam = classif$TCGAprediction
names(pam) = rownames(classif)
pamNoNa = lapply(pam, function(x) {
  if (is.na(x)) {
    return(classif[names(x),2])
  } else { return(x) }
})
pamNoNa = unlist(pamNoNa)
names(pamNoNa) = names(pam)
pamNoNa = pamNoNa[rownames(annotation_row)]

annotation_row$PAM50 = unlist(pamNoNa)
annotation_row$PAM50 = gsub("Lum", "Luminal ", annotation_row$PAM50)
annotation_row$PAM50 = gsub("Normal", "Normal-like ", annotation_row$PAM50)

annot_col = brewer.pal(6, "Dark2")
pam_col = c("red", "#808000", "#ffe119", "#4363d8", "black", "grey")
names(pam_col) = unique(annotation_row$PAM50)
names(annot_col) = unique(clusters)

ann_colors = list(
  'Consensus Cluster' = annot_col,
  'PAM50' = pam_col
)

mat_breaks <- quantile_breaks(as.matrix(data), n = 10)

p=pheatmap(t(data), fontsize=5, show_rownames=F, show_colnames=F, cluster_rows=T, cluster_cols = F,
           annotation_col=annotation_row, annotation_colors=ann_colors,
           annotation_names_row=T, annotation_names_col = T,
           breaks = mat_breaks, color = inferno(length(mat_breaks) - 1) ) 

save_pheatmap_pdf(p, "methylationHeatmap.pdf")