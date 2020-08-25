library("FactoMineR")
library("factoextra")
library("genefilter")
library("affy")
library("pca3d")
library("optparse")

option_list = list(
  make_option(c("-c", "--cancer"), type="character", default=NULL, 
              help="Cancer", metavar="character"),
  make_option(c("-o", "--out_plot"), type="character", default=NULL, 
              help="output file name", metavar="character"),
  make_option(c("-r", "--clust_results"), type="character", default=NULL, 
              help="Clustering result file", metavar="character")
) 

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

cancer = opt$cancer
clustering = opt$clust_results
out_plot = opt$out_plot

data_path <- paste("./raw_data", cancer, sep="/")
exp_data <- paste(data_path, "exp", sep="/")
mirna_data <- paste(data_path, "mirna", sep="/")
met_data <- paste(data_path, "methy", sep="/")

exp <- read.table(exp_data, header = T)
mirna <- read.table(mirna_data, header = T)
met <- read.table(met_data, header = T)

clustering = read.table(clustering, header=F, skip = 1, sep="\t")
clust = clustering$V2
names(clust) = clustering$V1

exp = exp[,which(colnames(exp) %in% names(clust))]
met = met[,which(colnames(met) %in% names(clust))]
mirna = mirna[,which(colnames(mirna) %in% names(clust))]

exp = log(1+exp)
exp = t(scale(t(exp)))
exp = exp[complete.cases(exp), ]
expSet = ExpressionSet(assayData = as.matrix(exp))
expSet_filt = varFilter(expSet, var.func=IQR, var.cutoff=0.9, filterByQuantile=TRUE)

mirna = log(1+mirna)
mirna = t(scale(t(mirna)))
mirna = mirna[complete.cases(mirna), ]
if (dim(mirna)[1]>2000) {
  vars = apply(mirna, 1, var)
  threshold = vars[order(vars, decreasing = T)][2000]
  mirna = mirna[vars >= threshold,]
}

met = t(scale(t(met)))
met = met[complete.cases(met),]
if (dim(met)[1]>2000) {
  vars = apply(met, 1, var)
  threshold = vars[order(vars, decreasing = T)][2000]
  met = met[vars >= threshold,]
}

#### EXP
res.pca.exp <- PCA(t(exprs(expSet_filt)), scale.unit = F, graph = FALSE)
fviz_eig(res.pca.exp, addlabels = TRUE, ylim = c(0, 50))
cl_pca_exp = factor(clust[colnames(exp)])
fviz_pca_ind (res.pca.exp, geom = "point", axes = c(1, 2), mean.point = FALSE, addEllipses = T,
              repel = TRUE, title = "PCA expression", col.ind =  cl_pca_exp, legend.title = "Clusters")
fviz_pca_ind (res.pca.exp, geom = "point", axes = c(2, 3), mean.point = FALSE, addEllipses = T,
              repel = TRUE, title = "PCA expression", col.ind =  cl_pca_exp, legend.title = "Clusters")

title_plot = paste(c("Dim1 (", "Dim2 (", "Dim3 ("), round(res.pca.exp$eig[c(1:3), 'percentage of variance'], digits=2), rep("%)", 3), sep="")
pca3d(res.pca.exp$ind$coord, group=cl_pca_exp,
      show.ellipses=TRUE, ellipse.ci=0.75, show.plane=FALSE, axe.titles = title_plot)
snapshotPCA3d(file=paste(out_plot, "_exp.png", sep=""))

#### MIRNA
res.pca.mirna <- PCA(t(mirna), scale.unit = F, graph = FALSE)
fviz_eig(res.pca.mirna, addlabels = TRUE, ylim = c(0, 50))
cl_pca_mirna = factor(clust[colnames(mirna)])
fviz_pca_ind (res.pca.mirna, geom = "point", axes = c(1, 2), mean.point = FALSE, addEllipses = T,
              repel = TRUE, title = "PCA mirna", col.ind =  cl_pca_mirna, legend.title = "Clusters")
fviz_pca_ind (res.pca.mirna, geom = "point", axes = c(2, 3), mean.point = FALSE, addEllipses = T,
              repel = TRUE, title = "PCA mirna", col.ind =  cl_pca_mirna, legend.title = "Clusters")

title_plot = paste(c("Dim1 (", "Dim2 (", "Dim3 ("), round(res.pca.mirna$eig[c(1:3), 'percentage of variance'], digits=2), rep("%)", 3), sep="")
pca3d(res.pca.mirna$ind$coord, group=cl_pca_mirna,
      show.ellipses=TRUE, ellipse.ci=0.75, show.plane=FALSE, axe.titles = title_plot)
snapshotPCA3d(file=paste(out_plot, "_mirna.png", sep=""))

#### MET
res.pca.met <- PCA(t(met), scale.unit = F, graph = FALSE)
fviz_eig(res.pca.met, addlabels = TRUE, ylim = c(0, 50))
cl_pca_met = factor(clust[colnames(met)])
fviz_pca_ind (res.pca.met, geom = "point", axes = c(1, 2), mean.point = FALSE, addEllipses = T,
              repel = TRUE, title = "PCA methylation", col.ind =  cl_pca_met, legend.title = "Clusters")
fviz_pca_ind (res.pca.met, geom = "point", axes = c(2, 3), mean.point = FALSE, addEllipses = T,
              repel = TRUE, title = "PCA methylation", col.ind =  cl_pca_met, legend.title = "Clusters")

title_plot = paste(c("Dim1 (", "Dim2 (", "Dim3 ("), round(res.pca.met$eig[c(1:3), 'percentage of variance'], digits=2), rep("%)", 3), sep="")
pca3d(res.pca.met$ind$coord, group=cl_pca_met,
      show.ellipses=TRUE, ellipse.ci=0.75, show.plane=FALSE, axe.titles = title_plot)
snapshotPCA3d(file=paste(out_plot, "_met.png", sep=""))
