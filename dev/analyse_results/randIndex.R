library(pdfCluster)
library(RColorBrewer)
library(gplots)
library(pheatmap)
library(ggplot2)
library(viridis)

save_pheatmap_png <- function(x, filename, width=1200, height=1000, res = 150) {
  png(filename, width = width, height = height, res = res)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

cancers = c("AML", "BIC", "COAD", "GBM", "KIRC", "LIHC", "LUSC", "OV", "SARC", "SKCM")

# SINGLE TO MULTI
methods = c("PINS", "NEMO", "SNF", "rMKL", "MCCA", "kmeans")
datatypes = c("expression", "mirna", "methylation")

clusterings = c()
for (datatype in datatypes) {
  for (method in methods) {
    if ((method == "MCCA") && (datatype != "multiomics")) {
      next
    } else {
      run = paste(datatype, method, sep="_")
      clusterings = c(clusterings, run)
    }
  }
}
clusterings = c(clusterings, "ClustOmics StoM", "COCA StoM")
combinations = expand.grid(clusterings, clusterings)

all_rands = data.frame()

for (cancer in cancers) {
  out_plot = paste("outAll/plots/", cancer, "/", cancer, "_heatmap_rand_index_tree_SingleToMulti.png", sep="")

  combinations$rands = apply(combinations, 1, function(x) {
    if (x[[1]] == "ClustOmics StoM") {
      file1 = paste("outAll/results/", cancer, "/", cancer, ".", cancer, "_EXP_MIRNA_MET_NEMO_PINS_SNF_rMKL_kmeans_all.ClustOmicsClustering", sep="")
    } else if (x[[1]] == "COCA StoM") {
      file1 = paste("dataAll/coca/", cancer, "_singleToMulti_COCA.clst", sep="")
    } else {
      file1 = paste("dataAll/", cancer, "/", cancer, "_", x[[1]], ".clst", sep="")
    }
    
    if (x[[2]] == "ClustOmics StoM") {
      file2 = paste("outAll/results/", cancer, "/", cancer, ".", cancer, "_EXP_MIRNA_MET_NEMO_PINS_SNF_rMKL_kmeans_all.ClustOmicsClustering", sep="")
    } else if (x[[2]] == "COCA StoM") {
      file2 = paste("dataAll/coca/", cancer, "_singleToMulti_COCA.clst", sep="")
    } else {
      file2 = paste("dataAll/", cancer, "/", cancer, "_", x[[2]], ".clst", sep="")
    }

    c1 = read.table(file1, header = T)
    c2 = read.table(file2, header = T)

    
    cl1 = c1$Cluster
    names(cl1) = c1$Patient
    cl2 = c2$Cluster
    names(cl2) = c2$Patient
    
    multi = intersect(names(cl1), names(cl2))
    cl1 = cl1[multi]
    cl2 = cl2[multi]
    
    return(adj.rand.index(cl1, cl2))
  
  })
  
  combinations = combinations[which(!is.na(combinations$rands)),]
  df = matrix(ncol = length(clusterings), nrow = length(clusterings), data = combinations$rands)
  colnames(df) = unique(combinations$Var2)
  rownames(df) = unique(combinations$Var1)
  
  rownames(df) = gsub("ClustOmics StoM", "multiomics_ClustOmics StoM", rownames(df))
  rownames(df) = gsub("COCA StoM", "multiomics_COCA StoM", rownames(df))
  colnames(df) = row.names(df)
  
  my_group = as.factor(gsub("_.*", "", rownames(df)))
  names(my_group) = rownames(df)
  my_group = as.data.frame(my_group)
  names(my_group) = "Omic"
  my_group$Omic = gsub("expression", "Expression", my_group$Omic)
  my_group$Omic = gsub("mirna", "miRNA", my_group$Omic)
  my_group$Omic = gsub("methylation", "Methylation", my_group$Omic)
  my_group$Omic = gsub("multiomics", "Multiomics", my_group$Omic)
  my_colour = list(Omic = c(Expression = "#F8766D", Methylation="#00BA38", miRNA="#619CFF", Multiomics="#C77CFF"))
  p = pheatmap(df, annotation_row = my_group, annotation_col = my_group, annotation_colors = my_colour, color=magma(10))

  save_pheatmap_png(p, out_plot)
  
  combinations$cancer = cancer
  all_rands = rbind(all_rands, combinations)
}

clusto = all_rands[which((all_rands$Var1 == "ClustOmics StoM") & (all_rands$Var2 != "ClustOmics StoM")),]
clusto$omic = gsub("_.*", "", clusto$Var2)
clusto$method = gsub(".*_", "", clusto$Var2)
clusto$omic = gsub("expression", "Expression", clusto$omic)
clusto$omic = gsub("mirna", "miRNA", clusto$omic)
clusto$omic = gsub("methylation", "Methylation", clusto$omic)

color_vals=c("#F8766D","#00BA38","#619CFF", "#774936", "black")
names(color_vals) = c("Expression", "Methylation", "miRNA", "COCA StoM", "ClustOmics StoM")
symb_vals = c(21, 21, 21, 23, 22)
names(symb_vals) = names(color_vals)

p <- ggplot(clusto, aes(x=cancer, y=rands, fill=omic)) + geom_boxplot(fill="white") + 
  geom_point(data=clusto, aes(shape=omic), size=3, alpha=0.8, position=position_dodge(width=0.5)) +
  scale_fill_manual(name = "Omic", values = color_vals) +
  scale_colour_manual(name = "Omic",values = color_vals) +   
  scale_shape_manual(name = "Omic", values = symb_vals) +  theme_bw() +
  labs(x ="Cancer", y = "Adjusted Rand Index", shape="Omic", fill = "Omic")
p

# svg(file = "./out/plots/overall_randIndex_SingleToMulti.svg")
# p
# dev.off()

## COCA
coca = all_rands[which((all_rands$Var1 == "COCA StoM") & (all_rands$Var2 != "COCA StoM")),]
coca$omic = gsub("_.*", "", coca$Var2)
coca$method = gsub(".*_", "", coca$Var2)
coca$omic = gsub("expression", "Expression", coca$omic)
coca$omic = gsub("mirna", "miRNA", coca$omic)
coca$omic = gsub("methylation", "Methylation", coca$omic)

color_vals=c("#F8766D","#00BA38","#619CFF", "#774936", "black")
names(color_vals) = c("Expression", "Methylation", "miRNA", "COCA StoM", "ClustOmics StoM")
symb_vals = c(21, 21, 21, 23, 22)
names(symb_vals) = names(color_vals)

q <- ggplot(coca, aes(x=cancer, y=rands, fill=omic)) + geom_boxplot(fill="white") + 
  geom_point(data=coca, aes(shape=omic), size=3, alpha=0.8, position=position_dodge(width=0.5)) +
  scale_fill_manual(name = "Omic", values = color_vals) +
  scale_colour_manual(name = "Omic",values = color_vals) +   
  scale_shape_manual(name = "Omic", values = symb_vals) +  theme_bw() +
  labs(x ="Cancer", y = "Adjusted Rand Index", shape="Omic", fill = "Omic") 
q

## Both
color_vals=c("#F8766D","#00BA38","#619CFF", "#774936", "black")
names(color_vals) = c("Expression", "Methylation", "miRNA", "COCA StoM", "ClustOmics StoM")
color_box = c("#F8766D", "#F8766D")
names(color_box) = c("COCA StoM", "ClustOmics StoM")
colorss = list(omic = color_vals, Var1=color_box)
symb_vals = c(21, 21, 21, 23, 22)
names(symb_vals) = names(color_vals)

cocaclusto = rbind(coca, clusto)

r <- ggplot(cocaclusto, aes(x=Var1, y=rands)) + 
  geom_boxplot(color="black", fill="grey", show.legend = F, alpha=0.5) + 
  theme(axis.text.y = element_text(angle = 90), axis.text.x = element_text(angle = 90)) +
  geom_point(data=cocaclusto, aes(shape=omic, fill=omic), size=3, alpha=0.8, position=position_dodge(width=0.5)) +
  facet_wrap(~ cancer, nrow=10, switch = "y") + coord_flip() + scale_x_discrete(position = "top") +
  scale_fill_manual(name = "Omic", values = color_vals) +
  scale_colour_manual(name = "Omic",values = color_vals) +   
  scale_shape_manual(name = "Omic", values = symb_vals) +  theme_bw() +
  labs(x="", y = "Adjusted Rand Index", shape="Omic", fill = "Omic") + theme(legend.position="bottom")
r

svg(file = "./outAll/plots/ClustOmics_COCA_randIndex_SingleToMulti.svg", width = 10, height = 10)
r
dev.off()


# MULTI TO MULTI
methods = c("PINS", "NEMO", "SNF", "rMKL", "MCCA")
datatypes = c("multiomics")

clusterings = c()
for (datatype in datatypes) {
  for (method in methods) {
    if ((method == "MCCA") && (datatype != "multiomics")) {
      next
    } else {
      run = paste(datatype, method, sep="_")
      clusterings = c(clusterings, run)
    }
  }
}
clusterings = c(clusterings, "ClustOmics", "COCA")
combinations = expand.grid(clusterings, clusterings)

all_rands = data.frame()

for (cancer in cancers) {
  out_plot = paste("outOnlyMulti/plots/", cancer, "/", cancer, "_heatmap_rand_index_tree_MultiToMulti.png", sep="")
  
  combinations$rands = apply(combinations, 1, function(x) {
    if ((x[[1]] %in% c("mirna_MCCA", "methylation_MCCA", "expression_MCCA")) | (x[[2]] %in% c("mirna_MCCA", "methylation_MCCA", "expression_MCCA"))) {
      return(NA)
    } else {
      if (x[[1]] == "ClustOmics") {
        file1 = paste("outOnlyMulti/results/", cancer, "/", cancer, ".", cancer, "_MULTI_MCCA_NEMO_PINS_SNF_rMKL.ClustOmicsClustering", sep="")
      } else if (x[[1]] == "COCA") {
        file1 = paste("dataOnlyMulti/coca/", cancer, "_multiomics_COCA.clst", sep="")
      } else {
        file1 = paste("dataOnlyMulti/", cancer, "/", cancer, "_", x[[1]], ".clst", sep="")
      }
      
      if (x[[2]] == "ClustOmics") {
        file2 = paste("outOnlyMulti/results/", cancer, "/", cancer, ".", cancer, "_MULTI_MCCA_NEMO_PINS_SNF_rMKL.ClustOmicsClustering", sep="")
      } else if (x[[2]] == "COCA") {
        file2 = paste("dataOnlyMulti/coca/", cancer, "_multiomics_COCA.clst", sep="")
      } else {
        file2 = paste("dataOnlyMulti/", cancer, "/", cancer, "_", x[[2]], ".clst", sep="")
      }

      c1 = read.table(file1, header = T)
      c2 = read.table(file2, header = T)
      
      cl1 = c1$Cluster
      names(cl1) = c1$Patient
      cl2 = c2$Cluster
      names(cl2) = c2$Patient
      
      multi = intersect(names(cl1), names(cl2))
      cl1 = cl1[multi]
      cl2 = cl2[multi]
      
      return(adj.rand.index(cl1, cl2))
    }
  })
  
  combinations = combinations[which(!is.na(combinations$rands)),]
  df = matrix(ncol = length(clusterings), nrow = length(clusterings), data = combinations$rands)
  colnames(df) = unique(combinations$Var2)
  rownames(df) = unique(combinations$Var1)
  
  rownames(df) = gsub("ClustOmics", "multiomics_ClustOmics", rownames(df))
  rownames(df) = gsub("COCA", "multiomics_COCA", rownames(df))
  colnames(df) = row.names(df)
  
  my_group = as.factor(gsub("_.*", "", rownames(df)))
  my_group = gsub("ClustOmics", "Multiomics", my_group)
  names(my_group) = rownames(df)
  my_group = as.data.frame(my_group)
  names(my_group) = "Omic"
  my_group$Omic = gsub("expression", "Expression", my_group$Omic)
  my_group$Omic = gsub("mirna", "miRNA", my_group$Omic)
  my_group$Omic = gsub("methylation", "Methylation", my_group$Omic)
  my_group$Omic = gsub("multiomics", "Multiomics", my_group$Omic)
  my_colour = list(Omic = c(Expression = "#F8766D", Methylation="#00BA38", miRNA="#619CFF", Multiomics="#C77CFF"))
  p = pheatmap(df, color=magma(10))

  save_pheatmap_png(p, out_plot)

  combinations$cancer = cancer
  all_rands = rbind(all_rands, combinations)
}

clusto = all_rands[which((all_rands$Var1 == "ClustOmics") & (all_rands$Var2 != "ClustOmics")),]
clusto$omic = gsub("_.*", "", clusto$Var2)
clusto$method = gsub(".*_", "", clusto$Var2)
clusto$omic = gsub("expression", "Expression", clusto$omic)
clusto$omic = gsub("mirna", "miRNA", clusto$omic)
clusto$omic = gsub("methylation", "Methylation", clusto$omic)

vals=c("#F8766D","#00BA38","#619CFF", "#E6AB02", "#A9A9A9","black", "#A6761D")
names(vals) = c("PINS", "SNF", "rMKL", "MCCA", "NEMO", "ClustOmics", "COCA")

p <- ggplot(clusto, aes(x=cancer, y=rands)) + geom_boxplot(fill="white") +
  geom_point(data=clusto, aes(shape=method, fill=method, color=method), size=3, position=position_dodge(width=0.5)) +
  scale_fill_manual(name = "Input Clustering",values = vals) +
  scale_colour_manual(name = "Input Clustering",values = vals) +   
  scale_shape_manual(name = "Input Clustering",values = c(PINS=21, SNF=21, NEMO=21, MCCA=21, rMKL=21, ClustOmics=24, COCA=25)) +
  labs(x ="Cancer", y = "Adjusted Rand Index") +
  theme_bw() 
p

svg(file = "./outOnlyMulti/plots/overall_randIndex_MultiToMulti.svg")
p
dev.off()
