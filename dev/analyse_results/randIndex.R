library(pdfCluster)
library(RColorBrewer)
library(gplots)
library(pheatmap)
library(ggplot2)
library(viridis)

setwd("~/ClustOmics")

save_pheatmap_png <- function(x, filename, width=1200, height=1000, res = 150) {
  png(filename, width = width, height = height, res = res)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

cancers = c("AML", "BIC", "COAD", "GBM", "KIRC", "LIHC", "LUSC", "OV", "SARC", "SKCM")

# SINGLE TO MULTI
methods = c("PINS", "NEMO", "SNF", "rMKL", "MCCA")
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
clusterings = c(clusterings, "ClustOmics")
combinations = expand.grid(clusterings, clusterings)

all_rands = data.frame()

for (cancer in cancers) {
  out_plot = paste("out/plots/", cancer, "/", cancer, "_heatmap_rand_index_tree_SingleToMulti.png", sep="")

  combinations$rands = apply(combinations, 1, function(x) {
    if (x[[1]] == "ClustOmics") {
      file1 = paste("out/results/", cancer, "/", cancer, ".", cancer, "_EXP_MIRNA_MET_NEMO_PINS_SNF_rMKL.ClustOmicsClustering", sep="")
    } else {
      file1 = paste("data/", cancer, "/", cancer, "_", x[[1]], ".clst", sep="")
    }
    
    if (x[[2]] == "ClustOmics") {
      file2 = paste("out/results/", cancer, "/", cancer, ".", cancer, "_EXP_MIRNA_MET_NEMO_PINS_SNF_rMKL.ClustOmicsClustering", sep="")
    } else {
      file2 = paste("data/", cancer, "/", cancer, "_", x[[2]], ".clst", sep="")
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
  
  rownames(df) = gsub("ClustOmics", "multiomics_ClustOmics", rownames(df))
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
  p = pheatmap(df, annotation_row = my_group, annotation_colors = my_colour, color=magma(5))

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

p <- ggplot(clusto, aes(x=cancer, y=rands, fill=omic)) + geom_boxplot(fill="white") + 
  geom_dotplot(binaxis='y', stackdir='center', dotsize = 0.5, position=position_dodge(width=0.2)) + theme_bw() +
  labs(x ="Cancer", y = "Adjusted Rand Index", fill = "Omic") #+ scale_fill_brewer(palette = "Paired")

svg(file = "./out/plots/overall_randIndex_SingleToMulti.svg")
p
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
clusterings = c(clusterings, "ClustOmics")
combinations = expand.grid(clusterings, clusterings)

all_rands = data.frame()

for (cancer in cancers) {
  out_plot = paste("out/plots/", cancer, "/", cancer, "_heatmap_rand_index_tree_MultiToMulti.png", sep="")
  
  combinations$rands = apply(combinations, 1, function(x) {
    if ((x[[1]] %in% c("mirna_MCCA", "methylation_MCCA", "expression_MCCA")) | (x[[2]] %in% c("mirna_MCCA", "methylation_MCCA", "expression_MCCA"))) {
      return(NA)
    } else {
      if (x[[1]] == "ClustOmics") {
        file1 = paste("out/results/", cancer, "/", cancer, ".", cancer, "_MULTI_MCCA_NEMO_PINS_SNF_rMKL.ClustOmicsClustering", sep="")
      } else {
        file1 = paste("data/", cancer, "/", cancer, "_", x[[1]], ".clst", sep="")
      }
      
      if (x[[2]] == "ClustOmics") {
        file2 = paste("out/results/", cancer, "/", cancer, ".", cancer, "_MULTI_MCCA_NEMO_PINS_SNF_rMKL.ClustOmicsClustering", sep="")
      } else {
        file2 = paste("data/", cancer, "/", cancer, "_", x[[2]], ".clst", sep="")
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
  p = pheatmap(df, color=magma(5))
  
  save_pheatmap_png(p, out_plot)
  
  # png(file = out_plot)
  # heatmap.2(df, margins=c(10,10), #scale="column", #Colv = NA, Rowv = NA,
  #           RowSideColors=colSide, ColSideColors=colSide, revC=F,
  #           trace = "none", density.info = "none")
  # dev.off()
  
  combinations$cancer = cancer
  all_rands = rbind(all_rands, combinations)
}

clusto = all_rands[which((all_rands$Var1 == "ClustOmics") & (all_rands$Var2 != "ClustOmics")),]
clusto$omic = gsub("_.*", "", clusto$Var2)
clusto$method = gsub(".*_", "", clusto$Var2)
clusto$omic = gsub("expression", "Expression", clusto$omic)
clusto$omic = gsub("mirna", "miRNA", clusto$omic)
clusto$omic = gsub("methylation", "Methylation", clusto$omic)


p <- ggplot(clusto, aes(x=cancer, y=rands, fill=method)) + geom_boxplot(fill="white") +
  geom_dotplot(binaxis='y', stackdir='center', dotsize = 0.5) + theme_bw() +
  labs(x ="Cancer", y = "Adjusted Rand Index", fill = "Input Clustering")

svg(file = "./out/plots/overall_randIndex_MultiToMulti.svg")
p
dev.off()
