library(chisq.posthoc.test)
library(FSA)
library(ggpubr)

setwd("describeBIC")

cancer = "BIC"
clustering = "../outAll/results/BIC/BIC.BIC_EXP_MIRNA_MET_NEMO_PINS_SNF_rMKL_kmeans_all.ClustOmicsClustering"
clinical = "../raw_data/clinical/BIC"
myPAM = "myclassifPAM.csv"

clinical_metadata = list(gender='DISCRETE', age_at_initial_pathologic_diagnosis='NUMERIC', pathologic_M='DISCRETE', 
                         pathologic_N='DISCRETE', pathologic_T='DISCRETE', pathologic_stage='DISCRETE', histological_type='DISCRETE',
                         new_neoplasm_event_type='DISCRETE', neoplasm_histologic_grade='DISCRETE')
sup_clin = list(PAM50Call_RNAseq='DISCRETE', breast_carcinoma_estrogen_receptor_status='DISCRETE', 
                breast_carcinoma_progesterone_receptor_status="DISCRETE", er_level_cell_percentage_category='DISCRETE',
                progesterone_receptor_level_cell_percent_category='DISCRETE')
clinical_metadata = append(clinical_metadata, sup_clin)


clust <- read.table(clustering, header=T, sep="\t")

clin = read.table(clinical, sep="\t", header=T, row.names = 1, stringsAsFactors = F, na.strings = c("NA",""))
clin = clin[which(names(clin) %in% names(clinical_metadata))] 
clin = clin[which(!duplicated(substring(rownames(clin), 1, 15))),]
row.names(clin) = gsub('-', '\\.', toupper(substring(rownames(clin), 1, 15)))
clin = clin[which(row.names(clin) %in% clust$Patient),]
clin$Cluster = lapply(row.names(clin), function(x) {
  return(clust[which(clust$Patient == x),2])
})
clin$Cluster = unlist(clin$Cluster) 

myPam = read.table(myPAM, header=T, sep=";")
myPam$pam = gsub("My", "", myPam$pam)
clin$PAM50Call_RNAseq = lapply(row.names(clin), function(x) {
  return(myPam[which(myPam$Patient == x),2])
})
clin$PAM50Call_RNAseq = unlist(clin$PAM50Call_RNAseq) 

clin$pathologic_M = gsub("c", "", clin$pathologic_M)
clin$pathologic_M = gsub(" \\(i\\+\\)", "", clin$pathologic_M)

clin$pathologic_N = gsub("a", "", clin$pathologic_N)
clin$pathologic_N = gsub("b", "", clin$pathologic_N)
clin$pathologic_N = gsub("c", "", clin$pathologic_N)
clin$pathologic_N = gsub("mi", "", clin$pathologic_N)
clin$pathologic_N = gsub(" \\(i\\+\\)", "", clin$pathologic_N)
clin$pathologic_N = gsub(" \\(i\\-\\)", "", clin$pathologic_N)
clin$pathologic_N = gsub(" \\(mol\\+\\)", "", clin$pathologic_N)

clin$pathologic_T = gsub("\\[Discrepancy\\]", NA, clin$pathologic_T)
clin$pathologic_T = gsub("a", "", clin$pathologic_T)
clin$pathologic_T = gsub("b", "", clin$pathologic_T)
clin$pathologic_T = gsub("c", "", clin$pathologic_T)
clin$pathologic_T = gsub("d", "", clin$pathologic_T)

clin$pathologic_stage = gsub("\\[Discrepancy\\]", NA, clin$pathologic_stage)
clin$pathologic_stage = gsub("A", "", clin$pathologic_stage)
clin$pathologic_stage = gsub("B", "", clin$pathologic_stage)
clin$pathologic_stage = gsub("C", "", clin$pathologic_stage)

for (clinical_param in intersect(names(clinical_metadata), names(clin))) {
  is_discrete = clinical_metadata[clinical_param] == 'DISCRETE'
  print(clinical_param)
  if (is_discrete) {
    tab = table(clin$Cluster, clin[,clinical_param])
    print(chisq.posthoc.test(tab, simulate.p.value = TRUE))
  }
}

dunn_age_diag = dunnTest(age_at_initial_pathologic_diagnosis ~ Cluster,
                         data=clin,
                         method="bh")
print(dunn_age_diag)
a=dunn_age_diag$res
a[which(a$P.adj<=0.001),]

my_comparisons <- list( c("D", "C"), c("D", "E"), c("D", "F"), c("A", "C"), c("A", "E"),
                        c("A", "F"), c("B", "C"), c("B", "E"), c("B", "F"))

clin$Cluster = gsub("1", "A", clin$Cluster)
clin$Cluster = gsub("2", "B", clin$Cluster)
clin$Cluster = gsub("3", "C", clin$Cluster)
clin$Cluster = gsub("4", "D", clin$Cluster)
clin$Cluster = gsub("6", "E", clin$Cluster)
clin$Cluster = gsub("8", "F", clin$Cluster)
clin$Cluster = factor(clin$Cluster, level = c('A', 'B', 'C', 'D', 'E', 'F'))
plot = ggboxplot(clin, x = "Cluster", y = "age_at_initial_pathologic_diagnosis", color='Cluster')+
  stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = 5) + xlab("Consensus Cluster") +
  ylab("Age at diagnosis")

pdf(file="ageDiagPlot.pdf")
plot
dev.off()

