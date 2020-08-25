#!/usr/bin/env Rscript

library(parallel)

check_enrichment <- function(clust_path, clin_path, cancer) {
  clinical_metadata = list(gender='DISCRETE', age_at_initial_pathologic_diagnosis='NUMERIC', pathologic_M='DISCRETE', 
                           pathologic_N='DISCRETE', pathologic_T='DISCRETE', pathologic_stage='DISCRETE', histological_type='DISCRETE',
                           new_neoplasm_event_type='DISCRETE', neoplasm_histologic_grade='DISCRETE')
  
  if (cancer == "BIC") {
    sup_clin = list(PAM50Call_RNAseq='DISCRETE', breast_carcinoma_estrogen_receptor_status='DISCRETE', 
                    breast_carcinoma_progesterone_receptor_status="DISCRETE", er_level_cell_percentage_category='DISCRETE',
                    progesterone_receptor_level_cell_percent_category='DISCRETE')
  } else if (cancer == "LUSC") {
    sup_clin = list(tobacco_smoking_history='DISCRETE', number_pack_years_smoked='NUMERIC')
  } else if (cancer == "SKCM") {
    sup_clin = list(melanoma_clark_level_value='DISCRETE', melanoma_ulceration_indicator='DISCRETE')
  } else if (cancer == "AML") {
    sup_clin = list(acute_myeloid_leukemia_calgb_cytogenetics_risk_category='DISCRETE', leukemia_french_american_british_morphology_code='DISCRETE')
  } else if (cancer == "COAD") {
    sup_clin = list(colon_polyps_present='DISCRETE', history_of_colon_polyps='DISCRETE')
  } else if (cancer == "GBM") {
    sup_clin = list(prior_glioma='DISCRETE')
  } else if (cancer == "KIRC") {
    sup_clin = list(hemoglobin_result='DISCRETE', platelet_qualitative_result='DISCRETE', serum_calcium_result='DISCRETE', white_cell_count_result='DISCRETE')
  } else if (cancer == "LIHC") {
    sup_clin = list(adjacent_hepatic_tissue_inflammation_extent_type='DISCRETE', albumin_result_specified_value='NUMERIC', creatinine_value_in_mg_dl='NUMERIC',
                    fetoprotein_outcome_value='NUMERIC', fibrosis_ishak_score='DISCRETE')
  } else if (cancer == "OV" | cancer == "SARC") {
    sup_clin = list()
  } 
  
  clinical_metadata = append(clinical_metadata, sup_clin)
  

  clust <- read.table(clust_path, header=T, sep="\t")

  clin <- read.table(clin_path, sep="\t", header=T, row.names = 1, stringsAsFactors = F, na.strings = c("NA",""))
  clin <- clin[which(names(clin) %in% names(clinical_metadata))] 
  clin <- clin[which(!duplicated(substring(rownames(clin), 1, 15))),]
  row.names(clin) <- gsub('-', '\\.', toupper(substring(rownames(clin), 1, 15)))
  clin <- clin[which(row.names(clin) %in% clust$Patient),]
  
  stopifnot(nrow(clin) == nrow(clust))
  
  indices <- match(clust$Patient, row.names(clin))
  clin <- clin[indices,]
  clin$cluster <- clust$Cluster
  
  pvalues = c()
  params_tested = c()
  
  for (clinical_param in names(clinical_metadata)) {
    
    print(clinical_param)
    
    if (!(clinical_param %in% colnames(clin))) {
      next
    }
    
    is_discrete_param <- clinical_metadata[clinical_param] == 'DISCRETE' #boolean
    is_numeric_param <- clinical_metadata[clinical_param] == 'NUMERIC'
    stopifnot(is_discrete_param | is_numeric_param)
    
    df = clin[,c(clinical_param, "cluster")]
    
    if (clinical_param == "pahtologic_M") {
      #df = clin[,c("pathologic_M", "cluster")] 
      df$pathologic_M = gsub("a", "", df$pathologic_M)
      df$pathologic_M = gsub("b", "", df$pathologic_M)
      df$pathologic_M = gsub("c", "", df$pathologic_M)
      df$pathologic_M = gsub(" \\(i\\+\\)", "", df$pathologic_M)
      df[which(!df$pathologic_M %in% c("MX", "M0", "M1")),1] = NA
    } else if (clinical_param == "pahtologic_T") {
      #df = clin[,c("pathologic_T", "cluster")]
      df[which(df$pathologic_T == "[Discrepancy]"),1] = NA
      df$pathologic_T = gsub("a", "", df$pathologic_T)
      df$pathologic_T = gsub("b", "", df$pathologic_T)
      df$pathologic_T = gsub("c", "", df$pathologic_T)
      df$pathologic_T = gsub("d", "", df$pathologic_T)
      df[which(!df$pathologic_T %in% c("TX", "T0", "T1", "T2", "T3", "T4", "Tis")),1] = NA
    } else if (clinical_param == "pathologic_N") {
      #df = clin[,c("pathologic_N", "cluster")] 
      df$pathologic_N = gsub("a", "", df$pathologic_N)
      df$pathologic_N = gsub("b", "", df$pathologic_N)
      df$pathologic_N = gsub("c", "", df$pathologic_N)
      df$pathologic_N = gsub("mi", "", df$pathologic_N)
      df$pathologic_N = gsub(" \\(i\\+\\)", "", df$pathologic_N)
      df$pathologic_N = gsub(" \\(i\\-\\)", "", df$pathologic_N)
      df$pathologic_N = gsub(" \\(mol\\+\\)", "", df$pathologic_N)
      df[which(!df$pathologic_N %in% c("NX", "N1", "N2", "N3", "N0")),1] = NA
    } else if (clinical_param == "pathologic_stage") {
      #df = clin[,c("pathologic_stage", "cluster")]
      df[which(df$pathologic_stage == "[Discrepancy]"),1] = NA
      df$pathologic_stage = gsub("A", "", df$pathologic_stage)
      df$pathologic_stage = gsub("B", "", df$pathologic_stage)
      df$pathologic_stage = gsub("C", "", df$pathologic_stage)
      df[which(!df$pathologic_stage %in% c("Stage X", "Stage 0", "Stage I", "Stage II", "Stage III", "Stage IV")),1] = NA
    } else if (clinical_param == "histological_type") { ## TROP VARIABLE
      #df = clin[,c("histological_type", "cluster")]
      df[which(df$histological_type == "[Discrepancy]"),1] = NA
      df[which(df$histological_type == "Other, specify"),1] = NA
      df[which(df$histological_type == "Mixed Histology (please specify)"),1] = NA
      # if (length(unique(df$histological_type)) > 10) {
      #   next
      # }
    } else if (clinical_param == "new_neoplasm_event_type") {
      #df = clin[,c("new_neoplasm_event_type", "cluster")]
      df[which(!is.na(df$new_neoplasm_event_type)),1] = "Yes"
      df[which(is.na(df$new_neoplasm_event_type)),1] = "None"
    } else if (clinical_param == "neoplasm_histologic_grade") {
      #df = clin[,c("neoplasm_histologic_grade", "cluster")]
      df[which(!df$neoplasm_histologic_grade %in% c("GX", "G1", "G2", "G3", "G4")),1] = NA
    } else if (clinical_param == "tobacco_smoking_history") {
      df[which(df$tobacco_smoking_history == "[Discrepancy]"),1] = NA
    } else if (clinical_param == "leukemia_french_american_british_morphology_code") {
      df[which(df$leukemia_french_american_british_morphology_code == "Not Classified"),1] = NA
    } 
    
    clinical_values <- df[,clinical_param]
    
    if (is_numeric_param) {
      numeric_entries = !is.na(as.numeric(clinical_values)) # Boolean list > TRUE if not NA and FALSE if is NA
      # if too many missing data (more than half), skip the param
      if (2 * sum(numeric_entries) < length(clinical_values)) {
        next
      }
    } else {
      not_na_entries = !is.na(clinical_values) 
      should_skip = F
      if (2 * sum(not_na_entries) < length(clinical_values)) {
        should_skip = T
      } else if (length(table(clinical_values[not_na_entries])) == 1) {
        should_skip = T
      }
      if (should_skip) {
        next
      }
    }
    
    params_tested = c(params_tested, clinical_param)
    
    clustering_list <- df$cluster
    names(clustering_list) <- row.names(df)
    
    if (is_discrete_param) {
      pvalue = get.empirical.clinical(clustering_list[!is.na(clinical_values)], clinical_values[!is.na(clinical_values)], T)
    } else if (is_numeric_param) {
      pvalue = get.empirical.clinical(clustering_list[numeric_entries], as.numeric(clinical_values[numeric_entries]), F)
    }
    
    pvalues = c(pvalues, pvalue)
    
  }
  
  names(pvalues) = params_tested
  
  sum(pvalues < 0.05)
  
  return(pvalues)
  
}

get.empirical.clinical <- function(clustering, clinical.values, is.chisq) {
  set.seed(42)
  if (is.chisq) {
    clustering.with.clinical = cbind(clustering, clinical.values)
    tbl = table(as.data.frame(clustering.with.clinical))
    test.res = chisq.test(tbl)
  } else {
    test.res = kruskal.test(as.numeric(clinical.values), clustering)
  }
  orig.pvalue = test.res$p.value
  num.iter = 1000
  total.num.iters = 0
  total.num.extreme = 0
  should.continue = T
  
  while (should.continue) {
    perm.pvalues = as.numeric(mclapply(1:num.iter, function(i) {
      cur.clustering = sample(clustering)
      names(cur.clustering) = names(clustering)
      
      if (is.chisq) {
        clustering.with.clinical = cbind(cur.clustering, clinical.values)
        tbl = table(as.data.frame(clustering.with.clinical))
        test.res = chisq.test(tbl)
      } else {
        test.res = kruskal.test(as.numeric(clinical.values), cur.clustering)
      }
      cur.pvalue = test.res$p.value
      return(cur.pvalue)
    }, mc.cores=5))
    total.num.iters = total.num.iters + num.iter
    total.num.extreme = total.num.extreme + sum(perm.pvalues <= orig.pvalue)
    
    binom.ret = binom.test(total.num.extreme, total.num.iters)
    cur.pvalue = binom.ret$estimate
    cur.conf.int = binom.ret$conf.int
    
    sig.threshold = 0.05
    is.threshold.in.conf = cur.conf.int[1] < sig.threshold & cur.conf.int[2] < sig.threshold
    if (!is.threshold.in.conf | total.num.iters > 1e5) {
      should.continue = F
    }
  }
  
  return(cur.pvalue)
}

library("optparse")

option_list = list(
  make_option(c("-c", "--clust"), type="character", default=NULL, 
              help="Clustering file name", metavar="character"),
  make_option(c("-f", "--clinical_file"), type="character", default=NULL, 
              help="Survival data", metavar="character"),
  make_option(c("-o", "--out_enrich"), type="character", default=NULL, 
              help="output file name (clinical labels enrichment pvalues)", metavar="character"),
  make_option(c("-s", "--cancer"), type="character", default=NULL, 
              help="Cancer", metavar="character")
) 

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

pvalues <- check_enrichment(opt$clust, opt$clinical_file, cancer = opt$cancer)
pvalues <- as.data.frame(pvalues)
write.table(pvalues, opt$out_enrich, col.names=F, row.names=T, sep="\t", quote=F)


       
