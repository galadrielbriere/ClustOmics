library(survival)
library(optparse)
library(parallel)
library(survminer)

check.survival <- function(groups, survival.file.path) {
  survival.data = read.table(survival.file.path, header = TRUE)
  patient.names = names(groups)
  patient.names.in.file = as.character(survival.data[, 1])
  patient.names.in.file = toupper(gsub('-', '\\.', substring(patient.names.in.file, 1, 15)))
  
  if (!all(patient.names %in% patient.names.in.file)) {
    patient.names <- toupper(gsub('-', '\\.', substring(patient.names, 1, 12)))
    stopifnot(all(patient.names %in% patient.names.in.file))
  }

  indices = match(patient.names, patient.names.in.file)
  ordered.survival.data = survival.data[indices,]
  ordered.survival.data["cluster"] <- groups
  #ordered.survival.data$Survival[is.na(ordered.survival.data$Survival)] = 0
  #ordered.survival.data$Death[is.na(ordered.survival.data$Death)] = 0
  return(survdiff(Surv(Survival, Death) ~ cluster, data=ordered.survival.data))
}

plot.survival <- function(groups, survival.file.path) {
  survival.data = read.table(survival.file.path, header = TRUE)
  patient.names = names(groups)
  patient.names.in.file = as.character(survival.data[, 1])
  patient.names.in.file = toupper(gsub('-', '\\.', substring(patient.names.in.file, 1, 15)))
  
  if (!all(patient.names %in% patient.names.in.file)) {
    patient.names <- toupper(gsub('-', '\\.', substring(patient.names, 1, 12)))
    stopifnot(all(patient.names %in% patient.names.in.file))
  }
  
  indices = match(patient.names, patient.names.in.file)
  ordered.survival.data = survival.data[indices,]
  ordered.survival.data["cluster"] <- groups
  #ordered.survival.data$Survival[is.na(ordered.survival.data$Survival)] = 0
  #ordered.survival.data$Death[is.na(ordered.survival.data$Death)] = 0
  res = survfit(Surv(Survival, Death) ~ cluster, data=ordered.survival.data)
  return(list(data = ordered.survival.data, res = res))
}

get.logrank.pvalue <- function(survdiff.res) {
  1 - pchisq(survdiff.res$chisq, length(survdiff.res$n) - 1)  
}

get.empirical.surv <- function(clustering, surv) {
  set.seed(42)
  surv.ret = check.survival(clustering, surv)
  orig.chisq = surv.ret$chisq
  orig.pvalue = get.logrank.pvalue(surv.ret)
  # The initial number of permutations to run
  num.perms = round(min(max(10 / orig.pvalue, 1000), 1e6))
  should.continue = T
  
  total.num.perms = 0
  total.num.extreme.chisq = 0
  
  while (should.continue) {
    print(num.perms)
    perm.chisq = as.numeric(mclapply(1:num.perms, function(i) {
      cur.clustering = sample(clustering)
      names(cur.clustering) = names(clustering)
      cur.chisq = check.survival(cur.clustering, surv)$chisq
      return(cur.chisq)
    }, mc.cores=4))
    
    total.num.perms = total.num.perms + num.perms
    total.num.extreme.chisq = total.num.extreme.chisq + sum(perm.chisq >= orig.chisq)
    
    binom.ret = binom.test(total.num.extreme.chisq, total.num.perms)
    cur.pvalue = binom.ret$estimate
    cur.conf.int = binom.ret$conf.int
    
    print(c(total.num.extreme.chisq, total.num.perms))
    print(cur.pvalue)
    print(cur.conf.int)
    
    sig.threshold = 0.05
    is.conf.small = ((cur.conf.int[2] - cur.pvalue) < min(cur.pvalue / 10, 0.01)) & ((cur.pvalue - cur.conf.int[1]) < min(cur.pvalue / 10, 0.01))
    is.threshold.in.conf = cur.conf.int[1] < sig.threshold & cur.conf.int[2] > sig.threshold
    if ((is.conf.small & !is.threshold.in.conf) | (total.num.perms > 1e5)) {
      #if (is.conf.small) {
      should.continue = F
    } else {
      num.perms = 1e5
    }
  }
  
  return(list(pvalue = cur.pvalue, conf.int = cur.conf.int, total.num.perms=total.num.perms, 
              total.num.extreme.chisq=total.num.extreme.chisq))
}


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
names(clustering) <- clust$Patient

result_surv <- get.empirical.surv(clustering, opt$surv_file)
pval <- result_surv$pvalue

write.table(pval, opt$out_survival, col.names=F, row.names=F, sep="\t", quote=F)

res = plot.survival(clustering, opt$surv_file)

svg(file = opt$out_fig)
plot = ggsurvplot(res$res, data = res$data, risk.table = F, pval = TRUE)
print(plot)
dev.off()

