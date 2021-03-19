configfile: "config.yaml"

NEO4J_ID = config['neo4j_id']
NEO4J_PASSWORD =  config['neo4j_password']
DATABASE = config['database']
CLUSTERINGS_FOLDER = config['clusterings_folder']
OBJECT_NODES_NAME = config['object_nodes_name'],
CLUSTER_NODES_NAME = config['cluster_nodes_name']
WRITE_ALL = config['write_all']
REASSIGN = config['reassign_unclassified']

# Run ClustOmics
rule All:
    input:
        clust="out/{subject}.{rel_name}.FuseClusterings.log",
        MQ_plot="out/plots/{subject}/wMQ_{subject}.{rel_name}.svg",
        survival="out/survival/{subject}/{subject}.{rel_name}.ClustOmicsClustering.pval",
        clinical="out/clinical/{subject}/{subject}.{rel_name}.ClustOmicsClustering.pval"
    output:
        "out/{subject}.{rel_name}.all.log"
    shell:
        "touch {output}"

# Build the graph from metadata file and raw clustering results
rule BuildGraph:
    input:
        metadata = "data/{subject}/{subject}_metadata.txt"
    output:
        "out/{subject}.BuildGraph.log"
    shell:
        "python dev/BuildGraph.py -id {NEO4J_ID} -pwd {NEO4J_PASSWORD} -host {DATABASE}\
        -cls {CLUSTERINGS_FOLDER} -subj {wildcards.subject} -out {output}"

# Create edges between each pair of objects clustered together in at least one raw clustering results
rule SupportEdges:
    input:
        "out/{subject}.BuildGraph.log"
    output:
        "out/{subject}.SupportEdges.log"
    shell:
        "python dev/SupportEdges.py -id {NEO4J_ID} -pwd {NEO4J_PASSWORD} -host {DATABASE}\
         -subj {wildcards.subject} -obj_name {OBJECT_NODES_NAME} -cls_name {CLUSTER_NODES_NAME} -out {output}"

# Create integration edges according to the specified methods and omics to integrate
rule IntegrationEdges:
    input:
        "out/{subject}.SupportEdges.log"
    output:
        "out/{subject}.{rel_name}.IntegrationEdges.log"
    params:
         datatypes=lambda wildcards: config['datatypes'][wildcards.rel_name],
         methods=lambda wildcards: config['methods'][wildcards.rel_name]
    shell:
        "python dev/IntegrationEdges.py -id {NEO4J_ID} -pwd {NEO4J_PASSWORD} -host {DATABASE} \
        -subj {wildcards.subject} -obj_name {OBJECT_NODES_NAME} -rel {wildcards.rel_name} \
         -dt '{params.datatypes}' -met '{params.methods}' -op 'OR' -out {output}"

rule FuseClusterings:
    input:
        "out/{subject}.{rel_name}.IntegrationEdges.log"
    output:
        out_log="out/{subject}.{rel_name}.FuseClusterings.log",
        out_clust="out/results/{subject}/{subject}.{rel_name}.ClustOmicsClustering"
    params:
         datatypes=lambda wildcards: config['datatypes'][wildcards.rel_name],
         methods=lambda wildcards: config['methods'][wildcards.rel_name],
         min_size_consensus=lambda wildcards: config['min_size_consensus'][wildcards.rel_name],
         min_size_clust=lambda wildcards: config['min_size_clust'][wildcards.rel_name]
    shell:
        "mkdir -p out/results/{wildcards.subject};"
        "python dev/FuseClusterings.py -out_cls {output.out_clust} -id {NEO4J_ID} -pwd {NEO4J_PASSWORD} -host {DATABASE} \
        -subj {wildcards.subject} -obj_name {OBJECT_NODES_NAME} -cls_name {CLUSTER_NODES_NAME} -rel {wildcards.rel_name} \
        -dt '{params.datatypes}' -met '{params.methods}' -min_clust {params.min_size_clust} \
        -min_nodes {params.min_size_consensus} -out {output.out_log} -nb_sup False --writeAll {WRITE_ALL} --reassign {REASSIGN}"

rule FuseClusteringsWithNbSupp:
    input:
        "out/{subject}.{rel_name}.IntegrationEdges.log"
    output:
        out_log="out/{subject}.{rel_name}.FuseClusterings.{nb_supports}_supports.log",
        out_clust="out/results/{subject}/{subject}.{rel_name}.ClustOmicsClustering.{nb_supports}_supports"
    params:
         datatypes=lambda wildcards: config['datatypes'][wildcards.rel_name],
         methods=lambda wildcards: config['methods'][wildcards.rel_name],
         min_size_consensus=lambda wildcards: config['min_size_consensus'][wildcards.rel_name],
         min_size_clust=lambda wildcards: config['min_size_clust'][wildcards.rel_name]
    shell:
        "mkdir -p out/results/{wildcards.subject};"
        "python dev/FuseClusterings.py -out_cls {output.out_clust} -id {NEO4J_ID} -pwd {NEO4J_PASSWORD} -host {DATABASE} \
        -subj {wildcards.subject} -obj_name {OBJECT_NODES_NAME} -cls_name {CLUSTER_NODES_NAME} -rel {wildcards.rel_name} \
        -dt '{params.datatypes}' -met '{params.methods}' -min_clust {params.min_size_clust} \
        -min_nodes {params.min_size_consensus} -out {output.out_log} -nb_sup {wildcards.nb_supports} --writeAll {WRITE_ALL} --reassign {REASSIGN}"

# Compute and plot the Modularization Quality for each possible number of supports
rule PlotMQ:
    input:
        "out/{subject}.{rel_name}.IntegrationEdges.log",
    output:
        "out/plots/{subject}/wMQ_{subject}.{rel_name}.svg"
    params:
         datatypes=lambda wildcards: config['datatypes'][wildcards.rel_name],
         methods=lambda wildcards: config['methods'][wildcards.rel_name],
         min_size_clust=lambda wildcards: config['min_size_clust'][wildcards.rel_name],
         min_size_consensus=lambda wildcards: config['min_size_consensus'][wildcards.rel_name]
    shell:
        "mkdir -p out/plots/{wildcards.subject};"
        "python dev/PlotMQ.py -out {output} -id {NEO4J_ID} -pwd {NEO4J_PASSWORD} -host {DATABASE} \
        -subj {wildcards.subject} -obj_name {OBJECT_NODES_NAME} -cls_name {CLUSTER_NODES_NAME} \
        -rel {wildcards.rel_name} -dt '{params.datatypes}' -met '{params.methods}' \
        -min_clust {params.min_size_clust} -min_nodes {params.min_size_consensus} "

# Compute Survival Analysis for ClustOmics clusterings
rule SurvivalClustOmics:
    input:
        surv="raw_data/{subject}/survival",
        clust="out/results/{subject}/{subject}.{rel_name}.ClustOmicsClustering"
    output:
        out="out/survival/{subject}/{subject}.{rel_name}.ClustOmicsClustering.pval",
        fig="out/survival/{subject}/{subject}.{rel_name}.ClustOmicsClustering.svg"
    shell:
        "mkdir -p out/survival/{wildcards.subject};"
        "Rscript dev/analyse_results/SurvivalHeinze.R -c {input.clust} -s {input.surv} -o {output.out} -f {output.fig}"

# Compute Survival Analysis for COCA clusterings
rule SurvivalCOCA:
    input:
        surv="raw_data/{subject}/survival",
        clust="data/coca/{subject}_{datatype}_COCA.clst"
    output:
        out="out/survival/coca/{subject}_{datatype}_COCA.pval",
        fig="out/survival/coca/{subject}_{datatype}_COCA.svg"
    shell:
        "mkdir -p out/survival/coca;"
        "Rscript dev/analyse_results/SurvivalHeinze.R -c {input.clust} -s {input.surv} -o {output.out} -f {output.fig}"

# Compute Survival Analysis for raw clusterings
rule SurvivalRaw:
    input:
        surv="raw_data/{subject}/survival",
        clust="data/{subject}/{subject}_{datatype}_{method}.clst"
    output:
        out="out/survival/{subject}/{subject}_{datatype}_{method}.pval",
        fig="out/survival/{subject}/{subject}_{datatype}_{method}.svg"
    shell:
        "mkdir -p out/survival/{wildcards.subject};"
        "Rscript dev/analyse_results/SurvivalHeinze.R -c {input.clust} -s {input.surv} -o {output.out} -f {output.fig}"

# Compute Clinical Labels Enrichment Analysis for ClustOmics clusterings
rule ClinicalClustOmics:
    input:
        clin="raw_data/clinical/{subject}",
        clust="out/results/{subject}/{subject}.{rel_name}.ClustOmicsClustering"
    output:
        "out/clinical/{subject}/{subject}.{rel_name}.ClustOmicsClustering.pval"
    shell:
        "mkdir -p out/clinical/{wildcards.subject};"
        "Rscript dev/analyse_results/ClinicalLabelsEnrichment.R -s {wildcards.subject} -c {input.clust} -f {input.clin} -o {output}"

# Compute Clinical Labels Enrichment Analysis for COCA clusterings
rule ClinicalCOCA:
    input:
        clin="raw_data/clinical/{subject}",
        clust="data/coca/{subject}_{datatype}_COCA.clst"
    output:
        "out/clinical/coca/{subject}_{datatype}_COCA.pval"
    shell:
        "mkdir -p out/clinical/{wildcards.subject};"
        "Rscript dev/analyse_results/ClinicalLabelsEnrichment.R -s {wildcards.subject} -c {input.clust} -f {input.clin} -o {output}"

# Compute Clinical Labels Enrichment Analysis for raw clusterings
rule ClinicalRaw:
    input:
        clin="raw_data/clinical/{subject}",
        clust="data/{subject}/{subject}_{datatype}_{method}.clst"
    output:
        "out/clinical/{subject}/{subject}_{datatype}_{method}.pval"
    shell:
        "mkdir -p out/clinical/{wildcards.subject};"
        "Rscript dev/analyse_results/ClinicalLabelsEnrichment.R -s {wildcards.subject} -c {input.clust} -f {input.clin} -o {output}"

rule AllSurvClinMulti:
    input:
        survival_PINS="out/survival/{subject}/{subject}_multiomics_PINS.pval",
        survival_MCCA="out/survival/{subject}/{subject}_multiomics_MCCA.pval",
        survival_SNF="out/survival/{subject}/{subject}_multiomics_SNF.pval",
        survival_rMKL="out/survival/{subject}/{subject}_multiomics_rMKL.pval",
        survival_NEMO="out/survival/{subject}/{subject}_multiomics_NEMO.pval",
        clinical_PINS="out/clinical/{subject}/{subject}_multiomics_PINS.pval",
        clinical_MCCA="out/clinical/{subject}/{subject}_multiomics_MCCA.pval",
        clinical_SNF="out/clinical/{subject}/{subject}_multiomics_SNF.pval",
        clinical_rMKL="out/clinical/{subject}/{subject}_multiomics_rMKL.pval",
        clinical_NEMO="out/clinical/{subject}/{subject}_multiomics_NEMO.pval",

    output:
        "out/{subject}.surv_clin_multi.log"
    shell:
        "touch {output}"

rule AllSurvClinSingle:
    input:
        survival_PINS="out/survival/{subject}/{subject}_{omic}_PINS.pval",
        survival_SNF="out/survival/{subject}/{subject}_{omic}_SNF.pval",
        survival_rMKL="out/survival/{subject}/{subject}_{omic}_rMKL.pval",
        survival_NEMO="out/survival/{subject}/{subject}_{omic}_NEMO.pval",
        survival_kmeans="out/survival/{subject}/{subject}_{omic}_kmeans.pval",
        clinical_PINS="out/clinical/{subject}/{subject}_{omic}_PINS.pval",
        clinical_SNF="out/clinical/{subject}/{subject}_{omic}_SNF.pval",
        clinical_rMKL="out/clinical/{subject}/{subject}_{omic}_rMKL.pval",
        clinical_NEMO="out/clinical/{subject}/{subject}_{omic}_NEMO.pval",
        clinical_kmeans="out/clinical/{subject}/{subject}_{omic}_kmeans.pval"
    output:
        "out/{subject}_{omic}.surv_clin_single.log"
    shell:
        "touch {output}"

# Run Survival & Clinical analysis for input clusterings
rule AllSurvClin:
    input:
        surv_clin_multi="out/{subject}.surv_clin_multi.log",
        surv_clin_exp="out/{subject}_expression.surv_clin_single.log",
        surv_clin_mirna="out/{subject}_mirna.surv_clin_single.log",
        surv_clin_met="out/{subject}_methylation.surv_clin_single.log"
    output:
        "out/{subject}.allSurvClin.log"
    shell:
        "touch {output}"
