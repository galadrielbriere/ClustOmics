

INSTALLATION

1) Install Neo4j, create an empty database and install "APOC" and "Graph Data Science" Neo4j libraries.

2) Create conda environment for ClustOmics:

	conda create -n MyClustEnv python=3.8
	conda activate MyClustEnv
	conda install -c bioconda snakemake 
 	conda install -c conda-forge neo4j-python-driver 
	conda install numpy networkx matplotlib

To run analysis on generated consensus clusterings (Survival analysis, Clinical label enrichment, ...):

	conda install r-essentials r-base r-survival r-optparse r-survminer r-FactoMineR r-pca3d
	conda install -c bioconda bioconductor-genefilter  bioconductor-affy 


3) Install Markov Clustering from: https://github.com/GuyAllard/markov_clustering
	python -m pip install markov_clustering

RUN CLUSTOMICS

1) DATA FOLDER
The ./data folder contains directories for each study case. For instance, a specific folder is created for each cancer type analysed. Input clusterings and metadata file are sored in ./data/<subject>. For instance, all input data for AML cancer type are stored in ./data/AML. 

Each data/<subject> folder must contain an object metadata file and input clusterings files. 

1.1) The metadata file
The metadata file stores informations on the objects considered. This file is used to instanciate the graph. All objects appearing in at least one input clustering must be described in this file. 
The metadata file must be named as: <subject>_metadata.txt. For instance, for AML cancer study case, patients are described in ./data/AML/AML_metadata.txt.
The metadata file is tab delimited. 
In its minimal form, the metadata file must be organised as followed:

main_node
<main_node_name>
<node1_id>
<node2_id>
...
<nodeN_id>

The first line always start by "main_node". The second line allows you to chose the name to give your main_nodes. For instance, for AML, main_nodes are "Patients". The 3rd line and the rest of the file must contain each object ids, as they appear in input clusterings. For instance, the TCGA sample ID form AML study case. 

To add information about main_nodes in the graph, additional columns can be added:

main_node	property	node	label
<main_node_name>	<property_name>	<node_name>	<label_name>
<node1_id>	<property_value>	<node_value>	<label_value>
...

For instance: 

main_node	property	node	label	label	label	label	label
Patient	age_at_initial_pathologic_diagnosis	gender	Tissue	expression	mirna	methylation	multiomics
TCGA.AB.2802.03	50	MALE	PrimaryBloodDerivedCancer 		mirna	methylation	
TCGA.AB.2803.03	61	FEMALE	PrimaryBloodDerivedCancer 	expression	mirna	methylation	multiomics

From this metadata file, ClustOmics will create Patient nodes for each provided TCGA id. Each Patient node will display a property "age_at_initial_pathologic_diagnosis" to indicate the age of the patient at diagnosis. Each Patient node will share a "gender" relationship with a "gender" node "MALE" or "FEMALE" node to indicate the gender of the Patient. Each Patient node will carry additional labels to indicate the Tissue and the omic they were measured for. For instance, Patient "TCGA.AB.2802.03" have been measured for miRNA and methylation but not for expression, but not for expression.

1.2) Input clustering files
Input clusterings must be named as follow: <subject>_<datatype>_<method>.clst. For instance, input clustering computed with NEMO from the expression dataset is named AML_expression_NEMO.clst. 
Clustering files are tab delimited and must contain a header to indicate the main_node name (as in the metadata file) and the name to give to cluster nodes.

For instance:
Patient Cluster
TCGA.AB.2803.03	3
TCGA.AB.2805.03	7

2) CONFIGURATION FILE
Fill the configuration file config.yaml according to the specifications indicated in it.

3) RUN CLUSTOMICS

	snakemake out/{subject}.{rel_name}.FuseClusterings.log --cores 1

Or:

	snakemake out/<subject>.BuildGraph.log --cores 1
	snakemake out/{subject}.SupportEdges.log --cores 1
	snakemake out/{subject}.{rel_name}.IntegrationEdges.log --cores 1
	snakemake out/{subject}.{rel_name}.FuseClusterings.log --cores 1

To run ClustOmics with a user-defined number of supports (/!\ parameter min_size_consensus will be ignored):

	snakemake out/{subject}.{rel_name}.FuseClusterings.{nb_supports}_supports.log

To run all analysis (including Survival analysis, clinical label enrichment, ...):

	snakemake out/{subject}.{rel_name}.all.log --cores 1

For it to work, you need to download raw datsets here: http://acgt.cs.tau.ac.il/multi_omic_benchmark/download.html. Files must be stored in ./raw_data, and organised as follow:
./raw_data/
├── AML
│   ├── exp
│   ├── methy
│   ├── mirna
│   └── survival
├── BIC
│   ├── exp
│   ├── methy
│   ├── mirna
│   └── survival
├── clinical
│   ├── AML
│   ├── BIC
│   ├── COAD
│   ├── GBM
│   ├── KIRC
│   ├── LIHC
│   ├── LUSC
│   ├── OV
│   ├── SARC
│   └── SKCM
├── COAD
│   ├── exp
│   ├── methy
│   ├── mirna
│   └── survival
├── GBM
│   ├── exp
│   ├── methy
│   ├── mirna
│   └── survival
├── KIRC
│   ├── exp
│   ├── methy
│   ├── mirna
│   └── survival
├── LIHC
│   ├── exp
│   ├── methy
│   ├── mirna
│   └── survival
├── LUSC
│   ├── exp
│   ├── methy
│   ├── mirna
│   └── survival
├── OV
│   ├── exp
│   ├── methy
│   ├── mirna
│   └── survival
├── SARC
│   ├── exp
│   ├── methy
│   ├── mirna
│   └── survival
└── SKCM
    ├── exp
    ├── methy
    ├── mirna
    └── survival




We don't recommend using multiple cores (defined by the --cores parameter).



