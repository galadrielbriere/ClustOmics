# Installation

1) Download and install Neo4j Desktop from https://neo4j.com/.
From the application, create a new project, and a local database  in the project through "Add - Add Local DBMS": set the database name, a password, and the Neo4j version of the database (version tested: 4.0.6). By default, the username to connect to the database is "neo4j". Once the database created, click on it, go to "Plugins", and install "APOC" and "Graph Data Science" Neo4j libraries. 
Open the database settings file ("..." track) and check if the following line is commented, if not, comment it.
```
#dbms.connector.bolt.listen_address=:7687
```

Start the database. Update the "neo4j_password" field in the Snakemake configuration file <i>config.yaml</i> according to the password you chose for the database.


2) Install Conda (https://docs.conda.io/projects/conda/en/latest/index.html). Create conda environment for ClustOmics from ClustOmicsCondaEnv.yml:

```
	git clone https://github.com/galadrielbriere/ClustOmics.git
	cd ClustOmics
	conda env create -f ClustOmicsCondaEnv.yml
	conda activate ClustOmics
	Rscript -e "devtools::install_github('Shamir-Lab/Logrank-Inaccuracies/logrankHeinze')"
```

# Run ClustOmics
ClustOmics uses **Snakemake** for a quick and easy execution. See https://snakemake.readthedocs.io/en/stable/ to get started with Snakemake.

## TLDR
SingleToMulti scenario for AML cancer type:
```
	cd ClustOmics
	conda activate ClustOmics
	mv dataAll data #or change clusterings_folder parameter to 'dataAll' in the Snakemake configuration file config.yaml
	snakemake out/AML.AML_EXP_MIRNA_MET_NEMO_PINS_SNF_rMKL_kmeans_all.FuseClusterings.log --cores 1
	mv data dataAll
```
MultiToMulti scenario for COAD cancer type:
```
	cd ClustOmics
	conda activate ClustOmics
	mv dataOnlyMulti data #or change clusterings_folder parameter to 'dataOnlyMulti' in the Snakemake configuration file config.yaml
	snakemake out/COAD.COAD_MULTI_MCCA_NEMO_PINS_SNF_rMKL.FuseClusterings.log --cores 1
	mv data dataOnlyMulti
```

## 1) Data folder organisation
The _./data_ folder contains directories for each study case. For instance, a specific folder is created for each cancer type analysed. Input clusterings and metadata file are sored in <i>./data/<subject></i>. For instance, all input data for AML cancer type are stored in <i>./data/AML</i>.

Each <i>data/<subject></i> folder must contain an **object metadata file** and **input clusterings files**.

### 1.1) The metadata file
The metadata file stores informations on the objects considered (for instance, Patients for cancer subtyping). This file is used to instanciate the graph. **All objects appearing in at least one input clustering must be described in this file.**
The metadata file must be named as: <i>subject_metadata.txt</i>. For instance, for AML cancer study case, patients are described in <i>./data/AML/AML_metadata.txt</i>.
The metadata file is tab delimited.
In its minimal form, the metadata file must be organised as followed:

```
main_node
<main_node_name>
<node1_id>
<node2_id>
...
<nodeN_id>
```

The first line always start by "main_node". The second line allows you to chose the name to give your main nodes. For instance, for AML, main nodes are "Patients". The 3rd line and the rest of the file must contain each object ids, as they appear in input clusterings. For instance, the TCGA sample ID form AML study case.

To add information about main nodes in the graph, additional columns can be added:

```
main_node	property	node	label
<main_node_name>	<property_name>	<node_name>	<label_name>
<node1_id>	<property_value>	<node_value>	<label_value>
...
```

For instance:

```
main_node	property	node	label	label	label	label	label
Patient	age_at_initial_pathologic_diagnosis	gender	Tissue	expression	mirna	methylation	multiomics
TCGA.AB.2802.03	50	MALE	PrimaryBloodDerivedCancer 		mirna	methylation
TCGA.AB.2803.03	61	FEMALE	PrimaryBloodDerivedCancer 	expression	mirna	methylation	multiomics
```

From this metadata file, ClustOmics will create _Patient_ nodes for each provided TCGA id. Each Patient node will display a property "age_at_initial_pathologic_diagnosis" to indicate the age of the patient at diagnosis. Each Patient node will share a "gender" relationship with a "gender" node "MALE" or "FEMALE" to indicate the gender of the patient. Each Patient node will carry additional labels to indicate the tissue and the omic(s) they were measured for. For instance, Patient "TCGA.AB.2802.03" have been measured for miRNA and methylation but not for expression.

### 1.2) Input clustering files
Input clusterings must be named as follow: <i>subject_datatype_method.clst</i>. For instance, input clustering computed with NEMO from the AML expression dataset is named <i>AML_expression_NEMO.clst</i>.
Clustering files are tab delimited and must contain a header to indicate the main node name (as set the metadata file) and the name to give to cluster nodes.

For instance:
```
Patient	Cluster
TCGA.AB.2803.03	3
TCGA.AB.2805.03	7
```

## 2) Snakemake Configuration file
Fill the configuration file <i>config.yaml</i> according to the specifications indicated in it. The integration scenarios are defined in this configuration file.

## 3) Run ClustOmics
```
	snakemake out/{subject}.{rel_name}.FuseClusterings.log --cores 1
```
Or:
```
	snakemake out/<subject>.BuildGraph.log --cores 1
	snakemake out/{subject}.SupportEdges.log --cores 1
	snakemake out/{subject}.{rel_name}.IntegrationEdges.log --cores 1
	snakemake out/{subject}.{rel_name}.FuseClusterings.log --cores 1
```

For instance, to compute consensus clustering for AML SingleToMulti:
```
	mv dataAll data #or change clusterings_folder parameter to 'dataAll' in the Snakemake configuration file
	snakemake out/AML.AML_EXP_MIRNA_MET_NEMO_PINS_SNF_rMKL_kmeans_all.FuseClusterings.log --cores 1
```

To run ClustOmics with a user-defined number of supports (/!\ parameter min_size_consensus will be ignored),
we recommend choosing the appropriate number of supports threshold by ploting wMQ evolution:

```
	snakemake out/plots/{subject}/wMQ_{subject}.{rel_name}.svg --cores 1
```

This results in, for instance:

![alt tag](https://raw.githubusercontent.com/galadrielbriere/ClustOmics/bf89187112d0ff7662e27623ae8ae4d6ef39bcae/outAll/plots/AML/wMQ_AML.AML_EXP_MIRNA_MET_NEMO_PINS_SNF_rMKL_kmeans_all.svg)

After selecting the number of supports threshold to use, run:
```
	snakemake out/{subject}.{rel_name}.FuseClusterings.{nb_supports}_supports.log --cores 1
```

For instance, on AML SingleToMulti with number of supports treshold of 10:
```
	snakemake out/AML.AML_EXP_MIRNA_MET_NEMO_PINS_SNF_rMKL_kmeans_all.FuseClusterings.10_supports.log  --cores 1
```

To run all analysis (including Survival analysis, clinical label enrichment, ...):
```
	snakemake out/{subject}.{rel_name}.all.log --cores 1
```

For latest command to work, you need to download raw datasets here: http://acgt.cs.tau.ac.il/multi_omic_benchmark/download.html. Files must be stored in <i>./raw_data</i>, and organised as follow:
```
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
├── Etc....

```

We don't recommend using multiple cores (defined by the --cores parameter), unless when running Snakemake rules AllSurvClinMulti, AllSurvClinSingle or AllSurvClin, for which we recommend using as many cores as possible.

## 4) Visualize results in Neo4j browser
From Neo4j Desktop, open the database with Neo4j browser (or enter 'http://localhost:7474/browser/' in your web browser).
From the "Settings" track, bottom left, uncheck "Connect results nodes".
By clicking on the "Database" icon, top left, you can access all node types and relationship types stored in the database. By clicking on each node type of relation type, you can visualize a subset of the corresponding objects.
You can also use the Cypher command box to enter your own visualization queries, using Cypher language.
When running ClustOmics, consensus clustering is created in rule FuseClusterings. At the end of this rule, ClustOmics gives you the Cypher query you need to visualize consensus results (directly on your terminal, or in the rule log file).
For instance:
```
Check the graph with the following Cypher query : 
 MATCH (o:Patient:AML)-[r:FROM_MARKOVCLUST]-(c:OptimalNbSupports:AML:MarkovCluster:AML_EXP_MIRNA_MET_NEMO_PINS_SNF_rMKL_kmeans_all) RETURN o, r, c
```
When visualizing the result of this query, you can choose to color nodes or relationship according to labels of interest. Juste click on the label or relationship you want to highlight (on the panel just above the graph), and select the color you want to display (on the pannel just under the graph).

