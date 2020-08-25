#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 16 09:25:04 2019

@author: briere
"""

import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import numpy as np
import argparse
import logging

from GraphMetrics import *
from neo4j import GraphDatabase


# neo_id = 'neo4j'
# neo_pwd = "4jneo"
# neo_localhost = "bolt://localhost:7687"
# obj = "Patient"
# clust = "Cluster"
# rel_name = "AML_EXP_MIRNA_MET_NEMO_PINS_SNF_rMKL"
# datatypes = "expression|methylation|mirna"
# methods = 'PINS|NEMO|SNF|rMKL'
# #datatypes='multiomics'
# #methods='MCCA|PINS|NEMO|SNF|rMKL'
# datatypes = list(datatypes.split('|'))
# methods = list(methods.split('|'))
# nb_nodes_clust = 8
# clust_algo = "Markov"
# opt=True

# graph = Graph(neo_localhost,auth=(neo_id, neo_pwd))
# host = neo_localhost.replace('://', '://' + neo_id + ":" + neo_pwd + "@") + "/db/data"

# obj = "Patient"
# clust = "Cluster"
# nb_nodes_clust = 8
# subject = "BIC"
# rel_name = "BIC_MULTI_MCCA_NEMO_PINS_SNF_rMKL"
# datatypes='multiomics'
# methods='MCCA|PINS|NEMO|SNF|rMKL'
# min_tot_nodes = 589
# nb_supports = 4
# algo = "Louvain"

# mq, nodes, clust = weighted_modularity(subject, obj, clust, rel_name, algo, nb_nodes_clust, nb_supports, graph, host)
# print(mq, nodes, clust)

descr = 'Describe the IntegrationGraph filtered with different nb_supports values.'

parser = argparse.ArgumentParser(description = descr)
parser.add_argument("-id", "--neo_id", help="Neo4j ID used to connect to the database")
parser.add_argument("-pwd", "--neo_pwd", help="Neo4j password used to connect to the database")
parser.add_argument("-host", "--neo_localhost", help="Neo4j graph localhost")
parser.add_argument("-subj", "--subject", help="Analysis subject")
parser.add_argument("-cls_name", "--cluster_nodes_name", help="Name of the Cluster nodes in the graph (defined in the metadata file)")
parser.add_argument("-obj_name", "--object_nodes_name", help=" Name of the Object nodes in the graph (main nodes, defined in the metadata file)")
parser.add_argument("-out", "--plot_out", help="Output plot file")
parser.add_argument("-rel", "--rel_name", help="Name of the IntegrationEdges")
parser.add_argument("-dt", "--datatypes", help="Datatypes to integrate (pipe | separated)")
parser.add_argument("-met", "--methods", help="Clustering methods to integrate (pipe | separated)")
parser.add_argument("-min_clust", "--min_nb_node_in_cluster", help="Minimum accepted number of nodes a cluster must contain to be returned")
parser.add_argument("-min_nodes", "--min_nodes_analyse", help="Minimum number of nodes for a clustering to be returned")

args = parser.parse_args()

out = args.plot_out
neo_id = args.neo_id
neo_pwd = args.neo_pwd
neo_localhost = args.neo_localhost
subject = args.subject
obj = args.object_nodes_name
clust = args.cluster_nodes_name
rel_name = args.rel_name
datatypes = args.datatypes
methods = args.methods
datatypes = list(datatypes.split('|'))
methods = list(methods.split('|'))
nb_nodes_clust = int(args.min_nb_node_in_cluster)
min_tot_nodes = int(args.min_nodes_analyse)


driver = GraphDatabase.driver(uri=neo_localhost, auth=(neo_id, neo_pwd))

nb_supp = []

louvain_MQ = []
louvain_nb_nodes_kept = []
louvain_nb_clusters = []

markov_MQ = []
markov_nb_nodes_kept = []
markov_nb_clusters = []

for nb_supports in range(1, len(methods)*len(datatypes)+1):
    nb_supp.append(nb_supports)

    louvain_mq, louvain_nodes, louvain_clust = weighted_modularity(subject, obj, clust, rel_name, 'Louvain', nb_nodes_clust, nb_supports, driver)
    markov_mq, markov_nodes, markov_clust = weighted_modularity(subject, obj, clust, rel_name, 'Markov', nb_nodes_clust, nb_supports, driver)

    louvain_MQ.append(louvain_mq)
    louvain_nb_nodes_kept.append(louvain_nodes)
    louvain_nb_clusters.append(louvain_clust)

    markov_MQ.append(markov_mq)
    markov_nb_nodes_kept.append(markov_nodes)
    markov_nb_clusters.append(markov_clust)

driver.close()

louvain_nb_nodes_kept = np.asarray(louvain_nb_nodes_kept)
markov_nb_nodes_kept = np.asarray(markov_nb_nodes_kept)
louvain_MQ = np.asarray(louvain_MQ)
markov_MQ = np.asarray(markov_MQ)
nb_supp = np.asarray(nb_supp)

louvain_opt, max_mq_louvain = find_opt_nb_supports(louvain_MQ, nb_supp, louvain_nb_nodes_kept, min_tot_nodes)
markov_opt, max_mq_markov = find_opt_nb_supports(markov_MQ, nb_supp, markov_nb_nodes_kept, min_tot_nodes)

fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.plot(nb_supp, louvain_MQ, 'c', linestyle = "-")
ax1.plot(nb_supp, markov_MQ, 'm', linestyle = "-")
#ax1.plot((louvain_opt, louvain_opt), (min_MQ, louvain_MQ[louvain_opt-1]), 'c', linestyle="-")
ax1.set_ylabel('Modularization Quality')

ax2 = ax1.twinx()
ax2.plot(nb_supp, louvain_nb_nodes_kept, 'c', linestyle = ":")
ax2.plot(nb_supp, markov_nb_nodes_kept, 'm', linestyle = ":")
ax2.set_ylabel('Number of nodes returned')

for i, txt in enumerate(markov_nb_clusters):
    ax1.annotate(txt, (nb_supp[i], markov_MQ[i]), color='m')

for i, txt in enumerate(louvain_nb_clusters):
    ax1.annotate(txt, (nb_supp[i], louvain_MQ[i]), color='c')

plt.axhline(y=min_tot_nodes, linewidth=1, color='k', linestyle="--")
plt.axvline(x=louvain_opt, linewidth=1, color='c', linestyle="--")
plt.axvline(x=markov_opt, linewidth=1, color='m', linestyle="--")

magenta_line = mlines.Line2D([], [], color='m',
                          markersize=5, label='MQ Markov & number of clusters')
cyan_line = mlines.Line2D([], [], color='c',
                          markersize=5, label='MQ Louvain & number of clusters')
magenta_line2 = mlines.Line2D([], [], color='m', linestyle=':',
                          markersize=5, label='Number of nodes Markov')
cyan_line2 = mlines.Line2D([], [], color='c', linestyle=':',
                          markersize=5, label='Number of nodes Louvain')
threshold_line = mlines.Line2D([], [], color='k', linestyle="--", label='Threshold population size')
nb_sup_louvain_line = mlines.Line2D([], [], color='c', linestyle="--", label='Optimal number of supports Louvain')
nb_sup_markov_line = mlines.Line2D([], [], color='m', linestyle="--", label='Optimal number of supports Markov')

legend = plt.legend(handles=[magenta_line, cyan_line, magenta_line2, cyan_line2, threshold_line, nb_sup_louvain_line, nb_sup_markov_line], loc='center left', bbox_to_anchor=(0., 0.5, 0., -1.7),
           ncol=1, fancybox=True)
plt.title('Modularization Quality according to the nb_supports threshold \n (' + subject + ' cancer)')
plt.savefig(out, bbox_inches='tight')
plt.close()

print('Graphic saved. See : %s' % out)

