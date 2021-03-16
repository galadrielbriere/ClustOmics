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

with driver.session() as session:
    clust_node = "(c:" + clust + ":" + subject + ")"
    datatypes_cond = "c:"
    for datatype in datatypes:
        if datatype == datatypes[0]:
            datatypes_cond = datatypes_cond + datatype
        else:
            datatypes_cond = datatypes_cond + " OR c:" + datatype
    methods_cond = "c:"
    for method in methods:
        if method == methods[0]:
            methods_cond = methods_cond + method
        else:
            methods_cond = methods_cond + " OR c:" + method
    max_pos_nb_sup = session.run("MATCH " + clust_node + " WHERE (" + datatypes_cond + ") AND (" + methods_cond + ") RETURN count(distinct labels(c)) as max")
    max_pos_nb_sup = int([record['max'] for record in max_pos_nb_sup][0])

for nb_supports in range(1, max_pos_nb_sup+1):
    nb_supp.append(nb_supports)

    louvain_mq, louvain_nodes, louvain_clust = wModularization(subject, obj, clust, rel_name, 'Louvain', nb_nodes_clust, nb_supports, max_pos_nb_sup, driver)
    markov_mq, markov_nodes, markov_clust = wModularization(subject, obj, clust, rel_name, 'Markov', nb_nodes_clust, nb_supports, max_pos_nb_sup, driver)

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

louvain_opt, max_mq_louvain = find_opt_nb_supports(louvain_MQ, nb_supp, louvain_nb_nodes_kept, min_tot_nodes, 0.05)
markov_opt, max_mq_markov = find_opt_nb_supports(markov_MQ, nb_supp, markov_nb_nodes_kept, min_tot_nodes, 0.05)

fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.plot(nb_supp, louvain_MQ, 'r', linestyle = "-")
ax1.plot(nb_supp, markov_MQ, 'b', linestyle = "-")
ax1.set_ylabel('Modularization Quality')

ax2 = ax1.twinx()
ax2.plot(nb_supp, louvain_nb_nodes_kept, 'r', linestyle = ":")
ax2.plot(nb_supp, markov_nb_nodes_kept, 'b', linestyle = ":")
ax2.set_ylabel('Number of nodes returned')

for i, txt in enumerate(markov_nb_clusters):
    ax1.annotate(txt, (nb_supp[i], markov_MQ[i]), color='b')

for i, txt in enumerate(louvain_nb_clusters):
    ax1.annotate(txt, (nb_supp[i], louvain_MQ[i]), color='r')

plt.axhline(y=min_tot_nodes, linewidth=1, color='k', linestyle="--")
plt.axvline(x=louvain_opt, linewidth=1, color='r', linestyle="--")
plt.axvline(x=markov_opt, linewidth=1, color='b', linestyle="--")

magenta_line = mlines.Line2D([], [], color='b',
                          markersize=5, label='MQ Markov & number of clusters')
cyan_line = mlines.Line2D([], [], color='r',
                          markersize=5, label='MQ Louvain & number of clusters')
magenta_line2 = mlines.Line2D([], [], color='b', linestyle=':',
                          markersize=5, label='Number of nodes Markov')
cyan_line2 = mlines.Line2D([], [], color='r', linestyle=':',
                          markersize=5, label='Number of nodes Louvain')
threshold_line = mlines.Line2D([], [], color='k', linestyle="--", label='Threshold population size')
nb_sup_louvain_line = mlines.Line2D([], [], color='r', linestyle="--", label='Optimal number of supports Louvain')
nb_sup_markov_line = mlines.Line2D([], [], color='b', linestyle="--", label='Optimal number of supports Markov')

legend = plt.legend(handles=[magenta_line, cyan_line, magenta_line2, cyan_line2, threshold_line, nb_sup_louvain_line, nb_sup_markov_line], loc='center left', bbox_to_anchor=(0., 0.5, 0., -1.7),
           ncol=1, fancybox=True)
plt.title('Modularization Quality according to the nb_supports threshold \n (' + subject + ')')
plt.savefig(out, bbox_inches='tight')
plt.close()

print('Graphic saved. See : %s' % out)
