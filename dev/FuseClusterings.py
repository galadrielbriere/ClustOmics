#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 15 10:17:32 2020

@author: briere
"""

import argparse
import logging
from GraphMetrics import *
from neo4j import GraphDatabase

descr = '''This script computes the optimal number of supports to cluster the IntegrationGraph.
Once filtered, the graph is clustered using either the Louvain community detection algorithm
or the Markov Clustering algorithm.
The resulting clustering is an integrated clustering and is stored in file.'''

parser = argparse.ArgumentParser(description=descr)
parser.add_argument("-id", "--neo_id",
                    help="Neo4j ID used to connect to the database")
parser.add_argument("-pwd", "--neo_pwd",
                    help="Neo4j password used to connect to the database")
parser.add_argument("-host", "--neo_localhost",
                    help="Neo4j graph localhost")
parser.add_argument("-subj", "--subject",
                    help="Analysis subject")
parser.add_argument("-cls_name", "--cluster_nodes_name",
                    help="Name of the Cluster nodes in the graph (defined in the metadata file)")
parser.add_argument("-obj_name", "--object_nodes_name",
                    help=" Name of the Object nodes in the graph\
                    (main nodes, defined in the metadata file)")
parser.add_argument("-out", "--log_out",
                    help="Log file")
parser.add_argument("-out_cls", "--clustering_out",
                    help="Output clustering file")
parser.add_argument("-rel", "--rel_name",
                    help="Name of the IntegrationEdges")
parser.add_argument("-dt", "--datatypes",
                    help="Datatypes to integrate (pipe | separated)")
parser.add_argument("-met", "--methods",
                    help="Clustering methods to integrate (pipe | separated)")
parser.add_argument("-min_clust", "--min_nodes_clust",
                    help="Minimum number of nodes for a cluster to be returned")
parser.add_argument("-min_nodes", "--min_nodes_analyse",
                    help="Minimum number of nodes for a clustering to be returned")
parser.add_argument("-nb_sup", "--user_nb_supports",
                    help="If 0, NeOmics coputes the optimal number of supports to filter the graph \
                    else, the provided nuber of supports is used.")

args = parser.parse_args()

out = args.clustering_out
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
nb_nodes_clust = int(args.min_nodes_clust)
min_tot_nodes = int(args.min_nodes_analyse)
user_nb_supports = args.user_nb_supports
log_file = args.log_out
logging.basicConfig(filename=log_file, level=logging.INFO, format='%(message)s')

logging.info('Integrating clusterings from : %s' % (subject))
logging.info('On datatypes : %s' % (str(datatypes)))
logging.info('For methods : %s' % (str(methods)))
logging.info('Relationship type : %s' % (rel_name))

def print_and_log(string):
    print(string)
    logging.info(string)

def fuse_clusterings(subject, obj, clust, rel_name, nb_nodes_clust, min_tot_nodes,
                     methods, datatypes, driver, out, opt=True, user_nb_sup=None):
    '''

        Parameters
        ----------
        subject : STRING
            Subject of analysis (e.g: AML to run the analysis form AML cancer).
        obj : STRING
            Main nodes name (e.g: Patient)
        clust : STRING
            Cluster nodes name.
        rel_name : STRING
            Name of integratio edges.
        nb_nodes_clust : INTEGER
            Minimum number of nodes in a cluster for it to be returned (min_size_clust).
        min_tot_nodes : INTEGER
            Minimum number of nodes to be partitioned in the consensus (min_size_consensus).
        methods : LIST OF STRINGS
            Input clusterings to consider (algorithm).
        datatypes : LIST OF STRINGS
            Input clusterings to consider (datatype).
        driver : NEO4J PYTHON DRIVER
            driver = GraphDatabase.driver(uri=neo_localhost, auth=(neo_id, neo_pwd)).
        out : STRING
            Output filename.
        opt : BOOL, optional
            Should the number of supports threshold be computed automatically ? The default is True.
            If set to False, user_nb_sup will be used to filter the graph
            and the min_tot_nodes_parameter will be ignored.
        user_nb_sup : INT, optional
            If opt=False, which number of supports threshold should be used ? The default is None.

        Returns
        -------
        Consensus clustering.

        '''

    clust_algos = ["Markov", "Louvain"]

    with driver.session() as session:
        max_nb_sup_in_graph = session.run("MATCH ()-[r:" + rel_name + "]->() RETURN max(r.nb_supports) as max")
        max_nb_sup_in_graph = int([record['max'] for record in max_nb_sup_in_graph][0])

    if opt == True:
        MQ_louvain_markov = []
        nb_supp_louvain_markov = []
        for clust_algo in clust_algos:
            nb_supp = []
            MQ = []
            nb_nodes_kept = []
            nb_clusters = []

            for nb_supports in range(1, max_nb_sup_in_graph+1):
                nb_supp.append(nb_supports)

                mq, nb_nodes, nb_clust = weighted_modularity(subject, obj, clust, rel_name, clust_algo, nb_nodes_clust, nb_supports, driver)

                MQ.append(mq)
                nb_nodes_kept.append(nb_nodes)
                nb_clusters.append(nb_clust)

            MQ = np.asarray(MQ)
            nb_nodes_kept = np.asarray(nb_nodes_kept)
            nb_clusters = np.asarray(nb_clusters)
            nb_supp = np.asarray(nb_supp)

            nb_supports, max_MQ = find_opt_nb_supports(MQ, nb_supp, nb_nodes_kept, min_tot_nodes)
            MQ_louvain_markov.append(max_MQ)
            nb_supp_louvain_markov.append(nb_supports)

        MQ_louvain_markov = np.asarray(MQ_louvain_markov)
        nb_supp_louvain_markov = np.asarray(nb_supp_louvain_markov)

        max_MQ_index = np.where(MQ_louvain_markov == max(MQ_louvain_markov))

        if len(max_MQ_index[0]) > 1:
            max_nb_supports = np.where(nb_supp_louvain_markov == max(nb_supp_louvain_markov))
            if len(max_nb_supports[0]) > 1:
                nb_supports = nb_supp_louvain_markov[1]
                clust_algo = 'Markov'
            else:
                nb_supports = nb_supp_louvain_markov[max_nb_supports[0][0]]
                clust_algo = clust_algos[max_nb_supports[0][0]]
        else:
            nb_supports = nb_supp_louvain_markov[max_MQ_index[0][0]]
            clust_algo = clust_algos[max_MQ_index[0][0]]

        print_and_log("Optimal number of supports is %s" % nb_supports)
        print_and_log("The graph will be clustered with %s clustering using only IntegrationEdges with nb_supports >= %s" % (clust_algo, nb_supports))

    else:
        nb_supports = user_nb_sup
        MQ_louvain_markov = []
        nb_nodes_kept = []
        nb_clusters = []
        nb_supp_louvain_markov = [user_nb_sup, user_nb_sup]
        for clust_algo in clust_algos:
            mq, nb_nodes, nb_clust = weighted_modularity(subject, obj, clust, rel_name, clust_algo, nb_nodes_clust, nb_supports, driver)
            MQ_louvain_markov.append(mq)
            nb_nodes_kept.append(nb_nodes)
            nb_clusters.append(nb_clust)
        MQ_louvain_markov = np.asarray(MQ_louvain_markov)
        nb_nodes_kept = np.asarray(nb_nodes_kept)
        nb_clusters = np.asarray(nb_clusters)
        nb_supp_louvain_markov = np.asarray(nb_supp_louvain_markov)

        max_MQ_index = np.where(MQ_louvain_markov == max(MQ_louvain_markov))

        if len(max_MQ_index[0]) > 1:
            nb_supports = nb_supp_louvain_markov[1]
            clust_algo = 'Markov'
        else:
            nb_supports = nb_supp_louvain_markov[max_MQ_index[0][0]]
            clust_algo = clust_algos[max_MQ_index[0][0]]

    if clust_algo == "Louvain":
        node_clust_name = "LouvainCommunity"
        rel_clust_name = "FROM_COMMUNITY"
    elif clust_algo == "Markov":
        node_clust_name = "MarkovCluster"
        rel_clust_name = "FROM_MARKOVCLUST"

    if opt:
        communities = "(c:OptimalNbSupports:"  + subject + ":" + node_clust_name + ":" + rel_name + ")"
    else:
        communities = "(c:UserNbSupports:"  + subject + ":" + node_clust_name + ":" + rel_name + "{nb_supports:" + str(nb_supports) + "})"

    results = get_main_results(subject, obj, clust, rel_name, clust_algo, nb_nodes_clust, nb_supports, driver)

    if results is not None:
        nb_nodes_kept = len(results)
        nb_clusters = len(np.unique(results[:, 1]))
        print_and_log('%s nodes classified in %s clusters after remooving too small clusters (min_accepted_nb_nodes_in_clusters set to %s in the configuration file)' % (nb_nodes_kept, nb_clusters, nb_nodes_clust))
        print_and_log("Weighted Modularization Quality for the clustering is %s" % str(max(MQ_louvain_markov)))
        print_and_log('Storing results into the graph')
        results_to_neo4j(results, subject, obj, rel_name, clust_algo, nb_supports, driver, opt)

        obj_nodes = "(o:"+ obj + ":" + subject + ")"
        check_graph_query = "MATCH " + obj_nodes + "-[r:" + rel_clust_name + "]-"  + communities + " RETURN o, r, c"
        print_and_log('Check the graph with the following Cypher query : \n %s' % check_graph_query)

        head = obj + "\t" + clust
        file = open(out, "w")
        np.savetxt(file, results, fmt='%s', header=head, comments='', delimiter='\t')
        file.close()

        print_and_log('Clustering results stored in file : %s' % out)
    else:
        exit("ERROR: Can not cluster the graph (only 1 cluster found). Please, use a lower number of supports.")


def results_to_neo4j(results, subject, obj, rel_name, clust_algo, nb_supports, driver, opt=False):

    if clust_algo == "Louvain":
        node_clust_name = "LouvainCommunity"
        rel_clust_name = "FROM_COMMUNITY"
    elif clust_algo == "Markov":
        node_clust_name = "MarkovCluster"
        rel_clust_name = "FROM_MARKOVCLUST"

    for community in np.unique(results[:, 1]):
        community_id = rel_name + "_" + community
        if opt == True:
            community_node = "(c:OptimalNbSupports:"  + subject + ":" + node_clust_name + ":" + rel_name + \
            " {id: '" +  community_id + "'" + ", nb_supports: " + str(nb_supports) + ", clust:'" + community + "'})"

        else:
            community_node = "(c:UserNbSupports:"  + subject + ":" + node_clust_name + ":" + rel_name + \
            " {id: '" +  community_id + "'" + ", nb_supports: " + str(nb_supports) + ", clust:'" + community + "'})"

        create_community_node = "MERGE " + community_node

        with driver.session() as session:
            session.run(create_community_node)

            # Get all nodes belonging to same community
            nodes = [node for node in results[np.where(results == str(community))[0], :][:, 0]]

            for node in nodes:
                # Make relationship between community node and patients nodes
                make_rel = "MATCH (o:" + obj + ":" + subject + " {id: '" + str(node) + "'}) " + \
                "MATCH " + community_node + " MERGE (o)-[r:" + rel_clust_name + "]-(c) RETURN o, r, c"
                session.run(make_rel)


def main():
    driver = GraphDatabase.driver(uri=neo_localhost, auth=(neo_id, neo_pwd))
    if user_nb_supports == 'False':
        fuse_clusterings(subject, obj, clust, rel_name, nb_nodes_clust, min_tot_nodes, methods, datatypes, driver, out, opt=True, user_nb_sup=None)
    else:
        fuse_clusterings(subject, obj, clust, rel_name, nb_nodes_clust, min_tot_nodes, methods, datatypes, driver, out, opt=False, user_nb_sup=int(user_nb_supports))
    driver.close()



if __name__ == "__main__":
    main()
