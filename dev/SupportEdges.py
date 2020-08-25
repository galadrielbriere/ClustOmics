#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 29 09:21:28 2019

@author: briere
"""

import logging
import argparse
from itertools import combinations
from neo4j import GraphDatabase

descr = ('This script creates SupportEdges for each pair of subject nodes.'
         'SupportEdges link nodes classified in the same cluster for the various'
         'clustering methods and summarize the datatypes that led to this prediction.')

parser = argparse.ArgumentParser(description=descr)
parser.add_argument("-id", "--neo_id", help="Neo4j ID used to connect to the database")
parser.add_argument("-pwd", "--neo_pwd", help="Neo4j password used to connect to the database")
parser.add_argument("-host", "--neo_localhost", help="Neo4j graph localhost")
parser.add_argument("-subj", "--subject", help="Analysis subject")
parser.add_argument("-cls_name", "--cluster_nodes_name",
                    help="Name of the Cluster nodes in the graph (defined in the metadata file)")
parser.add_argument("-obj_name", "--object_nodes_name", help=" Name of the Object nodes in"
                    "the graph (main nodes, defined in the metadata file)")
parser.add_argument("-out", "--log_out", help="Log file")
args = parser.parse_args()

neo_id = args.neo_id
neo_pwd = args.neo_pwd
neo_localhost = args.neo_localhost
subject = args.subject
object_nodes_name = args.object_nodes_name
cluster_nodes_name = args.cluster_nodes_name
log_file = args.log_out

logging.basicConfig(filename=log_file, level=logging.INFO, format='%(message)s')

logging.info('Computing support edges for %s analysis...' % (subject))

def compute_support_edges(subject, obj, clust, driver):
    '''
    Parameters
    ----------
    subject : STRING
        Subject of analysis (e.g: AML to run the analysis form AML cancer).
    obj : STRING
        Main node name (e.g: Patient)
    clust : STRING
        Cluster node name (e.g: Cluster)

    Returns
    -------
    Computes support edges for all raw clustering results.

    '''

    match_clust = "(c:"+ clust + ":" + subject + ")"
    get_all_clust = 'MATCH ' + match_clust + ' RETURN distinct(c.id) AS id'
    graph = driver.session()
    all_clust = graph.run(get_all_clust)
    all_clust = [cluster["id"] for cluster in all_clust]
    [get_supports(subject, obj, clust, cluster, graph) for cluster in all_clust]

    graph.close()

def get_supports(subject, obj, clust, cluster, graph):
    '''
    Parameters
    ----------
    subject : STRING
        Subject of analysis (e.g: AML to run the analysis form AML cancer).
    obj : STRING
        Main node name (e.g: Patient)
    clust : STRING
        Cluster node name (e.g: Cluster)
    cluster : STRING
        A cluster id.

    Returns
    -------
    None.

    '''

    print("Computing Support Edges for clustering " + "_".join(cluster.split("_")[:-1]))
    logging.info('Computing Support Edges for cluster %s.', cluster)
    datatype = cluster.split('_')[1]
    method = cluster.split('_')[2]


    clust_node = "(c1:"+ clust + ":" + subject + " {id:'" + cluster + "'})"
    obj_node = "(o:"+ obj + ")"
    rel = "[r:" + method + "]"

    query = "MATCH " + clust_node + "-[p:PART_OF]-" + obj_node + \
        "WITH collect(o) as objs " + \
        " WITH apoc.coll.combinations(objs, 2, 2) as combis " + \
        " UNWIND combis as comb WITH comb[0] as o1, comb[1] as o2 " + \
        " MERGE (o1)-" + rel + "-(o2) SET r." + datatype + " = 1"

    graph.run(query)


def main():
    '''
    Parameters
    ----------
    neo_id : STRING
        Neo4j database connection ID.
    neo_pwd : STRING
        Password to connect to Neo4j database.
    neo_localhost : STRING
        Database localhost adress (e.g: http://localhost:7474).
    subject : STRING
        Subject of analysis (e.g: AML to run the analysis form AML cancer).
    obj : STRING
        Main node name (e.g: Patient)
    clust : STRING
        Cluster node name (e.g: Cluster)

    Returns
    -------
    Compute support edges for each raw clustering.
    '''

    driver = GraphDatabase.driver(uri=neo_localhost, auth=(neo_id, neo_pwd))
    compute_support_edges(subject, object_nodes_name, cluster_nodes_name, driver)
    driver.close()


if __name__ == "__main__":
    main()
