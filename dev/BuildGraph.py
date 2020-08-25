#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 25 13:08:45 2019

@author: briere
"""

import logging
import argparse
import os
from neo4j import GraphDatabase

descr = ('This script instanciate a Neo4j graph database from the metadata '
         'file and the various clustering results to integrate.')

parser = argparse.ArgumentParser(description=descr)
parser.add_argument("-id", "--neo_id", help="Neo4j ID used to connect to the database")
parser.add_argument("-pwd", "--neo_pwd", help="Neo4j password used to connect to the database")
parser.add_argument("-host", "--neo_localhost", help="Neo4j graph localhost")
parser.add_argument("-cls", "--clustering_folder", help="Path to directory containing clustering "
                    "results files and metadata file \n Please note that cluster files must be "
                    "names following the pattern : {SUBJECT}_{DATATYPE}_{METHOD}.clst \n"
                    "The metadata file must be names following the pattern : "
                    "{SUBJECT}_metadata.txt \n Please, read the docs to make sure your "
                    "files respect the mandatory format.")
parser.add_argument("-subj", "--subject", help="Analysis subject")
parser.add_argument("-out", "--log_out", help="Log file")

args = parser.parse_args()

neo_id = args.neo_id
neo_pwd = args.neo_pwd
neo_localhost = args.neo_localhost
results_path = args.clustering_folder
subject = args.subject
log_file = args.log_out

logging.basicConfig(filename=log_file, level=logging.INFO, format='%(message)s')

logging.info('Loading results from %s in %s \n' % (subject, results_path))

def make_graph(results_path, subject, driver):
    '''
    Parameters
    ----------
    results_path : STRING
        Path to folder containing raw clustering results.
    subject : STRING
        Subject of analysis (e.g: AML to run the analysis form AML cancer).

    Returns
    -------
    Log file. Instanciates the graph database with data from metadata file
        and clustering results files.
    '''

    files = os.listdir(os.path.join(results_path, subject))
    metadata_file = [file for file in files if "metadata" in file][0]
    logging.info('Initiating graph from %s metadata.', metadata_file)
    create_object_nodes(os.path.join(results_path, subject, metadata_file), subject, driver)
    clustering_files = [file[0:-5] for file in files if file[-5:len(file)] == ".clst"]
    for file in clustering_files:
        datatype = file.split('_')[1]
        method = file.split('_')[2]
        clust_file = str(file) + '.clst'
        logging.info('Loading %s clustering results.',
                     os.path.join(results_path, subject, clust_file))
        logging.info('- computed on datatype : %s', datatype)
        logging.info('- computed with method : %s', method)
        create_cluster_nodes(os.path.join(results_path, subject, clust_file),
                             subject, datatype, method, driver)


def create_object_nodes(metadata_file, subject, driver):
    '''
    Parameters
    ----------
    metadata_file : STRING
        Path to metadata file. Must be named as follow:
            {SUBJECT}_metadata.txt
    subject : STRING
        Subject of analysis (e.g: AML to run the analysis form AML cancer).

    Returns
    -------
    Instanciate graph with data from metadata file.
    metadata_file : rows = "obects" nodes ids
                    cols = metadata nodes ids (with header as node names)
    On first line, for each column, specify how it should be represented in the graph:
              - "main_nodes" for the objects to cluster (e.g Patients TCGA id)
              - as a node linked to main_nodes: "node"
              - as a label of main_nodes: "label"
              - as a property of main nodes : "prop"
    '''

    print("Instanciating graph for " + subject + "from metadata file : " + metadata_file)
    print("Check the graph in your Neo4j browser by querying : \n" +
          "MATCH (o:" + subject + ") RETURN o LIMIT 100")

    graph = driver.session()

    with open(metadata_file) as file:
        line = file.readline()
        line = line.strip()
        col_types = line.split('\t')
        for col_type in col_types[1:len(col_types)]:
            if col_type not in ('label', 'node', 'property'):
                logging.error('Error : Wrong column type. Must be "label",'
                              '"node" or "property". Please, modify the metadata file. \n')

        line = file.readline()
        line = line.strip()
        col_names = line.split('\t')
        main_node = col_names[0]

        line = file.readline()

        while line:
            line = line.strip()
            data = line.split('\t')
            # Create object node
            node_id = data[0]

            node = "(o:" + subject + ":" + main_node + "{id:'" + node_id + "'})"

            graph.run("MERGE" + node)

            # Get metadata
            for i in range(1, len(data)):
                if col_types[i] == "property" and data[i] != '':
                    graph.run("MATCH" + node + "SET o." + col_names[i] + "='" + data[i] +"'")
                elif col_types[i] == "label" and data[i] != '':
                    graph.run("MATCH" + node + "SET o:" + data[i])
                elif col_types[i] == "node" and data[i] != '':
                    meta_node = "(m:" + col_names[i] + "{id:'" + data[i] + "'})"
                    graph.run("MERGE" + meta_node)
                    meta_rel = "[r:" + col_names[i] + "]"
                    graph.run("MATCH" + node + " MATCH" + meta_node +\
                              " CREATE (o)-" + meta_rel + "->(m)")
                elif col_types[i] not in ["property", "label", "node"]:
                    logging.info('"%s" is not a valid keyword.' % (col_types[i]))
                    logging.info('Please use "property", "label" or "node". Skipping.')
                    next
            line = file.readline()

    graph.close()

def create_cluster_nodes(clust_file, subject, datatype, method, driver):
    '''
    Parameters
    ----------
    clust_file : STRING
        Path to a raw clustering result file. This file must be named as follow:
            {SUBJECT}_{DATATYPE}_{METHOD}.clst
    subject : STRING
        Subject of analysis (e.g: AML to run the analysis form AML cancer).
    datatype : STRING
        Datatype used to compute this clustering.
    method : STRING
        Method used to compute this clustering.

    Returns
    -------
    Compute edges between main_nodes and their respective cluster nodes.

    '''

    graph = driver.session()

    with open(clust_file) as file:
        line = file.readline()
        line = line.strip()
        col_names = line.split('\t')
        main_node = col_names[0]
        clust_node = col_names[1]
        clust_labels = ":".join([clust_node, subject, datatype, method])
        id_base = "_".join([subject, datatype, method])

        print("Loading input clustering results for " + subject + "from file: " + clust_file)
        print("Check the graph in your Neo4j browser by querying : \n" +
          "MATCH (o:" + main_node + ":" + subject + ")-[r:PART_OF]-(c:"+\
              clust_labels + ") RETURN * LIMIT 100")

        line = file.readline()
        while line:
            node_id = line.split('\t')[0]
            clust_nb = str(line.split('\t')[1].strip())
            clust_id = id_base + "_" + str(clust_nb)

            node = "(o:" + subject + ":" + main_node + "{id:'" + node_id + "'})"
            clust = "(c:" + clust_labels + "{id:'" + clust_id + "'})"
            graph.run("MERGE" + clust)
            graph.run("MATCH" + clust + " MATCH" + node +\
                      " CREATE (o)-[:PART_OF]->(c)")

            line = file.readline()

    graph.close()

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
    Returns
    -------
    Instanciate graph from metadata file and raw clustering files.
    '''

    driver = GraphDatabase.driver(uri=neo_localhost, auth=(neo_id, neo_pwd))
    make_graph(results_path, subject, driver)
    driver.close()

if __name__ == "__main__":
    main()

