#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 29 13:56:49 2019

@author: briere
"""

import logging
import argparse
from neo4j import GraphDatabase


descr = ('This script creates IntegrationEdges summarizing methods and datatypes that led to',
         ' similar predictions for each pair of main nodes.',
         'IntegrationEdges are specific to one integration, defined with a list of datatypes to ',
         'integrate and a list of methods to integrate.')

parser = argparse.ArgumentParser(description=descr)
parser.add_argument("-id", "--neo_id", help="Neo4j ID used to connect to the database")
parser.add_argument("-pwd", "--neo_pwd", help="Neo4j password used to connect to the database")
parser.add_argument("-host", "--neo_localhost", help="Neo4j graph localhost")
parser.add_argument("-subj", "--subject", help="Analysis subject")
parser.add_argument("-obj_name", "--object_nodes_name",
                    help=" Name of the Object nodes in the graph\
                        (main nodes, defined in the metadata file)")
parser.add_argument("-out", "--log_out", help="Log file")
parser.add_argument("-op", "--operator",
                    help="'AND' or 'OR' (default). If set to 'AND', stringent integration : \
                    IntegrationEdges are created only if the different datatypes lead to the same prediction \
                    (pair classified in same cluster for each datatype). If set to 'OR' (default value), IntegrationEdges \
                    are create if a least one datatype classifies the pair in the same cluster.")
parser.add_argument("-rel", "--rel_name", help="Name of the IntegrationEdges")
parser.add_argument("-dt", "--datatypes", help="Datatypes to integrate (pipe | separated)")
parser.add_argument("-met", "--methods", help="Clustering methods to integrate (pipe | separated)")

args = parser.parse_args()

neo_id = args.neo_id
neo_pwd = args.neo_pwd
neo_localhost = args.neo_localhost
subject = args.subject
object_nodes_name = args.object_nodes_name
rel_name = args.rel_name

datatypes = args.datatypes
methods = args.methods
datatypes = list(datatypes.split('|'))
methods = list(methods.split('|'))

operator = args.operator
log_file = args.log_out
logging.basicConfig(filename=log_file, level=logging.INFO, format='%(message)s')


if operator not in ('AND', 'OR'):
    logging.error('\noperator parameter must be set to "OR" or "AND"')

logging.info('\nIntegrating clusterings from : %s' % (subject))
logging.info('\nOn datatypes : %s' % (str(datatypes)))
logging.info('\nFor methods : %s' % (str(methods)))
logging.info('\nRelationship type : %s' % (rel_name))

def create_rel_to_query(driver, subject, obj, datatypes, methods, rel_name, operator="OR"):
    '''
    Parameters
    ----------
    subject : STRING
        Subject of analysis (e.g: AML to run the analysis form AML cancer).
    obj : STRING
        Main node name (e.g: Patient)
    datatypes : LIST OF STRINGS
        List of datatypes to fuse\
            (e.g: ["methylation", "mirna"] for the integration of methylation and mirna clusterings).
    methods : LIST OF STRINGS
        List of clustering methods to fuse\
            (e.g: ["k-means", "hierarchical"] for the integration of k-means and hierarchical\
                produced clusterings).
    rel_name : STRING
        Name of the integration edges to compute.
    operator : STRING, optional
        DESCRIPTION. The default is "OR". It means an integration edge will be added\
            if co-clustering is observed in one of the given datatypes.\
            If set to "AND", an integration edge will be added only if all considered\
            datatypes co-cluster the pair of objects.

    Returns
    -------
    Compute integration edges for the given integration scenario\
        (set of datatypes and methods to fuse).

    '''

    graph = driver.session()

    # Methods
    method_rel = "[r"
    for method in methods:
        if method == methods[0]:
            method_rel = method_rel + ":" + method
        else:
            method_rel = method_rel + "|" + method
    method_rel = method_rel + "]"

    # Condition
    condition = "r."
    if operator == "OR":
        or_and = " OR "
    else:
        or_and = " AND "
    # Datatypes
    for datatype in datatypes:
        if datatype == datatypes[0]:
            condition = condition + datatype + "=1"
        else:
            condition = condition + or_and + "r." + datatype + "=1"

    q_obj1 = "(o1:"+ obj + ":" + subject + ")"
    q_obj2 = "(o2:"+ obj + ":" + subject + ")"

    return_values = ""
    data_values = ""
    datatypes_str = ""
    for datatype in datatypes:
        if datatype == datatypes[0]:
            return_values = return_values + "r." + datatype
            data_values = data_values + 'sum(r.' + datatype + ") as " + datatype
            datatypes_str = datatypes_str + datatype
        else:
            return_values = return_values = return_values + ", r." + datatype
            data_values = data_values + ', sum(r.' + datatype + ") as " + datatype
            datatypes_str = datatypes_str + "+" + datatype

    new_rel = "[nr:" + rel_name + "]"

    merge_new_rel_query = "MATCH " + q_obj1 + "-" + method_rel + "-" + q_obj2 + \
        " WHERE " + condition + " WITH o1.id as o1, o2.id as o2, " + data_values + \
        " WITH o1 as o1, o2 as o2, " + datatypes_str + " as nb_supports" + \
        " MATCH (n1), (n2) WHERE n1.id = o1 and n2.id = o2 " +\
        " MERGE (n1)-" + new_rel + "-(n2) SET nr.nb_supports = nb_supports"

    graph.run(merge_new_rel_query)

    graph.close()

def main():
    '''
    Parameters
    ----------
    subject : STRING
        Subject of analysis (e.g: AML to run the analysis form AML cancer).
    obj : STRING
        Main node name (e.g: Patient)
    datatypes : LIST OF STRINGS
        List of datatypes to fuse\
            (e.g: ["methylation", "mirna"] for the integration of methylation and mirna clusterings).
    methods : LIST OF STRINGS
        List of clustering methods to fuse\
            (e.g: ["k-means", "hierarchical"] for the integration of k-means and hierarchical\
                produced clusterings).
    rel_name : STRING
        Name of the integration edges to compute.
    operator : STRING, optional
        DESCRIPTION. The default is "OR". It means an integration edge will be added\
            if co-clustering is observed in one of the given datatypes.\
            If set to "AND", an integration edge will be added only if all considered\
            datatypes co-cluster the pair of objects.

    Returns
    -------
    Compute integration edges for the given integration scenario\
        (set of datatypes and methods to fuse).

    '''

    logging.info('Building the integration graph. Integration Edges will be named "%s".', rel_name)
    print('Building the integration graph, this could take a while.'+\
          ' Integration Edges will be named "' + rel_name + '".')

    driver = GraphDatabase.driver(uri=neo_localhost, auth=(neo_id, neo_pwd), max_connection_lifetime=14000)
    create_rel_to_query(driver, subject, object_nodes_name,\
                        datatypes, methods, rel_name, operator)
    driver.close()

if __name__ == "__main__":
    main()
