#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 13 19:13:22 2019

@author: briere
"""
import numpy as np
import itertools
import networkx as nx
import markov_clustering as mc
#pip install markov-clustering
#Not working with scipy 1.3 : must install 1.2 = pip install scipy==1.2.1
#import pkg_resources
#pkg_resources.require("scipy==1.2.1")


def call_louvain_algo(subject, obj, clust, rel_name, nb_supports, driver):
    with driver.session() as session:
        obj1_nodes = "(o1:"+ obj + ":" + subject + ")"
        obj2_nodes = "(o2:"+ obj + ":" + subject + ")"

        cat_graph_name = "_".join([rel_name, str(nb_supports)])

        objects =' MATCH ' + obj1_nodes + '-[r:' + rel_name + ']-' + obj2_nodes + \
        ' WHERE r.nb_supports >= ' + str(nb_supports) + ' RETURN DISTINCT id(o1) as id'

        relations = 'MATCH ' + obj1_nodes + '-[r:' + rel_name + ']->' + obj2_nodes + \
        ' WHERE r.nb_supports >= ' + str(nb_supports) + ' RETURN id(o1) as source, id(o2) as target, r.nb_supports as weight'

        get_integration_graph = "MATCH " + obj1_nodes + "-[r:" + rel_name + "]->" + obj2_nodes + \
            "WHERE r.nb_supports >= " + str(nb_supports) + " RETURN o1.id, r.nb_supports, o2.id"
        int_graph = session.run(get_integration_graph)
        int_graph = [(record['o1.id'], record['o2.id'], {"weight": record['r.nb_supports']}) for record in int_graph]

        G = nx.DiGraph()
        G.add_edges_from(int_graph)
        
        if (G.number_of_edges()==0):
            return(0)
        
        make_graph_cat = "CALL gds.graph.create.cypher('" + cat_graph_name +\
            "', '" + objects + "', '" + relations + "')"
        session.run(make_graph_cat)

        call_louvain = "CALL gds.louvain.stream('" + cat_graph_name +\
            "',{relationshipWeightProperty:'weight', maxIterations:10})" +\
            "YIELD nodeId, communityId RETURN gds.util.asNode(nodeId).id as id, communityId"
        results = session.run(call_louvain)
        results = [[record['id'], record['communityId']] for record in results]

        remove_graph_cat = "CALL gds.graph.drop('" + cat_graph_name + "') YIELD graphName;"
        session.run(remove_graph_cat)

    return(np.array(results))


def call_markov_algo(subject, obj, clust, rel_name, nb_supports, driver):
    with driver.session() as session:

        obj1_nodes = "(o1:"+ obj + ":" + subject + ")"
        obj2_nodes = "(o2:"+ obj + ":" + subject + ")"

        get_integration_graph = "MATCH " + obj1_nodes + "-[r:" + rel_name + "]->" + obj2_nodes + \
            "WHERE r.nb_supports >= " + str(nb_supports) + " RETURN o1.id, r.nb_supports, o2.id"
        int_graph = session.run(get_integration_graph)
        int_graph = [(record['o1.id'], record['o2.id'], {"weight": record['r.nb_supports']}) for record in int_graph]

    G = nx.DiGraph()
    G.add_edges_from(int_graph)
    
    if (G.number_of_edges()==0):
        return(0)
    
    matrix = nx.to_scipy_sparse_matrix(G, weight = "weight")

    result = mc.run_mcl(matrix)           # run MCL with default parameters
    clusters = mc.get_clusters(result)    # get clusters

    o_nodes = []

    for pair in int_graph:
        o1 = pair[0]
        o2 = pair[1]
        if o1 not in o_nodes:
            o_nodes.append(o1)
        if o2 not in o_nodes:
            o_nodes.append(o2)

    clustering = []
    marked = []
    for cluster in clusters:
        nodes = []
        for node in cluster:
            if o_nodes[node] not in marked:
                nodes.append(o_nodes[node])
                marked.append(o_nodes[node])
        clustering.append(nodes)

    result = []
    num_clust = 0
    for cluster in clustering:
        num_clust = num_clust + 1
        for node in cluster:
            result.append([node, num_clust])
    return(np.array(result))


def get_main_results(subject, obj, clust, rel_name, clust_algo, nb_nodes_clust, nb_supports, driver):

    if clust_algo == "Louvain":
        results = call_louvain_algo(subject, obj, clust, rel_name, nb_supports, driver)
    elif clust_algo == "Markov":
        results = call_markov_algo(subject, obj, clust, rel_name, nb_supports, driver)

    if type(results) != np.ndarray: 
        return(None, None, None)
     
    communities = np.unique(results.transpose()[1])

    main_communities = [[community, len(np.where(results == community)[0])] for community in communities if len(np.where(results == community)[0]) >= nb_nodes_clust]
    main_communities = [ [np.int(j) for j in i] for i in main_communities]
    main_communities = np.array(main_communities)
    small_communities = [[community, len(np.where(results == community)[0])] for community in communities if len(np.where(results == community)[0]) < nb_nodes_clust]
    small_communities = [ [np.int(j) for j in i] for i in small_communities]
    small_communities = np.array(small_communities)

    with driver.session() as session:
        obj1_nodes = "(o1:"+ obj + ":" + subject + ")"
        obj2_nodes = "(o2:"+ obj + ":" + subject + ")"
        rel = "-[r:" + rel_name + "]-"
        all_objects = "MATCH " + obj1_nodes + rel + obj2_nodes + " RETURN distinct(o1.id) as id"
        all_objects = session.run(all_objects)
        all_objects = [record['id'] for record in all_objects]

    if len(main_communities) > 1:
        nb_nodes_kept = np.sum(main_communities, axis=0)[1]

        main_results = [[result[0],int(result[1])] for result in results if int(result[1]) in main_communities[:,0]]
        main_results = np.array(main_results)

        if len(small_communities) > 1:
            small = [[result[0],int(result[1])] for result in results if int(result[1]) in small_communities[:,0]]
            small = np.array(small)
        else:
            small = None

        if small is None:
            sum_len = len(main_results)
        else:
            sum_len = len(main_results) + len(small)

        if sum_len < len(all_objects):
            if small is None:
                unclassified = np.setdiff1d(all_objects, main_results[:,0])
                unclassified = [[objid, "unclassified"] for objid in unclassified]
                unclassified = np.array(unclassified)
            else:
                unclassified = np.setdiff1d(all_objects, np.concatenate([main_results[:,0], small[:,0]]))
                unclassified = [[objid, "unclassified"] for objid in unclassified]
                unclassified = np.array(unclassified)

        else:
            unclassified = None

        return(main_results, small, unclassified)

    else:
        return(None, None, None)


def wModularization(subject, obj, clust, rel_name, clust_algo, nb_nodes_clust, nb_supports, max_pos_nb_sup, driver):
    obj_nodes = "(o:"+ obj + ":" + subject + ")"
    obj1_nodes = "(o1:"+ obj + ":" + subject + ")"
    obj2_nodes = "(o2:"+ obj + ":" + subject + ")"
    rel = "[r:" + rel_name + "]"
    rel_condition = "r.nb_supports >=" + str(nb_supports)

    main_results, small_clusters, unclassified = get_main_results(subject, obj, clust, rel_name, clust_algo, nb_nodes_clust, nb_supports, driver)

    if main_results is not None:
        with driver.session() as session:
            nb_nodes_kept = len(main_results)
            clusters = np.unique(main_results[:,1])
            mq_cl = []
            for i in range(0, len(clusters)):
                clust_id = clusters[i]
                other_clusts = np.delete(clusters, i)

                nodes_in_clust = [res[0] for res in main_results if res[1] == clust_id]

                somme_inter = 0
                for oth_clust in other_clusts:
                    nodes_out_clust = [res[0] for res in main_results if res[1] == oth_clust]
                    e_inter_query = "MATCH " + obj1_nodes + "-" + rel + "-" + obj2_nodes + " WHERE " + rel_condition + \
                    " AND o1.id in " + str(nodes_in_clust) + " AND o2.id in " + str(nodes_out_clust) + \
                    " WITH distinct r as dr RETURN sum(dr.nb_supports) as sum"

                    res_e_inter = session.run(e_inter_query)
                    e_inter = int([record['sum'] for record in res_e_inter][0])
                    e_inter = e_inter/(max_pos_nb_sup*len(nodes_in_clust)*len(nodes_out_clust))

                    somme_inter = somme_inter + e_inter

                if len(nodes_in_clust) <= 1:
                    e_intra = 0
                else:
                    # Compute number of edges within the cluster
                    e_intra_query = "MATCH " + obj1_nodes + "-" + rel + "-" + obj2_nodes + " WHERE " + rel_condition + \
                        " AND o1.id in " + str(nodes_in_clust) + " AND o2.id in " + str(nodes_in_clust) + \
                        " WITH distinct r as dr RETURN sum(dr.nb_supports) as sum"

                    res_e_intra = session.run(e_intra_query)
                    e_intra = int([record['sum'] for record in res_e_intra][0])
                    e_intra = 2*e_intra/(max_pos_nb_sup*len(nodes_in_clust)*(len(nodes_in_clust)-1))

                mq = e_intra - 1/(len(clusters)-1) * somme_inter
                mq_cl.append(mq)

        MQ = sum(mq_cl)/len(clusters)
        return(MQ, nb_nodes_kept, len(clusters))

    else:
        return(-float('Inf'), 0, 0)


# MQ = np array listing all ModulariztionQuality values for each possible number of supports
# nb_supp = np array listing each nb_supports used to compute MQ
# nb_nodes_kept = np array listing the number of nodes returned for each number of supports tested
# min_tot_nodes = minimum number of nodes accepted for a clustering
def find_opt_nb_supports(MQ, nb_supp, nb_nodes_kept, min_tot_nodes, epsilon):
    # Only clustering results with sufficient number of nodes
    min_nodes_ok = np.where(nb_nodes_kept >= min_tot_nodes)
    if len(min_nodes_ok[0]) > 0:
        reduce_MQ = MQ[min_nodes_ok]
        reduce_nb_supp = nb_supp[min_nodes_ok]

        max_MQ = max(reduce_MQ)
        interval_l = max_MQ - epsilon*max_MQ

        ok_MQ = np.where(reduce_MQ >= interval_l)
        # Keep MQ obtained with max nb_supports
        opt_MQ = reduce_MQ[ok_MQ[0][-1]]
        opt_nb_sup = reduce_nb_supp[ok_MQ[0][-1]]

        return(opt_nb_sup, opt_MQ)

    else:
        exit("ERROR: The 'min_size_consensus' parameter is too high (see config file)")


def reassign(subject, obj, rel_name, main_results, unclassified, driver):
    with driver.session() as session:
        results = []

        for ind in unclassified[:,0]:
            ind_node = "(o1:" + subject + ":" + obj + "{id:'" + ind + "'})"
            obj_node = "(o2:" + subject + ":" + obj + ")"
            rel = "-[r:" + rel_name + "]-"
            all_rels_querry = "MATCH " + ind_node + rel + obj_node + " RETURN o1.id, o2.id, r.nb_supports"
            all_rels = session.run(all_rels_querry)
            all_rels = [[record['o1.id'], record['o2.id'], record['r.nb_supports']] for record in all_rels]
            all_rels = np.array(all_rels)

            sum_clusts = []
            clusters = np.unique(main_results[:,1])
            for cluster in clusters:
                members = main_results[np.where(main_results[:,1]==cluster),0]
                members_rel = np.in1d(all_rels[:,1], members)
                scores = all_rels[members_rel,2]
                scores = [int(score) for score in scores]
                if len(scores) > 0:
                    sum_clusts.append(sum(scores)/len(scores))
                else:
                    sum_clusts.append(0)

            reassign_clust = clusters[sum_clusts.index(max(sum_clusts))]
            results.append([ind, reassign_clust])

    return(np.array(results))
