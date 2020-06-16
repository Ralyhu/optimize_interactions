from interaction_graph import *
from pqdict import pqdict
from collections import *
import constants
import random
import numpy as np
import math
import scipy
import time

def optimize_clustering_by_relocations(interaction_graph, estimate_p_plus, estimate_p_minus, clustering):
    """ 
    Hill Climbing algorithm for the MIL problem

    Arguments:
        interaction_graph {InteractionGraph} -- the input interaction graph object
        estimate_p_plus {list} -- list of expected values of interactions of P^+, one float value per each edge 
        estimate_p_minus {list} -- list of expected values of interactions of P^-, one float value for each edge 
        clustering {list} -- initial clustering

    Returns:
        clustering, cur_it {list, int} -- final clustering after the hill climbing optimization, number of completed iterations 
    """
    num_nodes = interaction_graph.get_number_nodes()
    max_iters = constants.maximum_number_iterations_relocation
    graph = interaction_graph.get_graph()
    # ref for avoiding dots
    incident = graph.incident
    edges = graph.es
    has_changed = True
    cur_it = 0
    next_id_cluster = max(clustering) + 1
    while cur_it < max_iters and has_changed:
        # random order for nodes evaluation
        permutation = np.arange(num_nodes)
        np.random.shuffle(permutation)
        permutation = permutation.tolist()
        has_changed = False
        for cur_node in permutation:
            gain_by_cluster = defaultdict(int)
            c_u = clustering[cur_node]
            # get neighbors of current node
            incident_edges = incident(cur_node)
            for edge_id in incident_edges:
                edge = edges[edge_id]
                if edge.source != cur_node:
                    neigh = edge.source
                else: 
                    neigh = edge.target
                puv_plus = estimate_p_plus[edge_id]
                puv_minus = estimate_p_minus[edge_id]
                c_neigh = clustering[neigh]
                # delta plus
                delta = puv_plus - puv_minus
                if c_u == c_neigh:
                    # need delta minus, change sign to delta
                    delta = - delta 
                gain_by_cluster[c_neigh] += delta
            # intra cluster gain by moving to a new cluster
            gain_by_cluster[next_id_cluster] = 0
            
            # compute delta f(cur_node) by taking the max over neighbor clusters
            # first remove the entry associated to current cluster (if cur_node is linked to another node in c_u)
            if c_u in gain_by_cluster:
                delta_cur = gain_by_cluster.pop(c_u)
            else:
                delta_cur = 0
            #assert len(gain_by_cluster.keys()) != 0
            cur_max_value = float("-inf")
            cur_max_cluster = -1
            for neigh_cluster, value in gain_by_cluster.items():
                #assert neigh_cluster != c_u
                if value > cur_max_value:
                    cur_max_value = value
                    cur_max_cluster = neigh_cluster
            #assert cur_max_cluster != -1
            delta_f = delta_cur + cur_max_value
            if delta_f > 0:
                # move cur_node from current cluster to the best one
                has_changed = True
                clustering[cur_node] = cur_max_cluster
                # next id for moving to new cluster
                if cur_max_cluster == next_id_cluster:
                    next_id_cluster += 1
        cur_it += 1
    return clustering, cur_it


def min_cc_degree_ailon_heap_lazy(interaction_graph, estimate_p_plus, estimate_p_minus, final_optimization=False, file_path_times=None):
    """
    D-MIL algorithm with updatable heap in order to manage the sampling of nodes.
    The updates on the heap data structure are made according to a lazy policy for efficiency. After a pivot node is sampled, priorities of neighbor 
    nodes need to be updated. However, since these priorities are non-increasing, their updates are deferred to the moment they become the nodes with
    highest priority. This is equivalent to update priorities as modifications occur but it requires a lower number of updates on the heap.

    Arguments:
        interaction_graph {InteractionGraph} -- the input interaction graph object
        estimate_p_plus {list} -- list of expected values of interactions of P^+, one float value per each edge 
        estimate_p_minus {list} -- list of expected values of interactions of P^-, one float value for each edge 

    Keyword Arguments:
        final_optimization {bool} -- if True HillClimbing step is performed as a postprocessing step (default: {False})
        file_path_times {str} -- if specifies the path where running time information is saved  (default: {None} -- no time info is memorized)

    Returns:
        cluster_membership {list} -- a clustering solution
    """
    start1 = time.time()
    estimate_p_plus = estimate_p_plus.tolist()
    estimate_p_minus = estimate_p_minus.tolist()
    M = interaction_graph.get_M()
    num_nodes = interaction_graph.get_number_nodes()
    cluster_membership = [-1] * num_nodes
    graph = interaction_graph.get_graph()
    degree = interaction_graph.get_degrees()
    # random numbers for each node
    rnds = np.random.uniform(0, 1, size=num_nodes)
    rnds = rnds.tolist()
    data = {}
    for node, degree_value in enumerate(degree):
        data[node] = degree_value * rnds[node]
    pq = pqdict(data, reverse=True)
    # lazy updates, degree to decrement for each node in order to get updated values, 0 means no update is required
    updates = [0] * num_nodes
    cluster_id = 0
    # ref for avoiding dots
    neighbors = graph.neighbors
    incident = graph.incident
    edges = graph.es
    while len(pq) > 0:
        while updates[pq.top()] > 0:
            to_update_node = pq.top()
            degree[to_update_node] -= updates[to_update_node]
            pq.updateitem(to_update_node, (degree[to_update_node]) * rnds[to_update_node])
            updates[to_update_node] = 0
        cur_node = pq.pop()
        cluster_membership[cur_node] = cluster_id
        # get neighbors of current node
        incident_edges = incident(cur_node)
        for edge_id in incident_edges:
            edge = edges[edge_id]
            if edge.source != cur_node:
                neigh = edge.source
            else: 
                neigh = edge.target
            if cluster_membership[neigh] == -1:
                puv_plus = estimate_p_plus[edge.index]
                puv_minus = estimate_p_minus[edge.index]
                w_plus = (1.0/M) * (puv_plus + (M - puv_plus - puv_minus)/2.0)
                w_minus = (1.0/M) * (puv_minus + (M - puv_plus - puv_minus)/2.0)
                #assert math.isclose(w_plus + w_minus, 1.0, rel_tol=1e-6)
                if w_plus > w_minus:
                    # remove element from heap
                    pq.pop(neigh)
                    cluster_membership[neigh] = cluster_id
                    # update degrees of neighbors of this neighbor node
                    for neigh_of_neigh in neighbors(neigh):
                        if cluster_membership[neigh_of_neigh] == -1:
                            updates[neigh_of_neigh] += 1
                else:
                    updates[neigh] += 1
        cluster_id += 1
    total_time1 = time.time() - start1
    if final_optimization:
        start2 = time.time()
        _, n_it = optimize_clustering_by_relocations(interaction_graph, estimate_p_plus, estimate_p_minus, cluster_membership)
        total_time2 = time.time() - start2
        if file_path_times != None:
            file2 = file_path_times + "times2.txt"
            with open(file2, "w+") as f:
                f.write(str(total_time2) + "\n")
                f.write(str(n_it))
    if file_path_times != None:
        file1 = file_path_times + "times1.txt"
        with open(file1, "w+") as f:
            f.write(str(total_time1))
        
    return cluster_membership

def min_cc_ailon(interaction_graph, estimate_p_plus, estimate_p_minus, final_optimization = False, file_path_times=None):
    """
    MIL algorithm

    Arguments:
        interaction_graph {InteractionGraph} -- the input interaction graph object
        estimate_p_plus {list} -- list of expected values of interactions of P^+, one float value per each edge 
        estimate_p_minus {list} -- list of expected values of interactions of P^-, one float value for each edge 

    Keyword Arguments:
        final_optimization {bool} -- if True HillClimbing step is performed as a postprocessing step (default: {False})
        file_path_times {str} -- it specifies the path where running time information is saved  (default: {None} -- no time info is memorized)

    Returns:
        cluster_membership {list} -- a clustering solution
    """
    start1 = time.time()
    estimate_p_plus = estimate_p_plus.tolist()
    estimate_p_minus = estimate_p_minus.tolist()
    M = interaction_graph.get_M()
    num_nodes = interaction_graph.get_number_nodes()
    permutation = np.arange(num_nodes)
    marked_nodes = [False] * num_nodes
    cluster_membership = [-1] * num_nodes
    graph = interaction_graph.get_graph()
    # generate random permutation of nodes (Fisher-Yates algorithm O(n))
    np.random.shuffle(permutation)
    permutation = permutation.tolist()
    cluster_id = 0
    # ref for avoiding dots
    incident = graph.incident
    edges = graph.es
    for cur_node in permutation:
        if not marked_nodes[cur_node]:
            cluster_membership[cur_node] = cluster_id
            marked_nodes[cur_node] = True
            # get neighbors of current node
            incident_edges = incident(cur_node)
            for edge_id in incident_edges:
                edge = edges[edge_id]
                if edge.source != cur_node:
                    neigh = edge.source
                else: 
                    neigh = edge.target
                if not marked_nodes[neigh]:
                    puv_plus = estimate_p_plus[edge.index]
                    puv_minus = estimate_p_minus[edge.index]
                    w_plus = (1.0/M) * (puv_plus + (M - puv_plus - puv_minus)/2.0)
                    w_minus = (1.0/M) * (puv_minus + (M - puv_plus - puv_minus)/2.0)
                    #assert math.isclose(w_plus + w_minus, 1.0, rel_tol=1e-6)
                    if w_plus > w_minus:
                        cluster_membership[neigh] = cluster_id
                        marked_nodes[neigh] = True
            cluster_id += 1
    total_time1 = time.time() - start1
    if final_optimization:
        start2 = time.time()
        _, n_it = optimize_clustering_by_relocations(interaction_graph, estimate_p_plus, estimate_p_minus, cluster_membership)
        total_time2 = time.time() - start2
        if file_path_times != None:
            file2 = file_path_times + "times2.txt"
            with open(file2, "w+") as f:
                f.write(str(total_time2) + "\n")
                f.write(str(n_it))
    if file_path_times != None:
        file1 = file_path_times + "times1.txt"
        with open(file1, "w+") as f:
            f.write(str(total_time1))
    return cluster_membership