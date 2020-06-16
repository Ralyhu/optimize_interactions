from interaction_graph import *
from igraph import *
import random
import constants
import numpy as np

def count_clusters(clustering):
    clusters = set()
    for c in clustering:
        clusters.add(c)
    return len(clusters)

def count_majority_positive(interaction_graph):
    graph = interaction_graph.get_graph()
    n_edges = interaction_graph.get_number_edges()
    cont = 0
    for edge in graph.es:
        p_plus = edge["p+"]
        p_minus = edge["p-"]
        if p_plus > p_minus:
            cont += 1
    return cont, n_edges

def count_internal_external_edges(graph, cluster_membership):
    n_edges = len(graph.es)
    same_cluster_mask = np.array([cluster_membership[edge.source] == cluster_membership[edge.target] for edge in graph.es])
    intra_community_edges = same_cluster_mask.sum()
    inter_community_edges = n_edges - intra_community_edges
    return intra_community_edges, inter_community_edges