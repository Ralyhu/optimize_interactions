from igraph import *
from scipy.stats import truncnorm 
from scipy.stats import *
import numpy as np
from numpy.random import binomial as binomial
import constants

class InteractionGraph:
    
    # main constructor which builds an interaction graph from an igraph graph object and maximum strength of interaction M
    def __init__(self, graph, M = 1):
        self.M = M
        self.graph = graph
        degrees = graph.degree()
        # plus one because they are multiplied by some probability in order to sample by degree
        self.degrees = [d+1 for d in degrees]
        self.p_plus = np.array(graph.es["p+"])
        self.p_minus = np.array(graph.es["p-"])

    @staticmethod
    def get_instance_from_file(path):
        with open(path) as f:
            header = f.readline()
            header_splitted = header.split(" ")
            n_nodes = int(header_splitted[0])
            n_edges = int(header_splitted[1])
            g = Graph()
            g.add_vertices(n_nodes)
            edges = [None] * n_edges
            p_plus_means = [None] * n_edges
            p_minus_means = [None] * n_edges
            id_edge = 0
            for line in f.readlines():
                line_splitted = line.split(" ")
                i = int(line_splitted[0])
                j = int(line_splitted[1])
                p_plus = float(line_splitted[2])
                p_minus = float(line_splitted[3])
                edges[id_edge] = (i, j)
                p_plus_means[id_edge] = p_plus
                p_minus_means[id_edge] = p_minus
                id_edge += 1
            g.add_edges(edges)
            g.es["p+"] = p_plus_means
            g.es["p-"] = p_minus_means
            return InteractionGraph(g)

    # compute expected cumulative loss 
    def analytically_expected_loss(self, cluster_membership, p_plus, p_minus):
        graph = self.graph
        n_edges = self.get_number_edges()
        same_cluster_mask = np.array([cluster_membership[edge.source] == cluster_membership[edge.target] for edge in graph.es])
        different_cluster_mask = (np.ones(shape=n_edges) - same_cluster_mask)
        interactions = same_cluster_mask * p_plus + different_cluster_mask * p_minus
        loss = self.M * n_edges -  np.sum(interactions)
        return loss   
    
    def get_graph(self):
        return self.graph

    def get_number_nodes(self):
        return len(self.graph.vs)

    def get_degrees(self):
        # return a copy of degrees
        return list(self.degrees)

    def get_number_edges(self):
        return len(self.graph.es)

    def get_M(self):
        return self.M

    def get_p_plus(self):
        return self.p_plus
    
    def get_p_minus(self):
        return self.p_minus

    def print(self):
        print(self.graph)

    def get_approximate_condition(self):
        sum_weights = 0.0
        sum_weights += np.sum(self.p_plus)
        sum_weights += np.sum(self.p_minus)
        return sum_weights
    
    def get_minimum_loss_possible(self):
        # discarding non-linked pairs losses
        return self.get_number_edges() - np.sum(np.maximum(self.p_minus, self.p_plus))
    
    
