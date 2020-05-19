import sys

import networkx as nx
import networkx.algorithms.approximation as nxa

from collections import deque
from random      import choice, shuffle

def nx_vercov(homology_graph, verbose):
    return nxa.vertex_cover.min_weighted_vertex_cover(homology_graph)

def get_max_degree_candidate_nodes(homology_graph, candidate_nodes):
    max_degree, max_degree_nodes = 0, []
    for node, degree in candidate_nodes:
        if degree == max_degree:
            max_degree_nodes.append(node)
        elif degree > max_degree:
            max_degree, max_degree_nodes = degree, [node]
    shuffle(max_degree_nodes)
    return max_degree_nodes

def get_max_degree_nodes(homology_graph, candidate_nodes=None):
    if candidate_nodes is None:
        return get_max_degree_candidate_nodes(homology_graph, candidate_nodes=homology_graph.degree())
    else:
        return get_max_degree_candidate_nodes(homology_graph, candidate_nodes=homology_graph.degree(candidate_nodes))

def get_pendants_from_node(homology_graph, candidate_node, verbose):
    pendant_list = []
    for neighbor in homology_graph[candidate_node]:
        if len(homology_graph[neighbor]) <= 2:
            pendant_list.append(neighbor)
    pendant_list.sort(key=lambda node: len(homology_graph[node])-1)
    return deque(pendant_list)

def remove_node_and_pendants(homology_graph, candidate_node, vercov_nodes, verbose):
    if candidate_node in homology_graph:
        pendant_deque = get_pendants_from_node(homology_graph, candidate_node, verbose)
        homology_graph.remove_node(candidate_node)
        if pendant_deque:
            eliminate_pendants(homology_graph, pendant_deque, vercov_nodes, verbose)

def remove_edge_and_pendants(homology_graph, candidate_edge, vercov_nodes, verbose):
    for node in candidate_edge:
        remove_node_and_pendants(homology_graph, node, vercov_nodes, verbose)

def eliminate_pendants(homology_graph, pendant_deque, vercov_nodes, verbose):
    while pendant_deque:
        candidate_node = pendant_deque.popleft()
        if candidate_node in homology_graph:
            if len(homology_graph[candidate_node]) == 1:
                if verbose:
                    print '  [x] Pendant node {} eliminated'.format(candidate_node)
                adj_node = homology_graph[candidate_node].iterkeys().next()
                homology_graph.remove_edge(adj_node, candidate_node)
                vercov_nodes.add(adj_node)
                new_candidates = [node for node in homology_graph[adj_node] if len(homology_graph[node]) <= 2]
                if new_candidates:
                    new_candidates.sort(key=lambda x: len(homology_graph[x])-1)
                    pendant_deque.extend(new_candidates)
                homology_graph.remove_nodes_from([candidate_node, adj_node])
            elif len(homology_graph[candidate_node]) == 0:
                if verbose:
                    print '  [x] Isolated node {} eliminated'.format(candidate_node)
                    homology_graph.remove_node(candidate_node)
        else:
            if verbose:
                print '  [+] Node {} covered'.format(candidate_node)

def get_pendant_deque(homology_graph):
    return deque(sorted([node for node in homology_graph if len(homology_graph[node]) <= 1], key=lambda x: len(homology_graph[x])))

def init_vercov(homology_graph, verbose):
    vercov_nodes = set()

    if verbose:
        print '\n Pendant checking is in progress...'
    pendant_deque = get_pendant_deque(homology_graph)

    if not pendant_deque:
        if verbose:
            print '  [x] No pendants found'
    else:
        if verbose:
            print '  [+] {} Pendants found\n\n Pendant elimination initiated...'.format(len(pendant_deque))
        eliminate_pendants(homology_graph, pendant_deque, vercov_nodes, verbose)
        if verbose:
            print
    
    return vercov_nodes

def get_arbitrary_neighbor(homology_graph, candidate_node, verbose):
    if len(homology_graph[candidate_node]) > 0:
        return homology_graph[candidate_node].iterkeys().next()
    return None

def get_sorted_nodes(homology_graph, verbose):
    return sorted(homology_graph.nodes(), key=lambda node: len(homology_graph[node]))

def std_vercov_approx(homology_graph, verbose):
    vercov_nodes = init_vercov(homology_graph, verbose)

    sorted_nodes = get_sorted_nodes(homology_graph, verbose)

    while homology_graph.number_of_edges():
        if verbose:
            print ' [Edges remaining = {}] [Nodes remaining = {}]'.format(homology_graph.number_of_edges(), homology_graph.number_of_nodes())

        node_i = sorted_nodes.pop()
        if node_i in homology_graph:
            node_j = get_arbitrary_neighbor(homology_graph, node_i, verbose)
            if not node_j is None:
                vercov_nodes.update([node_i, node_j])
                remove_edge_and_pendants(homology_graph, candidate_edge=(node_i, node_j), vercov_nodes=vercov_nodes, verbose=verbose)
            else:
                homology_graph.remove_node(node_i)

    return vercov_nodes

def nrp_vercov_approx(homology_graph, verbose):
    vercov_nodes = init_vercov(homology_graph, verbose)

    while homology_graph.number_of_edges():
        if verbose:
            print ' [Edges remaining = {}] [Nodes remaining = {}]'.format(homology_graph.number_of_edges(), homology_graph.number_of_nodes())

        max_degree_nodes = get_max_degree_nodes(homology_graph)

        for max_node in max_degree_nodes:
            if max_node in homology_graph:
                max_neighbor = get_max_degree_nodes(homology_graph, candidate_nodes=homology_graph.neighbors(max_node))
                if max_neighbor:
                    max_neighbor = max_neighbor.pop()
                    vercov_nodes.update([max_node, max_neighbor])
                    remove_edge_and_pendants(homology_graph, candidate_edge=(max_node, max_neighbor), vercov_nodes=vercov_nodes, verbose=verbose)
                else:
                    homology_graph.remove_node(max_node)

    return vercov_nodes

def nrp_vercov_greedy(homology_graph, verbose):
    vercov_nodes = init_vercov(homology_graph, verbose)

    while homology_graph.number_of_edges():
        if verbose:
            print ' [Edges remaining = {}] [Nodes remaining = {}]'.format(homology_graph.number_of_edges(), homology_graph.number_of_nodes())

        max_degree_nodes = get_max_degree_nodes(homology_graph)

        for max_node in max_degree_nodes:
            if max_node in homology_graph:
                vercov_nodes.add(max_node)
                remove_node_and_pendants(homology_graph, candidate_node=max_node, vercov_nodes=vercov_nodes, verbose=verbose)
    
    return vercov_nodes