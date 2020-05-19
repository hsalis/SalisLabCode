import sys

import vercov

import networkx as nx

from itertools import count
from time      import time
from math      import log10

def is_graph_empty(homology_graph):
    if not homology_graph.number_of_nodes():
        return True
    return False

def is_graph_complete(homology_graph):
    num_nodes = homology_graph.number_of_nodes()
    num_edges = homology_graph.number_of_edges()
    if float(num_edges) == num_nodes * (num_nodes - 1) * 0.5:
        return True
    return False

def completex_elimination(homology_graph, verbose):
    ### REMOVE THIS ###
    return set()
    ### REMOVE THIS ###
    
    indiset_nodes   = []
    eliminated_subs = 0
    retained_subs   = 0
    nodes_to_delete = []

    if verbose:
        clear_length = 0
        print '\nCompletex processing'

    for homology_sub_graph in nx.connected_component_subgraphs(homology_graph):
        if is_graph_complete(homology_sub_graph):
            eliminated_subs += 1
            indiset_nodes.append(homology_sub_graph.nodes().iterkeys().next())
            nodes_to_delete.append(homology_sub_graph.nodes())
        else:
            retained_subs += 1
        if verbose:
            printed = ' [Components processed = {}] [Eliminated = {}] [Retained = {}]\r'.format(eliminated_subs + retained_subs, eliminated_subs, retained_subs)
            sys.stdout.write(' '*clear_length+'\r')
            sys.stdout.write(printed)
            clear_length = len(printed)
            sys.stdout.flush()
            # print ' [Components processed = {}] [Eliminated = {}] [Retained = {}]\r'.format(eliminated_subs + retained_subs, eliminated_subs, retained_subs)
    for nodeset_to_delete in nodes_to_delete:
        homology_graph.remove_nodes_from(nodeset_to_delete)
    
    if verbose: print
    
    return indiset_nodes

def get_powertex_degree(homology_graph, powertex_coeff, verbose):
    max_degree = max(len(homology_graph[node]) for node in homology_graph)
    avg_degree = (homology_graph.number_of_edges() * 1.0) / homology_graph.number_of_nodes()
    ptx_degree = max_degree - max_degree*powertex_coeff + avg_degree*powertex_coeff # y = -(M - A)*x + M
    try:
        ptx_degree = ptx_degree if int(log10(max_degree)) - int(log10(avg_degree)) > 0 else 0
    except:
        ptx_degree = 0
    if verbose:
        print ' [Max degree = {}] [Avg degree = {}] [Powertex degree = {}]'.format(max_degree, avg_degree, ptx_degree)
    return ptx_degree

def powertex_elimination(homology_graph, verbose):
    if verbose:
        clear_length    = 0
        print '\nPowertex processing'

    powertex_coeff      = 0.3
    powertex_degree     = get_powertex_degree(homology_graph, powertex_coeff, verbose)
    eliminated_powertex = []
    processed_nodes     = 0

    if powertex_degree > 0:
        for processed_nodes, node in enumerate(homology_graph.nodes()):
            # In case we are dealing with a complete graph
            if homology_graph.number_of_nodes() - len(eliminated_powertex) == 1:
                break
            if len(homology_graph[node]) > powertex_degree:
                eliminated_powertex.append(node)
                if verbose:
                    printed = ' [Nodes processed = {}] [Powertex eliminated = {}] [Retained = {}]\r'.format(processed_nodes+1, len(eliminated_powertex), homology_graph.number_of_nodes()-len(eliminated_powertex))
                    sys.stdout.write(' '*clear_length+'\r')
                    sys.stdout.write(printed)
                    clear_length = len(printed)
                    sys.stdout.flush()
    
    if verbose:
        if not eliminated_powertex:
            print ' [Nodes processed = {}] [Powertex eliminated = 0] [Retained = {}]\r'.format(processed_nodes, homology_graph.number_of_nodes())
        else:
            print

    homology_graph.remove_nodes_from(eliminated_powertex)

def get_vercov_func(vercov_func, homology_graph):
    # TO DO: NEED TO USE DIFFERENT FUNCTIONS AT DIFFERENT EDGE COUNT SCALES!!
    if vercov_func == 'nrpG':
        return vercov.nrp_vercov_greedy, 'NRP Greedy'
    elif vercov_func == '2apx':
        # return vercov.std_vercov_approx, 'Standard 2-approximation'
        return vercov.nx_vercov, 'Standard 2-approximation'
    else:
        return vercov.nrp_vercov_approx, 'NRP 2-approximation'

def dump_homology_graph(homology_graph, graph_file):
    nx.write_adjlist(homology_graph, graph_file)

def load_homology_graph(graph_file):
    return nx.read_adjlist(graph_file, nodetype=int)

def get_recovered_non_homologs(homology_graph, graph_file, vercov_func=None, verbose=True):
    indiset_nodes   = set()
    # powertex_elimination(homology_graph, verbose)
    completex_nodes = [] #completex_elimination(homology_graph, verbose)
    # indiset_nodes.update(completex_nodes)
    
    possible_nodes  = set(homology_graph.nodes())   
    
    if verbose:
        print '\n [+] Initial independent set = {} [{} completex], computing vertex cover on remaining {} nodes.'.format(len(indiset_nodes), len(completex_nodes), len(possible_nodes))
    
    if is_graph_empty(homology_graph):
        if verbose:
            print ' [X] Graph is empty, further independent set expansion not possible, terminating.'
    else:
        vercov_func, vercov_func_name = get_vercov_func(vercov_func, homology_graph)
        iteration = -1

        if verbose:
            print ' [+] Vertex Cover Function: {}'.format(vercov_func_name)
            sys.stdout.write(' [+] Dumping graph into: {}'.format(graph_file))

        t0 = time()
        dump_homology_graph(homology_graph, graph_file)
        if verbose:
            sys.stdout.write(' in {} seconds\n'.format(time()-t0))

        while True:
            iteration += 1

            if verbose:
                print '\n----------------------'
                print 'Now running iteration: {}'.format(iteration)
                print '----------------------'

            t0 = time()

            if iteration > 0:
                homology_graph = nx.Graph(homology_graph.subgraph(possible_nodes))

            vercov_nodes = vercov_func(homology_graph, verbose)
            
            if verbose:
                print '\n [+] Computed vertex cover of size: {} (in {:.4} seconds)'.format(len(vercov_nodes), time()-t0)
                print ' [+] Loading graph from: {}'.format(graph_file)

            homology_graph = load_homology_graph(graph_file)

            new_indiset_nodes = possible_nodes - vercov_nodes
            indiset_nodes |= new_indiset_nodes

            possible_nodes = vercov_nodes
            prev_possibility_count = len(possible_nodes)
            for indi_node in new_indiset_nodes:
                possible_nodes.difference_update(homology_graph[indi_node])
            curr_possibility_count = len(possible_nodes)

            if verbose:
                print ' [+] Current independent set size:  {}'.format(len(indiset_nodes))
                print ' [+] Potential nodes for expansion: {} (projected independent set size: {})'.format(len(possible_nodes), len(indiset_nodes)+len(possible_nodes))
            if len(possible_nodes) == 0 or prev_possibility_count == curr_possibility_count:
                if verbose:
                    print ' [X] Cannot expand independent set, terminating.'
                break

    return indiset_nodes