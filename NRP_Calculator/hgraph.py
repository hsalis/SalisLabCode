import adder
import utils
import sys

import networkx as nx

from time             import time
from itertools        import islice
from collections      import deque, defaultdict

def uniquify_seq_list(seq_list):
    unique_seq_index = {}
    index = len(seq_list)-1
    while index > -1:
        seq = seq_list.pop()
        seq.upper().replace('U', 'T')
        unique_seq_index[seq] = index
        index -= 1
    unique_seq_list = []
    while unique_seq_index:
        unique_seq_list.append(unique_seq_index.popitem()[::-1])
    return unique_seq_list

def get_background_filter(homology, background_list, verbose):
    return adder.init_back_db(background_list, homology, seq_list=[], target_list=[0], verbose=verbose)[1]

def get_repeat_dict(seq_deque, background_list, homology, internal_repeats, verbose):
    repeat_dict  = {}
    back_bfilter = None

    if background_list:
        back_bfilter = get_background_filter(homology, background_list, verbose)

    if verbose:
        clear_length = 0

    if not back_bfilter is None:
        if verbose: print

    while seq_deque:
        
        if verbose:
            printed = ' [Sequence processing remaining] = {}\r'.format(len(seq_deque))
            sys.stdout.write(' '*clear_length+'\r')
            sys.stdout.write(printed)
            clear_length = len(printed)
            sys.stdout.flush()

        seq_id, seq = seq_deque.popleft()
        hmers  = list(utils.stream_kmers(seq, k=homology))
        rhmers = map(utils.get_revcomp, hmers)

        # Resolve internal repeats
        if not internal_repeats:
            hmers_set = set(hmers)

            if len(hmers_set) < len(hmers): # Resolve direct internal repeat
                continue

            inverted_repeats_found = False
            for rhmer in rhmers:
                if rhmer in hmers_set:
                    inverted_repeats_found = True
                    break

            if inverted_repeats_found:      # Resolve inverted internal repeat
                continue

        # Resolve background repeats
        background_repeats_found = False
        if not back_bfilter is None:
            for i in xrange(len(hmers)):
                fkmer = min(hmers[i], rhmers[i])
                if fkmer in back_bfilter:
                    # print 'Seq {} conflicts with background_list'.format(seq)
                    background_repeats_found = True
                    break
            if background_repeats_found:
                continue

        # Populate repeat_dict for sequences sharing repeats
        for i in xrange(len(hmers)):
            hmer  = hmers[i]
            rhmer = rhmers[i]
            mmer  = min(hmer, rhmer)
            if mmer in repeat_dict:
                repeat_dict[mmer].append(seq_id)
            else:
                repeat_dict[mmer] = [seq_id]

    if verbose:
        print

    return repeat_dict

def get_maximal_cliques(repeat_cliques_list):
    # Note: https://goo.gl/2NMKzh
    element_clique_index = defaultdict(set)
    for i,clique in enumerate(repeat_cliques_list):
        for element in clique:
            element_clique_index[element].add(i)
    maximal_cliques_set = set()
    for i, clique in enumerate(repeat_cliques_list):
        maximal_clique = set()
        for element in clique:
            if len(maximal_clique) == 0:
                maximal_clique |= element_clique_index[element]
            else:
                maximal_clique &= element_clique_index[element]
        if len(maximal_clique) == 1: # Current set has no superset
            maximal_cliques_set.add(i)
        else:                        # Current set has a superset
            maximal_cliques_set.update(maximal_clique - set([i]))
    return [repeat_cliques_list[i] for i in maximal_cliques_set]

def get_repeat_cliques(seq_deque, background_list, homology, internal_repeats, verbose):
    repeat_dict = get_repeat_dict(seq_deque, background_list, homology, internal_repeats, verbose)

    # Get unique cliques
    repeat_cliques_set = set()
    while repeat_dict:
        hmer, repeat_clique = repeat_dict.popitem()
        # Note: each repeat_clique is sorted by construction
        repeat_cliques_set.add(tuple(repeat_clique))
    
    # Transform clique tuples into sets
    repeat_cliques_list = list()
    while repeat_cliques_set:
        repeat_cliques_list.append(set(repeat_cliques_set.pop()))

    # Get maximal cliques
    repeat_cliques_list = get_maximal_cliques(repeat_cliques_list)

    # Compute new repeat cliques set
    repeat_cliques_set = set()
    while repeat_cliques_list:
        repeat_cliques_set.add(tuple(repeat_cliques_list.pop()))

    # Stream maximal repeat cliques
    while repeat_cliques_set:
        yield deque(repeat_cliques_set.pop())

def build_homology_graph(cliques, verbose):
    homology_graph = nx.Graph()

    if verbose:
        print
        clear_length = 0

    i = 0
    while True:
        # Stream Clique
        try:
            repeat_clique = cliques.next()
        except:
            break

        # Insert Clique
        if len(repeat_clique) == 1:
            homology_graph.add_node(repeat_clique.popleft())
        else:
            repeat_clique_set = set(repeat_clique)
            while repeat_clique:
                vi = repeat_clique.popleft()
                repeat_clique_set.remove(vi)
                homology_graph.add_node(vi)
                new_neighbors = repeat_clique_set - set(homology_graph[vi])
                for vj in new_neighbors:
                    homology_graph.add_edge(vi, vj)

        # Verbose Output
        if verbose:
            printed = ' [Cliques inserted] = {}\r'.format(i+1)
            sys.stdout.write(' '*clear_length+'\r')
            sys.stdout.write(printed)
            clear_length = len(printed)
            sys.stdout.flush()

        # Update Insertion Count
        i += 1

    if verbose:
        print

    return homology_graph

def get_homology_graph(seq_list, background_list, homology, internal_repeats, seq_file, verbose):
    t0 = time()
    num_seqs = len(seq_list)
    seq_list = uniquify_seq_list(seq_list)
    unq_seqs = len(seq_list)
    if verbose:
        print 'Extracted {} unique sequences out of {} sequences in {:2.4} seconds'.format(unq_seqs, num_seqs, time()-t0)

    t0 = time()
    with open(seq_file, 'w') as outfile:
        seq_deque = deque()
        while seq_list:
            seq_id,seq = seq_list.pop()
            seq_deque.append((seq_id, seq))
            outfile.write('{},{}\n'.format(seq_id, seq))
    if verbose:
        print '\nWritten {} unique sequences out to {} in {:2.4} seconds'.format(unq_seqs, seq_file, time()-t0)

    t0 = time()
    repeat_cliques = get_repeat_cliques(seq_deque, background_list, homology, internal_repeats, verbose)
    homology_graph = build_homology_graph(repeat_cliques, verbose)
    if verbose:
        print '\nBuilt homology graph in {:2.4} seconds. [Edges = {}] [Nodes = {}]'.format(time()-t0, homology_graph.number_of_edges(), homology_graph.number_of_nodes())
        print ' [Intital Nodes = {}] - [Repetitive Nodes = {}] = [Final Nodes = {}]'.format(unq_seqs, unq_seqs - homology_graph.number_of_nodes(), homology_graph.number_of_nodes())
    return homology_graph