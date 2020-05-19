import hgraph
import recovery
import utils

import os
import uuid
import atexit
import shutil

import time

current_uuid = None

def setup_proj_dir():
    global current_uuid
    directory = './{}/'.format(current_uuid)
    if not os.path.exists(directory):
        os.makedirs(directory)

def nrp_finder(seq_list, background_list, homology, internal_repeats=False, vercov_func=None, verbose=True):
    if verbose:
        print '\n[Non-Repetitive Parts Calculator - Finder Mode]\n'
    global current_uuid
    current_uuid = str(uuid.uuid4())
    setup_proj_dir()
    seq_file       = './{}/seq_list.txt'.format(current_uuid)
    graph_file     = './{}/repeat_graph.txt'.format(current_uuid)
    homology_graph = hgraph.get_homology_graph(seq_list, background_list, homology, internal_repeats, seq_file, verbose)
    indi_seqs_indx = recovery.get_recovered_non_homologs(homology_graph, graph_file, vercov_func, verbose)
    indi_seqs      = {}
    with open(seq_file, 'r') as infile:
        for line in infile:
            index,seq = line.strip().split(',')
            index = int(index)
            if index in indi_seqs_indx:
                indi_seqs[index] = seq
    remove_proj_dir()
    if verbose:
        print '\nFinal Non-Repetitive Toolbox size: {}'.format(len(indi_seqs))
    return indi_seqs

@atexit.register
def remove_proj_dir():
    global current_uuid
    directory = './{}/'.format(current_uuid)
    if not current_uuid is None and os.path.isdir(directory):
        shutil.rmtree(directory)

def check_validity(seq_list, homology):
    kmer_dict = {}
    for seq in seq_list:
        for kmer in utils.stream_min_kmers(seq, k=homology):
            if kmer in kmer_dict:
                print kmer
            else:
                kmer_dict[kmer] = seq

def main():
    t0              = time.time()
    homology        = 16
    fasta_filename  = 'riboz.fa' #'input.fa.bk104'
    seq_list        = utils.get_fasta_seq_list(fasta_filename)
    # check_validity(seq_list, homology)
    background_list = []#utils.get_fasta_seq_list(fasta_filename='input.fa')
    non_homologs    = nrp_finder(seq_list, background_list, homology, verbose=True)
    print '{}\t{}'.format(homology-1, len(non_homologs))
    print 'Wall Time {} sec'.format(time.time()-t0)
    check_validity(seq_list=non_homologs.values(), homology=homology)

if __name__ == '__main__':
    main()
