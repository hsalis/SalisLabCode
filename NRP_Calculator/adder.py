from glob         import glob
from psutil       import virtual_memory

from sqlitedict   import SqliteDict  as SD
from pybloom_live import BloomFilter as BF

import utils
import finder

import os
import atexit
import math
import shutil

import hyperloglog

current_uuid = None

def setup_proj_dir():
    global current_uuid
    directory = './{}/'.format(current_uuid)
    # print directory
    if not os.path.exists(directory):
        os.makedirs(directory)

def get_non_mem_kmer_db():
    global current_uuid
    with open('./{}/kmer.db'.format(current_uuid), 'w') as tmp_file:
        pass

    return SD('./{}/kmer.db'.format(current_uuid), autocommit=False)

def get_space_capacity(homology):
    return (4 ** homology) / 2

def infinite_parts(target_list):
    for target in target_list:
        if target == float('inf'):
            return True
    return False

def get_background_capacity(homology, background_list, verbose):
    total_kmers_background = hyperloglog.HyperLogLog(0.01) # An HLL object with 1% estimation error
    if background_list:
        if verbose:
            print 'Estimating k-mer Storage Requirements'
        for i,seq in enumerate(background_list):
            if verbose:
                if ((i+1) % 100) == 0:
                    print ' Estimating Background {}/{}'.format(i+1, len(background_list))
            for kmer in utils.stream_min_kmers(seq, k=homology):
                total_kmers_background.add(kmer)
    return int(math.ceil(len(total_kmers_background) * 1.1))

def get_toolbox_capacity(homology, seq_list, target_list):
    if infinite_parts(target_list):
        return float('inf') 
    else:
        total_kmers_parts = 0
        for i, seq in enumerate(seq_list):
            if len(seq) <= homology:
                total_kmers_parts += 1
            else:
                total_kmers_parts += (len(seq) - homology + 1) * target_list[i]
        return total_kmers_parts

def get_system_capacity():
    sys_bits = virtual_memory().available*8.0
    h = 1
    prev_capacity = 0
    cap_padding   = 1024
    while True:
        capacity   = (2**h) + cap_padding
        alloc_bits = int( math.ceil( capacity * math.log(capacity, 2) / ( math.ceil( math.log(capacity, 2) ) * math.log(2) ) ) * math.ceil( math.log(capacity, 2) ) )
        # print h, alloc_bits, sys_bits, alloc_bits < sys_bits
        if alloc_bits < sys_bits:
            h += 1
            prev_capacity = capacity
        else:
            return prev_capacity

def get_capacity(homology, background_list, seq_list, target_list, verbose):
    # Compute relevant capacities
    space_capacity      = get_space_capacity(homology)
    system_capacity     = get_system_capacity()
    background_capacity = get_background_capacity(homology, background_list, verbose)
    toolbox_capacity    = get_toolbox_capacity(homology, seq_list, target_list)
    problem_capacity    = background_capacity + toolbox_capacity
    bf_buffer_space     = 1024
    min_reqd_capacity   = min([space_capacity, problem_capacity, system_capacity])
    capacity            = min_reqd_capacity + bf_buffer_space
    mem                 = True

    # Print out capacity information
    if background_list: print
    print 'Estimated k-mer   Space  Capacity = {} bits'.format(space_capacity)
    print 'Estimated System  Memory Capacity = {} bits'.format(system_capacity)
    print 'Estimated Problem Space  Capacity = {} bits'.format(problem_capacity)
    if problem_capacity+bf_buffer_space > system_capacity:
        print ' >> System memory insufficient for toolbox and/or background ... switching to secondary storage'
        mem = False
    else:
        if problem_capacity > space_capacity:
            print ' >> k-mer Space is smaller than Problem Space ... fewer parts will be generated'
        print 'Allocated BloomFilter    Memory   = {} bits'.format(capacity)

    return mem, capacity

def get_mem_kmer_db(capacity):
    error_rate = 1.0 / capacity
    return BF(capacity=capacity, error_rate=error_rate)

def init_back_db(background_list, homology, seq_list, target_list, verbose):
    if background_list:
        background_list = utils.uniquify_background_list(background_list)

    mem, capacity = get_capacity(homology, background_list, seq_list, target_list, verbose)
    if mem:
        back_db = get_mem_kmer_db(capacity)
    else:
        setup_proj_dir()
        back_db = get_non_mem_kmer_db()
        
    # back_db = set()

    if background_list:
        if verbose:
            print '\nInitializing Background Data Structure'
        for i, seq in enumerate(background_list):
            for kmer in utils.stream_min_kmers(seq, k=homology):
                if mem:
                    back_db.add(kmer)
                else:
                    try:
                        back_db[kmer] = True
                    except:
                        raise MemoryError('Secondary Storage Exhausted.')
            if verbose:
                if ((i+1) % 100) == 0:
                    print ' Background Processed {}/{}'.format(i+1, len(background_list))
            if not mem:
                if ((i+1) % 100) == 0:
                    back_db.commit()
        if not mem:
            back_db.commit()
        if verbose:
            print 'Inserted {} {}-mers in Background'.format(len(back_db), homology)
    return mem, back_db

@atexit.register
def remove_proj_dir():
    global current_uuid
    directory = './{}/'.format(current_uuid)
    if not current_uuid is None and os.path.isdir(directory):
        shutil.rmtree(directory)

def nrp_adder(back_db, sel_list, homology, mem, vercov_func=None, verbose=True):
    if not sel_list:
        return sel_list
    if len(sel_list) > 1:
        nh_sel = finder.nrp_finder(list(sel_list), homology, internal_repeats=False, vercov_func=vercov_func, verbose=verbose)
    else:
        return sel_list
    # print ' nh-sel found = {}'.format(len(nh_sel))
    return [sel_list[i] for i in nh_sel]

def kmerdb_update(back_db, sel_list, homology, mem):
    for seq in sel_list:
        if not seq is None:
            for kmer in utils.stream_min_kmers(seq, k=homology):
                if not mem:
                    try:
                        back_db[kmer] = True
                    except:
                        return False # Secondary storage exhausted
                else:
                    try:
                        back_db.add(kmer)
                    except:
                        return False # Cannot add any more k-mers to BloomFilter
    return True
