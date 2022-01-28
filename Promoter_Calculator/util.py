from sklearn.preprocessing import OneHotEncoder, LabelEncoder
from Bio.Seq               import Seq
from scipy                 import stats
import numpy as np

import itertools
import collections
import operator

groove_access = {'CG' : 43,
                 'CA':  42, 'TG' : 42,  #revcomp
                 'GG':  42, 'CC' : 42,  #revcomp
                 'GC':  25,
                 'GA':  22, 'TC' : 22,  #revcomp
                 'TA':  14,
                 'AG':  9, 'CT' : 9,    #revcomp
                 'AA':  5, 'TT' : 5,    #revcomp
                 'AC':  4, 'GT' : 4,    #revcomp
                 'AT':  0
                }

stacking_dict = {'AA' :-2.41,
                 'AC':-2.06, 'CA' : -2.06, #reverse
                 'AG':-2.54, 'GA' : -2.54, #reverse
                 'AT':-2.25, 'TA' : -2.25, #reverse
                 'CC':-1.93,
                 'CG':-2.24, 'GC' : -2.24, #reverse
                 'CT':-2.01, 'TC' : -2.01, #reverse
                 'GG':-2.66,
                 'GT':-2.46, 'TG' : -2.46, #reverse
                 'TT':-1.82
                 }

average_twist = {'AA' : 35.6, 'TT' : 35.6, #revcomp
                 'AT':  32.1, 'AT' : 32.1, #revcomp
                 'TA':  35.3, 'TA' : 35.3, #revcomp
                 'GG':  33.65, 'CC' : 33.65, #revcomp
                 'GC':  40.2,
                 'CG':  29.9,
                 'AG':  27.9, 'CT' : 27.9, #revcomp
                 'GA':  36.8, 'TC' : 36.8, #revcomp
                 'AC':  34.0, 'GT' : 34.0, #revcomp
                 'CA':  34.4, 'TG' : 34.4 #revcomp
                 }

# Sequence dependence of DNA bending rigidity, Geggier and Vologodskii
persistence   = {'AA' : 50.4, 'TT' : 50.4, #revcomp
                 'AC':  55.4, 'GT' : 55.4, #revcomp
                 'AG':  51.0, 'CT' : 51.0, #revcomp
                 'AT':  40.9, 'AT' : 40.9, #revcomp
                 'CA':  46.7, 'TG' : 46.7, #revcomp
                 'CC':  41.7, 'GG' : 41.7, #revcomp
                 'CG':  56.0,
                 'GA':  54.4, 'TC' : 54.4, #revcomp
                 'GC':  44.6,
                 'TA':  44.7,
                }
RNA_DNA_hybrids = {'TT': -1.0,
                   'TG': -2.1,
                   'TC': -1.8,
                   'TA': -0.9,
                   'GT': -0.9,
                   'GG': -2.1,
                   'GC': -1.7,
                   'GA': -0.9,
                   'CT': -1.3,
                   'CG': -2.7,
                   'CC': -2.9,
                   'CA': -1.1,
                   'AT': -0.6,
                   'AG': -1.5,
                   'AC': -1.6,
                   'AA': -0.2
                   }

DNA_DNA_hybrids = { 'AA':   -1.00, 'TT' : -1.00, #revcomp
                    'AT':   -0.88,
                    'TA':   -0.58,
                    'CA':   -1.45, 'TG' : -1.45, #revcomp
                    'GT':   -1.44, 'AC' : -1.44, #revcomp
                    'CT':   -1.28, 'AG' : -1.28, #revcomp
                    'GA':   -1.30, 'TC' : -1.30, #revcomp
                    'CG':   -2.17,
                    'GC':   -2.24,
                    'GG':   -1.42, 'CC' : -1.42  #revcomp
                  }


### This function calcs the hamming dist between 2 sequences of = len
def hamming(seq1, seq2):
    distance = 0
    L = len(seq1)
    if L == len(seq2):
        for i in range(L):
            if seq1[i] != seq2[i]:
                distance += 1
    else:
        print('Sequences are not of equal length')
        print(seq1, seq2)
        distance = None
    return distance

def get_rsquared(TXpred,TXempirical):
    slope, intercept, r_value, p_value, std_err = stats.mstats.linregress(TXpred,TXempirical)
    return r_value**2

### RETURN THE REVERSE COMPLIMENT OF A SEQ
revcomp_letters_util = {'U' : 'A', 'A' : 'T', 'G' : 'C', 'T' : 'A', 'C' : 'G'}
def reverse_complement(seq):
    return "".join([revcomp_letters_util[letter] for letter in seq[::-1] ])
    ## seq = Seq(seq)
    ## return seq.reverse_complement()

def AT_content(seq):
    ATcont = (seq.count('A') + seq.count('T'))/float(len(seq))*100
    return ATcont

def get_dg_matrices(parameters, encoder):
    dg_dict = {}
    mer_lis = []
    for entry in encoder.categories_[0]:
        mer_lis.append(entry)
    reference = mer_lis
    my_dict = dict(zip(reference,parameters))
    my_dict2 = collections.OrderedDict(sorted(my_dict.items()))
    sorted_by_val = sorted(my_dict2.items(), key=operator.itemgetter(1))
    for entry in sorted_by_val:
        dg_dict[entry[0]] = entry[1]
    return dg_dict

# Functions to calculate DNA or RNA properties
def calc_stacking_and_twist_energy(seq):
    #find stacking free enegy of spacer
    spacer_sfe  = 0
    spacer_twist = 0
    for j in range(0,len(seq)-1,2):
        mer = seq[j:j+2]
        spacer_sfe += stacking_dict[mer]
        spacer_twist += average_twist[mer]
        # try forward and reverse or compliment
        ## try:
        ##     spacer_sfe += stacking_dict[mer]
        ## except:
        ##     spacer_sfe += stacking_dict[mer[::-1]]
        ## try:
        ##     spacer_twist += average_twist[mer]
        ## except:
        ##    spacer_twist += average_twist[reverse_complement(mer)]

    return spacer_sfe, spacer_twist

# local DNA element rigidity
def calc_rigidity(seq):
    rigidity = 0
    for m in range(0, len(seq), 2):
        mer = seq[m:m+2]
        rigidity += persistence[mer]
        ## try:
        ##     rigidity += persistence[mer]
        ## except:
        ##    rigidity += persistence[reverse_complement(mer)]

    rigidity = rigidity / len(seq)
    return rigidity

def calc_disc_melting(seq, downstream):
    dg_dna = 0
    for j in range(0,len(seq[0:len(seq)])-1,2):
        mer = seq[j:j+2]
        if len(mer) < 2:    mer = mer + downstream[0]
        if mer in DNA_DNA_hybrids:  dg_dna += DNA_DNA_hybrids[mer]
        dg_dna += DNA_DNA_hybrids[mer]
        ## try:
        ##     dg_dna += DNA_DNA_hybrids[mer]
        ## except:
        ##    dg_dna += DNA_DNA_hybrids[reverse_complement(mer)]

    return dg_dna

# find free enegy of DNA:DNA and RNA:DNA hybrids in ITR
def calc_DNA_RNA_hybrid_energy(seq):
    dg_dna = 0
    dg_rna = 0
    for j in range(0,len(seq[0:15])-1,2): # 15 is the length of long abortive transcripts
        mer = seq[j:j+2]
        dg_dna += DNA_DNA_hybrids[mer]
        # try forward and reverse or compliment
        ## try:
        ##     dg_dna += DNA_DNA_hybrids[mer]
        ## except:
        ##    dg_dna += DNA_DNA_hybrids[reverse_complement(mer)]
        dg_rna += RNA_DNA_hybrids[mer]

    dg_hybrid    = dg_dna - dg_rna
    return dg_dna,dg_rna,dg_hybrid

def calc_groove_width(seq):
    groove_width = 0
    for l in range(0, len(seq), 2):
        mer1 = seq[l:l+2]
        groove_width += groove_access[mer1]
        ## try:
        ##    groove_width+= groove_access[mer1]
        ## except:
        ##    groove_width+= groove_access[reverse_complement(mer1)]
    return groove_width

def length_encoders(lower_bound, upper_bound):
    length_array =[]
    for x in range(lower_bound, upper_bound+1):
        length_array.append(str(x))
    length_array = np.array(length_array)

    #one hot encode the array of length into bit vectors
    onehot_encoder = OneHotEncoder(sparse=False)
    onehot_encoder.fit(length_array.reshape(-1, 1))
    return onehot_encoder

def kmer_encoders(k):
    # k is the kmer size for DNA sequence
    # this will give bit vectors for 2mers or 3mers for us
    kmer_array =[]
    nts = ['A','T','G','C']
    for output in itertools.product(nts, repeat = k):
        kmer_array.append(''.join(output))
    kmer_array = np.array(kmer_array)

    #one hot encode the sequence
    onehot_encoder = OneHotEncoder(sparse=False)
    onehot_encoder.fit(kmer_array.reshape(-1, 1))
    return onehot_encoder