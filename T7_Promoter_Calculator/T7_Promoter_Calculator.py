import os, sys, itertools, pickle, random
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from sklearn.preprocessing import OneHotEncoder
from warnings import simplefilter
from sklearn.exceptions import InconsistentVersionWarning
simplefilter("ignore", category=InconsistentVersionWarning)

def _revcomp(seq):
    comp = {'A' : 'T', 'G' : 'C', 'T' : 'A', 'C' : 'G', 'U' : 'A',
            'W' : 'W', 'S' : 'S', 'M' : 'K', 'K' : 'M', 'R' : 'Y', 'Y' : 'R',
            'B' : 'V', 'V' : 'B', 'D' : 'H', 'H' : 'D', 'N' : 'N',
            }
    x = "".join([comp[letter] for letter in seq ])
    return x[::-1]  #reverse

class T7_PromoterCalculator_Predict(object):

    #V-11544 reference. T7 3.8 promoter from NGS in vitro transcription reads. Reference is ddG = 0.
    BETA = 1.63600
    K = 100000.00000
    kmer_length = 4
    PROMOTER_LENGTH = 39 #bp
    TSS_OFFSET_FROM_START = 27 #bp
    model_filename = 'Ridge_model_1000_thresh_4mers_1bcerr_3maebound.pkl'

    def __init__(self):

        kmer_array = [''.join(p) for p in itertools.product(['A', 'C', 'G', 'T'], repeat=self.kmer_length)]
        kmer_array = np.array(kmer_array).reshape(-1, 1)
        onehot_encoder = OneHotEncoder(sparse_output=False, handle_unknown='ignore')
        onehot_encoder.fit(kmer_array)
        self.kmer_encoder = onehot_encoder
        self.num_features = len(self.kmer_encoder.get_feature_names_out())

        self.directory = os.path.dirname(os.path.abspath(__file__))
        with open(os.path.join(self.directory, self.model_filename), 'rb') as fc:
            self.model = pickle.load(fc)

    def generate_kmers(self, sequence):
        return [sequence[i:i+self.kmer_length] for i in range(len(sequence) - self.kmer_length + 1)]

    def encode_sequences(self, sequenceList, VERBOSE=False):
        encoded_sequences = []
        for sequence in sequenceList:
            kmers = self.generate_kmers(sequence)
            num_kmers = len(kmers)

            # Initialize a matrix to store encoded k-mers
            encoded_matrix = np.zeros((num_kmers, self.num_features))

            for i, kmer in enumerate(kmers):
                encoded_kmer = self.kmer_encoder.transform([[kmer]])
                encoded_matrix[i] = encoded_kmer.flatten()

            # Flatten the matrix and append to the list
            encoded_sequences.append(encoded_matrix.flatten())

            # Print the encoded sequence
            if VERBOSE:
                print(f"Encoded sequence for {sequence}:")
                print(encoded_matrix)
                print("Flattened encoded sequence:")
                print(encoded_matrix.flatten())

        encoded_array = np.atleast_2d(np.array(encoded_sequences))
        if VERBOSE: print(f"Encoded sequences shape: {encoded_array.shape}")
        return encoded_array

    def scan_sequence(self, sequence, TSS_range):
        #Breaks the sequence into 39 bp segments, since the model is trained on 39 bp sequences.
        #pos is the position of the sequence that the model starts from (i.e. it is at -27 from the TSS. The +1nt TSS is expected to sequence[27]

        begin = max(0, TSS_range[0] - self.TSS_OFFSET_FROM_START)
        end = min(TSS_range[1] - self.TSS_OFFSET_FROM_START, len(sequence)-self.PROMOTER_LENGTH + 1)
        seq_scan_list = [(pos+self.TSS_OFFSET_FROM_START, sequence[pos:(pos+self.PROMOTER_LENGTH)]) for pos in range(begin, end)]
        return seq_scan_list

    def run(self, sequence, TSS_range = None, onlyStrand = None):

        seqlen = len(sequence)

        if TSS_range is None: TSS_range = [0, seqlen]

        #Initialize result dictionary
        output = {'beta' : self.BETA,
                  'k' : self.K,
                  'sequence' : sequence,
                  'TSS_range' : TSS_range,
                  'Forward_Predictions_per_TSS' : {}, 'Reverse_Predictions_per_TSS' : {},
                 }

        # Forward Strand
        if onlyStrand is None or onlyStrand == '+':
            TSS_list = []
            subsequence_list = []

            for (TSS, subsequence) in self.scan_sequence(sequence, TSS_range):
                TSS_list.append(TSS)
                subsequence_list.append(subsequence)

            if len(subsequence_list) > 0:
                encoded_subsequences = self.encode_sequences(subsequence_list)
                ddG_predicted_list = [float(x) for x in self.model.predict(encoded_subsequences)]
                TX_predicted_list = [ float(self.K * np.exp(-self.BETA * ddG_predicted)) for ddG_predicted in ddG_predicted_list ]

                for (TSS, subsequence, ddG_predicted, TX_predicted) in zip(TSS_list, subsequence_list, ddG_predicted_list, TX_predicted_list):
                    output['Forward_Predictions_per_TSS'][TSS] = {'Tx_rate' : TX_predicted,
                                                                  'promoter_sequence' : subsequence,
                                                                  'dG_total' : ddG_predicted,
                                                                  }
        # Reverse Strand
        if onlyStrand is None or onlyStrand == '-':
            rc_sequence = _revcomp(sequence)
            rc_TSS_range = [seqlen - TSS_range[1], seqlen - TSS_range[0]]

            TSS_list = []
            subsequence_list = []
            for (TSS, subsequence) in self.scan_sequence(rc_sequence, rc_TSS_range):
                TSS_list.append(TSS)
                subsequence_list.append(subsequence)

            if len(subsequence_list) > 0:
                encoded_subsequences = self.encode_sequences(subsequence_list)
                ddG_predicted_list = [float(x) for x in self.model.predict(encoded_subsequences)]
                TX_predicted_list = [ float(self.K * np.exp(-self.BETA * ddG_predicted)) for ddG_predicted in ddG_predicted_list ]

                for (TSS, subsequence, ddG_predicted, TX_predicted) in zip(TSS_list, subsequence_list, ddG_predicted_list, TX_predicted_list):

                    output['Reverse_Predictions_per_TSS'][seqlen - TSS] = {'Tx_rate' : TX_predicted,
                                                                           'promoter_sequence' : subsequence,
                                                                           'dG_total' : ddG_predicted,
                                                                           }
        return output

if __name__ == '__main__':

    test_sequence = 'CGAAGACGGTGGTGTTGTTGTCCCTAAATTAATACGACTCACTATAGGGAGATAGGGGTTTCACAGATTATTCATTCATACTACCCCATTACTGCCTACCGCACCGCCACCATGGCGCTCCTGCTGTGCTTCGTGCTCCTGTGCGGAGTAGTGGATTTCGCCAGAAGTTTGAGTATCACTACTCCTGAAGAGATGATTGAAAAAGCCAAAGGGGAAACTGCCTATCTGCCATGCAAATTTACGCTTAGTCCCGAAGACCAGGGACCGCTGGACATCGAGTGGCTGATATCACCAGCTGATAATCAGAAGGTGGATCAAGTGATTATTTTATATTCTGGAGACAAAATTTATGATGACTACTATCCAGATCTGAAAGGCCGAGTACATTTTACGAGTAATGATCTCAAATCTGGTGATGCATCAATAAATGTAACGAATTTACAACTGTCAGATATTGGCACATATCAGTGCAAAGTGAAAAAAGCTCCTGGTGTTGCAAATAAGAAGATTCATCTGGTAGTTCTTGTTAAGCCTTCAGGTGCGAGATGTTACGTTGATGGATCTGAAGAAATTGGAAGTGACTTTAAGATAAAATGTGAACCAAAAGAAGGTTCACTTCCATTACAGTATGAGTGGCAAAAATTGTCTGACTCACAGAAAATGCCCACTTCATGGTTAGCAGAAATGACTTCATCTGTTATATCTGTAAAAAATGCCTCTTCTGAGTACTCTGGGACATACAGCTGTACAGTCAGAAACAGAGTGGGCTCTGATCAGTGCCTGTTGCGTCTAAACGTTGTCCCTCCTTCAAATAAAGCTGGACTAATTGCAGGAGCCATTATAGGAACTTTGCTTGCTCTAGCGCTCATTGGTCTTATCATCTTTTGCTGTCGTAAAAAGCGCAGAGAAGAAAAATATGAAAAGGAAGTTCATCACGATATCAGGGAAGATGTGCCACCTCCAAAGAGCCGTACGTCCACTGCCAGAAGCTACATCGGCAGTAATCATTCATCCCTGGGGTCCATGTCTCCTTCCAACATGGAAGGATATTCCAAGACTCAGTATAACCAAGTACCAAGTGAAGACTTTGAACGCACTCCTCAGAGTCCGACTCTCCCACCTGCTAAGGTAGCTGCCCCTAATCTAAGTCGAATGGGTGCGATTCCTGTGATGATTCCAGCACAGAGCAAGGATGGGTCTATAGTATAG'

    model = T7_PromoterCalculator_Predict()
    output = model.run(test_sequence)

    for (TSS, data) in output['Forward_Predictions_per_TSS'].items():
        print('fwd TSS: {}. dG_total: {}. TX rate: {}. Sequence: {}'.format(TSS, data['dG_total'], data['Tx_rate'], data['promoter_sequence']))
    print('Max fwd TIR: ', max([data['Tx_rate'] for data in output['Forward_Predictions_per_TSS'].values()]))
    for (TSS, data) in output['Reverse_Predictions_per_TSS'].items():
        print('rev TSS: {}. dG_total: {}. TX rate: {}. Sequence: {}'.format(TSS, data['dG_total'], data['Tx_rate'], data['promoter_sequence']))
    print('Max rev TIR: ', max([data['Tx_rate'] for data in output['Reverse_Predictions_per_TSS'].values()]))
