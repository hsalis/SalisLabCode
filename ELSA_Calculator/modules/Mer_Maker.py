from Bio import SeqIO
from time import time
import re
import pickle
import random
import sys
import os

def mers(length): 
	"""Generates multimers for sorting through list of 10mers based on user 
	specification. Multimers generated act as the keys for generating a 
	hashtable to eliminate undesired sequence patterns from those 10mers not
	found in the genome. 
	
	Usage: mers(N) = 4^(N) unique Nmers 
	"""
	seq_list = ['']
	counter=0
	while counter < length:
		for seq in seq_list:
			if len(seq) == counter:
				for x in ['A', 'T', 'C', 'G']:
					seq_list.append(seq+x)
		counter += 1
	last_N_Mers = 4**length
	return seq_list[len(seq_list)-last_N_Mers :]
