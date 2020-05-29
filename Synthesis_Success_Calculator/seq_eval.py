import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), "../.."))  
from PyVRNA import PyVRNA
from Bio.Seq import Seq
from Bio import SeqIO
import re
import numpy as np
from bisect import bisect
import collections
from subprocess import Popen, PIPE, STDOUT
import cStringIO
import scipy.io as sio
from Bio.SeqUtils import MeltingTemp as TM
from FastFinder import FastFinder
import itertools
import time 
from string import maketrans
import json
from collections import OrderedDict, defaultdict
from sklearn.ensemble import RandomForestClassifier as RFC
import pickle as pkl
import pandas as pd
from treeinterpreter import treeinterpreter as ti
import random

class seq_assess(object):
    def __init__(self, sequence = None, begin=0, verbose = False):
        # defaults
        self.rc = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
        self.bad_seqs=[]
        self.score=[0,0,0,0]
        self.comp_table = maketrans('ATGC', 'TACG')
        self.verbose = verbose 
        self.begin = begin
        if sequence:
            self.load(sequence)

    def load(self,sequence):
        #input validation
        try:
            if 'U' in sequence.seq or 'u' in sequence.seq:
                print sequence.seq
                sequence=sequence.back_transcribe()
            sequence1 = str(sequence.seq).upper()  
            sequence2 = str(sequence.seq.reverse_complement()).upper()

        except:
            try:
                if 'U' in sequence or 'u' in sequence:
                    print sequence
                    sequence=sequence.back_transcribe()
                sequence1 = str(sequence).upper() 
                sequence2 = str(sequence.reverse_complement()).upper()
            except:
                try:
                     sequence1 = str(sequence)
                     sequence1 = sequence1.upper()
                     #print sequence1
                     bases = list(sequence1) 
                     bases = reversed([self.rc.get(base,base) for base in bases])
                     sequence2 = ''.join(bases)
                     #print sequence2
                except:
                    print "The sequence provided was neither a Biopython record, Biopython Seq class, nor a valid string with your DNA sequence. Check what you are passing to seq_assess."
                    self.bad_seqs.append("The sequence provided was neither a Biopython record, Biopython Seq class, nor a valid string with your DNA sequence. Check what you are passing to seq_assess.")
        # setting up internal sequence
        self.size=len(sequence)
        self.sequence=sequence1.upper()
        self.sequence1=sequence1.upper()
        self.sequence2=sequence2.upper()

    def run(self):
        ### quickmode for running as an evaluate method
        self.s_info_2={}
        self.s_info_2['size']=[]
        #t0=time.time()
        self.find_repeats(RL=9,Rmax=11)
        #t1=time.time()
        self.find_hairpins_windowed()
        #self.fast_hairpins()
        #t2=time.time()
        self.find_nucleotides()
        if self.verbose:
            for i in self.bad_seqs:
                print i

        self.r_metrics = OrderedDict([("size", self.size),
                               ("terminal_repeat", self.terminal_repeat),
                                        ('tandem', self.tandem), # number of tandem repeats
                                   ("freq_repeat", self.freq_repeat),  # count of the most frequent repeat
                                ("longest_repeat", self.longest_repeat), # length of the longest repeat
                          ("high_local_density_1", self.high_local_density_1), # number of windows that exceed local density screen (0.90 of sequence is in some way repetitive for window of 70bp)
                          ("high_local_density_2", self.high_local_density_2),  # number of windows that exceed local density screen (0.60 of sequence is in some way repetitive for window of 500bp)
                         ("high_specific_density", self.high_specific_density), # number of windows that exceed specific density screen (0.40 of sequence is repetitive to a given repeat)
                            ("high_total_density", self.high_total_density),
                                     ("repeat_10", self.repeat_10),
                                     ("repeat_15", self.repeat_15),
                                     ("repeat_20", self.repeat_20),         
                                     ("repeat_25", self.repeat_25),
                                     ("repeat_40", self.repeat_40), 
                                  ("repeat_large", self.repeat_large),
                                  ("repeat_9_metric", self.repeat9)]) # whether sequence exceeds repeat screen (0.69 of sequence is repetitive in some way)}

        self.h_metrics = OrderedDict([     ("long_hairpins", self.long_hairpins), # hairpins (type1)
                             ("strong_hairpins", self.strong_hairpins), # hairpins (type2)
                               ("wide_hairpins", self.wide_hairpins), # hairpins (type3)
                             ("longest_hairpin", self.longest_hairpin),
                             ("richest_hairpin", self.richest_hairpin), # GC content of GC richest hairpin
                                 ("palindromes", self.palindromes), # palindromes}
                           ("terminal_hairpins", self.terminal_hairpins)])
                          #"g_quad_predictions", self.g_quad_predictions}

        self.n_metrics = OrderedDict([   ("poly_runs", self.poly_runs), 
                         ("pattern_runs", self.pattern_runs),
                             ("i_motifs", self.i_motifs),
                        ("g_quad_motifs", self.g_quadruplexes), 
                           ("GC_short_h", self.GC_short_h), 
                           ("GC_short_l", self.GC_short_l), 
                            ("GC_long_h", self.GC_long_h), 
                            ("GC_long_l", self.GC_long_l), 
                            ("GC_term_h", self.GC_term_h), 
                            ("GC_term_l", self.GC_term_l), 
                                  ("dGC", self.dGC), 
                               ("Tm_low", self.Tm_low), 
                              ("Tm_high", self.Tm_high),
                                  ("dTm", self.dTm),
                             ("total_GC", self.total_GC)])

        self.s_count = self.r_count + self.h_count + self.n_count
        self.s_loc_list = self.r_loc_list + self.h_loc_list+self.n_loc_list 

        self.s_info = {}
        for info in [self.r_info,self.h_info,self.n_info]:
            for key in info.keys():
                
                try:
                    x = len(info[key]["locations"][0]) #this should be a tuple and have length, but sometimes it's not!
                except:
                    print info[key]
            
                try:
                    for location in info[key]["locations"]:
                        if location[0] > location[1]:
                            print "location is ", location
                            print "The begin position must be less than the end position."
                            print info[key]
                except:
                    print "Error"
                    print info[key]

                if key not in self.s_info:
                    self.s_info[key] = info[key]
                else:
                    self.s_info[key]["locations"] += info[key]["locations"]
                    self.s_info[key]["types"] += info[key]["types"]

    def find_repeats(self, RL = 8,
                           Rmax=11, 
                           Rthreshold = [0.69,0.40,0.90,0.60], 
                           Rwindow = [70,500,60], 
                           Rtandem = 5):
        ### code for testing the following six rules:
        # A) Any sequence that contains more that 69 % of  bases that are part of repeats 8 bases or longer *   "X= 69 Y=8" 10
        # B) Any single repeated sequence contains more that 40 % of  bases that are part of repeats 8 bases or longer *    "X=40 Y=8"  6
        # C) Any sequence that contains more that 90 % of  bases that are part of repeats 8 bases or longer within a window of 70 length *  "X=90 Y=8 Z=70" 10
        # D) Any sequence that contains more that 60 % of  bases that are part of repeats 8 bases or longer within a window of 500 length for sequences longer than 1000 *  "X=60 Y=8 Z=500 L=1000" 10
        # E) a tandem repeat of length 5 or greater beginning or ending at either termini of the sequence   X=5 10
        # F) repeats more than 12 bp long risk impede synthesis
        
       
        ## outputs for Operon Calculator
        self.r_count = 0
        self.r_loc_list = [] 
        self.r_info = {}

        # metrics for tracking repeats for ML analysis
        self.tandem = 0 # number of tandem repeats
        self.s_info_2['tandem']=[]
        self.freq_repeat = 0 # count of the most frequent repeat
        self.s_info_2['freq_repeat']=[]
        self.longest_repeat = 0 # length of the longest repeat
        self.s_info_2['longest_repeat']=[]
        self.high_local_density_1 = 0 # number of windows that exceed local density screen (0.90 of sequence is in some way repetitive for window of 70bp)
        self.s_info_2['high_local_density_1']=[]
        self.high_local_density_2 = 0 # number of windows that exceed local density screen (0.60 of sequence is in some way repetitive for window of 500bp)
        self.s_info_2['high_local_density_2']=[]
        self.high_specific_density = 0 # number of windows that exceed specific density screen (0.40 of sequence is repetitive to a given repeat)
        self.s_info_2['high_specific_density']=[]
        self.high_total_density = 0 # whether sequence exceeds repeat screen (0.69 of sequence is repetitive in some way)
        self.s_info_2['high_total_density']=[]
        self.repeat_10 = 0
        self.s_info_2['repeat_10']=[]
        self.repeat_15 = 0
        self.s_info_2['repeat_15']=[]
        self.repeat_20 = 0
        self.s_info_2['repeat_20']=[]
        self.repeat_25 = 0
        self.s_info_2['repeat_25']=[]
        self.repeat_40 = 0 
        self.s_info_2['repeat_40']=[]
        self.repeat_large = 0
        self.s_info_2['repeat_large']=[]
        self.terminal_repeat = 0
        self.s_info_2['terminal_repeat']=[]
        self.repeat9 = 0
        self.s_info_2['repeat_9_metric']=[]
        ## FastFinder, courtesy of Ayaan Hossain
        RP=FastFinder()

        self.RepeatDict=RP.get_repeat_dict([self.sequence1], RL, verbose = False)

        total=np.zeros(self.size, dtype=bool)
        total2=np.zeros(self.size)
        r_distribution={} # for gathering distributions of repeats for visualizations

        repeat9_numerator = 0

        location={}

        for i in self.RepeatDict.keys():
            self.r_info[i]={"locations":[],"types":[]}
            val=np.zeros(self.size, dtype=bool)

            if len(i) > Rmax: # F)
                self.bad_seqs.append("This repeat is too long: %s"%(i))
                if len(i)>=self.longest_repeat:
                    self.longest_repeat = len(i)
                    for j in self.RepeatDict[i].keys():
                        location_list = list(self.RepeatDict[i][j][0])
                        self.s_info_2['longest_repeat'].extend([((self.begin + k, self.begin + k+len(i)),1.0) for k in location_list])
                    
            for j in self.RepeatDict[i].keys():
                location_list = list(self.RepeatDict[i][j][0])
                self.r_info[i]["types"].extend([j])
                self.r_info[i]["locations"].extend([(self.begin + k, self.begin + k+len(i)) for k in location_list])
                if len(location_list)>self.freq_repeat:
                    self.freq_repeat=len(location_list)
                    self.s_info_2['freq_repeat']=[((self.begin + k, self.begin + k+len(i)),1.0) for k in location_list]
                for k in location_list:
                    if len(i) > Rmax:
                        self.score[0]+=len(i)
                    if len(i) > 8:
                        for l in range(8,len(i)):
                            repeat9_numerator +=9
                    if (self.begin + k)<60 or (self.begin + k) >=self.size-60:
                        self.terminal_repeat+=1
                        self.s_info_2['terminal_repeat'].extend([((self.begin + k, self.begin + k+len(i)),1.0) for k in location_list])

                    self.r_count+=1
                    self.r_loc_list.extend([( (self.begin + k, self.begin + k+len(i)), 1.0)])

                    for v in range(k,k+len(i)):
                        val[v]=True

                    if len(i) not in r_distribution:
                        #for l in range(RL,len(i)):
                        #    if l not in r_distribution:
                        #        r_distribution[l] = 0
                        r_distribution[len(i)] = 1
                    else:
                        r_distribution[len(i)]+=1
                if len(i)<=10:
                    self.s_info_2['repeat_10'].extend([((self.begin + k, self.begin + k+len(i)),1.0) for k in location_list])
                elif 10<len(i) and len(i)<=15:
                    self.s_info_2['repeat_15'].extend([((self.begin + k, self.begin + k+len(i)),1.0) for k in location_list])
                elif 15<len(i) and len(i)<=20:
                    self.s_info_2['repeat_20'].extend([((self.begin + k, self.begin + k+len(i)),1.0) for k in location_list])
                elif 20<len(i) and len(i)<=25:
                    self.s_info_2['repeat_25'].extend([((self.begin + k, self.begin + k+len(i)),1.0) for k in location_list])
                elif 25<len(i) and len(i)<=40:
                    self.s_info_2['repeat_40'].extend([((self.begin + k, self.begin + k+len(i)),1.0) for k in location_list])
                else:
                    self.s_info_2['repeat_large'].extend([((self.begin + k, self.begin + k+len(i)),1.0) for k in location_list])

    
            ind_fract=float(sum(val.astype(int)))/self.size
            if ind_fract> Rthreshold[1]: # B)
                self.bad_seqs.append("Too much of this repeat: %s. Redesign to remove "%(i))
                self.score[0]+=6
                self.high_specific_density+=1
                for j in self.RepeatDict[i].keys():
                    location_list = list(self.RepeatDict[i][j][0])
                    self.s_info_2['high_specific_density'].extend([((self.begin + k, self.begin + k+len(i)),1.0) for k in location_list])
            total+=val
            total2+=val.astype(int)

        total_fract=float(sum(total.astype(int)))/self.size
        if total_fract>Rthreshold[0]: # A)
            self.bad_seqs.append("Too many repeats. Redesign to remove ")
            self.score[0]+=10
            self.high_total_density+=1
            self.s_info_2['high_total_density']=self.r_loc_list

        self.repeat9 = repeat9_numerator / float(self.size)
        self.s_info_2['repeat_9_metric'] = self.r_loc_list

        windowed_repeats=[] # for more visuals
        for i in range(len(Rwindow)-1): # C) & D)

            for j in range(self.size-Rwindow[i]):
                peek = float(sum(total[j:j+Rwindow[i]].astype(int)))

                if peek/Rwindow[i] > Rthreshold[2+i]: 
                    #self.bad_seqs.append("%i bp window starting at %i is too repeat rich. Redesign to remove "%(Rwindow[i],j))
                    self.score[0]+=10
                    if i==0:
                        self.high_local_density_1+=1
                        self.s_info_2['high_local_density_1'].extend([((self.begin + j, self.begin + j+Rwindow[i]),1.0)])

                    if i==1:
                        self.high_local_density_2+=1
                        self.s_info_2['high_local_density_2'].extend([((self.begin + j, self.begin + j+Rwindow[i]),1.0)])
                if i==0:
                    windowed_repeats.append(peek/Rwindow[i])
        # E)
        term1=self.sequence1[:Rwindow[2]]
        self.tandem_subroutine(term1,Rtandem,"five")
        term2=self.sequence1[len(self.sequence1)-Rwindow[2]:]
        self.tandem_subroutine(term2,Rtandem,"three")

        self.repeat_density=windowed_repeats
        self.repeat_distribution=r_distribution
        for i in range(RL,self.longest_repeat):
            if i in r_distribution:
                if i<=10:
                    self.repeat_10+=r_distribution[i]
                elif 10<i and i<=15:
                    self.repeat_15+=r_distribution[i]
                elif 15<i and i<=20:
                    self.repeat_20+=r_distribution[i]
                elif 20<i and i<=25:
                    self.repeat_25+=r_distribution[i]
                elif 25<i and i<=40:
                    self.repeat_40+=r_distribution[i]
                else:
                    self.repeat_large+=r_distribution[i]

    def tandem_subroutine(self, term, Rtandem, end):
        # subroutine that finds tandem repeats using regex
       
        run =self.get_tandem_dict(term, Rtandem, 5, 0)
        if len(run)>0:
            for i in run.keys():
                for j in range(len(run[i])):
                    tandrep=term[i:run[i][j]+Rtandem]
                    self.r_count+=1
                    self.tandem+=1
                    self.r_info[tandrep]={"locations":[],"types":[]}

                    if end == "five":
                        self.r_info[tandrep]["locations"].append( (self.begin + i, self.begin + i + len(tandrep)) ) #Adding offset (begin) if sequence is a subsequence
                        self.r_info[tandrep]["types"].append("tandem")
                        self.r_loc_list += [( (self.begin + i, self.begin + i + len(tandrep)), 1.0)] #Adding offset (begin) if sequence is a subsequence
                        self.s_info_2['tandem'].extend([((self.begin + i, self.begin + i + len(tandrep)),1.0)])
                        self.bad_seqs.append("tandem repeat too close to 5' end: %s"%(tandrep))

                    elif end == "three":
                        self.r_info[tandrep]["locations"].append( (self.begin + self.size - 60 + i, self.begin + self.size - 60 + i+len(tandrep) ) ) #Adding offset (begin) if sequence is a subsequence
                        self.r_info[tandrep]["types"].append("tandem")
                        self.r_loc_list += [( (self.begin + self.size - 60 + i, self.begin + self.size - 60 + i+len(tandrep)), 1.0)] #Adding offset (begin) if sequence is a subsequence
                        self.s_info_2['tandem'].extend([((self.begin + i, self.begin + i + len(tandrep)),1.0)])
                        self.bad_seqs.append("tandem repeat too close to 3' end: %s"%(tandrep))
                self.score[0]+=10

    def is_tandem(self, seq, start, gap, rulen, max_x):
        x = 0
        for i in xrange(rulen):
            if seq[start+i] != seq[start+i+gap+rulen]:
                x += 1
            if x > max_x:
                return False
        return True

    def get_tandem_dict(self, seq, rulen, max_g, max_x):
        tr_chains = defaultdict(list)
        depth     = len(seq)-2*rulen+1
        for start in xrange(depth):
            for gap in xrange(min(depth-start, max_g+1)):
                if self.is_tandem(seq, start, gap, rulen, max_x):
                    tr_chains[start].append(start+rulen+gap)
        return tr_chains

    def hdist(self, seq, start, gap, rulen, max_x):
        x = 0
        for i in xrange(rulen):
            #print len(seq)
            #print start+i
            #print start+i+gap+rulen
            if seq[start+i] != seq[start+i+gap+rulen]:
                x += 1
            if x > max_x:
                return False
        self.tandem+=1
        return True

    def get_tandem_dict(self, seq, rulen, max_g, max_x):
        tr_chains = defaultdict(list)
        for i in xrange(len(seq) - rulen -rulen - 1):
            for g in xrange(min(len(seq)-i-rulen-max_g,max_g+1)):
                if self.hdist(seq, i, g, rulen, max_x):
                    tr_chains[i].append(i+rulen+g)
        return tr_chains

    def hairpin_processing(self, index, lengths, details, gquads, l_threshold, gc_threshold, mfe, rc=False):
        for i, l in enumerate(lengths):
            if l>l_threshold and details[i][3]<mfe and l_threshold<16:
                hpinx = details[i][0]
                structure = details[i][1]
                gc = self.GC_count(hpinx,len(hpinx))
                if rc:
                    location = self.size -(index + details[i][2] + len(hpinx))
                else:
                    location = index + details[i][2]
                location_tup = (self.begin + location, self.begin + location+len(hpinx) )
                if l>self.longest_hairpin:
                    self.longest_hairpin=l
                    self.s_info_2['longest_hairpin']=[((self.begin + location, self.begin + location+len(hpinx) ),1.0)]
                if gc>gc_threshold:
                    if hpinx not in self.h_info.keys():
                        if any(i<60 for i in location_tup)  or any(i>(self.size-60) for i in location_tup):
                            self.terminal_hairpins+=1
                            self.s_info_2['terminal_hairpins'].extend([((self.begin + location, self.begin + location+len(hpinx) ),1.0)])
                        self.h_info[hpinx] = {"locations":[(self.begin + location, self.begin + location+len(hpinx) )],"types":["strong hairpin"]} #Adding offset (begin) if sequence is a subsequence
                        self.h_loc_list += [((self.begin + location, self.begin + location+len(hpinx) ), 1.0)] #Adding offset (begin) if sequence is a subsequence
                        self.h_count +=1
                        self.strong_hairpins+=1
                        self.s_info_2['strong_hairpins'].extend([((self.begin + location, self.begin + location+len(hpinx) ),1.0)])
                        self.hairpins[hpinx]=(self.begin+location,structure)
                        self.score[1]+=10
                        #print "too long hairpin at %i: %s, %s"%(location,hpinx,structure)
                        self.bad_seqs.append("strong hairpin at %i: %s, %s"%(location,hpinx,structure))
                        if gc>self.richest_hairpin:
                                self.richest_hairpin=gc
                                self.s_info_2['richest_hairpin']=[((self.begin + location, self.begin + location+len(hpinx)),1.0)]
                        if gc==self.richest_hairpin:
                            self.s_info_2['richest_hairpin'].extend([((self.begin + location, self.begin + location+len(hpinx) ),1.0)])
                    else:
                        continue
            elif l>l_threshold and details[i][3]<mfe and l_threshold>16 and l_threshold<20:
                hpinx = details[i][0]
                structure = details[i][1]
                gc = self.GC_count(hpinx,len(hpinx))
                if rc:
                    location = self.size -(index + details[i][2] + len(hpinx))
                else:
                    location = index + details[i][2]
                location_tup = (self.begin + location, self.begin + location+len(hpinx) )
                if hpinx not in self.h_info.keys():
                    if any(i<60 for i in location_tup)  or any(i>(self.size-60) for i in location_tup):
                        self.terminal_hairpins+=1
                    self.h_info[hpinx] = {"locations":[(self.begin + location, self.begin + location+len(hpinx) )],"types":["long hairpin"]} #Adding offset (begin) if sequence is a subsequence
                    self.h_loc_list += [((self.begin + location, self.begin + location+len(hpinx) ), 1.0)] #Adding offset (begin) if sequence is a subsequence
                    self.h_count +=1
                    self.wide_hairpins+=1
                    self.s_info_2['wide_hairpins'].extend([((self.begin + location, self.begin + location+len(hpinx) ),1.0)])
                    self.hairpins[hpinx]=(self.begin+location,structure)
                    self.score[1]+=10
                    #print "too long hairpin at %i: %s, %s"%(location,hpinx,structure)
                    self.bad_seqs.append("long hairpin at %i: %s, %s"%(location,hpinx,structure))
                else:
                    continue
            elif l>l_threshold and details[i][3]<mfe and l_threshold>20:
                hpinx = details[i][0]
                structure = details[i][1]
                gc = self.GC_count(hpinx,len(hpinx))
                if rc:
                    location = self.size -(index + details[i][2] + len(hpinx))
                else:
                    location = index + details[i][2]
                location_tup = (self.begin + location, self.begin + location+len(hpinx) )
                if hpinx not in self.h_info.keys():
                    if any(i<60 for i in location_tup)  or any(i>(self.size-60) for i in location_tup):
                        self.terminal_hairpins+=1
                    self.h_info[hpinx] = {"locations":[(self.begin + location, self.begin + location+len(hpinx) )],"types":["long hairpin"]} #Adding offset (begin) if sequence is a subsequence
                    self.h_loc_list += [((self.begin + location, self.begin + location+len(hpinx) ), 1.0)] #Adding offset (begin) if sequence is a subsequence
                    self.h_count +=1
                    self.long_hairpins+=1
                    self.s_info_2['long_hairpins'].extend([((self.begin + location, self.begin + location+len(hpinx) ),1.0)])
                    self.hairpins[hpinx]=(self.begin+location,structure)
                    self.score[1]+=10
                    #print "too long hairpin at %i: %s, %s"%(location,hpinx,structure)
                    self.bad_seqs.append("long hairpin at %i: %s, %s"%(location,hpinx,structure))
            else:
                continue
        for gquad in gquads:
            hpinx = gquad[0]
            structure = gquad[1]
            #gc = self.GC_count(hpinx,len(hpinx))
            if rc:
                location = self.size -(index + gquad[2] + len(hpinx))
            else:
                location = index + gquad[2]
            location_tup = (self.begin + location, self.begin + location+len(hpinx) )
            if hpinx not in self.h_info.keys():
                if any(i<60 for i in location_tup)  or any(i>(self.size-60) for i in location_tup):
                    self.terminal_hairpins+=1
                self.h_info[hpinx] = {"locations":[(self.begin + location, self.begin + location+len(hpinx) )],"types":["gquadruplex"]} #Adding offset (begin) if sequence is a subsequence
                self.h_loc_list += [((self.begin + location, self.begin + location+len(hpinx) ), 1.0)] #Adding offset (begin) if sequence is a subsequence
                #self.h_count +=1
                self.g_quad_predictions+=1
                self.hairpins[hpinx]=(self.begin+location,structure)
                self.score[1]+=10
                #print "too long hairpin at %i: %s, %s"%(location,hpinx,structure)
                self.bad_seqs.append("g quadruplex predicted at %i: %s, %s"%(location,hpinx,structure))
        return

    def find_hairpins_windowed(self,Hwindow=[50,100,50], 
                               HL = [11,17,21],
                               HGC = 0.8,
                               HMFE = [-10,-15,-20],
                               pal=11):

        self.h_count = 0
        self.h_loc_list = [] 
        self.h_info = {}
        self.hairpins={}

        #metrics for ML
        self.long_hairpins = 0 # hairpins (type1)
        self.strong_hairpins = 0 # hairpins (type2)
        self.wide_hairpins = 0 # hairpins (type3)
        self.longest_hairpin = 0 # length of longest hairpin
        self.richest_hairpin = 0 # GC content of GC richest hairpin
        self.palindromes = 0 # palindromes
        self.terminal_hairpins = 0
        self.g_quad_predictions = 0

        self.s_info_2['long_hairpins']=[]
        self.s_info_2['strong_hairpins']=[]
        self.s_info_2['wide_hairpins']=[]
        self.s_info_2['longest_hairpin']=[]
        self.s_info_2['richest_hairpin']=[]
        self.s_info_2['palindromes']=[]
        self.s_info_2['terminal_hairpins']=[]


        self.DNAmodel = PyVRNA(parameter_file="dna_mathews2004.par",dangles=0,noGU = True, gquad=True) # initialize PyVRNA model object

        for j in range(self.size-100):
            seq1 = self.sequence1[j:j+100]
            seq2 = self.sequence2[self.size-(j+100):self.size-j]
            for i,h in enumerate(Hwindow):
                lengths1, details1, gquads1= self.get_hairpin_lengths(seq1,h)
                lengths2, details2, gquads2 = self.get_hairpin_lengths(seq2,h)
                self.hairpin_processing(j,lengths1,details1,gquads1,HL[i],HGC,HMFE[i])
                self.hairpin_processing(j,lengths2,details2,gquads2,HL[i],HGC,HMFE[i],rc=True)

        self.energies,self.max_bp_hairpin=self.get_windowed_DNA_structure_info(50)
        self.fast_hairpins(fast=False)

    def get_hairpin_lengths(self, seq, maxspan):
        # hairpin length subroutine
        # details = list of tuples (subsequence,substructure,start_pos of str)
        self.DNAmodel.settings.max_bp_span = maxspan
        fold = self.DNAmodel.RNAfold(seq)
        #print fold.structure
        parsed_fold = self.DNAmodel.vienna2bp(fold.structure)
        #print fold.structure
        bp_x = parsed_fold.bpx
        bp_y = parsed_fold.bpy
        bp_gquad = parsed_fold.gquad
        lengths = []
        details = []
        gquads = []
        while len(bp_x) > 0:
            indx = bisect(bp_x,bp_y[0])
            lengths.append(indx) # save number of bp in each hairpin (height)
            details.append((seq[bp_x[0]-1:bp_y[0]],fold.structure[bp_x[0]-1:bp_y[0]],bp_x[0]-1,fold.energy)) # save string and location of each hairpin
            bp_x = bp_x[indx:]
            bp_y = bp_y[indx:]
        for i in range(len(bp_gquad)/3):
            first = bp_gquad[3*i][0] 
            last =  bp_gquad[3*i][3]
            gquads.append((seq[first-1:last],fold.structure[first-1:last],first-1))
        return lengths,details,gquads

    def _is_hairpin_pass(self, seq, stem, loop, max_mismatch, gc_high):
        #quickmode Haripin subroutine
        i = 0
        while i < len(seq)-(stem+loop+stem)+1:
            j, mismatch, gc_count = 0, 0, 0.0

            while j < stem:

                if seq[i+j] in ['G', 'C']:
                    gc_count += 1.0

                if seq[i+(stem+loop+stem)-j-1] in ['G', 'C']:
                    gc_count += 1.0

                if seq[i+j] != seq[i+(stem+loop+stem)-j-1].translate(self.comp_table):
                    mismatch += 1

                if mismatch > max_mismatch:
                    break

                j += 1

            if j == stem:
                gc_content = (gc_count / 2.0) / stem

                if gc_content >= gc_high:
                    hairpin= seq[i:i+(stem+loop+stem)]

                    if hairpin not in self.h_info.keys():
                        if i<60 or i>len(seq)-61:
                            self.h_info[hairpin] = {"locations":[(self.begin + i,i+(stem))],"types":["terminal hairpin"]}
                            self.terminal_hairpins+=1
                            self.s_info_2['terminal_hairpins'].extend([((self.begin + i, self.begin + i+(stem) ),1.0)])
                        if gc_high > 0.50:
                            self.h_info[hairpin] = {"locations":[(self.begin + i,i+(stem))],"types":["strong hairpin"]}
                            self.strong_hairpins+=1
                            self.s_info_2['strong_hairpins'].extend([((self.begin + i, self.begin + i+(stem) ),1.0)])
                            if gc_content>self.richest_hairpin:
                                self.richest_hairpin=gc_content
                                self.s_info_2['richest_hairpin']=[((self.begin + i, self.begin + i+(stem) ),1.0)]
                            if gc_content==self.richest_hairpin:
                                self.s_info_2['richest_hairpin'].extend([((self.begin + i, self.begin + i+(stem) ),1.0)])

                        else:
                            self.h_info[hairpin] = {"locations":[(self.begin + i,i+(stem))],"types":["long hairpin"]}
                            if loop>0 and loop<49:
                                self.long_hairpins+=1
                                self.s_info_2['long_hairpins'].extend([((self.begin + i, self.begin + i+(stem) ),1.0)])
                            elif loop==0:
                                self.palindromes+=1
                                self.s_info_2['palindromes'].extend([((self.begin + i, self.begin + i+(stem) ),1.0)])
                            else:
                               self.wide_hairpins+=1
                               self.s_info_2['wide_hairpins'].extend([((self.begin + i, self.begin + i+(stem) ),1.0)])
                            if gc_content>self.richest_hairpin:
                                self.richest_hairpin=gc_content
                                self.s_info_2['richest_hairpin']=[((self.begin + i, self.begin + i+(stem) ),1.0)]
                            if gc_content==self.richest_hairpin:
                                self.s_info_2['richest_hairpin'].extend([((self.begin + i, self.begin + i+(stem) ),1.0)])

                    else:
                        self.h_info[hairpin]["locations"].append((self.begin + i,i+(stem)))
                        if i<60 or i>len(seq)-61:
                            self.h_info[hairpin]["types"].append("terminal hairpin")
                            self.terminal_hairpins+=1
                            self.s_info_2['terminal_hairpins'].extend([((self.begin + i, self.begin + i+(stem) ),1.0)])
                        if gc_high > 0.50:
                            self.h_info[hairpin]["types"].append("strong hairpin")
                            self.strong_hairpins+=1
                            if gc_content>self.richest_hairpin:
                                self.richest_hairpin=gc_content
                                self.s_info_2['richest_hairpin']=[((self.begin + i, self.begin + i+(stem) ),1.0)]
                            if gc_content==self.richest_hairpin:
                                self.s_info_2['richest_hairpin'].extend([((self.begin + i, self.begin + i+(stem) ),1.0)])
                        else:
                            self.h_info[hairpin]["types"].append("long hairpin")
                            if loop>0 and loop<49:
                                self.long_hairpins+=1
                                self.s_info_2['long_hairpins'].extend([((self.begin + i, self.begin + i+(stem) ),1.0)])
                            elif loop==0:
                                self.palindromes+=1
                                self.s_info_2['palindromes'].extend([((self.begin + i, self.begin + i+(stem) ),1.0)])
                            else:
                                self.wide_hairpins+=1
                                self.s_info_2['wide_hairpins'].extend([((self.begin + i, self.begin + i+(stem) ),1.0)])
                            if gc_content>self.richest_hairpin:
                                self.richest_hairpin=gc_content
                                self.s_info_2['richest_hairpin']=[((self.begin + i, self.begin + i+(stem) ),1.0)]
                            if gc_content==self.richest_hairpin:
                                self.s_info_2['richest_hairpin'].extend([((self.begin + i, self.begin + i+(stem) ),1.0)])

                    self.h_loc_list += [((self.begin + i,i+(stem)), 1.0)]
                    self.h_count +=1
                    self.score[1]+=10
                    self.bad_seqs.append("too long hairpin at position %i-%i:%s "%(i,i+(stem+loop+stem),seq[i:i+(stem+loop+stem)]))

            i += 1

    def fast_hairpins(self, fast = True):
        if fast:
            # hairpin quickmode
            self.h_count = 0
            self.h_loc_list = [] 
            self.h_info = {}

            #metrics for ML
            self.long_hairpins = 0 # hairpins (type1)
            self.strong_hairpins = 0 # hairpins (type2)
            self.wide_hairpins = 0 # hairpins (type3)
            self.longest_hairpin = 0 # length of longest hairpin
            self.richest_hairpin = 0 # GC content of GC richest hairpin
            self.palindromes = 0 # palindromes
            self.terminal_hairpins=0

            self.s_info_2['long_hairpins']=[]
            self.s_info_2['strong_hairpins']=[]
            self.s_info_2['wide_hairpins']=[]
            self.s_info_2['longest_hairpin']=[]
            self.s_info_2['richest_hairpin']=[]
            self.s_info_2['palindromes']=[]
            self.s_info_2['terminal_hairpins']=[]

        # Type 1: Stem = 11bp, Loop = 3bp to 48bp, Mismatches = 2, Stem must have GC-content >= 0.80
        seq=self.sequence1
        stem, max_mismatch, gc_high = 11, 2, 0.80

        for loop in xrange(3, 48+1):

            if len(seq) >= (stem+loop+stem):
                self._is_hairpin_pass(seq, stem, loop, max_mismatch, gc_high)

            else:
                break

        # Type 2: Stem = 17bp, Loop = 3bp to 100bp, Mismatches = 3
        stem, max_mismatch, gc_high = 17, 3, 0.0

        for loop in xrange(3, 100+1):

            if len(seq) >= (stem+loop+stem):
                self._is_hairpin_pass(seq, stem, loop, max_mismatch, gc_high)

            else:
                break

        if len(seq) > 500:
            # Type 3: Stem = 22bp, Loop = 100bp to 500bp, Mismatches = 5
            stem, max_mismatch, gc_high = 22, 5, 0.0

            for loop in xrange(100, 500+1):

                if len(seq) >= (stem+loop+stem):
                    self._is_hairpin_pass(seq, stem, loop, max_mismatch, gc_high)

                else:
                    break

        #Palindromes
        stem, loop, max_mismatch, gc_high = 11, 0, 1, 0.0

        if len(seq) >= (stem+loop+stem):
            self._is_hairpin_pass(seq, stem, loop, max_mismatch, gc_high)

    def get_windowed_DNA_structure_info(self, window_size):
        # hairpin visualization subroutine
        energies = []
        max_bp_hairpin = []   

        for i in range(self.size-window_size):
            seq=self.sequence1[i:i+window_size]
            fold = self.DNAmodel.RNAfold(seq)
            mfe=fold.energy
            #print fold.structure
            parsed_fold = self.DNAmodel.vienna2bp(fold.structure)
            bp_x = parsed_fold.bpx
            bp_y = parsed_fold.bpy
            lengths = []
            details = []

            while len(bp_x) > 0:
                indx = bisect(bp_x,bp_y[0])
                lengths.append(indx) # save number of bp in each hairpin (height)
                details.append((seq[bp_x[0]-1:bp_y[0]],fold.structure[bp_x[0]-1:bp_y[0]],bp_x[0]-1)) # save string and location of each hairpin
                bp_x = bp_x[indx:]
                bp_y = bp_y[indx:]

            lengths=sorted(lengths)
            #print lengths

            if lengths:
                max_bp_hairpin.append(lengths[-1])

            else:
                max_bp_hairpin.append(0)

            energies.append(mfe)

        return energies,max_bp_hairpin

    def find_nucleotides(self,   GCwindow=[20,100],
                                 GC_max=0.67,
                                 GC_low=[0.20,0.30],
                                 GC_high=[0.80,0.70],
                                 dGC=0.50,
                                 T_window=50,
                                 T_bounds=[0.30,0.70],
                                 cons={'A':13,'C':9,'G':9,'T':13},
                                 run_size=[3,2],
                                 run_count=[6,10],
                                 tm_bounds=[40,70],
                                 g_quad_motif="G{3,5}[ATGC]{2,7}",
                                 i_motif="C{2,5}[ATGC]{2,7}"):
        ### Finds sequecnes that violate the following rules:
        # Greater than or equal to X consecutive bases within the sequence   "A=13 C=9 G=9 T=13" 10
        # Overall GC% is equal to or greater than X  x=67    10
        # GC% in any window of X  is less than or equal to  Y %  "X=20, Y=20; X=100, Y=30"   "3; 10"
        # GC% in any window of X bp is equal to or greater than Y %  "X=20, Y=80; X=100, Y=70"   "3;10"
        # GC content of terminal 50 bp is X % GC or greater  X=70    3
        # GC content of terminal 50 bp is X % GC or less X=30    3
        # No more than X consecutive  copies of a  Y nucleotide repeat are present   "X=6, Y=3; X=10, Y=2"   10
        # change in GC content no greater than 50% in 100 bp windows

        self.n_count = 0
        self.n_count_2 = 0
        self.n_loc_list = [] 
        self.n_info = {}

        # Metrics for ML
        self.poly_runs = 0
        self.pattern_runs = 0 
        self.total_GC = 0
        self.GC_overage = 0
        self.GC_short_h = 0
        self.GC_short_l = 0
        self.GC_long_h = 0
        self.GC_long_l = 0
        self.GC_term_h = 0
        self.GC_term_l = 0
        self.dGC = 0
        self.Tm_low = 0
        self.Tm_high = 0
        self.dTm = 0
        self.s_info_2["poly_runs"]=[]
        self.s_info_2["pattern_runs"]=[]
        self.s_info_2["GC_short_l"]=[]
        self.s_info_2["GC_short_h"]=[]
        self.s_info_2["GC_long_l"]=[]
        self.s_info_2["GC_long_h"]=[]
        self.s_info_2["GC_term_h"]=[]
        self.s_info_2["GC_term_l"]=[]
        self.s_info_2["dGC"]=[]
        self.s_info_2["Tm_low"]=[]
        self.s_info_2["Tm_high"]=[]
        self.s_info_2["dTm"]=[]
        self.s_info_2["g_quad_motifs"]=[]
        self.s_info_2["i_motifs"]=[]
        self.s_info_2["total_GC"]=[]


        self.dGC_list = []
        self.dGC_highlights = [] 

        self.dTm_list = []
        self.dTm_highlights = [] 

        self.g_quadruplexes = 0 
        self.i_motifs = 0

        GCcount=self.GC_count(self.sequence1,self.size)
        #print GCcount
        self.total_GC = GCcount
        #GC_count=sum([1 for i in self.sequence1 if (i=='G' or i =='C')])/self.size

        if GCcount>=GC_max:
            #print "GC content of entire sequence too high"
            self.bad_seqs.append("GC content of entire sequence too high")
            self.score[2]+=10
            self.GC_overage=1

        for key in cons:
            rule="%s{%i,1000}"%(key,cons[key])
            self.poly_runs = self.reg_ex(rule, self.sequence1, self.poly_runs, "poly N run")

        for j in range(len(GCwindow)):

            if j==0:
                self.tm_list=[]
                self.GC_list=[]

            for i in range(self.size-GCwindow[j]):
                seq=self.sequence1[i:i+GCwindow[j]]
                val=self.GC_count(seq,GCwindow[j])

                if j==0:
                    tm=tm_calc(seq)
                    self.tm_list.append(tm)
                    self.GC_list.append(val)

                    if val<=GC_low[j]:
                        self.GC_report(i,GCwindow[j],seq,True,False)
                        self.GC_short_l+=1
                        self.score[2]+=3

                    if val>=GC_high[j]:

                        self.GC_short_h+=1
                        self.GC_report(i,GCwindow[j],seq,False,False)
                        self.score[2]+=3

                    if tm<=tm_bounds[0]:
                        self.Tm_low+=1
                        self.GC_report(i,GCwindow[j],seq,True,True)
                        self.score[3]+=3

                    if tm>=tm_bounds[1]:
                        self.Tm_high+=1
                        self.GC_report(i,GCwindow[j],seq,False,True)
                        self.score[3]+=3
                if val<=GC_low[j]:
                    self.GC_long_l+=1
                    self.score[2]+=10

                if val>=GC_high[j]:
                    self.GC_long_h+=1
                    self.score[2]+=10
        
        self.delta_GC_fxn(Tm_mode=False)
        self.delta_GC_fxn(Tm_mode=True)
        if len(self.sequence1)>T_window:
            for i in range(T_window-20):
                fiveterm=self.sequence1[i:i+20]
                fiveGC=self.GC_count(fiveterm,20)

                threeterm=self.sequence1[self.size-T_window+i:self.size-T_window+i+20]
                threeGC=self.GC_count(threeterm,20)
                    
                if fiveGC<T_bounds[0]:
                    seq=self.sequence1[i:i+20]
                    self.GC_term_l+=1
                    self.GC_report(i,20,seq,True,False,True)
                    self.score[2]+=10

                if fiveGC>T_bounds[1]:
                    seq=self.sequence1[i:i+20]
                    self.GC_term_h+=1
                    self.GC_report(i,20,seq,False,False,True)
                    self.score[2]+=10

                if threeGC<T_bounds[0]:
                    seq=self.sequence1[self.size-T_window+i:self.size-T_window+i+20]
                    self.GC_term_l+=1
                    self.GC_report(self.size-T_window+i,20,seq,True,False,True)
                    self.score[2]+=10

                if threeGC>T_bounds[1]:
                    seq=self.sequence1[self.size-T_window+i:self.size-T_window+i+20]
                    self.GC_term_h+=1
                    self.GC_report(self.size-T_window+i,20,seq,False,False,True)
                    self.score[2]+=10
        else:
            for i in range(len(self.sequence1)-20):
                fiveterm=self.sequence1[i:i+20]
                fiveGC=self.GC_count(fiveterm,20)

                threeterm=self.sequence1[i:i+20]
                threeGC=self.GC_count(threeterm,20)
                    
                if fiveGC<T_bounds[0]:
                    seq=self.sequence1[i:i+20]
                    self.GC_term_l+=1
                    self.GC_report(i,20,seq,True,False,True)
                    self.score[2]+=10

                if fiveGC>T_bounds[1]:
                    seq=self.sequence1[i:i+20]
                    self.GC_term_h+=1
                    self.GC_report(i,20,seq,False,False,True)
                    self.score[2]+=10

                if threeGC<T_bounds[0]:
                    seq=self.sequence1[i:i+20]
                    self.GC_term_l+=1
                    self.GC_report(i,20,seq,True,False,True)
                    self.score[2]+=10

                if threeGC>T_bounds[1]:
                    seq=self.sequence1[i:i+20]
                    self.GC_term_h+=1
                    self.GC_report(i,20,seq,False,False,True)
                    self.score[2]+=10    

        for i in range(len(run_size)):
            merlist=mers(run_size[i])

            for j in merlist:
                rule=j*run_count[i]
                self.pattern_runs = self.reg_ex(rule,self.sequence1,self.pattern_runs,"mer_run")
        self.g_quadruplexes = self.reg_ex(g_quad_motif*3,self.sequence1,self.g_quadruplexes, "g_quadruplex")
        self.i_motifs=self.reg_ex(i_motif*3,self.sequence1,self.i_motifs, "i_motifs")
        if self.total_GC>0.5:
            self.s_info_2['total_GC'].extend(self.s_info_2['GC_short_h'])
            self.s_info_2['total_GC'].extend(self.s_info_2['GC_long_h'])
            self.s_info_2['total_GC'].extend(self.s_info_2['GC_term_h'])
            #print self.s_info_2
        else:
            self.s_info_2['total_GC'].extend(self.s_info_2['GC_short_l'])
            self.s_info_2['total_GC'].extend(self.s_info_2['GC_long_l'])
            self.s_info_2['total_GC'].extend(self.s_info_2['GC_term_l'])
            #print self.s_info_2

        

    def GC_report(self,i, window, seq, low = None, Tm = None, terminal= None):
        if low:
            nucs = "AT"
            report1 = "low"  
        else:
            nucs="GC"
            report1 = "high"
        if Tm:
            report2 = "Tm"
        else: 
            report2 = "GC"
        if terminal:
            term = "terminal "
        else:
            term = ""
        
        sub_seqs=[]
        try: 
            for k in range(window):
                if seq[k] in nucs:
                    sub_seqs.append(i+k)
        except:
            print len(seq)
            print window
        #sub_seqs = [i+k for k in range(window) if seq[k] in nucs]
        new_seqs = []
        if len(sub_seqs)>0:
            five_p = sub_seqs[0]
            for n,s in enumerate(sub_seqs[1:]):
                if s == sub_seqs[n]+1 and s != sub_seqs[-1]:
                    continue
                elif s == sub_seqs[n]+1 and s == sub_seqs[-1]: 
                    new_seqs.append((self.begin + five_p,self.begin + sub_seqs[-1]))
                else: 
                    new_seqs.append((self.begin + five_p,self.begin + sub_seqs[n]))
                    five_p = s
                    if s == sub_seqs[-1]:
                        new_seqs.append((self.begin + s,self.begin + s))
            
            for k in new_seqs:
                if k[1]>k[0]+2:
                    if seq not in self.n_info.keys():
                        self.n_info[seq] = {"locations":[ k ],"types":["{}{} {}".format(term,report1,report2)]}

                    else:
                        self.n_info[seq]["locations"].append( k )
                        self.n_info[seq]["types"].append("{}{} {}".format(term,report1,report2))
                    if window==100:
                        if low:
                            self.s_info_2["GC_long_l"].extend([(k,1.0)])
                        else:
                            self.s_info_2["GC_long_h"].extend([(k,1.0)])
                    else:
                        if terminal:
                            if low:
                                self.s_info_2["GC_term_l"].extend([(k,1.0)])
                            else:
                                self.s_info_2["GC_term_h"].extend([(k,1.0)])
                        elif Tm:
                            if low:
                                self.s_info_2["Tm_low"].extend([(k,1.0)])
                            else:
                                self.s_info_2["Tm_high"].extend([(k,1.0)])
                        else:
                            if low:
                                self.s_info_2["GC_short_l"].extend([(k,1.0)])
                            else:
                                self.s_info_2["GC_short_h"].extend([(k,1.0)])
                    self.n_loc_list += [(k, 1.0)]
                    self.n_count +=1
            self.n_count_2 +=1
            self.bad_seqs.append("{} {} in a {} bp {}window at position {}: {}".format(report1,report2,window, term,i,seq))

    def delta_GC_fxn(self, Tm_mode=None, thresholds=[0.5, 30]):
        if Tm_mode:
            strtype="Tm"
            feed_list = self.tm_list
            d_list = self.dTm_list
            highlights= self.dTm_highlights
            th = thresholds[1]
        else:
            strtype = "GC"
            feed_list = self.GC_list
            d_list = self.dGC_list
            highlights = self.dGC_highlights
            th = thresholds[0]
        last_seq=[None,None]
        for i in range(len(feed_list)-80):
            local = feed_list[i:i+80]
            hi=max(local)
            lo=min(local)
            d_list.append(hi-lo)
            if hi-lo>th:
                hi_seq=self.sequence[i+local.index(hi):i+local.index(hi)+20]
                lo_seq=self.sequence[i+local.index(lo):i+local.index(lo)+20]
                if hi_seq != last_seq[1] or lo_seq != last_seq[0]:
                    if i+local.index(hi) not in highlights:
                        highlights.append(i+local.index(hi))
                    if i+local.index(lo) not in highlights:
                        highlights.append(i+local.index(lo))
                    last_seq=[lo_seq, hi_seq]
                    seq=self.sequence1[i:i+100]
                    if Tm_mode:
                        self.dTm+=1
                    else:
                        self.dGC+=1
                    sub_seqs_l = [i+local.index(lo)+k for k in range(20) if last_seq[0][k] in "AT"]
                    sub_seqs_h = [i+local.index(hi)+k for k in range(20) if last_seq[1][k] in "GC"]
                    new_seqs = []
                    five_p = sub_seqs_l[0]
                    for n,s in enumerate(sub_seqs_l[1:]):
                        #print s
                        #print sub_seqs[n]+1
                        if s == sub_seqs_l[n]+1 and s != sub_seqs_l[-1]:
                            continue
                        elif s == sub_seqs_l[n]+1 and s == sub_seqs_l[-1]: 
                            new_seqs.append((self.begin + five_p,self.begin + sub_seqs_l[-1]))
                        else: 
                            new_seqs.append((self.begin + five_p,self.begin + sub_seqs_l[n]))
                            five_p = s
                            if s == sub_seqs_l[-1]:
                                new_seqs.append((self.begin + s,self.begin + s))
                    five_p = sub_seqs_h[0]
                    for n,s in enumerate(sub_seqs_h[1:]):
                        #print s
                        #print sub_seqs[n]+1
                        if s == sub_seqs_h[n]+1 and s != sub_seqs_h[-1]:
                            continue
                        elif s == sub_seqs_h[n]+1 and s == sub_seqs_h[-1]: 
                            new_seqs.append((self.begin + five_p,self.begin + sub_seqs_h[-1]))
                        else: 
                            new_seqs.append((self.begin + five_p,self.begin + sub_seqs_h[n]))
                            five_p = s
                            if s == sub_seqs_h[-1]:
                                new_seqs.append((self.begin + s,self.begin + s))
                    #print sub_seqs
                    #print new_seqs
                    for k in new_seqs:
                        if k[1]>k[0]+2:
                            if Tm_mode:
                                self.s_info_2["dTm"].extend([(k,1.0)])
                            else:
                                self.s_info_2["dGC"].extend([(k,1.0)])
                            if seq not in self.n_info.keys():
                                self.n_info[seq] = {"locations":[ k ],"types":["delta {}".format(strtype)]}

                            else:
                                self.n_info[seq]["locations"].append( k )
                                self.n_info[seq]["types"].append("delta {}".format(strtype))

                            self.n_loc_list += [(k, 1.0)]
                            self.n_count +=1
                    self.n_count_2 +=1
                    self.bad_seqs.append("delta {} in a 100 bp window at position {}: {}, {}".format(strtype,i,lo_seq, hi_seq))
        if Tm_mode:
            self.dTm_list = d_list
            self.dTm_highlights = highlights
        else:
            self.dGC_list = d_list 
            self.dGC_highlights = highlights

    def reg_ex(self, rule, seq, metrictype, ruletype):

        run=re.finditer(r"%s"%(rule),seq)
        m = [[n.start(), n.end()] for n in run]

        if m>0:

            for r in m:
                run = self.sequence1[r[0]:r[1]]
                metrictype+=1
                if run not in self.n_info.keys():
                    self.n_info[run] = {"locations":[ (self.begin + r[0],self.begin + r[1]) ],"types":[ruletype]}

                else:
                    self.n_info[run]["locations"].append( (self.begin + r[0],self.begin + r[1]) )
                    self.n_info[run]["types"].append(ruletype)

                if ruletype=="poly N run":
                    self.s_info_2["poly_runs"].extend([((self.begin + r[0],self.begin + r[1]),1.0)])
                if ruletype=="mer_run":
                    self.s_info_2["pattern_runs"].extend([((self.begin + r[0],self.begin + r[1]),1.0)])
                if ruletype=="g_quadruplex":
                    self.s_info_2["g_quad_motifs"].extend([((self.begin + r[0],self.begin + r[1]),1.0)])
                if ruletype=="i_motifs":
                    self.s_info_2["i_motifs"].extend([((self.begin + r[0],self.begin + r[1]),1.0)])
                self.n_loc_list += [(self.begin + r[0],self.begin + r[1])]
                self.n_count +=1
                self.n_count_2 +=1
                # print "%s run at position %i"%(motif,r[0])
                self.bad_seqs.append("%s at position %i"%(ruletype,r[0]))
                self.score[2]+=10
        return metrictype

    def GC_count(self,seq,window):
        counts=dict(collections.Counter(seq))

        if 'G' in counts.keys() and 'C' in counts.keys():
            val=float(counts['G']+counts['C'])/window

        elif 'G' in counts.keys():
            val=float(counts['G'])/window

        elif 'C' in counts.keys():
            val=float(counts['C'])/window

        else:
            val=0

        return val

def tm_calc(primer_seq):
    ## calculates melting temperature
    tm=TM.Tm_NN(primer_seq,dnac1=250.0, saltcorr=7)

    return tm

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


def seq_evaluate(seq, verbose=False, begin=0, jsonify=False):
    assess=seq_assess(seq, begin=0, verbose = verbose)
    assess.run()

    results = OrderedDict([     ('scores', assess.score), # lumped scores based on specifications (deprecated)
                   ('r_count', assess.r_count), # repeat count
                   ('h_count', assess.h_count), # hairpin count
                   ('n_count', assess.n_count_2), # nucleotide complexity count
                  ('total_GC', assess.total_GC),
                      ('size', assess.size),
                   ('tm_list', assess.tm_list), # 20bp windowed measurement of melting temperature
                   ('GC_list', assess.GC_list), # 20bp windowed measurement of GC content
                  ('dGC_list', assess.dGC_list), #100bp windowed evaluate of dGC
                  ('dTm_list', assess.dTm_list), # 100bp windowed evaluate of dGC
                ('h_loc_list', assess.h_loc_list), # list of hairpin locations as tuples
          ('folding_energies', assess.energies), # folding energy of structures found in 50 bp windows across the sequence 
            ('max_bp_hairpin', assess.max_bp_hairpin), # longest hairpin stem found in 50 bp windows across the sequence   
            ('repeat_density', assess.repeat_density), # density of disallowed repeats(9+bp) found in 70 bp windows 
       ('repeat_distribution', assess.repeat_distribution), # binned repeats
              ('flag_strings', assess.bad_seqs), # warnings for each flagged motif
                 ('r_metrics', assess.r_metrics),
                 ('h_metrics', assess.h_metrics),
                 ('n_metrics', assess.n_metrics)])
    
    if jsonify:
        json.dumps(results)
    return results

def biggest_contributor_2(rf, X_test, number):
    predictions,bias,contributions = ti.predict(rf,X_test)

    cancel_p= [[contributions[i][j][0] for j in range(len(contributions[i]))]for i in range(len(contributions))]
    tally=[]
    for i in range(len(cancel_p)):
        tally.append([])
        c= sorted(zip(cancel_p[i], X_test.columns), key=lambda x: -x[0])
        if predictions[i][1]<.5:
            for j in range(number):
                tally[i].append((c[j][1],c[j][0],X_test[c[j][1]].iloc[i]))
            #print    
            
        if predictions[i][1]>.5:
            #tally.append([])
            for j in range(number):
                tally[i].append((c[j][1],c[j][0],X_test[c[j][1]].iloc[i]))
            #print
    #tally = pd.DataFrame(tally,columns=X_test.columns,index=[0])
    return tally

def contributors(rf, X_test):
    predictions,bias,contributions = ti.predict(rf,X_test)

    cancel_p= [[contributions[i][j][0] for j in range(len(contributions[i]))]for i in range(len(contributions))]
    feature_list=[i for i in X_test.columns]

    tally=[]
    for i in range(len(cancel_p)):
        tally.append({})
        for j in range(len(feature_list)):
            tally[i][feature_list[j]]=cancel_p[i][j]
    
    #tally = pd.DataFrame(tally, columns=feature_list)
    return tally


def SSC_evaluate(seq_list, verbose = False, jsonify=False):
    print "Sequence List"
    print seq_list
    
    all_results = []
    repeat_dict = OrderedDict()
    hairpin_dict = OrderedDict()
    nucleotide_dict = OrderedDict()
    s_info_dict=OrderedDict()
    
    for (i, seq) in enumerate(seq_list):
        print "Running sequence ", seq
        assess=seq_assess(seq,begin=0,verbose = verbose)
        assess.run()
        
        
        repeat_dict[seq] = assess.r_metrics
        hairpin_dict[seq] = assess.h_metrics
        nucleotide_dict[seq] = assess.n_metrics
        s_info_dict[seq] = assess.s_info_2
        
        
        result = {'total_GC'        : assess.total_GC,
                  'size'            : assess.size,
                  'tm_list'         : assess.tm_list, # 20bp windowed measurement of melting temperature
                  'GC_list'         : assess.GC_list, # 20bp windowed measurement of GC content
                  'dGC_list'        : assess.dGC_list, #100bp windowed evaluate of dGC
                  'dTm_list'        : assess.dTm_list, # 100bp windowed evaluate of dGC
                  'h_loc_list'      : assess.h_loc_list, # list of hairpin locations as tuples
                  'folding_energies': assess.energies, # folding energy of structures found in 50 bp windows across the sequence 
                  'max_bp_hairpin'  : assess.max_bp_hairpin, # longest hairpin stem found in 50 bp windows across the sequence   
                  'repeat_density'  : assess.repeat_density, # density of disallowed repeats(9+bp) found in 70 bp windows 
                  'repeat_distribution': assess.repeat_distribution,
                  'repeat_dictionary' : assess.RepeatDict
                }
        all_results.append(result)
        
    repeat_frame = pd.DataFrame(repeat_dict).T
    hairpin_frame = pd.DataFrame(hairpin_dict).T
    nucleotide_frame = pd.DataFrame(nucleotide_dict).T
    
    dirname, filename = os.path.split(os.path.abspath(__file__))
    forestModel = pkl.load(open(os.path.join(dirname,'SSC_forest.pkl'),'r'))
    SSC=forestModel["forest"]
    features=forestModel["features"]
    
    all_frame=pd.concat([repeat_frame, hairpin_frame, nucleotide_frame], axis=1)
    all_frame_test=all_frame[features]
    y_pred=SSC.predict(all_frame_test)
    y_proba=SSC.predict_proba(all_frame_test)[:,1]
    tally=contributors(SSC, all_frame_test)
    
    IMPORTANCE_THRESHOLD = 0.01
    
    for (i, seq) in enumerate(seq_list):
        current_features = tally[i].keys()
        importance_scores = {feature : importance for (feature, importance) in tally[i].items() if importance > IMPORTANCE_THRESHOLD}
        total_importance = sum(importance_scores.values())
        normalized_importance_scores = {feature : importance / total_importance for (feature, importance) in importance_scores.items()}
        
        target_locations = {f : s_info_dict[seq][f] for f in current_features}
        target_locations_with_importance_scores = []
        for (feature, items) in target_locations.items():
            if feature in normalized_importance_scores:
                importance = normalized_importance_scores[feature]
                target_locations_with_importance_scores.extend( [ ((begin, end), importance) for ((begin, end), one) in items] )
        target_locations_with_importance_scores.sort(key = lambda location: location[0][0], reverse=False)
        
        info_dict={}
        for k in target_locations.keys():
            for l in target_locations[k]:
                subseq=seq[l[0][0]:l[0][1]]
                if subseq not in info_dict.keys():
                    info_dict[subseq]={"types":[k],
                                        "locations": [l]}
                else:
                    info_dict[subseq]["locations"].extend([l])
                    info_dict[subseq]["types"].extend([k])
        
        all_results[i]["synthesis_failure_contributions"]=tally[i]
        all_results[i]["synthesis_probability"]=y_proba[i]
        if y_pred[i] == 'Synthesized':
            outcome = True
        else:
            outcome = False
        all_results[i]["synthesis_outcome"] = outcome
        all_results[i]["synthesis_contributor_loc_dict"]=target_locations
        all_results[i]["synthesis_contributor_locations_with_priorities"] = target_locations_with_importance_scores
        all_results[i]["synthesis_contributor_info_dict"]=info_dict
        
        important_features = []
        highlighted_regions = {}

        for (feature, score) in all_results[i]["synthesis_failure_contributions"].items():
            if score >= IMPORTANCE_THRESHOLD:
                important_features.append(feature)
                locationSnippetList = all_results[i]["synthesis_contributor_loc_dict"][feature]
                highlighted_regions[feature] = []
                for location in locationSnippetList:
                    begin = location[0][0]
                    end = location[0][1]
                    interior_point_list = [location for location in highlighted_regions[feature] if begin >= location[0] and end <= location[1]]
                    if len(interior_point_list) == 0:   highlighted_regions[feature].append( [begin, end] )
            
                highlighted_regions[feature].sort(key = lambda location: location[0], reverse=False)
        
        all_highlighted_regions = []
        for (feature, locationList) in highlighted_regions.items():
            for location in locationList:
                begin = location[0]
                end = location[1]
                interior_point_list = [location for location in all_highlighted_regions if begin >= location[0] and end <= location[1]]
                if len(interior_point_list) == 0:
                    all_highlighted_regions.append( [begin, end, [feature] ] )
                else:
                    indexList = [index for index, location in enumerate(all_highlighted_regions) if location[0] == interior_point_list[0][0]]
                    for index in indexList:
                        all_highlighted_regions[index][2].append(feature)

        all_highlighted_regions.sort(key = lambda location: location[0], reverse=False)
        all_results[i]['highlighted_regions_with_features'] = all_highlighted_regions

    if jsonify:
        return json.dumps(all_results)
    return all_results

if __name__=="__main__":
    
    sequences = pd.read_csv('SupplementaryData1.csv', header=0).loc[1:1, 'Sequence']
    results = SSC_evaluate(sequences)
    print results
    