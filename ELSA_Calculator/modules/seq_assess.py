import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), "../.."))  
#print sys.path
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
#from RepeatFinder import RepeatFinder
import itertools
import time 
from string import maketrans
import json
from Mer_Maker import mers


class seq_assess(object):
    def __init__(self,sequence=None, begin=0, verbose=None):
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
            sequence1 = str(sequence.seq).upper()    
            sequence2 = str(sequence.seq.reverse_complement()).upper()
        except:
            try:
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
        self.sequence=sequence1
        self.sequence1=sequence1
        self.sequence2=sequence2

    def run(self):
        ### quickmode for running as an evaluate method
        
        self.find_repeats(RL=9,Rmax=11)
        self.find_hairpins_windowed()
        self.find_nucleotides()
        if self.verbose:
            for i in self.bad_seqs:
                print i


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

    def run_ELSA_analysis(self):
        ## quickmode for running as an evaluate method
        
        self.find_repeats(RL=9,Rmax=11)
        self.fast_hairpins()
        self.find_nucleotides()

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
                           Rmax=20, 
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

        ## FastFinder, courtesy of Ayaan Hossain
        RP=FastFinder()

        self.RepeatDict=RP.get_repeat_dict([self.sequence1],RL,verbose=False)

        total=np.zeros(self.size, dtype=bool)
        total2=np.zeros(self.size)
        r_distribution={} # for gathering distributions of repeats for visualizations

        location={}

        for i in self.RepeatDict.keys():
            self.r_info[i]={"locations":[],"types":[]}
            val=np.zeros(self.size, dtype=bool)

            if len(i) > Rmax: # F)
                self.bad_seqs.append("This repeat is too long: %s"%(i))
                        
            for j in self.RepeatDict[i].keys():
                location_list = list(self.RepeatDict[i][j][0])
                self.r_info[i]["types"].extend([j])
                self.r_info[i]["locations"].extend([(self.begin + k, self.begin + k+len(i)) for k in location_list])

                for k in location_list:
                    if len(i) > Rmax:
                        self.score[0]+=len(i)
                    
                     
                    self.r_count+=1
                    self.r_loc_list += [( (self.begin + k, self.begin + k+len(i)), 1.0)]

                    for v in range(k,k+len(i)):
                        val[v]=True

                    if len(i) not in r_distribution:
                        #for l in range(RL,len(i)):
                        #    if l not in r_distribution:
                        #        r_distribution[l] = 0
                        r_distribution[len(i)] = 1
                    else:
                        r_distribution[len(i)]+=1
    
            ind_fract=float(sum(val.astype(int)))/self.size
            if ind_fract> Rthreshold[1]: # B)
                self.bad_seqs.append("Too much of this repeat: %s. Redesign to remove "%(i))
                self.score[0]+=6


            total+=val
            total2+=val.astype(int)

        total_fract=float(sum(total.astype(int)))/self.size
        if total_fract>Rthreshold[0]: # A)
            self.bad_seqs.append("Too many repeats. Redesign to remove ")
            self.score[0]+=10

        windowed_repeats=[] # for more visuals
        for i in range(len(Rwindow)-1): # C) & D)

            for j in range(self.size-Rwindow[i]):
                peek = float(sum(total[j:j+Rwindow[i]].astype(int)))

                if peek/Rwindow[i] > Rthreshold[2+i]: 
                    #self.bad_seqs.append("%i bp window starting at %i is too repeat rich. Redesign to remove "%(Rwindow[i],j))
                    self.score[0]+=10
                if i==0:
                    windowed_repeats.append(peek/Rwindow[i])
        # E)
        term1=self.sequence1[:Rwindow[2]]
        self.tandem_subroutine(term1,Rtandem,"five")

        term2=self.sequence1[len(self.sequence1)-Rwindow[2]:]
        self.tandem_subroutine(term2,Rtandem,"three")

        self.repeat_density=windowed_repeats
        rep_list=[]
        for k,v in r_distribution.items():
            rep_list.append((k,v))
        self.repeat_distribution=rep_list

    def tandem_subroutine(self,term,Rtandem,end):
        # subroutine that finds tandem repeats using regex
        i=int(Rtandem)
        L=int(Rtandem)

        while i <len(term):
            seq = term[(i - L):i]
            run = re.finditer(r"%s%s"%(seq,seq), term)
            m = [[n.start(), n.end()] for n in run]

            if len(m) > 0:
                tandrep=''.join([seq,seq])
                self.r_count+=1
                self.r_info[tandrep]={"locations":[],"types":[]}

                if end == "five":
                    self.r_info[tandrep]["locations"].append( (self.begin + i, self.begin + i + len(tandrep)) ) #Adding offset (begin) if sequence is a subsequence
                    self.r_info[tandrep]["types"].append("tandem")
                    self.r_loc_list += [( (self.begin + i, self.begin + i + len(tandrep)), 1.0)] #Adding offset (begin) if sequence is a subsequence
                    self.bad_seqs.append("tandem repeat too close to 5' end: %s"%(tandrep))

                elif end == "three":
                    self.r_info[tandrep]["locations"].append( (self.begin + self.size - 60 + i-L, self.begin + self.size - 60 + i-L+len(tandrep) ) ) #Adding offset (begin) if sequence is a subsequence
                    self.r_info[tandrep]["types"].append("tandem")
                    self.r_loc_list += [( (self.begin + self.size - 60 + i-L, self.begin + self.size - 60 + i-L+len(tandrep)), 1.0)] #Adding offset (begin) if sequence is a subsequence
                    self.bad_seqs.append("tandem repeat too close to 3' end: %s"%(tandrep))

                else:
                    raise InputError("internal error in tandem repeat finding")

                self.score[0]+=10

            i+=1


    def hairpin_processing(self,index,lengths,details,l_threshold, gc_threshold, mfe, rc=False):
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
                #print location_tup
                if gc>gc_threshold:
                    if hpinx not in self.h_info.keys():
                        self.h_info[hpinx] = {"locations":[(self.begin + location, self.begin + location+len(hpinx) )],"types":["strong hairpin"]} #Adding offset (begin) if sequence is a subsequence
                        self.h_loc_list += [((self.begin + location, self.begin + location+len(hpinx) ), 1.0)] #Adding offset (begin) if sequence is a subsequence
                        self.h_count +=1
                        self.hairpins[hpinx]=(self.begin+location,structure)
                        self.score[1]+=10
                        #print "too long hairpin at %i: %s, %s"%(location,hpinx,structure)
                        self.bad_seqs.append("strong hairpin at %i: %s, %s"%(location,hpinx,structure))
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
                    self.h_info[hpinx] = {"locations":[(self.begin + location, self.begin + location+len(hpinx) )],"types":["long hairpin"]} #Adding offset (begin) if sequence is a subsequence
                    self.h_loc_list += [((self.begin + location, self.begin + location+len(hpinx) ), 1.0)] #Adding offset (begin) if sequence is a subsequence
                    self.h_count +=1
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
                    self.h_info[hpinx] = {"locations":[(self.begin + location, self.begin + location+len(hpinx) )],"types":["long hairpin"]} #Adding offset (begin) if sequence is a subsequence
                    self.h_loc_list += [((self.begin + location, self.begin + location+len(hpinx) ), 1.0)] #Adding offset (begin) if sequence is a subsequence
                    self.h_count +=1
                    self.hairpins[hpinx]=(self.begin+location,structure)
                    self.score[1]+=10
                    #print "too long hairpin at %i: %s, %s"%(location,hpinx,structure)
                    self.bad_seqs.append("long hairpin at %i: %s, %s"%(location,hpinx,structure))
            else:
                continue
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

        self.DNAmodel = PyVRNA(parameter_file="dna_mathews2004.par",dangles=0,noGU = True) # initialize PyVRNA model object

        for j in range(self.size-100):
            seq1 = self.sequence1[j:j+100]
            seq2 = self.sequence2[self.size-(j+100):self.size-j]
            for i,h in enumerate(Hwindow):
                lengths1, details1= self.get_hairpin_lengths(seq1,h)
                lengths2, details2= self.get_hairpin_lengths(seq2,h)
                self.hairpin_processing(j,lengths1,details1,HL[i],HGC,HMFE[i])
                self.hairpin_processing(j,lengths2,details2,HL[i],HGC,HMFE[i],rc=True)

        self.energies,self.max_bp_hairpin=self.get_windowed_DNA_structure_info(50)
        self.fast_hairpins(fast=False)

    def get_hairpin_lengths(self,seq,maxspan):
        # hairpin length subroutine
        # details = list of tuples (subsequence,substructure,start_pos of str)
        self.DNAmodel.settings.max_bp_span = maxspan
        fold = self.DNAmodel.RNAfold(seq)
        parsed_fold = self.DNAmodel.vienna2bp(fold.structure)

        bp_x = parsed_fold.bpx
        bp_y = parsed_fold.bpy
        lengths = []
        details = []

        while len(bp_x) > 0:
            indx = bisect(bp_x,bp_y[0])
            lengths.append(indx) # save number of bp in each hairpin (height)
            details.append((seq[bp_x[0]-1:bp_y[0]],fold.structure[bp_x[0]-1:bp_y[0]],bp_x[0]-1,fold.energy)) # save string and location of each hairpin
            bp_x = bp_x[indx:]
            bp_y = bp_y[indx:]

        return lengths,details

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

                        if gc_high > 0.50:
                            self.h_info[hairpin] = {"locations":[(self.begin + i,self.begin + i+(stem))],"types":["strong hairpin"]}

                        else:
                            self.h_info[hairpin] = {"locations":[(self.begin + i,self.begin + i+(stem))],"types":["long hairpin"]}

                    else:
                        self.h_info[hairpin]["locations"].append((self.begin + i,self.begin + i+(stem)))

                        if gc_high > 0.50:
                            self.h_info[hairpin]["types"].append("strong hairpin")

                        else:
                            self.h_info[hairpin]["types"].append("long hairpin")

                    self.h_loc_list += [((self.begin + i,self.begin + i+(stem)), 1.0)]
                    self.h_count +=1
                    self.score[1]+=10
                    self.bad_seqs.append("too long hairpin at position %i-%i:%s "%(i,i+(stem+loop+stem),seq[i:i+(stem+loop+stem)]))

            i += 1

    def fast_hairpins(self,fast=True):
        if fast:
            # hairpin quickmode
            self.h_count = 0
            self.h_loc_list = [] 
            self.h_info = {}


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

    def get_windowed_DNA_structure_info(self,window_size):
        # hairpin visualization subroutine
        energies = []
        max_bp_hairpin = []   

        for i in range(self.size-window_size):
            seq=self.sequence1[i:i+window_size]
            fold = self.DNAmodel.RNAfold(seq)
            mfe=fold.energy
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

        self.dGC_list = []
        self.dGC_highlights = [] 

        self.dTm_list = []
        self.dTm_highlights = [] 

        GCcount=self.GC_count(self.sequence1,self.size)
        #print GCcount
        #GC_count=sum([1 for i in self.sequence1 if (i=='G' or i =='C')])/self.size

        if GCcount>=GC_max:
            #print "GC content of entire sequence too high"
            self.bad_seqs.append("GC content of entire sequence too high")
            self.score[2]+=10

        for key in cons:
            rule="%s{%i,1000}"%(key,cons[key])
            self.reg_ex(rule, self.sequence1, "poly N run")

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
                        self.score[2]+=3

                    if val>=GC_high[j]:
                        self.GC_report(i,GCwindow[j],seq,False,False)
                        self.score[2]+=3

                    if tm<=tm_bounds[0]:
                        self.GC_report(i,GCwindow[j],seq,True,True)
                        self.score[3]+=3

                    if tm>=tm_bounds[1]:
                        self.GC_report(i,GCwindow[j],seq,False,True)
                        self.score[3]+=3
                if val<=GC_low[j]:
                    self.score[2]+=10

                if val>=GC_high[j]:
                    self.score[2]+=10
        
        self.delta_GC_fxn(Tm_mode=False)
        self.delta_GC_fxn(Tm_mode=True)
        
        for i in range(T_window-20):
            fiveterm=self.sequence1[i:i+20]
            fiveGC=self.GC_count(fiveterm,20)

            threeterm=self.sequence1[self.size-T_window+i:self.size-T_window+i+20]
            threeGC=self.GC_count(threeterm,20)
                
            if fiveGC<T_bounds[0]:
                seq=self.sequence1[i:i+20]
                self.GC_report(i,20,fiveterm,True,False,True)
                self.score[2]+=10

            if fiveGC>T_bounds[1]:
                seq=self.sequence1[i:i+20]
                self.GC_report(i,20,fiveterm,False,False,True)
                self.score[2]+=10

            if threeGC<T_bounds[0]:
                seq=self.sequence1[self.size-T_window+i:self.size-T_window+i+20]
                self.GC_report(self.size-T_window+i,20,threeterm,True,False,True)
                self.score[2]+=10

            if threeGC>T_bounds[1]:
                seq=self.sequence1[self.size-T_window+i:self.size-T_window+i+20]
                self.GC_report(self.size-T_window+i,20,threeterm,False,False,True)
                self.score[2]+=10

        for i in range(len(run_size)):
            merlist=mers(run_size[i])

            for j in merlist:
                rule=j*run_count[i]
                self.reg_ex(rule,self.sequence1,"mer_run")
        self.reg_ex(g_quad_motif*3,self.sequence1, "g_quadruplex")
        self.reg_ex(i_motif*3,self.sequence1, "i_motifs")

        

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
        #print len(seq)
        #print type(window)
        #print 
        sub_seqs = [i+k for k in range(window) if seq[k] in nucs]
        if len(sub_seqs)==0:
            print seq, window, nucs
        new_seqs = []
        five_p = sub_seqs[0]
        for n,s in enumerate(sub_seqs[1:]):
            #print s
            #print sub_seqs[n]+1
            if s == sub_seqs[n]+1 and s != sub_seqs[-1]:
                continue
            elif s == sub_seqs[n]+1 and s == sub_seqs[-1]: 
                new_seqs.append((self.begin + five_p,self.begin + sub_seqs[-1]))
            else: 
                new_seqs.append((self.begin + five_p,self.begin + sub_seqs[n]))
                five_p = s
                if s == sub_seqs[-1]:
                    new_seqs.append((self.begin + s,self.begin + s))
        #print sub_seqs
        #print new_seqs
        for k in new_seqs:
            if seq not in self.n_info.keys():
                self.n_info[seq] = {"locations":[ k ],"types":["{}{} {}".format(term,report1,report2)]}

            else:
                self.n_info[seq]["locations"].append( k )
                self.n_info[seq]["types"].append("{}{} {}".format(term,report1,report2))

            self.n_loc_list += [(k, 1.0)]
            self.n_count +=1
        self.n_count_2 +=1
        self.bad_seqs.append("{} {} in a {} bp {}window at position {}: {}".format(report1,report2,window, term,i,seq))

    def delta_GC_fxn(self,Tm_mode=None, thresholds=[0.5,30]):
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
                    self.score[2]+=10
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

    def reg_ex(self, rule, seq, ruletype):

        run=re.finditer(r"%s"%(rule),seq)
        m = [[n.start(), n.end()] for n in run]

        if m>0:

            for r in m:
                run = self.sequence1[r[0]:r[1]]
                if run not in self.n_info.keys():
                    self.n_info[run] = {"locations":[ (self.begin + r[0],self.begin + r[1]) ],"types":[ruletype]}

                else:
                    self.n_info[run]["locations"].append( (self.begin + r[0],self.begin + r[1]) )
                    self.n_info[run]["types"].append(ruletype)

                self.n_loc_list += [((self.begin + r[0],self.begin + r[1]), 1.0)]
                self.n_count +=1
                self.n_count_2 +=1
                # print "%s run at position %i"%(motif,r[0])
                self.bad_seqs.append("%s at position %i"%(ruletype,r[0]))
                self.score[2]+=10

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



if __name__=="__main__":
    
    
    sequence = 'CAACCAGAAAAGCCAACCTGCGGGTTGGCTTTTTTATGCAAATAGGATCCTAGTTTATTCGCTCTATTGAGGTAGTCGTCAGAACCCTTATCTTGACATTTCGTCAAGAGTCGACTATAATATCGCGGCGATAGTTGATCCTCAGCGGTTTTAGATCACGAAAGTGAAAGTTAAAATAAGCCTAGCCCGTTACCAACTGGAAACAGTGACTTAAGACCGCCGGTCTTGTCCACTACCTTGCAGTAATGCGGTGGACAGGATCGGCGGTTTTCTTTTCTGAACCAGTCTCGTAGTTGTTACAGCGATAAGAATAGGTGTTGAAATACTCTTGACATGAGCTCGTCGTCAGGATATATAGCTTTGTACGCAAGTTCACGTAAAAGTTGTAGAGCTAGCAATAGCAGGTTACAATAAGGCTCGTCCGTTATAAACATGAAAATGTGTTCACAAATGCCGCCACTCAAACAGAGCGGCATTTTTCTTCCCCATCTCTTACCGAGTTTTACTTCAGTGTGCGAATAGACAACAATTGACAGAGGCAGTACTACCGTTTATAATTCGGACAATGCCTAAAGAGTTACCATGGAATAGAAAACAAAAGTTTAAGTTATTCTAAGGCCAGTCCGGAATCATCCTAAAAAGGAGTTATTGAACACCCGAAAGGGTGTTTTTTTGTTTTGTGAGACTTATTTATCCCGAAACTATTGTGTTACTGAAGCAACCGCAGATTGACATGCGTGATTTAACATTCTATAATTGCACAAACGCCTCCCATTCAGGGGAGATTTCGAGCTAGGCATAGCAAGTGAAATTAAGGCTGGTCCATTAACACCTTGAAAAAGGGAACAATAAGGCCTCCCTTTAGGGGGGGCCTTTTTTATTGATGAAAAGCAATCCCTCGTGAAGTAACTCAATAGTGTTCTCTGGTATCGTATTGACATAAGTCGTATTCAAAGATATAATATAGGTACAGTAAGTCGACCGACCGTTTTTCAGATTTGGAAACAAAACGTTGAAAAAAGGCAAGTCCGTTATGAACGCGAAAGCGTGCGAAAAAACCCGCTTCGGCGGGTTTTTTTATAGTTTCCAGAGGACCTTCACGGATAAAATAGATTACAGTTCTCGTCGTAGTATTGACAGTTGTGTTATCCGGCCATATAATATCTCTGTGAATACACTGCCCGTTTCGATGTAGATGTAGAAATACAAGGTTACATTAAGGCCCGTCCGTAATCAACTTGAAGAAGTGTTCCATCGGGTCCGAATTTTCGGACCTTTTCTCCGCATTACAATCAGCAGTCAGAACTTTTACGAAGAATAGTGGTCGCTCAACCTTTTGACAGTGTGCTAAAATTTGTCTATAATGAGTACCAGCGCGTCTTTTTGCATGAGTTGGAGAGCAAGACATTGCAAGTTCCAATAAGGCGTGTCCGATAAAAGCTTGAGAAAGCAAAGTAATACAAAACAGGCCCAGGCGGCCTGTTTTGTCTTTTTAATGTCCGTAGATAATAGAATAAGGTGCCCTCAGATTGTTGGAAGCGACTTTTATTGACAGCATCTGCTTTGTCACCTATAATTCAATGTGTGGTTGTGCTCTGGCAGAATCTGAGAGCCAAAAATGGCAAGTTCAGATAAGGCCAGACCGTTACCAGCTTAAATAAGCGATCCTAAAGCCCCGAATTTTTTATAAATTCGGGGCTTTTTTACTAGAGTATCGTGAAAACCTTTATTACCACACTCTGAACTGTAGGACGGGATTTTTGACAGACCTTATCTACATGGTTATAATCTGAATCAGGTTAGCAGTTCGAGAGTGCTTCAGATCCAGAAATGGAAAGTTGAAGTGAGGCAGGTCCGGTAGCAACTCGAAAGAGTGAGAAAAGAGGGGAGCGGGAAACCGCTCCCCTTTTTTCGTTTTATCGTATTCGTCACACCAGATTGGCGTAAGAAGTCGCTATTGAAACTATTTGACACTTTGCACATGTCCCGTTATAATCATGATCAGGCTAATCACTCGTAGAACATTTTGGCGTCGAAAGACGAAGTAAAATGAAGGCGAGACCGATATCAACTGGAAGCAGTGTCTGGTAGTCCTGGTAAGACGCGAACAGCGTCGCATCAGGCATATTGCCAACTAGCTGAATAAGCACTGTTGATAATCGCAATCTGTCTCTTCGTGAAAAGTAGCTTGACACGGATCTTCGCTGAACGTATAATGAGAAATACTGTACTAAAGTCACTTAGTTTTGGACCTAGAAATAGGAAGTCAAAATAAGGCTGGACCGACATGTAATCGAAAGATTTAGTCAAAAGCCTCCGGTCGGAGGCTTTTGACTTTCGTGAACGACACTACTATTTCTTACGAGATACTTATTCTGGAAGCAACGGTTTGACACAGCCCAGCCGGAGAGTATAATCCTATTATTAAACGCATCATAAAAATCTTACCGAACTAGGAATAGTAAGTGGTAAGAAGGCCTGACCGTAATAAGCCTGAAAAGGCGACCAAAAAGGGGGGATTTTATCTCCCCTTTAATTTTTCAAAGGTGGTATTTATTACGCAGACAACTCCCTGAGAACGGTTTTCAATCTAGAAATCATCCTTAGCGAAAGCTAAGGATTTTTTTTATCTGTTACACTGCGC'    #sequence = 'tgcatgatctacgtgcgtcacatgcaCgCGtacCAACCAGAAAAGCCAACCTGCGGGTTGGCTTTTTTATGCAAATAGGATCCATATCTCTTAGCGGCCGCAGTTTTACTTCAGTGTGCGAATAGACAACAATTGACAGAGGCAGTACTACCGTTTATAATTCGGACAGTCTCTTTTTTCTGTATCGGGAATAGAAAACAAAAGTTTAAGTTATTCTAAGGCCAGTCCGGAATCATCCTAAAAAGGAGTTATTGAACACCCGAAAGGGTGTTTTTTTGTTTTGTGAGACTTATTTATCCCGAAACTATTGTGTTACTGAAGCAACCGCAGATTGACATGCGTGATTTAACATTCTATAATTGCACATGGCTTTCCAATAATGCAGGGATTTCGAGCTAGGCATAGCAAGTGAAATTAAGGCTGGTCCATTAACACCTTGAAAAAGGGAACAATAAGGCCTCCCTTTAGGGGGGGCCTTTTTTATTGATGAAAAGCAATCCCTCGTGAAGTAACTCAATAGTGTTCTCTGGTATCGTATTGACATAAGTCGTATTCAAAGATATAATATAGGTTGAACAATGTTTATTCATCCTTTTCAGATTTGGAAACAAAACGTTGAAAAAAGGCAAGTCCGTTATGAACGCGAAAGCGTGCGAAAAAACCCGCTTCGGCGGGTTTTTTTATAGGAGATTACTTTACAGAAGACTCACTTATTTCACGGAACTGGTGCTGACAATTGACAGTTGTGTTATCCGGCCATATAATATCTCTGCACCCATTCCCGCGAAACGGATGTAGATGTAGAAATACAAGGTTACATTAAGGCCCGTCCGTAATCAACTTGAAGAAGTGTTCCATCGGGTCCGAATTTTCGGACCTTTTCTCCGCATTACAATCAGCAGTCAGAACTTTTACGAAGAATAGTGGTCGCTCAACCTTTTGACAGTGTGCTAAAATTTGTCTATAATGAGTACATCGTCGACATATTGCTGCCGTTGGAGAGCAAGACATTGCAAGTTCCAATAAGGCGTGTCCGATAAAAGCTTGAGAAAGCAAAGTAATACAAAACAGGCCCAGGCGGCCTGTTTTGTCTTTTTAATGTCCGTAGATAATAGAATAAGGTGCCCTCAGATTGTTGGAAGCGACTTTTATTGACAGCATCTGCTTTGTCACCTATAATTCAATGAAAAATCCATAATAAATATAATCTGAGAGCCAAAAATGGCAAGTTCAGATAAGGCCAGACCGTTACCAGCTTAAATAAGCGATCCTAAAGCCCCGAATTTTTTATAAATTCGGGGCTTTTTTACTAGAGTATCGTGAAAACCTTTATTACCACACTCTGAACTGTAGGACGGGATTTTTGACAGACCTTATCTACATGGTTATAATCTGAATGCTCTTGTTATAAATGGGAAGCTTCAGATCCAGAAATGGAAAGTTGAAGTGAGGCAGGTCCGGTAGCAACTCGAAAGAGTGAGAAAAGAGGGGAGCGGGAAACCGCTCCCCTTTTTTCGTTTTATCGTATTCGTCACACCAGATTGGCGTAAGAAGTCGCTATTGAAACTATTTGACACTTTGCACATGTCCCGTTATAATCATGATAAACAGATAAGCCTGAATTACATTTTGGCGTCGAAAGACGAAGTAAAATGAAGGCGAGACCGATATCAACTGGAAGCAGTGTCTGGTAGTCCTGGTAAGACGCGAACAGCGTCGCATCAGGCATATTGCCAACTAGCTGAATAAGCACTGTTGATAATCGCAATCTGTCTCTTCGTGAAAAGTAGCTTGACACGGATCTTCGCTGAACGTATAATGAGAAACTTTTTCATTAAAGCAATCAGTTTTGGACCTAGAAATAGGAAGTCAAAATAAGGCTGGACCGACATGTAATCGAAAGATTTAGTCAAAAGCCTCCGGTCGGAGGCTTTTGACTTTCGTGAACGACACTACTATTTCTTACGAGATACTTATTCTGGAAGCAACGGTTTGACACAGCCCAGCCGGAGAGTATAATCCTATTaatcattattttctcaaatgCTTACCGAACTAGGAATAGTAAGTGGTAAGAAGGCCTGACCGTAATAAGCCTGAAAAGGCGACCAAAAAGGGGGGATTTTATCTCCCCTTTAATTTTTCAAAGGTGGTATTTATTACGCAGACAACTCCCTGAGAACGGTTTTCAATCTAGAAATCATCCTTAGCGAAAGCTAAGGATTTTTTTTATCTGTTACACTGCGCcatatgctcagattcagtagaccgctgttg'

    ### main fxns
    assess = seq_assess(sequence)
    assess.load(sequence)
    assess.run()
    print assess.bad_seqs
    
