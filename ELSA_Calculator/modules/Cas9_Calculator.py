import scipy.io
import math
import operator
from time import time
from Bio import SeqIO
from Bio.Seq import *
import csv
import dill
from Mer_Maker import mers

import sys, os
cwd = os.path.dirname(os.path.abspath(__file__))

def identifyNucleotidePositionsOfMers(full_sequence, length=10):
    """Saves list of nucleotide positions in genome that all match a unique N-mer
    sequence. Counting begins at __ending_ of MER.
    
    Usage:   genomePositionsAtMers[mer_sequence] is a list of nucleotide positions 
    within the inputted fullSequence that match mer_sequence
             If mer_sequence ends in a PAM site, then this can be used to match 
             the first N-3 nt of a guide strand plus a PAM site sequence.
    """

    #Create a list of all possible N-mers
    all_possible_mers = mers(length)
    
    #Search through the genome and add nucleotide positions for match to an N-mer
    positions_at_mers = {}
    for mer in all_possible_mers:
        positions_at_mers[mer] = []
    
    # print "Number of Mers: ", len(positions_at_mers.keys())
    counter = 0
    while counter < (len(full_sequence)-length):
        word = full_sequence[counter : counter+length]
        positions_at_mers[word].append(counter+length)
        counter += 1
    return positions_at_mers

def identifyTargetSequencesMatchingPAM(PAM_seq, positions_at_mers, full_sequence, 
                                        target_sequence_length=20):
    """ Generates a list of target nucleotide sequences and corresponding nt 
    positions for an inputted sequence that matched the PAM_seq.
        Uses the positionsAtMers dictionary to accelerate the identification. 
        Good for large genomes.
        
    Usage:  listOfTargets = identifyTargetSequencesMatchingPAM('CGG', 
                                positionsAtMers, genome_sequence)
    """
    target_sequence_list = []
    all_mers = positions_at_mers.keys()
    mer_length = len(all_mers[0])
    list_of_mers_with_PAM = [mer + PAM_seq for mer in mers(mer_length - len(PAM_seq))]
    for mer_with_PAM in list_of_mers_with_PAM:
        nt_list = positions_at_mers[mer_with_PAM]
        for nt in nt_list:
            begin = nt-target_sequence_length - len(PAM_seq)
            end = nt - len(PAM_seq)
            if begin > 0 and end < len(full_sequence): #Does not account for circular DNAs
                target_sequence = full_sequence[begin : end]
                target_sequence_list.append((target_sequence, nt))
    return target_sequence_list
    
class sgRNA(object):

    def __init__(self, guideSequence, Cas9Calculator):
    
        self.guideSequence = guideSequence
        self.Cas9Calculator = Cas9Calculator
        
        self.partition_function = 1
        self.targetSequenceEnergetics = {}
        
        self.debug = False
        
    def run(self, sequence=None):
    
        begin_time = time()
        dG_best=float(1e12)
        dG_source=None
        dG_pos=None
        targetDictionary = self.Cas9Calculator.targetDictionary
        for source in targetDictionary.keys():
            for (fullPAM,targets) in targetDictionary[source].iteritems():
                self.targetSequenceEnergetics[source] = {}
                #print fullPAM
                dG_PAM = Cas9Calculator.calc_dG_PAM(fullPAM)
                dG_supercoiling = Cas9Calculator.calc_dG_supercoiling(sigmaInitial = -0.05, targetSequence = 20 * "N")  #only cares about length of sequence
                for (targetSequence,targetPosition) in targets:
                    dG_exchange = Cas9Calculator.calc_dG_exchange(self.guideSequence,targetSequence)
                    dG_target = dG_PAM + dG_supercoiling + dG_exchange
                    if sequence:
                        if targetSequence==sequence:
                            print "true sequence found"
                            dG_best=float(dG_target)
                            dG_source=str(source)
                            dG_pos=int(targetPosition)
                    else:
                        if dG_target<dG_best:
                            print "overtaken"
                            dG_best=float(dG_target)
                            dG_source=str(source)
                            dG_pos=int(targetPosition)
                    self.targetSequenceEnergetics[source][targetPosition] = {'sequence' : targetSequence, 'dG_PAM' : dG_PAM, 'full_PAM' : fullPAM, 'dG_exchange' : dG_exchange, 'dG_supercoiling' : dG_supercoiling, 'dG_target' : dG_target}
                    self.partition_function += math.exp(-dG_target / self.Cas9Calculator.RT)
                    
                    if self.debug:
                        print "targetSequence : ", targetSequence
                        print "fullPAM: " , fullPAM
                        print "dG_PAM: ", dG_PAM
                        print "dG_supercoiling: ", dG_supercoiling
                        print "dG_exchange: ", dG_exchange
                        print "dG_target: ", dG_target
                        print "Partition function (so far): ", self.partition_function
        
        end_time = time()
        print "Elapsed Time: ", end_time - begin_time
        bestDict=self.targetSequenceEnergetics[dG_source][dG_pos]
        return math.exp(-bestDict['dG_target'] / self.Cas9Calculator.RT)/self.partition_function, self.partition_function, bestDict['dG_target'], bestDict['dG_PAM'], bestDict['dG_supercoiling'], bestDict['dG_exchange']
    
    def exportAsDill(self):
    
        handle = open('sgRNA_%s.dill' % self.guideSequence,'wb')
        dill.dump(self, handle, -1)
        handle.close()
        
class clCas9Calculator(object):

    def __init__(self,seqRecords, quickmode=False, ModelName='InvitroModel.mat'):
        
        self.quickmode=quickmode
        self.ModelName=ModelName
        data=scipy.io.loadmat(os.path.join(cwd, self.ModelName))
        self.weights=data['w1']
        self.decNN=data['decNN']
        self.RT = 0.61597        
        
        # the PAMs with the highest dG, ignoring other PAM sequences by setting their dG to 0
        self.PAM_energy={'GGA':-9.8,'GGT':-10,'GGC':-10,'GGG':-9.9,'CGG':-8.1,'TGG':-7.8,'AGG':-8.1,'AGC':-8.1,'AGT':-8.1,'AGA':-7.9,'GCT':-7.1,'GCG':-6.9,'ATT':-7.2,'ATC':-6.4,'TTT':-7.6,'TTG':-6.8,'GTA':-7.4,'GTT':-7.9,'GTG':-7.7,'AAT':-7,'AAG':-7,'TAT':-7.2,'TAG':-7.2,'GAA':-7.2,'GAT':-7.3,'GAC':-7.2,'GAG':-7.3}
    
        self.initTargetFinder(seqRecords)
    
    def returnAllPAMs(self):
        
        for (PAMpart,energies) in sorted(self.PAM_energy.items(), key = lambda x: x[1]):  #PAMpart will be 'GGT'
            for nt in ('A','G','C','T'):        #nt + PAMpart will be all possible 'NGGT'
                yield nt + PAMpart
    
    def initTargetFinder(self, seqRecords):
    
        targetDictionary = {}        
        for record in seqRecords:
            fullSequence = str(record.seq)+str(record.seq.reverse_complement())
            positionsAtMers = identifyNucleotidePositionsOfMers(fullSequence, length = 10)
            targetDictionary[record.id] = {}
            targetSequenceList = []
            for fullPAM in self.returnAllPAMs():
                targetSequenceList = identifyTargetSequencesMatchingPAM(fullPAM, positionsAtMers, fullSequence)
                targetDictionary[record.id][fullPAM] = targetSequenceList
            print "Loaded %s genome into Cas9 Calculator. %s nucleotides long." % (record.id, len(record.seq))
        self.targetDictionary = targetDictionary
        
    def printModelInfo(self):
        m=0
        s=0
        negative_val=0
        for i,l in enumerate(self.decNN):
            for j,e in enumerate(l):
                if float(e)<0:
                    negative_val+=1
                if i!=j:
                    s+=float(e)
                    m+=1
        meanNN=float(s)/float(m)
        
        sw=0
        for w in self.weights:
            sw+=w
        
        meanw=sw/len(self.weights)
        print 'average mismatchc energy: ', meanNN
        print 'average weight:', meanw
        print 'number of negative energies: ', negative_val
        
    def Calc_Exchange_Energy(self,crRNA,targetSeq):
        nt_pos={'A':0,'T':1,'C':2,'G':3,'a':0,'t':1,'c':2,'g':3}
        dG=0
        RNA=''
        DNA=''
        for i in range(0,len(crRNA)):
            if i>0:
                RNA=crRNA[(i-1):(i+1)]
                DNA=targetSeq[(i-1):(i+1)]
                RNA_index=nt_pos[RNA[0]]+4*nt_pos[RNA[1]]
                DNA_index=nt_pos[DNA[0]]+4*nt_pos[DNA[1]]
                
                dG1=float(self.decNN[RNA_index][DNA_index])
                if abs(dG1-0.000015)<1e-6:
                    dG1=10000
                    dG1=2.3 # during model identification, I set the value of every unknown dG to 0.000015 (if I did not find a value for it)
                    
                    
                pos=20-i
                w1=float(self.weights[pos])
                #print 'b1',RNA[0],RNA[1],DNA[0],DNA[1],RNA_index, DNA_index, pos,dG1, w1
            else:
                w1=0
                dG1=0
            if i<(len(crRNA)-1):
                RNA2=crRNA[i:(i+2)]
                DNA2=targetSeq[i:(i+2)]         
                RNA_index=nt_pos[RNA2[0]]+4*nt_pos[RNA2[1]]
                DNA_index=nt_pos[DNA2[0]]+4*nt_pos[DNA2[1]]
                dG2=float(self.decNN[RNA_index][DNA_index]) 
                if abs(dG2-0.000015)<1e-6:
                    dG2=10000
                    dG2=2.3 # during model identification, I set the value of every unknown dG to 0.000015 (if I did not find a value for it)
                
                pos=20-i-1
                w2=float(self.weights[pos])
                #print 'b2',RNA2[0],RNA2[1],DNA2[0],DNA2[1],RNA_index, DNA_index, pos,dG2, w2
            else:
                w2=0
                dG2=0 
            dG+=w1*dG1+w2*dG2
        return float(dG)

    def QuickCalc_Exchange_Energy(self,crRNA,TargetSeq):
        nt_pos={'A':0,'T':1,'C':2,'G':3,'a':0,'t':1,'c':2,'g':3}
        dG=0
        RNA=''
        DNA=''
        self.nt_mismatch_in_first8=0
        for i in range(0,len(crRNA)):
            pos=20-i
            w1=self.weights[pos]
            if nt_pos[crRNA[i]]==nt_pos[TargetSeq[i]]:
                dG1=0
            else:
                # using a bioinformatics search approach to find sequences with up to x mismatches
                dG1=2.3 # kcal/mol
                if pos<=8:
                    self.nt_mismatch_in_first8=self.nt_mismatch_in_first8+1
            dG+=w1*dG1
        return float(dG)
        
    def calc_dG_PAM(self, PAM_full_seq):
    
        #PAM sequence is 5' - N xxx N - 3' where the energy of xxx is listed below. A normal PAM of 'NGG' with 'TC' afterwards would be listed as 'GGT'
        key = PAM_full_seq[1:4]
        if key in self.PAM_energy:
            return self.PAM_energy[key]
        else:
            return 0.0
        
        # acceptedPAMList=PAM_dic_energy.keys()
        # self.dG_PAM_List=[]
        # self.WarningPAM_List=[]
        # PAMsize=len(self.PAM)
        # for target in self.sequence_list:
            # tPAM=target[-(PAMsize):-1]+target[-1]
            # if tPAM in acceptedPAMList:
                # dGPAM=PAM_dic_energy[tPAM]
                # warning=''
            # else:
                # dGPAM=0
                # warning='N.B'             
            # self.dG_PAM_List.append(dGPAM)
            # self.WarningPAM_List.append(warning)
           
    def calc_dG_exchange(self, guideSequence, targetSequence):
        self.nt_mismatch_in_first8_list=[]
        if self.quickmode:
            solverfunc=self.QuickCalc_Exchange_Energy           
        else:
            solverfunc=self.Calc_Exchange_Energy
        
        dG_exchange = solverfunc(guideSequence,targetSequence)
        
        return dG_exchange
        
    def calc_dG_supercoiling(self, sigmaInitial, targetSequence):
    
        
        sigmaFinal = -0.08
        dG_supercoiling = 10.0 * len(targetSequence) * self.RT * (sigmaFinal**2 - sigmaInitial**2)
        return dG_supercoiling
                        
if __name__ == "__main__":

    # guideSequence = 'TACGTACACAAGAGCTCTAG'      
    # Cas9Calculator=clCas9Calculator(['Neural_Network/NC_000913_3.gb'])
    # sgRNA1 = sgRNA(guideSequence, Cas9Calculator)
    # sgRNA1.run()
    # sgRNA1.exportAsDill()
    
    #PAM='GGA' # NGGA
    #sequence_list=['AGTCCTCATCTCCCTCAAGCCGGA','AGTCCTCATCTCCCTCAAGTCGGA','AGTCCTCATCTCCCTCATGCCGGA']  # list of all potential on- and off-targets
    #Cas9Calculator=clCas9Calculator(quickmode=True) # using quick approach
    #Cas9Calculator=clCas9Calculator(quickmode=False) # using Invitro or complete model
    #Cas9Calculator=clCas9Calculator(quickmode=False,cModelName='All_dataModel.mat') # using Invitro or complete model
    #Cas9Calculator.loadData(sequence_list,crRNAseq,PAM,True)
    #Cas9Calculator.calcTarget_energy()
    #Cas9Calculator.export_dG() # in an excel file
    #print Cas9Calculator.dG_total_List

    # For the sgRNA project
    # guides = [
    #     'GTCGCTCGAGACACGAAAAG',
    #     'CTACAGAGCGATAGTTCGTA',
    #     'AGCTTCTGAGTAGTCTAAGT',
    #     'GCTCTCTAGTCATCTCACAA',
    #     'AGCACACACGCGCTAGAGAC',
    #     'ACGAGACTAGCTAGCGAATT',
    #     'CTAGTCACGTCTCTAAGCAC',
    #     'TAGCTCTGTCGACGTCTAGC',
    #     'TAGCAGCTTGCTGAGAGCTC',
    #     'ATAACACAGATTCGCTTAGA',
    #     'AAGACTAGATCTTAGCTAGT',
    #     'CTACGCACGAGCTAGCACAG',
    #     'TAACACTCGAAAGTCGCTAG',
    #     'GCTAGATTGCGTACGTCTAA',
    #     'GCACACTAGCAGACTTCTCT',
    #     'TGTATAGAGCGCAAGACTCT',
    #     'GACAAGACACTCTAGCGACG',
    #     'ATACTCGACTTCGAAGCTTG',
    #     'GAGTCGTTCGCGAGCTAGCT',
    #     'CTCTAATCGCGCTTCTAGAG',
    #     'CATATACTACATACTCTATT',
    #     'GTATAATAGAAAACTAGACA',
    #     'TCTAAGCTCGCTGTAGACTA',
    #     'TCACGACTAGAGCTAGCACG',
    #     'ACTAGTTGTCGTAACATCTA',
    #     'TCACTAGTTGTCTGAGCTCT',
    #     'CTAGATCGTTTAACTAAGTA',
    #     'AATTAGTACTTCTTATAGTT',
    #     'ATCTAGTTCGCTATACACGT',
    #     'GATACAGTCAGTCGCTTGAC',
    #     'AAGTCTAGTACGCTTACATA',
    #     'TCGAAGTATATACTAACGAG',
    #     'CTCATAGTAGTACATACGAA',
    #     'AGACAGAGTCGCTATCTAAA',
    #     'TGTCTAACTACACTTGACAG',
    #     'GAGTCGTTCGCGAGCTAGCT',
    #     'CTCTAATCGCGCTTCTAGAG',
    #     'CATATACTACATACTCTATT',
    #     'GTATAATAGAAAACTAGACA',
    #     'TCTAAGCTCGCTGTAGACTA',
    #     'TCACGACTAGAGCTAGCACG',
    #     'GTCGCTCGAGACACGAAAAG',
    #     'TCACTAGTTGTCTGAGCTCT',
    #     'ACTTACAGTATAAATTATAG',
    #     'CTAGATCGTTTAACTAAGTA',
    #     'AATTAGTACTTCTTATAGTT',
    #     'ATCTAGTTCGCTATACACGT',
    #     'GATACAGTCAGTCGCTTGAC',
    #     'AAGTCTAGTACGCTTACATA',
    #     'TCGAAGTATATACTAACGAG',
    #     'CTCATAGTAGTACATACGAA',
    #     'AGACAGAGTCGCTATCTAAA',
    #     'TGTCTAACTACACTTGACAG']
    # guides = ['AAGTACCTATAATTGATACG',
    #         'ATAACGCGCATCTTTCATGA',
    #         'GCACCCATTCCCGCGAAACG',
    #         'GTTTAGTTCATCTGACGGAG',
    #         'AGGGTGTTCCGGAGACCTGG',
    #         'GAAGAAGCTGCCCATAATCG']
    guides = [
    'CCATCTCCTGAATGTGATAA'
    ]
    f = open(os.path.join(cwd, 'check_guides_genome.txt'),'w')
    # for seq in guides:
    #     handle = open('psc101_target3_riboj_mrfp1.gb','r')
    #     records = SeqIO.parse(handle,"genbank")
    #     record = records.next()
    #     handle.close()
    #     print record[1005:1025]
    #     revcomp=str(Seq(seq).reverse_complement())
    #     new_record=record[0:1005]+revcomp+record[1025:-1]
    #     handle = open('psc101_target3_riboj_mrfp1.gb','w')
    #     new=[new_record]
    #     SeqIO.write(new,handle,"genbank")
    #     handle.close()
    #     Cas9Calculator=clCas9Calculator(['NC_000913.gb','psc101_target3_riboj_mrfp1.gb'],quickmode=True)
    Cas9Calculator=clCas9Calculator(['/home/ubuntu/Genomes/refseq/bacteria/GCF_000005845.2/GCF_000005845.2_ASM584v2_genomic.gbff.gz'],quickmode=True)
    for seq in guides:    
        print seq
        sgRNAobj = sgRNA(seq, Cas9Calculator)
        x = sgRNAobj.run(seq)
        f.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(x[0],x[1],x[2],x[3],x[4],x[5]))
    f.close()
    