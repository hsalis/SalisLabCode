import sys,traceback,types,time, cPickle, os, re, itertools, time, random, math, copy_reg, json, StringIO
from PyVRNA import PyVRNA
from Bio.Seq import Seq
import numpy as np
from Bio import pairwise2,SeqIO,SeqFeature
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.SeqFeature import SeqFeature,FeatureLocation
from deap import base
from deap import creator
from deap import tools
from seq_assess import seq_assess
from string import maketrans
from Cas9_Calculator import clCas9Calculator
from mpi4py import MPI
import cPickle as pickle
import pickle as pkl
from RepeatFinder import KMP

model = PyVRNA(parameter_file="rna_turner1999.par",dangles=0)
cwd = os.path.dirname(os.path.abspath(__file__))

#reads in part sequences for ELSAs from FASTA files
backbones=[str(record.seq) for record in SeqIO.parse(cwd + "/../parts/backbones.fa","fasta")]
terminators = [str(record.seq) for record in SeqIO.parse(cwd + "/../parts/terminators.fa","fasta")]
handles_all=[str(record.seq) for record in SeqIO.parse(cwd + "/../parts/handles.fa","fasta")]       
spacers=[str(record.seq) for record in SeqIO.parse(cwd + "/../parts/spacers.fa","fasta")]
promoters=[str(record.seq) for record in SeqIO.parse(cwd + "/../parts/promoters.fa","fasta")]

#gets names of parts from FASTA files
part_names = {"backbones":[record.name for record in SeqIO.parse(cwd + "/../parts/backbones.fa","fasta")],
              "terminators":[record.name for record in SeqIO.parse(cwd + "/../parts/terminators.fa","fasta")],
              "handles":[record.name for record in SeqIO.parse(cwd + "/../parts/handles.fa","fasta")],
              "spacers":[record.name for record in SeqIO.parse(cwd + "/../parts/spacers.fa","fasta")],
              "promoters":[record.name for record in SeqIO.parse(cwd + "/../parts/promoters.fa","fasta")]}

part_seqs = {"backbones": backbones,
              "terminators": terminators,
              "handles":handles_all,
              "spacers":spacers,
              "promoters":promoters}
# indexes for successfully synthesized ELSAs, used in quickmode
'''
DO NOT EDIT
'''
quick_small = {"promoters":[21,14,24,19,9,18,10,15],"handles":[7,18,1,16,21,19,13,16],"spacers":[49,43,9,6,5,23,52,3]}
quick_medium = {"promoters":[25,14,15,16,17,18,19,20,21,22],"handles":[23,10,17,18,20,21,22,24,19,26],"spacers":[14,15,33,17,18,19,20,21,22,58]}
quick_large = {"promoters":[1,2,25,14,15,16,17,18,19,20,21,22],"handles":[3,9,23,10,17,18,20,21,22,24,19,26],"spacers":[23,13,14,15,33,17,18,19,20,21,22,58]}
quick_max = {"promoters":[24,14,19,21,10,9,18,26,15,8,27,1,5,2,28,12,3,23,0,4],"handles":[1,18,6,7,13,16,19,20,21,25,15,26,24,23,3,22,10,11,9,17],"spacers":[9,43,49,6,52,23,5,12,3,2,13,45,16,34,21,1,46,18,56,44]}
'''
DO NOT EDIT
'''

class ELSA_Calculator(object):
    def __init__(self,Cas9Calculator):
        '''
        initializes main ELSA calculator class
        '''

        # initalizing Cas9Calculator
        self.Cas9Calculator=Cas9Calculator
        self.RT=Cas9Calculator.RT
        self.rc = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}

        #initializes PyVRNA
        self.model = model

        # initializes GA params
        self.POP_SIZE = 5
        self.NGEN = 5
        self.CXPB = 0.65
        self.MUTPB = 0.25
        self.NUM_ELITES = 2  

        # defaults for guide sequences, handle selection
        self.prefeed=None
        self.quickmode=False

    def assemble_handles(self,H, backbones, terminators):
        '''
        builds new handles using terminator and backbone libraries
        '''

        # converts backbone seqs to RNA for folding
        true_backbones=[]

        for i in range(len(backbones)):
            temp=backbones[i].replace("T","U")
            true_backbones.append(temp)

        # builds list of all possible backbone/terminator pairs
        tmp=list(itertools.product(true_backbones,terminators))
        sgRNA_list=[]

        for i in tmp:
            sgRNA_list.append(i[0]+i[1])
        
        # filters out those pairs that disrupt SL2 and R/AR
        tmp_list=[]
        
        constraint="((((((..((((....))))....))))))..................((((....))))."
        for i in sgRNA_list:
            structure=model.RNAfold(i).structure
            #print test
            test="".join([structure[j] for j in range(61)])

            for j in range(61):
                if constraint[j]=='(' and test[j]!="(":
                    break
                elif constraint[j]==')' and test[j]!=")": 
                    break
                elif j ==7 and test[j]!=".": 
                    break
                elif j in range(12,16) and test[j]!=".": 
                    break
                elif j in range(20,23) and test[j]!=".": 
                    break
                else:
                    tmp_list.append(i)
                    #print 'added'

        # converts back to DNA 
        new_list=[]
        
        for i in range(len(tmp_list)):
            temp=tmp_list[i].replace("U","T")
            new_list.append(temp)

        kmers_last=[]
        sgRNAs=[]
        for n in new_list:
            kmers = [n[i:i+H]for i in range(len(n)-H)]
            kmers_new=set(kmers_last+kmers)
            if len(kmers_new)<len(kmers)+len(kmers_last):
                continue
            else:
                kmers_last=list(kmers_new)
                sgRNAs.append(n)

        return sgRNAs

    def find_guides(self,target_regions,array_length,strand,guide_distribution,prefeed=None,pool=None):
        '''
        finds guides in genome to use in GA/ extracts data on preselected guides 
        '''

        # parallelization
        if pool is not None:
            pool.start()
            if pool.is_master():
                map_function = pool.map
                print 'INFO: Using MPI Pool'

        else:
            map_function = map
            print 'INFO: Using serial mode'

        if pool is None or pool.is_master():
            # set up Cas9Calculator
            guide_list = []
            total_guide_dict = {'warnings':[]}
            dist_error=False

            # processing and flags for preselected guides  
            if prefeed:
                true_dist=None
                g_flag=True
                guides = prefeed

                # processing by Cas9Calculator
                inputs = [(g, self.Cas9Calculator) for g in guides]
                outputs = map_function(partition_function, inputs)

                # populating information dictionaries
                for i in range(len(guides)):
                    total_guide_dict[guides[i]] = {}
                    dG_target, partition_fxn,on_targets,off_targets = outputs[i]
                    total_guide_dict[guides[i]]['dG_target'] = dG_target
                    total_guide_dict[guides[i]]['partition_fxn'] = partition_fxn
                    total_guide_dict[guides[i]]['position'] = 0
                    total_guide_dict[guides[i]]['off_targets'] = off_targets
                    total_guide_dict[guides[i]]['on_targets']=on_targets


                guide_list = guides

            # if no preselected guides, builds pool of guides for each locus
            else:
                true_dist=list(guide_distribution)
                g_flag=False
                print "target_regions: ", target_regions
                for r in range(len(target_regions)):
                    print "INFO: target region: %i"%(r)
                    
                    #guide selection based on strand selection
                    guide_dict = {}

                    if strand=="coding":
                        print target_regions[r]

                        # if targeting coding region, guide is reverse complement
                        for i in range(len(target_regions[r])-23):

                            if target_regions[r][i:i+2]=='CC':
                                guide = ''.join(reversed([self.rc.get(j.upper(),j.upper()) for j in target_regions[r][i+3:i+23]]))
                                print target_regions[r][i:i+23], guide
                                guide_dict[guide]={}
                                guide_dict[guide]['position']=i

                    elif strand=="template":
                        print target_regions[r]

                        # if targeting template region, guide is coding sequence
                        for i in range(len(target_regions[r])-23):
                            if target_regions[r][i+21:i+23]=='GG':
                                guide = target_regions[r][i:i+20]
                                print target_regions[r][i:i+23], guide
                                guide_dict[guide]={}
                                guide_dict[guide]['position']=i+23
                    else:
                        #if targeting both strands, guides from both strands are used
                        print target_regions[r]
                        for i in range(len(target_regions[r])-23):

                            if target_regions[r][i:i+2]=='CC':
                                guide = ''.join(reversed([self.rc.get(j.upper(),j.upper()) for j in target_regions[r][i+3:i+23]]))
                                print target_regions[r][i:i+23], guide
                                guide_dict[guide]={}
                                guide_dict[guide]['position']=i

                            if target_regions[r][i+21:i+23]=='GG':
                                guide = target_regions[r][i:i+20]
                                print target_regions[r][i:i+23], guide
                                guide_dict[guide]={}
                                guide_dict[guide]['position']=i+23

                    # isolates guides for direct handling by GA
                    guides=[i for i in guide_dict.keys()]

                    # filter if region has insufficient target regions
                    if len(guides)<guide_distribution[r]:
                        if len(guides)>0:
                            true_dist[r]=len(guides)
                            warning="Total number of guides detected ({}) in target_region {} less than number of guides requested ({}), ELSA design has been modified to compensate.".format(len(guides),r,guide_distribution[r])
                            total_guide_dict['warnings'].append(warning)
                            print warning
                        else: #len(guides)==0
                            true_dist[r]=len(guides)
                            dist_error=True
                            warning="no guides found in target region {}".format(r)
                            total_guide_dict['warnings'] = warning
                            return guide_list, total_guide_dict, g_flag, true_dist, dist_error

                    # runs guides on Cas9Calculator
                    inputs=[(g,self.Cas9Calculator) for g in guides]
                    outputs=map_function(partition_function,inputs)

                    # populating information dictionaries                    
                    for i in range(len(guides)):
                        total_guide_dict[guides[i]]={}
                        dG_target, partition_fxn, on_targets,off_targets = outputs[i]
                        total_guide_dict[guides[i]]['dG_target']=dG_target
                        total_guide_dict[guides[i]]['partition_fxn']=partition_fxn
                        total_guide_dict[guides[i]]['position']=guide_dict[guides[i]]['position']
                        total_guide_dict[guides[i]]['off_targets']=off_targets
                        total_guide_dict[guides[i]]['on_targets']=on_targets

                    guide_list.append(guides)
        else:
            pass #workers see this
  

        return guide_list, total_guide_dict, g_flag, true_dist, dist_error

def _reduce_method(m):
    '''
    helper function to prevent classes from gumming up MPI
    '''
    if m.im_self is None:
        return getattr, (m.im_class, m.im_func.func_name)

    else:
        return getattr, (m.im_self, m.im_func.func_name)

def init_de_novo(individual, VLSA):
    '''
    assembles initial array designs
    '''
    parts=VLSA.new_parts
    values=[]

    for i in parts['order']:
        if VLSA.quickmode:
            print "INFO: Selecting Quick Design Mode for ELSA Calculator"
            
            # if guides are not preselected
            if i == "guides" and VLSA.g_flag == False:
                array_parts=[]

                # randomly selects guide indices for each target region to achieve desired coverage 
                indexes=[random.sample(xrange(len(parts['guides'][j])),VLSA.guide_distribution[j]) for j in range(len(VLSA.guide_distribution))]
                
                # reorganizes index information
                for j in range(len(VLSA.guide_distribution)):
                    for k in range(VLSA.guide_distribution[j]):
                        guide_inds=[j,indexes[j][k]]
                        array_parts.append(guide_inds)

                # shuffles guides post selection to distrubute targeting
                random.shuffle(array_parts)
                values.append(array_parts)

            #if guides are preselected, guides are randomized
            elif i == "guides" and VLSA.g_flag == True:
                array_parts=random.sample(xrange(len(parts[i])),VLSA.array_length)
                values.append(array_parts)

            # if not guides, parts are fixed
            else:
                array_parts = [j for j in range(VLSA.array_length)]
                values.append(array_parts)

        else:
            
            # if guides are not preselected
            if i == "guides" and VLSA.g_flag == False:
                array_parts=[]
                
                # randomly selects guide indices for each target region to achieve desired coverage 
                indexes=[random.sample(xrange(len(parts['guides'][j])),VLSA.guide_distribution[j]) for j in range(len(VLSA.guide_distribution))]
                
                # reorganizes index information
                for j in range(len(VLSA.guide_distribution)):
                    for k in range(VLSA.guide_distribution[j]):
                        guide_inds=[j,indexes[j][k]]
                        array_parts.append(guide_inds)

                # shuffles guides post selection to distrubute targeting
                random.shuffle(array_parts)
                values.append(array_parts)

            # if guides preselected or any other part,randomized
            else:
                array_parts=random.sample(xrange(len(parts[i])),VLSA.array_length)
                values.append(array_parts)

    # initializes individuals
    ind=individual(values)
    
    # sets initial fitness really high so  will change to minimize down
    ind.fitness.values = (1e12,1e12,1e12,1e12,1e12,1e12)
    return ind

def mutate_de_novo(individual,VLSA,indpb1=0.15,indpb2=0.30,indpb3=0.30,indpb4=0.05):
    '''
    replaces ELSA parts with non duplicate parts in the GA
    '''
    for i in VLSA.new_parts['order']:
        # rolls for mutation event
        if random.random() < indpb1:
            
            # if quickmode, only can edit guides
            if VLSA.quickmode:
                if i=="guides": 
                    
                    # if no preselected guides
                    if VLSA.g_flag == False:
                        array_parts=[]

                        # randomly selects guide indices for each target region to achieve desired coverage
                        indexes=[random.sample(xrange(len(VLSA.new_parts['guides'][j])),VLSA.guide_distribution[j]) for j in range(len(VLSA.guide_distribution))]
                        
                        # reorganizes index information
                        for j in range(len(VLSA.guide_distribution)):
                            for k in range(VLSA.guide_distribution[j]):
                                guide_inds=[j,indexes[j][k]]
                                array_parts.append(guide_inds)
                        
                        # shuffles guides post selection to distrubute targeting
                        random.shuffle(array_parts)
                        individual[VLSA.new_parts['order'].index(i)]=array_parts

                    # if preselected guides, randomized
                    else:
                        individual[VLSA.new_parts['order'].index(i)]=random.sample(xrange(len(VLSA.new_parts[i])),VLSA.array_length)

                else: continue
            # if not qucikmode, replaces all parts of a given type
            else:
                if i=="guides" and VLSA.g_flag == False:
                    array_parts=[]
                    indexes=[random.sample(xrange(len(VLSA.new_parts['guides'][j])),VLSA.guide_distribution[j]) for j in range(len(VLSA.guide_distribution))]
                    for j in range(len(VLSA.guide_distribution)):
                        for k in range(VLSA.guide_distribution[j]):
                            guide_inds=[j,indexes[j][k]]
                            array_parts.append(guide_inds)
                    
                    random.shuffle(array_parts)
                    individual[VLSA.new_parts['order'].index(i)]=array_parts
                
                else:
                    individual[VLSA.new_parts['order'].index(i)]=random.sample(xrange(len(VLSA.new_parts[i])),VLSA.array_length)
            
        # replaces one part of a given type
        if random.random() < indpb2 and i!="guides" and VLSA.quickmode == False:
            
            for j in xrange(VLSA.array_length):

                # select specific part to remove 
                if random.random() < indpb3:
                    replacement=random.randint(0,len(VLSA.new_parts[i])-1)
                    
                    #replace with new part
                    if replacement not in individual[VLSA.new_parts['order'].index(i)]:
                        individual[VLSA.new_parts['order'].index(i)][j]= replacement

    return individual      


def customSelect(new_population, population,ELSA):
    """
    NSGA2 selection enforces pareto optimality
    """

    elite_population = tools.selBest(new_population + population, k = ELSA.NUM_ELITES)
    selected_population = tools.selNSGA2(new_population + population, k = (len(population) - ELSA.NUM_ELITES) )
    return elite_population + selected_population

def evaluate_construct((individual,part_dict,array_length,arms,total_guide_dict,guide_distribution,g_flag,RT)):
    '''
    Evaluates the synthesis complexity of an array
    '''
    
    pickled = cPickle.dumps((individual,part_dict,array_length,arms,total_guide_dict,guide_distribution,g_flag,RT))
    #print "Size of Pickled (individual,part_dict,array_length,arms,total_guide_dict,guide_distribution,g_flag,RT): %s bytes" % sys.getsizeof(pickled)
    new_parts=dict(part_dict)

    #starts sequence w/ 5' homology arm 
    sequence=str(arms[0])
    canonical_structure='xxxxxxxxxxxxxxxxxxxx'
    dG_t=[]
    pf=[]
    p_list=[]
    off_targets=[]
    for i in range(array_length):

        for j in range(len(new_parts["order"])):

            # assembles sequences
            if new_parts["order"][j]=="guides" and g_flag == False:

                #retrieves guide sequence, dG_target information and partition function information
                dna=new_parts[new_parts["order"][j]][individual[j][i][0]][individual[j][i][1]]
                sequence+=dna
                dG_t.append(total_guide_dict[dna]["dG_target"])
                pf.append(total_guide_dict[dna]["partition_fxn"])
                position = total_guide_dict[dna]["position"]

                adjacent_guide_list = [new_parts[new_parts["order"][j]][individual[j][g][0]][individual[j][g][1]] for g in range(len(individual[j])) if individual[j][g][0] == individual[j][i][0]]

                d_x = sum([abs(position - total_guide_dict[p]["position"]) for p in adjacent_guide_list])
                p_list.append(d_x)
                off_targets.append(total_guide_dict[dna]["off_targets"])

            elif new_parts["order"][j]=="guides" and g_flag == True:

                #retrieves guide sequence, dG_target information and partition function information
                dna=new_parts[new_parts["order"][j]][individual[j][i]]
                sequence+=dna
                dG_t.append(total_guide_dict[dna]["dG_target"])
                pf.append(total_guide_dict[dna]["partition_fxn"])
                off_targets.append(total_guide_dict[dna]["off_targets"])
                p_list.append(0)

            else:
                #retrieves part sequence
                dna=new_parts[new_parts["order"][j]][individual[j][i]]
                sequence+=dna
        #sgRNA=new_parts["guides"][individual[new_parts["order"].index("guides")][i][0]][[individual[new_parts["order"].index("guides")][i][1]]+new_parts["handles"][individual[new_parts["order"].index("handles")][i]]
        #print sgRNA
        #con="".join([canonical_structure[j] if j<len(canonical_structure)  else "." for j in range(len(sgRNA))]) #32<j<40
        #dG_native=model.RNAfold(sgRNA).energy
        #dG_final=model.RNAfold(sgRNA,con).energy
        #print dG_final, dG_native
        #ddGs.append(dG_final-dG_native)

    # ends sequence w/ 5' homology arm 
    sequence+=str(arms[1])

    # starts sequence assessment helpers
    
    assess=seq_assess()
    assess.load(sequence)
    #print "debug"
    assess.run_ELSA_analysis()
    #print "rebug"

    # repeat score from seq assess, scaled
    repeats=int(assess.score[0])

    # haripin score from seq assess
    hairpins=int(assess.score[1])

    # nucleotide score from seq assess
    GC=int(assess.score[2])

    # converts partition fxn to linear term 
    PartFxn=math.log(sum(pf))
    
    PF_off = 0
    for o in range(len(off_targets)):
        for source,position in off_targets[o].items():
            for targets, data in position.items(): 
                PF_off += math.exp(-data['dG_target']/RT) 
    PartFxn_off=math.log(PF_off)
    dx_t = sum(p_list)/len(p_list)
    
    pickled = cPickle.dumps((repeats,hairpins,GC,dx_t,PartFxn,PartFxn_off))
    #print "Size of Pickled (repeats,hairpins,GC,dx_t,PartFxn,PartFxn_off): %s bytes" % sys.getsizeof(pickled)
    return (repeats,hairpins,GC,dx_t,PartFxn,PartFxn_off)
    
def partition_function((g,Cas9Calculator)):
    '''
    Cas9 Calculator parallel function, runs each guide calculation
    '''
    pickled = cPickle.dumps({'g' : g, 'Cas9Calculator' : Cas9Calculator})
    #print "Size of Pickled Dictionary containing g and Cas9Calculator: %s bytes" % sys.getsizeof(pickled)
    
    print "..."
    print g
    print "..."
    print Cas9Calculator.targetDictionary.keys()
    #initial value of partition fxn
    partition_function=1

    all_off_targets={}
    all_on_targets={}
    dG_offs=[(0,0,0,0,0,0,0,100000)]
    # choose genbank file
    for source in Cas9Calculator.targetDictionary.keys():
        print "Source: ", source
        # for genome in genbank file
        #best_off_targets[source]={}
        #all_on_targets[source]={}

        #Cas9Calculator.targetDictionary[source].iteritems()
        for (fullPAM,targets) in Cas9Calculator.targetDictionary[source].iteritems():
            #print fullPAM
            #print targets.keys()
            dG_PAM = Cas9Calculator.calc_dG_PAM(fullPAM)
            dG_supercoiling = Cas9Calculator.calc_dG_supercoiling(sigmaInitial = -0.05, targetSequence = 20 * "N")  #only cares about length of sequence
            
            # for target sites with each PAM
            for (targetSequence, targetPosition) in targets:
                #best_off_targets[source][targetPosition]={}
                #all_on_targets[source][targetPosition]={}
                dG_exchange = Cas9Calculator.calc_dG_exchange(g,targetSequence)
                dG_target = dG_PAM + dG_supercoiling + dG_exchange
                found_dg_target = False
                # if guide sequence matches, save dG target
                if targetSequence==g:
                    print "Target Binding Site Found with dG_target = %s kcal/mol" % dG_target
                    found_dg_target=dG_target
                    if source not in all_on_targets.keys():
                        all_on_targets[source] = {targetPosition:{'sequence': targetSequence,
                                                            'dG_PAM': dG_PAM,
                                                            'full_PAM': fullPAM,
                                                            'dG_exchange': dG_exchange,
                                                            'dG_supercoiling': dG_supercoiling,
                                                            'dG_target': dG_target}}
                    else:
                        all_on_targets[source][targetPosition] = {'sequence': targetSequence,
                                                            'dG_PAM': dG_PAM,
                                                            'full_PAM': fullPAM,
                                                            'dG_exchange': dG_exchange,
                                                            'dG_supercoiling': dG_supercoiling,
                                                            'dG_target': dG_target}
                else:
                    if dG_target<-3:
                        dG_offs.append((source,targetPosition,targetSequence,dG_PAM,fullPAM,dG_exchange,dG_supercoiling,dG_target))
                # adds each dG_target to guide's partition fxn
                partition_function+=math.exp(-dG_target/Cas9Calculator.RT)
    dG_offs=sorted(dG_offs, reverse=True, key = lambda x: x[7])
    
    for i in dG_offs[:1000]:
        if i[0] not in all_off_targets.keys():
            all_off_targets[i[0]] = {i[1]:{'sequence': i[2],
                                            'dG_PAM': i[3],
                                            'full_PAM': i[4],
                                            'dG_exchange': i[5],
                                            'dG_supercoiling': i[6],
                                            'dG_target': i[7]}}
        else:
            all_off_targets[i[0]][i[1]] = {'sequence': i[2],
                                            'dG_PAM': i[3],
                                            'full_PAM': i[4],
                                            'dG_exchange': i[5],
                                            'dG_supercoiling': i[6],
                                            'dG_target': i[7]}

    pickled = cPickle.dumps({'all_on_targets' : all_on_targets, 'all_off_targets' : all_off_targets})
    #print "Size of Pickled Dictionary containing all_on_targets and all_off_targets: %s bytes" % sys.getsizeof(pickled)
    
    return found_dg_target, partition_function , all_on_targets , all_off_targets


# set up GA toolbox
toolbox = base.Toolbox()
toolbox.register("mate", tools.cxTwoPoint) # crossover between models
toolbox.register("init_full", init_de_novo) # initializing parameters
toolbox.register("init_quick", init_de_novo) # initializing parameters
toolbox.register("select", customSelect) # selection
toolbox.register("evaluate", evaluate_construct) # evaluation for model quality
    
def run_GA_full(outputData,name,Cas9Calculator,array_length,target_regions,guide_distribution,parts,arms,strand,H=None,prefeed=None, background = None,quickmode = False,terms=True,cp=True,verbose=False):
    '''
    primary function of ELSA Calculator, calls part selection functions and runs GA
    '''

    try:
        
        #sets up multiprocessing via MPI
        from mpi4py import MPI
        from MPI_pool import Pool as MPIPool
        import copy_reg
        import cPickle as pickle

        # Using copy_reg to allow the pickle module to instance methods
        # Replicates multiprocessing module's ForkingPickler

        #copy_reg.pickle(types.MethodType, _reduce_method)
        #print "test"
        pool = MPIPool(MPI.COMM_WORLD)
        print "MPI Pool Started"

    except:
        pool = None
        print traceback.format_exc()
        print "Could not start MPI Pool. Using serial map instead."

    #start time for benchmarking
    t0=time.time()

    # sets up classes for GA
    ELSA= ELSA_Calculator(Cas9Calculator)

    creator.create("FitnessMulti", base.Fitness, weights=(-4.0,-1.0,-2.0,1.0,-1.0,-1.0))
    creator.create("Individual", list, fitness=creator.FitnessMulti)

    t1=time.time()
    print "VLSA init:", t1-t0

    ELSA.H=H
    

    module_length=len(parts["order"]) #

    #if multiprocessing, sets up mapping to multiple cores
    if pool is not None:
        pool.start()
        if pool.is_master():
            map_function = pool.map
            print 'INFO: Using MPI Pool'
    else:
        map_function = map
        print 'INFO: Using Serial Mode'

    if pool is None or pool.is_master():
        
        # checkpoint(cp) timesaving measure for reruns
        if cp:

            try:
                #loads guide information from checkpoint pickle
                print "trying to open checkpoint pickle file"
                guide_storage = pkl.load(open('%s_checkpoint.pkl'%(name), "r"))
                guides_list = guide_storage["guides_list"]
                guide_dict = guide_storage["guides_dict"]
                g_flag_temp = guide_storage["g_flag"]
                true_dist = guide_storage["true_dist"]

                print "Opened checkpoint pickle file"
                #print "prefeed: ", prefeed
                #print "g_flag_temp: ", g_flag_temp
                if prefeed and g_flag_temp:
                    g_flag = g_flag_temp
                    #print "Prefeed and g_flag both True"
                elif not prefeed and not g_flag_temp:
                    g_flag = g_flag_temp
                    #print "Prefeed and g_flag both False"
                    #print g_flag
                else:
                    #print "Prefeed and g_flag mismatched"
                    print g_flag
            except:
                # checkpointing fails
                print "Could not open checkpoint pickle file"

                # calculates guide information
                guides_list,guide_dict,g_flag,true_dist,error = ELSA.find_guides(target_regions,array_length,strand,guide_distribution,prefeed,pool)
                
                if error:
                    if pool is not None: pool.stop()
                    outputData['error']=guide_dict['warnings']
                    return outputData

                guide_storage = {"guides_list":guides_list,
                    "guides_dict":guide_dict,
                    "g_flag":g_flag,
                    "true_dist":true_dist}

                # saves checkpoint
                pickle.dump(guide_storage,open("%s_checkpoint.pkl"%(name),"w"))

        else:
            print "Checkpointing is not enabled"

            # calculates guide information
            guides_list,guide_dict,g_flag,true_dist,error = ELSA.find_guides(target_regions,array_length,strand,guide_distribution, prefeed,pool)
            
            if error:
                if pool is not None: pool.stop()
                outputData['error']=guide_dict['warnings']
                return outputData

            guide_storage = {"guides_list":guides_list,
                "guides_dict":guide_dict,
                 "g_flag":g_flag,
                 "true_dist":true_dist}

            #saves checkpoint
            pickle.dump(guide_storage,open("%s_checkpoint.pkl"%(name),"w"))
        
        # adds guide info to parts dict
        parts["guides"]=guides_list
        if true_dist:
            array_length =sum(true_dist)
        t3=time.time()
        print 'guide generation:',t3-t1

        # if quickmode, loads part lists from previous designs
        if quickmode:
            ELSA.quickmode=True
            if array_length<=8:
                parts['promoters'] = [promoters[p] for p in quick_small['promoters']]
                parts['handles'] = [handles_all[h] for h in quick_small['handles']]
                parts['spacers'] = [spacers[s] for s in quick_small['spacers']]
            elif array_length<=10:
                parts['promoters'] = [promoters[p] for p in quick_medium['promoters']]
                parts['handles'] = [handles_all[h] for h in quick_medium['handles']]
                parts['spacers'] = [spacers[s] for s in quick_medium['spacers']]
            elif array_length<=12:
                parts['promoters'] = [promoters[p] for p in quick_large['promoters']]
                parts['handles'] = [handles_all[h] for h in quick_large['handles']]
                parts['spacers'] = [spacers[s] for s in quick_large['spacers']]
            elif array_length<=20:
                parts['promoters'] = [promoters[p] for p in quick_max['promoters']]
                parts['handles'] = [handles_all[h] for h in quick_max['handles']]
                parts['spacers'] = [spacers[s] for s in quick_max['spacers']]
            else:
                if pool is not None: pool.stop()
                outputData['error']="ELSA is too large to design with quickmode"
                return outputData
            new_parts=parts

        #if not quickmode, use full part lists to choose parts with GA
        else:
            parts['promoters']=[str(record.seq) for record in SeqIO.parse(cwd + "/../parts/promoters_{}.fa".format(H),"fasta")]
            parts['spacers']=[str(record.seq) for record in SeqIO.parse(cwd + "/../parts/spacers_{}.fa".format(H),"fasta")]
            
            # if designing new handles, use assemble handles
            if terms:
                    parts["handles"]=ELSA.assemble_handles(H,[str(record.seq) for record in SeqIO.parse(cwd + "/../parts/backbones_{}.fa".format(H),"fasta")],[str(record.seq) for record in SeqIO.parse(cwd + "/../parts/terminators_{}.fa".format(H),"fasta")])

            else:
                parts["handles"]=[str(record.seq) for record in SeqIO.parse(cwd + "/../parts/handles_{}.fa".format(H),"fasta")]

            # establishes kmers for background prefilter
            if background:
                MF=KMP()
                kmers = set([background[i:i+20]for i in range(len(background)-20)])


            # reassemble parts dict using non repetitive parts according to H
            new_parts={'order':list(parts["order"])}

            for i in parts["order"]:
                t4=time.time()

                if i !='guides':
                    temp_parts = parts[i]

                    # removes parts found in background
                    if background:
                        for k in kmers:
                            remaining_parts_list = []
                            for t in temp_parts:
                                match = False
                                try:
                                    match = MF.match_iter(t,k).next() >-1
                                except:
                                    pass
                                if not match:
                                    remaining_parts_list.append(t)
                                else:
                                    print k,t
                            temp_parts = remaining_parts_list

                    if array_length>len(temp_parts):
                        if pool is not None: pool.stop()
                        outputData['error'] = "Insufficient parts for provided number of guides/provided H. Rerun with higher H or fewer guides"
                        return outputData
                    new_parts[i]=temp_parts
                
            new_parts["guides"]=guides_list

        # start GA
        t5=time.time()
        ELSA.new_parts=new_parts
        ELSA.guide_distribution=true_dist
        if true_dist:
            ELSA.array_length=sum(true_dist)
        else:
            ELSA.array_length=array_length
        ELSA.guide_dict=guide_dict
        ELSA.g_flag=g_flag

        outputData['guides'] = guides_list
        outputData['total_guide_dict'] = guide_dict
        outputData['actual_sgRNAs_per_target']=true_dist

        print "total part validation:",t5-t3

        # initialize output_dict
        outputData['GA']={}
        outputData['POP_SIZE'] = ELSA.POP_SIZE
        outputData['NGEN'] = ELSA.NGEN
        outputData['CXPB'] = ELSA.CXPB
        outputData['MUTPB'] = ELSA.MUTPB
        outputData['NUM_ELITES'] = ELSA.NUM_ELITES

        # creates tools to create design model parameters for optimizing
        toolbox.register("construct", toolbox.init_full, creator.Individual,
                        (ELSA))
        toolbox.register("population", tools.initRepeat, list, toolbox.construct)

        # creates population of ELSA designs
        population = toolbox.population(n=ELSA.POP_SIZE)

        # establishes hall of fame for fittest members of the run
        best = tools.ParetoFront()

        # runs for NGEN generations
        print "init time:", t5-t0


        for g in range(ELSA.NGEN):

            #prints population members and their fitnesses
            print "*" * 80
            print "P GENERATION #%s" % g

            for i in population:
                print i,  i.fitness.values

            print "P GENERATION #%s" % g
            print "*" * 80

            t6=time.time()

            print "current run time:",t6-t0

            # Clone the selected individuals
            offspring = map(toolbox.clone, population)

            # Apply crossover on the offspring
            for child1, child2 in zip(offspring[::2], offspring[1::2]):

                if random.random() < ELSA.CXPB:
                    toolbox.mate(child1, child2)
                    del child1.fitness.values
                    del child2.fitness.values    

            # Apply mutation 1 (scaling) on the offspring        
            for mutant in offspring:

                if random.random() < ELSA.MUTPB:
                    mutant = mutate_de_novo(mutant, ELSA)
                    del mutant.fitness.values

            # Evaluate the individuals with an invalid fitness
            invalid = [i for i in offspring if not i.fitness.valid]
            evaluation_inputs = [(i,ELSA.new_parts,ELSA.array_length,arms,guide_dict,ELSA.guide_distribution,g_flag,ELSA.RT) for i in invalid]
            fitnesses = map_function(toolbox.evaluate,  evaluation_inputs)

            for i, fit in zip(invalid, fitnesses):
                i.fitness.values = fit

            # The population is entirely replaced by the offspring
            # Select the next generation individuals
            best.update(population)
            successors = toolbox.select(offspring, population,ELSA)
            population[:] = successors

            t7=time.time()
            print "last round time:", t7-t6

        #sorts and returns population of final generation, prints best member of the run
        sorted_population = [x for x in sorted(best, key = lambda x: x.fitness.values)]

        outputData['sorted_population']={}

        for i in range(len(sorted_population)):
            outputData['sorted_population'][i]=(generate_product_dict(name,i,list(sorted_population[i]),ELSA.array_length,ELSA.new_parts,g_flag),list(sorted_population[i].fitness.values)) # print i, i.fitness.values 

        t8=time.time()
        print "total run time:", t8-t0
        
        if pool is not None: pool.stop()
        return outputData

    else:
        pass #workers see this                  
def generate_product_dict(name,index,design,array_length,new_parts,g_flag):
    
    design_list =[]
    sequence=SeqRecord(Seq('',IUPAC.unambiguous_dna))
    sequence.name="{}_{}".format(name,index)
    
    for i in range(array_length):

        for j in range(len(new_parts["order"])):

            if new_parts["order"][j]=="guides" and g_flag == False:
                dna=SeqRecord(new_parts[new_parts["order"][j]][design[j][i][0]][design[j][i][1]])
                f = SeqFeature(
                    FeatureLocation(0,len(dna.seq)),
                    type="guide")
                f.qualifiers['label'] = ["guide {}".format(str(i))]
                dna.features.append(f)
                part = new_parts['order'][j]
                index = (design[j][i][0],design[j][i][1])
                design_list.append((part,index,str(dna.seq)))
                sequence+=dna

            elif new_parts["order"][j]=="guides" and g_flag == True:
                dna=SeqRecord(new_parts[new_parts["order"][j]][design[j][i]])
                f = SeqFeature(
                    FeatureLocation(0,len(dna.seq)),
                    type="guide")
                f.qualifiers['label'] = ["guide {}".format(str(i))]
                dna.features.append(f)
                part = new_parts['order'][j]
                index = design[j][i]
                design_list.append((part,index,str(dna.seq)))
                sequence+=dna
            elif new_parts["order"][j]=="handles":
                dna=SeqRecord(new_parts[new_parts["order"][j]][design[j][i]])
                backbone = dna[:61]
                terminator = dna[61:]
                f1= SeqFeature(
                    FeatureLocation(0,61),
                    type="backbone")
                f1.qualifiers['label'] = [part_names['backbones'][backbones.index(str(backbone.seq))]]
                f2= SeqFeature(
                    FeatureLocation(61,len(dna.seq)),
                    type="terminator")
                f2.qualifiers['label'] = [part_names['terminators'][terminators.index(str(terminator.seq))]]
                dna.features.append(f1)
                dna.features.append(f2)
                part = new_parts['order'][j]
                index1 = backbones.index(str(backbone.seq))
                index2 = terminators.index(str(terminator.seq))
                design_list.append(('backbone',index1,str(backbone.seq)))
                design_list.append(('terminator',index2,str(terminator.seq)))
                sequence+=dna
            else:
                dna=SeqRecord(new_parts[new_parts["order"][j]][design[j][i]])
                f = SeqFeature(
                    FeatureLocation(0,len(dna.seq)),
                    type=new_parts["order"][j])
                f.qualifiers['label'] = [part_names[new_parts["order"][j]][part_seqs[new_parts['order'][j]].index(new_parts[new_parts["order"][j]][design[j][i]])]]
                dna.features.append(f)
                part = new_parts['order'][j]
                index = part_seqs[new_parts['order'][j]].index(new_parts[new_parts["order"][j]][design[j][i]])
                design_list.append((part,index,str(dna.seq)))
                sequence+=dna

    assess=seq_assess()
    assess.load(sequence.seq)
    assess.run()
    s_dict = assess.s_info
    for (key, info) in s_dict.items():
                    
        location_list = info['locations']
        type_list = info['types']

        for (location, type) in zip(location_list, type_list):
            synthesis_complexity_color = '#CC3300'
            loc = FeatureLocation(location[0], location[1])
            f = SeqFeature(location = loc, type = 'bad_synth', qualifiers={"label" : type, "ApEinfo_label" : type,
                    "ApEinfo_fwdcolor" : synthesis_complexity_color, "ApEinfo_revcolor" : synthesis_complexity_color} )
            sequence.features.append(f)


    genbank = StringIO.StringIO() 

    handle_temp = StringIO.StringIO() 

    SeqIO.write(sequence, handle_temp, "gb")

    fixGenBankforAPE(handle_temp,genbank)

    handle_temp.close()


    result = {"sequence":str(sequence.seq),"composition": design_list,"s_dict": s_dict, "r_dist": assess.repeat_distribution, 'RepeatDict':assess.RepeatDict, 'genbank':genbank.getvalue()}
    genbank.close()
    return result

def run_ELSA(inputs):

    # passes input dicts to variables
    name = inputs["name"][:16]
    host_genomes=inputs["host_genomes"]
    mode = inputs["mode"]
    outputData = {}

    Cas9Calculator = clCas9Calculator(host_genomes,quickmode=True)
    if mode == "design":
        quick_mode = inputs["quickmode"]
        guide_mode = inputs["guide_mode"]
        strand = inputs["strand"]
        verbose = inputs['verbose']
        arms = inputs['arms']
        background = inputs['background']
        # initializes Cas9 Calculator
        
        #ELSA order default
        parts={"order":["promoters","guides","handles","spacers"]}

        # selects inputs to propogate based on flags
        if guide_mode == 'manual':
            # guides pre picked, regions not needed
            prefeed = inputs["guide_sequence_list"]
            target_regions = None
            guide_distribution = None
            number_of_guides = len(prefeed)
        elif guide_mode == 'auto':
            # guides chosen from targeted regions
            target_regions = inputs['target_sequence_list']
            guide_distribution = inputs["sgRNAs_per_target_list"]
            prefeed = None
            number_of_guides = sum(guide_distribution)
        else:
            raise Exception, "The input variable guide_mode must be either 'manual' or 'auto'."
        
        # sets quickmode rules
        if quick_mode:
            L=20
            outputData = run_GA_full(outputData,name,Cas9Calculator,number_of_guides,target_regions,guide_distribution,parts,arms,strand,L,prefeed, background, quick_mode,terms=False, cp=True,verbose = verbose)
        else:
            L= inputs["L"]
            outputData = run_GA_full(outputData,name,Cas9Calculator,number_of_guides,target_regions,guide_distribution,parts,arms,strand,L,prefeed, background, quick_mode,terms=True, cp=True,verbose = verbose)
    else:
        sequence = inputs["sequence"]
        outputData = ELSA_evaluate_mode(outputData,name,sequence, host_genomes, Cas9Calculator,cp = True)
    return outputData
    
def ELSA_evaluate_mode(outputData,name,sequence, host_genomes, Cas9Calculator,cp=False):
    #start MPI
    try:
        
        #sets up multiprocessing via MPI
        from mpi4py import MPI
        from MPI_pool import Pool as MPIPool
        import copy_reg
        import cPickle as pickle

        # Using copy_reg to allow the pickle module to instance methods
        # Replicates multiprocessing module's ForkingPickler

        copy_reg.pickle(types.MethodType, _reduce_method)
        pool = MPIPool(MPI.COMM_WORLD)
        print "MPI Pool Started"

    except:
        pool = None
        print traceback.format_exc()
        print "Could not start MPI Pool. Using serial map instead."

    if pool is not None:
        pool.start()
        if pool.is_master():
            map_function = pool.map
            #print 'yes'
    else:
        map_function = map
        #print 'no'

    #initialize itermatch for guide finding
    MF=KMP()

    if pool is None or pool.is_master():
        outputData['name'] = name
        if cp:
            try:
                print "trying to load checkpoint"
                checkpoint=pickle.load(open("%s_checkpoint.pkl"%(name),"r"))

                guides = checkpoint['guides']
                g_pos = checkpoint['g_pos']
                ELSA_handles = checkpoint['ELSA_handles']
                h_pos = checkpoint['h_pos']
                total_guide_dict=  checkpoint['tgt']

                pf = [total_guide_dict[g]["partition_fxn"] for g in guides]
                off_targets=[total_guide_dict[g]["off_targets"] for g in guides]
                outputData['guides'] = guides
                outputData['total_guide_dict'] = total_guide_dict

            except:
                print "unableto load checkpoint"
                g_input_set = [(sequence,b,MF,False, True) for b in backbones]
                g_l_input_set = [(sequence,b,MF,True, True) for b in backbones]
                h_input_set = [(sequence,b,MF,False, False) for b in backbones]
                h_l_input_set = [(sequence,b,MF,True, False) for b in backbones]
                g_output = map_function(part_scan,g_input_set)
                g_loc = map_function(part_scan,g_l_input_set)
                h_output = map_function(part_scan,h_input_set)
                h_loc = map_function(part_scan,h_l_input_set)
                guides = []
                g_pos = []
                ELSA_handles = []
                h_pos = [] 
                for i in g_output:
                    guides.extend(i)
                for i in g_loc:
                    g_pos.extend(i)
                for i in h_output:
                    ELSA_handles.extend(i)
                for i in h_loc:
                    h_pos.extend(i)
                # set up Cas9Calculator
                if len(guides)==0:
                    outputData['error'] = "No guides/handles detected in provided sequence"
                    if pool is not None: pool.stop()
                    return outputData
                total_guide_dict = {}

                # processing by Cas9Calculator
                inputs = [(g, Cas9Calculator) for g in guides]
                outputs_2 = map_function(partition_function, inputs)

                # populating information dictionaries
                for i in range(len(guides)):
                    total_guide_dict[guides[i]] = {}
                    dG_target, partition_fxn, on, off = outputs_2[i]
                    total_guide_dict[guides[i]]['dG_target'] = dG_target
                    total_guide_dict[guides[i]]['partition_fxn'] = partition_fxn
                    total_guide_dict[guides[i]]["off_targets"] = off
                    if 'position' not in total_guide_dict[guides[i]].keys():
                        total_guide_dict[guides[i]]['position'] = [g_pos[i]]
                    else:
                        total_guide_dict[guides[i]]['position'].append(g_pos[i])

                guide_list = guides

                outputData['guides'] = guides
                outputData['total_guide_dict'] = total_guide_dict

                pf = [total_guide_dict[g]["partition_fxn"] for g in guides]
                off_targets=[total_guide_dict[g]["off_targets"] for g in guides]
                checkpoint={'guides':guides,
                            'g_pos':g_pos,
                            'ELSA_handles': ELSA_handles,
                            'h_pos':h_pos,
                            'tgt':total_guide_dict}

                pickle.dump(checkpoint,open("%s_checkpoint.pkl"%(name),"w"))

        else:
            print "no checkpoint"
            g_input_set = [(sequence,b,MF,False, True) for b in backbones]
            g_l_input_set = [(sequence,b,MF,True, True) for b in backbones]
            h_input_set = [(sequence,b,MF,False, False) for b in backbones]
            h_l_input_set = [(sequence,b,MF,True, False) for b in backbones]
            g_output = map_function(part_scan,g_input_set)
            g_loc = map_function(part_scan,g_l_input_set)
            h_output = map_function(part_scan,h_input_set)
            h_loc = map_function(part_scan,h_l_input_set)
            guides = []
            g_pos = []
            ELSA_handles = []
            h_pos = [] 
            for i in g_output:
                guides.extend(i)
            for i in g_loc:
                g_pos.extend(i)
            for i in h_output:
                ELSA_handles.extend(i)
            for i in h_loc:
                h_pos.extend(i)
            # set up Cas9Calculator
            if len(guides)==0:
                outputData['error'] = "No guides/handles detected in provided sequence"
                if pool is not None: pool.stop()
                return outputData
            total_guide_dict = {}

            # processing by Cas9Calculator
            inputs = [(g, Cas9Calculator) for g in guides]
            outputs_2 = map_function(partition_function, inputs)

            # populating information dictionaries
            for i in range(len(guides)):
                total_guide_dict[guides[i]] = {}
                dG_target, partition_fxn, on, off = outputs_2[i]
                total_guide_dict[guides[i]]['dG_target'] = dG_target
                total_guide_dict[guides[i]]['partition_fxn'] = partition_fxn
                total_guide_dict[guides[i]]["off_targets"] = off
                if 'position' not in total_guide_dict[guides[i]].keys():
                    total_guide_dict[guides[i]]['position'] = [g_pos[i]]
                else:
                    total_guide_dict[guides[i]]['position'].append(g_pos[i])

            guide_list = guides

            outputData['guides'] = guides
            outputData['total_guide_dict'] = total_guide_dict

            pf = [total_guide_dict[g]["partition_fxn"] for g in guides]
            off_targets=[total_guide_dict[g]["off_targets"] for g in guides]
        # starts sequence assessment helpers
        assess=seq_assess()
        assess.load(sequence)
        assess.run_ELSA_analysis()

        s_dict = assess.s_info

        outputData['s_dict'] = s_dict
        #outputData['flagged'] = assess.bad_seqs

        # repeat score from seq assess, scaled
        repeats=int(assess.score[0])

        # haripin score from seq assess
        hairpins=int(assess.score[1])

        # nucleotide score from seq assess
        GC=int(assess.score[2])

        # converts partition fxn to linear term 
        PartFxn=math.log(sum(pf))
        
        PF_off = 0
        for o in range(len(off_targets)):
            for source,position in off_targets[o].items():
                for targets, data in position.items(): 
                    PF_off += math.exp(-data['dG_target']/Cas9Calculator.RT) 

        PartFxn_off=math.log(PF_off)

        dx_t=0

        outputData['score'] = (repeats,hairpins,GC,dx_t,PartFxn,PartFxn_off)
        outputData['score_dict']={   "repeats" : repeats,
                                     "hairpins" : hairpins,
                                     "GC" : GC,
                                     "p_fxn" : PartFxn,
                                     "p_off" : PartFxn_off,
                                     "dx_t" : dx_t}
        print "AT end"
        print name, len(name)
        genbank = StringIO.StringIO()

        sequence,genbank=export_sequence2(sequence, guides, total_guide_dict, ELSA_handles, h_pos, s_dict, name, genbank)  

        outputData['sequence']=str(sequence.seq)
        outputData["genbank"]=genbank.getvalue()

        genbank.close()
        
        if pool is not None: pool.stop()

        return outputData
    else:
        pass #workers see this


def export_sequence2(raw, guides,  guide_dict, handles_list, handles_loc, s_dict,name,handle):
    '''
    Converts sequence into a .gb file with anotations
    '''

    
    sequence=SeqRecord(Seq(raw,IUPAC.unambiguous_dna))
    sequence.name=name
    
    for g in range(len(guides)):
        #print g
        for p in guide_dict[guides[g]]['position']:
            #print p
            #print len(str(sequence))
            f = SeqFeature(
                    FeatureLocation(p[0],p[1]),
                    type="guide")
            f.qualifiers['label'] = ["guide {}".format(str(g))]
            sequence.features.append(f)

    for h in range(len(handles_list)):
        f = SeqFeature(
                FeatureLocation(handles_loc[h][0],handles_loc[h][0]),
                type="backbone")
        f.qualifiers['label'] = [part_names["backbones"][backbones.index(handles_list[h])]]
        sequence.features.append(f)

    for (key, info) in s_dict.items():
                    
        location_list = info['locations']
        type_list = info['types']

        for (location, type) in zip(location_list, type_list):
            synthesis_complexity_color = '#CC3300'
            loc = FeatureLocation(location[0], location[1])
                
            f = SeqFeature(location = loc, type = 'bad_synth', qualifiers={"label" : type, "ApEinfo_label" : type,
                    "ApEinfo_fwdcolor" : synthesis_complexity_color, "ApEinfo_revcolor" : synthesis_complexity_color} )
        sequence.features.append(f)

    handle_temp = StringIO.StringIO() 

    SeqIO.write(sequence, handle_temp, "gb")

    fixGenBankforAPE(handle_temp,handle)

    handle_temp.close()
    return sequence, handle

def part_scan((sequence,part,KMP,location,guide)):
    if guide:
        if location:
            return [(result-20,result) for result in KMP.match_iter(sequence,part)]
        else:
            return [sequence[result-20:result] for result in KMP.match_iter(sequence,part)]
    else:        
        if location:
            return [(result,result+len(part)) for result in KMP.match_iter(sequence,part)]
        else: 
            return [sequence[result:result+len(part)] for result in KMP.match_iter(sequence,part)]
     

def fixGenBankforAPE(handle_input, handle_output):

    handle_input.seek(0)
    
    for line in handle_input.readlines():
        words = line.split('=')
        if (words[0].strip() == "/ApEinfo_fwdcolor") or (words[0].strip() == "/ApEinfo_revcolor"):
            handle_output.write(words[0] + "=" + words[1].replace('"',''))
        else:
            handle_output.write(line)

if __name__ == "__main__":
    #
    inputs={}

    # mode: design or evaluate
    inputs["mode"] = "design"

    # name = name for identifying run
    inputs['name']="test_tr_quick"

    # host_genomes = DNA sequences that sgRNAs might bind to
    #inputs["host_genomes"]=['/home/ubuntu/Genomes/refseq/bacteria/GCF_000005845.2/GCF_000005845.2_ASM584v2_genomic.gbff.gz']

    inputs["host_genomes"]=[record for record in SeqIO.parse(cwd + "/NC_000913.gb","gb")]

    # verbose = printing to console
    inputs["verbose"] = True

    # design mode
    if inputs["mode"] == 'design':
        # L = longest permitted repeat in the design 
        inputs["L"]=17

        # arms = upstream or downstream required motifs
        inputs["arms"]=['','']#["ctcacttgtagaacggtgatcagcctgtgctctagagcctgatagttgagcgatacacacCTCGAGTCGTCATGAGGGCTTTATCTCATATTGTTCAAATCACCAGCAAACACCGACATATTTGCAACTCAATATTCACAACAACCAGAAAAGCCAACCTGCGGGTTGGCTTTTTTATGCAACTCG","AAcACGAAtAAAGGATCCCCTTTCGGGAGGCCTCTTTTCTTTACACTGCGCCACTATTTTCGCTATGGTTATGCGTAAGCATTGCTGTTGCTTCGTCGCGGCAATATAATGAGAATTATTACTAGTtgatcgttgaagtcgacctacatcgagtgcgcactatcaagagtgttccagtcacgcgat"]
        
        
        
        # quickmode = whether or not you are reusing an established ELSA design
        inputs["quickmode"]=False
        
        # cut_mode = "coding", "template", None: determines which strand is used for determining usable guides  
        inputs["strand"]="template"
        
        
        # target_regions = regions to be targeted by Cas9
        inputs["target_sequence_list"] = ["TTTTCCATCACAGCGTTAAGCGTCTTAATCTTAAACATGTATATTAGGGGGGGGGGGG",
                                          "GGACCATAACTTTAAAAGATGTGATCatgTTTTGGTGTGGCAATGTGAACCCACTATCTGTACTCTTGCAGCAAACAGCTGCTCA"]
        # guide_distribution = number of guides targeting each locus
        inputs['sgRNAs_per_target_list']=[4,4]
        
        
        # guides = list of preselected guides for use directly, instead of picking from target regions 
        
        #inputs['guide_sequence_list'] = ['TAAGGTTAAGACGCTTAACG','AGCCGTTTCCTGCCGGAGTA', 'AACAAGTTTAGTTCATCTGA']

        
        inputs['background'] = 'TTGACATGACTCTCCAGCTGTGCTATAATTGTACT'


        # if guides are provided, only those guides will be used
        if "guide_sequence_list" in inputs.keys():
            inputs["guide_mode"]="manual"
        else:
            inputs["guide_mode"]="auto"

    # evaluate mode
    else:

        #sequence: string corresponding to ELSA to be evaluated
        inputs["sequence"] = 'AGTTACACTTACCCTACTTTATCGGATTCTGAGGAACAGGAGACTGATTATTGACATTAGCACTTGAGCTGATTATAATGGGCCGCTCTTTCGTTACCGCCGATTGTTCTAGAGCTGGTAACAGCAAGTTAGAATAAGTCTAGTCCATTATCAACTGGAAACAGTGGCCAGAAAGGGTCCTGAATTTCAGGGCCCTTTTTTTACATTTCACGCCGATTGTCCTTCAGGTTTTACAGAATAAGATAACTACGGATAGTTGACATGCGTGATTTAACATTCTATAATTGCACATAAGGTTAAGACGCTTAACGGATGTAGATGTAGAAATACAAGGTTACATTAAGGCCCGTCCGTAATCAACTTGAAGAAGTGTTCCATCGGGTCCGAATTTTCGGACCTTTTCTCCGCATACTGATAATAGTTGATTGTCTGAAGTGTAAACCCTCCACCGAAAGGCTTTTGACAGACCTTATCTACATGGTTATAATCTGAATACCATTTACTGCATCGATGAGTTGTAGATCTAGAAATAGAATGTTACAATTAGGCTAGTCCGTTATGAACATGAAAATGTGAGAAAAGAGGCCGCGAAAGCGGCCTTTTTTCGTTTATACTGAGAGCGTTGACCAATAAGGATTACACTATCTTCTGTTGTGACACTTGACACGGATCTTCGCTGAACGTATAATGAGAAAGAAACCCTGCTGTTGCATCGGTTTTAGAGTGAGAAATCACAAGTTAAAATAAGGCTAGACCGTTATCAACTAGAAATAGTGTTATTGAACACCCAAATCGGGTGTTTTTTTGTTTGATTACACTGTGATACTTCTGACGCAGGTTTTCCAACGAGTTACGAATAATTGACAGGTGAACGCTCAGCTCTTATAATGCCTATAGCCGTTTCCTGCCGGAGTACTTTTAGAGATAGAAATATCAAGTTAAAAGAAGGCTAGTCCGTTACCAACTTGAAAAAGTGAGAAAAAAGGCACGTCATCTGACGTGCCTTTTTTATTTAAGTGAGGTGACTACTTCTCTGAAATCGTTTACAACCCTTCGGAATAAGATTTGACAACTGCTCAGCGAAATACTATAATGACTACTTTATCCTGAACAGTGATCCGCTTTAGAGCTAGAAATAGCAGGTTAAAGTAAGGCCAGTCCGTAATAAACTGGAAACAGTGAGAAAAGAGACGCTTTCGAGCGTCTTTTTTCGTTTGAACCAGTCTCGTAGTTGTTACAGCGATAAGAATAGGTGTTGAAATACTCTTGACAGCATCTGCTTTGTCACCTATAATTCAATGAGCTGCAACCGTTTGTTTCAGTTTTGGACCTAGAAATAGGAAGTCAAAATAAGGCTGGACCGACATGTAATCGAAAGATTTAGTCAAAAGCCTCCGGTCGGAGGCTTTTGACTTTCTGTAAGACAGAGTAGGGTATTATCACTATTCGCTGGAACTTCACCAATCTTTGACACCTCATCTTATAGTTCCTATAATTTCTATCGGAGACCTGGCGGCAGTATGTTGGAGAGCAAGACATTGCAAGTTCCAATAAGGCGTGTCCGATAAAAGCTTGAGAAAGCAAAGTAATACAAAACAGGCCCAGGCGGCCTGTTTTGTCTTTTTAATGTAATAATCCCAGACTCAGAATAGGAATCTTACCTGTCGTGTTGCGAAGTTTTGACATAAGTCGTATTCAAAGATATAATATAGGTAGCCAAAATCCATCATCATGATCTGAGAGCCAAAAATGGCAAGTTCAGATAAGGCCAGACCGTTACCAGCTTAAATAAGCGATCCTAAAGCCCCGAATTTTTTATAAATTCGGGGCTTTTTTACTAGTAGTTTATTCGCTCTATTGAGGTAGTCGTCAGAACCCTTATCAGGAAAAGTTGACATCCTCTCGTAGGACTCATATAATACTCAGAGTCTCTTTTTTCTGTATCGTTTTTAGAGGAAGGAATTCCAAGTTAAAAAAAGGCAGGACCGGGAACATGTTGAAAAACAGGCAAAAAAGCGCCTTTAGGGCGCTTTTTTACATTAGATAGTGCGAAGTCTTATTTGGATCCCTACACGAGAGTGATTTCAACCGTATTGACACCGGGTTGAATACTATCTATAATGTACGGTACGTGGCTAAAAAAACGTCGGAATAGAAAACAAAAGTTTAAGTTATTCTAAGGCCAGTCCGGAATCATCCTAAAAAGGAGTTATTGAACACCCGAAAGGGTGTTTTTTTGTTTAGGCAGTTATCTCTTACCGAGTTTTACTTCAGTGTGCGAATAGACAACAATTGACATTTCGTCAAGAGTCGACTATAATATCGCGGTCTGTAGGTCCAGATTAACGATTTCGAGCTAGGCATAGCAAGTGAAATTAAGGCTGGTCCATTAACACCTTGAAAAAGGGAACAATAAGGCCTCCCTTTAGGGGGGGCCTTTTTTATTGAAAAGTTCACCGTTATTATCCTGTAGGTAGTATTTTCAGCCACCAGAGTAGTTGACAGTCCTCGAACACCTCTATATAATAGTGTCGAAGAAGCTGCCCATAATCGTTTTCAGATTTGGAAACAAAACGTTGAAAAAAGGCAAGTCCGTTATGAACGCGAAAGCGTGCGAAAAAACCCGCTTCGGCGGGTTTTTTTATAGTTTCCAGAGGACCTTCACGGATAAAATAGATTACAGTTCTCGTCGTAGTATTGACATGAGCTCGTCGTCAGGATATATAGCTTTTGCTAAAGGATACCTGATAGGTTGTAGAGCTAGCAATAGCAGGTTACAATAAGGCTCGTCCGTTATAAACATGAAAATGTGACTAAAAAGGCCGCTCTGCGGCCTTTTTTCTTTTTTGATTGAAGGACCGTAGCAGTCACAGAGTGTAACTTTATTCCCAGTAATTTGACACATTAGGATGGACGTATTATAATATGCCCAGACAGTTCACGAACCATTGGTTTTAGATCACGAAAGTGAAAGTTAAAATAAGCCTAGCCCGTTACCAACTGGAAACAGTGACTTAAGACCGCCGGTCTTGTCCACTACCTTGCAGTAATGCGGTGGACAGGATCGGCGGTTTTCTTTTCTTCCGTAGATAATAGAATAAGGTGCCCTCAGATTGTTGGAAGCGACTTTTATTGACACTATGGTCCGCAAGCATTATAATGCTCTGAACAAGTTTAGTTCATCTGAGTATCTGAACTCGACAGAGTAAGTAGATATAAGGCCAGTCCGTTAGCAACTTGAAAAAGTCGGAGATAAAACCGACCACGGCACCAGGCAGTGACCATGTGGTTTCTTCATCCTCTCACTTCTCTGTTGGCACGAAAAGGGCAATAAGATTTACGGATTACTATCTTGACATGAAGTGTTAGACGTCATATAATCGTGGTGCGGGTGGACATCCACTCGAGCTTCAGATCCAGAAATGGAAAGTTGAAGTGAGGCAGGTCCGGTAGCAACTCGAAAGAGTGAGAAAAGAGGGGAGCGGGAAACCGCTCCCCTTTTTTCGTTTCTGAATAAGCACTGTTGATAATCGCAATCTGTCTCTTCGTGAAAAGTAGCTTGACAATCGCTGTCTACGTGAATATAATGAATTTCAATGGCAGTGTGGCACTCACATTTTGGCGTCGAAAGACGAAGTAAAATGAAGGCGAGACCGATATCAACTGGAAGCAGTGTCTGGTAGTCCTGGTAAGACGCGAACAGCGTCGCATCAGGCATATTGCCAACTAGAGACTACTATTGTCGTTTATTATCGCAACAGAGGGAAGTTCACTGACCTATTGACATGACTCTCCAGCTGTGCTATAATTGTACTATTTGCCGAAACGTGCAGCCGTTTAAGAGCTAGAAATAGCACGTTTAAATAAGGCTAGTCCGTTTTCAACTTGAAAAAGTGATAACAAAGCCGGGTAATTCCCGGCTTTGTTGTATCGTGAACGACACTACTATTTCTTACGAGATACTTATTCTGGAAGCAACGGTTTGACATAGGCAAGCCAGTATAGTATAATCACATACACCGCACCAGCGACTGGACCTTACCGAACTAGGAATAGTAAGTGGTAAGAAGGCCTGACCGTAATAAGCCTGAAAAGGCGACCAAAAAGGGGGGATTTTATCTCCCCTTTAATTTTTCAGACTCCCTGTATCGTTGAAAAGTTGGAACACTGTGAATCCTATTACTGA'

    # runs ELSA Calculator
    outputData = run_ELSA(inputs)
    if outputData:
        if 'error' not in outputData.keys():
            pickle.dump(outputData,open("%s_results.pkl"%(inputs['name']),"w"))
            json.dumps(outputData)
        else:
            print outputData['error']
    
    #design mode outputs
    #outputData['name']: name of run
    #outputData['GA']: population at each generation
    #outputData['GA'][g]: population at a given generation g
    #outputData['POP_SIZE']: population size
    #outputData['NGEN']: number of generations
    #outputData['CXPB']: crossover probability
    #outputData['MUTPB']: mutation probability 
    #outputData['NUM_ELITES']: number of elites for selection
    #outputData['s_dict']: dictionary of synthesis complexity sites
    #outputData['flagged']: list of descriptions of detected synthesis issues
    #outputData['sorted_population']: final generation, sorted by fitness
    #outputData['sorted_population'][i]: each ith member of the final generation
    #outputData['design']: best individual design found during the GA run
    #outputData['scores']: scores of best member for each evaluate 
    #outputData['best_sequence']: sequence producted by the best design
    #outputData["genbank"]: stringified annotated genbank file of best design
    #outputData['error']: empty until input validation fails

    #evaluate results
    #outputData['name']: name of run
    #outputData['guides']: detected guides in ELSA
    #outputData['total_guide_dict']: guide information
    #outputData['s_dict']: dictionary of synthesis complexity sites
    #outputData['flagged']: list of descriptions of detected synthesis issues
    #outputData['score']: tuple of scores according to evaluation
    #outputData['score_dict']: dictionary of scores according to evaluation
    #outputData['sequence']: sequence of ELSA to be evaluated
    #outputData["genbank"]: stringified annotated genbank file with best design
    #outputData['error']: empty until input validation fails