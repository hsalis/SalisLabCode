
import sys
sys.path.append('../')
sys.path.append('models')
# sys.path.append('../datasets')
# sys.path.append('/usr/local/lib/python2.7/site-packages/')

import synbiomts
import cPickle as pickle
# import TranslationRateModels as tl

# We're going to use shelve to store model predictions
# as a dictionary-like persistance object of pandas dataframes
import shelve

# Import private Salis Lab code (latest versions of RBS Calculator)
sys.path.append('/home/alex/Private-Code')
from DNAc import *

from PyVRNA import PyVRNA
RNAEnergyModel = PyVRNA(dangles=0)

handle = open('models/rRNA_16S_3p_ends.p','r')
rRNA_16S_3p_ends = pickle.load(handle)
handle.close()

# Required function in RBS Calculators and UTR Designer
def get_rRNA(organism):
    # temporary until I update database:
    if   organism=="Corynebacterium glutamicum B-2784":
        organism = 'Corynebacterium glutamicum R'

    elif organism=="Pseudomonas fluorescens A506":

        return 'ACCTCCTTT'
    elif organism=="Escherichia coli BL21(DE3)":

        return 'ACCTCCTTA'
    else:
        pass
    return rRNA_16S_3p_ends[organism]

def find_best_start(mRNA, start_pos, predictions):
    '''Finds the number of start codons, the most highly translated start codon and the start range associated with it.
        Disallows start codons that have stop codons following them.'''
    
    stop_codons = ["TAA", "TAG", "TGA", "UAA", "UAG", "UGA"]
    num_start_codons = 0
    dG_tot_list = []
    dG_best = None
    best_start_index = None
    for index in range(len(predictions.RBS_list)):
    
        RBS = predictions.RBS_list[index]
        in_frame = (start_pos - RBS.start_position) % 3
        stop_codon_present = any([ (mRNA[pos:pos+3] in stop_codons) for pos in range(RBS.start_position+3, start_pos, 3)]) # stop codon present?

        if (in_frame == 0) and (not stop_codon_present): # if start is in frame and no downstream stop codon
            num_start_codons += 1
            dG_tot_list.append(RBS.dG_total)

        if (dG_best == None or RBS.dG_total < dG_best) and (in_frame == 0) and (not stop_codon_present):
            best_start_index = index
            dG_best = RBS.dG_total
            best_start_pos = RBS.start_position
    
    if not best_start_index is None:
        RBS = predictions.RBS_list[best_start_index]
    else:
        RBS = None
        best_start_pos = None

    return (RBS,best_start_pos)

import RBS_Calculator_v2_1
def RBSCalc_v2_1(sequence,rRNA,temp,startpos):
    start_range = [0,startpos+1]
    # rRNA = get_rRNA(organism)
    model = RBS_Calculator_v2_1.RBS_Calculator(sequence,start_range,rRNA)
    model.temp = temp
    model.run()
    output = model.output()
    (RBS,best_start_pos) = find_best_start(sequence,startpos,output)

    if RBS is None:
        results = {}
    else:    

        results = {
        'TIR' : RBS.tir,
        'dG_total' : RBS.dG_total,
        'dG_mRNA_rRNA' : RBS.dG_mRNA_rRNA,
        'dG_mRNA' : RBS.dG_mRNA,
        'dG_start' : RBS.dG_start,
        'dG_standby' : RBS.dG_standby,
        'dG_spacing' : RBS.dG_spacing,
        'used_mRNA_sequence' : RBS.sequence,
        'dG_UTR_fold' : RBS.dG_UTR_folding,
        'dG_SD_16S_hybrid' : RBS.dG_SD_hybridization,
        'dG_after_footprint' : RBS.dG_after_footprint_folding,
        'spacing_length' : RBS.spacing_length,
        'final_post_cutoff' : RBS.adaptive_post_cutoff,
        'binding_fragment' : RBS.optimum_hairpin_fragment,
        'surface_area' : RBS.optimum_available_RNA_surface_area,
        'dG_distortion' : RBS.optimum_distortion_energy,
        'dG_sliding' : RBS.optimum_sliding_energy,
        'dG_unfolding' : RBS.optimum_unfolding_energy,  
        'most_5p_paired_SD_mRNA' : RBS.most_5p_paired_SD_mRNA,
        'most_3p_paired_SD_mRNA' : RBS.most_3p_paired_SD_mRNA,
        'aligned_most_5p_paired_SD_mRNA' : RBS.aligned_most_5p_SD_border_mRNA,
        'aligned_most_3p_paired_SD_mRNA' : RBS.aligned_most_3p_SD_border_mRNA,
        'predicted_5pUTR' : sequence[:best_start_pos],
        'predicted_CDS' : sequence[best_start_pos:],
        }

        # Save dot-parentheses structure for initial state mRNA structure
        bpx = RBS.initial_structure['bp_x']
        bpy = RBS.initial_structure['bp_y']
        viennafld = RNAEnergyModel.bp2vienna(length=[len(RBS.sequence)],bpx=bpx, bpy=bpy)
        results['initial_structure'] = viennafld

    return results

def StructureOnlyModel(sequence,startpos):

    standby_seq = sequence[0:startpos-14]
    post_footprint_seq = sequence[startpos+13:]

    initial_state = RNAEnergyModel.RNAfold(sequence[:startpos+35])
    dG_mRNA = initial_state.energy
    initial_structure = initial_state.structure
    dG_UTR_fold = RNAEnergyModel.RNAfold(standby_seq).energy
    dG_after_footprint = RNAEnergyModel.RNAfold(post_footprint_seq).energy
    dG_total = dG_UTR_fold + dG_after_footprint - dG_mRNA

    results = {
    'dG_total': dG_total,
    'dG_UTR_fold': dG_UTR_fold,
    'dG_after_footprint': dG_after_footprint,
    'initial_structure': initial_structure
    }

    return results


def test_Ecoli_hypotheses():

    # last 13 nt of Escherichia coli 16S rRNA sequence
    rRNA_16S_end = "GAUCACCUCCUUA"

    # Generate all possible anti-SD sequences of length 6-13 nt
    candidates = [(rRNA_16S_end[i:i+j],i,j) for i in xrange(10) for j in xrange(4,14) if i+j<(len(rRNA_16S_end)+1)]
    candidates = sorted(candidates, key=lambda x: (x[2], x[1]))

    # define database filters
    # we are only going to test on E. coli for now
    filters = {
                "DATASET": ['EspahBorujeni_NAR_2013',
                            'EspahBorujeni_NAR_2015',
                            'EspahBorujeni_JACS_2016',
                            'EspahBorujeni_Footprint',
                            'Salis_Nat_Biotech_2009',
                            'Farasat_MSB_2014',
                            'Tian_NAR_2015',
                            'Bonde_NatMethods_IC_2016',
                            'Egbert_PNAS_2012',
                            'Hecht_NAR_2017',
                            'Beck_PLoS_2016'],

                "ORGANISM": ['Escherichia coli str. K-12 substr. MG1655',
                             'Escherichia coli str. K-12 substr. DH10B',
                             'Escherichia coli BL21(DE3)']
                }

    # Provide the pickled database file name
    dbfilename = '../geneticsystems.db'

    # add models to interface.Models
    # model names are the tested 16S rRNA sequence
    models = synbiomts.interface.Container()
    modelNames = []
    for (rRNA,pos,length) in candidates:
        models.add(model=RBSCalc_v2_1,rRNA=rRNA)
        models.changeName('RBSCalc_v2_1',rRNA)
        modelNames.append(rRNA)

    models.setform(modelNames, x="dG_total", y="PROT.MEAN", yScale='ln', a1=-0.45)

    ModelTestSystem = synbiomts.analyze.ModelTest(models,dbfilename,filters,verbose=False)
    ModelTestSystem.run()

    for modelName in modelNames:
        ModelTestSystem.compare2models(modelNames=[modelName,'CCUCCU'])

    # Write model predictions and statistics to Excel
    with open("labels/labels1.txt","r") as f:
        predictLabels = [x.strip('\n') for x in f.readlines()]
    with open("labels/labels_stats.txt","r") as f:
        statsLabels = [x.strip('\n') for x in f.readlines()]

    ModelTestSystem.to_excel('RBS_Calc_16SrRNA_diff_lengths',predictLabels,statsLabels)

def test_Bthetaiotaomicron_hypotheses():

    ### === Bacteroides thetaiotaomicron hypothesis testing === ###

    filters = {"DATASET": ['Mimee_Cell_Sys_2015']}
    models = synbiomts.interface.Container()

    # default RBS Calculator model
    models.add(RBSCalc_v2_1,rRNA='ACCUCCUUU')
    models.changeName('RBSCalc_v2_1','RBSCalc_v2_1-9ntaSD')

    # RBS Calculator using 13 nt 16S rRNA sequence
    models.add(RBSCalc_v2_1,rRNA='GAACACCUCCUUU')
    models.changeName('RBSCalc_v2_1','RBSCalc_v2_1-13ntaSD')

    # Prediction using only ddG_structure
    models.add(StructureOnlyModel)

    modelNames = ['RBSCalc_v2_1-9ntaSD','RBSCalc_v2_1-13ntaSD','StructureOnlyModel']
    models.setform(modelNames, x="dG_total", y="PROT.MEAN", yScale='ln', a1=-0.45)

    dbfilename = '../geneticsystems.db'
    ModelTestSystem = synbiomts.analyze.ModelTest(models,dbfilename,filters,verbose=True)
    ModelTestSystem.run()    

    for modelName in modelNames:
        ModelTestSystem.compare2models(modelNames=[modelName,'RBSCalc_v2_1-9ntaSD'])

    # Write model predictions and statistics to Excel
    with open("labels/labels1.txt","r") as f:
        predictLabels = [x.strip('\n') for x in f.readlines()]
    with open("labels/labels_stats.txt","r") as f:
        statsLabels = [x.strip('\n') for x in f.readlines()]

    ModelTestSystem.to_excel('RBSCalc_16SrRNA_Btheta',predictLabels,statsLabels)

def test_Bthetaiotaomicron_hypotheses_043019():

    ### === Bacteroides thetaiotaomicron hypothesis testing === ###

    filters = {"DATASET": ['Mimee_Cell_Sys_2015']}
    models = synbiomts.interface.Container()

    # default RBS Calculator model
    models.add(RBSCalc_v2_1,rRNA='ACCUCCUUU')
    models.changeName('RBSCalc_v2_1','RBSCalc_v2_1_default')

    # RBS Calculator using 13 nt 16S rRNA sequence
    models.add(RBSCalc_v2_1,rRNA='ACCUCUUCC')
    models.changeName('RBSCalc_v2_1','RBSCalc_v2_1_alt')

    # Prediction using only ddG_structure
    models.add(StructureOnlyModel)

    modelNames = ['RBSCalc_v2_1_default','RBSCalc_v2_1_alt','StructureOnlyModel']
    models.setform(modelNames, x="dG_total", y="PROT.MEAN", std="PROT.STD", yScale='ln', a1=-0.45)

    dbfilename = 'geneticsystems.db'
    ModelTestSystem = synbiomts.analyze.ModelTest(models,dbfilename,filters,verbose=True)
    ModelTestSystem.run()    

    for modelName in modelNames:
        ModelTestSystem.compare2models(modelNames=[modelName,'RBSCalc_v2_1_default'])

    for modelName in modelNames:
        ModelTestSystem.compare2models(modelNames=[modelName,'StructureOnlyModel'])

    # Write model predictions and statistics to Excel
    with open("labels/labels1.txt","r") as f:
        predictLabels = [x.strip('\n') for x in f.readlines()]
    with open("labels/labels_stats.txt","r") as f:
        statsLabels = [x.strip('\n') for x in f.readlines()]

    ModelTestSystem.to_excel('RBSCalc_16SrRNA_Btheta_043019',predictLabels,statsLabels)

if __name__ == "__main__":
    # test_Ecoli_hypotheses()
    test_Bthetaiotaomicron_hypotheses_043019()