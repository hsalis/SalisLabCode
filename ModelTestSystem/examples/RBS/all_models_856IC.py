
import sys
sys.path.append('../')
sys.path.append('../models')
sys.path.append('../datasets')
sys.path.append('/usr/local/lib/python2.7/site-packages/')

import synbiomts
from synbiomts import dbms
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

handle = open('../models/rRNA_16S_3p_ends.p','r')
rRNA_16S_3p_ends = pickle.load(handle)
handle.close()

# Required function in RBS Calculators and UTR Designer
def get_rRNA(organism):
    # temporary until I update database:
    if   organism=="Corynebacterium glutamicum B-2784": organism = 'Corynebacterium glutamicum R'
    elif organism=="Pseudomonas fluorescens A506":    return 'ACCTCCTTT'
    elif organism=="Escherichia coli BL21(DE3)": return 'ACCTCCTTA'
    else: pass    
    return rRNA_16S_3p_ends[organism]

def find_best_start(mRNA, start_pos, predictions):
    '''Finds the number of start codons, the most highly translated start codon and the start range associated with it.
        Disallows start codons that have stop codons following them.'''
    
    stop_codons = ["TAA", "TAG", "TGA", "UAA", "UAG", "UGA"]
    num_start_codons = 0
    dG_tot_list = []
    dG_best = None
    
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
    
    RBS = predictions.RBS_list[best_start_index]
    return (RBS,best_start_pos)

# Load and wrap models
import RBS_Calculator_v1_0
def RBSCalc_v1(sequence,startpos,organism,temp):
    start_range = [0,startpos+1]
    rRNA = get_rRNA(organism)
    model = RBS_Calculator_v1_0.RBS_Calculator(sequence,start_range,rRNA)
    model.temp = temp
    model.run()
    output = model.output()
    (RBS,best_start_pos) = find_best_start(sequence,startpos,output)
    
    results = {
    'TIR' : RBS.tir,
    'dG_total' : RBS.dG_total,
    'dG_mRNA_rRNA' : RBS.dG_mRNA_rRNA,
    'dG_mRNA' : RBS.dG_mRNA,
    'dG_start' : RBS.dG_start,
    'dG_standby' : RBS.dG_standby,
    'dG_spacing' : RBS.dG_spacing,
    'used_mRNA_sequence' : RBS.sequence,
    'warnings' : RBS.warnings,
    '5pUTR' : sequence[:best_start_pos],
    'CDS' : sequence[best_start_pos:]
    }

    return results

import RBS_Calculator_v1_1
def RBSCalc_v1_1(sequence,startpos,organism,temp):
    start_range = [0,startpos+1]
    rRNA = get_rRNA(organism)
    model = RBS_Calculator_v1_1.RBS_Calculator(sequence,start_range,rRNA)
    model.temp = temp
    model.run()
    output = model.output()
    (RBS,best_start_pos) = find_best_start(sequence,startpos,output)

    results = {
    'TIR' : RBS.tir,
    'dG_total' : RBS.dG_total,
    'dG_mRNA_rRNA' : RBS.dG_mRNA_rRNA,
    'dG_mRNA' : RBS.dG_mRNA,
    'dG_start' : RBS.dG_start,
    'dG_standby' : RBS.dG_standby,
    'dG_spacing' : RBS.dG_spacing,
    'used_mRNA_sequence' : RBS.sequence,
    'warnings' : RBS.warnings,
    '5pUTR' : sequence[:best_start_pos],
    'CDS' : sequence[best_start_pos:]
    }

    return results

import RBS_Calculator_v2_0
def RBSCalc_v2(sequence,organism,temp,startpos):
    start_range = [0,startpos+1]
    rRNA = get_rRNA(organism)
    model = RBS_Calculator_v2_0.RBS_Calculator(sequence,start_range,rRNA)
    model.temp = temp
    model.run()
    output = model.output()
    (RBS,best_start_pos) = find_best_start(sequence,startpos,output)
    
    results = {
    'TIR' : RBS.tir,
    'dG_total' : RBS.dG_total,
    'dG_mRNA_rRNA' : RBS.dG_mRNA_rRNA,
    'dG_mRNA' : RBS.dG_mRNA,
    'dG_start' : RBS.dG_start,
    'dG_standby' : RBS.dG_standby,
    'dG_spacing' : RBS.dG_spacing,
    'used_mRNA_sequence' : RBS.sequence,
    'bracket_string_mRNA': RBS.bracket_string_mRNA,
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
    '5pUTR' : sequence[:best_start_pos],
    'CDS' : sequence[best_start_pos:],
    }

    return results

import UTR_Designer
def fUTR_Designer(sequence,startpos,organism,temp):
    start_range = [0,startpos+1]
    rRNA = get_rRNA(organism)
    model = UTR_Designer.RBS_Calculator(sequence,start_range,rRNA)
    model.temp = temp
    model.run()
    output = model.output()
    (RBS,best_start_pos) = find_best_start(sequence,startpos,output)
    
    results = {
    'TIR' : RBS.tir,
    'dG_UTR' : RBS.dG_total,
    'dG_mRNA_rRNA' : RBS.dG_mRNA_rRNA,
    'dG_mRNA' : RBS.dG_mRNA,
    'dG_start' : RBS.dG_start,
    'dG_standby' : RBS.dG_standby,
    'dG_spacing' : RBS.dG_spacing,
    'used_mRNA_sequence' : RBS.sequence,
    'warnings' : RBS.warnings,
    '5pUTR' : sequence[:best_start_pos],
    'CDS' : sequence[best_start_pos:],
    }

    return results

handle = open('../models/RBSDesigner_ALL.p','r')
RBS_Designer = pickle.load(handle)
handle.close()
def fRBS_Designer(sequence):
    
    keys = ['translation efficiency',
    'Expression',
    'exposure probability',
    'mRNA-ribosome probability',
    'SD',
    'dG_SD:aSD',
    'spacer length'
    ]

    if sequence in RBS_Designer: RBS = RBS_Designer[sequence]
    else: RBS = {k: None for k in keys}
    return RBS

from emopec._emopec import predict_spacing, get_expression
def EMOPEC(sequence,startpos):
    RBS = predict_spacing(sequence[:startpos])
    preseq,SD,spacer,Expression = RBS
    maxExpression = get_expression(sd_seq='AGGAGA',sd_dist=len(spacer))
    minExpression = get_expression(sd_seq='TTGGGC',sd_dist=len(spacer))
    Expression_percent = (Expression - minExpression)/(maxExpression - minExpression)

    results = {
    'preseq' : preseq,
    'SD_seq' : SD,
    'spacer_seq' : spacer,
    'spacer_length' : len(spacer),
    'Expression' : Expression,
    'Expression_percent': Expression_percent
    }

    return results

def _calc_warning_flags(RBS):
    # Calculate warning flags identified in Reis & Salis, 2017
    # RNA structure overlapping the RBS > 17 nt long
    
    warning_dict = {}

    # ----------------------------------------------------------------------
    # Check for long structure that at least partially occludes the RBS
    cutoff = 17
    bpx = RBS.initial_structure['bp_x']
    bpy = RBS.initial_structure['bp_y']
    RBSposi = RBS.most_3p_paired_SD_mRNA
    RBSposj = RBS.start_position + 13
    y0 = -1
    warning_dict["RBS_hairpin_lengths"] = []
    for x,y in zip(bpx,bpy):
        if x > y0:
            if x<RBSposi and y>RBSposi:
                warning_dict["RBS_hairpin_lengths"].append(y-x)
                y0 = y
    if any(length > cutoff for length in warning_dict["RBS_hairpin_lengths"]):
        warning_dict["LONG_RBS_HAIRPIN"] = True
    else:
        warning_dict["LONG_RBS_HAIRPIN"] = False

    # ----------------------------------------------------------------------
    # Check for stable refolding of standby site structure (ddG_standby < -2 kcal/mol)
    # Where ddG_standby_site = dG_pre_ribosome - dG_standby_site_initial
    cutoff = -2.0
    y0 = -1
    dG_initial_list = []
    for x,y in zip(bpx,bpy):
        if x > y0:
            PyVRNA_fold_result = RNAEnergyModel.RNAfold(RBS.sequence[x:y+1])
            dG_initial_list.append(PyVRNA_fold_result.energy)
            y0 = y
    dG_standby_site_initial = sum(dG_initial_list)
    dG_pre_ribosome = RNAEnergyModel.RNAfold(RBS.sequence[x:y+1]).energy
    ddG_standby_site = dG_pre_ribosome - dG_standby_site_initial
    warning_dict["ddG_standby_site"] = ddG_standby_site
    if ddG_standby_site < cutoff:
        warning_dict["STABLE_STANDBY_SITE_REFOLD"] = True
    else:
        warning_dict["STABLE_STANDBY_SITE_REFOLD"] = False

    # ----------------------------------------------------------------------
    # Check for sequences with minimal ssRNA regions (ssRNA < 5 nt)
    cutoff = 5
    warning_dict["ssRNA_lengths"] = []
    y0 = -1
    for x,y in zip(bpx,bpy):
        if x > y0:
            warning_dict["ssRNA_lengths"].append(x-y0)
            y0 = y
    if not any(length > cutoff for length in warning_dict["ssRNA_lengths"]):
        warning_dict["MINIMAL_ssRNA"] = True
    else:
        warning_dict["MINIMAL_ssRNA"] = False

    return warning_dict


def main():

    # Add models to interface.Container()
    models = synbiomts.interface.Container()
    models.add(RBSCalc_v1)
    models.add(RBSCalc_v1_1)
    models.add(RBSCalc_v2)
    models.add(fUTR_Designer)
    models.add(fRBS_Designer)
    models.add(EMOPEC)

    # Define linear form between predictor and outcome (genetic system function)
    models.setform(["RBSCalc_v1","RBSCalc_v1_1","RBS_Calc_v2"], x="dG_total", y="PROT.MEAN", yScale='ln', a1=-0.45)
    models.setform(["UTRDesigner"], x="dG_UTR", y="PROT.MEAN", yScale='ln', a1=-0.45)
    models.setform(["RBSDesigner", "EMOPEC"], x="Expression", y="PROT.MEAN")

    # Set database filters
    filters = { "DATASET": ['EspahBorujeni_NAR_2013',
                            'EspahBorujeni_NAR_2015',
                            'EspahBorujeni_JACS_2016',
                            'EspahBorujeni_Footprint',
                            'Salis_Nat_Biotech_2009',
                            'Farasat_MSB_2014',
                            'Tian_NAR_2015',
                            'Mimee_Cell_Sys_2015',
                            'Bonde_NatMethods_IC_2016',
                            'Egbert_Spacers_PNAS_2012']
                }

    # Provide the pickled database file name
    dbfilename = '../geneticsystems.db'

    # customtest = synbiomts.analyze.ModelTest(transl_rate_models,dbfilename,filters,nprocesses=1,verbose=True)
    customtest = synbiomts.analyze.ModelTest(transl_rate_models,dbfilename,filters,nprocesses=10, verbose=True)
    customtest.run(filename='model_calcs.db')

    # Initialize ModelTestSystem and run (automates model and stat calcs)
    ModelTestSystem = synbiomts.analyze.ModelTest(models,dbfilename,filters,verbose=True)
    ModelTestSystem.run(filename='model_calcs.db')

def data_for_RBSDesigner():

    # Export information for RBS Designer
    filters = { "DATASET": ['EspahBorujeni_NAR_2013',
                        'EspahBorujeni_NAR_2015',
                        'EspahBorujeni_JACS_2016',
                        'EspahBorujeni_Footprint',
                        'EspahBorujeni_Bsubtilis_2016',
                        'Salis_Nat_Biotech_2009',
                        'Farasat_MSB_2014',
                        'Tian_NAR_2015',
                        'Mimee_Cell_Sys_2015',
                        'Bonde_NatMethods_IC_2016',
                        'Egbert_PNAS_2012',
                        'Hecht_NAR_2017',
                        'Beck_PLoS_2016']
    }

    dbfilename = '../geneticsystems.db'
    handle = open(dbfilename,'r')
    database = pickle.load(handle)
    handle.close()
    database = dbms.filter(database,filters,False)

    UTR = database['5pUTR'].tolist()
    CDS = database['CDS'].tolist()
    rRNA = [get_rRNA(organism) for organism in database['ORGANISM']]
    assert len(UTR) == len(CDS) == len(rRNA)
    assert type(UTR) == type(CDS) == type(rRNA) == list
    pickle.dump((UTR,CDS,rRNA),open('data4RBSDesigner.p', 'wb'))

if __name__ == "__main__":

    # main()
    data_for_RBSDesigner()
