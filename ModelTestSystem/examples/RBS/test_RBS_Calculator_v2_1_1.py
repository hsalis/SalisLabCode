
import synbiomts
import cPickle as pickle

# We're going to use shelve to store model predictions
# as a dictionary-like persistance object of pandas dataframes
import shelve

# Import models
import sys
# sys.path.append('models')

import numpy as np

import copy

# Import private Salis Lab code (latest versions of RBS Calculator)
sys.path.append('/home/alex/Private-Code')
from DNAc import *

from PyVRNA import PyVRNA
RNAEnergyModel = PyVRNA(dangles=0)

handle = open('models/rRNA_16S_3p_ends.p','r')
rRNA_16S_3p_ends = pickle.load(handle)
handle.close()

# Required function in RBS Calculators
def get_rRNA(organism):
    # temporary until I update database:
    if   organism=="Corynebacterium glutamicum B-2784": organism = 'Corynebacterium glutamicum R'
    elif organism=="Pseudomonas fluorescens A506":      return 'ACCTCCTTT'
    elif organism=="Escherichia coli BL21(DE3)":        return 'ACCTCCTTA'
    else: pass
    return rRNA_16S_3p_ends[organism]

import RBS_Calculator_v2_1_1

def RBSCalc_v2_1(sequence,organism,temp,startpos):

    # instance RBS Calculator
    model = RBS_Calculator_v2_1_1.RBS_Calculator(canonicalStartsOnly=False,temp=temp,rRNA_seq=get_rRNA(organism),organism=organism)

    # calculate translation rates of mRNA
    result  = model.calc_transl_rates(sequence,start_range=[0,startpos])
    
    best = None
    for RBS in result.RBS_list:
        RBS.get_results_as_dict()
        if (startpos - RBS.position) % 3 != 0:
            continue
        if best is None or RBS.dG_total < best.dG_total:
            best = RBS

    results = best.get_results_as_dict()
    mRNA = result.sequence
    results['predicted_UTR'] = mRNA[:best.position]
    results['predicted_CDS'] = mRNA[best.position:]

    return results

if __name__ == "__main__":

    # add models to interface.Models
    models = synbiomts.interface.Container()
    models.add(model=RBSCalc_v2_1)
    models.setform(['RBSCalc_v2_1'], x='dG_total', y='PROT.MEAN', std='PROT.STD', yScale='ln', a1=-0.45)
    # models.setform(['RBSCalc_v2_1_1'], x='dG_total', y='PROT.MEAN', std='PROT.STD', yScale='ln')

    # Provide the pickled database file name
    dbfilename = 'geneticsystems.db'

    filters = {"DATASET":
        [
        'EspahBorujeni_NAR_2013',
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
        'Beck_PLoS_2016',
        ]
    }

    # ModelTestSystem = synbiomts.analyze.ModelTest(models,dbfilename,filters,nprocesses=1,add_data=True,verbose=True)
    ModelTestSystem = synbiomts.analyze.ModelTest(models,dbfilename,filters,add_data=True,verbose=True)
    ModelTestSystem.run()

    # Write model predictions and statistics to Excel
    with open("labels/labels3.txt","r") as f:
        predictLabels = [x.strip('\n') for x in f.readlines()]

    with open("labels/labels_stats.txt","r") as f:
        statsLabels = [x.strip('\n') for x in f.readlines()]

    ModelTestSystem.to_excel('RBS_Calculator_v2_1_1_IC1014',predictLabels,statsLabels)

