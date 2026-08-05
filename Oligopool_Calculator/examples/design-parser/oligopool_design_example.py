import sys
import os
import time as tt
import pandas as pd
import pprint  as pp

import collections.abc

#hyper needs the four following aliases to be done manually.
collections.Iterable = collections.abc.Iterable
collections.Mapping = collections.abc.Mapping
collections.MutableSet = collections.abc.MutableSet
collections.MutableMapping = collections.abc.MutableMapping

sys.path.append('../../oligopool')
from oligopool import designparser

def get_promoter_list(variant_file):
    with open(variant_file) as infile:
        plist = [x.strip() for x in infile.readlines()]
    return plist

def main(variant_file):

    plist = get_promoter_list(variant_file)

    results = designparser(

        pool_size=len(plist),

        element_names=[
            'Primer1',
            'Cut1',
            'Promoter_Variant',
            'Cut2',
            'spacer',
            'Cut3',
            'Barcode',
            'Cut4',
            'Primer2',
            'Filler'],

        elements_spec={
            'Primer1': {
                        'type': 'primer',
                  'oligolimit': 250,
                   'primerseq': 'WWNNNNNNNNNNNNNNSS',
                  'primertype': 0,
                    'mintmelt': 53,
                    'maxtmelt': 62,
                   'maxreplen': 15, 
                'pairedprimer': 'Primer2',
                 'leftcontext': None,
                'rightcontext': 'Cut1',
                    'exmotifs': ['GAGCTC', 'GTCGAC', 'GGATCC', 'GCTAGC'] + ['CCCCC', 'AAAAA', 'TTTTT', 'GGGGG']},

            'Cut1': {
                        'type': 'motif',
                  'oligolimit': 250,
                    'motifseq': 'GAGCTC',
                 'leftcontext': 'Primer1',
                'rightcontext': 'Promoter_Variant',
                    'exmotifs': ['GAGCTC', 'GTCGAC', 'GGATCC', 'GCTAGC'] + ['CCCCC', 'AAAAA', 'TTTTT', 'GGGGG']},

            'Promoter_Variant': {
                     'type': 'variant',
                'sequences': plist},

            'Cut2': {
                        'type': 'motif',
                  'oligolimit': 250,
                    'motifseq': 'GTCGAC',
                 'leftcontext': 'Promoter_Variant',
                'rightcontext': 'spacer',
                    'exmotifs': ['GAGCTC', 'GTCGAC', 'GGATCC', 'GCTAGC'] + ['CCCCC', 'AAAAA', 'TTTTT', 'GGGGG']},

            'spacer': {
                        'type': 'spacer',
                  'oligolimit': 250,
                   'spacerlen': 6, 
                 'leftcontext': 'Cut2',
                'rightcontext': 'Cut3',
                    'exmotifs': ['GAGCTC', 'GTCGAC', 'GGATCC', 'GCTAGC'] + ['CCCCC', 'AAAAA', 'TTTTT', 'GGGGG']},

            'Cut3': {
                        'type': 'motif',
                  'oligolimit': 250,
                    'motifseq': 'GGATCC',
                 'leftcontext': 'spacer',
                'rightcontext': 'Barcode',
                    'exmotifs': ['GAGCTC', 'GTCGAC', 'GGATCC', 'GCTAGC'] + ['CCCCC', 'AAAAA', 'TTTTT', 'GGGGG']},

            'Barcode': {
                        'type': 'barcode',
                  'oligolimit': 250,
                  'barcodelen': 15, 
                    'minhdist': 3,
                   'maxreplen': 5,
                 'barcodetype': 1,
                    'exmotifs': ['GAGCTC', 'GTCGAC', 'GGATCC', 'GCTAGC'] + ['CCCCC', 'AAAAA', 'TTTTT', 'GGGGG'],
                 'leftcontext': 'Cut3',
                'rightcontext': 'Cut4'},

            'Cut4': {
                        'type': 'motif',
                  'oligolimit': 250,
                    'motifseq': 'GCTAGC',
                 'leftcontext': 'Barcode',
                'rightcontext': 'Primer2',
                    'exmotifs': ['GAGCTC', 'GTCGAC', 'GGATCC', 'GCTAGC'] + ['CCCCC', 'AAAAA', 'TTTTT', 'GGGGG']},

            'Primer2': {
                        'type': 'primer',
                  'oligolimit': 250,
                   'primerseq': 'SSNNNNNNNNNNNNNNWW',
                  'primertype': 1,
                    'mintmelt': 53,
                    'maxtmelt': 62,
                   'maxreplen': 15,
                'pairedprimer': 'Primer1',
                 'leftcontext': None,
                'rightcontext': 'Cut1',
                    'exmotifs': ['GAGCTC', 'GTCGAC', 'GGATCC', 'GCTAGC'] + ['CCCCC', 'AAAAA', 'TTTTT', 'GGGGG']},

            'Filler': {
                        'type': 'spacer',
                  'oligolimit': 250,
                   'spacerlen': None,
                 'leftcontext': 'Primer2',
                'rightcontext': None,
                    'exmotifs': ['GAGCTC', 'GTCGAC', 'GGATCC', 'GCTAGC'] + ['CCCCC', 'AAAAA', 'TTTTT', 'GGGGG']},
        },

        background_spec={
               'indata': plist,
            'maxreplen': 12},

        split_spec= None, #{
#             'splitlimit':None, #100
#             'mintmelt':40,
#             'minhdist':1,
#             'minoverlap':20, #20
#             'maxoverlap':50}, #50

        padding_spec=None, #{
#               'typeIIS': 'BsaI',
#               'oligolimit': 250,
#               'mintmelt': 40,
#               'maxtmelt': 80,
#               'maxreplen': 15,}
     )
    
    print(results)
    if results['success']:
        output = results['output']
        stats_dict = results['stats_dict']
        print(output['annotated_complete'])
        
        oligos_df = pd.DataFrame.from_dict(data=output['annotated_complete'],orient='index')
    
        print(oligos_df.head())
    
        oligos_df.to_csv('./oligopool_full_sequences.csv')

    else:
        print('Design failed at Step #{}'.format(results['step']))
        print(results['step_name'])

if __name__ == '__main__':

  variant_file = './oligopool_design_variants.txt'
  main(variant_file)