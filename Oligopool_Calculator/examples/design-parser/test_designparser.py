import sys
import os

import time    as tt
import pandas  as pd

import designparser as dp

import pprint  as pp

def get_promoter_list():
    with open('promoters.txt') as infile:
        plist = [x.strip() for x in infile.readlines()]
    return plist

def main():

    plist = get_promoter_list()

    output = dp.designparser(

        pool_size=len(plist),

        element_names=[
            'Primer1',
            'Cut1',
            'Promoter',
            'Barcode',
            'Primer2',
            'Cut2',
            'Primer3',
            'Filler'],

        elements_spec={
            'Primer1': {
                        'type': 'primer',
                  'oligolimit': 170,
                   'primerseq': 'NNNNNNNNNNNNNNNNNNNNNN',
                  'primertype': 0,
                    'mintmelt': 53,
                    'maxtmelt': 55,
                   'maxreplen': 10,
                'pairedprimer': 'Primer3',
                 'leftcontext': None,
                'rightcontext': 'Cut1',
                    'exmotifs': ['GGATCC', 'TCTAGA'] + ['CCCCC', 'AAAAA', 'TTTTT', 'GGGGG']},

            'Cut1': {
                        'type': 'motif',
                  'oligolimit': 194,
                    'motifseq': 'NNNGGATCCNNN',
                 'leftcontext': 'Primer1',
                'rightcontext': 'Promoter',
                    'exmotifs': ['GGATCC', 'TCTAGA'] + ['CCCCC', 'AAAAA', 'TTTTT', 'GGGGG']},

            'Promoter': {
                     'type': 'variant',
                'sequences': plist},

            'Barcode': {
                        'type': 'barcode',
                  'oligolimit': 170,
                  'barcodelen': 15,
                    'minhdist': 3,
                   'maxreplen': 5,
                 'barcodetype': 1,
                    'exmotifs': ['GGATCC', 'TCTAGA'] + ['CCCCC', 'AAAAA', 'TTTTT', 'GGGGG'],
                 'leftcontext': 'Promoter',
                'rightcontext': 'Primer2'},

            'Primer2': {
                        'type': 'primer',
                  'oligolimit': 170,
                   'primerseq': 'NNNNNNNNNNNNNNNNNNNNNN',
                  'primertype': 0,
                    'mintmelt': 53,
                    'maxtmelt': 55,
                   'maxreplen': 10,
                'pairedprimer': 'Primer3',
                 'leftcontext': 'Barcode',
                'rightcontext': 'Cut2',
                    'exmotifs': ['GGATCC', 'TCTAGA'] + ['CCCCC', 'AAAAA', 'TTTTT', 'GGGGG']},

            'Cut2': {
                        'type': 'motif',
                  'oligolimit': 194,
                    'motifseq': 'NNNTCTAGANNN',
                 'leftcontext': 'Primer2',
                'rightcontext': 'Primer3',
                    'exmotifs': ['GGATCC', 'TCTAGA'] + ['CCCCC', 'AAAAA', 'TTTTT', 'GGGGG']},

            'Primer3': {
                        'type': 'primer',
                  'oligolimit': 170,
                   'primerseq': 'NNNNNNNNNNNNNNNNNNNNNN',
                  'primertype': 1,
                    'mintmelt': 53,
                    'maxtmelt': 55,
                   'maxreplen': 10,
                'pairedprimer': 'Primer2',
                 'leftcontext': 'Cut2',
                'rightcontext': 'Filler',
                    'exmotifs': ['GGATCC', 'TCTAGA'] + ['CCCCC', 'AAAAA', 'TTTTT', 'GGGGG']},

            'Filler': {
                        'type': 'spacer',
                  'oligolimit': 170,
                   'spacerlen': None,
                 'leftcontext': 'Primer3',
                'rightcontext': None,
                    'exmotifs': ['GGATCC', 'TCTAGA'] + ['CCCCC', 'AAAAA', 'TTTTT', 'GGGGG']},
        },

        background_spec={
               'indata': plist,
            'maxreplen': 12},

        split_spec={
            'splitlimit':100,
              'mintmelt':40,
              'minhdist':1,
            'minoverlap':20,
            'maxoverlap':50},

        padding_spec={
               'typeIIS': 'bsaI',
            'oligolimit': 170,
              'mintmelt': 40,
              'maxtmelt': 80,
             'maxreplen': 15,}
    )

    pp.pprint(output['stats_dict'])

    # for df in output:
    #     print(output[df])

if __name__ == '__main__':
    main()