import sys, math, time, random
sys.path.append('/code/SalisLabCode/Oligopool_Calculator/oligopool')

import pandas as pd
import numpy as np

import barcode, primer
from utils import liner_engine
import collections as cx

def generateBarcodeEngine(barcodeLength, minHdist, numberBarcodes):
    maximumSharedRepeat = barcodeLength - minHdist

    emptyContext = ['' for n in range(barcodeLength)]

    liner = liner_engine()
    stats = {}
    oligorepeats = {}
    for n in range(numberBarcodes):
        oligorepeats[n] = ''

    stats = {   'vars'    : {   'targetcount': numberBarcodes,   # Required Number of Barcodes
                                'barcodecount': 0,             # Barcode Design Count
                                'distancefail': 0,             # Hamming Distance Fail Count
                                'repeatfail': 0,             # Repeat Fail Count
                                'exmotiffail': 0,             # Exmotif Elimination Fail Count
                                'exmotifcounter': {}, # Exmotif Encounter Counter
                                'edgefail': 0,             # Edge Effect Fail Count
                                'typefail' : 0,
                                'distancedistro': None,          # Hamming Distance Distribution
                            }
            }
    (barcodeList, store, stats) = barcode.barcode_engine(barcodelen = barcodeLength, minhdist = minHdist, maxreplen = maximumSharedRepeat, barcodetype = 0, oligorepeats = oligorepeats, leftcontext = emptyContext, rightcontext = emptyContext, exmotifs = [], targetcount = numberBarcodes, stats = stats, liner = liner)
    return (barcodeList, stats)

def designBarcodes(barcodeLength, minHdist):

    t1 = time.time()
    numberBarcodes = int(math.pow(10,math.log10(4) * (barcodeLength - minHdist - 2)))
    saveFile = 'barcodes_L{}_H{}.csv'.format(barcodeLength, minHdist)

    (barcodes, stats) = generateBarcodeEngine(barcodeLength, minHdist, numberBarcodes)

    counter = 1
    barcodeList = []
    for barcode in barcodes:
        barcodeList.append({'sequence' : barcode, 'id' : 'barcode_L{}_H{}_{}'.format(barcodeLength, minHdist, counter)})
        counter +=1

    df = pd.DataFrame(barcodeList)
    df.to_csv(saveFile)

    t2 = time.time()
    print("Elapsed Time: ", t2 - t1, " seconds.")
    

def designPrimers():
    
    # Create random DNA sequence variants within the oligopool
    maximum_oligonucleotide_length = 300
    variant_length = 200
    number_variants_list = [100, 1000, 10000, 100000, 1000000]
    
    for number_variants in number_variants_list:
        variant_list = []
        for n in range(number_variants):
            variant = {'ID' : n, 'sequence' : "".join([random.choice(['A','G','C','T']) for l in range(variant_length)]) }
            variant_list.append(variant)
        variant_df = pd.DataFrame(variant_list)
        
        t1 = time.time()
        
        # Forward Primer Design
        primerseq = 'N' * 16 + 'WW'
        primertype = 0  #forward primer
        mintmelt = 50.0 #degC
        maxtmelt = 62.0 #degC
        maxreplen = 15  #nt
        primercol = 'fwd_primer'
        
        (primer_df, stats) = primer.primer(variant_df, maximum_oligonucleotide_length, primerseq, primertype, mintmelt, maxtmelt, maxreplen, primercol,
                                    outfile=None, pairedcol=None, leftcontext=None, rightcontext=None, exmotifs=None, background=None, verbose=True)
        
        print(stats)
        
        # Reverse Primer Design
        primerseq = 'N' * 16 + 'WW'
        primertype = 1  #reverse primer
        mintmelt = 50.0 #degC
        maxtmelt = 62.0 #degC
        maxreplen = 15  #nt
        primercol = 'rev_primer'
        
        (primer_df, stats) = primer.primer(primer_df, maximum_oligonucleotide_length, primerseq, primertype, mintmelt, maxtmelt, maxreplen, primercol,
                                    outfile=None, pairedcol='fwd_primer', leftcontext=None, rightcontext=None, exmotifs=None, background=None, verbose=True)
        
        t2 = time.time()
        print(primer_df)
        print(stats)
        print('Elapsed Time [num_variants = {}]: {} seconds'.format(number_variants, t2 - t1))

if __name__ == "__main__":
    
    ## Primer Benchmarks
    designPrimers()

    
    ## Barcode Benchmarks
    for barcodeLength in range(6, 16):
        t1 = time.time()
        designBarcodes(barcodeLength, 2)
        t2 = time.time()
        print('Barcode Design: Hmin = 2. Length = {} nt. Elapsed Time: {} seconds'.format(barcodeLength, t2 - t1))

    for barcodeLength in range(6, 17):
        t1 = time.time()
        designBarcodes(barcodeLength, 3)
        t2 = time.time()
        print('Barcode Design: Hmin = 3. Length = {} nt. Elapsed Time: {} seconds'.format(barcodeLength, t2 - t1))
    
    
    
    