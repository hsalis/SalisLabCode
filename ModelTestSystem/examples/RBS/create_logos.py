
import cPickle as pickle
from synbiomts import (dbms,graphics)

handle = open('../geneticsystems.db','r')
database = pickle.load(handle)
handle.close()

# 1014IC dataset
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

df = dbms.filter(database,filters,ordered=False) # filter database
sequences = [seq.upper().replace('T','U') for seq in list(df['SEQUENCE'])] # get sequences
starts = list(df['STARTPOS']) # get starts
sequences = [seq[:pos+35] for seq,pos in zip(sequences,starts)] # remove CDS[35:]
maxlen = max([len(seq) for seq in sequences]) # get max sequence length
sequences = ['X'*(maxlen-len(seq))+seq for seq in sequences[:]] # buffer
graphics.create_logo(sequences,'IC1014.eps') # go!

# 6062FS dataset
filters = { "DATASET": ['Kosuri_PNAS_2013',
                        'Goodman_Science_2013']
}

df = dbms.filter(database,filters,ordered=False)
sequences = [seq.upper().replace('T','U') for seq in list(df['SEQUENCE'])]
starts = list(df['STARTPOS'])
sequences = [seq[:pos+35] for seq,pos in zip(sequences,starts)]
maxlen = max([len(seq) for seq in sequences])
sequences = ['X'*(maxlen-len(seq))+seq for seq in sequences[:]]
graphics.create_logo(sequences,'FS6062.eps')