import sys
import os

sys.path.append(
    os.path.abspath(
        os.path.join('../../', 'oligopool')))

import index
import pack
import xcount


def main():

    # Custom Analysis Function
    def myfunc(read, ID, count, coreid):
        '''
        A custom read processing function that is
        executed concurrently with counting, which
        performs additional analysis, and returns
        a boolean indicating whether a read should
        be accepted or not.

        These functions are useful when one needs
        to extract additional information from the
        reads being counted, or specify additional
        criteria for accepting a read for counting.

        All custom counting functions must at least
        accept the following arguments.

        :: read
           type - string
           desc - the read being counted / analyzed
        :: ID
           type - tuple
           desc - a tuple of IDs, one for each element
                  mapped from each of the given indexes
                  e.g. ('Index-1-Barcode-47',
                        'Index-2-Barcode-12',
                        '-',
                        'Index-4-Barcode-77)
                  is a potential ID for a given read
                  when four indexes are specified for
                  counting, and none of the elements
                  from the third index got mapped.
                  The '-' indicates a missing value.
        :: count
           type - integer
           desc - the associated read / ID frequency
                  count in a given read pack
        :: coreid
           type - integer
           desc - the coreid integer identifier that
                  processed this read
        '''
        return True # Accept everything, for now.

if __name__ == '__main__':
    
    index_file_name = 'analysis_barcode_index'
    '''   Step 1: Index barcodes   '''
    index_stats = index.index(
        barcodedata='./oligopool_analysis_example_barcodes.csv', # File Name Storing Barcode Information [Required]
        barcodecol='Barcode',            # Column Name for Barcode Sequence [Reqd]
        barcodeprefix='Prefix',       # Column Name for Barcode Prefix Constant [Provide Either]
        barcodesuffix='Suffix',          # Column Name for Barcode Suffix Constant [or Both of These]
        indexfile=index_file_name,     # File Name for Resulting Index File [Reqd]
        barcodepregap=0,             # Gap Length between Prefix and Barcode [Optional]
        barcodepostgap=0,            # Gap Length between Barcode and Suffix [Opt]
        associatedata=None,          # File Name Storing Associate Information [Opt]
        associatecol=None,           # Column Name for Associate Sequence [Provide these]
        associateprefix=None,        # Column Name for Assocaite Prefix Constant [if Assicate]
        associatesuffix=None,        # Column Name for Associate Suffix Constant [Data Present]
        verbose=True)                # Verbosity Control

    indexfile = index_file_name
    
            
    '''   Step 2: Pack FASTQ Files to be counted   '''

    # Build out File Names
    r1file   = f'./analysis_100k_reads_R1.fastq.gz'
    r2file   = f'./analysis_100k_reads_R2.fastq.gz'
    
    packfile = f'Pack/{index_file_name}'
    countfile = f'Count/{index_file_name}'

    
    pack.pack(
        r1file=r1file,     # File Name Storing R1 Reads [Reqd]
        r1type=0,          # Orientation of R1 Reads, 0=Forward [Reqd]
        r2file=r2file,     # File Name Storing R2 Reads [Opt]
        r2type=1,          # Orientation of R2 Reads, 1=Reverse [Provide if R2 File Present]
        packtype=0,        # 0=Store Reads, 1=Merge Reads [Reqd]
        packfile=packfile, # File Name to Store Packed Reads [Reqd]
        r1length=140,      # Min. Acceptable R1 Length [Reqd]
        r1qual=30,         # Min. Acceptable R1 Q-Score [Reqd]
        r2length=140,      # Min. Acceptable R2 Length [Opt]
        r2qual=30,         # Min. Acceptable R2 Q-Score [Provide if R2 File Present]
        packsize=2.0,      # No. of Reads in Each Pack in Millions [Opt]
        ncores=0,          # No. of Cores to be Used, 0=Auto [Opt]
        memlimit=0,        # Max. Amount of RAM per Core in GBs, 0=Auto [Opt]
        verbose=True)      # Verbosity Control
        




    xcount.xcount(
    indexfiles=indexfile, # Tuple of Index Files for Mapping [Reqd]
    packfile=packfile,     # Current Packfile to be Counted [Reqd]
    countfile=countfile,   # File Name to Store Counts [Reqd]
    maptype=1,             # 0=Fast/Approx Mapping, 1=Exact Mapping [Opt]
    barcodeerrors=3,       # Max. No. Seq Errors in Barcodes to Tolerate [Opt]
    callback=None,       # Custom Callback Analysis Function [Opt]
    ncores=0,             # No. of Cores to be Used, 0=Auto [Opt]
    memlimit=0,            # Max. Amount of RAM per Core in GBs, 0=Auto [Opt]
    verbose=True)          # Verbosity Control