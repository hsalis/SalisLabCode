import sys
import os
​
sys.path.append(
    os.path.abspath(
        os.path.join('..', 'oligopool')))
​
import index
import pack
import xcount
​
​
def main():
​
    '''
    Step 1: Index all Barcodes to be Mapped
    '''
​
    # Some Vars we Reuse in this Example
    indexfiles = ('Ribozymes-BC1','Ribozymes-BC2')
    conditions = ('MB', 'CB')
    replicates = (1, 2)
​
​
​
    # Index the First Barcode - BC1
    index.index(
        barcodedata='Ribozymes.csv', # File Name Storing Barcode Information [Required]
        barcodecol='BC1',            # Column Name for Barcode Sequence [Reqd]
        barcodeprefix='BC1-Left',    # Column Name for Barcode Prefix Constant [Provide Either]
        barcodesuffix=None,          # Column Name for Barcode Suffix Constant [or Both of These]
        indexfile=indexfiles[0],     # File Name for Resulting Index File [Reqd]
        barcodepregap=0,             # Gap Length between Prefix and Barcode [Optional]
        barcodepostgap=0,            # Gap Length between Barcode and Suffix [Opt]
        associatedata=None,          # File Name Storing Associate Information [Opt]
        associatecol=None,           # Column Name for Associate Sequence [Provide these]
        associateprefix=None,        # Column Name for Assocaite Prefix Constant [if Assicate]
        associatesuffix=None,        # Column Name for Associate Suffix Constant [Data Present]
        verbose=True)                # Verbosity Control
​
    # Index the Second Barcode - BC2
    index.index(
        barcodedata='Ribozymes.csv', # File Name Storing Barcode Information [Required]
        barcodecol='BC2',            # Column Name for Barcode Sequence [Reqd]
        barcodeprefix='BC2-Left',    # Column Name for Barcode Prefix Constant [Provide Either]
        barcodesuffix='BC2-Right',   # Column Name for Barcode Suffix Constant [or Both of These]
        indexfile=indexfiles[1],     # File Name for Resulting Index File [Reqd]
        barcodepregap=0,             # Gap Length between Prefix and Barcode [Optional]
        barcodepostgap=0,            # Gap Length between Barcode and Suffix [Opt]
        associatedata=None,          # File Name Storing Associate Information [Opt]
        associatecol=None,           # Column Name for Associate Sequence [Provide these]
        associateprefix=None,        # Column Name for Assocaite Prefix Constant [if Assicate]
        associatesuffix=None,        # Column Name for Associate Suffix Constant [Data Present]
        verbose=True)                # Verbosity Control
​
    '''
    Step 2: Pack all FASTQ Files to be Counted
    '''
​
    # Serially Pack Each Experiment
    for condition in conditions:
        for replicate in replicates:
​
            # Build out File Names
            r1file   = f'Reads/{condition}{replicate}_R1_001.fastq.gz'
            r2file   = f'Reads/{condition}{replicate}_R2_001.fastq.gz'
            packfile = f'Packs/{condition}{replicate}'
​
            # Execute Packing
            pack.pack(
                r1file=r1file,     # File Name Storing R1 Reads [Reqd]
                r1type=0,          # Orientation of R1 Reads, 0=Forward [Reqd]
                r2file=r2file,     # File Name Storing R2 Reads [Opt]
                r2type=1,          # Orientation of R2 Reads, 1=Reverse [Provide if R2 File Present]
                packtype=1,        # 0=Store Reads, 2=Merge Reads [Reqd]
                packfile=packfile, # File Name to Store Packed Reads [Reqd]
                r1length=140,      # Min. Acceptable R1 Length [Reqd]
                r1qual=30,         # Min. Acceptable R1 Q-Score [Reqd]
                r2length=140,      # Min. Acceptable R2 Length [Opt]
                r2qual=30,         # Min. Acceptable R2 Q-Score [Provide if R2 File Present]
                packsize=2.0,      # No. of Reads in Each Pack in Millions [Opt]
                ncores=0,          # No. of Cores to be Used, 0=Auto [Opt]
                memlimit=0,        # Max. Amount of RAM per Core in GBs, 0=Auto [Opt]
                verbose=True)      # Verbosity Control
​
    '''
    Step 3: Count all Packed Reads (Step 2) using Index (Step 1)
    '''
​
    # Custom Analysis Function
    def myfunc(read, ID, count, coreid):
        '''
        A custom read processing function that is
        executed concurrently with counting, which
        performs additional analysis, and returns
        a boolean indicating whether a read should
        be accepted or not.
​
        These functions are useful when one needs
        to extract additional information from the
        reads being counted, or specify additional
        criteria for accepting a read for counting.
​
        All custom counting functions must at least
        accept the following arguments.
​
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
​
    # Serially Count Each Experiment
    for condition in conditions:
        for replicate in replicates:
​
            # Build File Names
            packfile  = f'Packs/{condition}{replicate}'
            countfile = f'Count/{condition}{replicate}'
​
            # Execute Counting
            xcount.xcount(
                indexfiles=indexfiles, # Tuple of Index Files for Mapping [Reqd]
                packfile=packfile,     # Current Packfile to be Counted [Reqd]
                countfile=countfile,   # File Name to Store Counts [Reqd]
                maptype=1,             # 0=Fast/Approx Mapping, 1=Exact Mapping [Opt]
                barcodeerrors=3,       # Max. No. Seq Errors in Barcodes to Tolerate [Opt]
                callback=myfunc,       # Custom Callback Analysis Function [Opt]
                ncores=10,             # No. of Cores to be Used, 0=Auto [Opt]
                memlimit=0,            # Max. Amount of RAM per Core in GBs, 0=Auto [Opt]
                verbose=True)          # Verbosity Control
​
if __name__ == '__main__':
    main()
