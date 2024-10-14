# Oligopool Calculator v1.0
### Code written by Ayaan Hossain and Howard Salis

## Relevant Article: Oligopool Calculator: automated design of oligopools and rapid analysis of massively parallel barcoded measurements
### Ayaan Hossain, Daniel P. Cetnar, Travis L. La Fleur, James McLellan, and Howard M. Salis

Oligopool synthesis and next-generation sequencing enable the construction and characterization of large libraries of designed genetic parts and systems. As library sizes grow, it becomes computationally challenging to optimally design large numbers of primer binding sites, barcode sequences, and overlap regions to obtain efficient assemblies and precise measurements. We present the Oligopool Calculator, an end-to-end suite of algorithms and data structures, that rapidly designs many thousands of oligonucleotides within an oligopool and rapidly analyzes many billions of barcoded sequencing reads. We introduce several novel concepts that greatly increase the design and analysis throughput, including orthogonally symmetric barcode design, adaptive decision trees for primer design, a Scry barcode classifier, and efficient read packing. We demonstrate the Oligopool Calculator’s capabilities across computational benchmarks and real-data projects, including the design of over four million highly unique and compact barcodes in 1.2 hours, the design of universal primer binding sites for a million 200-mer oligos in 15 minutes, and the analysis of about 500 million sequencing reads per hour, all on an 8-core desktop computer. Overall, the Oligopool Calculator enables the creative use of massively parallel experiments by eliminating the computational complexity of their design and analysis.


## Installation

Install Python 3, pip, and more:

`apt-get update && apt-get install -y python3 python3-pip build-essential zip wget pypy graphviz nano apt-transport-https ca-certificates curl software-properties-common intel-mkl git-core cmake pandoc memcached mpich mpich-doc libgfortran5 apt-get clean && rm -rf /var/lib/apt/lists/*`

Install the Python module dependencies:

`pip3 install numpy pandas matplotlib seaborn scipy scikit-learn biopython dm-sonnet sonnet graphs mpi4py openpyxl leveldb==0.201 networkx nrpcalc statsmodels msgpack pyfastx edlib dinopy bounter numba primer3-py`

## Design Examples

This example uses the following Python files:
`oligopool_design_example.py
oligopool_design_variants.py`

and generates the following output CSV files:
`oligopool_annotated_sequences.csv
oligopool_full_sequences.csv`

We will used the Oligopool Calculator’s design mode to generate an oligopool containing 1000 sequence variants. Our objective is to design the primers, cut sites, barcodes, and padding around the defined sequence variants to ensure uniform PCR amplification, library-based DNA assembly, and quantification using next-generation sequencing. The example oligopool has the following architecture:

5’ -- [Forward Primer] [Cut1] [Sequence Variant] [Cut2] [spacer] [Cut3] [Barcode] [Cut4] [Filler] [Reverse Primer] -- 3’

Here, a general cloning and characterization strategy is to design the oligopool so that each oligonucleotide contains a unique barcode associated with a desired sequence variant. There are four restriction sites used to carry out a 2-step DNA assembly workflow. The first step introduces the oligopool mixture into the vector, while the second step inserts a large constant region (e.g. a RBS and reporter protein CDS) into the plasmid. This DNA assembly strategy enables each genetic system variant to be trackable and quantifiable using a unique barcode sequence, although the barcode sequence here is located in the 3' UTR region (far away from the designed sequence variant). Padding sequence (Filler) is added to ensure that all PCR amplicons have the same length. Highly specific primers are designed to ensure uniform PCR amplification. Their corresponding primer binding site sequences are inserted into the oligos as flanking regions. During the design of these elements, motif sequences are excluded to facilitate cloning and to prevent homologous recombination (e.g. excluded restriction sites).

The Python code ‘oligopool_design_example.py’ demonstrates the workflow for designing the oligopool described above. The asset file ‘oligopool_design_variants.txt’ contains a list of 1000 example promoter variants.

The oligopool design pipeline takes place in the designpasrser.designparser() function (Appendix A). In this script, we will go through each item and explain their role in this context. The ‘result’ variable stores the series of nested dictionaries that comprise the designparser() output. Here, designparser() takes the parameters pool_size, element_names, element_spec, and background_spec.

pool_size defines the size of the oligopool. This will be the number of variants we have, so  we simply store the text file that contains the promoter sequences as a list, and use its length as the pool_size.  

• element_names is a list that defines the names of each element in each sequence as dictionaries in the order of appearance. Note that in the provided script, this list matches the sequence design above. The defined elements are ‘Primer1’, ‘Cut1’, ‘Promoter_Variant’, ‘Cut2’, ‘spacer’, ‘Cut3’, ‘Barcode’, ‘Cut4’, ‘Filler’. 

• element_spec is a dictionary of element dictionaries where each key is a name from element_names. For each element, there is a key ‘type’ that will dictate what other keys are required within the element dictionary. See Appendix B for details on the valid keys compatible for each type in element_spec. For this example, we will go through the each element defined in the script, and note its function and important design decisions. The element types are: barcode, primer, motif, variant, and spacer. 

•	‘Barcode’ element type barcode, which defines parameters to design the unique barcodes for each promoter variant, such that no oligo two sequences will have the same barcode. This will be necessary for the analysis pipeline, where barcodes are counted in NGS data and mapped them to each unique promoter variant. Here, we’ve defined barcode length as 15, the max repetitive length (maxrepel) as 5 to ensure the barcodes are non-repetitive, and the minimum pairwise hamming distance (minhdist) as 3, which ensures that the barcodes are sufficiently different. 

•	‘Primer1’ and ‘Primer2’ are of type primer, which includes specifications for the melting temperature range and sequence constraints, including specifications for degenerate nucleotides. Importantly, ‘pairedprimer’ will specify which primers are paired for PCR together, which is accounted for in their design to ensure optimal melting temperature matching. At the 3’ end of each primer, we’ve specified two W nts to minimize dimer amplification.

•	The ‘Cut’ elements are defined as type motif.  The motif type allows a defined sequence (motifseq) to be defined and placed in each sequence the oligopool. Using the motif type will ensure the desired restriction enzyme cutsites are present in all designed sequences.

•	‘Promoter_Variant’ is of type variant. The variant type is where you will define the variants desired in your oligopool. All sequence variants you want included in your design should be stored as a list in the sequence key. Here, this is the list of promoters that were imported from the text file. 

•	‘spacer’ and ‘Filler’ elements are of type spacer, which are sequences generated by the Oligopool Calculator according to the constraints supplied by the user and inherent to the oligopool design, such as checking for excluded motifs and undesired primer binding sites. Note that ‘Filler’ has length ‘None’, which allows designparser() to fill out the remaining space for each sequence based on the defined ‘oligolimit’.

There are several common keys for each element type (except for variant): 
•	exmotifs are sequences the user defines to be excluded from the designed element. If you want to exclude all these sequences from your design oligopool, they will need to be specified in each element. Since we want to ensure that our four restriction cutsites in our sequence designs are unique, we’ve included each of these in the exmotifs of all elements.
•	oligolimit is the length of the oligos in your oligopool. This should be the same for all elements. 
•	leftcontext and rightcontext are the elements directly upstream (leftcontext) and downstream (rightcontext) of each element. 

The output of designparser() is a dictionary (in this case, titled result), with the keys, step, step_name, stats_dict, and output. The former 3 provide information on the design process and may be used for troubleshooting; refer to the documentation for a detailed explanation. The output key stores two dictionaries with your oligopool sequences stored in two keys, (1) presplit_complete is a dataframe which includes includes the full sequences indexed by a unique sequence ID, and (2) annotated_complete contains the same sequences as a dataframe, with each element separated into its own column. Our script saves these dataframes as separate files ‘oligopool_annotated_sequences.csv’ and ‘oligopool_full_sequences.csv’, so they may be accessed for ordering and further inspection. 

## Analysis Examples
This example uses the following Python files:
`oligopool_analysis_example.py
analysis_100k_reads_R1.fastq.gzip
analysis_100k_reads_R2.fastq.gzip
oligopool_analysis_example_barcodes.csv`

This example generates the following output files:
`oligopool_analysis_counts.csv`

For this example, we provide compressed .fastq files containing 100,000 paired-end sequencing reads. We will analyze these files and return counts of each barcode occurrence associated with each ID. This process has 3 steps: indexing, packing, and counting. 

In the indexing step, the function index.index() is called to read properly index the barcodes to their sequences IDs and save them to a .index file. This file is used for anchoring and mapping downstream in the counting step. In the packing step, pack.pack() processes the reads in the fastq.gzip files and packs them to be counted efficiently. In the counting step, xcount.xcount() anchors and maps reads to indexed sequences. Additional information on each of these is provided in Appendix A. 

In our script, in the index() function we define the parameters for barcodedata, barcodecol, barcodeprefix, barcodesuffix, indexfile. The barcodedata parameter takes our ‘oligopool_analysis_example_barcodes.csv’ file, which contains the information needed for indexing sequencing IDs, and anchoring and mapping reads to them later. In this .csv file, it is required that the sequence IDs are under a column header ‘ID’. The rest of the column headers may be appropriately labeled by the user. Note that you will need to have a labeled column for at least one of your barcodeprefix or barcodesuffix parameters. Either or both parameters may be defined. Our barcodes have known upstream and downstream sequences under the headers ‘Prefix’ and ‘Suffix’, which we specify in each of these parameters. The parameters barcodepregap and barcodepostgap may be defined as well, which specify the distance between the prefix and suffix sequences and the barcode. This can be useful if the sequences immediately upstream and downstream of your barcode contain variability, but here we define them as 0 since our constant prefix and suffix sequences are immediately adjacent to the barcodes. Finally, indexfile contains the filename for saving the output upon successful completion. Defining the filetype is not necessary here. The output of index() is a dictionary that contains information on the design process, such as completion status and warnings. 

Next, we call the pack.pack() to process our sequencing reads. In this example, we have paired-end reads stored analysis_100k_reads_R1.fastq.gzip and analysis_100k_reads_R2.fastq.gzip files. We have defined the parameters r1file, r1type, packtype, packfile, r1length, r1qual, r2file, r2type, r2length, r2qual, and packsize. Additional default settings can be viewed in Appendix A. The parameters r1file and r2file take in the file paths to their respective compressed .fastq files. r1type and r2type define the expected orientation of the reads in each file and take 0 or 1 integers as inputs, where 0 indicates forwards reads and 1 indicates reverse reads. Note that when counting later on, if the barcode sequence is not found, the xcount.xcount() function will attempt to find the barcode in the alternative orientation. r1length and r2length define the minimum acceptable length of the reads for the reads to be stored, which we set to 140. r1qual and r2qual define the minimum acceptable Phred-quality score for the reads to be stored, which we have set to 30 to ensure high confidence reads are stored. packtype defines whether to or not to join or merge reads, where 0 will store concatenated reads and 1 will merge the reads. Amplicons were about 400nt long, so the reads were not expected to overlap. Therefore, we defined this as 0, which will join reads. If overlapping reads are expected, setting packtype to 1 will merge reads. packsize defines the number of reads to be stored in each pack, in millions. Finally, packfile defines the filename to store the packfile. Since this process is computationally intensive, ncores and memlimit parameters will allow specifying the number of cores to run the process on, and the maximum amount of RAM (in GBs) per core. We have set these parameters to 0 for clarity in the script, but these are the default settings which will automatically determine how to allocate resources. Once again, this function returns a stats dictionary containing information about packing completion.

Now that we have indexed the barcodes and packed the sequencing reads, we can begin mapping and counting the barcode reads to the indexed sequences with xcount.xcount(). indexfiles, packfile, and countfile are the three parameters that are required to call xcount(). The former two simply contain the paths to their respective files that were generated by index() and pack(). The countfile parameter is the save filename where the counts for each associated ID are stored. The file will be saved in .csv format, specifying the file type is not necessary. We have also set the maptype to 1, which is slower but more accurate (for faster mapping, maptype may be set to 0). barcodeerrors specifies the maximum number of errors in a barcode to tolerate, which we’ve set to 3. Like index() and pack(), the output of xcount() is a dictionary containing information about counting completion.

Here, we’ve analyzed the data from a single sequencing reaction. When applying this pipeline to RT-qPCR, sequencing results for multiple reactions, such as DNA-seq and RNA-seq replicates can be consolidated into a single script and mapped to the same indexfile, to avoid repeating the indexing and packing steps. 

The filed saved by the script is ‘oligopool_analysis_counts.csv’ with two columns. The first column header is in the format ‘{your_indexfile}.ID’, which specifies the indexfile containing the IDs your reads were mapped to. The second column header is ‘Counts’, which contains the counts for each barcode associated with the respective ID, which can now be viewed.

## Appendix A: Package Functions 
A basic oligopool design pipeline will take place in the designparser.designparser() function using user-input design specifications. For analyzing barcode data, index.index() will index barcodes to IDs and generate .index files, pack.pack() will store .fastq files in packs according to user-input parameters, and xcount.xcount() will map and count barcodes in the pack files to the indexes and generate a .csv file containing ID’d barcode counts.

### function designparser
def designparser.designparser(pool_size, element_names, elements_spec, background_spec, padding_spec, split_spec=None)

Parameters
pool_size : int
The integer should be equal to the number of unique variants in your oligopool. If your variant sequences are stored in a list, you can simply use len(‘your_list’) to retrieve the pool size.

element_names : list
A list of strings for the elements in your oligo sequences, in the order you wish them to be arranged in each oligo sequence. The following example will name each element in your oligo sequence and specify their order, starting with ‘Primer1’ as the first element, and ending with ‘Padding’. 

Example: 
element_names = [‘Primer1’, ‘Cutsite2’, ‘Variant_Sequence’, ‘Cutsite2’,  ‘Barcode’, ‘Primer2’, ’Cutsite3’,  ‘Padding’]

elements_spec : dict
A nested dictionary of dictionaries. The keys of the element_spec dictionary should be the names in element_names. Each key in element_spec is also a nested dictionary, where each element’s type and specifications (which are unique to their type) are defined by the user. For a list of element types and specifications, see Appendix B.

Example:
elements_spec = {‘Primer1’: {‘type’: ‘Primer’, ‘oligolimit’: 150, …}, … ‘Padding’: {‘type’: ‘spacer’, … }}

background_spec : dict = {‘indata’ : list, ‘maxreplen’ : int }  / None
A dictionary containing the background sequences of your genetic system outside of your oligopool design. This ensures optimal oligopool design by preventing undesired sequences being introduced, for example, repetitive sequences. Including your sequence variants will ensure additional copies are not introduced in your designed sequences (for example, in padding sequence designs). ‘indata’ is a list of the background sequences. ‘maxreplen’ is the maximum repetitive sequence length allowed in the design. 

padding_spec : dict = {‘typeIIS’ : str / None, ‘oligolimit’ : int, ‘mintmelt’: int, ‘maxtmelt’: int, ‘maxreplen’: int}
A dictionary containing parameters for designing padding/filler sequences in your oligos with a specified TypeIIS recognition site, so the designed sequences are equal to the ‘oligolimit’. 'mintmelt’, ‘maxtmelt’ define the minimum and maximum melting temperatures of the design. ‘maxreplen’ defines the maximum repetitive sequence length allowed. 

split_spec : dict = {‘splitlimit’ “: int / None, ‘mintmelt’ : int, ‘minhdist’:  int,  ‘minoverlap’: int , ‘maxoverlap’: int } / None
split_spec defines parameters for splitting based on typeIIS recognition site given in padding_spec. ‘splitlimit’ defines maximum allowed oligolimit after splitting. ‘mintmelt’ defines the minimum melting temperature of the split regions. ‘minhdist’ defines the pairwise hamming distance between split regions at a given distance. ‘minoverlap’ and ‘maxoverlap’ define the minimum and maximum overlap allowed by the split. 

Output : dict = {‘output’: dict, ‘step’, ‘step_name’, ‘stats_dict’ }
designparser.designparser() returns a dictionary with information about the design. The ‘output’ key stores a dictionary with two dataframes. ‘annotated_complete’ stores your designed oligopool sequences indexed by their sequence ID, where each element type is a column (e.g. ‘Primer1’, ‘Cutsite1’, ‘Barcode’, etc.). ‘presplit_complete’ stores your designed full oligopool sequences under header ‘CompleteOligo’, indexed by their sequence ID. 
‘step’ stores the last step completed by designparser. 
‘step_name’ returns a description of the last step completed by designparser.
‘stats_dict’ stores nested dictionary with information about the design of the oligopool and each element, including ‘status’ and ‘basis’ describing success or failure, ‘step’ and ‘step_name’ store the last step completed and its name.

### function index
def index.index(barcodedata, barcodecol, barcodeprefix, barcodesuffix, indexfile, barcpde[rega[, barcpdepostgap, associatedata, associatecol, associateprefix, associatesuffix, verbose = True)

Parameters
barcodedata : str
A path to a .csv file storing barcode information. At minimum, column headers should include ‘ID’ for sequence IDs, and headers for barcodes, and barcode prefixes and/or suffixes. 
barcodecol : str
	Column name for the barcode sequences.

barcodeprefix : str
Column name for the constant upstream prefix sequence to the barcode. At least one of barcodeprefix or barcodesuffix must be provided.

barcodesuffix : str
Column name for the constant downstream suffix sequence to the barcode. At least one of barcodeprefix or barcodesuffix must be provided.

indexfile : str
Path to file for saving the index file. The suffix ‘oligopool.index’ will be automatically appended to the end of the file name.

barcodepregap : int
	Gap length between prefix sequence (if provided) and barcode. 

barcodpostgap : int
	Gap length between suffix sequence (if provided) and barcode. 

associatedata : str / None
	Filename for storing associate information. 

associatecol : str / None
	Column name for associated sequence.

associateprefix : str / None
	Column name for the constant upstream associate sequence prefix. 

associatesuffix : str / None
	Column name for the constant upstream associate sequence suffix.

verbose : bool
	If True, print processing statements. 

Output : dict = {‘output’: dict, ‘step’, ‘step_name’, ‘stats_dict’ }
A dictionary containing information on the indexing solved status. Keys of the dictionary are ‘status’ (True or False, depending on successful completion). A file of .index type is created which stores indexed reads to be used by xcount.xcount() for mapping and counting barcode sequences. 

### function pack
def pack.pack(r1file, r1type, packtype, packfile, r1length, r1qual, r2file, r2type, r2length, r2qual, packsize, ncores, memlimit, verbose)

r1file : str
	A path to the .fastq file containing R1 reads

 r1type : int
Defines orientation of reads. 0 specifies forward orientation. 1 specifies reverse orientation

packtype : int
A value of 0 will concatenate and store reads. A value of 1 will merge and store reads.

packfile : str
Path to filename to store packed reads. The suffix ‘oligopool.pack’ will be automatically appended to the end of the file name.

r1length : int
Minimum length of R1 reads to be stored. Reads shorter than defined length will not be stored. 

r1qual : int
	Minimum Phred-quality score for R1 reads to be stored. 

r2file : str
	A path to the .fastq file containing R2 reads

r2type : int 
Defines orientation of reads. 0 specifies forward orientation. 1 specifies reverse orientation

r2length : int
Minimum length of R2 reads to be stored. Reads shorter than defined length will not be stored. 

r2qual : int
Minimum Phred-quality score for R2 reads to be stored.

packsize : int
	Maximum number of reads (in millions) stored per pack.

ncores : int
	Total number of packers concurrently initiated. 

memlimit : nu.Real
Total amount of memory allowed per core. 

verbose : bool
	If True, print processing statements. 

xcount
def xcount.xcount(indexfiles, packfile, countfile, maptype, barcodeerrors, callback, ncores, memlimit, verbose)

indexfiles : str
	Path to .index file containing indexed barcodes generated by index.index().

packfile : str
	Path to .pack file containing packed reads generated by pack.pack().

countfile : str
Path to the filename to save counts found associated with each barcode ID. The suffix ‘oligopool.xcount.csv will be automatically append to the end of the file name.

maptype : int
A value of 0 will perform fast/approximate mapping. A value of 1 will perform slower, exact mapping. 

barcodeerrors : int
	Maximum number of sequencing errors allowed in barcode mapping to tolerate. 

callback : function(read, ID, count, coreid)
A custom callback function that may be passed for additional analysis. All callback functions must contain at least the parameters read, ID, count, and coreid.
read : str
	The read being counted/analyzed
ID : tuple
A tuple of IDs, one for each element mapped from each of the indexes provided. A value of ‘-‘ is given when an element was not mapped to any indices. 
count : int
	The associated read / ID frequency.
coreid : int
	The coreid integer identifier that processed this read.

ncores : int
	Total number of packers concurrently initiated. 

memlimit : nu.Real
Total amount of memory allowed per core. 
verbose : bool
	If True, print processing statements. 
 
## Appendix B: element_spec Types and Specifications
For each of the following types, included in the ‘type’ key of each element in the element_spec dictionary, the corresponding keywords must also be included in the dictionary unless specified as None.

### primer
‘oligolimit’ : int
		      Length of designed oligopool sequences. 

‘primerseq’ : str
       		An IUAC degenerate nucleotide primer sequence constraint

‘primertype’ : bool
          0 = forward primer design
          1 = reverse primer design

‘mintmelt’ : float
          Primer melting temperature lower bound

‘maxtmelt’ : float
          Primer melting temperature upper bound

‘maxreplen’ : int <= 6
          Maximum shared repeat length between the primers and flanking regions. 

‘pairedprimer’ : str
          The paired primer in the primer design. Ensures optimal design, including melting temperature matching and dimer minimization. 

‘leftcontext’ : str / None
          Upstream element. Checks design constraints are applied when elements are combined.

‘rightcontext’ : str / None
        Downstream element. Checks design constraints are applied when elements are combined.

‘exmotifs’ : list / None
	      Sequences to be excluded from element design. 

### ‘motif’ 
‘oligolimit’ : int
		    Length of designed oligopool sequences. 

‘motifseq’ : str
		    An IUAC degenerate nucleotide motif sequence. 

‘leftcontext’ : str / None
        Upstream element. Checks design constraints are applied when elements are combined.

‘rightcontext’ : str / None
        Downstream element. Checks design constraints are applied when elements are combined.

‘exmotifs’ : list / None
	      Sequences to be excluded from element design. 

### ‘variant’
	‘sequences’: list
        A list of all variants to be designed on your oligopool. The length of this list should inherently be equal to the pool_size parameter in desgingparser()

### ‘barcode’
‘oligolimit’ : int
		    Length of designed oligopool sequences

‘barcodelen’ : int
		    Length of designed barcodes

‘minhdist’ : int
		    Minimum pairwise hamming distance between barcode sequences

‘barcodetype’ : Bool int
        Indicates optimization parameter. ‘0’ for terminus optimized barcode. ‘1’ for spectrum optimized barcode

‘leftcontext’ : str / None
        Upstream element. Checks design constraints are applied when elements are combined

‘rightcontext’ : str / None
        Downstream element. Checks design constraints are applied when elements are combined

‘exmotifs’ : list / None
	      Sequences to be excluded from element design

### ‘spacer’
‘oligolimit’ : int
		    Length of designed oligopool sequences

‘spacerlen’ : int / None
        If int, generates a sequence according to design constraints of length int. If None, generates sequence according to design constraints for remaining sequence length up to ‘oligolimit’ considering lengths of all other elements.
 
‘leftcontext’ : str / None
        Upstream element. Checks design constraints are applied when elements are combined

‘rightcontext’ : str / None
        Downstream element. Checks design constraints are applied when elements are combined

‘exmotifs’ : list / None
	      Sequences to be excluded from element design

