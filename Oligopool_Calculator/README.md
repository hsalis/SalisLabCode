# Oligopool Calculator v1.0
### Code written by Ayaan Hossain and Howard Salis

## Relevant Article: Oligopool Calculator: automated design of oligopools and rapid analysis of massively parallel barcoded measurements
### Ayaan Hossain, Daniel P. Cetnar, Travis L. La Fleur, James McLellan, and Howard M. Salis

Oligopool synthesis and next-generation sequencing enable the construction and characterization of large libraries of designed genetic parts and systems. As library sizes grow, it becomes computationally challenging to optimally design large numbers of primer binding sites, barcode sequences, and overlap regions to obtain efficient assemblies and precise measurements. We present the Oligopool Calculator, an end-to-end suite of algorithms and data structures, that rapidly designs the many thousands of oligonucleotides within an oligopool and rapidly analyzes many billions of barcoded sequencing reads. We introduce several novel concepts that greatly increase the design and analysis throughput, including orthogonally symmetric barcode design, adaptive decision trees for primer design, a Scry barcode classifier, and efficient read packing. We demonstrate the Oligopool Calculator’s capabilities across computational benchmarks as well as real-data projects, illustrating how the flexible grammar is used to build and characterize thousands of genetic parts and systems, achieving design throughputs of XX primer binding sites and XX barcodes per hour and analysis throughputs of XX billion sequencing reads per hour. Overall, the Oligopool Calculator enables the creative design of massively parallel experiments, unburdening researchers of the computational complexities.

## Installation

Install Python 3, pip, and more:
`   apt-get update \
    && apt-get install -y python3 python3-pip build-essential zip wget pypy graphviz nano apt-transport-https ca-certificates curl software-properties-common intel-mkl git-core cmake pandoc memcached mpich mpich-doc libgfortran5 \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*`

Install the Python module dependencies:
`   pip3 install numpy pandas matplotlib seaborn scipy scikit-learn biopython dm-sonnet sonnet graphs mpi4py openpyxl leveldb==0.201 networkx nrpcalc statsmodels msgpack pyfastx edlib dinopy bounter numba primer3-py`

## Design Examples

In this example, we will used the Oligopool Calculator’s design functions to design an oligopool containing 1000 promoter sequence variants. The Oligopool Calculator will generate the needed primers, barcodes, and padding sequences using defined sequence constraints. 
Each oligonucleotide in the oligopool has the following architecture:

5’ -- [Forward Primer Site] [Cut1] [Promoter Variant] [Cut2] [spacer] [Cut3] [Barcode] [Cut4] [Reverse Primer Site] [Padding]-- 3’

where the Forward Primer Site and Reverse Primer Site are used to selectively amplify a subset of oligonucleotides, creating double-stranded DNA fragments;
the Cut1 and Cut4 restriction sites are used to clone the oligopool into a desired vector via digestion and ligation; 
the Promoter Variant sequence contains the designed promoter sequence to be constructed and characterized;
the Barcode is a uniquely identifying DNA sequence that is used to label each promoter sequence and quantitatively measure the DNA and mRNA composition within the library;
the Cut2 and Cut3 restriction sites are used to insert a DNA fragment into the middle of the oligopool, which we use to introduce a constant RBS-CDS sequence;
and the Padding is a padding sequence used to ensure that all oligonucleotides are the same size, which improves oligopool synthesis.

Provided are two files, ‘oligopool_design_examplel.py’ that demonstrates the workflow for designing the oligopool described above. The file ‘oligopool_design_variants.txt’ contains a list of 1000 example promoter variants.

The oligopool design pipeline takes place in the designpasrser.designparser() function (Appendix A). In this script, we will go through each item and explain their role in this context. The ‘result’ variable stores the series of nested dictionaries that comprise the designparser() output. Here, designparser() takes the parameters pool_size, element_names, element_spec, and background_spec.

pool_size defines the size of the oligopool. This will be the number of variants we have, so  we simply store the text file that contains the promoter sequences as a list, and use its length as the pool_size.  

element_names is a list that defines the names of each element in each sequence as 
dictionaries in the order of appearance. Note that in the provided script, this list matches
the sequence design above. The defined elements are ‘Primer1’, ‘Cut1’, ‘Promoter_Variant’, ‘Cut2’, ‘spacer’, ‘Cut3’, ‘Barcode’, ‘Cut4’, ‘Filler’. 

element_spec is a dictionary of element dictionaries where each key is a name from
element_names. For each element, there is a key ‘type’ that will dictate what other keys are
required within the element dictionary. See Appendix B for details on the valid keys 
compatible for each type in element_spec. For this example, we will go through the each element defined in the script, and note its function and important design decisions. The element types are: barcode, primer, motif, variant, and spacer. 

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


To ensure proper code execution, it is recommended you use the following python package versions:

| Package | Version |
| ------- | ------- |
| absl-py |  1.0.0  |
| anyio	4.0.0
| argon2-cffi	23.1.0
| argon2-cffi-bindings	21.2.0
| arrow	1.2.3
|asttokens	2.4.0
| astunparse	1.6.3
|async-lru	2.0.4
| attrs	23.1.0
Babel	2.12.1
backcall	0.2.0
beautifulsoup4	4.12.2
biopython	1.79
bleach	6.0.0
blinker	1.4
bounter	1.2.0
cachetools	5.0.0
certifi	2021.10.8
cffi	1.16.0rc2
charset-normalizer	2.0.12
cloudpickle	2.0.0
comm	0.1.4
cryptography	3.4.8
cycler	0.11.0
Cython	3.0.10
dbus-python	1.2.18
debugpy	1.8.0
decorator	5.1.1
defusedxml	0.7.1
dinopy	3.0.0
distro	1.7.0
distro-info	1.1build1
dm-sonnet	2.0.0
dm-tree	0.1.7
edlib	1.3.9
et-xmlfile	1.1.0
exceptiongroup	1.1.3
executing	1.2.0
fastjsonschema	2.18.0
flatbuffers	2
fonttools	4.32.0
fqdn	1.5.1
gast	0.5.3
google-auth	2.6.5
google-auth-oauthlib	0.4.6
google-pasta	0.2.0
graphs	0.1.3
grpcio	1.44.0
h5py	3.6.0
httplib2	0.20.2
idna	3.3
importlib-metadata	4.6.4
ipykernel	6.25.2
ipython	8.15.0
ipython-genutils	0.2.0
ipywidgets	8.1.1
isoduration	20.11.0
jedi	0.19.0
jeepney	0.7.1
Jinja2	3.1.2
joblib	1.1.0
json5	0.9.14
jsonpointer	2.4
jsonschema	4.19.1
jsonschema-specifications	2023.7.1
jupyter	1.0.0
jupyter_client	8.3.1
jupyter-console	6.6.3
jupyter_core	5.3.1
jupyter-events	0.7.0
jupyter-lsp	2.2.0
jupyter_server	2.7.3
jupyter_server_terminals	0.4.4
jupyterlab	4.1.0a1
jupyterlab-pygments	0.2.2
jupyterlab_server	2.25.0
jupyterlab-widgets	3.0.9
keras	2.8.0
Keras-Preprocessing	1.1.2
keyring	23.5.0
kiwisolver	1.4.2
launchpadlib	1.10.16
lazr.restfulclient	0.14.4
lazr.uri	1.0.6
leveldb	0.201
libclang	13.0.0
lightgbm	3.3.2
llvmlite	0.41.0
Markdown	3.3.6
MarkupSafe	2.1.3
matplotlib	3.5.1
matplotlib-inline	0.1.6
mistune	3.0.1
more-itertools	8.10.0
msgpack	1.0.6
nbclient	0.8.0
nbconvert	7.8.0
nbformat	5.9.2
nest-asyncio	1.5.8
networkx	3.1
notebook	7.0.4
notebook_shim	0.2.3
nrpcalc	1.6.3
numba	0.58.0
numpy	1.22.3
oauthlib	3.2.0
openpyxl	3.1.2
opt-einsum	3.3.0
overrides	7.4.0
packaging	21.3
pandas	1.4.2
pandocfilters	1.5.0
parasail	1.3.4
parso	0.8.3
patsy	0.5.3
pexpect	4.8.0
pickleshare	0.7.5
Pillow	9.1.0
pip	22.0.4
platformdirs	3.10.0
primer3-py	2.0.1
prometheus-client	0.17.1
prompt-toolkit	3.0.39
protobuf	3.20.0
psutil	5.9.5
ptyprocess	0.7.0
pure-eval	0.2.2
pyasn1	0.4.8
pyasn1-modules	0.2.8
pycparser	2.21
pyfastx	2.0.1
Pygments	2.16.1
PyGObject	3.42.0
PyJWT	2.3.0
pyparsing	2.4.7
python-apt	2.3.0+ubuntu2
python-dateutil	2.8.2
python-json-logger	2.0.7
pytz	2022.1
PyYAML	6.0.1
pyzmq	25.1.1
qtconsole	5.4.4
QtPy	2.4.0
referencing	0.30.2
requests	2.27.1
requests-oauthlib	1.3.1
rfc3339-validator	0.1.4
rfc3986-validator	0.1.1
rpds-py	0.10.3
rsa	4.8
scipy	1.8.0
seaborn	0.11.2
Send2Trash	1.8.0
setuptools	59.6.0
simplejson	3.18.1
six	1.16.0
sniffio	1.3.0
soupsieve	2.3.1
stack-data	0.6.2
statsmodels	0.13.2
tabulate	0.9.0

A suite of tools to design and analyze oligopool libraries for massively parallel assays.
