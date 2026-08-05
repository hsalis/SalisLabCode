import sys
import os
import gc
import shutil  as su
import time    as tt

import uuid    as uu
import pickle  as pk
import zipfile as zf
import atexit  as ae

import contextlib  as cl
import collections as cx

import ctypes          as ct
import multiprocessing as mp
import multiprocessing.sharedctypes as mpsct

import numpy   as np
import numba   as nb
import primer3 as p3
import psutil  as pu

import msgpack as mg
import pyfastx as pf
import edlib   as ed
import dinopy  as dp
import bounter as bt

import scry as sy


# Global Lookups

complement_table = str.maketrans(
    '-ACGTURYSWKMBVDHN',
    '-TGCATKYWSRMBDHVN')

dna_alpha = set(
    '-ATGC')

ddna_space = {
    '-': {'-'},
    'A': {'A'},
    'C': {'C'},
    'G': {'G'},
    'T': {'T'},
    'R': {'A', 'G'},
    'Y': {'C', 'T'},
    'S': {'G', 'C'},
    'W': {'A', 'T'},
    'K': {'G', 'T'},
    'M': {'A', 'C'},
    'B': {'C', 'G', 'T'},
    'V': {'A', 'C', 'G'},
    'D': {'A', 'G', 'T'},
    'H': {'A', 'C', 'T'},
    'N': {'A', 'T', 'G', 'C'}}

ddna_alpha = set(ddna_space.keys())

typeIIS_dict = {
    'acui'    : ('AcuI',     'CTGAAG', 16),
    'alwi'    : ('AlwI',     'GGATC',   5),
    'bbsi'    : ('BbsI',     'GAAGAC',  6),
    'bcci'    : ('BccI',     'CCATC',   5),
    'bceai'   : ('BceAI',    'ACGGC',  14),
    'bcivi'   : ('BciVI',    'GTATCC',  6),
    'bcodi'   : ('BcoDI',    'GTCTC',   5),
    'bmri'    : ('BmrI',     'ACTGGG',  5),
    'bpuei'   : ('BpuEI',    'CTTGAG', 16),
    'bsai'    : ('BsaI',     'GGTCTC',  5),
    'bseri'   : ('BseRI',    'GAGGAG', 10),
    'bsmai'   : ('BsmAI',    'GTCTC',   5),
    'bsmbi'   : ('BsmBI',    'CGTCTC',  5),
    'bsmfi'   : ('BsmFI',    'GGGAC',  14),
    'bsmi'    : ('BsmI',     'GAATGC',  1),
    'bspcni'  : ('BspCNI',   'CTCAG',   9),
    'bspqi'   : ('BspQI',    'GCTCTTC', 4),
    'bsrdi'   : ('BsrDI',    'GCAATG',  2),
    'bsri'    : ('BsrI',     'ACTGG',   1),
    'btgzi'   : ('BtgZI',    'GCGATG', 14),
    'btsci'   : ('BtsCI',    'GGATG',   2),
    'btsi'    : ('BtsI',     'GCAGTG',  2),
    'btsimuti': ('BtsIMutI', 'CAGTG',   2),
    'eari'    : ('EarI',     'CTCTTC',  4),
    'ecii'    : ('EciI',     'GGCGGA', 11),
    'esp3i'   : ('Esp3I',    'CGTCTC',  5),
    'faui'    : ('FauI',     'CCCGC',   6),
    'hgai'    : ('HgaI',     'GACGC',  10),
    'hphi'    : ('HphI',     'GGTGA',   8),
    'hpyav'   : ('HpyAV',    'CCTTC',   6),
    'mlyi'    : ('MlyI',     'GAGTC',   5),
    'mnli'    : ('MnlI',     'CCTC',    7),
    'sapi'    : ('SapI',     'GCTCTTC', 4),
    'sfani'   : ('SfaNI',    'GCATC',   9)}

trimflag = sys.platform.lower(
    ).startswith('linux')

# Shared Classes

class SafeCounter(object):
    '''
    Because mp.Value is a lie.
    Internal use only.
    '''

    def __init__(self, initval=0):
        '''
        Initialize SafeCounter instance with
        initval for storing integer counters.

        :: initval
           type - integer
           desc - initialization value
        '''
        self.counter = mpsct.RawValue(
            ct.c_longlong, initval)
        self.lock = mp.Lock()

    def increment(self, incr=1):
        '''
        Increment instance value by incr.

        :: incr
           type - integer
           desc - increment value
                  (default=1)
        '''
        with self.lock:
            self.counter.value += incr
        return self

    def decrement(self, decr=1):
        '''
        Decrease instance value by decr.

        :: decr
           type - integer
           desc - decrement value
                  (default=1)
        '''
        with self.lock:
            self.counter.value -= decr
        return self

    def value(self):
        '''
        Return the current value of instance.
        '''
        with self.lock:
            val = self.counter.value
        return val

class SafeQueue(object):
    '''
    Because mp.Queue is a lie.
    Internal use only.
    '''

    def __init__(
        self,
        maxsize=None):
        '''
        Initialize SafeQueue instance
        with locking for data sharing.
        '''
        if not maxsize is None:
            # Finite Queue
            self.queue = mp.Queue(
                maxsize=maxsize)
        else:
            # Infinite Queue
            self.queue = mp.Queue()

        # Setup Sentinels
        self.counter = SafeCounter(initval=0)
        self.lock    = mp.Lock()
        self.alive   = mpsct.RawValue(
            ct.c_bool, True)

    def alivemethod(func):
        '''
        Decorator to ensure SafeQueue instance
        is alive. Raise Error otherwise.
        '''
        def wrapper(self, *args, **kwargs):
            if self.alive.value is True:
                return func(self, *args, **kwargs)
            else:
                raise RuntimeError(
                    'SafeQueue instance is closed')
        return wrapper

    @alivemethod
    def put(self, item):
        '''
        Put item in queue.

        :: item
           type - PyObject
           desc - object to be placed
        '''
        with self.lock:
            self.queue.put(item)
            self.counter.increment()
        return self

    @alivemethod
    def multiput(self, iterable):
        '''
        Put multiple items in queue.

        :: iterable
           type - iterable
           desc - iterable of items
        '''
        with self.lock:
            itemcount = 0
            for item in iterable:
                self.queue.put(item)
                itemcount += 1
            self.counter.increment(
                incr=itemcount)

    @alivemethod
    def get(self):
        '''
        Remove and return an item
        from the queue.
        '''
        with self.lock:
            if self.counter.value() > 0:
                item = self.queue.get()
                self.counter.decrement()
            else:
                item = None
        return item

    @alivemethod
    def multiget(self, numitems=None):
        '''
        Get multiple items in queue.

        :: numitems
           type - integer
           desc - total number of items
                  to get from queue
                  (default=None)
        '''
        with self.lock:
            if numitems is None:
                numitems = float('inf')
            itemcount = 0
            threshold = min(
                self.counter.value(),
                numitems)
            while itemcount < threshold:
                yield self.queue.get()
                itemcount += 1
            self.counter.decrement(
                decr=itemcount)

    def length(self):
        '''
        Return the number of items
        in queue.
        '''
        with self.lock:
            length = self.counter.value()
        return length

    @alivemethod
    def __len__(self):
        '''
        Dunder function for length.
        '''
        return self.length()

    @alivemethod
    def empty(self):
        '''
        Return True if the queue is empty,
        False otherwise.
        '''
        with self.lock:
            status = self.queue.empty()
        return status

    @alivemethod
    def full(self):
        '''
        Return True if the queue is full,
        False otherwise.
        '''
        with self.lock:
            status = self.queue.full()
        return status

    @alivemethod
    def close(self):
        '''
        Close the queue and release internal
        resources.
        '''
        for x in self.multiget():
            pass
        with self.lock:
            self.queue.close()
            self.queue.join_thread()
            self.queue   = None
            self.counter = None
            self.alive.value = False

# Decorators and Contexts

def coroutine(func):
    '''
    Decorator used to prime a coroutine.
    Internal use only.

    :: func
       type - coroutine
       desc - a coroutine we want online
    '''

    def online(*args, **kwargs):
        _coroutine = func(*args, **kwargs)
        next(_coroutine)
        return _coroutine

    return online

@cl.contextmanager
def ignored(*exceptions):
    '''
    Raymond Hettinger's ignored
    context. Internal use only.
    '''

    try:
        yield
    except exceptions:
        pass

# Printing / Logging

def removestarfix(string, fix, loc):
    '''
    Remove prefix or suffix from string.
    Because Manifold Bio doesn't believe
    in Python 3.9 :'). Internal use only.

    :: string
       type - string
       desc - some string of interest
    :: fix
       type - string
       desc - some prefix or suffix of
              interest to be trimmed
              from given string
    :: loc
       type - integer
       desc - prefix / suffix identifier
              <= 0 <-- remove prefix
              >  0 <-- remove suffix

    '''

    # Remove prefix: str.removeprefix(...)
    if loc <= 0:
        if string.startswith(fix):
            return string[len(fix):]
    # Remove suffix: str.removesuffix(...)
    elif loc > 0:
        if string.endswith(fix):
            return string[:-len(fix)]
    # Default case
    return string

@coroutine
def liner_engine(online=True):
    '''
    Return a coroutine to print stuff
    to stdout and log run information.
    Internal use only.

    :: online
       type - boolean
       desc - if True will print received
              strings to stdout
              (default=True)
    '''

    # Book-keeping
    clrlen = 0

    # Online
    try:
        while True:

            # Receive String
            printstr = (yield)

            # Line Preservation
            if online and printstr.startswith('\*'):
                sys.stdout.write('\n')
                clrlen   = 0
                printstr = removestarfix(
                    string=printstr,
                    fix='\*',
                    loc=-1)

            # Line String
            if online:
                sys.stdout.write('\r' + ' '*clrlen)
                sys.stdout.write('\r' + printstr)
                clrlen = len(printstr)
                sys.stdout.flush()

    # Closure
    except GeneratorExit:
        if online:
            sys.stdout.write('\n')

def get_printlen(value):
    '''
    Return printing length.
    Internal use only.

    :: value
       type - nu.Real
       desc - value to evaluate
    '''

    svalue = str(value)
    if len(svalue) <= 3:
        return len(svalue)
    else:
        return len(str(value)) + \
            int(safelog10(value) / 3)

def get_notelen(printlen):
    '''
    Return notation length.
    Internal use only.

    :: value
       type - integer
       desc - printlen to evaluate
    '''

    if printlen <= 15:
        return 'd', printlen
    else:
        return 'e', 0

# Numeric Functions

def get_trials(prob):
    '''
    Return the number of required trials
    given a probability of success.
    Internal use only.

    :: prob
       type - float
       desc - probability of a rate event
    '''
    growth_rate = -100.0 # Lower is Rarer
    # print(-10 / np.log(1. - prob))
    return 2 * int(np.ceil(
        growth_rate / np.log(1. - prob)))

def get_prob(success, trials):
    '''
    Return the probability of success given
    the number of successes in some finite
    trials. Internal use only.

    :: success
       type - integer
       desc - total number of successes
    :: trials
       type - integer
       desc - total number of attempts made
    '''
    return (success + 1.) / (trials + 2.)

def safediv(A, B):
    '''
    Return A / B, except when B is zero return
    zero (no more pesky Div by Zero Errors).
    Internal use only.

    :: A
       type - nu.Real
       desc - numerator
    :: B
       type - nu.Real
       desc - denominator
    '''

    return 0. if B == 0. else float(A) / B

def safelog10(A):
    '''
    Return the log10 of A, except when A is less
    than or equal to zero return zero (no more
    pesky Math Domain Errors).
    Internal use only.

    :: A
       type - nu.Real
       desc - quantity
    '''

    return np.log10(A) if A > 0. else 0.

def safelog(A, n=10):
    '''
    Return the log10 of A, except when A is less
    than or equal to zero return zero (no more
    pesky Math Domain Errors).
    Internal use only.

    :: A
       type - nu.Real
       desc - quantity
    :: n
       type - nu.Real
       desc - base (default=10)
    '''

    return np.log(A) / np.log(n) if (A > 0. and n > 0.) else 0.

def get_sample(value, lf, uf):
    '''
    Return a random integer between
    (lf*value, uf*value).
    Internal use only.

    :: value
       type - nu.Real
       desc - numeric value
    :: lf
       type - nu.Real
       desc - sample lower bound fraction
              of value
    :: uf
       type - nu.Real
       desc - sample upper bound fraction
              of value
    '''

    if value <= 10:
        return round(value / 2)
    return np.random.randint(
        *map(np.round, (value * lf,
                        value * uf)))

# DataFrame Functions

def get_unique_count(iterable):
    '''
    Return the count of unique elements
    in iterable. Internal use only.

    :: iterable
       type - iterable
       desc - iterable to uniquify
    '''

    seen = set()
    for element in iterable:
        if not element in seen:
            seen.add(element)
    return len(seen)

def get_uniques(
    iterable,
    typer):
    '''
    Return the unique elements in iterable
    packed via typer. Internal use only.

    :: iterable
       type - iterable
       desc - iterable to uniquify
    :: typer
       type - function
       desc - factory function to wrap
              uniques extracted
    '''

    uniques = []
    seen    = set()
    for element in iterable:
        if not element in seen:
            uniques.append(element)
            seen.add(element)
    return typer(uniques)

def get_col_exist_idx(
    col,
    df):
    '''
    Determine if col exists as a column in df,
    and return its index. Internal use only.

    :: col
       type - string
       desc - column name to check for existence
    ::df
       type - pd.DataFrame
       desc - DataFrame to check existence in
    '''

    col = col
    collist = tuple(df.columns)
    if col == df.index.name or \
       col in collist:
        return True, collist.index(col)
    else:
        return False, None

def get_df_concat(df):
    '''
    Concatenate all columns in df.
    Internal use only.

    :: df
       type - pd.DataFrame
       desc - DataFrame to concatenate
              columns from
    '''

    return tuple(df.sum(
        axis=1).str.replace(
            '-', '').values)

def update_df(
    indf,
    lcname,
    rcname,
    out,
    outcol):
    '''
    Insert out values in outcol between
    lcname and rcname cols if provided.
    Internal use only.

    :: indf
       type - pd.DataFrame
       desc - DataFrame to update
    :: lcname
       type - string / None
       desc - left context column name
    :: rcname
       type - string
       desc - right context column name
    :: out
       type - list / string
       desc - value(s) to insert
    :: outcol
       type - string
       desc - column name for values
    '''

    # Book-keeping
    insidx    = len(indf.columns)
    insstatus = False

    # Do we have a left context?
    if not insstatus and \
       not lcname is None:

        # Determine existence and index
        (lcexists,
        lcidx) = get_col_exist_idx(
            col=lcname,
            df=indf)

        # Update Book-keeping
        if lcexists:
            insidx    = lcidx + 1
            insstatus = True

    # Do we have a right context?
    if not insstatus and \
       not rcname is None:

        # Determine existence and index
        (rcexists,
        rcidx) = get_col_exist_idx(
            col=rcname,
            df=indf)

        # Update Book-keeping
        if rcexists:
            insidx    = rcidx
            insstatus = True

    # Update indf
    indf.insert(
        loc=insidx,
        column=outcol,
        value=out)

# Oligo Functions

@nb.njit
def get_store_hdist(
    store,
    idx,
    i=0,
    j=None,
    t=None,
    direction=0):
    '''
    Return the minimum pairwise hamming
    distance between store[idx, i:j] and
    store[*: i:j] depending on direction.
    Internal use only.

    :: store
       type - np.array
       desc - numeric sequence array
    :: idx
       type - integer
       desc - location index of sequence
              being compared
    :: i
       type - integer
       desc - starting col index of comparison
    :: j
       type - integer
       desc - ending   col index of comparison
    :: t
       type - np.array / None
       desc - fancy row index for taking
    :: direction
       type - integer
       desc - direction flag;
              0 = upward   comparison
              1 = downward comparison
              2 = all-pair comparison
    '''

    # Default Result
    hdist = store.shape[1]
    if i >= 0 and not j is None:
        hdist = j-i

    # Specific Coordinates
    if not t is None:
        return min(
            hdist,
            (store[t, i:j] != store[idx, i:j]).sum(1).min())

    # Something to compare against?
    if idx > 0:

        # Upward / All-Pair Comparison
        if direction == 0 or direction == 2:
            hdist = min(
                hdist,
                (store[:idx, i:j] != store[idx, i:j]).sum(1).min())

    # Something to compare against?
    if idx < store.shape[0] - 1:

        # Downward / All-Pair Comparison
        if direction == 1 or direction == 2:
            hdist = min(
                hdist,
                (store[idx+1:, i:j] != store[idx, i:j]).sum(1).min())

    # Return Result
    return hdist

@nb.njit
def get_pair_hdist(
    store,
    idx1,
    idx2,
    i=0,
    j=None):
    '''
    Return the pairwise hamming distance between
    store[idx1, i:j] and store[idx2, i:j].
    Internal use only.

    :: store
       type - np.array
       desc - numeric sequence array
    :: idx1
       type - integer
       desc - location index of the first
              sequence being compared
    :: idx2
       type - integer
       desc - location index of the second
              sequence being compared
    :: i
       type - integer
       desc - starting index of comparison
    :: j
       type - integer
       desc - ending   index of comparison
    '''

    return (store[idx1, i:j] != store[idx2, i:j]).sum()

def get_tmelt(
    seq,
    i=0,
    j=None,
    mvc=50.,
    dvc=0.,
    ntc=0.8,
    olc=50.):
    '''
    Return the melting temperature of seq.
    Internal use only.

    Note: primer3-py API docs

    :: seq
       type - string
       desc - DNA string in context
    :: i
       type - integer
       desc - starting index of comparison
              (default=)
    :: j
       type - integer
       desc - ending index of comparison
    :: mvc
       type - float
       desc - monovalent cation conc. (mM)
              (default=50.0 nM)
    :: dvc
       type - float
       desc - divalent cation conc. (mM)
              (default=0.0 nM)
    :: ntc
       type - float
       desc - dNTP conc. (mM)
              (default=0.8 mM)
    :: olc
       type - float
       desc - oligo conc. (nM)
              (default=50.0 nM)
    '''
    return p3.bindings.calcTm(
        seq=seq[i:j],
        mv_conc=mvc,
        dv_conc=dvc,
        dntp_conc=ntc,
        dna_conc=olc)

def get_constant_regions(seqconstr):
    '''
    Extract all constant, non-degenerate
    regions embedded in seqconstr.
    Internal use only.

    :: seqconstr
       type - string
       desc - element sequence constraint
    '''

    # Make a local copy
    cseqconstr = str(seqconstr)

    # Remove dgenerate nucleotides
    for nt in set(seqconstr) - set('ATGC'):
        cseqconstr = cseqconstr.replace(nt, '-')

    # Split and extract defined regions
    regions = set(cseqconstr.split('-'))
    if '' in regions:
        regions.remove('')

    # Return results
    return regions

# Oligopool Functions

def get_variantlens(indf):
    '''
    Return the length of all variants
    in the input DataFrame. Internal
    use only.

    :: indf
       type - pd.DataFrame
       desc - input DataFrame containing
              all designed variants
    '''

    return np.array(list(map(
        len, get_df_concat(
            df=indf))))

def get_parsed_oligolimit(
    indf,
    variantlens,
    oligolimit,
    minelementlen,
    maxelementlen,
    element,
    liner):
    '''
    Check if enough space available
    for all variants in pool for the
    designed element.

    :: indf
       type - pd.DataFrame
       desc - input DataFrame containing
              all designed variants
    :: variantlens
       type - np.array / None
       desc - length of all variants stored
              in input DataFrame, if None
              then compute it explicitly
    :: oligolimit
       type - integer
       desc - maximum allowed oligo length
    :: maxelementlen
       type - integer
       desc - maximum allowed designed element
              length
    :: minelementlen
       type - integer
       desc - maximum allowed designed element
              length
    :: element
       type - string
       desc - element being designed
    :: liner
       type - coroutine
       desc - dynamic printing
    '''

    # Book-keeping
    parsestatus   = True
    minvariantlen = None
    maxvariantlen = None
    minspaceavail = None
    maxspaceavail = None
    t0 = tt.time()

    # Compute Variant Lengths
    liner.send(' Parsing Variant Lengths ...')

    if variantlens is None:
        variantlens = get_variantlens(indf=indf)

    minvariantlen = np.min(variantlens)
    maxvariantlen = np.max(variantlens)
    minspaceavail = oligolimit - maxvariantlen
    maxspaceavail = oligolimit - minvariantlen

    parsestatus = (0 <= minelementlen <= minspaceavail) and \
                  (0 <= maxelementlen <= maxspaceavail)

    # Show Updates
    plen = get_printlen(
        value=max(np.abs([
            oligolimit,
            minvariantlen,
            maxvariantlen,
            minelementlen,
            maxelementlen,
            minspaceavail,
            maxspaceavail]))) + 1

    liner.send(
        ' Maximum Oligo Length: {:{},d} Base Pair(s)\n'.format(
            oligolimit,
            plen))

    # How much space occupied by variants?
    if minvariantlen == maxvariantlen:
        liner.send(
            ' Input Variant Length: {:{},d} Base Pair(s)\n'.format(
                minvariantlen,
                plen))
    else:
        liner.send(
            ' Input Variant Length: {:{},d} to {:,} Base Pair(s)\n'.format(
                minvariantlen,
                plen,
                maxvariantlen))

    # How much space required by elements?
    if minelementlen == maxelementlen:
        liner.send(
            ' {:>13} Length: {:{},d} Base Pair(s)\n'.format(
                'Reqd. ' + element,
                minelementlen,
                plen))
    else:
        liner.send(
            ' {:>13} Length: {:{},d} to {:,} Base Pair(s)\n'.format(
                'Reqd. ' + element,
                minelementlen,
                plen,
                maxelementlen))

    # How much space available?
    if not parsestatus:
        msg = ' [INFEASIBLE] (Insufficient Oligo Space Available)'
    else:
        msg = ''

    if minspaceavail == maxspaceavail:
        liner.send(
            ' Free Space Available: {:{},d} Base Pair(s){}\n'.format(
                minspaceavail,
                plen,
                msg))
    else:
        liner.send(
            ' Free Space Available: {:{},d} to {:,} Base Pair(s){}\n'.format(
                minspaceavail,
                plen,
                maxspaceavail,
                msg))

    # Show Time Elapsed
    liner.send(
        ' Time Elapsed: {:.2f} sec\n'.format(
            tt.time()-t0))

    # Show Verdict
    if not parsestatus:
        liner.send(
            ' Verdict: {} Design Infeasible due to Oligo Limit Constraints\n'.format(
                element))
    else:
        liner.send(
            ' Verdict: {} Design Possibly Feasible\n'.format(
                element))

    # Return Results
    return (parsestatus,
        minvariantlen,
        maxvariantlen,
        minelementlen,
        maxelementlen,
        minspaceavail,
        maxspaceavail)

def get_parsed_oligopool_repeats(
    df,
    maxreplen,
    element,
    merge,
    liner):
    '''
    Check if oligopool repeats are
    feasible. Internal use only.

    :: df
       type - pd.DataFrame
       desc - DataFrame to extract
              oligopool repeats from
    :: maxreplen
       type - integer
       desc - maximum shared repeat
              length with oligopool
    :: element
       type - string
       desc - element being designed
    :: liner
       type - coroutine
       desc - dyanamic printing
    '''

    # Book-keeping
    if not merge:
        oligorepeats = {}
    else:
        oligorepeats = bt.bounter(
            need_iteration=False, # Static Checks Only
            size_mb=4096)         # 4 GB Spectrum Cap
    t0           = tt.time()
    kmerspace    = ((4**(maxreplen+1)) // 2)
    fillcount    = None
    freecount    = None
    repeatcount  = 0

    # Verbage Stuff
    plen = get_printlen(
        value=kmerspace)
    if plen > 15:
        sntn = 'e'
        plen = len('{:e}'.format(kmerspace))
    else:
        sntn = 'd'

    # Show Update
    if merge:
        sourcecontext = False
        liner.send(' Repeat Source: From Entire Oligopool\n')
    else:
        sourcecontext = True
        liner.send(' Repeat Source: Contextually per Variant\n')

    liner.send('  k-mer Space : {:{},{}} Unique {:,}-mers\n'.format(
        kmerspace,
        plen,
        sntn,
        maxreplen+1))

    # Extract Repeats
    liner.send(' Extracting Repeats ...')

    for idx,oligo in enumerate(
        get_df_concat(df=df)):

        # Show Update
        if (idx+1) % 1000 == 0:
            liner.send(
                ' Extracting Repeats: Processed {:,} Sequences'.format(
                    idx+1))

        # Store Repeats
        if not merge:
            transformation = set(map(
                lambda x: dp.conversion.encode_twobit(
                    seq=x,
                    sentinel=False),
                stream_canon_spectrum(
                    seq=oligo,
                    k=maxreplen+1)))
            oligorepeats[idx] = tuple(
                sorted(transformation))
            repeatcount = max(
                repeatcount,
                len(oligorepeats[idx]))
        else:
            transformation = stream_canon_spectrum(
                seq=oligo,
                k=maxreplen+1)
            oligorepeats.update(
                transformation)
    if merge:
        repeatcount = oligorepeats.cardinality()

    # Compute Feasibility
    fillcount   = min(kmerspace, repeatcount)
    freecount   = kmerspace - fillcount
    parsestatus = (freecount * 1.) / kmerspace > .01
    if parsestatus:
        statusmsg = ''
    else:
        statusmsg = ' [INFEASIBLE]'

    # Show Final Updates
    liner.send('   Fill Count : {:{},{}} Unique {:,}-mers ({:6.2f} %)\n'.format(
        fillcount,
        plen,
        sntn,
        maxreplen+1,
        safediv(
            A=fillcount*100.,
            B=kmerspace)))
    liner.send('   Free Count : {:{},{}} Unique {:,}-mers ({:6.2f} %){}\n'.format(
        freecount,
        plen,
        sntn,
        maxreplen+1,
        safediv(
            A=freecount*100.,
            B=kmerspace),
        statusmsg))

    if not parsestatus:
        liner.send(
            ' Verdict: {} Design Infeasible due to Repeat Length Constraint\n'.format(
                element))
    else:
        liner.send(
            ' Verdict: {} Design Possibly Feasible\n'.format(
                element))

    liner.send(
        ' Time Elapsed: {:.2f} sec\n'.format(
            tt.time()-t0))

    # Return Results
    return (parsestatus,
        sourcecontext,
        kmerspace,
        fillcount,
        freecount,
        oligorepeats)

# Exmotif Functions

def stream_exmotif_splits(exmotif):
    '''
    Stream all prefix-suffix splits
    for given motif. This is also
    left partition of motif.
    Internal use only.

    :: motif
       type - string
       desc - motif to split
    '''

    return ((exmotif[:i], exmotif[i:]) \
        for i in range(1, len(exmotif)))

def get_exmotif_partition(exmotifs):
    '''
    Compute the left partition of
    exmotifs. Internal use only.

    :: exmotifs
       type - cx.deque
       desc - deque of all motifs
              to be excluded
    '''

    partition = cx.defaultdict(list)
    for exmotif in exmotifs:
        for u,v in stream_exmotif_splits(
            exmotif=exmotif):
            partition[u].append(v)
    return partition

def get_inverted_exmotif_partition(partition):
    '''
    Compute the right partition of
    left partitioned exmotifs by
    inverting left partition.
    Internal use only.

    :: partition
       type - dict
       desc - dictionary of left
              partitioned exmotifs,
              prefix -> suffix
    '''

    inv_partition = cx.defaultdict(list)
    for u,v in partition.items():
        for vi in v:
            inv_partition[vi].append(u)
    return inv_partition

def is_right_blocked(suffixes):
    '''
    Determine if all suffixes are blocked
    for a given left partition.
    Internal use only.

    :: suffixes
       type - list
       desc - left partition motif suffixes
              that share a common prefix
    '''

    return len(set(vi[0] for vi in suffixes)) == 4

def is_left_blocked(prefixes):
    '''
    Determine if all suffixes are blocked
    for a given left partition.
    Internal use only.

    :: prefixes
       type - list
       desc - right partition motif prefixes
              that share a common suffix
    '''

    return len(set(vi[-1] for vi in prefixes)) == 4

def get_parsed_exmotifs(
    exmotifs,
    typer,
    element,
    leftcontext,
    rightcontext,
    warn,
    liner):
    '''
    Check feasibility and return all
    exmotifs sorted by length.
    Internal use only.

    :: exmotifs
       type - list
       desc - list of all motifs
              to be excluded
    :: typer
       type - function
       desc - factory function for
              wrapping exmotifs
    :: element
       type - string
       desc - element being designed
    :: leftcontext
       type - tuple
       desc - tuple of all left context
              sequences
    :: rightcontext
       type - tuple
       desc - tuple of all right context
              sequences
    :: warn
       type - dict
       desc - warning dictionary entry
    :: liner
       type - coroutine
       desc - dynamic printing
    '''

    # Time-keeping
    t0 = tt.time()

    # Sort Enque all motifs by length
    liner.send(' Sorting and Enqueing Excluded Motif(s) ...')

    exmotifs = sorted(
        exmotifs,
        key=len)
    dq = typer(exmotifs)

    liner.send(' Sorted and Enqued: {:,} Unique Excluded Motif(s)\n'.format(
        len(exmotifs)))

    # Check motif feasibility
    liner.send(' Computing Excluded Motif Length Distribution ...')

    # Compute length distribution
    cr = cx.Counter(len(m) for m in dq)

    liner.send(' Excluded Motif Length Distribution\n')

    klen = get_printlen(
        value=max(cr.keys()))
    vlen = get_printlen(
        value=max(cr.values()))

    # Check for infeasible lengths
    parsestatus = True
    problens    = []

    for mlen in sorted(cr.keys()):

        # Compute infeasibility
        parsemsg = ''
        if cr[mlen] == 4**(mlen):
            parsemsg    = ' [INFEASIBLE] (All {}-mers Excluded)'.format(mlen)
            parsestatus = False
            problens.append(mlen)

        # Show update
        liner.send('   - {:{},d} Motif(s) of Length {:{},d}{}\n'.format(
            cr[mlen],
            vlen,
            mlen,
            klen,
            parsemsg))

    # Finalize problens
    problens = None if len(problens) == 0 else problens

    # Partition Exmotifs
    if (not  leftcontext is None) or \
       (not rightcontext is None):

        # Show Update
        liner.send(' Left Partitioning Excluded Motif(s) ...')

        # Right Partition exmotifs
        leftpartition = get_exmotif_partition(
            exmotifs=exmotifs)

        # Show Update
        liner.send(' Right Partitioning Excluded Motif(s) ...')

        # Right Partition exmotifs
        rightpartition = get_inverted_exmotif_partition(
            partition=leftpartition)

        # Setup Warning Variables
        warn['vars'] = {
            'exmotifprefixgroup': set(),
            'exmotifsuffixgroup': set()}

        # Detect Potentially Blocked Left Exmotifs
        if not leftcontext is None:
            liner.send(' Detecting Blocked Excluded Motif Prefix(es) ...')
            for prefix,suffixes in leftpartition.items():
                if is_right_blocked(suffixes=suffixes):
                    warn['vars']['exmotifprefixgroup'].add(prefix)

        # Detect Potentially Blocked Right Exmotifs
        if not rightcontext is None:
            liner.send(' Detecting Blocked Excluded Motif Suffix(es) ...')
            for suffix,prefixes in rightpartition.items():
                if is_left_blocked(prefixes=prefixes):
                    warn['vars']['exmotifsuffixgroup'].add(suffix)

        # Update Warning Counts
        warn['warncount'] += len(warn['vars']['exmotifprefixgroup'])
        warn['warncount'] += len(warn['vars']['exmotifsuffixgroup'])

        # Show Blocked Motifs
        if warn['warncount']:
            liner.send(
                ' Found {:,} Problematic Excluded Motif Group(s)\n'.format(
                    warn['warncount']))
            plen = max(map(
                len,
                warn['vars']['exmotifprefixgroup'] | \
                warn['vars']['exmotifsuffixgroup'])) + 2 + 1

            # Show Prefix Updates
            for prefix in sorted(warn['vars']['exmotifprefixgroup'], key=len):
                prefix = '\'' + prefix + '.'*(plen-len(prefix)-2) + '\''
                liner.send(
                    '   - Motif(s) starting with {:<{}} [WARNING] (Prefix Prevents All 4 Bases After It)\n'.format(
                        prefix,
                        plen))

            # Show Suffix Updates
            for suffix in sorted(warn['vars']['exmotifsuffixgroup'], key=len):
                suffix = '\'' + '.'*(plen-len(suffix)-2) + suffix + '\''
                liner.send(
                    '   - Motif(s)   ending with {:>{}} [WARNING] (Suffix Prevents All 4 Bases Before It)\n'.format(
                        suffix,
                        plen))

            # Finalize Left Partition Warning
            if not warn['vars']['exmotifprefixgroup']:
                warn['vars']['exmotifprefixgroup'] = None

            # Finalize Left Partition Warning
            if not warn['vars']['exmotifsuffixgroup']:
                warn['vars']['exmotifsuffixgroup'] = None

    # No Partition Required
    else:
        leftpartition  = None
        rightpartition = None

    # Show feasibility verdict
    liner.send(
        ' Time Elapsed: {:.2f} sec\n'.format(
            tt.time()-t0))

    if not parsestatus:
        liner.send(
            ' Verdict: {} Design Infeasible due to Excluded Motif Constraints\n'.format(
                element))
    else:

        # Warnings Incurred
        if warn['warncount']:
            liner.send(
                ' Verdict: {} Design Potentially Infeasible\n'.format(
                    element))

        # No Warnings
        else:
            liner.send(
                ' Verdict: {} Design Possibly Feasible\n'.format(
                    element))

    # Return Results
    return (parsestatus,
        dq,
        problens,
        leftpartition,
        rightpartition)

def get_exmotif_conflict(
    seq,
    seqlen,
    exmotifs,
    partial=False,
    checkall=False):
    '''
    Determine if the sequence does not contain or is
    contained in one of the excluded motifs (motif
    feasibility). Internal use only.

    :: seq
       type - string
       desc - sequence for conflict checking
    :: seqlen
       type - integer
       desc - length of sequence
    :: exmotifs
       type - iterable / None
       desc - list of motifs to exclude
              in designed barcodes,
              None otherwise
    :: partial
       type - boolean
       desc - if True will not check for
              conflicts for motifs longer
              than seqlen, otherwise run
              all checks
              (default=False)
    :: checkall
       type - boolean
       desc - if True will check for all
              possible conflicts due to
              exmotifs
              (default=False)
    '''

    # Book-keeping
    status = True # No Conflict
    pmotif = set() if checkall else None

    # Do we have anything to exclude?
    if exmotifs:

        # Loop through all motifs
        for motif in exmotifs:

            # Is sequence partial?
            if partial:
                if len(motif) > seqlen:
                    break # No need to check further

            # Embedding and Embedded Conflict
            if (len(motif) <= seqlen and motif in seq) or \
               (len(motif)  > seqlen and seq   in motif):

                # We got Conflict!
                status = False

                # What is the Conflict?
                if checkall:
                    pmotif.add(motif)
                else:
                    pmotif = motif
                    break

    # Return result
    if status:
        return (True, None)
    else:
        return (False, pmotif)

def get_exmotif_conflict_index(
    seq,
    conflicts):
    '''
    Compute the location of all
    motif conflicts in seq.
    Internal use only.

    :: seq
       type - string
       desc - sequence for conflict
              checking
    :: conflicts
       type - iterable
       desc - iterable of
    '''

    # Book-keeping
    index = cx.defaultdict(set)

    # Search Loop
    for motif in conflicts:
        start = 0 # Initial Start Index
        while True:

            # Locate Motif
            idx = seq.find(motif, start)

            # None Found
            if idx == -1:
                break

            # Found Instance
            else:
                index[motif].add(idx)
                start = idx+1

    # Return Results
    return index

def is_local_exmotif_feasible(
    seq,
    exmotifs,
    exmotifindex):
    '''
    Determine if sequence devoid of exmotifs.
    Internal use only.

    :: seq
       type - string
       desc - a partially explored sequence path
    :: exmotifs
       type - set / None
       desc - set of all excluded motifs
    :: exmotifindex
       type - set / None
       desc - set of constraint embedded
              exmotif ending indices
    '''

    # Ignore Current Path?
    if not exmotifindex is None:
        if len(seq) in exmotifindex:
            return True, None

    # Exmotifs to Prevent
    if not exmotifs is None:

        # Verify Absence
        for mlen in exmotifs:

            # Path is Partial
            if len(seq) < mlen:
                # No Additional Checks Needed
                break

            # Check Trace
            trace = seq[-mlen:]
            if trace in exmotifs[mlen]:
                # Found Conflict
                return False, trace

    # No Conflict
    return True, None

# Context Functions

def get_edgeeffectlength(exmotifs):
    '''
    Get the context length to search in for exmotif
    edge effects. Internal use only.

    :: exmotifs
       type - iterable / None
       desc - iterable of motifs to exlude,
              None otherwise
    '''

    if not exmotifs is None:
        return len(exmotifs[-1])
    return 0

def get_grouped_sequences(sequences):
    '''
    Group all selected strings by length.
    Internal use only.

    :: sequences
       type - iterable
       desc - set of all sequences / motifs
              to group by length
    '''

    groupdict = cx.OrderedDict()
    for sequence in sorted(sequences, key=len):
        if not len(sequence) in groupdict:
            groupdict[len(sequence)] = set()
        groupdict[len(sequence)].add(sequence)
    return groupdict

def get_extracted_edge(
    contextseq,
    position,
    edgeeffectlength):
    '''
    Extract context sequence edge based on
    context position. Internal use only.

    :: contextseq
       type - string
       desc - context sequence to evaluate
    :: position
       type - integer
       desc - context position identifier
              0 =  left context
              1 = right context
    :: edgeeffectlength
       type - integer
       desc - length of context sequence to
              extract for edge-effect eval
    '''
    # Do we have edges to consider?
    if not edgeeffectlength is None:

        # Extract Left Context
        if position == 0:
            return contextseq[-edgeeffectlength:]

        # Extract Right Context
        else:
            return contextseq[:+edgeeffectlength]

    # No edges to consider
    return None

def get_extracted_context(
    leftcontext,
    rightcontext,
    edgeeffectlength,
    reduce,
    liner):
    '''
    Return parsed context sequences.
    Internal use only.

    :: leftcontext
       type - pd.Series
       desc - left context sequences
    :: rightcontext
       type - pd.Series
       desc - right context sequences
    :: edgeeffectlength
       type - integer / None
       desc - length of context sequence to
              extract for edge-effect eval
    :: reduce
       type - boolean
       desc - if True will reduce multi-unique
              context sequences to uniques;
              otherwise return full context
    :: liner
       type - coroutine
       desc - dynamic printing
    '''

    # Start Timer
    t0 = tt.time()

    # Setup Contexts
    contexts         = ((leftcontext, 0),  (rightcontext, 1))
    uniquecount      = []
    extractedcontext = []

    # Parse Context
    for context,position in contexts:

        # No Context
        if context is None:
            uniquecount.append(0)
            extractedcontext.append(None)

        # Context Present
        else:

            # Book-keeping
            uniques  = set()
            extracts = list()

            # Analyze Context
            for idx,seq in enumerate(context):

                # Extract Edge Sequence
                edgeseq = get_extracted_edge(
                    contextseq=seq,
                    position=position,
                    edgeeffectlength=edgeeffectlength)

                # Show Update
                if (idx+1) % 1000 == 0:
                    liner.send(
                        ' Extracting {} Context: Processed {:,} Sequences'.format(
                            ['Left', 'Right'][position],
                            idx+1))

                # Store Edge Sequence
                uniques.add(edgeseq)
                extracts.append(edgeseq)

            # Update uniquecount
            uniquecount.append(len(uniques))

            # Context is Constant
            if len(uniques) == 1:
                extractedcontext.append(
                    [extracts[-1]])

            # Context is Variable
            else:
                # Only Uniques
                if reduce:
                    extractedcontext.append(
                        list(uniques))
                else:
                    extractedcontext.append(
                        extracts)

    # Compute plen
    plen = get_printlen(
        value=max(uniquecount))

    # Show Updates
    liner.send(
        '  Left Context: {:{},d} Unique Sequence(s)\n'.format(
            uniquecount[0],
            plen))
    liner.send(
        ' Right Context: {:{},d} Unique Sequence(s)\n'.format(
            uniquecount[1],
            plen))
    liner.send('  Time Elapsed: {:.2f} sec\n'.format(tt.time()-t0))

    # Return Parsed Context
    return (extractedcontext[0],
        extractedcontext[1])

def get_context_type_selector(context):
    '''
    Return selector functions for building
    element context. Internal use only.

    :: context
       type - list / None
       desc - list of context sequence,
              or None
    '''

    # Multiple contexts
    if isinstance(context, list):

        # Compute uniquecount
        uniquecount = len(set(context))

        # Variable context
        if uniquecount > 1:
            return (1, lambda x: context[x])

        # Constant context
        else:
            return (2, lambda x: context[0])

    # No context
    elif context is None:
        return (3, lambda x: None)

    # Unknown packing
    else:
        raise TypeError(
            'Invalid Context Type {}'.format(
                type(context)))

def get_parsed_edgeeffects(
    sequence,
    element,
    leftcontext,
    rightcontext,
    leftpartition,
    rightpartition,
    exmotifs,
    merge,
    warn,
    liner):
    '''
    Cluster left and right context sequences
    and record forbidden prefix and suffixes.
    Internal use only.

    :: sequence
       type - string
       desc - element sequence constraint
    :: element
       type - string
       desc - designed element name
    :: leftcontext
       type - tuple
       desc - tuple of all left context
              sequences
    :: rightcontext
       type - tuple
       desc - tuple of all right context
              sequences
    :: leftpartition
       type - dict / None
       desc - left partition of excluded
              motifs if leftcontext present
    :: rightpartition
       type - dict / None
       desc - right partition of excluded
              motifs if rightcontext present
    :: exmotifs
       type - cx.deque
       desc - deque of all motifs
              to be excluded
    :: merge
       type - boolean
       desc - if True, merge all entries in
              prefixdict and suffixdict
    :: warn
       type - dict
       desc - warning dictionary entry
    :: liner
       type - coroutine
       desc - dynamic printing
    '''

    # Book-keeping
    prefixdict  = None
    suffixdict  = None
    lcwarncount = set()
    rcwarncount = set()
    constantprefix = None
    constantsuffix = None

    # Time-keeping
    t0 = tt.time()

    # Define plen and qlen
    plen = get_printlen(
        value=max(
            len(leftcontext)  if not  leftcontext is None else 0,
            len(rightcontext) if not rightcontext is None else 0))
    qlen = len(exmotifs[-1]) + 3

    # Assess Left Edge Effects
    if leftcontext:

        # Show Updates
        liner.send(' Parsing Left Context ...')

        # Define Prefix Dict
        prefixdict = cx.defaultdict(lambda: cx.defaultdict(set))

        # Compute Prefix Conficts
        for lcseq in leftcontext:

            # Have we seen this prefix already?
            if lcseq in prefixdict:
                continue

            # Crunch Constraints
            for prefix in leftpartition:

                # Does lcseq impose constraint?
                if lcseq.endswith(prefix):

                    # print((lcseq, prefix, leftpartition[prefix]))

                    # Identify Forbidden Motif Prefixes
                    prefixforbidden = set(leftpartition[prefix])

                    # Compute Tolerated Prefixes due to Motif Sequence
                    toleratedprefix = set()
                    for fp in prefixforbidden:
                        # RIP .. Have to Tolerate
                        if sequence.startswith(fp):
                            toleratedprefix.add(fp)
                            # Record Longest Constant Prefix
                            if constantprefix is None:
                                constantprefix = fp
                            else:
                                if len(constantprefix) < len(fp):
                                    constantprefix = fp
                    if toleratedprefix:
                        lcwarncount.add(lcseq)

                    # Update Forbidden Motif Prefixes
                    prefixforbidden -= toleratedprefix
                    # prefixdict[lcseq].update(prefixforbidden)
                    for fp in prefixforbidden:
                        prefixdict[lcseq][fp].add(prefix+fp)

            # Reduce Forbidden Prefixes
            fpx = sorted(prefixdict[lcseq].keys(), key=len)
            # idxo = index outer/shorter prefix
            # idxi = index inner/longer  prefix
            for idxo in range(len(fpx)-1):
                for idxi in range(idxo+1, len(fpx)):
                    if  fpx[idxi] in prefixdict[lcseq] and \
                        fpx[idxi].startswith(fpx[idxo]): # fpo is a prefix for fpx[idxi]
                        prefixdict[lcseq].pop(fpx[idxi])

            # Freeze Forbidden Prefixes
            for fp in prefixdict[lcseq]:
                prefixdict[lcseq][fp] = tuple(prefixdict[lcseq][fp])

            # Group Forbidden Prefixes
            prefixdict[lcseq]['keys'] = get_grouped_sequences(
                sequences=prefixdict[lcseq])

            # Show Update
            liner.send(
                '  Left Context {:{},d}: {:.>{}} Prevents {:,} Prefix(es)'.format(
                    len(prefixdict),
                    plen,
                    lcseq,
                    qlen,
                    len(prefixdict[lcseq])))

        # Retain Last Update
        liner.send('\*')

    # Assess Right Edge Effects
    if rightcontext:

        # Show Updates
        liner.send(' Parsing Right Context ...')

        # Define Suffix Dict
        suffixdict = cx.defaultdict(lambda: cx.defaultdict(set))

        # Compute Suffix Conficts
        for rcseq in rightcontext:

            # Have we seen this Suffix already?
            if rcseq in suffixdict:
                continue

            # Crunch Constraints
            for suffix in rightpartition:

                # Does rcseq impose constraint?
                if rcseq.startswith(suffix):

                    # print((rcseq, suffix, rightpartition[suffix]))

                    # Identify Forbidden Motif Suffixes
                    suffixforbidden = set(rightpartition[suffix])

                    # Compute Tolerated Suffixes due to Motif Sequence
                    toleratedsuffix = set()
                    for fs in suffixforbidden:
                        # RIP .. Have to Tolerate
                        if sequence.endswith(fs):
                            toleratedsuffix.add(fs)
                            # Record Longest Constant Suffix
                            if constantsuffix is None:
                                constantsuffix = fs
                            else:
                                if len(constantsuffix) < len(fs):
                                    constantsuffix = fs
                    if toleratedsuffix:
                        rcwarncount.add(rcseq)

                    # Update Forbidden Motif Suffixes
                    suffixforbidden -= toleratedsuffix
                    # suffixdict[rcseq].update(suffixforbidden)
                    for fs in suffixforbidden:
                        suffixdict[rcseq][fs].add(fs+suffix)

            # Reduce Forbidden Suffixes
            fsx = sorted(suffixdict[rcseq].keys(), key=len)
            # idxo = index outer/shorter suffix
            # idxi = index inner/longer  suffix
            for idxo in range(len(fsx)-1):
                for idxi in range(idxo+1, len(fsx)):
                    if  fsx[idxi] in suffixdict[rcseq] and \
                        fsx[idxi].endswith(fsx[idxo]): # fso is a suffix for fsx[idxi]
                        suffixdict[rcseq].pop(fsx[idxi])

            # Freeze Forbidden Suffixes
            for fs in suffixdict[rcseq]:
                suffixdict[rcseq][fs] = tuple(suffixdict[rcseq][fs])

            # Group Forbidden Suffixes
            suffixdict[rcseq]['keys'] = get_grouped_sequences(
                sequences=suffixdict[rcseq])

            # Show Update
            liner.send(
                ' Right Context {:{},d}: {:.<{}} Prevents {:,} Suffix(es)'.format(
                    len(suffixdict),
                    plen,
                    rcseq,
                    qlen,
                    len(suffixdict[rcseq])))

        # Retain Last Update
        liner.send('\*')

    # Merge Dictionaries
    if merge:

        # Show Update
        liner.send(' Merging Forbidden Prefixes ...')

        # Merge prefixdict
        if not prefixdict is None:

            # Define Merged Dictionary
            _prefixdict = {}

            # Merging Loop
            for d in prefixdict.values():
                for k in d:
                    if not k in _prefixdict:

                        # Special Keys Arrangement
                        if k == 'keys':
                            _prefixdict[k] = {}
                        # Regular Arrangement
                        else:
                            _prefixdict[k] = set()

                    # Special Keys Arrangement
                    if k == 'keys':
                        for l in d[k]:
                            if not l in _prefixdict[k]:
                                _prefixdict[k][l] = set()
                            else:
                                _prefixdict[k][l].update(d[k][l])
                    # Regular Arrangement
                    else:
                        _prefixdict[k].update(d[k])

            # Restore Order
            prefixdict = {}
            for k in _prefixdict:
                if k == 'keys':
                    prefixdict[k] = cx.OrderedDict()
                    for l in sorted(_prefixdict[k]):
                        prefixdict[k][l] = set(_prefixdict[k][l])
                else:
                    prefixdict[k] = tuple(
                        sorted(_prefixdict[k], key=len))

        # Show Update
        liner.send(' Merging Forbidden Suffixes ...')

        # Merge suffixdict
        if not suffixdict is None:

            # Define Merged Dictionary
            _suffixdict = {}

            # Merging Loop
            for d in suffixdict.values():
                for k in d:
                    if not k in _suffixdict:

                        # Special Keys Arrangement
                        if k == 'keys':
                            _suffixdict[k] = {}
                        # Regular Arrangement
                        else:
                            _suffixdict[k] = set()

                    # Special Keys Arrangement
                    if k == 'keys':
                        for l in d[k]:
                            if not l in _suffixdict[k]:
                                _suffixdict[k][l] = set()
                            else:
                                _suffixdict[k][l].update(d[k][l])
                    # Regular Arrangement
                    else:
                        _suffixdict[k].update(d[k])

            # Restore Order
            suffixdict = {}
            for k in _suffixdict:
                if k == 'keys':
                    suffixdict[k] = cx.OrderedDict()
                    for l in sorted(_suffixdict[k]):
                        suffixdict[k][l] = set(_suffixdict[k][l])
                else:
                    suffixdict[k] = tuple(
                        sorted(_suffixdict[k], key=len))

    # Show Final Updates
    if len(lcwarncount) or len(rcwarncount):

        # Compute plen
        maxlen = 0
        if lcwarncount:
            maxlen = max(maxlen, len(constantprefix))
        if rcwarncount:
            maxlen = max(maxlen, len(constantsuffix))
        plen = get_printlen(value=maxlen) + 2 + 3

        # Show Major Category
        liner.send(' Found Terminal Constant Base(s) in Sequence Constraint\n')

        # Show Responsible Constants
        if constantprefix:
            cp = '\'' + constantprefix + '.'*(plen-len(constantprefix)-2) + '\''
            liner.send(
                '   - {} Design starts with {:<{}} [WARNING] (Constant Prefix)\n'.format(
                    element, cp, maxlen))

        if constantsuffix:
            cs = '\'' + '.'*(plen-len(constantsuffix)-2) + constantsuffix + '\''
            liner.send(
                '   - {} Design   ends with {:>{}} [WARNING] (Constant Suffix)\n'.format(
                    element, cs, maxlen))

        # Optimization Status
        suboptimal = True

        # Compute plen
        plen = get_printlen(
            value=max(
                len(lcwarncount),
                len(rcwarncount)))

        # Show Major Category
        liner.send(
            ' Found {:,} Resulting Inexorable Edge-Effects\n'.format(
                len(lcwarncount) + len(rcwarncount)))

        # Left Context
        if len(lcwarncount):
            liner.send(
                '   - {:{},d} Unique  Left Context(s) Impacted [WARNING]\n'.format(
                    len(lcwarncount),
                    plen))

        # Right Context
        if len(rcwarncount):
            liner.send(
                '   - {:{},d} Unique Right Context(s) Impacted [WARNING]\n'.format(
                    len(rcwarncount),
                    plen))

    else:
        # Optimization Status
        suboptimal = False

    # Show Time Elapsed
    liner.send(
        ' Time Elapsed: {:.2f} sec\n'.format(
            tt.time()-t0))

    # Show Verdict
    if suboptimal:
        liner.send(
            ' Verdict: {} Design Potentially With Edge Effects\n'.format(
                element))
    else:
        liner.send(
            ' Verdict: {} Design Possibly Feasible\n'.format(
                element))

    # Update Warning
    warn['warncount'] = len(lcwarncount) + \
                        len(rcwarncount) + \
                        (1 if not constantprefix is None else 0) + \
                        (1 if not constantsuffix is None else 0)
    warn['vars'] = {
        'constantprefix'    : constantprefix,
        'constantsuffix'    : constantsuffix,
        'leftcontextimpact' : len(lcwarncount),
        'rightcontextimpact': len(rcwarncount)}

    # Return Results
    return (prefixdict, suffixdict)

def is_local_edge_feasible(
    seq,
    seqlen,
    lcseq,
    rcseq,
    edgeeffectlength,
    prefixforbidden,
    suffixforbidden):
    '''
    Determine if sequence prefix and suffix
    is forbidden. Internal use only.

    :: seq
       type - string
       desc - a paritally explored sequence path
    :: lcseq
       type - string / None
       desc - left context sequence
    :: rcseq
       type - string / None
       desc - right context sequence
    :: seqlen
       type - integer
       desc - full sequence length
    :: edgeeffectlength
       type - integer
       desc - context length for edge effects
    :: prefixforbidden
       type - dict / None
       desc - dictionary of forbidden primer
              prefix sequences
    :: suffixforbidden
       type - dict / None
       desc - dictionary of forbidden primer
              suffix sequences
    '''

    # Sequence Prefix Selected Forbidden?
    if not prefixforbidden is None:

        # In Exploration
        if len(seq) in prefixforbidden['keys']:

            # Prefix Conflict Found
            if seq in prefixforbidden['keys'][len(seq)]:

                # Detect Implicated Exmotifs
                dxmotifs = prefixforbidden[seq]

                # Return Results
                return False, dxmotifs, len(seq) - 1

        # Post Exploration
        if len(seq) == seqlen:

            # Build In Context Sequence
            if (not rcseq is None) and \
               (edgeeffectlength > seqlen):
                diff = edgeeffectlength - seqlen
                incontext = seq + rcseq[:diff]
            else:
                incontext = seq

            # Check for Conflicts
            for mlen in prefixforbidden['keys']:

                # Too Small within Context?
                if mlen > len(incontext):
                    break

                # Extract Prefix
                seqprefix = incontext[:mlen]

                # Prefix Conflict Found
                if seqprefix in prefixforbidden['keys'][mlen]:

                    # Detect Implicated Exmotifs
                    dxmotifs = prefixforbidden[seqprefix]

                    # Return Results
                    return False, dxmotifs, min(len(seq)-1, mlen-1)

    # Sequence Prefix Selected Forbidden?
    if not suffixforbidden is None:

        # Post Exploration
        if len(seq) == seqlen:

            # Build In Context Sequence
            if (not lcseq is None) and \
               (edgeeffectlength > seqlen):
                diff = edgeeffectlength - seqlen
                incontext = lcseq[-diff:] + seq
            else:
                incontext = seq

            # Check for Conflicts
            for mlen in suffixforbidden['keys']:

                # Too Small within Context?
                if mlen > len(incontext):
                    break

                # Extract Suffix
                seqsuffix = incontext[-mlen:]
                # mlen = len(seq)

                # Suffix Conflict Found
                if seqsuffix in suffixforbidden['keys'][mlen]:

                    # Detect Implicated Exmotifs
                    dxmotifs = suffixforbidden[seqsuffix]

                    # Return Results
                    return False, dxmotifs, len(seq) - 1

    # Everything is OK
    return True, None, None

# Sequence Analysis Functions

def is_DNA(seq, dna_alpha=dna_alpha):
    '''
    Determine if seq is a DNA string.
    Internal use only.

    :: seq
       type - string
       desc - a candidate sequence
    :: dna_alpha
       type - set
       desc - alphabet set considered
              valid for seq
              (defaul={A, T, G, C})
    '''

    if isinstance(seq, str) and \
       set(seq.upper()) <= dna_alpha:
        return True
    return False

def get_comp(seq):
    '''
    Return the complement of seq.
    Internal use only.

    :: seq
       type - string
       desc - a string in alphabet
              {A, T, G, C}
    '''

    return seq.translate(complement_table)

def get_revcomp(seq):
    '''
    Return the reverse complement of seq.
    Internal use only.

    :: seq
       type - string
       desc - a string in alphabet
              {A, T, G, C}
    '''

    return get_comp(seq)[::-1]

def stream_spectrum(seq, k):
    '''
    Stream the k-mer spectrum of seq.
    Internal use only.

    :: seq
       type - string
       desc - a string in alphabet
              {A, T, G, C}
    :: k
       type - integer
       desc - k-mer length / k-value
    '''

    return (seq[i:i+k] for i in range(len(seq)-k+1))

def stream_canon_spectrum(seq, k):
    '''
    Stream the canonical k-mer spectrum
    of seq. Internal use only.

    :: seq
       type - string
       desc - a string in alphabet
              {A, T, G, C}
    :: k
       type - integer
       desc - k-mer length / k-value
    '''

    return map(
        lambda x: min(x, get_revcomp(x)),
        stream_spectrum(seq=seq, k=k))

def stream_contigs(seq, scheme):
    '''
    Stream the k-mer contigs of seq.
    Internal use only.

    :: seq
       type - string
       desc - a string in alphabet
              {A, T, G, C}
    :: scheme
       type - integer
       desc - contig break scheme
    '''

    i = 0
    j = 0
    while i < len(seq):
        k = scheme[j]
        yield (seq[i:i+k],j)
        i += k
        j += 1

def get_hdist(seq1, seq2, max_hd=None):
    '''
    Compute the hamming distance between
    seq1 and seq2, within max_hd bound.
    Internal use only.

    :: seq1
       type - string
       desc - a string in alphabet
              {A, T, G, C}
    :: seq2
       type - string
       desc - a string in alphabet
              {A, T, G, C}
    :: max_hd
       type - integer / None
       desc - compute hamming distance
              between seq1 and seq2,
              within max_hd bound,
              return None otherwise
              (default=None)
    '''

    seq_len = min(len(seq1), len(seq2))
    i = 0
    hdist = 0
    while i < seq_len:
        hdist += seq1[i] != seq2[i]
        i += 1
        if not max_hd is None and \
           hdist > max_hd:
            return None
    return hdist

def get_edist(seq1, seq2, max_ed=None, mode='NW'):
    '''
    Compute the mode edit distance between
    seq1 and seq2, within max_ed bound.
    Internal use only.

    :: seq1
       type - string
       desc - a string in alphabet
              {A, T, G, C}
    :: seq2
       type - string
       desc - a string in alphabet
              {A, T, G, C}
    :: max_ed
       type - integer / None
       desc - compute edit distance
              between seq1 and seq2,
              within max_ed bound,
              return None otherwise
    :: mode
       type - string
       desc - 'NW' / 'HW' for global
              / infix alignments
    '''

    if len(seq2) > len(seq2):
        seq1, seq2 = seq2, seq1
    i = 0
    if max_ed is None:
        max_ed = -1
    ed_align = ed.align(
        seq2,
        seq1,
        mode=mode,
        task='distance',
        k=max_ed)
    edist = ed_align['editDistance']
    if edist == -1:
        return None
    return edist

def get_tvalue(elementlen):
    '''
    Return the mismatch tolerance
    value based on elementlen.
    Internal use only.

    :: elementlen
       type - integer
       desc - length of element
    '''

    if   elementlen <=  8:
        return 0
    elif elementlen <= 15:
        return 1
    return ((elementlen - 1) // 5) - 1

# System and IO Functions

def stream_fastq_engine(
    filepath,
    invert=False,
    qualvec=False,
    skipcount=0,
    coreid=0,
    ncores=1,
    fastqid=None,
    liner=None):
    '''
    Stream reads, and associated quality
    information stored in a FastQ file.
    Internal use only.

    :: filepath
       type - string
       desc - FastQ file storing reads
    :: invert
       type - boolean
       desc - if True will reverse complement
              read, and reverse the quality
              vector
              (default=False)
    :: qualvec
       type - boolean
       desc - if True will compute and return
              the quality vector instead of
              just the quality string
    :: skipcount
       type - integer
       desc - total number of reads from
              the beginning that may be
              (default=0)
    :: coreid
       type - integer
       desc - current core integer id
              (default=0)
    :: ncores
       type - integer
       desc - total number of readers
              concurrently initiated
              (default=1)
    :: liner
       type - coroutine / None
       desc - dynamic printing
    '''

    # Try reading file
    try:

        # Setup reading
        entry = iter(pf.Fastq(
            file_name=filepath,
            build_index=False,
            full_index=False))

        # Skip Book-keeping
        clen     = len(str(ncores))
        scount   = 0
        screach  = 0
        sctarget = round(
            skipcount * np.random.uniform(
                0.01, 0.02, 1)[0])
        skiphase = scount < skipcount

        # Parsing Book-keeping
        rcount = 0
        tcount = coreid

        # Scan reads
        while True:

            # Initial skipping
            try:
                # Skip initial entries
                if skiphase:

                    # Fetch Read
                    next(entry)

                    # Update Trackers
                    scount  += 1
                    screach += 1
                    skiphase = scount < skipcount

                    # Time to show Updates?
                    if screach == sctarget:
                        if not liner is None:
                            # Show Updates
                            liner.send(' Core {:{},d}: Reached R{} - {:6.2f} % '.format(
                                coreid,
                                clen,
                                fastqid,
                                (scount * 100.) / skipcount))

                        # Reset Tracker
                        screach = 0

                    # Continue Skipping ..
                    continue

            # EOF handling
            except:
                # We're done!
                break

            # Skip/parse normal entries
            try:
                # Skip entires
                if rcount < tcount:
                    next(entry)
                    rcount += 1
                    continue

                # Parse and yield (read, qual)
                fqe  = next(entry)
                read = fqe[1]
                qual = fqe[2]

                # Do we correct Orientation?
                if invert:
                    read = get_revcomp(
                        seq=read)
                    qual = qual[::-1]

                # Vectorize Quality
                if qualvec:
                    qual = np.frombuffer(
                        qual.encode(),
                        dtype=np.int8) - 33

                # Return Results
                yield (read, qual)

                # Update counters
                rcount += 1
                tcount += ncores

            # EOF handling
            except Exception as E:
                # We're done!
                break

    # Handle exceptions
    except Exception as E:
        raise E

    # End-of-file Token
    yield None, None

def needs_restart(memlimit):
    '''
    Determine if process needs be
    restart given memory profile.
    Internal use only.

    :: memlimit
       type - float
       desc - maximum memory allowed
              for process execution
              before restarting
    '''

    # Check Memory Consumed So Far (in GBs)
    memused = pu.Process(
        os.getpid()).memory_info().rss / 10**9
    # Check Memory Available (in GBs)
    memleft = pu.virtual_memory().available / 10**9
    # Compute Policy Status
    policystatus = any((
        memused >= memlimit,
        memleft <= 2,))
    # Return Results
    return policystatus

# Workspace Functions

def get_adjusted_path(path, suffix):
    '''
    Adjust outfile with suffix.
    Internal use only.

    :: path
       type - string
       desc - some path which must
              have given suffix
    :: suffix
       type - string / None
       desc - suffix to concatenate
              with given path
    '''

    # Do we have a valid string path?
    if not isinstance(path, str):
        return None

    # Adjust path
    path = path.strip()
    path = removestarfix(
        string=path,
        fix='/',
        loc=1)
    if not suffix is None and \
       not path.endswith(suffix):
        path += str(suffix)

    # Return adjusted path
    return path

def get_path_status(
    path,
    suffix=None,
    readable=True,
    writable=False,
    creatable=False):
    '''
    Arbitrary path evaluation.
    Internal use only.

    :: path
       type - string
       desc - path to store or output
              information
    :: suffix
       type - string / None
       desc - required path suffix
              (default=None)
    :: readable
       type - boolean
       desc - if True check is path is
              readable  (file/dir   read)
              (default=True)
    :: writable
       type - boolean
       desc - if True check if path is
              writable  (file/dir  write)
              (default=False)
    :: creatable
       type - boolean
       desc - if True check if path is
              creatable (file/dir create)
              (default=False)
    '''

    # Status:
    # 0  = Path is invalid      input
    # 1  = Path is unreadable   file
    # 2  = Path is unwritable   file
    # 3  = Path is empty        file
    # 4  = Path is non-empty    file
    # 5  = Path is unreadable   directory
    # 6  = Path is unwritable   directory
    # 7  = Path is empty        directory
    # 8  = Path is non-empty    directory
    # 9  = Path is non-existent
    # 10 = Path is creatable
    # 11 = Path is uncreatable
    # X  = Path is existent     special file

    # Is the path valid?
    if not isinstance(path, str) or not path:
        return 0 # only non-empty strings are paths

    # Adjust Path for Checks
    path = get_adjusted_path(
        path=path,
        suffix=suffix)

    # Is the path in existence?
    if os.path.exists(path):

        # Path points to a file?
        if os.path.isfile(path):

            # Read accessible?
            if  readable and not \
                os.access(path, os.R_OK):
                    return 1 # Unreadable file

            # Is the file empty?
            if not os.stat(path).st_size:

                # Can we write to this file?
                if  writable and not \
                    os.access(path, os.W_OK):
                        return 2 # Unwritable file

                return 3 # Empty file

            else:
                return 4 # Non-empty file

        # Path points to a directory?
        elif os.path.isdir(path):

            # What's the accessibility like?
            if  readable and not \
                os.access(path, os.R_OK):
                    return 5 # Unreadable directory

            # Directory empty?
            if not any(os.scandir(path)):

                # Can we write to this directory?
                if  writable and not \
                    os.access(path, os.W_OK):
                        return 6 # Unwritable directory

                return 7 # Empty Directory

            else:
                return 8 # Non-empty Directory

        # Path points to a special file such as
        # socket, FIFO, device file etc.
        else:
            return 'X' # Special File Identifier

    # Is the path creatable?
    if creatable:

        try:
            # Build parent directory path
            parentpath = os.path.dirname(path) or \
                         os.getcwd()
            # Try building parent path
            os.makedirs(
                name=parentpath,
                exist_ok=True)

        except:
            return 11 # Cannot create path

        else:
            return 10 # Can create path

    # Path is non-existent
    return 9

def setup_directory(dirpath):
    '''
    Build dirpath. Internal use only.

    :: dirpath
       type - string
       desc - a directory built to
              store information
    '''

    dirpath = get_adjusted_path(
        path=dirpath,
        suffix=None)
    if not os.path.exists(dirpath):
        os.makedirs(dirpath)

def setup_workspace(
    outfile,
    outfile_suffix):
    '''
    Adjust given output file and setup
    output directory, and schedule its
    removal at exit. Internal use only.

    :: outfile
       type - string
       desc - output file
    :: outfile_suffix
       type - string
       desc - output file suffix
    '''

    # Adjust outfile with suffix
    outfile = get_adjusted_path(
        path=outfile,
        suffix=outfile_suffix)

    # Define outdir
    outdir = '.'.join([
        outfile,
        str(uu.uuid4()),
        'dir'])

    # Schedule outdir deletion
    ae.register(
        remove_directory,
        outdir)

    # Setup outdir
    setup_directory(
        dirpath=outdir)

    # Return Workspace
    return outfile, outdir

def remove_file(filepath):
    '''
    Delete file stored at filepath.
    Internal use only.

    :: filepath
       type - string
       desc - a file to be removed
    '''

    if not filepath is None:
        filepath = get_adjusted_path(
            path=filepath,
            suffix=None)
        with ignored(OSError):
            os.remove(filepath)

def remove_directory(dirpath):
    '''
    Delete dirpath and its contents.
    Internal use only.

    :: dirpath
       type - string
       desc - a directory (with content)
              to be removed
    '''

    if not dirpath is None:
        dirpath = get_adjusted_path(
            path=dirpath,
            suffix=None)
        with ignored(OSError):
            su.rmtree(dirpath)

# Data Saving / Loading Functions

def msgpacksave(obj, filepath):
    '''
    MSGPACK dump an object to filepath.
    Internal use only.

    :: obj
       type - object
       desc - object to persist to disk
    :: filepath
       type - string
       desc - filepath to save object
    '''

    # Pack Data
    outcontent = mg.packb(
        o=obj,
        use_single_float=False,
        autoreset=True,
        use_bin_type=True,
        strict_types=False)

    # Check Previous File Non-Existent
    if os.path.exists(filepath):
        raise RuntimeError(
            '{} exists!'.format(filepath))

    # Write Pack to File
    with open(filepath, 'wb') as outfile:
        outfile.write(outcontent)

    # Cleanup
    del outcontent

def savedict(dobj, filepath):
    '''
    MSGPACK dump a dictionary to filepath.
    Internal use only.

    :: dobj
       type - dict
       desc - dictionary to persist to disk
    :: filepath
       type - string
       desc - filepath to save dictionary
    '''

    return msgpacksave(
        obj=dobj,
        filepath=filepath)

def sorteditems(dobj):
    '''
    Return a list of dicionary items
    sorted by keys. Internal use only.

    :: dobj
       type - dict
       desc - a dictionary to sort by
              keys
    '''

    return sorted(dobj.items())

def sortedcounts(cobj):
    '''
    Return a list of count dictionary items
    sorted by keys. Internal use only.

    :: cobj
       type - dict
       desc - a count dictionary to sort
    '''

    return sorted(map(
        lambda x: ([str(_x) for _x in x[0]], x[1]),
        cobj.items()))

def saveitems(iobj, filepath, sortfn):
    '''
    MSGPACK dump sorted items to filepath,
    for downstream consumption.
    Internal use only.

    :: iobj
       type - dict / cx.Counter
       desc - read pack to persist to disk
    :: filepath
       type - string
       desc - filepath to save read pack
    :: sortfn
       type - function
       desc - sorting function used
    '''

    items = sortfn(iobj)
    msgpacksave(
        obj=items,
        filepath=filepath)
    del items

def savepack(pobj, filepath):
    '''
    MSGPACK dump a read pack to filepath.
    Internal use only.

    :: pobj
       type - dict
       desc - read pack to persist to disk
    :: filepath
       type - string
       desc - filepath to save read pack
    '''

    return saveitems(
        iobj=pobj,
        filepath=filepath,
        sortfn=sorteditems)

def savecount(cobj, filepath):
    '''
    MSGPACK dump a count dict to filepath.
    Internal use only.

    :: cobj
       type - dict
       desc - count dict to persist to disk
    :: filepath
       type - string
       desc - filepath to save read pack
    '''

    return saveitems(
        iobj=cobj,
        filepath=filepath,
        sortfn=sortedcounts)

def picklesave(obj, filepath):
    '''
    Pickle dump an object to filepath.
    Internal use only.

    :: obj
       type - object
       desc - object to persist to disk
    :: filepath
       type - string
       desc - filepath to save object
    '''

    # Check Previous FIle Non-Existent
    if os.path.exists(filepath):
        raise RuntimeError(
            '{} exists!'.format(filepath))

    # Pickle Dump Object!
    with open(filepath, 'wb') as outfile:
        pk.dump(
            obj=obj,
            file=outfile,
            protocol=5,
            fix_imports=True)

def savemeta(pobj, filepath):
    '''
    Pickle dump a meta read pack to filepath.
    Internal use only.

    :: pobj
       type - dict
       desc - meta read pack to persist to disk
    :: filepath
       type - string
       desc - filepath to save meta read pack
    '''

    return picklesave(
        obj=pobj,
        filepath=filepath)

def savemodel(mobj, filepath):
    '''
    Pickle dump a TrueSignt model to filepath.
    Internal use only.

    :: mobj
       type - Scry
       desc - TrueSignt model to persist
              to disk
    :: filepath
       type - string
       desc - filepath to save model
    '''

    return picklesave(
        obj=mobj,
        filepath=filepath)

def savedump(dobj, filepath):
    '''
    Pickle dump a dictionary of counting
    callback inputs. Internal use only.

    :: dobj
       type - dict
       desc - dictionary to persist to disk
    :: filepath
       type - string
       desc - filepath to save dictionary
    '''

    return picklesave(
        obj=dobj,
        filepath=filepath)

def msgpackload(incontent):
    '''
    MSGPACK load object from content
    stream. Internal use only.

    :: incontent
       type - bytes
       desc - input content stream
    '''

    return mg.unpackb(
        packed=incontent,
        use_list=True,
        raw=False,
        strict_map_key=False)

def loaddict(archive, dfile):
    '''
    MSGPACK load a dictionary from archive.
    Internal use only.

    :: archive
       type - ZipFile
       desc - archive storing files
    :: dfile
       type - string
       desc - filepath to archived
              dictionary
    '''

    incontent = archive.read(dfile)
    obj = msgpackload(
        incontent=incontent)
    del incontent
    return obj

def loadpack(archive, pfile):
    '''
    MSGPACK load a read pack from archive.
    Internal use only.

    :: archive
       type - ZipFile
       desc - archive storing files
    :: pfile
       type - string
       desc - filepath to archived
              read pack
    '''

    incontent = archive.read(pfile)
    obj = msgpackload(
        incontent=incontent)
    del incontent
    return obj

def loadcount(cfile):
    '''
    MSGPACK load count items from file.
    Internal use only.

    :: cfile
       type - string
       desc - filepath to stored
              count items
    '''

    with open(cfile, 'rb') as infile:
        obj = msgpackload(
            incontent=infile.read())
    return obj

def pickleload(filepath):
    '''
    Pickle load content from filepath.
    Internal use only.

    :: filepath
       type - string
       desc - filepath to saved content
    '''

    with open(filepath, 'rb') as infile:
        obj = pk.load(
            file=infile,
            fix_imports=True)
    return obj

def loadmeta(mfile):
    '''
    Pickle load a meta read pack from
    filepath. Internal use only.

    :: mfile
       type - string
       desc - filepath to saved meta
              read pack
    '''

    return pickleload(
        filepath=mfile)

def loadmodel(archive, mfile):
    '''
    Pickle load a Scry instance from
    archive. Internal use only.

    :: archive
       type - ZipFile
       desc - archive storing files
    :: pfile
       type - string
       desc - filepath to archived
              Scry instance
    '''

    incontent = archive.read(mfile)
    obj = pk.loads(
        incontent,
        fix_imports=True)
    del incontent
    return obj

def loaddump(dfile):
    '''
    Pickle load callback input dump
    from file. Internal use only.

    :: dfile
       type - string
       desc - filepath to stored
              input dump
    '''

    return pickleload(
        filepath=dfile)

# Archive Functions

def get_archive(
    arcfile):
    '''
    Open and return ZipFile archive.
    Internal use only.

    :: arcfile
       type - string
       desc - archive storing all objects
    '''

    return zf.ZipFile(arcfile)

@coroutine
def archive_engine(arcfile, mode='x'):
    '''
    Return a coroutine to archive sent files.
    Internal use only.

    :: arcfile
       type - string
       desc - archive storing all objects
    :: mode
       type - 'x' / 'a'
       desc - archive operation mode
              (default='x')
    '''

    # Delete Previous Archive File
    if mode == 'x':
        remove_file(
            filepath=arcfile)

    # Setup Archive File
    arcfile = zf.ZipFile(
        file=arcfile,
        mode=mode,
        compression=zf.ZIP_DEFLATED,
        compresslevel=2)

    # Arcname Extractor
    arcname = lambda file: file[file.rfind('/')+1:]

    # Online
    try:
        while True:

            # Release Control
            tt.sleep(0)

            # Receive Objects
            obj,liner = (yield)

            # Remnant None Object?
            if obj is None:
                continue

            # Object hasn't Materialized?
            while not os.path.isfile(obj):
                tt.sleep(0)
                continue

            # Decide Archive Entry Name
            arcentry = arcname(obj)

            # Update Archive Progress
            if not liner is None:
                liner.send(
                    ' Archiving {} in Progress'.format(
                        arcentry))

            # Archive Object
            arcfile.write(
                filename=obj,
                arcname=arcentry)

    # Closure
    except GeneratorExit:
        arcfile.close()

def archive(
    objqueue,
    arcfile,
    mode,
    prodcount,
    prodactive,
    liner):
    '''
    Archive objects from paths stored in
    objqueue. Internal use only.

    :: objqueue
       type - SimpleQueue
       desc - queue storing object paths
              to be archived
    :: arcfile
       type - string
       desc - archive storing all objects
    :: mode
       type - 'x' / 'a'
       desc - archive operation mode
              (default='x')
    :: prodcount
       type - integer
       desc - total number of producers
              scheduled to add to objqueue
    :: prodactive
       type - SafeCounter
       desc - total number of producers
              actively adding to objqueue
    :: liner
       type - coroutine
       desc - dynamic printing and
              logging
    '''

    # Define Archiver
    archiver = archive_engine(
        arcfile=arcfile,
        mode=mode)

    # Archive In-Progress
    while prodcount:

        # Release Control
        tt.sleep(0)

        # Waiting ...
        if objqueue.empty():
            continue

        # Do we have a new product?
        obj = objqueue.get()

        # Did a producer finish?
        if obj is None:
            prodcount -= 1
            continue

        # Send Object to Archiver
        if prodactive and \
           prodactive.value() == 0:
            archiver.send((obj, liner))
        else:
            archiver.send((obj, None))

        # Remove Archived File
        os.remove(obj)

    # Archive Completed
    archiver.close()

# Memory Management Functions

def trim_malloc():
    '''
    Trim the heap.
    Internal use only.
    '''

    if trimflag:
        try:
            ct.CDLL('libc.so.6').malloc_trim(0)
        except:
            pass

def free_mem():
    '''
    Free memory used by a process.
    Internal use only.
    '''

    gc.collect()
    trim_malloc()
