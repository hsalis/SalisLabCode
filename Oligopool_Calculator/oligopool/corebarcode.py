import time  as tt

import collections as cx

import bisect as bs
import dinopy as dp
import numpy  as np
from numpy.lib.function_base import percentile

import utils  as ut


# Barcode Conversion Dictionary
decoder = dict(zip((0., 1., 2., 3.), 'AGTC'))

# Parser and Setup Functions

def get_parsed_barcode_length(
    barcodelen,
    indf,
    liner):
    '''
    Check feasibility of barcodelen
    for given indf, i.e. reachability.
    Internal use only.

    :: barcodelen
       type - integer
       desc - required barcode length
    :: indf
       type - pd.DataFrame
       desc - input pandas DataFrame
    :: liner
       type - coroutine
       desc - dynamic printing
    '''

    # Start Timer
    t0 = tt.time()

    # Compute barcodelen feasibility
    designspace = 4 ** barcodelen
    targetcount = len(indf.index)

    lenspace = max(designspace, targetcount)
    plen = ut.get_printlen(
        value=lenspace)
    sntn, plen = ut.get_notelen(
        printlen=plen)

    parsestatus = designspace >= targetcount

    if not parsestatus:
        parsemsg = ' [INFEASIBLE] (Design Space is Smaller than Target Space)'
    else:
        parsemsg = ''

    # Show update
    liner.send(
        ' Required Length: {:,} Base Pair(s)\n'.format(
            barcodelen))
    liner.send(
        '   Design Space : {:{},{}} Barcode(s)\n'.format(
            designspace,
            plen,
            sntn))
    liner.send(
        '   Target Count : {:{},{}} Barcode(s){}\n'.format(
            targetcount,
            plen,
            sntn,
            parsemsg))
    liner.send(
        ' Time Elapsed: {:.2f} sec\n'.format(
            tt.time() - t0))

    if parsestatus:
        liner.send(
            ' Verdict: Barcode Design Possibly Feasible\n')
    else:
        liner.send(
            ' Verdict: Barcode Infeasible due to Barcode Length Constraint\n')

    # Return results
    return (parsestatus,
        designspace,
        targetcount)

def get_finite_jumper(upb):
    '''
    Return finite space jumper.
    Internal use only.

    :: upb
       type - integer
       desc - space uppberbound
    '''

    # Setup Permutation
    hfb = upb // 2
    qtb = upb // 4
    fld = np.random.permutation(qtb)

    # Split Field Index
    idx = 0

    # Stream Coordinates
    while idx < len(fld):
        coo = fld[idx]
        yield coo
        yield upb - 1 - coo
        yield hfb + coo
        yield hfb - 1 - coo
        idx += 1

def get_infinite_jumper(upb):
    '''
    Return infinite space jumper.
    Internal use only.

    :: upb
       type - integer
       desc - space upperbound
    '''

    # Setup RNG Permutation
    rng = np.random.default_rng()

    # Stream Coordinates
    while True:
        yield rng.integers(0, upb)

def get_jumper(barcodelen):
    '''
    Return jumper based on jumper type.
    Internal use only.

    :: barcodelen
       type - integer
       desc - barcode length
    '''

    # Compute Upperbound
    upb = 4**barcodelen

    # Return Jumper
    if barcodelen <= 12:
        # Finite   Random Jumps
        return 1, get_finite_jumper(upb=upb)
    else:
        # Infinite Random Jumps
        return 2, get_infinite_jumper(upb=upb)

# Engine Objective and Helper Functions

def stream_barcodes(
    barcodelen,
    jumper):
    '''
    Randomly jump through sequence space and
    stream barcodes. Internal use only.

    :: barcodelen
       type - integer
       desc - required barcode length
    :: jumper
       type - generator
       desc - generator streaming barcodes
    '''

    # Transformer Setup
    txf = lambda x: float(ord(x)) - 48.

    # Stream Jumps!
    while True:
        # Fetch a Jump Coordinate
        try:
            idx = next(jumper)
            # print(idx)
        except StopIteration:
            yield None # Finite Iterator Exhaused

        # Stream Coordinate Transformed to Barcode
        yield np.array(tuple(map(
                txf,
                np.base_repr(
                    idx,
                    base=4).zfill(
                        barcodelen))),
            dtype=np.float64)

def show_update(
    count,
    plen,
    barcodeseq,
    optstatus,
    optstate,
    inittime,
    terminal,
    liner):
    '''
    Display the current progress in barcode
    generation. Internal use only.

    :: count
       type - integer
       desc - current storage fill count
    :: plen
       type - integer
       desc - barcode index width
    :: barcodeseq
       type - string
       desc - decoded barcode sequence
    :: optstatus
       type - boolean
       desc - barcode feasibility status
    :: optstate
       type - integer
       desc - feasibility failure state marker
    :: inittime
       type - tt.time
       desc - initial time stamp
    :: terminal
       type - boolean
       desc - if True will terminate update to newline
              otherwise, rewrite previous update
    :: liner
       type - coroutine
       desc - dynamic printing
    '''

    if len(barcodeseq) >= 30:
        design = barcodeseq[:14] + '..' + barcodeseq[-14:]
    else:
        design = barcodeseq

    liner.send(' Candidate {:{},d}: Barcode {} is {}{}'.format(
        count,
        plen,
        design,
        ['Rejected', 'Accepted'][optstatus],
        ['',
        ' due to Hamming Distance',
        ' due to Excluded Motif',
        ' due to Edge Effect',
        ' due to Oligopool Repeat',
        ' due to Composition'][optstate]))

    if terminal:
        liner.send('\* Time Elapsed: {:.2f} sec\n'.format(
            tt.time() - inittime))

def get_barcodeseq(barcode):
    '''
    Return the decoded sequence from
    barcode coordinate array.
    Internal use only.

    :: barcode
       type - np.array
       desc - encoded barcode array
    '''

    return ''.join(map(
        lambda x: decoder[x],
        barcode))

def get_contigsize_scheme(
    barcodelen,
    minhdist):
    '''
    Return optimal barcode contig parition
    scheme for coordinate checking.
    Internal use only.

    :: barcodelen
       type - integer
       desc - barcode length
    :: minhdist
       type - integer
       desc - minimum pairiwise hamming
              distance between barcodes
    '''

    contigsize = np.zeros(minhdist, dtype=np.int32) + (barcodelen // minhdist)
    covered = np.sum(contigsize)
    if covered < barcodelen:
        contigsize[:barcodelen-covered] += 1
    return contigsize

def is_hamming_feasible(
    barcodeseq,
    store,
    count,
    minhdist,
    contigsize,
    coocache):
    '''
    Determine if the minimum distance between
    barcode and store is greater or equal to
    given minhdist (hamming feasibility).
    Internal use only.

    :: barcodeseq
       type - string
       desc - decoded barcode sequence
    :: store
       type - np.array
       desc - vectorized storage of numeric
              encoding of all previous barcodes
    :: count
       type - integer
       desc - current storage fill count
    :: minhdist
       type - integer
       desc - minimum pairwise hamming
              distance between a pair
              of barcodes
    :: contigsize
       type - array
       desc - optimized barcode contig
              size scheme
    :: coocache
       type - dict
       desc - dictionary of k-mer contig
              coordinates
    '''

    # Bound Targets
    targets = set()
    coordinates = []
    maxsize = store.shape[1]
    bounded = True
    for contig,index in ut.stream_contigs(
        seq=barcodeseq,
        scheme=contigsize):
        if contig in coocache:
            targets.update(
                coocache[contig][index])
        if len(targets) == maxsize:
            bounded = False
        coordinates.append((contig,index))
    targets = tuple(targets)

    # Did we match coordinates?
    if not targets:
        return (True,coordinates)
    else:
        if bounded:
            targets = np.array(targets)
        else:
            targets = None

    # Compute Distance
    return (ut.get_store_hdist(
        store=store,
        idx=count,
        t=targets,
        direction=0) >= minhdist,
        coordinates)

def is_exmotif_feasible(
    barcodeseq,
    barcodelen,
    exmotifs):
    '''
    Determine if the barcode does not contain any
    of the exluded motifs. Internal use only.

    :: barcodeseq
       type - string
       desc - decoded barcode sequence
    :: barcodelen
       type - string
       desc - barcode length
    :: exmotifs
       type - list / None
       desc - list of motifs to exclude
              in designed barcodes,
              None otherwise
    '''

    return ut.get_exmotif_conflict(
        seq=barcodeseq,
        seqlen=barcodelen,
        exmotifs=exmotifs,
        partial=True,
        checkall=False)

def is_edge_feasible(
    barcodeseq,
    exmotifs,
    contextarray,
    leftselector,
    leftcontexttype,
    rightselector,
    rightcontexttype):
    '''
    Determine the context assignment index
    for designed sequence, and see if the
    assignment is devoid of edge effects.
    Internal use only.

    :: barcodeseq
       type - string
       desc - decoded barcode sequence
    :: exmotifs
       type - list / None
       desc - list of motifs to exclude
              in designed barcodes,
              None otherwise
    :: contextarray
       type - np.array
       desc - context assignment array
    :: leftcontexttype
       type - integer
       desc - inferred left context type
              1 = list
              2 = string
              3 = None
    :: leftselector
       type - lambda
       desc - selector for the left
              sequence context
    :: rightcontexttype
       type - integer
       desc - inferred right context type
              1 = list
              2 = string
              3 = None
    :: rightselector
       type - lambda
       desc - selector for the right
              sequence context
    '''

    # Book-keeping
    i = 0
    t = len(contextarray)
    cfails = cx.Counter()

    # Loop through contexts for assignment
    while i < t:

        # Fetch Context
        aidx = contextarray.popleft()

        # Define in-context sequence
        incntxseq = barcodeseq

        # Fetch left and right contexts
        llen  = 0
        lcntx =  leftselector(aidx)
        rcntx = rightselector(aidx)

        # Build out in-context sequence
        if not lcntx is None:
            llen = len(lcntx)
            incntxseq = lcntx + incntxseq
        if not rcntx is None:
            incntxseq = incntxseq + rcntx

        # Compute Feasibility
        mcond, motif = True, None
        if (not  exmotifs is None) and \
           ((not lcntx    is None) or \
            (not rcntx    is None)):

            # Book-keeping
            p = 0
            q = llen
            r = llen + len(barcodeseq)
            s = len(incntxseq)

            # Loop through exmotifs
            for motif in exmotifs:

                # Locate motif start in in-context
                start = incntxseq.find(motif)

                # Found a match!
                if start > -1:

                    # Compute motif end in in-context
                    end = start + len(motif)

                    # Ref:
                    # p      q                r       s
                    # |------|----------------|-------|
                    if ((p <= start < q) and (q+1 <= end < r)) or \
                       ((q <= start < r) and (r+1 <= end < s)):
                            mcond = False
                            break

        # We have a winner folks!
        if mcond:
            return (True, aidx, None)

        # We lost an assignment, record cause
        else:
            cfails[motif] += 1

        # Update Iteration
        i += 1

        # We will try this context again
        contextarray.append(aidx)

        # Did we deal with constant contexts
        # on both sides of the barcode?
        if  leftcontexttype  in [2, 3] and \
            rightcontexttype in [2, 3]:
            break

    # Failed Assignment
    return (False, None, cfails)

def is_nonrepetitive(
    barcodeseq,
    maxreplen,
    index,
    oligorepeats):
    '''
    Determine if barcode is contextually
    non-repetitive. Internal use only.

    :: barcodeseq
       type - string
       desc - decoded barcode sequence
    :: maxreplen
       type - integer
       desc - maximum shared repeat length
              between barcode and context
    :: index
       type - integer
       desc - context assignment index
              for barcode sequences
    :: oligorepeats
       type - dict
       desc - set of all repeats indexed
              by context id
    '''

    corpus = oligorepeats[index]
    for kmer in ut.stream_canon_spectrum(
        seq=barcodeseq,
        k=maxreplen+1):
        kmer = dp.conversion.encode_twobit(
            seq=kmer,
            sentinel=False)
        idxk = bs.bisect_left(
            corpus,
            kmer)
        if (idxk < len(corpus)) and \
           (corpus[idxk] == kmer):
            return False
    return True

def is_composition_feasible(
    barcodeseq,
    substringlen,
    substringcache,):
    '''
    Determine if barcode is terminally
    non-repetitive. Internal use only.

    :: barcodeseq
       type - string
       desc - decoded barcode sequence
    :: substringlen
       type - integer
       desc - barcode substring extraction
              length
    :: substringcache
       type - integer
       desc - barcode substring storage
    '''

    subset = []
    for substring in ut.stream_spectrum(
        seq=barcodeseq,
        k=substringlen):
        if substringcache[substring] > 0:
            return False, None
        subset.append(substring)
    return True, subset

def barcode_objectives(
    store,
    count,
    barcodeseq,
    barcodelen,
    substringlen,
    substringcache,
    minhdist,
    contigsize,
    coocache,
    maxreplen,
    oligorepeats,
    exmotifs,
    contextarray,
    leftcontexttype,
    leftselector,
    rightcontexttype,
    rightselector):
    '''
    Determine if a barcode satisfies all
    global objectives. Internal use only.

    :: store
       type - np.array
       desc - vectorized storage of numeric
              encoding of all previous
              barcodes
    :: count
       type - integer
       desc - current storage fill count
    :: barcodeseq
       type - string
       desc - decoded barcode sequence
    :: barcodelen
       type - integer
       desc - barcode length
    :: substringlen
       type - integer
       desc - barcode substring length
    :: substringcache
       type - set
       desc - barcode substring storage
    :: minhdist
       type - integer
       desc - minimum pairwise hamming
              distance between a pair
              of barcodes
    :: contigsize
       type - array
       desc - optimized barcode contig
              size scheme
    :: coocache
       type - dict
       desc - dictionary of k-mer contig
              coordinates
    :: maxreplen
       type - integer
       desc - maximum shared repeat length
              between barcode and context
    :: oligorepeats
       type - dict
       desc - set of all repeats indexed
              by context id
    :: exmotifs
       type - list / None
       desc - list of motifs to exclude
              in designed barcodes,
              None otherwise
    :: contextarray
       type - np.array
       desc - context assignment array
    :: leftcontexttype
       type - integer
       desc - inferred left context type
              1 = list
              2 = string
              3 = None
    :: leftselector
       type - lambda
       desc - selector for the left
              sequence context
    :: rightcontexttype
       type - integer
       desc - inferred right context type
              1 = list
              2 = string
              3 = None
    :: rightselector
       type - lambda
       desc - selector for the right
              sequence context
    '''

    # Objective 1: Hamming Distance
    obj1,coordinates = is_hamming_feasible(
        barcodeseq=barcodeseq,
        store=store,
        count=count,
        minhdist=minhdist,
        contigsize=contigsize,
        coocache=coocache)

    # Objective 1 Failed
    if not obj1:
        return (False, 1, None, None, None, None)         # Hamming Failure

    # Objective 2: Motif Embedding
    obj2, exmotif = is_exmotif_feasible(
        barcodeseq=barcodeseq,
        barcodelen=barcodelen,
        exmotifs=exmotifs)

    # Objective 2 Failed
    if not obj2:
        return (False, 2, None, {exmotif: 1}, None, None) # Motif Failure

    # Objective 3: Edge Feasibility (Edge-Effects)
    obj3, aidx, afails = is_edge_feasible(
        barcodeseq=barcodeseq,
        exmotifs=exmotifs,
        contextarray=contextarray,
        leftcontexttype=leftcontexttype,
        leftselector=leftselector,
        rightcontexttype=rightcontexttype,
        rightselector=rightselector)

    # Objective 3 Failed
    if not obj3:
        return (False, 3, None, afails, None, None)       # Edge Failure

    # Objective 4: Contextual Non-Repetitiveness
    obj4 = is_nonrepetitive(
        barcodeseq=barcodeseq,
        maxreplen=maxreplen,
        index=aidx,
        oligorepeats=oligorepeats)

    # Objective 4 Failed
    if not obj4:
        # Re-insert Assignment Back in Context
        # Otherwise, Bug: Assigned Index is Lost
        #                 Followed by Infinite Loop
        #                 (Critical Bug)
        contextarray.append(aidx)
        return (False, 4, None, None, None, None)         # Repeat Failure

    # Objective 5: Composition Objective
    subset = None
    if barcodelen - substringlen > 0: # We have Composition to Optimize
        obj5, subset = is_composition_feasible(
            barcodeseq=barcodeseq,
            substringlen=substringlen,
            substringcache=substringcache)

        # Objective 5 Failed
        if not obj5:
            # Re-insert Assignment Back in Context
            # Otherwise, Bug: Assigned Index is Lost
            #                 Followed by Infinite Loop
            #                 (Critical Bug)
            contextarray.append(aidx)
            return (False, 5, None, None, None, None)     # Repeat Failure

    # All conditions met!
    return (True, 0, aidx, None, subset, coordinates)

def get_distro(
    store,
    liner):
    '''
    Compute and return the hamming
    distance distribution of store.
    Internal use only.

    :: store
       type - np.array
       desc - vectorized storage of numeric
              encoding of all previous
              barcodes
    :: liner
       type - coroutine
       desc - dynamic printing
    '''

    # Book-keeping
    t0 = tt.time()
    hammingdistro = cx.Counter()
    count     = store.shape[0]
    threshold = min(min(50000, count // 10), count)
    plen      = ut.get_printlen(
        value=count)

    # Show Update
    liner.send(' Selecting Samples ...')

    # Select Samples
    if threshold < count:
        # store = store[np.lexsort(np.rot90(store))]
        idxsample = sorted(np.random.permutation(
            count)[:threshold])
        store = store[idxsample, :]
    else:
        idxsample = np.arange(count)

    # Compute DIstribution - Half per Edge
    for idx in np.arange(store.shape[0]):

        # Upward Distance
        hdist = ut.get_store_hdist(
            store=store,
            idx=idx,
            direction=0)

        # Show Update
        liner.send(
            ' Candidate {:{},d}: Pairwise Distance ≥ {:,} Mismatches'.format(
                idxsample[idx]+1, plen, hdist))

        # Add to Distribution
        hammingdistro[hdist] += 1

    # Normalize Distribution
    liner.send('\* Normalizing Distribution ...')
    normhammingdistro = []
    for hdist in hammingdistro:
        percentage = (100. * hammingdistro[hdist]) / threshold
        normhammingdistro.append((percentage, hdist))
    normhammingdistro.sort()

    # Final Update
    liner.send(
        ' Time Elapsed: {:.2f} sec\n'.format(
            tt.time()-t0))

    # Return Distribution
    return normhammingdistro
