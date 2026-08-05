import time  as tt

import collections as cx

import numpy as np
import numba as nb

import utils as ut


# Parser and Setup Functions

def get_parsed_splitlimit(
    indf,
    splitlimit,
    liner):
    '''
    Determine if splitting is feasible
    based on oligo synthesis limit.
    Internal use only.

    :: indf
       type - pd.DataFrame
       desc - input DataFrame containing
              all designed variants
    :: splitlimit
       type - integer
       desc - maximum allowed oligo length
              after splitting
    :: liner
       type - coroutine
       desc - dynamic printing
    '''

    # Book-keeping
    t0 = tt.time()
    lengthdiff    = round(splitlimit * 0.25)
    maxvariantlen = None
    minsplitcount = None
    maxsplitcount = None
    lengthcond    = False
    splitcond     = False

    # Finalize Oligopool
    liner.send(' Finalizing Oligopool ...')
    seqlist = ut.get_df_concat(df=indf)

    # _seqlist = []
    # for idx in range(len(seqlist)):
    #     seq = seqlist[idx]
    #     seq = seq[:np.random.randint(182, 326)]
    #     # print(len(seq))
    #     _seqlist.append(seq)
    # seqlist = _seqlist

    # Compute Variant Lengths
    liner.send(' Computing Variant Lengths ...')

    variantlens   = list(map(len, seqlist))
    minvariantlen = min(variantlens)
    maxvariantlen = max(variantlens)

    # variantlen = 190
    # maxvariantlen = 190

    # Compute Length Feasibility
    if minvariantlen < splitlimit:
        parsemsg = ' [INFEASIBLE] (Variants Shorter than Split Limit)'
    elif splitlimit - (maxvariantlen - minvariantlen) < lengthdiff:
        # Note: Length Differential is more than 1/4 th of Split Limit
        parsemsg = ' [INFEASIBLE] (Lengths Differ by More than {:,} bp)'.format(
            lengthdiff)
    else:
        parsemsg   = ''
        lengthcond = True

    # Show Updates
    if minvariantlen == maxvariantlen:
        liner.send(
            ' Input Variant Length: {:,} Base Pair(s){}\n'.format(
                minvariantlen,
                parsemsg))
    else:
        liner.send(
            ' Input Variant Length: {:,} to {:,} Base Pair(s){}\n'.format(
                minvariantlen,
                maxvariantlen,
                parsemsg))

    # Computing Split Counts
    minsplitcount = int(np.ceil(minvariantlen / (splitlimit * 1.)))
    maxsplitcount = int(np.ceil(maxvariantlen / (splitlimit * 1.)))

    # minsplitcount = 2
    # maxsplitcount = 3

    # Compute Feasibility
    splitcond = minsplitcount == maxsplitcount
    if not splitcond:
        parsemsg = ' [INFEASIBLE] (Uneven Number of Splits)'
    else:
        parsemsg = ''

    # Show Updates
    if minsplitcount == maxsplitcount:
        liner.send(
            ' Split Fragment Count: At least {:,} Split(s) per Variant{}\n'.format(
                minsplitcount,
                parsemsg))
    else:
        liner.send(
            ' Split Fragment Count: At least {:,} to {:,} Split(s) per Variant{}\n'.format(
                minsplitcount,
                maxsplitcount,
                parsemsg))

    # Show Time Elapsed
    liner.send(
        ' Time Elapsed: {:.2f} sec\n'.format(
            tt.time()-t0))

    # Compute Verdict
    parsestatus = lengthcond and splitcond
    if not parsestatus:
        liner.send(
            ' Verdict: Splitting Infeasible due to Split Limit Constraints\n')
    else:
        liner.send(
            ' Verdict: Splitting Possibly Feasible\n')

    # Return Results
    return (parsestatus,
        seqlist,
        not lengthcond,
        not splitcond,
        minvariantlen,
        maxvariantlen,
        minsplitcount,
        maxsplitcount)

def get_seqvec(
    seq):
    '''
    Return the numeric vector for seq.
    Internal use only.

    :: seq
       type - string
       desc - a sequence to split
    '''

    return np.array(
        tuple(float(ord(nt)) for nt in seq),
        dtype=np.float64)

def get_padded_seq(
    seq,
    diff):
    '''
    Return 3' padded sequence.
    Internal use only.

    :: seq
       type - string
       desc - sequence to pad
    :: diff
       type - integer
       desc - amount of 3' padding
    '''

    return seq + ''.join(np.random.choice(
        list('ATGC')) for _ in range(diff))

def get_seqmat_padvec(
    seqlist,
    maxvariantlen,
    liner):
    '''
    Return the numeric representation
    of seqlist and the vector of added
    padding. Internal use only.

    :: seqlist
       type - iterable
       desc - list of sequences to split
    :: maxvariantlen
       type - integer
       desc - length of the longest variant
              in the pool
    :: liner
       type - coroutine
       desc - dynamic printing
    '''

    # Setup Stores
    seqmat = np.zeros(
        (len(seqlist), maxvariantlen),
        dtype=np.float64)
    padvec = np.array([maxvariantlen - len(seq) \
        for seq in seqlist])

    # Update-keeping
    short = True if maxvariantlen > 10 else False

    # Time-keeping
    t0 = tt.time()

    # Fill Store
    for idx,seq in enumerate(seqlist):
        diff   = padvec[idx]
        padseq = get_padded_seq(seq=seq, diff=diff)
        seqvec = get_seqvec(seq=padseq)
        seqmat[idx, :] = seqvec

        if short:
            liner.send(
                ' Storing Vectorized Sequence {}: {}..{}'.format(
                    idx,
                    seqvec[:5],
                    seqvec[-5:]))
        else:
            liner.send(
                ' Storing Vectorized Sequence {}: {}'.format(
                    idx,
                    seqvec))

    # Final Updates
    plen = ut.get_printlen(
        value=max(idx+1, padvec.min()))
    liner.send(
        '   Vectorized: {:{},d} Sequences\n'.format(
            idx+1,
            plen))
    if padvec.min() == padvec.max():
        liner.send(
            '   3\' Padding: {:{},d} Base Pair(s)\n'.format(
                padvec.min(),
                plen))
    else:
        liner.send(
            '   3\' Padding: {:{},d} to {:,} Base Pair(s)\n'.format(
                padvec.min(),
                plen,
                padvec.max()))
    liner.send(' Time Elapsed: {:.2f} sec\n'.format(
        tt.time()-t0))

    # Return Results
    seqmat = np.array(seqmat, dtype=np.float64)
    return padvec, seqmat

def get_entvec(
    seqmat,
    maxvariantlen,
    liner):
    '''
    Return entropy of sequence matrix.
    Internal use only.

    :: seqmat
       type - np.array
       desc - numeric sequence matrix
    :: maxvariantlen
       type - integer
       desc - length of the longest variant
              in the pool
    :: liner
       type - coroutine
       desc - dynamic printing
    '''

    # Update-keeping
    short = True if maxvariantlen > 10 else False

    # Time-keeping
    t0  = tt.time()

    # Show Update
    liner.send(' Computing Count Vector ...')

    # Setup Data Structures
    d = np.zeros(
        (4, seqmat.shape[1]),
        dtype=np.float64)
    p = np.zeros(4)

    # Compute Count Frequency
    for idx in range(seqmat.shape[1]):

        # Count Symbol Counts at i-th Index
        m = np.unique(
            seqmat[:, idx],
            return_counts=True)[1]

        # Update Data Structure and Show Updates
        p[:m.shape[0]] = m # Absorb Local

        liner.send(' Index {} Unique Count: {}'.format(
            idx, p))

        d[:, idx] += p     # Update Global
        p *= 0.            # Reset  Local

    # Normalize Counts
    d = d / d.sum(0)

    # Show Updates
    liner.send('   Count Vector: Normalized\n')

    # Show Updates
    liner.send(' Computing Entropy Vector ...')

    # Compute Entropy Vector
    entvec = cx.deque(np.abs(
        (d*(np.log(d, where=d>0.) / np.log(4))).sum(0)))

    # Show Updates
    liner.send(' Entropy Vector: Computed\n')
    liner.send(' Time Elapsed: {:.2f} sec\n'.format(
        tt.time() - t0))

    # Return Results
    return entvec

def get_base_varcont(
    entvec):
    '''
    Return all base variable region span indices.
    Internal use only.

    :: entvec
       type - cx.deque
       desc - positional entropy
    '''

    # Setup Parsing
    varcont = cx.deque()
    start   = None
    end     = None
    idx     = -1

    # Time-keeping
    t0 = tt.time()

    # Parse Contigs
    while entvec:

        # Update index
        idx += 1

        # Extract entropy
        ent = entvec.popleft()

        # Constant Region
        if ent <= 0.25: # Constant Upper Bound

            # A contig built?
            if not start is None and \
               not end   is None:
                varcont.append((start, end))

            # Reset for next contig
            start = None
            end   = None

            # Next
            continue

        # Variable Region
        else:

            # Contig continues?
            if start is None:
                start = idx
            if end is None:
                end = idx

            # Update ending index
            end += 1

            # Next
            continue

    # Remnant Update
    if not end is None:
        varcont.append((start, end))

    # Return Results
    return varcont

def get_merged_varcont(
    varcont,
    mergegap):
    '''
    Return a merged varcont, merging variable regions
    separated by at most mergegap constant bases.
    Internal use only.

    :: varcont
       type - cx.deque
       desc - all variable region span indices
    :: mergegap
       type - integer
       desc - maximum gap length between two
              variable regions to be merged
              into a single contig
    '''

    # Do we merge?
    if not varcont:
        return varcont

    # Setup Data Structures
    merged   = cx.deque()
    previous = list(varcont.popleft())

    # Time-keeping
    t0 = tt.time()

    # Merge Contigs
    while varcont:
        current = varcont.popleft()
        if current[0]-previous[1] <= mergegap:
            previous[1] = current[1]
        else:
            merged.append(tuple(previous))
            previous = list(current)
    merged.append(tuple(previous))

    # Return Result
    return merged

def is_spannable(
    p,
    q,
    spanlen):
    '''
    Determine if a contig span satisfies spanlen.
    Internal use only.

    :: p
       type - integer
       desc - span start
    :: q
       type - integer
       desc - span end
    :: spanlen
       type - integer
       desc - minimum required split span length
    '''

    if q - p < spanlen:
        return False
    return True

def get_filtered_varcont(
    varcont,
    spanlen):
    '''
    Filter all variable regions shorter than
    the required spanlen. Internal use only.

    :: varcont
       type - cx.deque
       desc - a deque of tuple with start and
              end coordinates of variable regions
    :: spanlen
       type - integer
       desc - minimum required split span length
    '''

    # Do we filter varcont?
    if not varcont:
        return varcont

    # Setup Data Structure
    filtered = cx.deque()

    # Time-keeping
    t0 = tt.time()

    # Filter Contigs
    for p,q in varcont:
        if is_spannable(p, q, spanlen):
            filtered.append((p, q))

    # Return Result
    return filtered

def get_varcont(
    entvec,
    minhdist,
    spanlen,
    liner):
    '''
    Return all valid variable contig span regions.
    Internal use only.

    :: entvec
       type - cx.deque
       desc - positional entropy
    :: minhdist
       type - integer
       desc - minimum pairwise hamming distance
    :: spanlen
       type - integer
       desc - minimum required split span length
    :: liner
       type - coroutine
       desc - dynamic printing
    '''

    # Extract varcont
    liner.send(' Extracting Variable Contigs ...')
    t0 = tt.time()
    varcont = get_base_varcont(
        entvec=entvec)
    varcontcount = len(varcont)
    liner.send(
        ' Extracted Variable Contigs: {}\n'.format(
            varcontcount))

    # Merge varcont
    liner.send('   Merging Variable Contigs ...')
    t0 = tt.time()
    varcont = get_merged_varcont(
        varcont=varcont,
        mergegap=min(3, minhdist // 4))
    mergedcontcount = len(varcont)
    liner.send(
        '    Merged Variable Contigs: {}\n'.format(
            mergedcontcount))

    # Filter varcont
    liner.send(' Filtering Variable Contigs ...')
    t0 = tt.time()
    varcont = get_filtered_varcont(
        varcont=varcont,
        spanlen=spanlen)
    filtercontcount = len(varcont)
    liner.send(
        '  Filtered Variable Contigs: {}\n'.format(
            filtercontcount))

    # Compute Feasibility
    liner.send(' Time Elapsed: {:.2f} sec\n'.format(
        tt.time()-t0))

    parsestatus = len(varcont) > 0
    if not parsestatus:
        liner.send(
            ' Verdict: Splitting Infeasible due to Lack of Sequence Diversity\n')
    else:
        liner.send(
            ' Verdict: Splitting Possibly Feasible\n')

    # Return Results
    return (parsestatus,
        varcont,
        varcontcount,
        mergedcontcount,
        filtercontcount)

# Engine Objective and Helper Functions

def is_varcont_feasible(
    varcont,
    seqlen,
    splitlen,
    spanlen,
    liner):
    '''
    Determine if the variable contigs are
    within splitlen range, otherwise there
    is no solution to problem instance.
    Internal use only.

    :: varcont
       type - cx.deque
       desc - all variable region span indices
    :: seqlen
       type - integer
       desc - length of sequences to split
    :: splitlen
       type - integer
       desc - maximum split length
    :: spanlen
       type - integer
       desc - minimum required split span length
    :: liner
       type - coroutine
       desc - dynamic printing
    '''

    # Book-keeping
    ci     = iter(varcont)
    pp, qq = next(ci)

    # Check for the first fragment
    if pp + spanlen > splitlen:
        liner.send(
            ' Verdict: Infeasible (First Contig (Start={}, End={}) too Far from Beginning)\n'.format(
                pp, qq))
        return False

    # Check between contigs
    for p,q in varcont:
        if p - qq + 2*spanlen > splitlen:
            liner.send(
                ' Verdict: Infeasible (Adjacent Contigs (Start={}, End={}) and (Start={}, End={}) Far Apart)\n'.format(
                    pp, qq, p, q))
            return False
        qq = q
        pp = p

    # Check for the last fragment
    if qq - spanlen + splitlen < seqlen:
        liner.send(
            ' Verdict: Infeasible (Last Contig (Start={}, End={}) too far from Ending)\n'.format(
                pp, qq))
        return False

    # No problems found
    return True

def get_splitend(
    fstart,
    splitlimit,
    variantlen):
    '''
    Return the fragment end coordinates.
    Internal use only.

    :: fstart
       type - integer
       desc - current fragment starting point
    :: splitlimit
       type - integer
       desc - maximum length of split oligo
    :: variantlen
       type - integer
       desc - variant length considered
    '''

    return min(fstart+splitlimit, variantlen)

def get_splitqueue(
    varcont,
    fstart,
    sstart,
    splitlimit,
    variantlen):
    '''
    Return all variable contigs splitpoints given
    current oligo starting point and splitlen.
    Internal use only.

    :: varcont
       type - cx.deque
       desc - all variable region span indices
    :: fstart
       type - integer
       desc - current fragment starting point
    :: sstart
       type - integer
       desc - current split starting point
    :: splitlimit
       type - integer
       desc - maximum length of split oligo
    :: variantlen
       type - integer
       desc - variant length considered
    '''

    # Do we compute splitpoints?
    if not varcont:
        return varcont

    # Setup Endpoint Queue
    spq = cx.deque()
    end = get_splitend(
        fstart=fstart,
        splitlimit=splitlimit,
        variantlen=variantlen)

    # Determine Potential Splitpoints from Contigs
    for p,q in varcont:

        # Skip Condition:
        # Contig Ends before Start
        if q <= sstart:
            continue # Skip!

        # Absorb Condition:
        # Contig Starts before End
        if p < end:

            # Add Splitpoint
            spq.appendleft((
                max(p, sstart),
                min(q, end)))

            # Contig before End
            if q <  end:
                continue # Next!

            # Contig encompasses End
            if end <= q:
                break # We're done!

        # Exhaustion Condition:
        # Contig Starts after End
        if p >= end:
            break # We're done!

    # Return Results
    return spq

def continue_splitting(
    fstart,
    splitlimit,
    variantlen):
    '''
    Determine if the current split engenders the last
    frament from splitting. Internal use only.

    :: fstart
       type - integer
       desc - current fragment starting point
    :: splitlimit
       type - integer
       desc - maximum length of split oligo
    :: variantlen
       type - integer
       desc - variant length considered
    '''

    # Starting at fstart, we'll
    # reach end of sequence
    if get_splitend(
        fstart=fstart,
        splitlimit=splitlimit,
        variantlen=variantlen) == variantlen:
        return False # No more splitting required

    return True # Splitting may be required

def get_split(
    seqlist,
    seqmat,
    spq,
    mintmelt,
    minhdist,
    spanlen,
    maxoverlap,
    liner):
    '''
    Return an integer r (p < r < q) from spq intervals such
    that r < q - maxoverlap, Tm(runmat[r:q]) > mintmelt and
    HD(seqmat) > minhdist. Internal use only.

    :: seqlist
       type - list
       desc - list of sequences to split
    :: seqmat
       type - np.array
       desc - numeric sequence matrix
    :: spq
       type - deque
       desc - endpoint variable contigs
    :: mintmelt
       type - float
       desc - melting temperature lower bound
    :: minhdist
       type - integer
       desc - minimum pairwise hamming distance
    :: spanlen
       type - integer
       desc - minimum required split span length
    :: maxoverlap
       type - integer
       desc - maximum span of a split region
    :: liner
       type - coroutine
       desc - dynamic printing
    '''

    # Book-keeping
    r, q  = None, None # Split Coordinates
    state = False      # Solution State

    # Splitting Loop
    while spq:

        # Fetch Current Splitpoints
        p, q = spq.popleft()

        # Show Update
        liner.send('    Attempting Split in Region: (Start={}, End={})\n'.format(
            p, q))

        # Is current splitpoint feasigle?
        if not is_spannable(p, q, spanlen):
            liner.send('      Current Split Region of {} bp Infeasible ... Skipping\n'.format(
                q-p))
            state = False # No Solution (Yet)
            continue
        else:
            liner.send('      Current Split Region of {} bp maybe Feasible ... Optimizing\n'.format(
                q-p))

        # Adjusted r Value
        r = q - spanlen # Takes care of spanlen

        # Constraint Match Tracers
        condol = True
        condtm = False
        condhd = False

        # Split Adjustment Loop
        idx = 0 # Current Sequence for Verification
        while idx < seqmat.shape[0]: # Adjust/Verify wrt all Sequences

            # tt.sleep(1)

            # Optimize Overlap
            if r < q - maxoverlap:
                # Show Update
                liner.send(
                    '      Sequence {:,}: Region (Start={:,}, End={:,}) has LONGER Overlap ({:,} bp)'.format(
                        idx, r, q, q-r))

                # Too big of an overlap
                condol = False

            else:
                condol = True

                # Show Update
                liner.send(
                    '      Sequence {:,}: Region (Start={:,}, End={:,}) is Overlap ({:,} bp) Optimal'.format(
                        idx, r, q, q-r))


            # Optimize Tm (after Overlap)
            if condol and not condtm:

                # Compute Tm for current split
                tmelt = ut.get_tmelt(
                    seq=seqlist[idx],
                    i=r,
                    j=q)

                # Tm was Lower ..
                if tmelt < mintmelt:

                    # Show Update
                    liner.send(
                        '      Sequence {:,}: Region (Start={:,}, End={:,}) has LOWER Tm ({:.2f} C)'.format(
                            idx, r, q, tmelt))

                    # Minimize r Value
                    r = r - 1

                # Tm was OK!
                elif tmelt >= mintmelt:
                    condtm = True

                    # Show Update
                    liner.send(
                        '      Sequence {:,}: Region (Start={:,}, End={:,}) is Tm ({:.2f} C) Optimal'.format(
                            idx, r, q, tmelt))

            # Optimize Hamming Distance (after Tm)
            if condol and condtm and not condhd:

                # Compute Hamming Distance for current split
                hdist = ut.get_store_hdist(
                    store=seqmat,
                    idx=idx,
                    i=r,
                    j=q,
                    direction=0)

                # HDist was Lower ..
                if hdist < minhdist:

                    # Show Update
                    liner.send(
                        '      Sequence {:,}: Region (Start={:,}, End={:,}) has LOWER HDist ({:,})'.format(
                            idx, r, q, hdist))

                    # Minimize r Value
                    r = r - minhdist + hdist

                # HDist was OK!
                elif hdist >= minhdist:
                    condhd = True

                    # Show Update
                    liner.send(
                        '      Sequence {:,}: Region (Start={:,}, End={:,}) is HDist ({:,}) Optimal'.format(
                            idx, r, q, hdist))

            # Both Conditions Met!
            if condol and condtm and condhd:
                idx   += 1     # Move to next sequence
                condtm = False # Reset Tm    Tracer
                condhd = False # Reset HDist Tracer

                # Show Update
                liner.send(
                    '      Sequence {:,}: Region (Start={:,}, End={:,}) is Overlap, Tm and HDist Optimal'.format(
                        idx-1, r, q))

                # Analyze Next Sequence!
                continue

            else:
                # Unresolvable
                if (not condol) or (r < p):
                    liner.send(
                        '\*      Current Split Region Infeasible for Sequence {:,} ... Skipping\n'.format(idx))
                    r = None # No solution to current split
                    break    # Try Next Split Region ..

                # Try again ..
                else:
                    continue

        # Do we have a solution?
        if not r is None:
            liner.send(
                '\*      Current Split Region Optimized for All {:,} Sequences\n'.format(idx))
            state  = True  # Solution Found!
            break # We're done!

    # Return Results
    if not state:
        return None # No Solution
    else:
        return r,q

def aggregate_stats(
    seqlist,
    seqmat,
    split,
    overlap,
    stats,
    liner):
    '''
    Aggregate Melting Temperature, Hamming
    Distance distribution and other metrics
    for split regions. Internal use  only.

    :: seqlist
       type - list
       desc - list of sequences to split
    :: seqmat
       type - np.array
       desc - numeric sequence matrix
    :: split
       type - list
       desc - list of all (start, end)
              split fragment coordinates
    :: overlap
       type - list
       desc - list of all (start, end)
              split overlap cooridnates
    :: stats
       type - dict
       desc - split design stats storage
    :: liner
       type - coroutine
       desc - dynamic printing
    '''

    # Book-keeping
    t0 = tt.time()
    splitstore = cx.deque([] for _ in range(len(split)))
    overlapTmtotals = np.zeros(len(overlap))
    overlapHDtotals = np.zeros(len(overlap))
    TmDenom = np.zeros(len(overlap))
    HDDenom = np.zeros(len(overlap))
    edges   = min(10, max(1, 10**6 / len(seqlist)))
    stats['vars']['splitlens'] = [float('inf') for _ in range(len(split))]

    plen = ut.get_printlen(value=len(overlap))
    qlen = ut.get_printlen(value=len(seqlist))

    # Loop through Strings
    for jdx in range(len(seqlist)):

        # Store Split Strings
        seq = seqlist[jdx]
        for idx in range(len(split)):
            splitseq = seq[split[idx][0]:split[idx][1]]
            stats['vars']['splitlens'][idx] = min(
                len(splitseq),
                stats['vars']['splitlens'][idx])
            if idx % 2 == 1:
                splitseq = ut.get_revcomp(
                    seq=splitseq)
            splitstore[idx].append(splitseq)

        # Loop through Overlaps
        for idx,ol in enumerate(overlap):

            # Compute and Store Tm
            tmelt = ut.get_tmelt(
                seq=seqlist[jdx],
                i=overlap[idx][0],
                j=overlap[idx][1])
            overlapTmtotals[idx] += tmelt
            TmDenom[idx] += 1

            # Compute and Store HDist
            if jdx < len(seqlist)-1:

                # Fetch Subsample Coordinates
                sxr = np.random.randint(
                    low=jdx+1,
                    high=len(seqlist),
                    size=min(edges, len(seqlist)-jdx-1))

                # Compute HDist
                for kdx in sxr:
                    hdist = ut.get_pair_hdist(
                        store=seqmat,
                        idx1=jdx,
                        idx2=kdx,
                        i=overlap[idx][0],
                        j=overlap[idx][1])
                    overlapHDtotals[idx] += hdist
                    HDDenom[idx] += 1

            # Show Update
            liner.send(' Analyzing: Sequence {:{},d} - Overlap {:{},d}'.format(
                jdx+1,
                qlen,
                idx+1,
                plen))

    # Average Statistics
    stats['vars']['numsplits']    = len(split)
    stats['vars']['overlaplens']  = [y-x for x,y in overlap]
    stats['vars']['meanTmdistro'] = list(map(
        int, np.round(overlapTmtotals / TmDenom)))
    stats['vars']['meandistancedistro'] = list(map(
        int, np.round(overlapHDtotals / HDDenom)))

    # Show Time Elapsed
    liner.send('\* Time Elapsed: {:.2f} sec\n'.format(
        tt.time()-t0))

    # Return Results
    return (splitstore,
        stats)
