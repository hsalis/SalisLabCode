import time  as tt

import collections as cx
import atexit      as ae

import numpy  as np
import pandas as pd

import utils     as ut
import valparse  as vp
import coresplit as cs


def split_engine(
    seqlist,
    splitlimit,
    mintmelt,
    minhdist,
    maxoverlap,
    minvariantlen,
    maxvariantlen,
    spanlen,
    seqmat,
    varcont,
    stats,
    liner):
    '''
    Compute and return splitting coordinates.
    Internal use only.

    :: seqlist
       type - list
       desc - list of sequences to split

    :: splitlimit
       type - integer
       desc - maximum allowed oligo length
              after splitting
    :: mintmelt
       type - float
       desc - minimum melting temperature of
              split regions
    :: minhdist
       type - integer
       desc - minimum pairwise hamming distance
              between all split regions at a
              given index
    :: maxoverlap
       type - integer
       desc - maximum allowed split overlap
              length
    :: minvariantlen
       type - integer
       desc - length of the shortest variant
              in the oligopool
    :: maxvariantlen
       type - integer
       desc - length of the longest variant
              in the oligopool
    :: spanlen
       type - integer
       desc - minimum required split overlap
              length
    :: seqmat
       type - np.array
       desc - numerically encoded array for
              all variants in oligopool with
              additional padding to make all
              sequences have length equal to
              maxvariantlen
    :: varcont
       type - cx.deque
       desc - a deque of tuple with start and
              end coordinates of variable
              regions
    :: stats
       type - dict
       desc - split design stats storage
    :: liner
       type - coroutine
       desc - dynamic printing
    '''

    # Book-keeping
    split   = []    # Split Coordinate Storage
    overlap = []    # Overlap Coordinate Storage
    fstart  = 0     # Current Fragment  Start
    sstart  = 0     # Current Split     Start
    status  = False # Solution status
    state   = None  # Failure State

    # Compute Split Coordinates
    mt0 = tt.time() # Total Split Time-keeping
    while True:

        # Show Update
        liner.send('\n  Now Computing Split for Fragment: {}\n'.format(
            len(split)+1))
        liner.send('    Initial Fragment Coordinates: (Start={}, End={})\n'.format(
            fstart,
            cs.get_splitend(
                fstart=fstart,
                splitlimit=splitlimit,
                variantlen=maxvariantlen)))

        # Do we need to split any more?
        if not cs.continue_splitting(
            fstart=fstart,
            splitlimit=splitlimit,
            variantlen=maxvariantlen):

            # Store Final Fragment Coordinates
            split.append((fstart, maxvariantlen))

            # Show Updates
            liner.send('    Split Required? No\n')
            liner.send('    Final Fragment Coordinates: (Start={}, End={})\n'.format(
                    *split[-1]))

            # Book-keeping Update
            status = True # Problem Solved!
            break # No more splitting required

        else:

            # Show Updates
            liner.send('    Split Required? Yes\n')

            # Instance Time-keeping
            t0 = tt.time()

            # Get Splitpoints for Current Fragment
            liner.send('    Finding Splittable Regions ...')
            spq = cs.get_splitqueue(
                varcont=varcont,
                fstart=fstart,
                sstart=sstart,
                splitlimit=splitlimit,
                variantlen=maxvariantlen)

            # Did we find split regions?
            if not spq: # No Split Regions Found
                liner.send('    No Splittable Regions Found ... Terminating\n')
                status = False # No Solution
                state  = 0
                break
            else:       # Split Regions Found
                liner.send('    Splittable Regions Found: {} (in {:.2f} sec)\n'.format(
                    len(spq),
                    tt.time() - t0))

            # Get the Tm and HDist based split
            rq = cs.get_split(
                seqlist=seqlist,
                seqmat=seqmat,
                spq=spq,
                mintmelt=mintmelt,
                minhdist=minhdist,
                spanlen=spanlen,
                maxoverlap=maxoverlap,
                liner=liner)

            # Did we find feasible split regions?
            if rq is None: # No Feasible Split Found
                liner.send('    No Feasible Splits Found ... Terminating\n')
                status = False # No Solution
                state  = 0
                break

            else:          # Feasigle Split Found
                r,q = rq   # Parse Split
                # Store Current Fragment Coordinates
                split.append((fstart, q))
                overlap.append((r, q))

                # Book-keeping Update
                fstart = r # Next Fragment Start Coordinate
                sstart = q # Next Split    Start Coordinate

                # Show Updates
                liner.send('    Split Region Selected: (Start={}, End={}) (in {:.2f} sec)\n'.format(
                    r, q, tt.time()-t0))
                liner.send('    Final Fragment {} Coordinates: (Start={}, End={})\n'.format(
                    len(split), *split[-1]))

    # Did the shorter variants get split as well?
    if status:
        for u,v in overlap:
            if u > minvariantlen:
                status = False
                state  = 1
                break

    # Compute Verdict
    if status:
        # Update Stats
        stats['status'] = True
        stats['basis']  = 'solved'
        # Sh
        liner.send(
            '\n  Solution Status: Splitting Completed\n'.format())
    else:
        if state == 0:
            # Update Stats
            stats['vars']['infeasiblecontigs'] = True

            liner.send(
                '\n  Solution Status: Unsolved due to Infeasible Variable Contigs\n')
        else:
            # Update Stats
            stats['vars']['unevensplits'] = True

            liner.send(
                '\n  Solution Status: Unsolved due to Uneven Number of Splits\n')

    # Final Updates
    liner.send('  Time Elapsed: {:.2f} sec\n'.format(tt.time() - mt0))

    # Return Results
    return (split,
        overlap,
        stats)

def split(
    indata,
    splitlimit,
    mintmelt,
    minhdist,
    minoverlap,
    maxoverlap,
    outfile=None,
    verbose=True):
    '''
    TBW
    '''

    # Start Liner
    liner = ut.liner_engine(verbose)

    # Splitting Verbage Print
    liner.send('\n[Oligopool Calculator: Design Mode - Split]\n')

    # Required Argument Parsing
    liner.send('\n Required Arguments\n')

    # First Pass indata Parsing and Validation
    (indf,
    indata_valid) = vp.get_parsed_indata_info(
        indata=indata,
        indata_field='   Input Data       ',
        required_fields=('ID',),
        precheck=False,
        liner=liner)

    # Full splitlimit Validation
    splitlimit_valid = vp.get_numeric_validity(
        numeric=splitlimit,
        numeric_field='   Split Limit      ',
        numeric_pre_desc=' Split Fragments at most ',
        numeric_post_desc=' Base Pair(s) Each',
        minval=4,
        maxval=float('inf'),
        precheck=False,
        liner=liner)

    # Full mintmelt and maxtmelt Validation
    tmelt_valid = vp.get_numeric_validity(
        numeric=mintmelt,
        numeric_field=' Melting Temperature',
        numeric_pre_desc=' At least ',
        numeric_post_desc=' °C b/w on-target Overlaps',
        minval=4,
        maxval=float('inf'),
        precheck=False,
        liner=liner)

    # Full minhdist Validation
    minhdist_valid = vp.get_numeric_validity(
        numeric=minhdist,
        numeric_field=' Hamming Distance   ',
        numeric_pre_desc=' At least ',
        numeric_post_desc=' Mismatch(es) per Off-Target Overlap Pair ',
        minval=1,
        maxval=splitlimit if splitlimit_valid else float('inf'),
        precheck=False,
        liner=liner)

    # Full minoverlap and maxoverlap Validation
    (minoverlap,
    maxoverlap,
    overlap_valid) = vp.get_parsed_range_info(
        minval=minoverlap,
        maxval=maxoverlap,
        range_field=' Overlap Length     ',
        range_unit='Base Pair(s) Fragment Overlap(s)',
        range_min=15 if not minhdist_valid else max(15, minhdist),
        range_max=splitlimit if splitlimit_valid else float('inf'),
        liner=liner)

    # Full outfile Validation
    outfile_valid = vp.get_outdf_validity(
        outdf=outfile,
        outdf_suffix='.oligopool.split.csv',
        outdf_field='  Output File       ',
        liner=liner)

    # Adjust outfile Suffix
    if not outfile is None:
        outfile = ut.get_adjusted_path(
            path=outfile,
            suffix='.oligopool.split.csv')

    # First Pass Validation
    if not all([
        indata_valid,
        splitlimit_valid,
        tmelt_valid,
        minhdist_valid,
        overlap_valid,
        outfile_valid]):
        liner.send('\n')
        raise RuntimeError(
            'Invalid Argument Input(s).')

    # Start Timer
    t0 = tt.time()

    # Adjust Numeric Paramters
    splitlimit = round(splitlimit)
    minhdist   = round(minhdist)
    minoverlap = round(minoverlap)
    maxoverlap = round(maxoverlap)

    # Primer Design Book-keeping
    outdf = None
    stats = None
    warns = {}

    # Parse Oligopool Split Limit Feasibility
    liner.send('\n[Step 1: Parsing Split Limit]\n')

    # Parse splitlimit
    (parsestatus,
    seqlist,
    oligounderflow,
    unevensplit,
    minvariantlen,
    maxvariantlen,
    minsplitcount,
    maxsplitcount) = cs.get_parsed_splitlimit(
        indf=indf,
        splitlimit=splitlimit,
        liner=liner)

    # splitlimit infeasible
    if not parsestatus:

        # Prepare stats
        stats = {
            'status'  : False,
            'basis'   : 'infeasible',
            'step'    : 1,
            'stepname': 'parsing-split-limit',
            'vars'    : {
                    'splitlimit': splitlimit,
                'oligounderflow': oligounderflow,
                   'unevensplit': unevensplit,
                 'minvariantlen': minvariantlen,
                 'maxvariantlen': maxvariantlen,
                 'minsplitcount': minsplitcount,
                 'maxsplitcount': maxsplitcount,},
            'warns'   : warns}

        # Return results
        liner.close()
        return (outdf, stats)

    # Define spanlen
    # Note: spanlen is the minimum span
    #       for a split region
    spanlen = max(minhdist, minoverlap)

    # Compute Sequence Matrix
    liner.send('\n[Step 2: Computing Sequence Matrix]\n')

    # Compute padvec and seqmat
    (padvec,
    seqmat) = cs.get_seqmat_padvec(
        seqlist=seqlist,
        maxvariantlen=maxvariantlen,
        liner=liner)

    # Compute Sequence Matrix
    liner.send('\n[Step 3: Computing Entropy Vector]\n')

    # Compute padvec and seqmat
    entvec = cs.get_entvec(
        seqmat=seqmat,
        maxvariantlen=maxvariantlen,
        liner=liner)

    # Parse Oligopool Limit Feasibility
    liner.send('\n[Step 4: Parsing Variable Contigs]\n')

    # Parse splitlimit
    (parsestatus,
    varcont,
    varcontcount,
    mergedcontcount,
    filtercontcount) = cs.get_varcont(
        entvec=entvec,
        minhdist=minhdist,
        spanlen=spanlen,
        liner=liner)

    # splitlimit infeasible
    if not parsestatus:

        # Prepare stats
        stats = {
            'status'  : False,
            'basis'   : 'infeasible',
            'step'    : 4,
            'stepname': 'parsing-variable-contig',
            'vars'    : {
                   'varcontcount': varcontcount,
                'mergedcontcount': mergedcontcount,
                'filtercontcount': filtercontcount},
            'warns'   : warns}

        # Return results
        liner.close()
        return (outdf, stats)

    # Launching Split Design
    liner.send('\n[Step 5: Computing Split]\n')

    # Define Split Design Stats
    stats = {
        'status'  : False,
        'basis'   : 'unsolved',
        'step'    : 8,
        'stepname': 'computing-split',
        'vars'    : {
                     'numsplits': 0,   # Total Number of Splits
                     'splitlens': [],  # Split Oligo Lengths
                   'overlaplens': [],  # Split Overlap Lengths
                  'meanTmdistro': [],  # Mean Tm of Each Split
            'meandistancedistro': [],  # Mean HDist of Each Split
             'infeasiblecontigs': False,  # Infeasible Contigs Flag
                  'unevensplits': False}, # Uneven Splits Flag
        'warns'   : warns}

    # Schedule outfile deletion
    ofdeletion = ae.register(
        ut.remove_file,
        outfile)

    # Design Split
    (split,
    overlap,
    stats) = split_engine(
        seqlist=seqlist,
        splitlimit=splitlimit,
        mintmelt=mintmelt,
        minhdist=minhdist,
        maxoverlap=maxoverlap,
        minvariantlen=minvariantlen,
        maxvariantlen=maxvariantlen,
        spanlen=spanlen,
        seqmat=seqmat,
        varcont=varcont,
        stats=stats,
        liner=liner)

    # Success Relevant Stats
    if stats['status']:

        # Launching Stats Aggregation
        liner.send('\n[Step 6: Aggregating Stats]\n')

        # Compute Tm and HDist Distribution
        # and Finalize Split Sequences
        (splitstore,
        stats) = cs.aggregate_stats(
            seqlist=seqlist,
            seqmat=seqmat,
            split=split,
            overlap=overlap,
            stats=stats,
            liner=liner)

    # Split Status
    if stats['status']:
        splitstatus = 'Successful'
    else:
        splitstatus = 'Failed'

    # Insert split into outdf
    if stats['status']:

        # Prepare outdf
        outdf = pd.DataFrame(
            index=indf.index)

        # Insert Splits into Columns
        for sidx in range(len(split)):
            outdf['Split-{}'.format(sidx+1)] = splitstore[sidx]

        # Write outdf to file
        if not outfile is None:
            outdf.to_csv(
                path_or_buf=outfile,
                sep=',')

    # Split Design Statistics
    liner.send('\n[Split Design Statistics]\n')

    liner.send(
        '     Design Status  : {}\n'.format(
            splitstatus))

    # Success Relevant Stats
    if stats['status']:

        maxval = max(max(min(stats['vars'][field]) for field in (
            'splitlens',
            'overlaplens',
            'meanTmdistro',
            'meandistancedistro')),
            stats['vars']['numsplits'])

        sntn, plen = ut.get_notelen(
            printlen=ut.get_printlen(
                value=maxval))

        liner.send(
            '     No. of Splits  : {:{},d} Fragments per Variant\n'.format(
                stats['vars']['numsplits'],
                plen))
        liner.send(
            '      Split Length  : {:{},d} {}Base Pair(s)\n'.format(
                min(stats['vars']['splitlens']),
                plen,
                ['', 'to {:,d} '.format(max(stats['vars']['splitlens']))][
                    max(stats['vars']['splitlens']) != min(stats['vars']['splitlens'])
                ]))
        liner.send(
            '    Overlap Length  : {:{},d} {}Base Pair(s)\n'.format(
                min(stats['vars']['overlaplens']),
                plen,
                ['', 'to {:,d} '.format(max(stats['vars']['overlaplens']))][
                    max(stats['vars']['overlaplens']) != min(stats['vars']['overlaplens'])
                ]))
        liner.send(
            '    Overlap Tm      : {:{},d} {}°C\n'.format(
                min(stats['vars']['meanTmdistro']),
                plen,
                ['', 'to {:,d} '.format(max(stats['vars']['meanTmdistro']))][
                    max(stats['vars']['meanTmdistro']) != min(stats['vars']['meanTmdistro'])
                ]))
        liner.send(
            '    Overlap Distance: {:{},d} {}Mismatch(es)\n'.format(
               min(stats['vars']['meandistancedistro']),
                plen,
                ['', 'to {:,d} '.format(max(stats['vars']['meandistancedistro']))][
                    max(stats['vars']['meandistancedistro']) != min(stats['vars']['meandistancedistro'])
                ]))

    # Failure Relavant Stats
    else:
        liner.send(
            ' Infeasible Contigs : {}\n'.format(
                ['No', 'Yes'][stats['vars']['infeasiblecontigs']]))
        liner.send(
            '      Unven Split   : {}\n'.format(
                ['No', 'Yes'][stats['vars']['unevensplits']]))

    # Show Time Elapsed
    liner.send(
        ' Time Elapsed: {:.2f} sec\n'.format(
            tt.time()-t0))

    # Unschedule outfile deletion
    if splitstatus == 'Successful':
        ae.unregister(ofdeletion)

    # Close Liner
    liner.close()

    # Return Solution and Statistics
    return (outdf, stats)
