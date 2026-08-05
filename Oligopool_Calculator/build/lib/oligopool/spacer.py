import time  as tt

import collections as cx
import numpy       as np
import atexit      as ae

import nrpcalc     as nr

import utils       as ut
import valparse    as vp
import coremotif   as cm


def spacer_engine(
    spacergroup,
    leftcontext,
    rightcontext,
    exmotifs,
    edgeeffectlength,
    prefixdict,
    suffixdict,
    targetcount,
    stats,
    liner):
    '''
    TBW
    '''

    # Book-keeping
    t0      = tt.time()          # Start Timer
    spacers = [None]*targetcount # Store Spacers
    plen    = ut.get_printlen(   # Target Print Length
        value=targetcount)
    abort   = False              # Abortion Flag

    # Context Setup
    contextarray   = cx.deque(range(targetcount))  # Context Array
    (_,
    leftselector)  = ut.get_context_type_selector( # Left  Context Selector
        context=leftcontext)
    (_,
    rightselector) = ut.get_context_type_selector( # Right Context Selector
        context=rightcontext)

    # Optimize exmotifs
    if not exmotifs is None:
        exmotifs = ut.get_grouped_sequences(
            sequences=exmotifs)

    # Define Maker Instance
    maker = nr.base.maker.NRPMaker(
        part_type='DNA',
        seed=None)

    # Core Design Outer Loop
    while spacergroup:

        # Fetch Context Array
        (spacerlen,
        contextarray) = spacergroup.popitem()

        # Core Design Inner Loop
        while contextarray:

            # Case 1: Zero Length Spacer
            if spacerlen == 0:

                # Assign Gap to Entire Group
                while contextarray:

                    # Fetch Context Index
                    idx = contextarray.popleft()

                    # Assign Gap
                    spacers[idx] = '-'
                    stats['vars']['spacercount'] += 1

                    # Show Update
                    cm.show_update(
                        idx=idx+1,
                        plen=plen,
                        element='Spacer',
                        motif='*GAP*',
                        optstatus=2,
                        optstate=0,
                        inittime=None,
                        terminal=False,
                        liner=liner)

                # Zero Spacers Completed
                break # Jump to Next Group

            # Case 2: Non-Zero Length Spacer
            #         (this requires computation)
            else:

                # Fetch Context Index
                idx = contextarray.popleft()

                # Fetch Context Sequences
                lcseq =  leftselector(idx)
                rcseq = rightselector(idx)

                # Define Forbidden Prefix and Suffix
                prefixforbidden = prefixdict[lcseq] if not lcseq is None else None
                suffixforbidden = suffixdict[rcseq] if not rcseq is None else None

                # Define Objective Function
                objectivefunction = lambda spacer: cm.motif_objectives(
                    motif=spacer,
                    motiflen=spacerlen,
                    exmotifs=exmotifs,
                    exmotifindex=None,
                    lcseq=lcseq,
                    rcseq=rcseq,
                    edgeeffectlength=edgeeffectlength,
                    prefixforbidden=prefixforbidden,
                    suffixforbidden=suffixforbidden,
                    inittime=t0,
                    stats=stats,
                    idx=idx+1,
                    plen=plen,
                    element='Spacer',
                    liner=liner)

                # Design Spacer via Maker
                spacer = maker.nrp_maker(
                    homology=min(spacerlen, 6),
                    seq_constr='N'*spacerlen,
                    struct_constr='.'*spacerlen,
                    target_size=1,
                    background=None,
                    struct_type=None,
                    synth_opt=False,
                    local_model_fn=objectivefunction,
                    jump_count=100,
                    fail_count=100,
                    output_file=None,
                    verbose=False,
                    abortion=True,
                    allow_internal_repeat=True,
                    check_constraints=False)

                # Did we succeed? No ..
                if len(spacer) == 0:

                    # Terminate Design Loop
                    abort = True
                    break # RIP .. We failed!

                # A spacer was designed!
                else:

                    # Extract Spacer
                    spacer = spacer[0]

                    # Record Designed Motif
                    spacers[idx] = spacer
                    stats['vars']['spacercount'] += 1

                    # Show Update
                    cm.show_update(
                        idx=idx+1,
                        plen=plen,
                        element='Spacer',
                        motif=spacers[idx],
                        optstatus=2,
                        optstate=0,
                        inittime=None,
                        terminal=False,
                        liner=liner)

                    cm.extra_assign_motif(
                        motif=spacer,
                        contextarray=contextarray,
                        leftselector=leftselector,
                        rightselector=rightselector,
                        edgeeffectlength=edgeeffectlength,
                        prefixdict=prefixdict,
                        suffixdict=suffixdict,
                        storage=spacers,
                        stats=stats,
                        plen=plen,
                        element='Spacer',
                        element_key='spacercount',
                        liner=liner)

        # Continue Outer Loop?
        if abort:
            break

    # Check Status and Return Solution
    if not abort and \
       stats['vars']['spacercount'] == targetcount:

        # We solved this!
        stats['status'] = True
        stats['basis']  = 'solved'

         # Determine Last Known Motif
        if spacers[-1] != '-':
            lastmotif = spacers[-1]
        else:
            lastmotif = '*GAP*'

        # Final Update
        cm.show_update(
            idx=targetcount,
            plen=plen,
            element='Spacer',
            motif=lastmotif,
            optstatus=2,
            optstate=0,
            inittime=t0,
            terminal=True,
            liner=liner)

        # Return Results
        return (spacers, stats)

    # Design Unsuccessful
    else:

        # This was a miscarriage
        stats['status'] = False
        stats['basis']  = 'unsolved'

        # Final Update
        liner.send('\* Time Elapsed: {:.2f} sec\n'.format(
            tt.time()-t0))

        # Return Results
        return (None, stats)

def spacer(
    indata,
    oligolimit,
    spacercol,
    outfile,
    spacerlen=None,
    leftcontext=None,
    rightcontext=None,
    exmotifs=None,
    verbose=True):
    '''
    TBW
    '''

    # Start Liner
    liner = ut.liner_engine(verbose)

    # Spacer Verbage Print
    liner.send('\n[Oligopool Calculator: Design Mode - Spacer]\n')

    # Required Argument Parsing
    liner.send('\n Required Arguments\n')

    # First Pass indata Parsing and Validation
    (indf,
    indata_valid) = vp.get_parsed_indata_info(
        indata=indata,
        indata_field='    Input Data   ',
        required_fields=('ID',),
        precheck=False,
        liner=liner)

    # Full oligolimit Validation
    oligolimit_valid = vp.get_numeric_validity(
        numeric=oligolimit,
        numeric_field='    Oligo Limit  ',
        numeric_pre_desc=' At most ',
        numeric_post_desc=' Base Pair(s)',
        minval=4,
        maxval=float('inf'),
        precheck=False,
        liner=liner)

    # Full spacercol Validation
    spacercol_valid = vp.get_parsed_column_info(
        col=spacercol,
        df=indf,
        col_field='   Spacer Column ',
        col_desc='Output in Column',
        col_type=1,
        adjcol=None,
        adjval=None,
        iscontext=False,
        typecontext=None,
        liner=liner)

    # Full outfile Validation
    outfile_valid = vp.get_outdf_validity(
        outdf=outfile,
        outdf_suffix='.oligopool.spacer.csv',
        outdf_field='   Output File   ',
        liner=liner)

    # Adjust outfile Suffix
    if not outfile is None:
        outfile = ut.get_adjusted_path(
            path=outfile,
            suffix='.oligopool.spacer.csv')

    # Optional Argument Parsing
    liner.send('\n Optional Arguments\n')

    # Full spacerlen Validation
    (spacerlen,
    spacerlen_valid) = vp.get_parsed_spacerlen_info(
        spacerlen=spacerlen,
        spacerlen_field='   Spacer Length ',
        df_field='Length',
        oligolimit=oligolimit,
        oligolimit_valid=oligolimit_valid,
        indf=indf,
        indata_valid=indata_valid,
        liner=liner)

    # Store Context Names
    leftcontextname  = leftcontext
    rightcontextname = rightcontext

    # Full leftcontext Parsing and Validation
    (leftcontext,
    leftcontext_valid) = vp.get_parsed_column_info(
        col=leftcontext,
        df=indf,
        col_field='     Left Context',
        col_desc='Input from Column',
        col_type=0,
        adjcol=rightcontextname,
        adjval=+1,
        iscontext=True,
        typecontext=0,
        liner=liner)

    # Full leftcontext Parsing and Validation
    (rightcontext,
    rightcontext_valid) = vp.get_parsed_column_info(
        col=rightcontext,
        df=indf,
        col_field='    Right Context',
        col_desc='Input from Column',
        col_type=0,
        adjcol=leftcontextname,
        adjval=-1,
        iscontext=True,
        typecontext=1,
        liner=liner)

    # Full exmotifs Parsing and Validation
    (exmotifs,
    exmotifs_valid) = vp.get_parsed_exseqs_info(
        exseqs=exmotifs,
        exseqs_field=' Excluded Motifs ',
        exseqs_desc='Unique Motif(s)',
        df_field='Exmotif',
        required=False,
        liner=liner)

    # First Pass Validation
    if not all([
        indata_valid,
        oligolimit_valid,
        spacercol_valid,
        outfile_valid,
        spacerlen_valid,
        leftcontext_valid,
        rightcontext_valid,
        exmotifs_valid]):
        liner.send('\n')
        raise RuntimeError(
            'Invalid Argument Input(s).')

    # Start Timer
    t0 = tt.time()

    # Adjust Numeric Paramters
    oligolimit = round(oligolimit)

    # Define Edge Effect Length
    edgeeffectlength = None

    # Spacer Design Book-keeping
    targetcount = len(indf.index)
    has_context = False
    variantlens = None
    outdf = None
    stats = None
    warns = {}

    # Extract spacerlen (Auto-Inference)
    if spacerlen is None:
        liner.send('\n[Step 1: Extracting Spacer Length]\n')
        (spacerlen,
        variantlens) = cm.get_extracted_spacerlen(
            indf=indf,
            oligolimit=oligolimit,
            liner=liner)

    # Parse Oligopool Limit Feasibility
    liner.send('\n[Step 2: Parsing Oligo Limit]\n')

    # Parse oligolimit
    (parsestatus,
    minvariantlen,
    maxvariantlen,
    minelementlen,
    maxelementlen,
    minspaceavail,
    maxspaceavail) = ut.get_parsed_oligolimit(
        indf=indf,
        variantlens=variantlens,
        oligolimit=oligolimit,
        minelementlen=np.min(spacerlen),
        maxelementlen=np.max(spacerlen),
        element='Spacer',
        liner=liner)

    # oligolimit infeasible
    if not parsestatus:

        # Prepare stats
        stats = {
            'status'  : False,
            'basis'   : 'infeasible',
            'step'    : 2,
            'stepname': 'parsing-oligo-limit',
            'vars'    : {
                   'oligolimit': oligolimit,
                'limitoverflow': True,
                'minvariantlen': minvariantlen,
                'maxvariantlen': maxvariantlen,
                'minelementlen': minelementlen,
                'maxelementlen': maxelementlen,
                'minspaceavail': minspaceavail,
                'maxspaceavail': maxspaceavail},
            'warns'   : warns}

        # Return results
        liner.close()
        return (outdf, stats)

    # Extract Spacer Length Groups
    liner.send('\n[Step 3: Extracting Spacer Groups]\n')

    # Group spacerlen
    spacergroup = cm.get_grouped_spacerlen(
        spacerlen=spacerlen,
        liner=liner)

    # Parse Excluded Motifs
    if not exmotifs is None:

        # Show update
        liner.send('\n[Step 4: Parsing Excluded Motifs]\n')

        # Update Step 4 Warning
        warns[4] = {
            'warncount': 0,
            'stepname' : 'parsing-excluded-motifs',
            'vars': None}

        # Parse exmotifs
        (parsestatus,
        exmotifs,
        problens,
        leftpartition,
        rightpartition) = ut.get_parsed_exmotifs(
            exmotifs=exmotifs,
            typer=tuple,
            element='Spacer',
            leftcontext=leftcontext,
            rightcontext=rightcontext,
            warn=warns[4],
            liner=liner)

        # Remove Step 4 Warning
        if not warns[4]['warncount']:
            warns.pop(4)

        # exmotifs infeasible
        if not parsestatus:

            # Prepare stats
            stats = {
                'status'  : False,
                'basis'   : 'infeasible',
                'step'    : 4,
                'stepname': 'parsing-excluded-motifs',
                'vars'    : {
                     'problens': problens,
                    'probcount': tuple(list(
                        4**pl for pl in problens))},
                'warns'   : warns}

            # Return results
            liner.close()
            return (outdf, stats)

        # Update Edge-Effect Length
        edgeeffectlength = ut.get_edgeeffectlength(
            exmotifs=exmotifs)

        # Parse Edge Effects
        if ((not leftcontext  is None) or \
            (not rightcontext is None)):

            # Set Context Flag
            has_context = True

            # Show Update
            liner.send('\n[Step 5: Extracting Context Sequences]\n')

            # Extract Both Contexts
            (leftcontext,
            rightcontext) = ut.get_extracted_context(
                leftcontext=leftcontext,
                rightcontext=rightcontext,
                edgeeffectlength=edgeeffectlength,
                reduce=False,
                liner=liner)

            # Show update
            liner.send('\n[Step 6: Parsing Edge Effects]\n')

            # Update Step 6 Warning
            warns[6] = {
                'warncount': 0,
                'stepname' : 'parsing-edge-effects',
                'vars': None}

            # Compute Forbidden Prefixes and Suffixes
            (prefixdict,
            suffixdict) = cm.get_parsed_edgeeffects(
                motifseq='NNNN',
                leftcontext=leftcontext,
                rightcontext=rightcontext,
                leftpartition=leftpartition,
                rightpartition=rightpartition,
                exmotifs=exmotifs,
                warn=warns[6],
                liner=liner)

            # Remove Step 6 Warning
            if not warns[6]['warncount']:
                warns.pop(6)

    # Finalize Context
    if not has_context:
        (leftcontext,
        rightcontext,
        prefixdict,
        suffixdict) = (None, None, None, None)

    # Launching Motif Design
    liner.send('\n[Step 7: Computing Spacers]\n')

    # Define Motif Design Stats
    stats = {
        'status'  : False,
        'basis'   : 'unsolved',
        'step'    : 7,
        'stepname': 'computing-spacers',
        'vars'    : {
               'targetcount': targetcount,   # Required Number of Spacers
               'spacercount': 0,             # Spacer Design Count
               'exmotiffail': 0,             # Exmotif Elimination Fail Count
                  'edgefail': 0,             # Edge Effect Fail Count
            'exmotifcounter': cx.Counter()}, # Exmotif Encounter Counter
        'warns'   : warns}

    # Schedule outfile deletion
    ofdeletion = ae.register(
        ut.remove_file,
        outfile)

    # Design Spacers
    (spacers,
    stats) = spacer_engine(
        spacergroup=spacergroup,
        leftcontext=leftcontext,
        rightcontext=rightcontext,
        exmotifs=exmotifs,
        edgeeffectlength=edgeeffectlength,
        prefixdict=prefixdict,
        suffixdict=suffixdict,
        targetcount=targetcount,
        stats=stats,
        liner=liner)

    # Spacer Status
    if stats['status']:
        spacerstatus = 'Successful'
    else:
        spacerstatus = 'Failed'

    # Insert spacer into indf
    if stats['status']:

        # Update indf
        ut.update_df(
            indf=indf,
            lcname=leftcontextname,
            rcname=rightcontextname,
            out=spacers,
            outcol=spacercol)

        # Prepare outdf
        outdf = indf

        # Write outdf to file
        if not outfile is None:
            outdf.to_csv(
                path_or_buf=outfile,
                sep=',')

    # Spacer Design Statistics
    liner.send('\n[Spacer Design Statistics]\n')

    plen = ut.get_printlen(
        value=max(stats['vars'][field] for field in (
            'targetcount',
            'spacercount')))

    liner.send(
        '  Design Status   : {}\n'.format(
            spacerstatus))
    liner.send(
        '  Target Count    : {:{},d} Spacer(s)\n'.format(
            stats['vars']['targetcount'],
            plen))
    liner.send(
        '  Spacer Count    : {:{},d} Spacer(s) ({:6.2f} %)\n'.format(
            stats['vars']['spacercount'],
            plen,
            ut.safediv(
                A=stats['vars']['spacercount'] * 100.,
                B=targetcount)))

    # Failure Relavant Stats
    if not stats['status']:
        maxval = max(stats['vars'][field] for field in (
            'exmotiffail',
            'edgefail'))

        sntn, plen = ut.get_notelen(
            printlen=ut.get_printlen(
                value=maxval))

        total_conflicts = stats['vars']['exmotiffail']  + \
                          stats['vars']['edgefail']
        liner.send(
            ' Exmotif Conflicts: {:{},{}} Event(s) ({:6.2f} %)\n'.format(
                stats['vars']['exmotiffail'],
                plen,
                sntn,
                ut.safediv(
                    A=stats['vars']['exmotiffail'] * 100.,
                    B=total_conflicts)))
        liner.send(
            '    Edge Conflicts: {:{},{}} Event(s) ({:6.2f} %)\n'.format(
                stats['vars']['edgefail'],
                plen,
                sntn,
                ut.safediv(
                    A=stats['vars']['edgefail'] * 100.,
                    B=total_conflicts)))

        # Enumerate Motif-wise Fail Counts
        if stats['vars']['exmotifcounter']:

            qlen = max(len(motif) \
                for motif in stats['vars']['exmotifcounter'].keys()) + 2

            sntn, vlen = ut.get_notelen(
                printlen=ut.get_printlen(
                    value=max(
                        stats['vars']['exmotifcounter'].values())))

            liner.send('   Exmotif-wise Conflict Distribution\n')

            for exmotif,count in stats['vars']['exmotifcounter'].most_common():
                exmotif = '\'{}\''.format(exmotif)
                liner.send(
                    '     - Motif {:>{}} Triggered {:{},{}} Event(s)\n'.format(
                        exmotif, qlen, count, vlen, sntn))

    # Show Time Elapsed
    liner.send(
        ' Time Elapsed: {:.2f} sec\n'.format(
            tt.time()-t0))

    # Unschedule outfile deletion
    if spacerstatus == 'Successful':
        ae.unregister(ofdeletion)

    # Close Liner
    liner.close()

    # Return Solution and Statistics
    return (outdf, stats)
