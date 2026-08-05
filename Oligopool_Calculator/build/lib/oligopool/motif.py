import time  as tt

import collections as cx
import numpy       as np
import atexit      as ae

import nrpcalc     as nr

import utils       as ut
import valparse    as vp
import coremotif   as cm


def motif_engine(
    motifseq,
    homology,
    optrequired,
    leftcontext,
    rightcontext,
    exmotifs,
    exmotifindex,
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
    t0     = tt.time()          # Start Timer
    motifs = [None]*targetcount # Store Motifs
    plen   = ut.get_printlen(   # Target Print Length
        value=targetcount)

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

    # Motif Design Required
    if optrequired:

        # Define Maker Instance
        maker = nr.base.maker.NRPMaker(
            part_type='DNA',
            seed=None)

        # Core Design Loop
        while contextarray:

            # Fetch Context Index
            idx = contextarray.popleft()

            # Fetch Context Sequences
            lcseq =  leftselector(idx)
            rcseq = rightselector(idx)

            # Define Forbidden Prefix and Suffix
            prefixforbidden = prefixdict[lcseq] if not lcseq is None else None
            suffixforbidden = suffixdict[rcseq] if not rcseq is None else None

            # Define Objective Function
            objectivefunction = lambda motif: cm.motif_objectives(
                motif=motif,
                motiflen=len(motifseq),
                exmotifs=exmotifs,
                exmotifindex=exmotifindex,
                lcseq=lcseq,
                rcseq=rcseq,
                edgeeffectlength=edgeeffectlength,
                prefixforbidden=prefixforbidden,
                suffixforbidden=suffixforbidden,
                inittime=t0,
                stats=stats,
                idx=idx+1,
                plen=plen,
                element='Motif',
                liner=liner)

            # Design Motif via Maker
            motif = maker.nrp_maker(
                homology=min(len(motifseq), homology),
                seq_constr=motifseq,
                struct_constr='.'*len(motifseq),
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
            if len(motif) == 0:

                # Terminate Design Loop
                break # RIP .. We failed!

            # A motif was designed!
            else:

                # Extract Motif
                motif = motif[0]

                # Record Designed Motif
                motifs[idx] = motif
                stats['vars']['motifcount'] += 1

                # Show Update
                cm.show_update(
                    idx=idx+1,
                    plen=plen,
                    element='Motif',
                    motif=motifs[idx],
                    optstatus=2,
                    optstate=0,
                    inittime=None,
                    terminal=False,
                    liner=liner)

                cm.extra_assign_motif(
                    motif=motif,
                    contextarray=contextarray,
                    leftselector=leftselector,
                    rightselector=rightselector,
                    edgeeffectlength=edgeeffectlength,
                    prefixdict=prefixdict,
                    suffixdict=suffixdict,
                    storage=motifs,
                    stats=stats,
                    plen=plen,
                    element='Motif',
                    element_key='motifcount',
                    liner=liner)

    # Constant Motif
    else:

        # Constant Solution
        motifs = motifseq
        stats['vars']['motifcount'] = targetcount

    # Check Status and Return Solution
    if not optrequired or \
       stats['vars']['motifcount'] == targetcount:

        # We solved this!
        stats['status'] = True
        stats['basis']  = 'solved'

        # Determine Last Known Motif
        if optrequired:
            lastmotif = motifs[-1]
        else:
            lastmotif = motifseq

        # Final Update
        cm.show_update(
            idx=targetcount,
            plen=plen,
            element='Motif',
            motif=lastmotif,
            optstatus=2,
            optstate=0,
            inittime=t0,
            terminal=True,
            liner=liner)

        # Return Results
        return (motifs, stats)

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

def motif(
    indata,
    oligolimit,
    motifseq,
    motifcol,
    outfile,
    leftcontext=None,
    rightcontext=None,
    exmotifs=None,
    verbose=True):
    '''
    TBW
    '''

    # Start Liner
    liner = ut.liner_engine(verbose)

    # Motif Verbage Print
    liner.send('\n[Oligopool Calculator: Design Mode - Motif]\n')

    # Required Argument Parsing
    liner.send('\n Required Arguments\n')

    # First Pass indata Parsing and Validation
    (indf,
    indata_valid) = vp.get_parsed_indata_info(
        indata=indata,
        indata_field='    Input Data    ',
        required_fields=('ID',),
        precheck=False,
        liner=liner)

    # Full oligolimit Validation
    oligolimit_valid = vp.get_numeric_validity(
        numeric=oligolimit,
        numeric_field='    Oligo Limit   ',
        numeric_pre_desc=' At most ',
        numeric_post_desc=' Base Pair(s)',
        minval=4,
        maxval=float('inf'),
        precheck=False,
        liner=liner)

    # First Pass motifseq Validation
    motifseq_valid = vp.get_seqconstr_validity(
        seqconstr=motifseq,
        seqconstr_field='    Motif Sequence',
        minlenval=1,
        element='MOTIF',
        liner=liner)

    # Full motifcol Validation
    motifcol_valid = vp.get_parsed_column_info(
        col=motifcol,
        df=indf,
        col_field='    Motif Column  ',
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
        outdf_suffix='.oligopool.motif.csv',
        outdf_field='   Output File    ',
        liner=liner)

    # Adjust outfile Suffix
    if not outfile is None:
        outfile = ut.get_adjusted_path(
            path=outfile,
            suffix='.oligopool.motif.csv')

    # Optional Argument Parsing
    liner.send('\n Optional Arguments\n')

    # Store Context Names
    leftcontextname  = leftcontext
    rightcontextname = rightcontext

    # Full leftcontext Parsing and Validation
    (leftcontext,
    leftcontext_valid) = vp.get_parsed_column_info(
        col=leftcontext,
        df=indf,
        col_field='     Left Context ',
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
        col_field='    Right Context ',
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
        exseqs_field=' Excluded Motifs  ',
        exseqs_desc='Unique Motif(s)',
        df_field='Exmotif',
        required=False,
        liner=liner)

    # First Pass Validation
    if not all([
        indata_valid,
        oligolimit_valid,
        motifseq_valid,
        motifcol_valid,
        outfile_valid,
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

    # Motif Design Book-keeping
    targetcount = len(indf.index)
    has_context = False
    outdf = None
    stats = None
    warns = {}

    # Parse Oligopool Limit Feasibility
    liner.send('\n[Step 1: Parsing Oligo Limit]\n')

    # Parse oligolimit
    (parsestatus,
    minvariantlen,
    maxvariantlen,
    minelementlen,
    maxelementlen,
    minspaceavail,
    maxspaceavail) = ut.get_parsed_oligolimit(
        indf=indf,
        variantlens=None,
        oligolimit=oligolimit,
        minelementlen=len(motifseq),
        maxelementlen=len(motifseq),
        element='Motif',
        liner=liner)

    # oligolimit infeasible
    if not parsestatus:

        # Prepare stats
        stats = {
            'status'  : False,
            'basis'   : 'infeasible',
            'step'    : 1,
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

    # Parse Excluded Motifs
    if not exmotifs is None:

        # Show update
        liner.send('\n[Step 2: Parsing Excluded Motifs]\n')

        # Update Step 2 Warning
        warns[2] = {
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
            element='Motif',
            leftcontext=leftcontext,
            rightcontext=rightcontext,
            warn=warns[2],
            liner=liner)

        # Remove Step 2 Warning
        if not warns[2]['warncount']:
            warns.pop(2)

        # exmotifs infeasible
        if not parsestatus:

            # Prepare stats
            stats = {
                'status'  : False,
                'basis'   : 'infeasible',
                'step'    : 2,
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

    # Parsing Sequence Constraint Feasibility
    liner.send('\n[Step 3: Parsing Motif Sequence]\n')

    # Update Step 2 Warning
    warns[3] = {
        'warncount': 0,
        'stepname' : 'parsing-motif-sequence',
        'vars': None}

    # Parse primerseq
    (optrequired,
    homology,
    exmotifindex) = cm.get_parsed_sequence_constraint(
        motifseq=motifseq,
        exmotifs=exmotifs,
        warn=warns[3],
        liner=liner)

    # Remove Step 3 Warning
    if not warns[3]['warncount']:
        warns.pop(3)

    # Parse Edge Effects
    if ((not optrequired is False) and \
        (not exmotifs    is None)) and \
       ((not leftcontext  is None) or \
        (not rightcontext is None)):

        # Set Context Flag
        has_context = True

        # Show Update
        liner.send('\n[Step 4: Extracting Context Sequences]\n')

        # Extract Both Contexts
        (leftcontext,
        rightcontext) = ut.get_extracted_context(
            leftcontext=leftcontext,
            rightcontext=rightcontext,
            edgeeffectlength=edgeeffectlength,
            reduce=False,
            liner=liner)

        # Show update
        liner.send('\n[Step 5: Parsing Edge Effects]\n')

        # Update Step 5 Warning
        warns[5] = {
            'warncount': 0,
            'stepname' : 'parsing-edge-effects',
            'vars': None}

        # Compute Forbidden Prefixes and Suffixes
        (prefixdict,
        suffixdict) = cm.get_parsed_edgeeffects(
            motifseq=motifseq,
            leftcontext=leftcontext,
            rightcontext=rightcontext,
            leftpartition=leftpartition,
            rightpartition=rightpartition,
            exmotifs=exmotifs,
            warn=warns[5],
            liner=liner)

        # Remove Step 5 Warning
        if not warns[5]['warncount']:
            warns.pop(5)

    # Finalize Context
    if not has_context:
        (leftcontext,
        rightcontext,
        prefixdict,
        suffixdict) = (None, None, None, None)

    # Launching Motif Design
    liner.send('\n[Step 6: Computing Motifs]\n')

    # Define Motif Design Stats
    stats = {
        'status'  : False,
        'basis'   : 'unsolved',
        'step'    : 6,
        'stepname': 'computing-motifs',
        'vars'    : {
               'targetcount': targetcount,   # Required Number of Motifs
                'motifcount': 0,             # Motif Design Count
               'exmotiffail': 0,             # Exmotif Elimination Fail Count
                  'edgefail': 0,             # Edge Effect Fail Count
            'exmotifcounter': cx.Counter()}, # Exmotif Encounter Counter
        'warns'   : warns}

    # Schedule outfile deletion
    ofdeletion = ae.register(
        ut.remove_file,
        outfile)

    # Design Motifs
    (motifs,
    stats) = motif_engine(
        motifseq=motifseq,
        homology=homology,
        optrequired=optrequired,
        leftcontext=leftcontext,
        rightcontext=rightcontext,
        exmotifs=exmotifs,
        exmotifindex=exmotifindex,
        edgeeffectlength=edgeeffectlength,
        prefixdict=prefixdict,
        suffixdict=suffixdict,
        targetcount=targetcount,
        stats=stats,
        liner=liner)

    # Motif Status
    if stats['status']:
        motifstatus = 'Successful'
    else:
        motifstatus = 'Failed'

    # Insert motif into indf
    if stats['status']:

        # Update indf
        ut.update_df(
            indf=indf,
            lcname=leftcontextname,
            rcname=rightcontextname,
            out=motifs,
            outcol=motifcol)

        # Prepare outdf
        outdf = indf

        # Write outdf to file
        if not outfile is None:
            outdf.to_csv(
                path_or_buf=outfile,
                sep=',')

    # Motif Design Statistics
    liner.send('\n[Motif Design Statistics]\n')

    plen = ut.get_printlen(
        value=max(stats['vars'][field] for field in (
            'targetcount',
            'motifcount')))

    liner.send(
        '  Design Status   : {}\n'.format(
            motifstatus))
    liner.send(
        '  Target Count    : {:{},d} Motif(s)\n'.format(
            stats['vars']['targetcount'],
            plen))
    liner.send(
        '   Motif Count    : {:{},d} Motif(s) ({:6.2f} %)\n'.format(
            stats['vars']['motifcount'],
            plen,
            ut.safediv(
                A=stats['vars']['motifcount'] * 100.,
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
    if motifstatus == 'Successful':
        ae.unregister(ofdeletion)

    # Close Liner
    liner.close()

    # Return Solution and Statistics
    return (outdf, stats)
