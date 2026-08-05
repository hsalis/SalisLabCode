import time  as tt

import collections as cx
import atexit      as ae

import nrpcalc     as nr

import utils       as ut
import valparse    as vp
import coreprimer  as cp
import primer      as pr


def pad(
    indata,
    splitcol,
    typeIIS,
    oligolimit,
    mintmelt,
    maxtmelt,
    maxreplen,
    outfile=None,
    verbose=True):
    '''

    Note 2. 34 unique TypeIIS systems available for padding
            AcuI,  AlwI,  BbsI,  BccI,   BceAI,    BciVI,
            BcoDI, BmrI,  BpuEI, BsaI,   BseRI,    BsmAI,
            BsmBI, BsmFI, BsmI,  BspCNI, BspQI,    BsrDI,
            BsrI,  BtgZI, BtsCI, BtsI,   BtsIMutI, EarI,
            EciI,  Esp3I, FauI,  HgaI,   HphI,     HpyAV,
            MlyI,  MnlI,  SapI,  SfaNI
    '''

    # Start Liner
    liner = ut.liner_engine(verbose)

    # Barcoding Verbage Print
    liner.send('\n[Oligopool Calculator: Design Mode - Pad]\n')

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

    # Store Split Column Name
    splitcolname = splitcol

    # Full splitcol Validation
    (splitcol,
    splitcol_valid) = vp.get_parsed_column_info(
        col=splitcol,
        df=indf,
        col_field='   Split Column     ',
        col_desc='Input from Column',
        col_type=0,
        adjcol=None,
        adjval=None,
        iscontext=False,
        typecontext=None,
        liner=liner)

    # Full typeIIS Validation
    (typeIIS,
    typeIISname,
    typeIIS_valid) = vp.get_parsed_typeIIS_info(
        typeIIS=typeIIS,
        typeIIS_field=' TypeIIS System     ',
        liner=liner)

    # Full maxreplen Validation
    oligolimit_valid = vp.get_numeric_validity(
        numeric=oligolimit,
        numeric_field='   Oligo Limit      ',
        numeric_pre_desc=' Design ',
        numeric_post_desc=' Base Pair(s) Padded Oligos',
        minval=60,
        maxval=float('+inf'),
        precheck=False,
        liner=liner)

    # Full mintmelt and maxtmelt Validation
    (mintmelt,
    maxtmelt,
    tmelt_valid) = vp.get_parsed_range_info(
        minval=mintmelt,
        maxval=maxtmelt,
        range_field=' Melting Temperature',
        range_unit='°C',
        range_min=25,
        range_max=95,
        liner=liner)

    # Full maxreplen Validation
    maxreplen_valid = vp.get_numeric_validity(
        numeric=maxreplen,
        numeric_field='  Repeat Length     ',
        numeric_pre_desc=' Up to ',
        numeric_post_desc=' Base Pair(s) Oligopool Repeats',
        minval=6,
        maxval=20,
        precheck=False,
        liner=liner)

    # Full outfile Validation
    outfile_valid = vp.get_outdf_validity(
        outdf=outfile,
        outdf_suffix='.oligopool.pad.csv',
        outdf_field='  Output File       ',
        liner=liner)

    # Adjust outfile Suffix
    if not outfile is None:
        outfile = ut.get_adjusted_path(
            path=outfile,
            suffix='.oligopool.pad.csv')

    # First Pass Validation
    if not all([
        indata_valid,
        splitcol_valid,
        typeIIS_valid,
        oligolimit_valid,
        tmelt_valid,
        maxreplen_valid,
        outfile_valid]):
        liner.send('\n')
        raise RuntimeError(
            'Invalid Argument Input(s).')

    # Start Timer
    t0 = tt.time()

    # Adjust Numeric Paramters
    oligolimit = round(oligolimit)
    maxreplen  = round(maxreplen)

    # Define Additional Variables
    typeIISmotif = typeIIS.replace('N', '')
    homology     = len(typeIISmotif) + 1
    background   = None

    # Primer Design Book-keeping
    outdf = None
    stats = None
    warns = {}

    # Parse Split Column
    liner.send('\n[Step 1: Parsing Split Column]\n')

    # Parse splitcol
    (parsestatus,
    minfragmentlen,
    maxfragmentlen,
    maxallowedlen,
    paddingbalance) = cp.get_parsed_splitcol(
        splitcol=splitcol,
        oligolimit=oligolimit,
        liner=liner)

    # splitcol infeasible
    if not parsestatus:

        # Prepare stats
        stats = {
            'status'  : False,
            'basis'   : 'infeasible',
            'step'    : 1,
            'stepname': 'parsing-split-column',
            'vars'    : {
                'maxfragmentlen': maxfragmentlen,
                 'maxallowedlen': maxallowedlen},
            'warns'   : warns}

        # Return results
        liner.close()
        return (outdf, stats)

    # Parse TypeIIS Constraint
    liner.send('\n[Step 2: Parsing TypeIIS System]\n')

    # Parse typeIIS
    (parsestatus,
    fwdcore,
    revcore,
    fwdseq,
    revseq,
    minpadlen,
    maxpadlen,
    typeIISfree) = cp.get_parsed_typeIIS_constraint(
        typeIIS=typeIIS,
        typeIISname=typeIISname,
        minfragmentlen=minfragmentlen,
        maxfragmentlen=maxfragmentlen,
        oligolimit=oligolimit,
        liner=liner)

    # typeIIS infeasible
    if not parsestatus:

        # Prepare stats
        stats = {
            'status'  : False,
            'basis'   : 'infeasible',
            'step'    : 2,
            'stepname': 'parsing-typeIIS-system',
            'vars'    : {
                  'minpadlen': minpadlen,
                  'maxpadlen': maxpadlen,
                'typeIISfree': typeIISfree},
            'warns'   : warns}

        # Return results
        liner.close()
        return (outdf, stats)

    # Parse Melting Temperature
    liner.send('\n[Step 3: Parsing Melting Temperature]\n')

    # Parse mintmelt and maxtmelt
    (parsestatus,
    estimatedminTm,
    estimatedmaxTm,
    higherminTm,
    lowermaxTm,
    mintmelt,
    maxtmelt) = cp.get_parsed_primer_tmelt_constraint(
        primerseq=revseq[:revcore],
        pairedprimer=None,
        mintmelt=mintmelt,
        maxtmelt=maxtmelt,
        element='Pad',
        liner=liner)

    # mintmelt and maxtmelt infeasible
    if not parsestatus:

        # Prepare stats
        stats = {
            'status'  : False,
            'basis'   : 'infeasible',
            'step'    : 3,
            'stepname': 'parsing-melting-temperature',
            'vars'    : {
                'estimatedminTm': estimatedminTm,
                'estimatedmaxTm': estimatedmaxTm,
                   'higherminTm': higherminTm,
                    'lowermaxTm': lowermaxTm},
            'warns'   : warns}

        # Return results
        liner.close()
        return (outdf, stats)

    # Parse Excluded Motifs
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
        exmotifs=(typeIISmotif,ut.get_revcomp(typeIISmotif)),
        typer=tuple,
        element='Pad',
        leftcontext=splitcol,
        rightcontext=splitcol,
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

    # Show update
    liner.send('\n[Step 5: Extracting Context Sequences]\n')

    # Extract Pad Contexts
    (leftcontext,
    rightcontext) = ut.get_extracted_context(
        leftcontext=splitcol,
        rightcontext=splitcol,
        edgeeffectlength=edgeeffectlength,
        reduce=True,
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
    suffixdict) = cp.get_parsed_edgeeffects(
        primerseq=fwdseq[-fwdcore:],
        leftcontext=leftcontext,
        rightcontext=rightcontext,
        leftpartition=leftpartition,
        rightpartition=rightpartition,
        exmotifs=exmotifs,
        element='Pad',
        warn=warns[6],
        liner=liner)

    # Remove Step 6 Warning
    if not warns[6]['warncount']:
        warns.pop(6)

    # Parse Oligopool Repeats
    liner.send('\n[Step 7: Parsing Oligopool Repeats]\n')

    # Parse Repeats from indf
    (parsestatus,
    sourcecontext,
    kmerspace,
    fillcount,
    freecount,
    oligorepeats) = ut.get_parsed_oligopool_repeats(
        df=indf,
        maxreplen=maxreplen,
        element='Pad',
        merge=True,
        liner=liner)

    # Repeat Length infeasible
    if not parsestatus:

        # Prepare stats
        stats = {
            'status'  : False,
            'basis'   : 'infeasible',
            'step'    : 7,
            'stepname': 'parsing-oligopool-repeats',
            'vars'    : {
                'sourcecontext': sourcecontext,
                'kmerspace'    : kmerspace,
                'fillcount'    : fillcount,
                'freecount'    : freecount},
            'warns'   : warns}

        # Return results
        liner.close()
        return (outdf, stats)

    # Define Pad Design Stats
    stats = {
        'status'  : False,
        'basis'   : 'unsolved',
        'step'    : 8,
        'stepname': 'computing-primer',
        'vars'    : {
                'fwdpadprimerTm': None,          # Forward Pad Melting Temperature
                'revpadprimerTm': None,          # Reverse Pad Melting Temperature
                'fwdpadprimerGC': None,          # Forward Pad GC Content
                'revpadprimerGC': None,          # Reverse Pad GC Content
              'fwdpadhairpinMFE': None,          # Forward Pad Hairpin Free Energy
              'revpadhairpinMFE': None,          # Reverse Pad Hairpin Free Energy
            'fwdpadhomodimerMFE': None,          # Forward Pad Homodimer Free Energy
            'revpadhomodimerMFE': None,          # Reverse Pad Homodimer Free Energy
                'heterodimerMFE': None,          # Heterodimer Free Energy
                        'Tmfail': 0,             # Melting Temperature Fail Count
                    'repeatfail': 0,             # Repeat Fail Count
                 'homodimerfail': 0,             # Homodimer Fail Count
               'heterodimerfail': 0,             # Heterodimer Fail Count
                   'exmotiffail': 0,             # Exmotif Elimination Fail Count
                      'edgefail': 0,             # Edge Effect Fail Count
                'exmotifcounter': cx.Counter()}, # Exmotif Encounter Counter
        'warns'   : warns}

    # Schedule outfile deletion
    ofdeletion = ae.register(
        ut.remove_file,
        outfile)

    # Define Forward Primer-Pad Design Stats
    fwdstats = {
        'status'  : False,
        'basis'   : 'unsolved',
        'step'    : 8,
        'stepname': 'computing-forward-pad',
        'vars'    : {
                   'primerTm': None,          # Primer Melting Temperature
                   'primerGC': None,          # Primer GC Content
                 'hairpinMFE': None,          # Primer Hairpin Free Energy
               'homodimerMFE': None,          # Homodimer Free Energy
             'heterodimerMFE': None,          # Heterodimer Free Energy
                     'Tmfail': 0,             # Melting Temperature Fail Count
                 'repeatfail': 0,             # Repeat Fail Count
              'homodimerfail': 0,             # Homodimer Fail Count
            'heterodimerfail': 0,             # Heterodimer Fail Count
                'exmotiffail': 0,             # Exmotif Elimination Fail Count
                   'edgefail': 0,             # Edge Effect Fail Count
             'exmotifcounter': cx.Counter()}, # Exmotif Encounter Counter
        'warns'   : warns}

    # Define Forward Primer-Pad Attributes
    pairedrepeats = set()
    exmotifindex  = set([fwdseq.index(
        typeIISmotif) + len(typeIISmotif)])

    # Launching Forward Primer-Pad Design
    liner.send('\n[Step 8: Computing Forward Pad]\n')

    # Design Forward Primer-Pad
    (fwdpad,
    fwdstats) = pr.primer_engine(
        primerseq=fwdseq,
        primerspan=fwdcore,
        homology=homology,
        primertype=0,
        fixedbaseindex=cp.get_fixedbaseindex(seq=fwdseq),
        mintmelt=mintmelt,
        maxtmelt=maxtmelt,
        maxreplen=maxreplen,
        oligorepeats=oligorepeats,
        pairedprimer=None,
        pairedspan=None,
        pairedrepeats=pairedrepeats,
        exmotifs=exmotifs,
        exmotifindex=exmotifindex,
        edgeeffectlength=edgeeffectlength,
        prefixdict=None,
        suffixdict=suffixdict,
        background=background,
        stats=fwdstats,
        liner=liner)

    # Define Reverse Primer-Pad Design Stats
    revstats = {
        'status'  : False,
        'basis'   : 'unsolved',
        'step'    : 9,
        'stepname': 'computing-reverse-pad',
        'vars'    : {
                   'primerTm': None,          # Primer Melting Temperature
                   'primerGC': None,          # Primer GC Content
                 'hairpinMFE': None,          # Primer Hairpin Free Energy
               'homodimerMFE': None,          # Homodimer Free Energy
             'heterodimerMFE': None,          # Heterodimer Free Energy
                     'Tmfail': 0,             # Melting Temperature Fail Count
                 'repeatfail': 0,             # Repeat Fail Count
              'homodimerfail': 0,             # Homodimer Fail Count
            'heterodimerfail': 0,             # Heterodimer Fail Count
                'exmotiffail': 0,             # Exmotif Elimination Fail Count
                   'edgefail': 0,             # Edge Effect Fail Count
             'exmotifcounter': cx.Counter()}, # Exmotif Encounter Counter
        'warns'   : warns}

    # Do we Continue?
    if fwdstats['status']:

        # Define Reverse Primer-Pad Attributes
        pairedrepeats = set(ut.stream_canon_spectrum(
            seq=fwdpad[-fwdcore:],
            k=len(typeIIS)))
        exmotifindex  = set([revseq.index(
            ut.get_revcomp(typeIISmotif)) + len(typeIISmotif)])

        # Launching Reverse Primer-Pad Design
        liner.send('\n[Step 9: Computing Reverse Pad]\n')

        # Design Reverse Primer-Pad
        (revpad,
        revstats) = pr.primer_engine(
            primerseq=revseq,
            primerspan=revcore,
            homology=homology,
            primertype=1,
            fixedbaseindex=cp.get_fixedbaseindex(seq=revseq),
            mintmelt=fwdstats['vars']['primerTm']-1,
            maxtmelt=fwdstats['vars']['primerTm']+1,
            maxreplen=maxreplen,
            oligorepeats=oligorepeats,
            pairedprimer=fwdpad,
            pairedspan=fwdcore,
            pairedrepeats=pairedrepeats,
            exmotifs=exmotifs,
            exmotifindex=exmotifindex,
            edgeeffectlength=edgeeffectlength,
            prefixdict=prefixdict,
            suffixdict=None,
            background=background,
            stats=revstats,
            liner=liner)

    # Meta Merge
    stats['status']   = fwdstats['status'] and \
                        revstats['status']
    stats['basis']    = 'solved' if stats['status'] else 'unsolved'
    stats['step']     = fwdstats['step'] if not revstats['status'] \
                                         else   revstats['step']
    stats['stepname'] = fwdstats['stepname'] if not revstats['status'] \
                                             else   revstats['stepname']

    # Forward Stats Merge
    stats['vars']['fwdpadprimerTm']     = fwdstats['vars']['primerTm']
    stats['vars']['fwdpadprimerGC']     = fwdstats['vars']['primerGC']
    stats['vars']['fwdpadhairpinMFE']   = fwdstats['vars']['hairpinMFE']
    stats['vars']['fwdpadhomodimerMFE'] = fwdstats['vars']['homodimerMFE']
    stats['vars']['Tmfail']             = fwdstats['vars']['Tmfail']
    stats['vars']['repeatfail']         = fwdstats['vars']['repeatfail']
    stats['vars']['homodimerfail']      = fwdstats['vars']['homodimerfail']
    stats['vars']['heterodimerfail']    = fwdstats['vars']['heterodimerfail']
    stats['vars']['exmotiffail']        = fwdstats['vars']['exmotiffail']
    stats['vars']['edgefail']           = fwdstats['vars']['edgefail']
    stats['vars']['exmotifcounter']     = fwdstats['vars']['exmotifcounter']

    # Reverse Stats Merge
    stats['vars']['revpadprimerTm']     = revstats['vars']['primerTm']
    stats['vars']['revpadprimerGC']     = revstats['vars']['primerGC']
    stats['vars']['revpadhairpinMFE']   = revstats['vars']['hairpinMFE']
    stats['vars']['revpadhomodimerMFE'] = revstats['vars']['homodimerMFE']
    stats['vars']['heterodimerMFE']     = revstats['vars']['heterodimerMFE']
    stats['vars']['Tmfail']            += revstats['vars']['Tmfail']
    stats['vars']['repeatfail']        += revstats['vars']['repeatfail']
    stats['vars']['homodimerfail']     += revstats['vars']['homodimerfail']
    stats['vars']['heterodimerfail']   += revstats['vars']['heterodimerfail']
    stats['vars']['exmotiffail']       += revstats['vars']['exmotiffail']
    stats['vars']['edgefail']          += revstats['vars']['edgefail']
    stats['vars']['exmotifcounter']    += revstats['vars']['exmotifcounter']

    # Primer Status
    if stats['status']:
        padstatus = 'Successful'
    else:
        padstatus = 'Failed'

    # Insert primer into indf
    if stats['status']:

        # Extract indf
        indf = indf[[splitcolname]]

        # Compute columns
        LeftSpacer    = []
        ForwardPrimer = []
        RightSpacer   = []
        ReversePrimer = []

        # Decompose Balance
        for balance in paddingbalance:

            # Get the Current Balance
            p,q = balance

            # Left Pad Extration from Balance
            xfwdpad = fwdpad[-p:]
            s = p - fwdcore

            # Do we have Left Spacer?
            if s > 0:
                leftspacer = xfwdpad[:s]
            else:
                leftspacer = '-'
            fwdprimer = xfwdpad[-fwdcore:]

            # Right Pad Extration from Balance
            xrevpad = revpad[:q]
            s = q - revcore

            # Do we have Right Spacer?
            if s > 0:
                rightspacer = xrevpad[-s:]
            else:
                rightspacer = '-'
            revprimer = xrevpad[:+revcore]

            # Add Elements to Columns
            LeftSpacer.append(leftspacer)
            RightSpacer.append(rightspacer)
            ForwardPrimer.append(fwdprimer)
            ReversePrimer.append(revprimer)

        # Add columns
        indf['5\'Spacer']     = LeftSpacer
        indf['3\'Spacer']     = RightSpacer
        indf['ForwardPrimer'] = ForwardPrimer
        indf['ReversePrimer'] = ReversePrimer

        # Prepare outdf
        outdf = indf
        outdf = outdf[['5\'Spacer',
            'ForwardPrimer',
            splitcolname,
            'ReversePrimer',
            '3\'Spacer']]

        # Write outdf to file
        if not outfile is None:
            outdf.to_csv(
                path_or_buf=outfile,
                sep=',')

    # Primer Design Statistics
    liner.send('\n[Pad Design Statistics]\n')

    liner.send(
        '       Design Status    : {}\n'.format(
            padstatus))

    # Success Relevant Stats
    if stats['status']:

        liner.send(
            '     Melting Temperature: {:6.2f} °C and {:6.2f} °C\n'.format(
                stats['vars']['fwdpadprimerTm'],
                stats['vars']['revpadprimerTm']))
        liner.send(
            '          GC Content    : {:6.2f} %  and {:6.2f} %\n'.format(
                stats['vars']['fwdpadprimerGC'],
                stats['vars']['revpadprimerGC']))
        liner.send(
            '     Hairpin MFE        : {:6.2f} kcal/mol and {:6.2f} kcal/mol\n'.format(
                stats['vars']['fwdpadhairpinMFE'],
                stats['vars']['revpadhairpinMFE']))
        liner.send(
            '   Homodimer MFE        : {:6.2f} kcal/mol and {:6.2f} kcal/mol\n'.format(
                stats['vars']['fwdpadhomodimerMFE'],
                stats['vars']['revpadhomodimerMFE']))
        liner.send(
            ' Heterodimer MFE        : {:6.2f} kcal/mol\n'.format(
                stats['vars']['heterodimerMFE']))

    # Failure Relavant Stats
    else:
        maxval = max(stats['vars'][field] for field in (
            'Tmfail',
            'repeatfail',
            'homodimerfail',
            'heterodimerfail',
            'exmotiffail',
            'edgefail'))

        sntn, plen = ut.get_notelen(
            printlen=ut.get_printlen(
                value=maxval))

        total_conflicts = stats['vars']['Tmfail']        + \
                          stats['vars']['repeatfail']      + \
                          stats['vars']['homodimerfail']   + \
                          stats['vars']['heterodimerfail'] + \
                          stats['vars']['exmotiffail']     + \
                          stats['vars']['edgefail']

        liner.send(
            ' Melt. Temp. Conflicts  : {:{},{}} Event(s) ({:6.2f} %)\n'.format(
                stats['vars']['Tmfail'],
                plen,
                sntn,
                ut.safediv(
                    A=stats['vars']['Tmfail'] * 100.,
                    B=total_conflicts)))
        liner.send(
            '      Repeat Conflicts  : {:{},{}} Event(s) ({:6.2f} %)\n'.format(
                stats['vars']['repeatfail'],
                plen,
                sntn,
                ut.safediv(
                    A=stats['vars']['repeatfail'] * 100.,
                    B=total_conflicts)))
        liner.send(
            '   Homodimer Conflicts  : {:{},{}} Event(s) ({:6.2f} %)\n'.format(
                stats['vars']['homodimerfail'],
                plen,
                sntn,
                ut.safediv(
                    A=stats['vars']['homodimerfail'] * 100.,
                    B=total_conflicts)))
        liner.send(
            ' Heterodimer Conflicts  : {:{},{}} Event(s) ({:6.2f} %)\n'.format(
                stats['vars']['heterodimerfail'],
                plen,
                sntn,
                ut.safediv(
                    A=stats['vars']['heterodimerfail'] * 100.,
                    B=total_conflicts)))
        liner.send(
            '     Exmotif Conflicts  : {:{},{}} Event(s) ({:6.2f} %)\n'.format(
                stats['vars']['exmotiffail'],
                plen,
                sntn,
                ut.safediv(
                    A=stats['vars']['exmotiffail'] * 100.,
                    B=total_conflicts)))
        liner.send(
            '        Edge Conflicts  : {:{},{}} Event(s) ({:6.2f} %)\n'.format(
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
                    '     - Exmotif {:>{}} Triggered {:{},{}} Event(s)\n'.format(
                        exmotif, qlen, count, vlen, sntn))

    # Show Time Elapsed
    liner.send(
        ' Time Elapsed: {:.2f} sec\n'.format(
            tt.time()-t0))

    # Unschedule outfile deletion
    if padstatus == 'Successful':
        ae.unregister(ofdeletion)

    # Close Liner
    liner.close()

    # Return Solution and Statistics
    return (outdf, stats)