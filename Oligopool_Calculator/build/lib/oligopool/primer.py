import time  as tt

import collections as cx
import atexit      as ae

import nrpcalc     as nr

import utils       as ut
import valparse    as vp
import coreprimer  as cp


def primer_engine(
    primerseq,
    primerspan,
    homology,
    primertype,
    fixedbaseindex,
    mintmelt,
    maxtmelt,
    maxreplen,
    oligorepeats,
    pairedprimer,
    pairedspan,
    pairedrepeats,
    exmotifs,
    exmotifindex,
    edgeeffectlength,
    prefixdict,
    suffixdict,
    background,
    stats,
    liner):
    '''
    Return a primer fulfilling all constraints.
    Internal use only.

    :: primerseq
       type - string
       desc - degenerate primer design sequence
              constraint
    :: primerspan
       type - integer / None
       desc - extraction span for primer sequence
              for dimer evaluation
    :: homology
       type - integer
       desc - maximum allowed internal repeat
              length for default traceback
    :: primertype
       type - integer
       desc - primer design type identifier
              0 = forward primer design
              1 = reverse primer design
    :: fixedbaseindex
       type - set
       desc - base indices lacking more than one
              nucleotide choices
    :: mintmelt
       type - float
       desc - primer melting temperature lower
              bound
    :: maxtmelt
       type - float
       desc - primer melting temperature upper
              bound
    :: maxreplen
       type - integer
       desc - maximum shared repeat length
    :: oligorepeats
       type - set
       desc - set storing oligopool repeats
    :: pairedprimer
       type - string / None
       desc - paired primer sequence
    :: pairedspan
       type - integer / None
       desc - extraction span for paired primer
              sequence for dimer evaluation
    :: pairedrepeats
       type - set / None
       desc - set storing paired primer repeats
    :: exmotifs
       type - cx.deque / None
       desc - deque of all excluded motifs
    :: exmotifindex
       type - set
       desc - ending index locations of all
              embedded excluded motifs
    :: edgeeffectlength
       type - integer
       desc - context length for edge effects
    :: prefixdict
       type - dict / None
       desc - dictionary of forbidden primer
              prefix sequences
    :: suffixdict
       type - dict / None
       desc - dictionary of forbidden primer
              suffix sequences
    :: background
       type - db.vectorDB / None
       desc - vectorDB instance storing
              background repeats
    :: stats
       type - dict
       desc - primer design stats storage
    :: liner
       type - coroutine
       desc - dynamic printing
    '''

    # Book-keeping
    t0 = tt.time() # Start Timer
    # Flexible Structure Constraint
    primerstruct = ''.join('x.'[idx in fixedbaseindex] \
        for idx in fixedbaseindex)

    # Extract Paired Primer Span
    if (not pairedspan   is None) and \
       (not pairedprimer is None):

        # Forward-Pair = Reverse Extraction
        if primertype == 0:
            pairedprimer = pairedprimer[:+pairedspan]
        # Reverse-Pair = Forward Extraction
        else:
            pairedprimer = pairedprimer[-pairedspan:]

    # Update Paired Primer Orientation
    # Since Forward Primer Design,
    # interpret Paired Primer as
    # Reverse Primer specified in
    # terms of Forward Strand
    if primertype == 0:
        if not pairedprimer is None:
            pairedprimer = ut.get_revcomp(
                seq=pairedprimer)
    # Correct Free Base Index
    else:
        fixedbaseindex = set(len(primerseq)-1-idx \
            for idx in fixedbaseindex)

    # Optimize exmotifs
    if not exmotifs is None:
        exmotifs = ut.get_grouped_sequences(
            sequences=exmotifs)

    # Define Objective Function
    objectivefunction = lambda primer: cp.primer_objectives(
        primer=primer,
        primerlen=len(primerseq),
        primerspan=primerspan,
        primertype=primertype,
        fixedbaseindex=fixedbaseindex,
        mintmelt=mintmelt,
        maxtmelt=maxtmelt,
        maxreplen=maxreplen,
        oligorepeats=oligorepeats,
        pairedprimer=pairedprimer,
        pairedrepeats=pairedrepeats,
        exmotifs=exmotifs,
        exmotifindex=exmotifindex,
        edgeeffectlength=edgeeffectlength,
        prefixforbidden=prefixdict,
        suffixforbidden=suffixdict,
        background=background,
        inittime=t0,
        stats=stats,
        liner=liner)

    # Define Maker Instance
    maker = nr.base.maker.NRPMaker(
        part_type='DNA',
        seed=None)

    # Design Primer via Maker
    primer = maker.nrp_maker(
        homology=min(len(primerseq), homology),
        seq_constr=primerseq,
        struct_constr=primerstruct,
        target_size=1,
        background=None,
        struct_type='both',
        synth_opt=False,
        local_model_fn=objectivefunction,
        jump_count=1000,
        fail_count=100000,
        output_file=None,
        verbose=False,
        abortion=True,
        allow_internal_repeat=False,
        check_constraints=False)

    # Check Status and Return Solution
    if len(primer) > 0:

        # Final Update
        cp.show_update(
            primer=primer[0],
            element='Primer',
            optstatus=2,
            optstate=0,
            inittime=t0,
            terminal=True,
            liner=liner)

        # Extract Primer
        primer  = primer[0]
        xprimer = primer[:]

        # Adjust Extraction
        if not primerspan is None:
            if primertype == 0:
                primer = primer[-primerspan:]
            else:
                primer = primer[:+primerspan]

        # We solved this!
        stats['status'] = True
        stats['basis']  = 'solved'

        # Correct Primer Orientation
        cprimer = ut.get_revcomp(
            seq=primer) if primertype == 1 else primer

        # Update Tm
        stats['vars']['primerTm'] = ut.get_tmelt(
            seq=cprimer)

        # Update GC Percentage
        stats['vars']['primerGC'] = (cprimer.count('G') + \
                                     cprimer.count('C')) \
                                        / (len(cprimer) * 0.01)

        # Update Hairpin Free Energy
        stats['vars']['hairpinMFE'] = cp.folder.evaluate_mfe(
            seq=cprimer,
            dg=True)[-1]

        # Update Heterodimer Free Energy
        stats['vars']['homodimerMFE'] = cp.folder.evaluate_mfe_dimer(
            seq1=cprimer,
            seq2=cprimer)[-1]

        # Update Homodimer Free Energy
        if not pairedprimer is None:
            stats['vars']['heterodimerMFE'] = cp.folder.evaluate_mfe_dimer(
                seq1=cprimer,
                seq2=pairedprimer)[-1]

        # Return Results
        return (xprimer, stats)

    # Design Unsuccessful
    else:

        # Final Update
        liner.send('\* Time Elapsed: {:.2f} sec\n'.format(
            tt.time() - t0))

        # Return Results
        return (None, stats)

def primer(
    indata,
    oligolimit,
    primerseq,
    primertype,
    mintmelt,
    maxtmelt,
    maxreplen,
    primercol,
    outfile=None,
    pairedcol=None,
    leftcontext=None,
    rightcontext=None,
    exmotifs=None,
    background=None,
    verbose=True):
    '''
    The primer function designs constrained primers, at desired
    melting temperature, with desired non-repetitiveness, that
    work optimally for all variants in the oligopool, without
    any excluded motifs inside or introducing new ones at the
    primer edges. Additional constraints are enforced to ensure
    compatibility with a paired primer, and minimal formation
    of dimers of duplexes. The generated DataFrame containing
    designed primer is returned, and optionally also written
    out to <outfile> (CSV) if specified.

    :: indata
       type - string / pd.DataFrame
       desc - path to CSV file or a pandas DataFrame storing
              annotated oligopool variants and their parts
    :: oligolimit
       type - integer
       desc - maximum oligo length allowed in the oligopool,
              must be 4 or greater
    :: primerseq
       type - integer
       desc - an IUPAC degenerate primer sequence constraint
    :: primertype
       type - integer
       desc - primer design type identifier
              0 = a forward primer is designed
              1 = a reverse primer is designed
    :: mintmelt
       type - float
       desc - minimum allowed primer melting temperature in
              degree celsius, must be 25 or greater
    :: maxtmelt
       type - float
       desc - maximum allowed primer melting temperature in
              degree celsius, must be 95 or lesser
    :: maxreplen
       type - integer
       desc - maximum shared repeat length between the
              primers and flanking regions, must be 6 or
              greater
    :: primercol
       type - string
       desc - name of the column to store designed primer
    :: outfile
       type - string
       desc - filename to save ouput DataFrame with primers
              (suffix='.oligopool.primer.csv')
              (default=None)
    :: pairedcol
       type - string / None
              name of the column containing the full primer
              sequence to be paired with the designed primer
              (default=None)
    :: leftcontext
       type - string / None
       desc - name of the column containing DNA sequences
              to the left of the primer sequence
              (default=None)
    :: rightcontext
       type - string / None
       desc - name of the column containing DNA sequences
              placed right of the primer sequence
              (default=None)
    :: exmotifs
       type - iterable / string / pd.DataFrame / None
       desc - iterable of DNA string motifs to be excluded
              within and at the edges of the primer when
              placed between context sequences; optionally,
              this can be a path to a CSV file containing
              uniquely identified excluded motifs or an
              equivalent pandas DataFrame
              (default=None)
    :: background
       type - string / None
       desc - path to directory storing the background
              sequence k-mers, against which the designed
              primer is to be non-repetitive
              (suffix='.oligopool.background')
              (default=None)
    :: verbose
       type - boolean
       desc - if True will log updates to stdout
              (default=True)

    Output: A file named <outfile> with '.oligoool.primer.csv'
            suffix if specified; otherwise a pandas DataFrame is
            returned, along with design or warning statistics.

    Note 1. Specified <indata> must contain a column named 'ID',
            that uniquely identifies variants in a pool. Values
            in <indata> except 'ID' must be DNA strings. All
            rows and columns in <indata> must be non-empty,
            i.e. none of the cells must be empty.

    Note 2. Column names in <indata> must be unique, without
            <primercol> as a pre-existing column name.

    Note 3. Either <leftcontext> or <rightcontext> or both must
            be specified. If both are specified then they must
            be adjacent to each other and in order. Designed
            primers would be inserted next to or between them.

    Note 4. The paired primer type is automatically inferred
            based on current primer type, i.e. if a forward
            primer is being designed, the paired primer is
            assumed to be a reverse primer, and optimization
            parameters are adjusted accordingly.

    Note 5. If a paired primer is specified, then the melting
            temperature of the designed primer is considered
            within 1 °C of the paired primer Tm. E.g. If input
            Tm range is 53 to 57 °C, and the paired primer has
            a Tm of 59 °C, the designed primer is optimized
            to have a Tm between 58 to 60 °C.

    Note 6. The <maxreplen> parameter here controls the level
            of non-repetitiveness in designed primers with
            respect to sequences in <indata> and as such is
            independent of the <maxreplen> used to specify
            the background.

    Note 7. If <exmotifs> points to a CSV file or DataFrame,
            it must contain both an 'ID' and an 'Exmotif'
            column, with 'Exmotif' containing all of the
            excluded motif sequences.

    Note 8. Constant motifs or bases in input primer sequence
            constraint may sometimes make it impossible to
            optimize for excluded motifs, edge-effects or prevent
            favorable thermodynamic properties. In such cases,
            a sub-optimal primer is designed and returned,
            along with any warnings.
    '''

    # Start Liner
    liner = ut.liner_engine(verbose)

    # Primer Verbage Print
    liner.send('\n[Oligopool Calculator: Design Mode - Primer]\n')

    # Required Argument Parsing
    liner.send('\n Required Arguments\n')

    # First Pass indata Parsing and Validation
    (indf,
    indata_valid) = vp.get_parsed_indata_info(
        indata=indata,
        indata_field='      Input Data       ',
        required_fields=('ID',),
        precheck=False,
        liner=liner)

    # Full oligolimit Validation
    oligolimit_valid = vp.get_numeric_validity(
        numeric=oligolimit,
        numeric_field='      Oligo Limit      ',
        numeric_pre_desc=' At most ',
        numeric_post_desc=' Base Pair(s)',
        minval=4,
        maxval=float('inf'),
        precheck=False,
        liner=liner)

    # First Pass primerseq Validation
    primerseq_valid = vp.get_seqconstr_validity(
        seqconstr=primerseq,
        seqconstr_field='     Primer Sequence   ',
        minlenval=15,
        element='PRIMER',
        liner=liner)

    # Full primertype Validation
    primertype_valid = vp.get_categorical_validity(
        category=primertype,
        category_field='     Primer Type       ',
        category_pre_desc=' ',
        category_post_desc=' Primer Design',
        category_dict={
            0: 'Forward',
            1: 'Reverse'},
        liner=liner)

    # Full mintmelt and maxtmelt Validation
    (mintmelt,
    maxtmelt,
    tmelt_valid) = vp.get_parsed_range_info(
        minval=mintmelt,
        maxval=maxtmelt,
        range_field='    Melting Temperature',
        range_unit='°C',
        range_min=25,
        range_max=95,
        liner=liner)

    # Full maxreplen Validation
    maxreplen_valid = vp.get_numeric_validity(
        numeric=maxreplen,
        numeric_field='     Repeat Length     ',
        numeric_pre_desc=' Up to ',
        numeric_post_desc=' Base Pair(s) Oligopool Repeats',
        minval=6,
        maxval=len(primerseq) if primerseq_valid else float('inf'),
        precheck=False,
        liner=liner)

    # Full primercol Validation
    primercol_valid = vp.get_parsed_column_info(
        col=primercol,
        df=indf,
        col_field='     Primer Column     ',
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
        outdf_suffix='.oligopool.primer.csv',
        outdf_field='     Output File       ',
        liner=liner)

    # Adjust outfile Suffix
    if not outfile is None:
        outfile = ut.get_adjusted_path(
            path=outfile,
            suffix='.oligopool.primer.csv')

    # Optional Argument Parsing
    liner.send('\n Optional Arguments\n')

    # Full pairedprimer Validation
    (pairedprimer,
    pairedprimer_valid) = vp.get_constantcol_validity(
        constantcol=pairedcol,
        constantcol_field='     Paired Primer     ',
        df=indf,
        liner=liner)

    # Store Context Names
    leftcontextname  = leftcontext
    rightcontextname = rightcontext

    # Full leftcontext Parsing and Validation
    (leftcontext,
    leftcontext_valid) = vp.get_parsed_column_info(
        col=leftcontext,
        df=indf,
        col_field='       Left Context    ',
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
        col_field='      Right Context    ',
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
        exseqs_field='   Excluded Motifs     ',
        exseqs_desc='Unique Motif(s)',
        df_field='Exmotif',
        required=False,
        liner=liner)

    # Full background Parsing and Validation
    (background,
    background_valid) = vp.get_parsed_background(
        background=background,
        background_field=' Background Database   ',
        liner=liner)

    # First Pass Validation
    if not all([
        indata_valid,
        oligolimit_valid,
        primerseq_valid,
        primertype_valid,
        tmelt_valid,
        maxreplen_valid,
        primercol_valid,
        outfile_valid,
        pairedprimer_valid,
        leftcontext_valid,
        rightcontext_valid,
        exmotifs_valid,
        background_valid]):
        liner.send('\n')
        raise RuntimeError(
            'Invalid Argument Input(s).')

    # Start Timer
    t0 = tt.time()

    # Adjust Numeric Paramters
    oligolimit = round(oligolimit)
    maxreplen  = round(maxreplen)

    # Define Edge Effect Length
    edgeeffectlength = None

    # Primer Design Book-keeping
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
        minelementlen=len(primerseq),
        maxelementlen=len(primerseq),
        element='Primer',
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
            element='Primer',
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
    liner.send('\n[Step 3: Parsing Primer Sequence]\n')

    # Update Step 3 Warning
    warns[3] = {
        'warncount': 0,
        'stepname' : 'parsing-primer-sequence',
        'vars': None}

    # Parse primerseq
    (parsestatus,
    homology,
    fixedbaseindex,
    exmotifindex,
    designspace,
    internalrepeats,
    pairedrepeats) = cp.get_parsed_sequence_constraint(
        primerseq=primerseq,
        primertype=primertype,
        exmotifs=exmotifs,
        pairedprimer=pairedprimer,
        warn=warns[3],
        liner=liner)

    # Remove Step 3 Warning
    if not warns[3]['warncount']:
        warns.pop(3)

    # primerseq infeasible
    if not parsestatus:

        # Prepare stats
        stats = {
            'status'  : False,
            'basis'   : 'infeasible',
            'step'    : 3,
            'stepname': 'parsing-primer-sequence',
            'vars'    : {
                    'designspace': designspace,
                'internalrepeats': internalrepeats,
                  'pairedrepeats': pairedrepeats},
            'warns'   : warns}

        # Return results
        liner.close()
        return (outdf, stats)

    # Define pairedrepeats
    if not pairedprimer is None:
        pairedrepeats = set(ut.stream_canon_spectrum(
            seq=pairedprimer,
            k=6))
    else:
        pairedrepeats = None

    # Parse Melting Temperature
    liner.send('\n[Step 4: Parsing Melting Temperature]\n')

    # Parse mintmelt and maxtmelt
    (parsestatus,
    estimatedminTm,
    estimatedmaxTm,
    higherminTm,
    lowermaxTm,
    mintmelt,
    maxtmelt) = cp.get_parsed_primer_tmelt_constraint(
        primerseq=primerseq,
        pairedprimer=pairedprimer,
        mintmelt=mintmelt,
        maxtmelt=maxtmelt,
        element='Primer',
        liner=liner)

    # mintmelt and maxtmelt infeasible
    if not parsestatus:

        # Prepare stats
        stats = {
            'status'  : False,
            'basis'   : 'infeasible',
            'step'    : 4,
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

    # Parse Edge Effects
    if not exmotifs is None and \
       ((not leftcontext  is None) or \
        (not rightcontext is None)):

        # Show update
        liner.send('\n[Step 5: Extracting Context Sequences]\n')

        # Extract Primer Contexts
        (leftcontext,
        rightcontext) = ut.get_extracted_context(
            leftcontext=leftcontext,
            rightcontext=rightcontext,
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
            primerseq=primerseq,
            leftcontext=leftcontext,
            rightcontext=rightcontext,
            leftpartition=leftpartition,
            rightpartition=rightpartition,
            exmotifs=exmotifs,
            element='Primer',
            warn=warns[6],
            liner=liner)

        # Remove Step 6 Warning
        if not warns[6]['warncount']:
            warns.pop(6)

    else:
        (prefixdict,
        suffixdict) = None, None

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
        element='Primer',
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

    # Launching Primer Design
    liner.send('\n[Step 8: Computing Primer]\n')

    # Define Primer Design Stats
    stats = {
        'status'  : False,
        'basis'   : 'unsolved',
        'step'    : 8,
        'stepname': 'computing-primer',
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

    # Schedule outfile deletion
    ofdeletion = ae.register(
        ut.remove_file,
        outfile)

    # Design Primer
    (primer,
    stats) = primer_engine(
        primerseq=primerseq,
        primerspan=None,
        homology=homology,
        primertype=primertype,
        fixedbaseindex=fixedbaseindex,
        mintmelt=mintmelt,
        maxtmelt=maxtmelt,
        maxreplen=maxreplen,
        oligorepeats=oligorepeats,
        pairedprimer=pairedprimer,
        pairedspan=None,
        pairedrepeats=pairedrepeats,
        exmotifs=exmotifs,
        exmotifindex=exmotifindex,
        edgeeffectlength=edgeeffectlength,
        prefixdict=prefixdict,
        suffixdict=suffixdict,
        background=background,
        stats=stats,
        liner=liner)

    # Primer Status
    if stats['status']:
        primerstatus = 'Successful'
    else:
        primerstatus = 'Failed'

    # Insert primer into indf
    if stats['status']:

        # Update indf
        ut.update_df(
            indf=indf,
            lcname=leftcontextname,
            rcname=rightcontextname,
            out=primer,
            outcol=primercol)

        # Prepare outdf
        outdf = indf

        # Write outdf to file
        if not outfile is None:
            outdf.to_csv(
                path_or_buf=outfile,
                sep=',')

    # Primer Design Statistics
    liner.send('\n[Primer Design Statistics]\n')

    liner.send(
        '      Design Status     : {}\n'.format(
            primerstatus))

    # Success Relevant Stats
    if stats['status']:

        liner.send(
            '     Melting Temperature: {:6.2f} °C\n'.format(
                stats['vars']['primerTm']))
        liner.send(
            '          GC Content    : {:6.2f} %\n'.format(
                stats['vars']['primerGC']))
        liner.send(
            '     Hairpin MFE        : {:6.2f} kcal/mol\n'.format(
                stats['vars']['hairpinMFE']))
        liner.send(
            '   Homodimer MFE        : {:6.2f} kcal/mol\n'.format(
                stats['vars']['homodimerMFE']))

        if not pairedprimer is None:
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

        total_conflicts = stats['vars']['Tmfail']          + \
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
    if primerstatus == 'Successful':
        ae.unregister(ofdeletion)

    # Close Liner
    liner.close()

    # Return Solution and Statistics
    return (outdf, stats)
