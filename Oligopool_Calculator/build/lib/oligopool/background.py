import sys
import time     as tt
import shutil   as sh

import leveldb  as lv
import nrpcalc  as nr
import atexit   as ae

import utils    as ut
import vectordb as db
import valparse as vp


def background_engine(
    background,
    maxreplen,
    outdir,
    stats,
    liner):
    '''
    Extract and populate background.
    Internal use only.

    :: background
       type - list
       desc - list of background sequences
    :: maxreplen
       type - integer
       desc - maximum shared repeat length
    :: outdir
       type - string
       desc - output directory
    :: stats
       type - dict
       desc - background stats storage
    :: liner
       type - coroutine
       desc - dynamic printing
    '''

    # Open vectorDB instance
    vDB = db.vectorDB(
        path=outdir,
        maxreplen=maxreplen,
        mode=0)

    # Book-keeping
    t0   = tt.time()
    plen = ut.get_printlen(
        value=len(background))

    # Loop and insert background
    for idx,seq in enumerate(background):

        # Format sequence
        if len(seq) > 20:
            pseq = seq[:20]
            pbuf = '...'
        else:
            pseq = seq
            pbuf = ''

        # Insert sequence
        vDB.add(
            seq=seq,
            rna=False)

        # Show updates
        liner.send(
            ' Sequence {:{},d}: {}{} Inserted'.format(
                idx+1, plen, pseq, pbuf))

    # Final Update
    liner.send(
        ' Sequence {:{},d}: {}{} Inserted\n'.format(
            idx+1, plen, pseq, pbuf))
    liner.send(' Time Elapsed: {:.2f} sec\n'.format(
        tt.time()-t0))

    # Populate Stats
    kmerspace = ((4**(maxreplen+1)) // 2)
    fillcount = min(kmerspace, len(vDB))
    leftcount = kmerspace - fillcount

    stats['status'] = (leftcount * 1.) / kmerspace > .01
    stats['basis']  = 'solved' if stats['status'] else 'infeasible'
    stats['vars']['kmerspace'] = kmerspace
    stats['vars']['fillcount'] = fillcount
    stats['vars']['leftcount'] = leftcount

    # If Successful Update and Close DB
    if stats['status']:
        vDB.DB.Put(
            b'LEN',
            str(fillcount).encode())
        vDB.close()

    # Otherwise Drop DB
    else:
        vDB.drop()

    # Return Stats
    return stats

def background(
    indata,
    maxreplen,
    outdir,
    verbose=True):
    '''
    The background function builds a database of k-mers from
    a list or a CSV file of background sequences. This database
    is then specified during primer design, to generate primers
    which are non-repetitive to background thus minimizing any
    off-target amplification. Non-repetitiveness is controlled
    via the maximum shared repeat length (replength) parameter.
    Generated database is stored in <outdir>.

    :: indata
       type - iterable / string / pd.DataFrame
       desc - iterable of DNA strings against which designed
              primers are ensured to be non-repetitive;
              optionally, this can be a path to a CSV file
              containing uniquely identified background
              sequences or an equivalent pandas DataFrame
    :: maxreplen
       type - integer
       desc - maximum shared repeat length between the primers
              and the background sequences, must be between
              6 to 20
    :: outdir
       type - string
       desc - path to store the generated background k-mer
              database
    :: verbose
       type - boolean
       desc - if True will log updates to stdout
              (default=True)

    Output: A directory <outdir> with '.oligoool.background'
            suffix.

    Note 1. If <indata> points to a CSV file or a DataFrame,it
            must contain a column named 'ID', that uniquely
            identifies each background sequence, listed in a
            'Sequence' column. Values in <indata> except 'ID'
            must be DNA strings. All rows and columns in the
            <indata> must be non-empty, i.e. none of the cells
            must be empty.

    Note 2. The <maxreplen> parameter here controls the level
            of non-repetitiveness in designed primers with
            respect to a background sequences such as a genome
            or a plasmid, and as such is independent of the
            <maxreplen> used in primer design which controls
            the non-repetitiveness of the primers against the
            core oligopool variants.
    '''

    # Start Liner
    liner = ut.liner_engine(verbose)

    # Barcoding Verbage Print
    liner.send('\n[Oligopool Calculator: Design Mode - Background]\n')

    # Required Argument Parsing
    liner.send('\n Required Arguments\n')

    # First Pass indata Parsing and Validation
    (background,
    background_valid) = vp.get_parsed_exseqs_info(
        exseqs=indata,
        exseqs_field=' Background      Data',
        exseqs_desc='Unique Sequence(s)',
        df_field='Sequence',
        required=True,
        liner=liner)

    # Full maxreplen Validation
    maxreplen_valid = vp.get_numeric_validity(
        numeric=maxreplen,
        numeric_field='    Maximum    Repeat',
        numeric_pre_desc=' Up to ',
        numeric_post_desc=' Base Pair(s) Background Repeats',
        minval=6,
        maxval=20,
        precheck=False,
        liner=liner)

    # Full outdir Validation
    outdir_valid = vp.get_outdir_validity(
        outdir=outdir,
        outdir_suffix='.oligopool.background',
        outdir_field='     Output Directory',
        liner=liner)

    # First Pass Validation
    if not all([
        background_valid,
        maxreplen_valid,
        outdir_valid]):
        liner.send('\n')
        raise RuntimeError(
            'Invalid Argument Input(s).')

    # Adjust Numeric Paramters
    maxreplen = round(maxreplen)

    # Adjust outdir Suffix
    outdir = ut.get_adjusted_path(
        path=outdir,
        suffix='.oligopool.background')

    # Schedule outdir deletion
    oddeletion = ae.register(
        ut.remove_directory,
        outdir)

    # Launching Background Extraction
    liner.send('\n[Computing Background]\n')

    # Define Background Stats
    stats = {
        'status'  : False,
        'basis'   : 'infeasible',
        'step'    : 1,
        'stepname': 'computing-background',
        'vars'    : {
            'kmerspace': 0,  # kmer Space
            'fillcount': 0,  # kmer Fill Count
            'leftcount': 0}, # kmer Left Count
        'warns'   : {}}

    # Start Timer
    t0 = tt.time()

    # Extract Background
    stats = background_engine(
        background=background,
        maxreplen=maxreplen,
        outdir=outdir,
        stats=stats,
        liner=liner)

    # Counting Status
    if stats['status']:
        backgroundstatus = 'Successful'
    else:
        backgroundstatus = 'Failed'

    # Background Statistics
    liner.send('\n[Background Statistics]\n')

    plen = ut.get_printlen(
        value=max(stats['vars'][field] for field in (
            'kmerspace',
            'fillcount',
            'leftcount')))

    sntn = 'e' if plen > 15 else 'd'

    liner.send(
        ' Background Status: {}\n'.format(
            backgroundstatus))
    liner.send(
        '      k-mer  Space: {:{},{}} Unique {:,}-mers\n'.format(
            stats['vars']['kmerspace'],
            plen,
            sntn,
            maxreplen+1))
    liner.send(
        '       Fill  Count: {:{},{}} Unique {:,}-mers ({:6.2f} %)\n'.format(
            stats['vars']['fillcount'],
            plen,
            sntn,
            maxreplen+1,
            ut.safediv(
                A=stats['vars']['fillcount'] * 100.,
                B=stats['vars']['kmerspace'])))
    liner.send(
        '       Left  Count: {:{},{}} Unique {:,}-mers ({:6.2f} %)\n'.format(
            stats['vars']['leftcount'],
            plen,
            sntn,
            maxreplen+1,
            ut.safediv(
                A=stats['vars']['leftcount'] * 100.,
                B=stats['vars']['kmerspace'])))

    liner.send(' Time Elapsed: {:.2f} sec\n'.format(tt.time()-t0))

    # Unschedule outfile deletion
    if backgroundstatus == 'Successful':
        ae.unregister(oddeletion)

    # Close Liner
    liner.close()

    # Return Statistics
    return stats
