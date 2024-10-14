import time  as tt

import numpy    as np
import pandas   as pd
import atexit   as ae

import utils    as ut
import valparse as vp


def final(
    indata,
    outfile=None,
    verbose=True):
    '''
    TBW
    '''

    # Start Liner
    liner = ut.liner_engine(verbose)

    # Finalization Verbage Print
    liner.send('\n[Oligopool Calculator: Design Mode - Final]\n')

    # Required Argument Parsing
    liner.send('\n Required Arguments\n')

    # First Pass indata Parsing and Validation
    (indf,
    indata_valid) = vp.get_parsed_indata_info(
        indata=indata,
        indata_field='  Input Data',
        required_fields=('ID',),
        precheck=False,
        liner=liner)

    # Full outfile Validation
    outfile_valid = vp.get_outdf_validity(
        outdf=outfile,
        outdf_suffix='.oligopool.final.csv',
        outdf_field=' Output File',
        liner=liner)

    # Adjust outfile Suffix
    if not outfile is None:
        outfile = ut.get_adjusted_path(
            path=outfile,
            suffix='.oligopool.final.csv')

    # First Pass Validation
    if not all([
        indata_valid,
        outfile_valid]):
        liner.send('\n')
        raise RuntimeError(
            'Invalid Argument Input(s).')

    # Start Timer
    t0 = tt.time()

    # Schedule outfile deletion
    ofdeletion = ae.register(
        ut.remove_file,
        outfile)

    # Show Update
    liner.send('\n[Step 1: Finalizing Oligopool]\n')

    # Compute Final DataFrame
    outdf = pd.DataFrame(index=indf.index)
    outdf['CompleteOligo'] = ut.get_df_concat(df=indf)
    outdf['OligoLength']   = list(map(len, outdf['CompleteOligo']))

    # Show Update
    liner.send(' Finalization Completed\n')
    liner.send(
        ' Time Elapsed: {:.2f} sec\n'.format(
            tt.time()-t0))

    # Write outdf to file
    if not outfile is None:
        outdf.to_csv(
            path_or_buf=outfile,
            sep=',')

    # Compute Stats
    minoligolen = np.min(outdf.OligoLength)
    maxoligolen = np.max(outdf.OligoLength)

    # Build Stats Dictionary
    stats = {
        'status'  : True,
        'basis'   : 'solved',
        'step'    : 1,
        'stepname': 'finalizing-oligopool',
        'vars'    : {
            'minoligolen': minoligolen,
            'maxoligolen': maxoligolen},
        'warns'   : {}}

    # Finalization Statistics
    liner.send('\n[Finalization Statistics]\n')

    plen = ut.get_printlen(
        value=max(stats['vars'][field] for field in (
            'minoligolen',
            'maxoligolen')))

    liner.send(
        ' Final Status: Successful\n')

    if minoligolen == maxoligolen:
        liner.send(
            ' Oligo Length: {:{},d} Base Pair(s)\n'.format(
                stats['vars']['minoligolen'],
                plen))
    else:
        liner.send(
            ' Oligo Length: {:{},d} to {:{},d} Base Pair(s)\n'.format(
                stats['vars']['minoligolen'],
                plen,
                stats['vars']['maxoligolen'],
                plen))

    # Show Time Elapsed
    liner.send(
        ' Time Elapsed: {:.2f} sec\n'.format(
            tt.time()-t0))

    # Unschedule outfile deletion
    ae.unregister(ofdeletion)

    # Close Liner
    liner.close()

    # Return Solution and Statistics
    return (outdf, stats)