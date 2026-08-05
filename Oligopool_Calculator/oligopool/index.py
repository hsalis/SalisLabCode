import time as tt

import atexit    as ae

import utils     as ut
import valparse  as vp
import scry      as sy
import coreindex as ci


def index_engine(
    IDdict,
    barcodedict,
    barcodename,
    barcodecount,
    barcodelen,
    barcodeprefix,
    barcodesuffix,
    barcodepregap,
    barcodepostgap,
    associatedict,
    associateprefix,
    associatesuffix,
    indexdir,
    indexqueue,
    liner):
    '''
    Prepare final indexes and models for
    counting, based on parsed objects.
    Internal use only.

    :: IDdict
       type - dict
       desc - barcode/associate ID dictionary
    :: barcodedict
       type - dict
       desc - complete barcode sequence
              dictionary
    :: barcodename
       type - string
       desc - barcode column name
    :: barcodecount
       type - integer
       desc - total number of barcodes
    :: barcodelen
       type - integer
       desc - length of barcodes
    :: barcodeprefix
       type - string / None
       desc - constant barcode prefix
    :: barcodesuffix
       type - string / None
       desc - constant barcode suffix
    :: barcodepregap
       type - integer
       desc - gap between prefix and barcode
    :: barcodepostgap
       type - integer
       desc - gap between barcode and suffix
    :: associatedict
       type - dict
       desc - complete associate sequence
              dictionary
    :: associateprefix
       type - string / None
       desc - constant associate prefix
    :: associatesuffix
       type - string / None
       desc - constant associate suffix
    :: indexdir
       type - string
       desc - path to directory temporarily
              storing indexed structures
              and data models
    :: indexqueue
       type - mp.SimpleQueue
       desc - queue storing indexed object
              paths once saved into indexdir
    :: liner
       type - coroutine
       desc - dynamic printing
    '''

    # Start Timer
    t0 = tt.time()

    # Base Index Directory Path
    indexdir += '/'

    # Counting Meta Data Initialization
    maxtval = 2
    metamap = {
        'variantcount'   : barcodecount,
        'barcodename'    : barcodename,
        'barcodelen'     : barcodelen,
        'barcodeprefix'  : barcodeprefix,
        'barcodesuffix'  : barcodesuffix,
        'barcodepregap'  : barcodepregap,
        'barcodepostgap' : barcodepostgap,
        'barcodegapped'  : any((barcodepregap, barcodepostgap)),
        'association'    : len(associatedict) > 0,
        'associateprefix': associateprefix,
        'associatesuffix': associatesuffix
    }

    # Compute Barcode Set t-value
    liner.send(' Inferring Barcode t-value ...')
    barcodetval = ci.infer_tvalue(
        elementlen=barcodelen,
        maximum=2)
    liner.send(
        '   Barcode t-value: {:,} Mismatch(es)\n'.format(
            barcodetval))
    metamap['barcodetval'] = barcodetval

    # Compute Barcode Set k-value
    liner.send(' Inferring Barcode k-value ...')
    barcodekval = ci.infer_kvalue(
        elementlen=barcodelen,
        tvalue=barcodetval,
        minimum=3)
    liner.send(
        '   Barcode k-value: {:,} Base Pair(s)\n'.format(
            barcodekval))
    metamap['barcodekval'] = barcodekval

    # Compute Associate Set t-value
    if associatedict:
        liner.send(' Inferring Associate t-value ...')
        (associatetvaldict,
        mintval,
        maxtval) = ci.get_associate_tvalues(
            associatedict=associatedict,
            maximum=4,
            liner=liner)
        if mintval == maxtval:
            liner.send(
                ' Associate t-value: {:,} Mismatch(es)\n'.format(
                    mintval))
        else:
            liner.send(
                ' Associate t-value: {:,} to {:,} Mismatch(es)\n'.format(
                    mintval,
                    maxtval))
        metamap['associatetvalmin'] = mintval
        metamap['associatetvalmax'] = maxtval
    else:
        metamap['associatetvalmin'] = None
        metamap['associatetvalmax'] = None

    # Compute Barcode Prefix t-value
    if not barcodeprefix is None:
        bpxtval = ci.infer_tvalue(
            elementlen=len(barcodeprefix),
            maximum=2)
        metamap['bpxtval'] = bpxtval
    else:
        metamap['bpxtval'] = None

    # Compute Barcode Suffix t-value
    if not barcodesuffix is None:
        bsxtval = ci.infer_tvalue(
            elementlen=len(barcodesuffix),
            maximum=2)
        metamap['bsxtval'] = bsxtval
    else:
        metamap['bsxtval'] = None

    # Compute Trimming Policy
    metamap['trimpolicy'] = None
    if not barcodeprefix is None and \
       not barcodesuffix is None:
        metamap['trimpolicy'] = 1.5
    elif not barcodeprefix is None:
        metamap['trimpolicy'] = 1
    elif not barcodesuffix is None:
        metamap['trimpolicy'] = 2

    # Compute Read Anchor Motif
    if not barcodeprefix is None:
        metamap['anchor']     = barcodeprefix
        metamap['revanchor']  = ut.get_revcomp(
            seq=barcodeprefix)
        metamap['anchortval'] = bpxtval
    else:
        metamap['anchor']     = barcodesuffix
        metamap['revanchor']  = ut.get_revcomp(
            seq=barcodesuffix)
        metamap['anchortval'] = bsxtval

    # Compute Associate Prefix t-value
    if not associateprefix is None:
        apxtval = ci.infer_tvalue(
            elementlen=len(associateprefix),
            maximum=2)
        metamap['apxtval'] = apxtval
    else:
        metamap['apxtval'] = None

    # Compute Associate Suffix t-value
    if not associatesuffix is None:
        asxtval = ci.infer_tvalue(
            elementlen=len(associatesuffix),
            maximum=2)
        metamap['asxtval'] = asxtval
    else:
        metamap['asxtval'] = None

    # Save and Delete ID Dict
    liner.send(' Writing ID Map ...')
    opath = indexdir+'ID.map'
    ut.savedict(
        dobj=IDdict,
        filepath=opath)
    indexqueue.put(opath)
    del IDdict
    liner.send(' Writing        ID Map  : Done\n')

    # Save and Delete metamap
    liner.send(' Writing Meta Map ...')
    opath = indexdir+'meta.map'
    ut.savedict(
        dobj=metamap,
        filepath=opath)
    indexqueue.put(opath)
    del metamap
    liner.send(' Writing      Meta Map  : Done\n')

    # Train Barcode Model
    liner.send(' Training Barcode Model ...')
    X,Y = ci.split_map(
        maptosplit=barcodedict)
    t0  = tt.time()
    barcodemodel = sy.Scry().train(
        X=X,
        Y=Y,
        n=barcodelen,
        k=barcodekval,
        t=barcodetval,
        liner=liner)

    # Save Barcode Model
    opath = indexdir+'barcode.model'
    ut.savemodel(
        mobj=barcodemodel,
        filepath=opath)
    indexqueue.put(opath)
    del barcodedict
    del barcodemodel
    liner.send(' Writing   Barcode Model: Done\n')

    # Update Associate Dictionary
    if associatedict:
        liner.send(' Updating Associate Map ...')
        for k in associatedict:
            associatedict[k] = (
                associatedict[k],
                associatetvaldict.pop(k))

    # Save and Delete associatedict
    if associatedict:
        liner.send(' Writing Associate Map ...')
        opath = indexdir+'associate.map'
        ut.savedict(
            dobj=associatedict,
            filepath=opath)
        indexqueue.put(opath)
        del associatedict
        liner.send(' Writing Associate Map  : Done\n')

    # Show Time Elapsed
    liner.send(
        ' Time Elapsed: {:.2f} sec\n'.format(
            tt.time()-t0))

    # Indexing Completed
    indexqueue.put(None)

def index(
    barcodedata,
    barcodecol,
    indexfile,
    barcodeprefix=None,
    barcodesuffix=None,
    barcodepregap=0,
    barcodepostgap=0,
    associatedata=None,
    associatecol=None,
    associateprefix=None,
    associatesuffix=None,
    verbose=True):
    '''
    TBW
    '''

    # Start Liner
    liner = ut.liner_engine(verbose)

    # Indexing Verbage Print
    liner.send('\n[Oligopool Calculator: Analysis Mode - Index]\n')

    # Required Argument Parsing
    liner.send('\n Required Arguments\n')

    # First Pass barcodedata Parsing and Validation
    (bardf,
    barcodedata_valid) = vp.get_parsed_indata_info(
        indata=barcodedata,
        indata_field='   Barcode Data   ',
        required_fields=('ID',),
        precheck=False,
        liner=liner)

    # Sort bardf
    if barcodedata_valid:
        bardf.sort_index(
            inplace=True)

    # Full barcodecol Validation
    (barcodes,
    barcodecol_valid) = vp.get_parsed_column_info(
        col=barcodecol,
        df=bardf,
        col_field='   Barcode Column ',
        col_desc='Input from Column',
        col_type=0,
        adjcol=None,
        adjval=None,
        iscontext=False,
        typecontext=None,
        liner=liner)

    # Full indexfile Validation
    indexfile_valid = vp.get_outfile_validity(
        outfile=indexfile,
        outfile_suffix='.oligopool.index',
        outfile_field='     Index File   ',
        liner=liner)

    # Adjust indexfile Suffix
    indexfile = ut.get_adjusted_path(
        path=indexfile,
        suffix='.oligopool.index')

    # Optional Argument Parsing
    liner.send('\n Optional Arguments\n')

    # Full barcodeprefix Validation
    (barcodeprefix,
    barcodeprefix_valid) = vp.get_parsed_column_info(
        col=barcodeprefix,
        df=bardf,
        col_field='   Barcode Prefix ',
        col_desc='Input from Column',
        col_type=0,
        adjcol=barcodecol,
        adjval=+1,
        iscontext=False,
        typecontext=None,
        liner=liner)

    # Full barcodesuffix Validation
    (barcodesuffix,
    barcodesuffix_valid) = vp.get_parsed_column_info(
        col=barcodesuffix,
        df=bardf,
        col_field='   Barcode Suffix ',
        col_desc='Input from Column',
        col_type=0,
        adjcol=barcodecol,
        adjval=-1,
        iscontext=False,
        typecontext=None,
        liner=liner)

    # Full barcodepregap Validation
    barcodepregap_valid = vp.get_numeric_validity(
        numeric=barcodepregap,
        numeric_field='   Barcode Pregap ',
        numeric_pre_desc=' Allow Exactly ',
        numeric_post_desc=' bp Gap b/w Prefix and Barcode',
        minval=0,
        maxval=float('inf'),
        precheck=False,
        liner=liner)

    # Full barcodepostgap Validation
    barcodepostgap_valid = vp.get_numeric_validity(
        numeric=barcodepostgap,
        numeric_field='   Barcode Postgap',
        numeric_pre_desc=' Allow Exactly ',
        numeric_post_desc=' bp Gap b/w Barcode and Suffix',
        minval=0,
        maxval=float('inf'),
        precheck=False,
        liner=liner)

    # First Pass associatedata Parsing and Validation
    (assdf,
    associatedata_valid) = vp.get_parsed_associatedata_info(
        associatedata=associatedata,
        associatedata_field=' Associate Data   ',
        required_fields=('ID',),
        bardf=bardf,
        barcodedata_valid=barcodedata_valid,
        liner=liner)

    # Full associatecol Validation
    (associates,
    associatecol_valid) = vp.get_parsed_column_info(
        col=associatecol,
        df=assdf,
        col_field=' Associate Column ',
        col_desc='Input from Column',
        col_type=0,
        adjcol=None,
        adjval=None,
        iscontext=False,
        typecontext=None,
        liner=liner)

    # Full associateprefix Validation
    (associateprefix,
    associateprefix_valid) = vp.get_parsed_column_info(
        col=associateprefix,
        df=assdf,
        col_field=' Associate Prefix ',
        col_desc='Input from Column',
        col_type=0,
        adjcol=associatecol,
        adjval=+1,
        iscontext=False,
        typecontext=None,
        liner=liner)

    # Full associatesuffix Validation
    (associatesuffix,
    associatesuffix_valid) = vp.get_parsed_column_info(
        col=associatesuffix,
        df=assdf,
        col_field=' Associate Suffix ',
        col_desc='Input from Column',
        col_type=0,
        adjcol=associatecol,
        adjval=-1,
        iscontext=False,
        typecontext=None,
        liner=liner)

    # First Pass Validation
    if not all([
        barcodedata_valid,
        barcodecol_valid,
        indexfile_valid,
        barcodeprefix_valid,
        barcodesuffix_valid,
        barcodepregap_valid,
        barcodepostgap_valid,
        associatedata_valid,
        associatecol_valid,
        associateprefix_valid,
        associatesuffix_valid]):
        liner.send('\n')
        raise RuntimeError(
            'Invalid Argument Input(s).')

    # Start Timer
    t0 = tt.time()

    # Setup Warning Dictionary
    warns = {}

    # Parse Barcode Data Feasibility
    liner.send('\n[Step 1: Parsing Barcode Data]\n')

    # Parse bardf
    (parsestatus,
    IDdict,
    barcodedict,
    duplicates,
    barcodecount,
    barcodesuniq,
    barcodelen,
    barcodelenuniq,
    minbarcodelen,
    maxbarcodelen) = ci.get_parsed_barcodes(
        bardf=bardf,
        barcodes=barcodes,
        liner=liner)

    # Barcode Data infeasible
    if not parsestatus:

        # Prepare stats
        stats = {
            'status'  : False,
            'basis'   : 'infeasible',
            'step'    : 1,
            'stepname': 'parsing-barcode-data',
            'vars'    : {
                      'duplicates': duplicates,
                    'barcodecount': barcodecount,
                  'barcodesunique': barcodesuniq,
                   'minbarcodelen': minbarcodelen,
                   'maxbarcodelen': maxbarcodelen,
                'barcodelenunique': barcodelenuniq},
            'warns'   : warns}

        # Return results
        liner.close()
        return stats

    # Parse Barcode Constant Feasibility
    liner.send('\n[Step 2: Parsing Barcode Constants]\n')

    # Parse barcodeprefix and barcodesuffix
    (parsestatus,
    constantsextracted,
    constantsuniq,
    longconstants,
    barcodeprefix,
    barcodesuffix,
    prefixuniq,
    suffixuniq,
    prefixlen,
    suffixlen) = ci.get_parsed_constants(
        prefixconstant=barcodeprefix,
        suffixconstant=barcodesuffix,
        attachetype=0,
        liner=liner)

    # Barcode Constants infeasible
    if not parsestatus:

        # Prepare stats
        stats = {
            'status'  : False,
            'basis'   : 'infeasible',
            'step'    : 2,
            'stepname': 'parsing-barcode-constants',
            'vars'    : {
                'constantsextracted': constantsextracted,
                   'constantsunique': constantsuniq,
                     'longconstants': longconstants,
                      'prefixunique': prefixuniq,
                      'suffixunique': suffixuniq,
                         'prefixlen': prefixlen,
                         'suffixlen': suffixlen},
            'warns'   : warns}

        # Return results
        liner.close()
        return stats

    # Do we have associated variants?
    if not associatedata is None:

        # Extract Associate Data
        liner.send('\n[Step 3: Extracting Associate Data]\n')

        # Update Step 3 Warning
        warns[3] = {
            'warncount': 0,
            'stepname' : 'parsing-associate-data',
            'vars': None}

        # Extract associatedata
        associatedict = ci.get_extracted_associates(
            associates=associates,
            expcount=barcodecount,
            IDdict=IDdict,
            warn=warns[3],
            liner=liner)

        # Parse Associate Constant Feasibility
        liner.send('\n[Step 4: Parsing Associate Constants]\n')

        # Parse associateprefix and associatesuffix
        (parsestatus,
        constantsextracted,
        constantsuniq,
        longconstants,
        associateprefix,
        associatesuffix,
        prefixuniq,
        suffixuniq,
        prefixlen,
        suffixlen) = ci.get_parsed_constants(
            prefixconstant=associateprefix,
            suffixconstant=associatesuffix,
            attachetype=1,
            liner=liner)

        # Associate Constants infeasible
        if not parsestatus:

            # Prepare stats
            stats = {
                'status'  : False,
                'basis'   : 'infeasible',
                'step'    : 4,
                'stepname': 'parsing-associate-constants',
                'vars'    : {
                    'constantsextracted': constantsextracted,
                       'constantsunique': constantsuniq,
                         'longconstants': longconstants,
                          'prefixunique': prefixuniq,
                          'suffixunique': suffixuniq,
                             'prefixlen': prefixlen,
                             'suffixlen': suffixlen},
                'warns'   : warns}

            # Return results
            liner.close()
            return stats

    else:
        associatedict   = {}
        associateprefix = None
        associatesuffix = None

    # Setup Workspace
    (indexfile,
    indexdir) = ut.setup_workspace(
        outfile=indexfile,
        outfile_suffix='.oligopool.index')

    # Schedule indexfile deletion
    ifdeletion = ae.register(
        ut.remove_file,
        indexfile)

    # Prepared Objects Queue
    indexqueue = ut.SafeQueue()

    # Launching Indexing
    liner.send('\n[Step 5: Computing Index]\n')

    # Define Index Stats
    stats = {
        'status'  : True,
        'basis'   : 'solved',
        'step'    : 5,
        'stepname': 'computing-index',
        'vars'    : None,
        'warns'   : warns}

    # Compute Index Objects
    index_engine(
        IDdict=IDdict,
        barcodedict=barcodedict,
        barcodename=barcodecol,
        barcodecount=barcodecount,
        barcodelen=barcodelen,
        barcodeprefix=barcodeprefix,
        barcodesuffix=barcodesuffix,
        barcodepregap=barcodepregap,
        barcodepostgap=barcodepostgap,
        associatedict=associatedict,
        associateprefix=associateprefix,
        associatesuffix=associatesuffix,
        indexdir=indexdir,
        indexqueue=indexqueue,
        liner=liner)

    # Archive All Indexed Objects
    ut.archive(
        objqueue=indexqueue,
        arcfile=indexfile,
        mode='x',
        prodcount=1,
        prodactive=0,
        liner=liner)

    # Indexing Status?
    indexstatus = 'Successful'

    # Indexing Stats
    liner.send('\n[Indexing Statistics]\n')

    liner.send(
        ' Indexing Status : {}\n'.format(
            indexstatus))
    liner.send(
        '     Time Elapsed: {:.2f} sec\n'.format(
            tt.time() - t0))

    # Remove Workspace
    ut.remove_directory(
        dirpath=indexdir)

    # Unschedule packfile deletion
    if indexstatus == 'Successful':
        ae.unregister(ifdeletion)

    # Close Liner
    liner.close()

    return stats