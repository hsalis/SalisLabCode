'''
Automated Design Mode Parser.

Look at examples/oligopool_design_example.py to show how to use this parser.
'''

import sys, os, collections, uuid, atexit, pandas

# Oligopool Calculator modules

import background
import motif
import primer
import barcode
import spacer
import split
import pad
import final
import utils

def setup_workspace(workspacename):

    # Schedule outdir deletion
    atexit.register(
        utils.remove_directory,
        workspacename)

    # Setup outdir
    utils.setup_directory(
        dirpath=workspacename)

def setup_dataframe(
    pool_size,
    element_names):

    # Blank DataFrame
    dataframe = pandas.DataFrame()

    # Build IDs
    dataframe.insert(
        0, 'ID',  ['V-{}'.format(idx+1) for idx in range(pool_size)])

    # Build Element Columns
    for idx,colname in enumerate(element_names):
        dataframe.insert(idx, colname, '-')

    # Re-Index
    dataframe.set_index('ID', inplace=True)

    # Return Results
    return dataframe

def setup_background(
    workspacename,
    background_spec,
    stats_dict):

    # Build Background
    status = False
    if (not (background_spec is None)) and \
       (len(background_spec) > 0):

        # Define Background Name
        backgroundname = workspacename + '/oligopool_background'

        # Build Background
        try:
            stats = background.background(
                indata=background_spec['indata'],
                maxreplen=background_spec['maxreplen'],
                outdir=backgroundname)
            stats_dict['background'] = stats
            status = True
        except:
            stats = {
                'status'  : False,
                'basis'   : 'infeasible',
                'step'    : 0,
                'stepname': 'parsing-input',
                'vars'    : {},
                'warns'   : {}}
            status = False
    else:
        print()
        backgroundname = None

    # Return Results
    return status, backgroundname, stats_dict

def insert_variants(
    dataframe,
    elements_spec,
    stats_dict):

    # Process All Variants
    status = False
    for element_name in elements_spec:
        if elements_spec[element_name]['type'] == 'variant':
            try:
                dataframe[element_name] = elements_spec[element_name]['sequences']
                stats_dict['oligopool'][element_name] = {
                    'status'  : True,
                    'basis'   : 'solved',
                    'step'    : 1,
                    'stepname': 'inserting-variants',
                    'vars'    : {},
                    'warns'   : {}
                }
                status = True
            except:
                dataframe = None
                stats_dict['oligopool'][element_name] = {
                    'status'  : False,
                    'basis'   : 'infeasible',
                    'step'    : 0,
                    'stepname': 'parsing-input',
                    'vars'    : {},
                    'warns'   : {}
                }
                status = False
                break

    # Show Update
    print('[Oligopool Calculator: Design Mode - Variants]\n')
    print(dataframe)

    # Return Results
    return status, dataframe, stats_dict

def insert_motifs(
    dataframe,
    elements_spec,
    stats_dict):

    # Process All Variants
    status = False
    incomplete = False
    for element_name in elements_spec:
        if elements_spec[element_name]['type'] == 'motif':
            try:
                if element_name in dataframe.columns:
                    dataframe.drop(element_name, inplace=True, axis=1)
                dataframe, stats = motif.motif(
                    indata=dataframe,
                    oligolimit=elements_spec[element_name]['oligolimit'],
                    motifseq=elements_spec[element_name]['motifseq'],
                    motifcol=element_name,
                    outfile=None,
                    leftcontext=elements_spec[element_name]['leftcontext'],
                    rightcontext=elements_spec[element_name]['rightcontext'],
                    exmotifs=elements_spec[element_name]['exmotifs'])
                stats_dict['oligopool'][element_name] = stats
                status = True
                if stats['status'] is False:
                    incomplete = True
                    break
            except:
                dataframe = None
                stats_dict['oligopool'][element_name] = {
                    'status'  : False,
                    'basis'   : 'infeasible',
                    'step'    : 0,
                    'stepname': 'parsing-input',
                    'vars'    : {},
                    'warns'   : {}}
                status = False
                break

    # Return Results
    return status, incomplete, dataframe, stats_dict

def get_primer_order(
    elements_spec):

    order = collections.deque()
    for element_name in elements_spec:
        if elements_spec[element_name]['type'] == 'primer':
            if not element_name in order:
                order.append(element_name)
            pairedprimer = elements_spec[element_name]['pairedprimer']
            if not pairedprimer is None and \
               not pairedprimer in order:
                index = order.index(element_name)
                order.insert(index, pairedprimer)
    return order

def insert_primers(
    element_names,
    background,
    dataframe,
    elements_spec,
    stats_dict):

    # Process All Variants
    status = False
    incomplete = False
    order = get_primer_order(elements_spec=elements_spec)
    # print(order)
    # print(dataframe.columns)
    # print(element_names)
    while order:
        element_name = order.popleft()
        try:
            if element_name in dataframe.columns:
                dataframe.drop(element_name, inplace=True, axis=1)
            pairedcol = elements_spec[element_name]['pairedprimer']
            if pairedcol in order:
                pairedcol = None # element_name is a Predecessor
            # Left Context is Sensitive
            if elements_spec[element_name]['leftcontext'] in stats_dict:
                leftcontext = elements_spec[element_name]['leftcontext']
            else:
                leftcontext = None
            # Right Context is Sensitive
            if elements_spec[element_name]['rightcontext'] in stats_dict:
                rightcontext = elements_spec[element_name]['rightcontext']
            else:
                rightcontext = None
            dataframe, stats = primer.primer(
                indata=dataframe,
                oligolimit=elements_spec[element_name]['oligolimit'],
                primerseq=elements_spec[element_name]['primerseq'],
                primertype=elements_spec[element_name]['primertype'],
                mintmelt=elements_spec[element_name]['mintmelt'],
                maxtmelt=elements_spec[element_name]['maxtmelt'],
                maxreplen=elements_spec[element_name]['maxreplen'],
                primercol=element_name,
                outfile=None,
                pairedcol=pairedcol,
                leftcontext=leftcontext,
                rightcontext=rightcontext,
                exmotifs=elements_spec[element_name]['exmotifs'],
                background=background,
                verbose=True)
            # print(dataframe)
            # print(dataframe.columns)
            for col in element_names:
                if not col in dataframe.columns:
                    # print(col)
                    dataframe[col] = '-'
            # print(dataframe.columns)
            dataframe = dataframe[element_names]
            # print(dataframe)
            stats_dict['oligopool'][element_name] = stats
            status = True
            if stats['status'] is False:
                incomplete = True
                break
        except Exception as E:
            raise E
            dataframe = None
            stats_dict['oligopool'][element_name] = {
                'status'  : False,
                'basis'   : 'infeasible',
                'step'    : 0,
                'stepname': 'parsing-input',
                'vars'    : {},
                'warns'   : {}}
            status = False
            break

    # Return Results
    return status, incomplete, dataframe, stats_dict

def insert_barcodes(
    dataframe,
    elements_spec,
    stats_dict):

    # Process All Variants
    status = False
    incomplete = False
    for element_name in elements_spec:
        if elements_spec[element_name]['type'] == 'barcode':
            try:
                if element_name in dataframe.columns:
                    dataframe.drop(element_name, inplace=True, axis=1)
                dataframe, stats = barcode.barcode(
                    indata=dataframe,
                    oligolimit=elements_spec[element_name]['oligolimit'],
                    barcodelen=elements_spec[element_name]['barcodelen'],
                    minhdist=elements_spec[element_name]['minhdist'],
                    maxreplen=elements_spec[element_name]['maxreplen'],
                    barcodecol=element_name,
                    outfile=None,
                    barcodetype=elements_spec[element_name]['barcodetype'],
                    leftcontext=elements_spec[element_name]['leftcontext'],
                    rightcontext=elements_spec[element_name]['rightcontext'],
                    exmotifs=elements_spec[element_name]['exmotifs'])
                stats_dict['oligopool'][element_name] = stats
                status = True
                if stats['status'] is False:
                    incomplete = True
                    break
            except Exception as E:
                raise E
                dataframe = None
                stats_dict['oligopool'][element_name] = {
                    'status'  : False,
                    'basis'   : 'infeasible',
                    'step'    : 0,
                    'stepname': 'parsing-input',
                    'vars'    : {},
                    'warns'   : {}}
                status = False
                break

    # Return Results
    return status, incomplete, dataframe, stats_dict

def insert_spacers(
    dataframe,
    elements_spec,
    stats_dict):

    # Process All Variants
    status = False
    incomplete = False
    for element_name in elements_spec:
        if elements_spec[element_name]['type'] == 'spacer':
            try:
                if element_name in dataframe.columns:
                    dataframe.drop(element_name, inplace=True, axis=1)
                dataframe, stats = spacer.spacer(
                    indata=dataframe,
                    oligolimit=elements_spec[element_name]['oligolimit'],
                    spacercol=element_name,
                    outfile=None,
                    spacerlen=elements_spec[element_name]['spacerlen'],
                    leftcontext=elements_spec[element_name]['leftcontext'],
                    rightcontext=elements_spec[element_name]['rightcontext'],
                    exmotifs=elements_spec[element_name]['exmotifs'])
                stats_dict['oligopool'][element_name] = stats
                status = True
                if stats['status'] is False:
                    incomplete = True
                    break
            except:
                dataframe = None
                stats_dict['oligopool'][element_name] = {
                    'status'  : False,
                    'basis'   : 'infeasible',
                    'step'    : 0,
                    'stepname': 'parsing-input',
                    'vars'    : {},
                    'warns'   : {}}
                status = False
                break

    # Return Results
    return status, incomplete, dataframe, stats_dict

def split_oligos(
    dataframe,
    split_spec,
    stats_dict):

    # Process All Variants
    status = False
    incomplete = False
    try:
        dataframe, stats = split.split(
            indata=dataframe,
            splitlimit=split_spec['splitlimit'],
            mintmelt=split_spec['mintmelt'],
            minhdist=split_spec['minhdist'],
            minoverlap=split_spec['minoverlap'],
            maxoverlap=split_spec['maxoverlap'],
            outfile=None)
        stats_dict['split'] = stats
        status = True
        if stats['status'] is False:
            incomplete = True
    except:
        dataframe = None
        stats_dict['split'] = {
            'status'  : False,
            'basis'   : 'infeasible',
            'step'    : 0,
            'stepname': 'parsing-input',
            'vars'    : {},
            'warns'   : {}}
        status = False

    # Return Results
    return status, incomplete, dataframe, stats_dict

def pad_oligos(
    dataframe,
    padding_spec,
    stats_dict):

    # Process All Variants
    status = False
    incomplete = False
    paddedframes = []
    stats_dict['padding'] = {}
    for splitcol in dataframe.columns:
        try:
            paddedframe, stats = pad.pad(
                indata=dataframe,
                splitcol=splitcol,
                typeIIS=padding_spec['typeIIS'],
                oligolimit=padding_spec['oligolimit'],
                mintmelt=padding_spec['mintmelt'],
                maxtmelt=padding_spec['maxtmelt'],
                maxreplen=padding_spec['maxreplen'],
                outfile=None)
            stats_dict['padding'][splitcol] = stats
            status = True
            if stats['status'] is False:
                incomplete = True
                break
        except:
            paddedframe = None
            stats_dict['padding'][splitcol] = {
                'status'  : False,
                'basis'   : 'infeasible',
                'step'    : 0,
                'stepname': 'parsing-input',
                'vars'    : {},
                'warns'   : {}}
            status = False
            break
        paddedframes.append(paddedframe)

    # Return Results
    return status, incomplete, paddedframes, stats_dict

def designparser(
    pool_size,
    element_names,
    elements_spec,
    background_spec,
    padding_spec,
    split_spec):
    '''
    Oligopool Calculator Design Mode Web Interface
    Parser Logic.

    :: poolsize
       type - integer
       desc - an integer equal to the number of unique
              variants in the oligopool
        e.g.  pool_size = 4350

    :: element_names
       type - list
       desc - a list of names of each oligopool element
              in the order of desired adjacency
        e.g.  element_names = [
            'Primer1',
            'Cut1',
            'Sequence',
            'Barcode',
            'Primer2',
            'Cut2',
            'Primer3',
            'Filler']

    :: elements_spec
       type - dict
       desc - a dictionary of various element design params
       e.g.   elements_spec={
            'Primer1': {
                        'type': 'primer',
                  'oligolimit': 170,
                   'primerseq': 'NNNNNNNNNNNNNNNNNNNNNN',
                  'primertype': 0,
                    'mintmelt': 53,
                    'maxtmelt': 55,
                   'maxreplen': 10,
                'pairedprimer': 'Primer3',
                 'leftcontext': None,
                'rightcontext': 'Cut1',
                    'exmotifs': ['GGATCC', 'TCTAGA'] + ['CCCCC', 'AAAAA', 'TTTTT', 'GGGGG']},

            'Cut1': {
                        'type': 'motif',
                  'oligolimit': 194,
                    'motifseq': 'NNNGGATCCNNN',
                 'leftcontext': 'Primer1',
                'rightcontext': 'Promoter',
                    'exmotifs': ['GGATCC', 'TCTAGA'] + ['CCCCC', 'AAAAA', 'TTTTT', 'GGGGG']},

            'Promoter': {
                     'type': 'variant',
                'sequences': plist},

            'Barcode': {
                        'type': 'barcode',
                  'oligolimit': 170,
                  'barcodelen': 15,
                    'minhdist': 3,
                   'maxreplen': 5,
                 'barcodetype': 1,
                    'exmotifs': ['GGATCC', 'TCTAGA'] + ['CCCCC', 'AAAAA', 'TTTTT', 'GGGGG'],
                 'leftcontext': 'Promoter',
                'rightcontext': 'Primer2'},

            'Primer2': {
                        'type': 'primer',
                  'oligolimit': 170,
                   'primerseq': 'NNNNNNNNNNNNNNNNNNNNNN',
                  'primertype': 0,
                    'mintmelt': 53,
                    'maxtmelt': 55,
                   'maxreplen': 10,
                'pairedprimer': 'Primer3',
                 'leftcontext': 'Barcode',
                'rightcontext': 'Cut2',
                    'exmotifs': ['GGATCC', 'TCTAGA'] + ['CCCCC', 'AAAAA', 'TTTTT', 'GGGGG']},

            'Cut2': {
                        'type': 'motif',
                  'oligolimit': 194,
                    'motifseq': 'NNNTCTAGANNN',
                 'leftcontext': 'Primer2',
                'rightcontext': 'Primer3',
                    'exmotifs': ['GGATCC', 'TCTAGA'] + ['CCCCC', 'AAAAA', 'TTTTT', 'GGGGG']},

            'Primer3': {
                        'type': 'primer',
                  'oligolimit': 170,
                   'primerseq': 'NNNNNNNNNNNNNNNNNNNNNN',
                  'primertype': 1,
                    'mintmelt': 53,
                    'maxtmelt': 55,
                   'maxreplen': 10,
                'pairedprimer': 'Primer2',
                 'leftcontext': 'Cut2',
                'rightcontext': 'Filler',
                    'exmotifs': ['GGATCC', 'TCTAGA'] + ['CCCCC', 'AAAAA', 'TTTTT', 'GGGGG']},

            'Filler': {
                        'type': 'spacer',
                  'oligolimit': 170,
                   'spacerlen': None,
                 'leftcontext': 'Primer3',
                'rightcontext': None,
                    'exmotifs': ['GGATCC', 'TCTAGA'] + ['CCCCC', 'AAAAA', 'TTTTT', 'GGGGG']},
        },

    :: background_spec
       type - dict / None
       desc - a dictionary with background params
       e.g. background_spec = {
               'indata': [seq1, seq2, ..., seqN],
            'maxreplen': 15}
    '''

    # Book-keeping
    stats_dict = {
        'background': None,
         'oligopool': {},
             'split': None,
           'padding': None}
    output = None

    # Setup Temporary Workspace
    workspacename = str(uuid.uuid4())
    setup_workspace(
        workspacename=workspacename)

    # Build Initial DataFrame
    dataframe = setup_dataframe(
        pool_size=pool_size,
        element_names=element_names)
    print('\n[Oligopool Calculator: Design Mode - Initialized DataFrame]\n')
    print(dataframe)

    # Setup Background
    if ((background_spec is None) or \
        (len(background_spec) <= 0)):
        backgroundname = None
        print()
    elif 'backgroundname' in background_spec:
        backgroundname = background_spec['backgroundname']
    else:
        (status,
        backgroundname,
        stats_dict) = setup_background(
            workspacename=workspacename,
            background_spec=background_spec,
            stats_dict=stats_dict)
        if status is False:
            return {
                'step': 1,
                'success' : False,
                'desc': 'Step 1: Background Input Error',
                'output': output,
            'stats_dict': stats_dict}
        elif stats_dict['background']['status'] is False:
            return {
                'step': 1,
                'success' : False,
            'step_name': 'Step 1: Background Unsuccessful',
                'output': output,
            'stats_dict': stats_dict}

    # Insert Core Variants
    (status,
    dataframe,
    stats_dict) = insert_variants(
        dataframe=dataframe,
        elements_spec=elements_spec,
        stats_dict=stats_dict)
    if status is False:
        return {
              'step': 2,
              'success' : False,
         'step_name': 'Step 2: Variants Input Error',
            'output': output,
        'stats_dict': stats_dict}

    # Insert All Motifs
    (status,
    incomplete,
    dataframe,
    stats_dict) = insert_motifs(
        dataframe=dataframe,
        elements_spec=elements_spec,
        stats_dict=stats_dict)
    if status is False:
        return {
              'step': 3,
              'success' : False,
         'step_name': 'Step 3: Motif Input Error',
            'output': output,
        'stats_dict': stats_dict}
    if incomplete:
        return {
              'step': 3,
              'success' : False,
         'step_name': 'Step 3: Motif Infeasible',
            'output': output,
        'stats_dict': stats_dict}

    # Insert All Primers
    (status,
    incomplete,
    dataframe,
    stats_dict) = insert_primers(
        element_names=element_names,
        background=backgroundname,
        dataframe=dataframe,
        elements_spec=elements_spec,
        stats_dict=stats_dict)
    if status is False:
        return {
              'step': 4,
              'success' : False,
         'step_name': 'Step 4: Primer Input Error',
            'output': output,
        'stats_dict': stats_dict}
    if incomplete:
        return {
              'step': 4,
              'success' : False,
         'step_name': 'Step 4: Primer Infeasible',
            'output': output,
        'stats_dict': stats_dict}

    # Insert All Barcodes
    (status,
    incomplete,
    dataframe,
    stats_dict) = insert_barcodes(
        dataframe=dataframe,
        elements_spec=elements_spec,
        stats_dict=stats_dict)
    if status is False:
        return {
              'step': 5,
              'success' : False,
         'step_name': 'Step 5: Barcode Input Error',
            'output': output,
        'stats_dict': stats_dict}
    if incomplete:
        return {
              'step': 5,
              'success' : False,
         'step_name': 'Step 5: Barcode Infeasible',
            'output': output,
        'stats_dict': stats_dict}

    # Insert All Spacers
    (status,
    incomplete,
    dataframe,
    stats_dict) = insert_spacers(
        dataframe=dataframe,
        elements_spec=elements_spec,
        stats_dict=stats_dict)
    if status is False:
        return {
              'step': 6,
              'success' : False,
         'step_name': 'Step 6: Spacer Input Error',
            'output': output,
        'stats_dict': stats_dict}
    if incomplete:
        return {
              'step': 6,
              'success' : False,
         'step_name': 'Step 6: Spacer Infeasible',
            'output': output,
        'stats_dict': stats_dict}

    # Finalize DataFrame
    final_dataframe, _ = final.final(
        indata=dataframe,
        outfile=None,
        verbose=True)
    final_dataframe.drop(
        'OligoLength', inplace=True, axis=1)

    # Split Oligos
    if split_spec is not None and len(split_spec) > 0:
        (status,
        incomplete,
        split_dataframe,
        stats_dict) = split_oligos(
            dataframe=final_dataframe,
            split_spec=split_spec,
            stats_dict=stats_dict)
        if status is False:
            return {
                'step': 8,
                'success' : False,
             'step_name': 'Step 8: Split Input Error',
                'output': output,
            'stats_dict': stats_dict}
        if incomplete:
            return {
                'step': 8,
                'success' : False,
             'step_name': 'Step 8: Split Infeasible',
                'output': output,
            'stats_dict': stats_dict}
    else:
        split_dataframe = None

    # Pad Oligos
    if padding_spec is not None:
        if split_dataframe is None:
            df = final_dataframe
        else:
            df = split_dataframe
            
        (status,
        incomplete,
        padded_dataframes,
        stats_dict) = pad_oligos(
            dataframe=df,
            padding_spec=padding_spec,
            stats_dict=stats_dict)
        if status is False:
            return {
                'step': 9,
             'step_name': 'Step 9: Padding Input Error',
                'output': output,
            'stats_dict': stats_dict}
        if incomplete:
            return {
                'step': 9,
             'step_name': 'Step 9: Padding Infeasible',
                'output': output,
            'stats_dict': stats_dict}
        
        padded_complete = [df.to_dict(orient='index') for df in padded_dataframes]
        
    else:
        padded_complete = None

    # Remove Workspace
    utils.remove_directory(
        dirpath=workspacename)

    # Compute Ouput
    output = {
        'annotated_complete': final_dataframe.to_dict(orient='index'),
           'padded_complete': padded_complete,          
    }

    # Return Results
    return {  'success' : True,
              'step': 10,
         'step_name': 'Step 10: Design Complete',
            'output': output,
        'stats_dict': stats_dict}