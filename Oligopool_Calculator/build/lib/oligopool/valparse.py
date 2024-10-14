import os

import numbers         as nu
import multiprocessing as mp
import zipfile         as zf
import math            as mt
import collections     as cx

import numpy    as np
import pandas   as pd
import psutil   as pu

import vectordb as db
import utils    as ut


def get_infile_validity(
    infile,
    infile_suffix,
    infile_field,
    liner):
    '''
    Determine if a given infile exists and
    non-empty. Internal use only.

    :: infile
       type - string
       desc - an input file to check for
              existence and emptiness
    :: infile_suffix
       type - string
       desc - required infile suffix
    :: infile_field
       type - string
       desc - infile fieldname used in
              printing
    :: liner
       type - coroutine
       desc - dynamic printing
    '''

    # infile exists and non-empty?
    infile_status = ut.get_path_status(
        path=infile,
        suffix=infile_suffix,
        readable=True,
        writable=False,
        creatable=False)

    # infile is invalid
    if   infile_status == 0:
        liner.send(
            '{}: {} [INPUT TYPE IS INVALID]\n'.format(
                infile_field, infile))

    else:

        infile = ut.get_adjusted_path(
            path=infile,
            suffix=infile_suffix)

        if   infile_status == 1:
            liner.send(
                '{}: {} [READ PERMISSION DENIED]\n'.format(
                    infile_field, infile))

        elif infile_status == 3:
            liner.send(
                '{}: {} [FILE IS EMPTY]\n'.format(
                    infile_field, infile))

        elif infile_status == 'X':
            liner.send(
                '{}: {} [FILE IS SPECIAL]\n'.format(
                    infile_field, infile))

        elif 5 <= infile_status <= 8:
            liner.send(
                '{}: {} [FILE IS DIRECTORY]\n'.format(
                    infile_field, infile))

        elif infile_status == 9:
            liner.send(
                '{}: {} [FILE DOES NOT EXIST]\n'.format(
                    infile_field, infile))

    # infile valid
    return infile_status == 4

def get_indir_validity(
    indir,
    indir_suffix,
    indir_field,
    liner):
    '''
    Determine if a given indir exists and
    non-empty. Internal use only.

    :: indir
       type - string
       desc - an input directory to check for
              existence and emptiness
    :: indir_suffix
       type - string
       desc - required indir suffix
    :: indir_field
       type - string
       desc - indir fieldname used in
              printing
    :: liner
       type - coroutine
       desc - dynamic printing
    '''

    # indir exists and non-empty?
    indir_status = ut.get_path_status(
        path=indir,
        suffix=indir_suffix,
        readable=True,
        writable=False,
        creatable=False)

    # indir is invalid
    if   indir_status == 0:
        liner.send(
            '{}: {} [INPUT TYPE IS INVALID]\n'.format(
                indir_field, indir))

    else:

        indir = ut.get_adjusted_path(
            path=indir,
            suffix=indir_suffix)

        if   indir_status == 5:
            liner.send(
                '{}: {} [READ PERMISSION DENIED]\n'.format(
                    indir_field, indir))

        elif indir_status == 7:
            liner.send(
                '{}: {} [DIRECTORY IS EMPTY]\n'.format(
                    indir_field, indir))

        elif indir_status == 'X':
            liner.send(
                '{}: {} [DIRECTORY IS SPECIAL]\n'.format(
                    indir_field, indir))

        elif 1 <= indir_status <= 4:
            liner.send(
                '{}: {} [DIRECTORY IS FILE]\n'.format(
                    indir_field, indir))

        elif indir_status == 9:
            liner.send(
                '{}: {} [DIRECTORY DOES NOT EXIST]\n'.format(
                    indir_field, indir))

    # indir valid
    return indir_status == 8

def get_readfile_validity(
    readfile,
    readfile_field,
    paired_readfile,
    liner):
    '''
    Determine if a given readfile exists, is
    non-empty and of FastQ type? If provided,
    check if given readfile is a duplicate of
    paired_readfile? Internal use only.

    :: readfile
       type - string
       desc - path to FastQ file storing
              reads
    :: readfile_field
       type - string
       desc - readfile fieldname used in
              printing
    :: paired_readfile
       type - string / None
       desc - path to paired FastQ file
              storing reads
    :: liner
       type - coroutine
       desc - dynamic printing
    '''

    # readfile exists and non-empty?
    readfile_exists = get_infile_validity(
        infile=readfile,
        infile_suffix=None,
        infile_field=readfile_field,
        liner=liner)

    # readfile is of FastQ type?
    readfile_is_fastq = False
    if readfile_exists:

        # FastQ file read attempt
        try:
            next(ut.stream_fastq_engine(
                filepath=readfile))

        # Read unsuccesful
        except:
            liner.send(
                '{}: {} [INVALID FASTQ FILE]\n'.format(
                    readfile_field, readfile))

        # Read successful
        else:
            readfile_is_fastq = True

    # Pairs duplicate?
    readfile_duplicate = False
    if readfile_is_fastq:

        # Pair to compare?
        if not paired_readfile is None:
            if os.path.samefile(
                readfile,
                paired_readfile):
                liner.send(
                    '{}: {} [DUPLICATE OF R1 FILE]\n'.format(
                        readfile_field, readfile))
                readfile_duplicate = True

        # Pair comparison successful
        # or unnecessary
        if not readfile_duplicate:
            liner.send('{}: {}\n'.format(
                readfile_field, readfile))

    # Return readfile validity
    return all([
        readfile_exists,
        readfile_is_fastq,
        not readfile_duplicate])

def get_optional_readfile_validity(
    readfile,
    readfile_field,
    paired_readfile,
    liner):
    '''
    Determine if an optional readfile exists,
    is non-empty and of FastQ type? If provided,
    check if given readfile is a duplicate of
    paired_readfile? Internal use only.

    :: readfile
       type - string
       desc - path to FastQ file storing
              reads
    :: readfile_field
       type - string
       desc - readfile fieldname used in
              printing
    :: paired_readfile
       type - string / None
       desc - path to paired FastQ file
              storing reads
    :: liner
       type - coroutine
       desc - dynamic printing
    '''

    # readfile is None
    if readfile is None:
        liner.send(
            '{}: None Specified\n'.format(
                readfile_field))
        return True

    # Regular readfile Validation
    return get_readfile_validity(
        readfile=readfile,
        readfile_field=readfile_field,
        paired_readfile=paired_readfile,
        liner=liner)

def get_inzip_validity(
    inzip,
    inzip_suffix,
    inzip_field,
    inzip_type,
    liner):
    '''
    Determine if zipfile is valid.
    Internal use only.

    :: inzip
       type - string
       desc - path to compressed zipfile
    :: inzip_suffix
       type - string
       desc - zipfile suffix adjustment
    :: inzip_field
       type - string
       desc - zipfile fieldname used
              in printing
    :: inzip_type
       type - string
       desc - zipfile type name string
    :: liner
       type - coroutine
       desc - dynamic printing
    '''

    fname = inzip.split('/')[-1]
    liner.send(
        '{}: Loading {} ...'.format(
            inzip_field, fname))

    # zipfile exists and non-empty?
    inzip_exists = get_infile_validity(
        infile=inzip,
        infile_suffix=inzip_suffix,
        infile_field=inzip_field,
        liner=liner)

    # inzip is zipfile?
    inzip_is_zipfile = False
    if inzip_exists:
        inzip = ut.get_adjusted_path(
            path=inzip,
            suffix=inzip_suffix)
        inzip_is_zipfile = zf.is_zipfile(
            filename=inzip)
        if not inzip_is_zipfile:
            liner.send(
                '{}: {} [INVALID {} FILE]\n'.format(
                    inzip_field, inzip, inzip_type))

    # inzip is good?
    inzip_is_good = False
    archive = None
    if inzip_is_zipfile:
        archive = ut.get_archive(
            arcfile=inzip)
        inzip_is_good = True
        # inzip_is_good = archive.testzip() is None
        # if not inzip_is_good:
        #     liner.send(
        #         '{}: {} [CORRUPT {} FILE]\n'.format(
        #             inzip_field, inzip, inzip_type))

    # inzip valid?
    return (all([
        inzip_exists,
        inzip_is_zipfile,
        inzip_is_good]),
        inzip,
        archive)

def get_indexfile_validity(
    indexfile,
    indexfile_field,
    associated,
    liner):
    '''
    Determine if indexfile points to a
    valid zipfile with indexed objects.
    Internal use only.

    :: indexfile
       type - string
       desc - path to compressed zipfile
              storing prepared structures
              and data models
              (suffix='.oligopool.index')
    :: indexfile_field
       type - string
       desc - indexfile fieldname used in
              printing
    :: associated
       type - boolean
       desc - if True then indexfile is
              expected to have associates
    :: liner
       type - coroutine
       desc - dynamic printing
    '''

    # indexfile exists and non-empty?
    (indexfile_ok_format,
    indexfile,
    archive) = get_inzip_validity(
        inzip=indexfile,
        inzip_suffix='.oligopool.index',
        inzip_field=indexfile_field,
        inzip_type='INDEX',
        liner=liner)

    # indexfile content ok?
    indexfile_ok_content = False
    if indexfile_ok_format:
        try:
            indexed = set([
                'ID.map',
                'barcode.model',
                'meta.map'])
            if associated:
                indexed.add('associate.map')
            if not (indexed <= set(archive.namelist())):
                raise
            variantcount = ut.loaddict(
                archive=archive,
                dfile='meta.map')['variantcount']
        except:
            liner.send(
                '{}: {} [INVALID INDEX FILE]\n'.format(
                    indexfile_field, indexfile))
        else:
            liner.send(
                '{}: {} w/ {:,} Variant(s)\n'.format(
                    indexfile_field, indexfile, variantcount))
            indexfile_ok_content = True
        finally:
            archive.close()

    # indexfile valid?
    return all([
        indexfile_ok_format,
        indexfile_ok_content])

def get_indexfiles_validity(
    indexfiles,
    indexfiles_field,
    associated,
    liner):
    '''
    Determine if all indexfiles are valid.
    Internal use only.

    :: indexfiles
       type - string / list
       desc - path to compressed zipfile
              storing prepared structures
              and data models
              (suffix='.oligopool.index')
    :: indexfiles_field
       type - string
       desc - indexfile fieldname used in
              printing
    :: associated
       type - boolean
       desc - if True then indexfiles are
              expected to have associates
    :: liner
       type - coroutine
       desc - dynamic printing
    '''

    # Do we have a single indexfile?
    if isinstance(indexfiles, str):
        return get_indexfile_validity(
            indexfile=indexfiles,
            indexfile_field=indexfiles_field,
            associated=associated,
            liner=liner)

    # Are indexfiles iterable?
    indexfiles_iterable = False
    indexfile_store = cx.deque()
    try:
        for indexfile in indexfiles:
            indexfile_store.append(indexfile)
        else:
            indexfiles_iterable = True
    except:
        liner.send(
            '{}: {} [INPUT TYPE IS INVALID]\n'.format(
                indexfiles_field, indexfiles))

    # Are indexfiles valid?
    indexfiles_ok = True
    seen_indexfiles  = set()
    if indexfiles_iterable:

        # Prepare Spacing
        altspacing = ' '*len(indexfiles_field)

        # Show Header Update
        liner.send(
            '{}: {:,} Index Input(s)\n'.format(
                indexfiles_field, len(indexfile_store)))

        # Core Validation Loop
        while indexfile_store:

            # Fetch indexfile
            indexfile = indexfile_store.popleft()

            # Duplicate indexfile?
            if indexfile in seen_indexfiles:
                # Show Update
                liner.send(
                    '{}: {} [DUPLICATE INDEX FILE]\n'.format(
                        altspacing, ut.get_adjusted_path(
                            path=indexfile,
                            suffix='oligopool.index')))
                # Update Global Validity
                indexfiles_ok = indexfiles_ok and False
                # Next!
                continue

            # Record indexfile
            else:
                seen_indexfiles.add(indexfile)

            # Get Current Validity
            idxfile_valid = get_indexfile_validity(
                indexfile=indexfile,
                indexfile_field=altspacing,
                associated=associated,
                liner=liner)

            # Update Global Validity
            indexfiles_ok = indexfiles_ok and idxfile_valid

    # Return Results
    return indexfiles_iterable and indexfiles_ok

def get_parsed_packfile(
    packfile,
    packfile_field,
    liner):
    '''
    Determine if packfile points to a
    valid zipfile with .pack files.
    Also, return packcount.
    Internal use only.

    :: packfile
       type - string
       desc - path to compressed zipfile
              storing read packs
              (suffix='.oligopool.pack')
    :: packfile_field
       type - string
       desc - packfile fieldname used in
              printing
    :: liner
       type - coroutine
       desc - dynamic printing
    '''

    # packfile exists and non-empty?
    (packfile_ok_format,
    packfile,
    archive) = get_inzip_validity(
        inzip=packfile,
        inzip_suffix='.oligopool.pack',
        inzip_field=packfile_field,
        inzip_type='PACK',
        liner=liner)

    # packfile content ok and non-empty?
    packfile_ok_content = False
    packtype  = None
    packcount = 0
    packfile_nonempty = False
    if packfile_ok_format:
        try:

            for cpath in archive.namelist():
                if not (cpath.endswith('.pack') or \
                        cpath.endswith('.stat')):
                    raise
                packcount += 1

            packcount -= 1 # Stat File isn't Pack

            packstat = ut.loaddict(
                archive=archive,
                dfile='packing.stat')

            if packcount != packstat['packcount']:
                raise

            packfile_nonempty = packcount > 0

        except:
            # Invalid / corrupt packfile
            liner.send(
                '{}: {} [INVALID PACK FILE]\n'.format(
                    packfile_field, packfile))
        else:
            # Empty packfile
            if not packfile_nonempty:
                liner.send(
                    '{}: {} w/ {} Read Packs [EMPTY PACK FILE]\n'.format(
                        packfile_field, packfile, packcount))
            # Packfile has packs
            else:
                liner.send(
                    '{}: {} w/ {} Read Packs\n'.format(
                        packfile_field, packfile, packcount))
                packfile_ok_content = True
        finally:
            archive.close()

    # packfile valid?
    return (all([
        packfile_ok_format,
        packfile_ok_content,
        packfile_nonempty]),
        packcount)

def get_parsed_data_info(
    data,
    data_field,
    required_fields,
    liner):
    '''
    Determine if given data is a valid,
    non-empty CSV file or a DataFrame.
    Internal use only.

    :: data
       type - string / pd.DataFrame
       desc - path to CSV file or a pandas
              DataFrame storing information
    :: data_field
       type - string
       desc - data fieldname used in
              printing
    :: required_fields
       type - list / None
       desc - list of column names which
              must be present in data
    :: liner
       type - coroutine
       desc - dynamic printing
    '''

    # data flags
    data_type   = None
    data_is_df  = False
    data_is_csv = False

    # df flags
    df = None
    df_nonempty = False

    # data is a DataFrame?
    if isinstance(data, pd.DataFrame):
        df = data.copy().reset_index(
            drop=data.index.name is None)
        data_is_df = True

    # data is a valid CSV file?
    else:

        # data file exists and non-empty?
        data_is_file = get_infile_validity(
            infile=data,
            infile_suffix=None,
            infile_field=data_field,
            liner=liner)

        # data is a CSV file?
        if data_is_file:

            # CSV file read attempt
            try:
                # Parse full df
                df = pd.read_csv(
                    filepath_or_buffer=data,
                    sep=',',
                    header=0,
                    engine='c')

            # Read and/or Indexing Unsuccesful
            except:
                liner.send('{}: {} [INVALID CSV FILE]\n'.format(
                    data_field,
                    data))

            # Read succesful
            else:
                data_is_csv = True

    # df was extracted?
    df_extracted = data_is_df or data_is_csv

    # Update data and data_type
    if   data_is_df:
        data      = 'DataFrame'
        data_type = 'DATAFRAME'
    elif data_is_csv:
        data_type = 'CSV FILE'

    # df is non-empty?
    df_nonempty = False

    if df_extracted:

        # Compute emptiness
        df_nonempty = not df.empty and \
                      len(df.columns) > 1

        if not df_nonempty:
            liner.send(
                '{}: {} w/ {:,} Record(s) [{} IS EMPTY]\n'.format(
                    data_field,
                    data,
                    len(df.index),
                    data_type))

    # df columns have missing values?
    df_no_missing_vals = False

    if df_nonempty:

        # Check column-wise emptiness
        for col in df.columns:
            if df[col].isnull().values.any():
                liner.send(
                    '{}: {} w/ {:,} Record(s) [MISSING VALUES IN COLUMN=\'{}\']\n'.format(
                        data_field,
                        data,
                        len(df.index),
                        col))
                break

        # No missing values
        else:
            df_no_missing_vals = True

    # df indexible?
    df_indexible = False

    if df_no_missing_vals:

        # Try indexing on unique ID
        try:

            # Reset and (Re-)Index by ID
            # print(df.index.is_unique)
            df.set_index(
                keys='ID',
                inplace=True)

            # Assert ID keys are unique
            if not df.index.is_unique:
                raise Exception

            # Everything checked out
            df_indexible = True

        # Indexing unsuccessful
        except:

            # Unindexible df
            liner.send(
                '{}: {} w/ {:,} Record(s) [NON-UNIQUE OR MISSING COLUMN=\'ID\']\n'.format(
                    data_field,
                    data,
                    len(df.index)))

            # Indexing failed
            df_indexible = False

    # df contains required columns?
    df_contains_required_cols = False

    if df_indexible:

        # Required fields present?
        if not required_fields is None:

            # Track requirement fulfillment
            required_fields_found = 0

            # Get df columns
            present_cols = set([c.lower() for c in df.columns])
            present_cols.update((df.index.name.lower(),))

            # Loop through required fields
            for required_field in required_fields:

                # Required field in present df columns
                if required_field.lower() in present_cols:
                    required_fields_found += 1
                # Missing field!
                else:
                    break

            # Requirements fulfilled
            if required_fields_found == len(required_fields):
                df_contains_required_cols = True

            # Requirement failed
            else:
                liner.send(
                    '{}: {} w/ {:,} Record(s) [MISSING COLUMN=\'{}\']\n'.format(
                        data_field,
                        data,
                        len(df.index),
                        required_field))

        # No requirements specified
        else:
            df_contains_required_cols = True

    # Compute final validity
    df_valid = all([
        df_extracted,
        df_nonempty,
        df_no_missing_vals,
        df_indexible,
        df_contains_required_cols])

    # Return Results
    data_name = data
    return (df, data_name, df_valid)

def get_parsed_indata_info(
    indata,
    indata_field,
    required_fields,
    precheck,
    liner):
    '''
    Determine if indata consisting of DNA columns
    only is valid. Internal use only.

    :: indata
       type - string / pd.DataFrame
       desc - path to CSV file or a pandas
              DataFrame storing information
    :: indata_field
       type - string
       desc - indata fieldname used in
              printing
    :: required_fields
       type - list / None
       desc - list of column names which
              must be present in data
    :: precheck
       type - boolean
       desc - if False prints content description
              when validation successful too,
              otherwise this is a pre-check
    :: liner
       type - coroutine
       desc - dynamic printing
    '''

    # Is indata valid CSV or DataFrame?
    (df,
    data_name,
    df_valid) = get_parsed_data_info(
        data=indata,
        data_field=indata_field,
        required_fields=required_fields,
        liner=liner)

    # df columns contain DNA strings only?
    df_contains_DNA_only = False

    if df_valid:

        # Are all entries DNA strings?

        non_DNA_found = False

        # Loop through all columns
        for column in df.columns:

            # Loop through all entires in column
            for value in df[column]:

                # A Non-DNA entry?
                if not ut.is_DNA(seq=value):

                    # Set flag
                    non_DNA_found = True

                    # Terminate inner loop
                    break

            # We encountered Non-DNA entry
            if non_DNA_found:
                # Terminate outer loop
                break

        df_contains_DNA_only = not non_DNA_found

        # SHow update
        if not df_contains_DNA_only:
            liner.send(
                '{}: {} w/ {:,} Record(s) [NON-DNA VALUE=\'{}\' IN COLUMN=\'{}\']\n'.format(
                    indata_field,
                    data_name,
                    len(df.index),
                    value,
                    column))

    # Compute final validity
    df_valid = all([
        df_valid,
        df_contains_DNA_only])

    # Is df valid?
    if df_valid:

        # Uppercase all DNA Strings
        for col in df.columns:
            df[col] = df[col].str.upper()

        # Show update?
        if not precheck:
            liner.send(
                '{}: {} w/ {:,} Record(s)\n'.format(
                    indata_field,
                    data_name,
                    len(df.index)))

    else:
        # Erase df
        df = None

    # Return data validity
    if precheck:
        return (df, data_name, df_valid)
    return (df, df_valid)

def get_outfile_validity(
    outfile,
    outfile_suffix,
    outfile_field,
    liner):
    '''
    Determine if outfile points to an existing
    empty file or non-existent but creatable path?
    Internal use only.

    :: outfile
       type - string
       desc - an output file storing
              computed information
    :: outfile_suffix
       type - string
       desc - required outfile suffix
    :: outfile_field
       type - string
       desc - outfile fieldname used in
              printing
    :: liner
       type - coroutine
       desc - dynamic printing
    '''

    # Append suffix to path?

    # outfile is an existing empty file or
    # non-existent but creatable path?
    outfile_status = ut.get_path_status(
        path=outfile,
        suffix=outfile_suffix,
        readable=False,
        writable=True,
        creatable=True)

    outfile_valid = False

    # outfile is non-string type?
    if outfile_status == 0:
        liner.send(
            '{}: {} [INPUT TYPE IS INVALID]\n'.format(
                outfile_field, outfile))

    else:

        # Adjust outfile with suffix
        outfile = ut.get_adjusted_path(
            path=outfile,
            suffix=outfile_suffix)

        # outfile is invalid
        if outfile_status in (2, 11):
            liner.send(
                '{}: {} [WRITE PERMISSION DENIED]\n'.format(
                    outfile_field, outfile))

        elif outfile_status == 4:
            liner.send(
                '{}: {} [FILE ALREADY EXISTS]\n'.format(
                    outfile_field, outfile))

        elif outfile_status == 'X':
            liner.send(
                '{}: {} [FILE IS SPECIAL]\n'.format(
                    outfile_field, outfile))

        elif 6 <= outfile_status <= 8:
            liner.send(
                '{}: {} [FILE IS DIRECTORY]\n'.format(
                    outfile_field, outfile))

        # outfile is valid
        elif outfile_status in (3, 10):
            liner.send('{}: {}\n'.format(
                outfile_field, outfile))
            outfile_valid = True

    # Return outfile validity
    return outfile_valid

def get_outdir_validity(
    outdir,
    outdir_suffix,
    outdir_field,
    liner):
    '''
    Determine if outdir points to an existing
    empty directory or non-existent but creatable
    path? Internal use only.

    :: outdir
       type - string
       desc - an output directory storing
              computed information
    :: outdir_suffix
       type - string
       desc - required outdir suffix
    :: outdir_field
       type - string
       desc - outdir fieldname used in
              printing
    :: liner
       type - coroutine
       desc - dynamic printing
    '''

    # Append suffix to path?

    # outdir is an existing empty directory
    # or non-existent but creatable path?
    outdir_status = ut.get_path_status(
        path=outdir,
        suffix=outdir_suffix,
        readable=False,
        writable=True,
        creatable=True)

    outdir_valid = False

    # outdir is non-string type?
    if outdir_status == 0:
        liner.send(
            '{}: {} [INPUT TYPE IS INVALID]\n'.format(
                outdir_field, outdir))

    else:

        # Adjust outdir with suffix
        outdir = ut.get_adjusted_path(
            path=outdir,
            suffix=outdir_suffix)

        # outdir is invalid
        if outdir_status in (6, 11):
            liner.send(
                '{}: {} [WRITE PERMISSION DENIED]\n'.format(
                    outdir_field, outdir))

        elif outdir_status == 8:
            liner.send(
                '{}: {} [DIRECTORY ALREADY EXISTS]\n'.format(
                    outdir_field, outdir))

        elif outdir_status == 'X':
            liner.send(
                '{}: {} [DIRECTORY IS SPECIAL]\n'.format(
                    outdir_field, outdir))

        elif 1 <= outdir_status <= 4:
            liner.send(
                '{}: {} [DIRECTORY IS FILE]\n'.format(
                    outdir_field, outdir))

        # outdir is valid
        elif outdir_status in (7, 10):
            liner.send('{}: {}\n'.format(
                outdir_field, outdir))
            outdir_valid = True

    # Return outdir validity
    return outdir_valid

def get_outdf_validity(
    outdf,
    outdf_suffix,
    outdf_field,
    liner):
    '''
    Determine if outdf points to a valid outfile
    if specified. Internal use only.

    :: outdf
       type - string / None
       desc - output file storing
              computed information
    :: outdf_suffix
       type - string
       desc - required outdf suffix
    :: outfile_field
       type - string
       desc - dipath fieldname used in
              printing
    :: liner
       type - coroutine
       desc - dynamic printing
    '''

    # outdf is None?
    if outdf is None:
        liner.send(
            '{}: In-Memory DataFrame\n'.format(
                outdf_field))
        return True

    # outdf is valid outfile?
    return get_outfile_validity(
        outfile=outdf,
        outfile_suffix=outdf_suffix,
        outfile_field=outdf_field,
        liner=liner)

def get_parsed_column_info(
    col,
    df,
    col_field,
    col_desc,
    col_type,
    adjcol,
    adjval,
    iscontext,
    typecontext,
    liner):
    '''
    Determine if col is a valid column in/for df
    depending on col_type. Internal use only.

    :: col
       type - string / None
       desc - column name to parse
    :: df
       type - pd.DataFrame
       desc - DataFrame to store output in
              or get input from
    :: col_field
       type - string
       desc - col fieldname used in
              printing
    :: col_desc
       type - string
       desc - col description used in
              printing
    :: col_type
       type - string
       desc - if 0, treat the column as input,
              otherwise, treat as output column
    :: adjcol
       type - string / None
       desc - ensure col is adjacent to adjcol
    :: adjval
       type - integer / None
       desc - ensure adjcol adjacent by adjval
    :: iscontext
       type - boolean
       desc - if True treat column as context,
              otherwise, column is a singular
              and context free element
    :: typecontext
       type - integer / None
       desc - if 0, treat context as left,
              otherwse, treat as right
    :: liner
       type - coroutine
       desc - dynamic printing
    '''

    # Is col None?
    if col is None:
        # OK for Input Column
        if col_type == 0:
            liner.send(
                '{}: None Specified\n'.format(
                    col_field,
                    col_desc,
                    col))
            return (None, True)

    # Is col string?
    col_is_string = False

    if isinstance(col, str):
        col_is_string = True

    else:
        liner.send(
            '{}: {} \'{}\' [INPUT TYPE IS INVALID]\n'.format(
                col_field,
                col_desc,
                col))

    # Is col non-existent?
    col_existence_valid = False

    if col_is_string:

        # Does df exist? No
        if not isinstance(df, pd.DataFrame):

            # There's no df to check col in
            liner.send(
                '{}: {} \'{}\' [COLUMN EXISTENCE UNKNOWN]\n'.format(
                    col_field,
                    col_desc,
                    col))

            # Guess, that means we don't know
            # the status of column specified
            col_existence_valid = False

        # Does df exist? Yes
        else:

            # Get col existence
            (col_existence,
            col_idx) = ut.get_col_exist_idx(
                col=col,
                df=df)

            # Is col_existence valid?
            if  ((col_type == 0 and    #  Input Column
                  col_existence) or    #  Input Column
                 (col_type == 1 and    # Output Column
                  not col_existence)): # Output Column
                col_existence_valid = True
            else:
                col_existence_valid = False

            # Show update for invalidation
            if not col_existence_valid:

                # Input Column
                if col_type == 0:
                    liner.send(
                    '{}: {} \'{}\' [COLUMN DOES NOT EXIST]\n'.format(
                        col_field,
                        col_desc,
                        col))

                # Output Column
                else:
                    liner.send(
                    '{}: {} \'{}\' [COLUMN ALREADY EXISTS]\n'.format(
                        col_field,
                        col_desc,
                        col))

    # Is col adjacent?
    col_is_adjacent = False

    if col_type == 0       and \
       col_existence_valid and \
       not adjcol is None  and \
       not adjval is None:

        # Get adjcol existence and index
        (adj_existence,
        adj_idx) = ut.get_col_exist_idx(
            col=adjcol,
            df=df)

        # adjcol non-existent?
        if not adj_existence:
            liner.send(
                '{}: {} \'{}\' [ADJACENT COLUMN \'{}\' DOES NOT EXIST]\n'.format(
                    col_field,
                    col_desc,
                    col,
                    adjcol))

        # adjcol more than adjval away?
        elif adj_idx - col_idx != adjval:
            liner.send(
                '{}: {} \'{}\' [COLUMN \'{}\' NON-ADJACENT]\n'.format(
                    col_field,
                    col_desc,
                    col,
                    adjcol))

        # All conditions met!
        else:
            col_is_adjacent = True

    else:
        col_is_adjacent = True


    # Is col valid?
    col_valid = all([col_is_string,
        col_existence_valid,
        col_is_adjacent])

    if col_valid:
        liner.send(
            '{}: {} \'{}\'\n'.format(
                col_field,
                col_desc,
                col))

        # Output Column
        if col_type == 1:
            parsedcol = None

        # Input Column
        else:
            # Non-Context Extraction
            if not iscontext:
                parsedcol = df[col].str.upper()
            # Context Extraction
            else:
                if typecontext == 0: #  Left Context
                    parsedcol = df.loc[:, :col]
                else:                # Right Context
                    parsedcol = df.loc[:, col:]
                parsedcol = parsedcol.sum(
                    axis=1).str.upper().str.replace(
                        '-', '')
    else:
        parsedcol = None

    # Return based on col_type
    if col_type == 1:
        return col_valid
    else:
        return (parsedcol,
            col_valid)

def get_parsed_exseqs_info(
    exseqs,
    exseqs_field,
    exseqs_desc,
    df_field,
    required,
    liner):
    '''
    Determine if given excluded sequences are valid.
    Internal use only.

    :: exseqs
       type - iterable / string / pd.DataFrame / None
       desc - iterable of DNA strings to be excluded
              from barcodes and around edges; optionally
              a DataFrame or a path to a CSV file with
              such sequences
    :: exseqs_field
       type - string
       desc - exseqs fieldname used in printing
    :: exseqs_desc
       type - string
       desc - exseqs description used in printing
    :: df_field
       type - string
       desc - name of the column in DataFrame storing
              excluded sequences
    :: required
       type - boolean
       desc - if True then exseqs cannot be None
    :: liner
       type - coroutine
       desc - dynamic printing
    '''

    # Is exseqs None?
    if exseqs is None:
        liner.send(
            '{}: 0 {}\n'.format(
                exseqs_field,
                exseqs_desc,))
        # Required but missing
        if required:
            return (None, False)
        # Non required so OK
        else:
            return (None, True)

    # Is exseqs string?
    if isinstance(exseqs, str):

        # Is exseqs a single DNA string?
        if ut.is_DNA(seq=exseqs.upper()):
            liner.send(
                '{}: 1 {}\n'.format(
                    exseqs_field,
                    exseqs_desc,))
            return ([exseqs.upper()], True)

    # Is exseqs iterable?
    if not isinstance(exseqs, pd.DataFrame) and \
       not isinstance(exseqs, str):

        # Try extracting exseqs
        try:
            exseqs = list(map(lambda x: x.upper(), ut.get_uniques(
                iterable=exseqs,
                typer=list)))

        # Error during extraction
        except:
            liner.send(
                '{}: {} [INPUT TYPE IS INVALID]\n'.format(
                    exseqs_field,
                    exseqs))
            return (None, False)

        # Ensure all motifs are DNA strings
        for seq in exseqs:

            # Non-DNA element found!
            if not ut.is_DNA(seq=seq):
                liner.send(
                    '{}: {:,} {} [NON-DNA SEQUENCE=\'{}\']\n'.format(
                        exseqs_field,
                        len(exseqs),
                        exseqs_desc,
                        seq))
                return (None, False)

        # No error in extraction
        # All sequences are DNA strings
        liner.send(
            '{}: {:,} {}\n'.format(
                exseqs_field,
                len(exseqs),
                exseqs_desc))

        return (exseqs, True)

     # Is exseqs a CSV file or DataFrame?
    (df, _,
    df_valid) = get_parsed_indata_info(
        indata=exseqs,
        indata_field=exseqs_field,
        required_fields=('ID', df_field,),
        precheck=True,
        liner=liner)

    # Is exseqs df valid?
    if df_valid:
        exseqs = list(map(lambda x: x.upper(), ut.get_uniques(
            iterable=df[df_field].to_list(),
            typer=list)))
        liner.send(
            '{}: {:,} {}\n'.format(
                exseqs_field,
                len(exseqs),
                exseqs_desc,))
    else:
        exseqs = None

    # Return Results
    return (exseqs, df_valid)

def get_numeric_validity(
    numeric,
    numeric_field,
    numeric_pre_desc,
    numeric_post_desc,
    minval,
    maxval,
    precheck,
    liner):
    '''
    Determine if numeric is a Real number.
    Internal use only.

    :: numeric
       type - Real
       desc - number to validate
    :: numeric_field
       type - string
       desc - numeric fieldname used in
              printing
    :: numeric_pre_desc
       type - string
       desc - numeric pre-description
              used in printing
    :: numeric_post_desc
       type - string
       desc - numeric post-description
              used in printing
    :: minval
       type - Real
       desc - minimum allowed value
    :: maxval
       type - Real
       desc - maximum allowed value
    :: precheck
       type - boolean
       desc - if False prints numeric_desc
              when validation successful too,
              otherwise this is a pre-check
    :: liner
       type - coroutine
       desc - dynamic printing
    '''

    # numeric is real?
    numeric_valid = False
    if not isinstance(numeric, nu.Real):
        liner.send('{}:{}{}{} [INPUT TYPE IS INVALID]\n'.format(
            numeric_field, numeric_pre_desc, numeric, numeric_post_desc))

    # numeric less than lowerbound?
    elif numeric < minval:
        liner.send('{}:{}{:,}{} [INPUT VALUE SMALLER THAN {:,}]\n'.format(
            numeric_field, numeric_pre_desc, numeric, numeric_post_desc, minval))

    # numeric greater than upperbound?
    elif numeric > maxval:
        liner.send('{}:{}{:,}{} [INPUT VALUE LARGER THAN {:,}]\n'.format(
            numeric_field, numeric_pre_desc, numeric, numeric_post_desc, maxval))

    # All conditions met!
    else:
        if not precheck:
            liner.send('{}:{}{:,}{}\n'.format(
                numeric_field, numeric_pre_desc, numeric, numeric_post_desc))
        numeric_valid = True

    # Return numeric validity
    return numeric_valid

def get_optional_numeric_validity(
    numeric,
    numeric_field,
    numeric_pre_desc,
    numeric_post_desc,
    minval,
    maxval,
    precheck,
    liner):
    '''
    Determine if optionally provided
    numeric is a Real number.
    Internal use only.

    :: numeric
       type - Real
       desc - number to validate
    :: numeric_field
       type - string
       desc - numeric fieldname used in
              printing
    :: numeric_pre_desc
       type - string
       desc - numeric pre-description
              used in printing
    :: numeric_post_desc
       type - string
       desc - numeric post-description
              used in printing
    :: minval
       type - Real
       desc - minimum allowed value
    :: maxval
       type - Real
       desc - maximum allowed value
    :: precheck
       type - boolean
       desc - if False prints numeric_desc
              when validation successful too,
              otherwise this is a pre-check
    :: liner
       type - coroutine
       desc - dynamic printing
    '''

    # numeric is None
    if numeric is None:
        liner.send(
            '{}: None Specified\n'.format(
                numeric_field))
        return True

    # Regular numeric Validation
    return get_numeric_validity(
        numeric=numeric,
        numeric_field=numeric_field,
        numeric_pre_desc=numeric_pre_desc,
        numeric_post_desc=numeric_post_desc,
        minval=minval,
        maxval=maxval,
        precheck=precheck,
        liner=liner)

def get_parsed_spacerlen_info(
    spacerlen,
    spacerlen_field,
    df_field,
    oligolimit,
    oligolimit_valid,
    indf,
    indata_valid,
    liner):
    '''
    Determine if given spacer length(s) input are valid.
    Internal use only.

    :: spacerlen
       type - iterable / integer / pd.DataFrame / None
       desc - iterable of integers denoting the length
              of spacers; optionally a DataFrame or a
              path to a CSV file with such lengths
    :: spacerlen_field
       type - string
       desc - spacerlen fieldname used in printing
    :: df_field
       type - string
       desc - name of the column in DataFrame storing
              spacer lengths
    :: oligolimit
       type - integer
       desc - maximum oligo length after inserting
              designed spacers
    :: oligolimit_valid
       type - boolean
       desc - oligolimit parsing status
    :: indf
       type - pd.DataFrame / None
       desc - Associated input DataFrame
    :: indata_valid
       type - boolean
       desc - Associated input DataFrame validity
    :: liner
       type - coroutine
       desc - dynamic printing
    '''

    # Is spacerlen None?
    if spacerlen is None:
        liner.send(
            '{}: Computed from Oligo Length (Auto-Inferred)\n'.format(
                spacerlen_field))
        return (None, True)

    # Quick compute indexlen
    if indata_valid:
        lenindex = len(indf.index)
    else:
        lenindex = 1

    # Correct oligolimit
    if not oligolimit_valid:
        oligolimit = float('+inf')
    else:
        oligolimit = round(oligolimit)

    # Is spacerlen numeric?
    if isinstance(spacerlen, nu.Real):

        # Ensure spacerlen valid
        spacerlen_valid = get_numeric_validity(
            numeric=spacerlen,
            numeric_field=spacerlen_field,
            numeric_pre_desc=' Exact1y ',
            numeric_post_desc=' Base Pair(s)',
            minval=1,
            maxval=oligolimit,
            precheck=False,
            liner=liner)

        # Build spacerlen
        if spacerlen_valid:
            spacerlen = np.zeros(
                lenindex, dtype=np.int64) + round(spacerlen)
        else:
            spacerlen = None

        # Return Results
        return (spacerlen, spacerlen_valid)

    # Is spacerlen extracted?
    spacerlen_extracted = False

    # Is spacerlen iterable?
    if not isinstance(spacerlen, pd.DataFrame) and \
       not isinstance(spacerlen, str) and \
       not isinstance(spacerlen, nu.Real):

        # Try extracting spacerlen
        try:
            spacerlen = list(n for n in spacerlen)

        # Error during extraction
        except:
            # Handled in sink below:
            # 'if not spacerlen_extracted: ...'
            spacerlen_extracted = False

        # Iteration successful
        else:
            spacerlen_extracted = True

        # Do we have enough spacers?
        if indata_valid:
            if len(spacerlen) != len(indf.index):

                # Compute category of mismatch
                catm = ['FEWER', 'MORE'][len(spacerlen) > len(indf.index)]

                # Show Update
                liner.send(
                    '{}: Iterable of {:,} Record(s) [{} THAN {:,} SPACERS SPECIFIED]\n'.format(
                        spacerlen_field,
                        len(spacerlen),
                        catm,
                        len(indf.index)))

                # Return Results
                return (None, False)

    # Is spacerlen a CSV file or DataFrame?
    else:

        # Compute spacerlen data validity
        (df,
        data_name,
        df_valid) = get_parsed_data_info(
            data=spacerlen,
            data_field=spacerlen_field,
            required_fields=('ID', df_field,),
            liner=liner)

        # Is spacerlen df valid?
        if df_valid:

            # Are the indexes matching?
            idx_match = False
            if indata_valid:

                # Matching IDs?
                idx_match = len(df.index)    == len(indf.index) and \
                            sorted(df.index) == sorted(indf.index)

                # Everything OK!
                if idx_match:
                    df = df.reindex(indf.index)
                    spacerlen = list(n for n in df[df_field])
                    spacerlen_extracted = True

            # No basis for matching, so show indifference
            # here, return and fail later on
            else:
                liner.send(
                    '{}: {} w/ {:,} Record(s)\n'.format(
                        spacerlen_field,
                        data_name,
                        len(df.index)))
                return (None, False)

            # Indexes don't match
            if not idx_match:
                liner.send(
                    '{}: {} w/ {:,} Record(s) [COLUMN=\'ID\' DOES NOT MATCH INPUT DATA]\n'.format(
                        spacerlen_field,
                        data_name,
                        len(df.index)))
                return (None, False)

        else:
            # Note: We've already shown how
            # invalidity impacts spacerlen,
            # so we're directly returning
            return (None, False)

    # spacerlen not extracted?
    if not spacerlen_extracted:
        # Note: if spacerlen is neither numeric, pd.DataFrame
        # iterable, nor a string path, then it should sink here
        liner.send(
            '{}: {} [INPUT TYPE IS INVALID]\n'.format(
                spacerlen_field,
                spacerlen))
        return (None, False)

    # Finalize spacerlen
    minimum = float('+inf')
    maximum = float('-inf')
    finalspacerlen = []

    for sl in spacerlen:

        # Ensure sl valid
        sl_valid = get_numeric_validity(
            numeric=sl,
            numeric_field=spacerlen_field,
            numeric_pre_desc=' Found an Entry for ',
            numeric_post_desc=' Base Pair(s)',
            minval=0,
            maxval=oligolimit,
            precheck=True,
            liner=liner)

        # Invalid sl
        if not sl_valid:
            return (None, False)
        # So far valid
        else:
            sl = round(sl)
            finalspacerlen.append(sl)
            minimum = min(minimum, sl)
            maximum = max(maximum, sl)

    # No errors whatsoever
    spacerlen = finalspacerlen
    spacerlen_valid = True

    if minimum == maximum:

        # We don't tolerate all zeros
        if minimum == 0:
            liner.send(
                '{}: Exactly 0 Base Pair(s) [ALL SPACERS HAVE ZERO LENGTH]\n'.format(
                    spacerlen_field,
                    minimum))
            return (None, False)

        # At least some entries require spacers
        else:
            liner.send(
                '{}: Exactly {:,} Base Pair(s)\n'.format(
                    spacerlen_field,
                    minimum))
    else:
        liner.send(
            '{}: {:,} to {:,} Base Pair(s)\n'.format(
                spacerlen_field,
                minimum,
                maximum))

    # Pack spacerlen
    spacerlen = np.array(spacerlen)

    # Return Results
    return (spacerlen, True)

def get_parsed_associatedata_info(
    associatedata,
    associatedata_field,
    required_fields,
    bardf,
    barcodedata_valid,
    liner):
    '''
    Determine if associatedata is valid and the
    ID is consistant with bardf.
    Internal use only.

    :: associatedata
       type - string / pd.DataFrame
       desc - path to CSV file or a pandas
              DataFrame storing variant
              information
    :: associatedata_field
       type - string
       desc - associatedata fieldname used in
              printing
    :: required_fields
       type - list / None
       desc - list of column names which
              must be present in data
    :: bardf
       type - pd.DataFrame / None
       desc - pandas DataFrame containing
              barcode information
    :: barcodedata_valid
       type - boolean
       desc - if True indicates that shared
              ID is to be computed, otherwise
              ID matching is skipped
    :: liner
       type - coroutine
       desc - dynamic printing
    '''

    # Is associatedata None?
    if associatedata is None:
        liner.send(
            '{}: None Specified\n'.format(
                associatedata_field))
        return (None, True) # Because variant is optional

    # Is associatedata valid?
    (assdf,
    data_name,
    associatedata_valid) = get_parsed_indata_info(
        indata=associatedata,
        indata_field=associatedata_field,
        required_fields=required_fields,
        precheck=True,
        liner=liner)

    # Does assdf share ID with bardf?
    associatedata_idx_match = False
    if associatedata_valid:

        # Are the indexes matching?
        idx_match = False
        if barcodedata_valid:

            # Matching IDs?
            idx_match = len(assdf.index)    == len(bardf.index) and \
                        sorted(assdf.index) == sorted(bardf.index)

            # Everything OK!
            if idx_match:
                assdf = assdf.reindex(bardf.index)
                associatedata_idx_match = True # same as idx_match

            # Indexes don't match
            else:
                liner.send(
                    '{}: {} w/ {:,} Record(s) [COLUMN=\'ID\' DOES NOT MATCH BARCODE DATA]\n'.format(
                        associatedata_field,
                        data_name,
                        len(assdf.index)))

        # No basis for matching, so show indifference
        # here, return and fail later on
        else:
            associatedata_idx_match = True # Technically, we can't complain!

    # Compute final validity
    assdf_valid = all([
        associatedata_valid,
        associatedata_idx_match])

    # Is df valid?
    if assdf_valid:
        liner.send(
            '{}: {} w/ {:,} Record(s)\n'.format(
                associatedata_field,
                data_name,
                len(assdf.index)))
    else:
        # Erase df
        assdf = None

    # Return data validity
    return (assdf, assdf_valid)

def get_categorical_validity(
    category,
    category_field,
    category_pre_desc,
    category_post_desc,
    category_dict,
    liner):
    '''
    Determine if category is valid with
    respect to category_dict.
    Internal use only.

    :: category
       type - Real
       desc - category to validate
    :: category_field
       type - string
       desc - category fieldname used in
              printing
    :: category_pre_desc
       type - string
       desc - category pre-description
              used in printing
    :: category_post_desc
       type - string
       desc - category post-description
              used in printing
    :: category_dict
       type - dict
       desc - category description used
              in printing
    :: liner
       type - coroutine
       desc - dynamic printing
    '''

    # Is category numeric?
    cat_is_numeric = get_numeric_validity(
        numeric=category,
        numeric_field=category_field,
        numeric_pre_desc=category_pre_desc,
        numeric_post_desc=category_post_desc,
        minval=min(category_dict.keys()),
        maxval=max(category_dict.keys()),
        precheck=True,
        liner=liner)

    # Is category present?
    cat_is_present = False
    if cat_is_numeric:
        if not category in category_dict:
            liner.send('{}:{}{}{} [INPUT VALUE IS INVALID]\n'.format(
                category_field,
                category_pre_desc,
                category,
                category_post_desc))
        else:
            cat_is_present = True

    # Compute validity
    cat_valid = cat_is_numeric and cat_is_present

    # Show update
    if cat_valid:
        liner.send('{}:{}{}{}\n'.format(
            category_field,
            category_pre_desc,
            category_dict[category],
            category_post_desc))

    # Return catgory validity
    return cat_valid

def get_optional_categorical_validity(
    category,
    category_field,
    category_pre_desc,
    category_post_desc,
    category_dict,
    liner):
    '''
    Determine if optional category is valid
    with respect to category_dict.
    Internal use only.

    :: category
       type - Real
       desc - category to validate
    :: category_field
       type - string
       desc - category fieldname used in
              printing
    :: category_pre_desc
       type - string
       desc - category pre-description
              used in printing
    :: category_post_desc
       type - string
       desc - category post-description
              used in printing
    :: category_dict
       type - dict
       desc - category description used
              in printing
    :: liner
       type - coroutine
       desc - dynamic printing
    '''

    # category is None
    if category is None:
        liner.send(
            '{}: None Specified\n'.format(
                category_field))
        return True

    # Regular numeric Validation
    return get_categorical_validity(
        category=category,
        category_field=category_field,
        category_pre_desc=category_pre_desc,
        category_post_desc=category_post_desc,
        category_dict=category_dict,
        liner=liner)

def get_parsed_typeIIS_info(
    typeIIS,
    typeIIS_field,
    liner):
    '''
    Determine if typeIIS system selected
    is valid. Internal use only.

    :: typeIIS
       type - string
       desc - name of the TypeIIS enzyme
              selected for excision
    :: typeIIS_field
       type - string
       desc - typeIIS fieldname used in
              printing
    :: liner
       type - coroutine
       desc - dynamic printing
    '''

    # Is typeIIs string?
    typeIIS_is_string = False
    if not isinstance(typeIIS, str):
        liner.send('{}: Enzyme \'{}\' [INPUT TYPE IS INVALID]\n'.format(
                typeIIS_field,
                typeIIS))
    else:
        typeIIS_is_string = True

    # Is typeIIS known?
    typeIIS_is_known = False
    if typeIIS_is_string:
        typeIIS_ = typeIIS.lower()
        if not typeIIS_ in ut.typeIIS_dict:
            liner.send('{}: Enzyme \'{}\' [UNSUPPORTED ENZYME]\n'.format(
                    typeIIS_field,
                    typeIIS))
        else:
            typeIIS_is_known = True

    # Compute validity
    typeIIS_valid = typeIIS_is_string and typeIIS_is_known

    # Show update
    if typeIIS_valid:
        liner.send('{}: Enzyme \'{}\' Recognizing Motif \'{}\'\n'.format(
            typeIIS_field,
            ut.typeIIS_dict[typeIIS_][0],
            ut.typeIIS_dict[typeIIS_][1]))

    # Return typeIIS validity
    if typeIIS_valid:
        typeIISname = ut.typeIIS_dict[typeIIS_][0]
        typeIIS = ut.typeIIS_dict[typeIIS_][1] + \
                 ('N' * ut.typeIIS_dict[typeIIS_][2])
        return typeIIS, typeIISname, typeIIS_valid
    else:
        return None, None, typeIIS_valid

def get_seqconstr_validity(
    seqconstr,
    seqconstr_field,
    minlenval,
    element,
    liner):
    '''
    Determine if seqconstr is a valid
    degenerate sequence constraint.
    Internal use only.

    :: seqconstr
       type - string
       desc - primer sequence constraint
    :: seqconstr_field
       type - string
       desc - seqconstr fieldname used
              in printing
    :: minlenval
       type - integer
       desc - minimum length allowed for
              sequence constraint
    :: element
       type - string
       desc - element being designed
    :: liner
       type - coroutine
       desc - dynamic printing
    '''

    # Is seqconstr string?
    seqconstr_is_string = False
    if not isinstance(seqconstr, str):
        liner.send(
            '{}: {} [INPUT TYPE IS INVALID]\n'.format(
                seqconstr_field,
                seqconstr))
    else:
        seqconstr_is_string = True

    # Is seqconstr degenerate DNA string?
    seqconstr_is_ddna = False
    if seqconstr_is_string:
        if not ut.is_DNA(
            seq=seqconstr,
            dna_alpha=ut.ddna_alpha):
            liner.send(
                '{}: {} [NON-IUPAC VALUE]\n'.format(
                    seqconstr_field,
                    seqconstr))
        else:
            seqconstr_is_ddna = True

    # Is seqconstr long enough?
    seqconstr_is_long = False
    if seqconstr_is_ddna:
        if len(seqconstr) < minlenval:
            liner.send(
                '{}: A {:,} Base Pair IUPAC Constraint [{} SHORTER THAN {:,} BASE PAIRS]\n'.format(
                    seqconstr_field,
                    len(seqconstr),
                    element,
                    minlenval))
        else:
            seqconstr_is_long = True

    # Compute final validity
    seqconstr_valid = all([
        seqconstr_is_string,
        seqconstr_is_ddna,
        seqconstr_is_long])

    # Show update
    if seqconstr_valid:
        liner.send(
            '{}: A {:,} Base Pair IUPAC Constraint\n'.format(
                seqconstr_field,
                len(seqconstr)))

    # Return validity
    return seqconstr_valid

def get_constantcol_validity(
    constantcol,
    constantcol_field,
    df,
    liner):
    '''
    Determine if constantcol is valid.
    Internal use only.

    :: constantcol
       type - string / None
       desc - the name of column in df where
              constantcol is stored
    :: constantcol_field
       type - string
       desc - constantcol fieldname used
              in printing
    :: df
       type - pd.DataFrame
       desc - input DataFrame where constant
              is present
    :: liner
       type - coroutine
       desc - dynamic printing
    '''

    # Is constantcol None?
    if constantcol is None:
        liner.send(
            '{}: None Specified\n'.format(
                constantcol_field))
        return None, True

    # Is constantcol a string?
    if isinstance(constantcol, str):

        # Is constantcol a column in df?
        if constantcol in df.columns:

            # Extract unique candidates
            uniques = ut.get_uniques(
                iterable=df[constantcol],
                typer=tuple)

            # Too many constant primer candidates?
            if len(uniques) > 1:
                liner.send(
                    '{}: Input from Column \'{}\' [NON-UNIQUE COLUMN=\'{}\']\n'.format(
                        constantcol_field,
                        constantcol,
                        constantcol))
                return None, False

            # Unique constant primer
            else:
                liner.send(
                    '{}: A {:,} Base Pair DNA Sequence\n'.format(
                        constantcol_field,
                        len(uniques[0])))
                return uniques[0], True

        # Nothing matches
        else:
            liner.send(
                '{}: Input from Column \'{}\' [MISSING COLUMN=\'{}\']\n'.format(
                    constantcol_field,
                    constantcol,
                    constantcol))
            return None, False

    # Non-string constantcol
    else:
        liner.send(
            '{}: Input from Column \'{}\' [INPUT TYPE IS INVALID]\n'.format(
                constantcol_field,
                constantcol))
        return None, False

def get_parsed_range_info(
    minval,
    maxval,
    range_field,
    range_unit,
    range_min,
    range_max,
    liner):
    '''
    Determine if range information
    is valid. Internal use only.

    :: minval
       type - nu.Real
       desc - minimum range value
    :: maxval
       type - nu.Real
       desc - maximum range value
    :: range_field
       type - string
       desc - range fieldname used in
              printing
    :: range_unit
       type - string
       desc - range value unit used in
              printing
    :: range_min
       type - nu.Real
       desc - allowed range lowerbound
    :: range_max
       type - nu.Real
       desc - allowed range upperbound
    :: liner
       type - coroutine
       desc - dynamic printing
    '''

    # Define value range
    vrange = (minval, maxval)

    # Is vrange numeric?
    vrange_is_numeric = False
    if  not isinstance(minval, nu.Real) or \
        not isinstance(maxval, nu.Real):
        liner.send(
            '{}: {} to {} {} [INPUT TYPE IS INVALID]\n'.format(
                range_field,
                minval,
                maxval,
                range_unit))
    else:
        vrange_is_numeric = True

    # Sort vrange
    if vrange_is_numeric:
        vrange = tuple([
            min(minval, maxval),
            max(minval, maxval)])

    # Is vrange in practical range?
    vrange_min_is_practical = False
    if vrange_is_numeric:
        if vrange[0] < range_min:
            liner.send(
                '{}: {} to {} {} [MINIMUM VALUE SMALLER THAN {} {}]\n'.format(
                    range_field,
                    vrange[0],
                    vrange[1],
                    range_unit,
                    range_min,
                    range_unit))
        else:
            vrange_min_is_practical = True

    # Is maxtmelt in practical range?
    vrange_max_is_practical = False
    if vrange_min_is_practical:
        if vrange[1] > range_max:
            liner.send(
                '{}: {} to {} {} [MAXIMUM VALUE LARGER THAN {} {}]\n'.format(
                    range_field,
                    vrange[0],
                    vrange[1],
                    range_unit,
                    range_max,
                    range_unit))
        else:
            vrange_max_is_practical = True

    # Compute validity
    vrange_valid = all([
        vrange_is_numeric,
        vrange_min_is_practical,
        vrange_max_is_practical])

    # Show update
    if vrange_valid:
        liner.send(
            '{}: {} to {} {}\n'.format(
                range_field,
                vrange[0],
                vrange[1],
                range_unit))

    # Return validity
    return (
        vrange[0],
        vrange[1],
        vrange_valid)

def get_parsed_background(
    background,
    background_field,
    liner):
    '''
    Determine if background is valid.
    Internal use only.

    :: background
       type - string / db.vectorDB / None
       desc - path to background storage,
              or a vectorDB instance
    :: background_field
       type - string
       desc - background fieldname used in
              printing
    :: liner
       type - coroutine
       desc - dynamic printing
    '''

    # Is background None?
    if background is None:
        liner.send(
            '{}: None Specified\n'.format(
                background_field))
        return None, True

    # Is background a vectorDB instance?
    if isinstance(background, db.vectorDB):
        liner.send(
            '{}: Contains {:,} Unique {}-mers\n'.format(
                background_field,
                len(background),
                background.K))
        return background, True

    # Is background a string?
    if isinstance(background, str):

        # Adjust background path
        indir = ut.get_adjusted_path(
            path=background,
            suffix='.oligopool.background')

        # Does background exist?
        background_exists = get_indir_validity(
            indir=indir,
            indir_suffix=None,
            indir_field=background_field,
            liner=liner)

        # Non-existent background
        if not background_exists:
            return None, False

        # Is background a valid vectorDB storage?
        else:

            # Open path as vectorDB instance
            try:
                vDB = db.vectorDB(
                    path=indir,
                    maxreplen=None,
                    mode=1)

            # Invalid attempt
            except Exception as E:
                liner.send(
                    '{}: {} [INVALID OR PRE-OPENED BACKGROUND OBJECT]\n'.format(
                        background_field,
                        indir))
                return None, False

            # Valid attempt
            else:
                liner.send(
                    '{}: Contains {:,} Unique {}-mers\n'.format(
                        background_field,
                        len(vDB),
                        vDB.K))
                return vDB, True

    # Invalid input type
    else:
        liner.send(
            '{}: {} [INPUT TYPE IS INVALID]\n'.format(
                background_field,
                background))
        return None, False

def get_callback_validity(
    callback,
    callback_field,
    liner):
    '''
    Determine if callback function is valid.
    Internal use only.

    :: callback
       type - function / None
       desc - callback function specified
    :: background_field
       type - string
       desc - callback function fieldname
              used in printing
    :: liner
       type - coroutine
       desc - dynamic printing
    '''

    # Is callback None?
    if callback is None:
        liner.send(
            '{}: None Specified\n'.format(
                callback_field))
        return True

    # Is callback callable?
    if not callable(callback):
        liner.send(
            '{}: {} [INVALID FUNCTION]\n'.format(
                callback_field,
                callback))
        return False
    else:
        liner.send(
            '{}: Function w/ ID={}\n'.format(
                callback_field,
                id(callback)))
        return True

def get_errors_validity(
    errors,
    errors_field,
    errors_pre_desc,
    errors_post_desc,
    errors_base,
    indexfiles_valid,
    indexfiles,
    liner):
    '''
    Determine if numeric is a None or a
    positive Real. Internal use only.

    :: errors
       type - Real / None
       desc - error value to validate
    :: errors_field
       type - string
       desc - errors fieldname used in
              printing
    :: errors_pre_desc
       type - string
       desc - errors pre-description used
              in printing
    :: errors_post_desc
       type - string
       desc - errors post-description used
              in printing
    :: errors_base
       type - string
       desc - either 'A' or 'B'
    :: indexfile_valid
       type - boolean
       desc - True if indexfiles are valid,
              False otherwise
    :: indexfile
       type - iterable
       desc - filenames storing prepared
              indexes and models
    :: liner
       type - coroutine
       desc - dynamic printing
    '''

    # errors is positive real?
    error_numeric = get_numeric_validity(
        numeric=errors,
        numeric_field=errors_field,
        numeric_pre_desc=errors_pre_desc,
        numeric_post_desc=errors_post_desc,
        minval=float('-inf'),
        maxval=float('+inf'),
        precheck=True,
        liner=liner)

    # errors status?
    if error_numeric:
        errors = round(errors)
        # errors is default
        if errors < 0.:
            if indexfiles_valid:
                errors = float('-inf')
                for indexfile in indexfiles:
                    archive = zf.ZipFile(
                        file=indexfile)
                    metamap = ut.loaddict(
                        archive=archive,
                        dfile='meta.map')
                    if errors_base == 'A':
                        errors = max(
                            errors,
                            metamap['associatetvalmax'])
                    else:
                        errors = max(
                            errors,
                            metamap['barcodetval'])
                    errors = round(errors)
                    archive.close()
                liner.send(
                    '{}:{}{}{} (Auto-Inferred)\n'.format(
                        errors_field,
                        errors_pre_desc,
                        errors,
                        errors_post_desc))
            else:
                liner.send(
                    '{}: Indeterminable\n'.format(
                        errors_field))
        # errors is specified
        else:
            liner.send('{}:{}{}{}\n'.format(
                errors_field,
                errors_pre_desc,
                errors,
                errors_post_desc))

    # Return errors validity
    return (errors,
        error_numeric)

def get_parsed_core_info(
    ncores,
    core_field,
    default,
    offset,
    liner):
    '''
    Determine if ncores is a positive Real and
    valid for given function. Internal use only.

    :: ncores
       type - Real
       desc - number of cores to use for function
    :: core_field
       type - string
       desc - ncores fieldname used in printing
    :: default
       type - None / integer
       desc - default value to use when ncores
              is Real but logically invalid
    :: offset
       type - integer
       desc - maximum numbers of core to withold
              from allocation if all cores usable
    :: liner
       type - coroutine
       desc - dynamic printing
    '''

    # How many cores available in total?
    sys_cores = mp.cpu_count()

    # ncores is non-integer?
    if not isinstance(ncores, nu.Real):
        liner.send(
            '{}: Use {} out of {:,} Cores [INPUT TYPE IS INVALID]\n'.format(
                core_field, ncores, sys_cores))
        ncores_valid = False

    # ncores must be >= 0
    elif ncores < 0:
        liner.send(
            '{}: Use {:,} out of {:,} Cores [INPUT VALUE IS INVALID]\n'.format(
                core_field, ncores, sys_cores))
        ncores_valid = False

    # ncores is valid
    elif ncores >= 0:
        autoinferred = ''

        # ncores adjusted
        if  ncores > sys_cores:
            ncores = sys_cores
            autoinferred = '(Auto-Inferred)'

        # ncores defaulting
        elif ncores == 0:
            autoinferred = '(Auto-Inferred)'

            # Amdahl's Ironclad Law
            optncores = max(
                round(mt.sqrt(sys_cores)),
                round(mt.log2(sys_cores)))

            if default is None:
                ncores = optncores
            else:
                ncores = min(default, sys_cores)

        liner.send(
            '{}: Use {:,} out of {:,} Cores {}\n'.format(
                core_field,
                ncores,
                sys_cores,
                autoinferred))

        ncores = round(ncores)
        ncores_valid = True

    # Return adjusted ncores and validity
    adjcores = max(1, min(ncores, sys_cores-offset))
    return adjcores, ncores_valid

def get_parsed_memory_info(
    memlimit,
    memlimit_field,
    ncores,
    ncores_valid,
    liner):
    '''
    Determine if memlimit is a positive Real and
    valid for given function. Internal use only.

    :: memlimit
       type - Real
       desc - amount of memory to be allocated
              per core for a function
    :: core_field
       type - string
       desc - memlimit fieldname used in printing
    :: ncores
       type - integer
       desc - total number of cores used in the
              function
    :: ncores_valid
       type - boolean
       desc - if True ncores was parsed to be
              valid
    :: liner
       type - coroutine
       desc - dynamic printing
    '''

    # Is ncores invalid?
    if not ncores_valid:
        liner.send(
            '{}: Use {} GB RAM per Core\n'.format(
                memlimit_field, memlimit))
        return memlimit, False

    # How much memory available in total?
    sys_mem = np.floor(
        pu.virtual_memory().total / (10**9)) - 2.0

    # memlimit is non-integer?
    if not isinstance(memlimit, nu.Real):
        liner.send(
            '{}: Use {} GB RAM per Core [INPUT TYPE IS INVALID]\n'.format(
                memlimit_field, memlimit))
        memlimit_valid = False

    # memlimit must be >= 0
    elif memlimit < 0:
        liner.send(
            '{}: Use {:.2f} GB RAM per Core [INPUT VALUE IS INVALID]\n'.format(
                memlimit_field, memlimit))
        memlimit_valid = False

    # memlimit is valid
    elif memlimit >= 0:
        autoinferred = ''

        # Normalize to Core Count
        sys_mem /= ncores

        # memlimit adjusted
        if  memlimit > sys_mem or \
            memlimit == 0.:
            memlimit = sys_mem
            autoinferred = '(Auto-Inferred)'

        liner.send(
            '{}: Use {:.2f} GB RAM per Core {}\n'.format(
                memlimit_field,
                memlimit,
                autoinferred))

        memlimit_valid = True

    # Return adjusted memlimit and validity
    return memlimit, memlimit_valid
