import time as tt

import collections as cx
import random      as rn
import atexit      as ae

import multiprocessing as mp

import utils    as ut
import valparse as vp
import corepack as cp


def pack_engine(
    r1file,
    r2file,
    r1type,
    r2type,
    r1length,
    r2length,
    r1qual,
    r2qual,
    packdir,
    metaqueue,
    packqueue,
    packtype,
    assemblyparams,
    packsize,
    r1truncfile,
    r2truncfile,
    previousreads,
    scannedreads,
    ambiguousreads,
    shortreads,
    survivedreads,
    packedreads,
    packsbuilt,
    coreid,
    batchid,
    ncores,
    nactive,
    memlimit,
    launchtime,
    restart,
    shutdown,
    liner):
    '''
    Process r1file and r2file and pack read
    pairs by frequency. Internal use only.

    Scan, process and stream reads from
    FastQ files. Internal use only.

    :: r1file
       type - string
       desc - path to FastQ file storing
              R1 reads
    :: r2file
       type - string / None
       desc - path to FastQ file storing
              R2 reads
    :: r1type
       type - integer
       desc - R1 read type identifier
              0 = Forward Reads
              1 = Reverse Reads
    :: r2type
       type - integer / None
       desc - R2 read type identifier
              0 = Forward Reads
              1 = Reverse Reads
    :: r1length
       type - integer
       desc - minimum required R1 read
              length for packing
    :: r2length
       type - integer / None
       desc - minimum required R2 read
              length for packing
    :: r1qual
       type - integer
       desc - minimum required R1 read
              Q-Score for packing
    :: r2qual
       type - integer / None
       desc - minimum required R2 read
              Q-Score for packing
    :: packdir
       type - string
       desc - path to directory temporarily
              storing packed reads
    :: metaqueue
       type - SafeQueue
       desc - queue storing meta read pack
              file paths once saved into
              packdir
    :: packqueue
       type - SafeQueue
       desc - queue storing read pack file
              paths once saved into packdir
    :: packtype
       type - integer
       desc - packing operation identifier
              0 = concatenated / joined reads
              1 = assembled / merged reads
    :: assemblyparams
       type - dict / None
       desc - dictionary storing read overlap
              merging parameters
    :: packsize
       type - integer
       desc - maximum number of reads (in
              millions) stored per pack
    :: r1truncfile
       type - mp.Event
       desc - multiprocessing Event triggered
              on truncated r1file
    :: r2truncfile
       type - mp.Event
       desc - multiprocessing Event triggered
              on truncated r2file
    :: previousreads
       type - SafeCounter
       desc - total number of read pairs scanned
              previously from r1file and r2file
              by current core
    :: scannedreads
       type - SafeCounter
       desc - total number of read pairs scanned
              from r1file and r2file
    :: ambiguousreads
       type - SafeCounter
       desc - total number of rejected read
              pairs in r1file and r2file due
              to ambiguous bases in the reads
    :: shortreads
       type - SafeCounter
       desc - total number of rejected read
              pairs in r1file and r2file due
              to shorter than expected lengths
    :: survivedreads
       type - SafeCounter
       desc - total number of filtered read
              pairs in r1file and r2file
    :: packedreads
       type - SafeCounter
       desc - total number of uniquely
              packed reads in batches
    :: packsbuilt
       type - SafeCounter
       desc - total number of read packs
              built from r1file and r2file
    :: coreid
       type - integer
       desc - current core integer id
    :: batchid
       type - integer
       desc - current batch integer id
    :: ncores
       type - integer
       desc - total number of packers
              concurrently initiated
    :: nactive
       type - SafeCounter
       desc - total number of packers
              concurrently active
    :: memlimit
       type - nu.Real
       desc - total amount of memory
              allowed per core
    :: launchtime
       type - time
       desc - initial launch timestamp
    :: restart
       type - mp.Event
       desc - multiprocessing event when
              process needs to restart
    :: shutdown
       type - mp.Event
       desc - multiprocessing event when
              process needs to shutdown
    :: liner
       type - coroutine
       desc - dynamic printing
    '''

    # Book-keeping Variables
    t0 = tt.time()
    truncated = False

    c_packsbuilt  = 0
    c_packedreads = 0

    min_dump_reach  = 0
    min_dump_target = packsize // 2

    max_dump_reach  = 0
    max_dump_target = packsize

    # Setup Verbage Variables
    if packtype:
        verbagefactor = 0.2
    else:
        verbagefactor = 1
    verbagetarget = rn.randint(
        *map(round, (min_dump_target * 0.080 * verbagefactor,
                     min_dump_target * 0.120 * verbagefactor)))

    clen = ut.get_printlen(value=ncores)
    plen = ut.get_printlen(value=packsize)

    # Current and Meta Pack Storage
    cpack = cx.Counter()
    mpack = {}
    mpath = None

    # Read Packing Loop
    for reads in cp.stream_processed_fastq(
        r1file=r1file,
        r2file=r2file,
        r1type=r1type,
        r2type=r2type,
        r1length=r1length,
        r2length=r2length,
        r1qual=r1qual,
        r2qual=r2qual,
        packtype=packtype,
        assemblyparams=assemblyparams,
        r1truncfile=r1truncfile,
        r2truncfile=r2truncfile,
        previousreads=previousreads,
        scannedreads=scannedreads,
        ambiguousreads=ambiguousreads,
        shortreads=shortreads,
        survivedreads=survivedreads,
        verbagetarget=verbagetarget,
        coreid=coreid,
        ncores=ncores,
        memlimit=memlimit,
        launchtime=launchtime,
        restart=restart,
        shutdown=shutdown,
        liner=liner):

        # Truncated Files?
        if  r1truncfile.is_set() or \
            r2truncfile.is_set():
            truncated = True
            shutdown.set()
            break

        # Pack Read
        if reads in mpack:
            mpack[reads] += 1
        else:
            cpack[reads] += 1

        # Update Book-keeping
        min_dump_reach  = len(cpack)
        max_dump_reach += 1

        # Time to dump a pack?
        # Dump packs per packsize / 2 (e.g. 1.5M) unique
        # entries or packsize (3.0M) total reads processed,
        # whichever is reached earlier
        if  max_dump_reach >= max_dump_target or \
            min_dump_reach >= min_dump_target:

            # Show Updates
            liner.send(
                ' Core {:{},d}: Built Pack {}.{}.{} w/ {:{},d} Reads in {:.2f} sec\n'.format(
                    coreid,
                    clen,
                    coreid,
                    c_packsbuilt,
                    batchid,
                    len(cpack),
                    plen,
                    tt.time()-t0))

            # Define Pack Path
            cpath = '{}/{}.{}.{}.pack'.format(
                packdir,
                coreid,
                c_packsbuilt,
                batchid)

            # Update Book-keeping
            t0 = tt.time()
            c_packedreads  += len(cpack)
            c_packsbuilt   += 1
            max_dump_reach  = 0
            min_dump_reach  = 0

            # Dump Current Pack
            if mpack:

                # Persist Current Pack to Disk
                ut.savepack(
                    pobj=cpack,
                    filepath=cpath)

                # Enqueue Current Pack Path
                packqueue.put(cpath)

                # Release Control
                tt.sleep(0)

                # Cleanup Current Pack
                cpack.clear()
                del cpack
                cpack = cx.Counter()
                ut.free_mem()

            # Setup Meta Pack
            else:
                mpack = cpack
                mpath = cpath
                cpack = cx.Counter()
                ut.free_mem()

    # Final Dumping
    if not truncated:

        # Dump Final Pack?
        if cpack:

            # Show Updates
            liner.send(
                ' Core {:{},d}: Built Pack {}.{}.{} w/ {:{},d} Reads in {:05.2f} sec\n'.format(
                    coreid,
                    clen,
                    coreid,
                    c_packsbuilt,
                    batchid,
                    len(cpack),
                    plen,
                    tt.time()-t0))

            # Dump Read Pack Appropriately

            # Meta pack
            if c_packsbuilt == 0:

                # Define Pack Path
                cpath = '{}/{}.{}.{}.pack.meta'.format(
                    packdir,
                    coreid,
                    c_packsbuilt,
                    batchid,)

                # Dump Meta Pack
                ut.savemeta(
                    pobj=cpack,
                    filepath=cpath)

                # Enqueue Meta Pack Path
                metaqueue.put(cpath)

                # Release Control
                tt.sleep(0)

            # Non-Meta Pack
            else:

                # Define Pack Path
                cpath = '{}/{}.{}.{}.pack'.format(
                    packdir,
                    coreid,
                    c_packsbuilt,
                    batchid)

                # Dump Pack
                ut.savepack(
                    pobj=cpack,
                    filepath=cpath)

                # Enqueue Pack Path
                packqueue.put(cpath)

                # Release Control
                tt.sleep(0)

            # Update Book-keeping
            c_packedreads += len(cpack)
            c_packsbuilt  += 1

            # Cleanup Final Pack
            cpack.clear()
            del cpack
            ut.free_mem()

        # Dump Meta Pack?
        if mpack:

            # Define Pack Path
            mpath = '{}.meta'.format(
                mpath)

            # Dump Meta Pack
            ut.savemeta(
                pobj=mpack,
                filepath=mpath)

            # Enqueue Meta Pack Path
            metaqueue.put(mpath)

            # Release Control
            tt.sleep(0)

            # Cleanup Meta Pack
            mpack.clear()
            del mpack
            ut.free_mem()

        # Update Read Packing Book-keeping
        packedreads.increment(incr=c_packedreads)
        packsbuilt.increment(incr=c_packsbuilt)

    # Packing Completed
    nactive.decrement()
    if shutdown.is_set():
        packqueue.put(None)

    # Release Control
    tt.sleep(0)

def pack(
    r1file,
    r1type,
    packtype,
    packfile,
    r1length=1,
    r1qual=20,
    r2file=None,
    r2type=None,
    r2length=None,
    r2qual=None,
    packsize=3.0,
    ncores=0,
    memlimit=0,
    verbose=True):
    '''
    TBD
    '''

    # Start Liner
    liner = ut.liner_engine(online=verbose)

    # Packing Verbage Print
    liner.send('\n[Oligopool Calculator: Analysis Mode - Pack]\n')

    # Required Argument Parsing
    liner.send('\n Required Arguments\n')

    # Full r1file Validation
    r1valid = vp.get_readfile_validity(
        readfile=r1file,
        readfile_field='   R1 File   ',
        paired_readfile=None,
        liner=liner)

    # Full r1type Validation
    t1valid = vp.get_categorical_validity(
        category=r1type,
        category_field='   R1 Type   ',
        category_pre_desc=' R1 has ',
        category_post_desc='',
        category_dict={
            0: 'Forward Reads 5\' ---F--> 3\'',
            1: 'Reverse Reads 3\' <--R--- 5\''},
        liner=liner)

    # Full packtype Validation
    packtype_valid = vp.get_categorical_validity(
        category=packtype,
        category_field=' Pack Type   ',
        category_pre_desc=' ',
        category_post_desc='',
        category_dict={
            0: 'Store Concatenated / Joined Reads',
            1: 'Store Assembled / Merged Reads'},
        liner=liner)

    # Full packfile Validation
    packfile_valid = vp.get_outfile_validity(
        outfile=packfile,
        outfile_suffix='.oligopool.pack',
        outfile_field=' Pack File   ',
        liner=liner)

    # Adjust packfile Suffix
    packfile = ut.get_adjusted_path(
        path=packfile,
        suffix='.oligopool.pack')

    # Optional Argument Parsing
    liner.send('\n Optional Arguments\n')

    # Full r1length Validation
    l1valid = vp.get_numeric_validity(
        numeric=r1length,
        numeric_field='   R1 Length ',
        numeric_pre_desc=' Use Reads of Length ',
        numeric_post_desc=' bp or Longer',
        minval=1,
        maxval=float('inf'),
        precheck=False,
        liner=liner)

    # Full r1qual Validation
    q1valid = vp.get_numeric_validity(
        numeric=r1qual,
        numeric_field='   R1 Quality',
        numeric_pre_desc=' Use Reads w/ Mean Q-Score of ',
        numeric_post_desc=' or Higher',
        minval=0,
        maxval=float('inf'),
        precheck=False,
        liner=liner)

    # Full r2file Validation
    r2valid = vp.get_optional_readfile_validity(
        readfile=r2file,
        readfile_field='   R2 File   ',
        paired_readfile=r1file,
        liner=liner)

    # Full r2type Validation
    if r2valid and (not r2file is None):
        validfn = vp.get_categorical_validity
    else:
        validfn = vp.get_optional_categorical_validity
    t2valid = validfn(
        category=r2type,
        category_field='   R2 Type   ',
        category_pre_desc=' R2 has ',
        category_post_desc='',
        category_dict={
            0: 'Forward Reads 5\' ---F--> 3\'',
            1: 'Reverse Reads 3\' <--R--- 5\''},
        liner=liner)

    # Full r2length Validation
    if r2valid and (not r2file is None):
        validfn = vp.get_numeric_validity
    else:
        validfn = vp.get_optional_numeric_validity
    l2valid = validfn(
        numeric=r2length,
        numeric_field='   R2 Length ',
        numeric_pre_desc=' Use Reads of Length ',
        numeric_post_desc=' bp or Longer',
        minval=1,
        maxval=float('inf'),
        precheck=False,
        liner=liner)

    # Full r2qual Validation
    if r2valid and (not r2file is None):
        validfn = vp.get_numeric_validity
    else:
        validfn = vp.get_optional_numeric_validity
    q2valid = validfn(
        numeric=r2qual,
        numeric_field='   R2 Quality',
        numeric_pre_desc=' Use Reads w/ Mean Q-Score of ',
        numeric_post_desc=' or Higher',
        minval=0,
        maxval=float('inf'),
        precheck=False,
        liner=liner)

    # Full packsize Validation
    packsize_valid = vp.get_numeric_validity(
        numeric=packsize,
        numeric_field=' Pack Size   ',
        numeric_pre_desc=' Store up to ',
        numeric_post_desc=' Million Reads per Pack',
        minval=0.10,
        maxval=5.00,
        precheck=False,
        liner=liner)

    # Full num_core Parsing and Validation
    (ncores,
    ncores_valid) = vp.get_parsed_core_info(
        ncores=ncores,
        core_field='  Num Cores  ',
        default=None if not packtype else mp.cpu_count() // 3,
        offset=2,
        liner=liner)

    # Full num_core Parsing and Validation
    (memlimit,
    memlimit_valid) = vp.get_parsed_memory_info(
        memlimit=memlimit,
        memlimit_field='  Mem Limit  ',
        ncores=ncores,
        ncores_valid=ncores_valid,
        liner=liner)

    # First Pass Validation
    if not all([
        r1valid,
        t1valid,
        packfile_valid,
        l1valid,
        q1valid,
        r2valid,
        t2valid,
        l2valid,
        q2valid,
        packtype_valid,
        packsize_valid,
        ncores_valid,
        memlimit_valid]):
        liner.send('\n')
        raise RuntimeError(
            'Invalid Argument Input(s).')

    # Start Timer
    t0 = tt.time()

    # Adjust Numeric Parameters
    if r2file is None:
        packtype = 0

    # Setup Warning Dictionary
    warns = {}

    # Some Assembly Required?
    if packtype:

        # Compute Assembly Parameters
        liner.send('\n[Step 1: Extracting Assembly Parameters]\n')

        # Extract Overlap Parameters
        (mergefnhigh,
        mergefnlow) = cp.get_extracted_overlap_parameters(
            r1file=r1file,
            r2file=r2file,
            r1type=r1type,
            r2type=r2type,
            liner=liner)

        # Store Parameters
        assemblyparams = {
            'mergefnhigh': mergefnhigh,
             'mergefnlow': mergefnlow}

    # Vanilla Storage
    else:
        # Nothing to see here ...
        assemblyparams = None

    # Free Memory
    ut.free_mem()

    # Setup Workspace
    (packfile,
    packdir) = ut.setup_workspace(
        outfile=packfile,
        outfile_suffix='.oligopool.pack')

    # Schedule packfile deletion
    pkdeletion = ae.register(
        ut.remove_file,
        packfile)

    # Expand packsize
    packsize = int(packsize * (10.**6))

    # Read Pack Queues
    metaqueue = ut.SafeQueue()
    packqueue = ut.SafeQueue()

    # Read File Processing Book-keeping
    r1truncfile  = mp.Event()
    r2truncfile  = mp.Event()
    restarts  = [
        mp.Event() for _ in range(ncores)]
    shutdowns = [
        mp.Event() for _ in range(ncores)]

    # Read Packing Book-keeping
    nactive        = ut.SafeCounter(initval=ncores)
    scannedreads   = ut.SafeCounter()
    ambiguousreads = ut.SafeCounter()
    shortreads     = ut.SafeCounter()
    survivedreads  = ut.SafeCounter()
    packedreads    = ut.SafeCounter()
    packsbuilt     = ut.SafeCounter()
    batchids       = [0] * ncores
    previousreads  = [
        ut.SafeCounter() for _ in range(ncores)]

    # Launching Read Packing
    liner.send('\n[Step 2: Computing Read Packs]\n')

    # Define Packing Stats
    stats = {
        'status'  : False,
        'basis'   : 'unsolved',
        'step'    : 2,
        'stepname': 'computing-read-packs',
        'vars'    : {
            'metaaggcount': 0,
             'r1truncated': False,
             'r2truncated': False,
                'packsize': packsize,
               'packcount': int(packsbuilt.value()),
                 'scanned': int(scannedreads.value()),
               'ambiguous': int(ambiguousreads.value()),
                   'short': int(shortreads.value()),
                'survived': int(survivedreads.value()),
                  'packed': int(packedreads.value())},
        'warns'   : warns}

    # Engine Timer
    et = tt.time()

    # Define Archiver (Non-Meta Read Packs)
    archiver = mp.Process(
        target=ut.archive,
        args=(packqueue,
            packfile,
            'x',
            ncores,
            nactive,
            liner,))

    # Fire-off Archiver
    archiver.start()

    # Define Packer Process Store
    readpackers = []

    # Fire-off Initial Read Packers
    coreid = 0
    clen = ut.get_printlen(value=ncores)
    while coreid < ncores:

        # Define Packer
        readpacker = mp.Process(
            target=pack_engine,
            args=(r1file,
                r2file,
                r1type,
                r2type,
                r1length,
                r2length,
                r1qual,
                r2qual,
                packdir,
                metaqueue,
                packqueue,
                packtype,
                assemblyparams,
                packsize,
                r1truncfile,
                r2truncfile,
                previousreads[coreid],
                scannedreads,
                ambiguousreads,
                shortreads,
                survivedreads,
                packedreads,
                packsbuilt,
                coreid,
                batchids[coreid],
                ncores,
                nactive,
                memlimit,
                et,
                restarts[coreid],
                shutdowns[coreid],
                liner,))

        # Show Update
        liner.send(
            ' Core {:{},d}: Starting Up\n'.format(
                coreid,
                clen))

        # Start Packer
        readpacker.start()

        # Update Book-keeping
        readpackers.append(readpacker)
        coreid += 1

    # Packer Management
    coreid = 0
    activepackers = ncores
    while activepackers:

        # Had Packer Finished?
        if readpackers[coreid] is None:
            pass

        # Has Packer Shutdown?
        elif shutdowns[coreid].is_set():
            # Cleanup
            readpackers[coreid].join()
            readpackers[coreid].close()
            # Update
            readpackers[coreid] = None
            activepackers -= 1
            ut.free_mem()
            # Reset
            restarts[coreid].clear()
            shutdowns[coreid].clear()

        # Must Packer Restart?
        elif restarts[coreid].is_set():
            # Cleanup
            readpackers[coreid].join()
            readpackers[coreid].close()
            # Update
            readpackers[coreid] = None
            batchids[coreid] += 1
            ut.free_mem()
            # Reset
            restarts[coreid].clear()
            shutdowns[coreid].clear()
            readpackers[coreid] = mp.Process(
                target=pack_engine,
                args=(r1file,
                    r2file,
                    r1type,
                    r2type,
                    r1length,
                    r2length,
                    r1qual,
                    r2qual,
                    packdir,
                    metaqueue,
                    packqueue,
                    packtype,
                    assemblyparams,
                    packsize,
                    r1truncfile,
                    r2truncfile,
                    previousreads[coreid],
                    scannedreads,
                    ambiguousreads,
                    shortreads,
                    survivedreads,
                    packedreads,
                    packsbuilt,
                    coreid,
                    batchids[coreid],
                    ncores,
                    nactive,
                    memlimit,
                    et,
                    restarts[coreid],
                    shutdowns[coreid],
                    liner,))
            readpackers[coreid].start()

        # Next Iteration
        coreid = (coreid + 1) % ncores
        tt.sleep(0)

    # Join Archiver
    archiver.join()
    archiver.close()

    # Free Memory
    ut.free_mem()

    # Handle Truncated Read Files
    if r1truncfile.is_set():
        liner.send(
            ' R1 File Truncated or Incompatible with R2 File\n')
        ut.remove_file(
            filepath=packfile)
        stats['vars']['r1truncated'] = True
    elif r2truncfile.is_set():
        liner.send(
            ' R2 File Truncated or Incompatible with R1 File\n')
        ut.remove_file(
            filepath=packfile)
        stats['vars']['r1truncated'] = True

    # Packing Successful
    else:
        stats['status'] = True
        stats['basis']  = 'solved'

    # Show Time Elapsed
    liner.send(
        ' Time Elapsed: {:.2f} sec\n'.format(
            tt.time() - et))

    # Did we succeed?
    if stats['status']:

        # Launching Read Packing
        liner.send('\n[Step 3: Aggregating Meta Packs]\n')

        # Update Stats
        stats['step'] = 3
        stats['stepname'] = 'aggregating-meta-packs'
        stats['vars']['metaaggcount'] = len(metaqueue)

        # Aggregation Timer
        at = tt.time()

        # Aggregate Meta Read Packs
        # Define Meta Pack Aggregator
        # on a Separate Core
        aggregator = mp.Process(
            target=cp.pack_aggregator,
            args=(metaqueue,
                packqueue,
                packedreads,
                packsbuilt,
                liner,))

        # Start Aggregator
        aggregator.start()

        # Update Producer Count
        nactive.increment()

        # Archive Meta Read Packs
        ut.archive(
            objqueue=packqueue,
            arcfile=packfile,
            mode='a',
            prodcount=1,
            prodactive=nactive,
            liner=liner)

        # Update Producer Count
        nactive.decrement()

        # Join and Close Aggregator
        aggregator.join()
        aggregator.close()

        # Show Time Elapsed
        liner.send(
            ' Time Elapsed: {:.2f} sec\n'.format(
                tt.time() - at))

    # Update Stats
    stats['vars']['packcount'] = int(packsbuilt.value())
    stats['vars']['scanned']   = int(scannedreads.value())
    stats['vars']['ambiguous'] = int(ambiguousreads.value())
    stats['vars']['short']     = int(shortreads.value())
    stats['vars']['survived']  = int(survivedreads.value())
    stats['vars']['packed']    = int(packedreads.value())

    # Packing Status
    if stats['status']:
        packstatus = 'Successful'
    else:
        packstatus = 'Failed'

    # Read Packing Stats
    liner.send('\n[Packing Statistics]\n')

    plen = ut.get_printlen(
        value=scannedreads.value())

    liner.send(
        '   Packing Status: {}\n'.format(
            packstatus))
    liner.send(
        '     R1 Truncated: {}\n'.format(
            ['No', 'Yes'][r1truncfile.is_set()]))
    liner.send(
        '     R2 Truncated: {}\n'.format(
            ['No', 'Yes'][r2truncfile.is_set()]))
    liner.send(
        '   Scanned Reads : {:{},d}\n'.format(
            scannedreads.value(),
            plen))
    liner.send(
        ' Ambiguous Reads : {:{},d} ({:6.2f} %)\n'.format(
            ambiguousreads.value(),
            plen,
            ut.safediv(
                A=(100. * ambiguousreads.value()),
                B=scannedreads.value())))
    liner.send(
        '     Short Reads : {:{},d} ({:6.2f} %)\n'.format(
            shortreads.value(),
            plen,
            ut.safediv(
                A=(100. * shortreads.value()),
                B=scannedreads.value())))
    liner.send(
        '  Survived Reads : {:{},d} ({:6.2f} %)\n'.format(
             survivedreads.value(),
            plen,
            ut.safediv(
                A=(100. * survivedreads.value()),
                B=scannedreads.value())))
    liner.send(
        '    Packed Reads : {:{},d}\n'.format(
            packedreads.value(),
            plen))
    liner.send(
        '   Packing Ratio : {:.2f} to 1\n'.format(
            ut.safediv(
                A=(1. * survivedreads.value()),
                B=packedreads.value())))
    liner.send(
        '   Packing Order : {:.2f} %\n'.format(
            (100. * (
                1. - ut.safediv(
                    A=packedreads.value(),
                    B=survivedreads.value())))))
    liner.send(
        '     Packs Built : {} Read Packs\n'.format(
            packsbuilt.value()))
    liner.send(
        '  Time Elapsed: {:.2f} sec\n'.format(
            tt.time() - t0))

    # Did we succeed?
    if stats['status']:

        # Archive Packing Stats
        packstat = {
            'packtype' : packtype,
            'packsize' : packsize,
            'packcount': int(packsbuilt.value()),
            'scanned'  : int(scannedreads.value()),
            'survived' : int(survivedreads.value()),
            'packed'   : int(packedreads.value())}

        # Prepare Archiving Sentinels
        statpath = '{}/packing.stat'.format(
            packdir)
        packqueue.put(statpath)
        packqueue.put(None)

        # Dump Pack Stat
        ut.savedict(
            dobj=packstat,
            filepath=statpath)

        # Update Producer Count
        nactive.increment()

        # Archive Pack Stat
        ut.archive(
            objqueue=packqueue,
            arcfile=packfile,
            mode='a',
            prodcount=1,
            prodactive=nactive,
            liner=liner)

        # Update Producer Count
        nactive.decrement()

        # Close Queues
        metaqueue.close()
        packqueue.close()

    # Remove Workspace
    ut.remove_directory(
        dirpath=packdir)

    # Unschedule packfile deletion
    if packstatus == 'Successful':
        ae.unregister(pkdeletion)

    # Close Liner
    liner.close()

    # Return Statistics
    return stats
