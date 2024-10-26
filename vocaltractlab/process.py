import multiprocessing
import tqdm


def _worker( args ):
    function, arg = args
    return function( **arg )

def OLD_multiprocess(
    function,
    args,
    initializer,
    initargs,
    return_data: bool = False,
    verbose = True,
    workers: int = None,
    ):
    if workers == None:
        workers = multiprocessing.cpu_count()
    pool = multiprocessing.Pool( workers )
    tasks = ( ( function, x ) for x in args)
    data = None
    if return_data:
        data = []
        if verbose:
            for x in tqdm.tqdm( pool.imap( _worker, tasks ), total=len( args ) ):
                data.append( x )
        else:
            for x in pool.imap( _worker, tasks ):
                data.append( x )
    else:
        if verbose:
            for x in tqdm.tqdm( pool.imap( _worker, tasks ), total=len( args ) ):
                pass
        else:
            pool.imap( _worker, tasks )
    pool.close()
    pool.join()
    return data

def multiprocess(
    function,
    args,
    initializer,
    initargs,
    return_data: bool = False,
    verbose = True,
    workers: int = None,
    ):
    if workers == None:
        workers = multiprocessing.cpu_count()

    # Set up the pool and initialize each worker with the config path
    with multiprocessing.Pool(
        workers,
        initializer=initializer,
        initargs=initargs,
    ) as pool:
        tasks = ((function, x) for x in args)
        if return_data:
            data = []
            if verbose :
                for x in tqdm.tqdm(pool.imap(_worker, tasks), total=len(args)):
                    data.append(x)
            else:
                for x in pool.imap(_worker, tasks):
                    data.append(x)
        else:
            data = None
            if verbose:
                for x in tqdm.tqdm(pool.imap(_worker, tasks), total=len(args)):
                    pass
            else:
                pool.imap(_worker, tasks)
    return data

def process(
        function,
        args,
        initializer,
        initargs,
        return_data: bool = False,
        verbose = True,
        workers: int = None,
        mp_threshold: int = 4,
        ):
    if len( args ) < mp_threshold:
        return single_process(
            function = function,
            args = args,
            return_data = return_data,
            verbose = verbose,
            )
    else:
        return multiprocess(
            function = function,
            args = args,
            initializer = initializer,
            initargs = initargs,
            return_data = return_data,
            verbose = verbose,
            workers = workers,
            )

def single_process(
        function,
        args,
        return_data: bool = False,
        verbose = True,
        ):
    y = tqdm.tqdm( args ) if verbose else args
    if return_data:
        data = [ function( **x ) for x in y ]
    else:
        [ function( **x ) for x in y ]
        data = None
    return data