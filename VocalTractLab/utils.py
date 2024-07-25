


from typing import Union, List, Tuple, Dict, Any, Optional, Callable, Iterable, Sequence
from numpy.typing import ArrayLike



def make_iterable( x ):
    if isinstance( x, str ) or not isinstance( x, Iterable ):
        return [ x ]
    return x