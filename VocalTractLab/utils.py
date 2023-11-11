
'''
Requires:
- numpy
- pandas
- torchaudio
- vocaltractlab_cython
- tools_mp
- whatever
'''


from typing import Union, List, Tuple, Dict, Any, Optional, Callable, Iterable, Sequence
from numpy.typing import ArrayLike

from tools_mp import multiprocess


def df_glottis_params_from_motor_file( index ):
    if (index > 7) and (index % 2 == 0):
        return False
    else:
        return True

def df_tract_params_from_motor_file( index ):
    if (index > 7) and ((index-1) % 2 == 0):
        return False
    else:
        return True

def make_iterable( x ):
    if isinstance( x, str ) or not isinstance( x, Iterable ):
        return [ x ]
    return x





if __name__ == '__main__':
    gesture_to_motor(
        gesture_files='test_1.csv',

        motor_files= 'deine_mudda',
    )