
'''
Requires:
- numpy
- pandas
- torchaudio
- vocaltractlab_cython
- tools_mp
- whatever
'''
import os
import numpy as np
import pandas as pd

from TargetApproximationModel import TargetSeries

#import torchaudio
from vocaltractlab_cython import get_constants
from vocaltractlab_cython import get_param_info


from typing import Union, List, Tuple, Dict, Any, Optional, Callable, Iterable, Sequence
from numpy.typing import ArrayLike

from .utils import df_tract_params_from_motor_file
from .utils import df_glottis_params_from_motor_file




class GestureScore:
    pass
class MotorScore:
    pass
class MotorSeries:
    pass



class TractSeries( TargetSeries ):
    
    def __init__(
            self,
            series: np.ndarray,
            sr: float,
            ):
        tiers = [
            x[ 'name' ]
            for x in get_param_info( 'tract' )
            ]
        super().__init__( series, sr, tiers )
        return
    
    @classmethod
    def from_motor_file(
            cls,
            path: str,
            sr = None,
            ):
        if sr is None:
            sr = get_constants()[ 'sr_internal' ]
        df = pd.read_csv(
            path,
            delim_whitespace = True,
            skiprows= lambda x: df_tract_params_from_motor_file(x),
            header = None,
            )
        return cls( df.to_numpy(), sr = sr )



class GlottalSeries( TargetSeries ):
    
    def __init__(
            self,
            series: np.ndarray,
            sr: float,
            ):
        tiers = [
            x[ 'name' ]
            for x in get_param_info( 'glottis' )
            ]
        super().__init__( series, sr, tiers )
        return
    
    @classmethod
    def from_motor_file(
            cls,
            path: str,
            sr = None,
            ):
        if sr is None:
            sr = get_constants()[ 'sr_internal' ]
        df = pd.read_csv(
            path,
            delim_whitespace = True,
            skiprows= lambda x: df_glottis_params_from_motor_file(x),
            header = None,
            )
        return cls( df.to_numpy(), sr = sr )



class MotorSeries( TargetSeries ):
    
    def __init__(
            self,
            tract_series: TractSeries,
            glottal_series: GlottalSeries,
            ):
        if not isinstance( tract_series, TractSeries ):
            raise TypeError(
                f'Argument "tract_series" must be of type TractSeries, not {type( tract_series )}!'
                )
        if not isinstance( glottal_series, GlottalSeries ):
            raise TypeError(
                f'Argument "glottal_series" must be of type GlottalSeries, not {type( glottal_series )}!'
                )
        if tract_series.sr != glottal_series.sr:
            raise ValueError(
                f'The sampling rate of the tract series: {tract_series.sr} does not match the sampling rate of the glottal series: {glottal_series.sr}!'
                )
        if len( tract_series ) != len( glottal_series ):
            raise ValueError(
                f'The length of the tract series: {len( tract_series )} does not match the length of the glottal series: {len( glottal_series )}!'
                )
        self.tract_tiers = tract_series.tiers
        self.glottal_tiers = glottal_series.tiers
        tiers = self.tract_tiers + self.glottal_tiers
        sr = tract_series.sr
        series = np.concatenate( [ tract_series.series, glottal_series.series ], axis = 1 )
        super().__init__( series, sr, tiers )
        return

    @classmethod
    def from_motor_file(
            cls,
            path: str,
            sr = None,
            ):
        if sr is None:
            sr = get_constants()[ 'sr_internal' ]
        df_vtp = pd.read_csv(
            path,
            delim_whitespace = True,
            skiprows= lambda x: df_tract_params_from_motor_file(x),
            header = None,
            )
        df_glp = pd.read_csv(
            path,
            delim_whitespace = True,
            skiprows= lambda x: df_glottis_params_from_motor_file(x),
            header = None,
            )
        return cls(
            TractSeries( df_vtp.to_numpy(), sr = sr ),
            GlottalSeries( df_glp.to_numpy(), sr = sr ),
            )
    
    def to_numpy(
            self,
            part = 'all',
            ):
        if part == 'all':
            x = super().to_numpy()
        elif part == 'glottis':
            x = self[ self.glottal_tiers ].to_numpy()
        elif part == 'tract':
            x = self[ self.tract_tiers ].to_numpy()
        else:
            raise ValueError(
                f'Argument "part" must be one of the following: "all", "glottis", "tract", not {part}!'
                )
        return x

    def to_glottal_series( self ):
        return GlottalSeries(
            self[ self.glottal_tiers ].to_numpy(),
            self.sr,
            )

    def to_tract_series( self ):
        return TractSeries(
            self[ self.tract_tiers ].to_numpy(),
            self.sr,
            )