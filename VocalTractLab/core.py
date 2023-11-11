
'''
Requires:
- numpy
- pandas
- torch
- torchaudio
- soundfile
- vocaltractlab_cython
- tools_mp
- whatever
'''
import os
import numpy as np
import pandas as pd
import torch
import torchaudio
from vocaltractlab_cython import get_constants
from vocaltractlab_cython import gesture_file_to_audio
from vocaltractlab_cython.VocalTractLabApi import _synth_block

from typing import Union, List, Tuple, Dict, Any, Optional, Callable, Iterable, Sequence
from numpy.typing import ArrayLike

from tools_mp import multiprocess

from .utils import make_iterable
from .dynamics import GestureScore
from .dynamics import MotorScore
from .dynamics import MotorSeries
from .audioprocessing import normalize_audio_amplitude
from .audioprocessing import resample_like_librosa




#def make_temporary_file():
#    

#gestural_score_to_audio(
#        ges_file_path: str,
#        audio_file_path: str = None,

def gesture_to_audio(
        gesture_data: Union[ Iterable[ str ], str ],
        audio_files: Optional[ Union[ Iterable[ str ], str ] ],
        normalize_audio: int = -1,
        sr: int = None,
        return_data: bool = False,
        workers: int = None,
        verbose: bool = True,
        ) -> None:

    gesture_data = make_iterable( gesture_data )
    audio_files = make_iterable( audio_files )
    if len( gesture_data ) != len( audio_files ):
        raise ValueError(
            f"""
            The number of gesture file paths: {len(gesture_data)}
            does not match the number of audio file paths: {len(audio_files)}.
            """
            )
    
    args = [
        dict(
            gesture_data = gd,
            audio_file_path = audio_file_path,
            verbose_api = False,
            )
        for gd, audio_file_path in zip(
            gesture_data,
            audio_files,
            )
        ]
    audio_data = multiprocess(
        _gesture_to_audio,
        args = args,
        return_data = return_data,
        workers = workers,
        verbose = verbose,
        mp_threshold = 4,
        )
    return

def _gesture_to_audio(
        gesture_data,
        audio_file_path,
        verbose_api,
        normalize_audio: int = -1,
        sr: int = None,
        ) -> np.ndarray:
    if isinstance( gesture_data, GestureScore ):
        gesture_file = gesture_data.to_gesture_file( file_path = None )
    elif isinstance( gesture_data, str ):
        gesture_file = gesture_data
    else:
        raise TypeError(
            f"""
            The specified motor data type: '{type(motor_data)}'
            is not supported. Type must be one of the following:
            - str
            - MotorScore
            - MotorSeries
            """
            )
    audio = gestural_score_to_audio(
        ges_file_path = gesture_file,
        audio_file_path = None,
        verbose_api = verbose_api,
    )
    
    if sr is not None:
        vtl_constants = get_constants()
        audio = resample( audio, sr, vtl_constants[ 'sampling_rate' ] )
    
    if normalize_audio is not None:
        audio = normalize( audio, normalize_audio )

    if audio_file_path is not None:
        torchaudio.save( audio_file_path, audio, sr )
    
    return audio
    

def gesture_to_motor(
        gesture_files: Union[ Iterable[ str ], str ],
        motor_files: Optional[ Union[ Iterable[ str ], str ] ],
        return_data: bool = False,
        workers: int = None,
        verbose: bool = True,
        ) -> None:

    gesture_files = make_iterable( gesture_files )
    motor_files = make_iterable( motor_files )
    if len( gesture_files ) != len( motor_files ):
        raise ValueError(
            f"""
            The number of gesture file paths: {len(gesture_files)}
            does not match the number of motor file paths: {len(motor_files)}.
            """
            )
    
    args = [
        dict(
            motor_data = md,
            audio_file_path = audio_file_path,
            normalize_audio = normalize_audio,
            sr = sr,
            )
        for md, audio_file_path in zip(
            motor_data,
            audio_file_path_list,
            )
        ]
    audio_data = multiprocess(
        _motor_to_audio,
        args = args,
        return_data = return_data,
        workers = workers,
        verbose = verbose,
        mp_threshold = 4,
        )
    return

def motor_to_audio(
        motor_data: Union[ MotorScore, MotorSeries, str ],
        audio_files: Optional[ Union[ Iterable[str], str ] ] = None,
        normalize_audio: int = -1,
        sr: int = None,
        return_data: bool = False,
        workers: int = None,
        verbose: bool = True,
        ) -> np.ndarray:
    motor_data = make_iterable( motor_data )
    if audio_files is None:
        audio_files = [ None ] * len( motor_data )
    else:
        audio_files = make_iterable( audio_files )
    if len( audio_files ) != len( motor_data ):
        raise ValueError(
            f"""
            The number of audio file paths: {len(audio_files)}
            does not match the number of motor data: {len(motor_data)}.
            """
            )
    
    args = [
        dict(
            motor_data = md,
            audio_file_path = audio_file_path,
            normalize_audio = normalize_audio,
            sr = sr,
            )
        for md, audio_file_path in zip(
            motor_data,
            audio_files,
            )
        ]
    audio_data = multiprocess(
        _motor_to_audio,
        args = args,
        return_data = return_data,
        workers = workers,
        verbose = verbose,
        #mp_threshold = 4,
        )
    return audio_data

def _motor_to_audio(
        motor_data,
        audio_file_path,
        normalize_audio,
        sr,
        state_samples = None,
        ):
    if isinstance( motor_data, str ):
        if not os.path.exists( motor_data ):
            raise FileNotFoundError( 
                f"""
                The specified motor file path: '{motor_data}'
                does not exist.
                """
            )
        motor_series = MotorSeries.from_motor_file( motor_data )
    elif isinstance( motor_data, MotorScore ):
        motor_series = motor_data.to_motor_series()
    elif isinstance( motor_data, MotorSeries ):
        motor_series = motor_data
    else:
        raise TypeError(
            f"""
            The specified motor data type: '{type(motor_data)}'
            is not supported. Type must be one of the following:
            - str
            - MotorScore
            - MotorSeries
            """
            )
    vtl_constants = get_constants()
    if state_samples is None:
        state_samples = vtl_constants[ 'n_samples_per_state' ]

    print( motor_series.to_numpy( part='tract' ) )

    tract_params = motor_series.to_numpy( part='tract' )
    glottal_params = motor_series.to_numpy( part='glottis' )
    print( tract_params.shape )
    print( glottal_params.shape )
    print( state_samples )

    
    audio = _synth_block(
        tract_parameters = tract_params,
        glottis_parameters = glottal_params,
        state_samples = state_samples,
        verbose_api = False,
        )
    
    audio = torch.tensor( audio ).unsqueeze( 0 )
    print( 'audio shape: ', audio.shape )
    
    if sr is None:
        sr = vtl_constants[ 'sr_audio' ]
    elif sr != vtl_constants[ 'sr_audio' ]:
        audio = resample_like_librosa( audio, sr, vtl_constants[ 'sr_audio' ] )
    
    if normalize_audio is not None:
        audio = normalize_audio_amplitude( audio, normalize_audio )

    if audio_file_path is not None:
        torchaudio.save( audio_file_path, audio, sr )
    
    return audio

def phoneme_to_gesture(
        phoneme: str,
        ) -> np.ndarray:
    vtl_constants = get_constants()
    return vtl_constants[ 'phonemes' ][ phoneme ]

def _phoneme_to_gesture(
        phoneme: str,
        ) -> np.ndarray:
    vtl_constants = get_constants()
    return vtl_constants[ 'phonemes' ][ phoneme ]

def phoneme_to_motor():
    pass

def _phoneme_to_motor():
    pass


#def _supra_glottal_state_to_svg_str( args ):
#    supra_glottal_state = args
#    svgStr = ( ' ' * 10000 ).encode()
#    constants = get_constants()
#    cdef np.ndarray[np.float64_t, ndim = 1] tractParams = np.zeros(
#        constants['n_tract_params'],
#        dtype = 'float64',
#        )
#    tractParams = supra_glottal_state.ravel()
#    vtlExportTractSvgToStr(
#        &tractParams[0],
#        svgStr,
#        )
#    return svgStr.decode()



if __name__ == '__main__':
    gesture_to_motor(
        gesture_files='test_1.csv',

        motor_files= 'deine_mudda',
    )