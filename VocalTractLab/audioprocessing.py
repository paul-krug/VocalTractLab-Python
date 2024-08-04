


import os
import torch
import torchaudio
import torchaudio.functional as F
#import torchaudio.transforms as T
import numpy as np
#import librosa
#from scipy import interpolate as ip
#import parselmouth
#import matplotlib.pyplot as plt
from vocaltractlab_cython import get_constants

from typing import Tuple
from typing import Union
from typing import Optional
from typing import Dict
from numpy.typing import ArrayLike



torch.set_num_threads(1)
torch.multiprocessing.set_sharing_strategy('file_system')

MAX_WAV_VALUE = 32768.0

def to_float(
          x: Union[torch.Tensor, ArrayLike],
    ) -> torch.Tensor:
    """Converts a tensor of ints into floats in the range [-1, 1].
    Args:
        x (Union[torch.Tensor, ArrayLike]): Tensor of ints with arbitrary shape.
    Returns:
        torch.Tensor: Tensor of floats with same shape as x.
    """
    # Convert to torch.Tensor if needed.
    if not isinstance(x, torch.Tensor):
        x = torch.tensor(x)
    # Convert to float if needed.
    if not x.dtype == torch.float:
        input_dtype = x.dtype
        x = x.float()
        if input_dtype == torch.int:
            x /= MAX_WAV_VALUE
    return x

def to_int(
        waveform,
    ):
    """
    Convert any audio array to Torch int tensor.
    :param waveform: Audio to convert.
    :return: Audio as int16.
    """
    # Convert to torch tensor
    if not isinstance(waveform, torch.Tensor):
        waveform = torch.tensor(waveform)
    # Conver to int   
    if not waveform.dtype == torch.int16:
        waveform = waveform * 32768
        waveform = waveform.short()
    return waveform

def normalize_audio_amplitude(
        x,
        dBFS = -1,
        ): #normalisation in dB
	norm_factor = 10**( -1 * dBFS * 0.05 ) -1
	norm_max = torch.max( torch.abs( x ) )#, axis=0)
	x /= ( norm_max + ( norm_max * norm_factor ) )
	return x

def postprocess(
        x: ArrayLike,
        sr_out: int,
        dBFS: int = -1,
        file_path: str = None,
        to_numpy: bool = False,
        ) -> np.ndarray:
    
    vtl_constants = get_constants()

    x = torch.tensor( x ).unsqueeze( 0 )

    if sr_out is None:
        sr_out = vtl_constants[ 'sr_audio' ]
    elif sr_out != vtl_constants[ 'sr_audio' ]:
        x = resample_like_librosa(
            x = x,
            sr_in = vtl_constants[ 'sr_audio' ],
            sr_out = sr_out,
            )
    
    if dBFS is not None:
        x = normalize_audio_amplitude(
            x = x,
            dBFS = dBFS,
            )
        
    if file_path is not None:
        if not os.path.exists(
            os.path.dirname( file_path )
            ):
            os.makedirs(
                os.path.dirname( file_path ),
                exist_ok = True,
                )
        torchaudio.save(
            file_path,
            x,
            sr_out,
            )
        
    if to_numpy:
        x = x.numpy()
    
    return x

def resample_like_librosa(
        x: Union[torch.Tensor, ArrayLike],
        sr_in: int,
        sr_out: int,
    ) -> torch.Tensor:
    """
    Resample a time series, similar to librosa
    with 'kaiser_best' resampling method.
    Args:
        x (Union[torch.Tensor, ArrayLike]): Tensor of ints with arbitrary shape.
        sr_in (int): Input sampling rate.
        sr_out (int): Output sampling rate.
    Returns:
        torch.Tensor: Tensor of floats with same shape as x.
    """
    # Convert to torch.Tensor if needed.
    x = to_float(x)
    if sr_in != sr_out:
        x = F.resample(
            waveform = x,
            orig_freq = sr_in,
            new_freq = sr_out,
            lowpass_filter_width = 64,
            rolloff = 0.9475937167399596,
            resampling_method = 'sinc_interp_kaiser',
            beta = 14.769656459379492,
            )
    return x

def hz_to_st(
        frequency_hz,
        reference = 1.0,
    ):
    return 12.0*np.log( frequency_hz / reference ) / np.log(2.0)

def st_to_hz(
        frequency_st,
        reference = 1.0,
    ):
    return reference*pow( 2, frequency_st / 12.0 )

def power_to_db(
    x: np.ndarray,  
    ref: float = 1.0,
    eps: float = 1e-10,
    ) -> np.ndarray:

    x = np.asarray(x)

    if eps <= 0:
        raise ValueError(
            "Arg eps must be positive"
            )

    ref_value = np.abs(ref)

    db = 10.0 * np.log10(np.maximum(eps, x))
    db -= 10.0 * np.log10(np.maximum(eps, ref_value))

    return db

def amplitude_to_db(
    x: np.ndarray,
    **kwargs,
    ) -> np.ndarray:
    return power_to_db(
        np.abs(x) ** 2,
        **kwargs,
        )