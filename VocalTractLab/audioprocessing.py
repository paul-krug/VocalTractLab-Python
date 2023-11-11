import os
import torch
import torchaudio
import torchaudio.functional as F
import torchaudio.transforms as T
import numpy as np
#import librosa
#from scipy import interpolate as ip
#import parselmouth
import matplotlib.pyplot as plt

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



def dynamic_range_compression(x, C=1, clip_val=1e-5):
    return np.log(np.clip(x, a_min=clip_val, a_max=None) * C)



def dynamic_range_decompression(x, C=1):
    return np.exp(x) / C



def dynamic_range_compression_torch(x, C=1, clip_val=1e-5):
    return torch.log(torch.clamp(x, min=clip_val) * C)



def dynamic_range_decompression_torch(x, C=1):
    return torch.exp(x) / C



def spectral_normalize_torch(magnitudes):
    output = dynamic_range_compression_torch(magnitudes)
    return output



def spectral_de_normalize_torch(magnitudes):
    output = dynamic_range_decompression_torch(magnitudes)
    return output



def normalize_audio_amplitude( waveform, normalization = -1 ): #normalisation in dB
	norm_factor = 10**( -1 * normalization * 0.05 ) -1
	norm_max = torch.max( torch.abs( waveform ) )#, axis=0)
	waveform /= ( norm_max + ( norm_max * norm_factor ) )
	return waveform



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



def remove_outlier_data(
        data,
    ):
    # finding the 1st quartile
    #try:
    q1 = np.quantile(data[:,1], 0.25)
    #except Exception:
    #	print( data[:,1] )
    #	#plt.show()
    #	stop
    # finding the 3rd quartile
    q3 = np.quantile(data[:,1], 0.75)
    med = np.median(data[:,1])
    
    # finding the iqr region
    iqr = q3-q1
    
    # finding upper and lower whiskers
    upper_bound = q3+(1.5*iqr)
    lower_bound = q1-(1.5*iqr)
    #upper_bound = np.quantile(data[:,1], 0.99)
    #lower_bound = np.quantile(data[:,1], 0.01)
    data = data[ (data[:,1] > lower_bound) & (data[:,1] < upper_bound) ]
    return data



class Wav2Mel(torch.nn.Module):
    def __init__(
            self,
            sr_out: int = 22000,
            n_fft: int = 1024,
            win_length: int = 1024,
            hop_length: int = 220,
            f_min: int = 0,
            f_max: int = 8000,
            n_mels: int = 80,
            normalize_amplitude: bool = True,
            transpose: bool = True,
            ):
        super().__init__()
        self.sr_out = sr_out
        self.normalize_amplitude = normalize_amplitude
        self.transpose = transpose
        self.mel = T.MelSpectrogram(
            sample_rate = sr_out,
            n_fft=n_fft,
            win_length=win_length,
            hop_length=hop_length,
            f_min = f_min,
            f_max = f_max,
            center=True,
            pad_mode="reflect",
            power=1.0,
            norm="slaney",
            #onesided=True,
            n_mels=n_mels,
            mel_scale="slaney",
            )
        return
    
    def forward(
            self,
            x: Union[str, torch.Tensor, ArrayLike],
            sr_in: int = None,
            path_out: str = None,
            ):
        # Check if x is a path to a file.
        if isinstance(x, str):
            x, sr_in = torchaudio.load(x)
        elif sr_in is None:
            raise ValueError("Must provide sr_in if x is not a path.")
        # Resample to sr_out
        x = resample_like_librosa(x, sr_in, self.sr_out)
        if self.normalize_amplitude:
            x = normalize_audio_amplitude( x )
        x = self.mel( x )
        x = spectral_normalize_torch( x ).squeeze()
        if self.transpose:
            x = x.T
        if path_out is not None:
            os.makedirs( os.path.dirname( path_out ), exist_ok = True )
            torch.save(x, path_out)
        return x



class Wav2PitchPYIN(torch.nn.Module):
    def __init__(
            self,
            sr_out: int = 100,
            upper_f0_limit: int = 400,
            lower_f0_limit: int = 50,
            log: bool = False,
            ):
        super().__init__()
        self.sr_out = sr_out
        self.upper_f0_limit = upper_f0_limit
        self.lower_f0_limit = lower_f0_limit
        self.log = log

        return
    
    def forward(
            self,
            x: Union[str, torch.Tensor, ArrayLike],
            sr_in: int = None,
            path_out: str = None,
        ) -> torch.Tensor:
        # Check if x is a path to a file.
        if isinstance(x, str):
            x, sr_in = torchaudio.load(x)
        elif sr_in is None:
            raise ValueError("Must provide sr_in if x is not a path.")

        sr_pyin = 22000
        # Resample to sr_pyin
        x = resample_like_librosa(x, sr_in, sr_pyin)



        f0, voiced_flag, voiced_probs = librosa.pyin(
            x.squeeze().numpy(),
            sr = sr_pyin,
            fmin = self.lower_f0_limit,
            fmax = self.upper_f0_limit,
            frame_length = 2048,
            hop_length = 220,#220,
            fill_na = None,
            center = True,
            )
        if self.log:
            f0 = np.log(f0)
        #print(f0)
        # Put f0 andvoiced_probs in a tensor together.
        f0 = torch.tensor( np.array( [f0, voiced_flag, voiced_probs] ).T )
        #print( f0.shape )
        #voiced_probs = torch.tensor(voiced_probs)
        #f0 = torch.stack([f0, voiced_probs], dim=1)
        if path_out is not None:
            os.makedirs( os.path.dirname( path_out ), exist_ok = True )
            torch.save(x, path_out)

        return f0
    
def strictly_increasing(L):
    return all(x<y for x, y in zip(L, L[1:]))

class Wav2Pitch(torch.nn.Module):
    def __init__(
            self,
            sr_out: int = 100,
            upper_f0_limit: int = 400,
            lower_f0_limit: int = 50,
            log: bool = True,
            remove_outliers: bool = False,

            ):
        super().__init__()
        self.sr_out = sr_out
        self.upper_f0_limit = upper_f0_limit
        self.lower_f0_limit = lower_f0_limit
        self.log = log
        self.remove_outliers = remove_outliers

        return
    
    def forward(
            self,
            x: Union[str, torch.Tensor, ArrayLike],
            sr_in: int = None,
            path_out: str = None,
        ) -> torch.Tensor:
        # Check if x is a path to a file.
        if isinstance(x, str):
            x, sr_in = torchaudio.load(x)
        elif sr_in is None:
            raise ValueError("Must provide sr_in if x is not a path.")

        #x_max = np.ceil(x.shape[-1] / sr_in * 100)/100

        pitch_pm = parselmouth.Sound(
            values = x.numpy(),
            sampling_frequency = sr_in
            ).to_pitch()
        pitch_pm_times  = pitch_pm.xs()
        pitch_pm_values = pitch_pm.selected_array[ 'frequency' ]

        #f0 = np.arange( 0, x_max , 0.01)
        pitch_pm_samples = [ x + int( pitch_pm_times[0]*100 ) for x in range( 0, len( pitch_pm_values ) ) ]
        if not strictly_increasing( pitch_pm_samples ):
            print( pitch_pm_samples )
            print( pitch_pm_times )
            raise ValueError("Pitch samples are not strictly increasing.")
        f0 = np.ones( int( np.ceil(x.shape[-1] / sr_in * 100) ) )
        f0[ pitch_pm_samples ] = pitch_pm_values
        f0[ f0 >= self.upper_f0_limit ] = 1
        f0[ f0 <= self.lower_f0_limit ] = 1
        #0 /= self.upper_f0_limit
        if self.log:
            f0 = np.log(f0) / np.log( self.upper_f0_limit )

        f0 = torch.tensor( f0, dtype=torch.float).unsqueeze(1)
        if path_out is not None:
            os.makedirs( os.path.dirname( path_out ), exist_ok = True )
            torch.save( f0, path_out )

        return f0
    


        
        
    

if __name__ == '__main__':
    from asr import AutomaticSpeechRecognizer
    recognizer = AutomaticSpeechRecognizer()
    audio, sr = torchaudio.load( 'data/output.wav' )
    print( audio.shape )
    print( sr )
    result = recognizer( audio, sr_in = sr, graphemes = 'Aber sehen will sie ihn doch', language ='de' )
    print( result )
    stop
    import pandas as pd
    wav2mel = Wav2Mel()
    wav2pitch = Wav2Pitch()
    #wav2pitch_pm = Wav2PitchPM()

    df = pd.read_csv( "datasets/librispeech_dev.csv", sep = "\t" )
    # load first audio file from df
    sample_id = 724
    audio, sr = torchaudio.load( df.audio[sample_id] )
    print( df.audio[sample_id]  )
    print(sr)
    #stop
    # convert to mel spectrogram
    x = wav2mel( audio, sr_in = sr )
    print( x.shape )

    f0 = wav2pitch( audio, sr_in = sr )
    print( f0.shape )
    stop

    f0_pm, data_orig = wav2pitch_pm( df.audio[sample_id] )
    print( f0_pm.shape )
    print( data_orig.shape )
    print( data_orig )

    import matplotlib.pyplot as plt
    plt.plot( [ x * 0.01 for x in range(0, len(f0[:,0])) ], f0[:,0])
    plt.plot( data_orig[:,0], data_orig[:,1] )
    plt.show()
    stop
    fig, axs = plt.subplots(3, 1)
    axs[0].imshow( x.T, origin='lower', aspect='auto' )
    axs[1].plot( f0[:,0] )
    axs[2].plot( f0[:,1] )
    # set x limits to match the mel spectrogram
    axs[1].set_xlim( 0, x.shape[0] )
    axs[2].set_xlim( 0, x.shape[0] )
    # plot f0_pm onto the same axis as f0
    axs[1].plot( f0_pm )
    # plot original data
    axs[1].plot( [ int(x/0.01) for x in data_orig[:,0] ], data_orig[:,1] )

    plt.show()