


import os
import numpy as np

from vocaltractlab_cython.VocalTractLabApi import _close
from vocaltractlab_cython.VocalTractLabApi import _initialize
from vocaltractlab_cython import get_constants
from vocaltractlab_cython import gesture_file_to_audio
from vocaltractlab_cython import gesture_file_to_motor_file
from vocaltractlab_cython import phoneme_file_to_gesture_file
from vocaltractlab_cython import synth_block
from vocaltractlab_cython import tract_state_to_limited_tract_state
from vocaltractlab_cython import tract_state_to_transfer_function
from vocaltractlab_cython import tract_state_to_tube_state
#from vocaltractlab_cython.exceptions import VTLAPIError
from target_approximation import TargetSeries
from target_approximation.vocaltractlab import MotorSequence
from target_approximation.vocaltractlab import MotorSeries
from target_approximation.vocaltractlab import SupraGlottalSequence
from target_approximation.vocaltractlab import SupraGlottalSeries

from typing import Union, List, Tuple, Dict, Any, Optional, Callable, Iterable, Sequence
from numpy.typing import ArrayLike

from tools_mp import multiprocess

from .utils import make_iterable
from .audioprocessing import audio_to_f0
from .audioprocessing import postprocess
from .frequency_domain import TransferFunction
from .tube_state import TubeState



def limit(
        x: Union[
            #MotorSequence,
            MotorSeries,
            #SupraGlottalSequence,
            SupraGlottalSeries,
            str,
            ],
        workers: int = None,
        verbose: bool = True,
        ):
    if isinstance( x, MotorSequence ):
        ms = x.to_series()
        sgs = ms.tract()
    elif isinstance( x, MotorSeries ):
        sgs = x.tract()
    elif isinstance( x, SupraGlottalSequence ):
        sgs = x.to_series()
    elif isinstance( x, str ):
        sgs = SupraGlottalSeries.load( x )
    elif isinstance( x, SupraGlottalSeries ):
        sgs = x
    else:
        raise TypeError(
            f"""
            The specified data type: '{type(x)}'
            is not supported. Type must be one of the following:
            - MotorSequence
            - MotorSeries
            - SupraGlottalSequence
            - SupraGlottalSeries
            - str
            """
            )
    
    args = [
        dict(
            tract_state = ts,
            )
        for ts in sgs.to_numpy( transpose = False )
        ]
    
    states = multiprocess(
        tract_state_to_limited_tract_state,
        args = args,
        return_data = True,
        workers = workers,
        verbose = verbose,
        #mp_threshold = 4, # TODO: not implemented yet
        )
    
    states = np.array( states )
    lim = SupraGlottalSeries( states )
    if isinstance( x, MotorSeries ):
        lim = MotorSeries( lim & states.glottis() )
    
    return lim
    
    return

def load_speaker(
        speaker: str,
        ) -> None:
    if not speaker.endswith( '.speaker' ):
        speaker = f"{speaker}.speaker"
    _close()
    # check if speaker is a valid file path
    if os.path.exists( speaker ):
        speaker_path = speaker
    else:
        speaker_path = os.path.join(
            os.path.dirname( __file__ ),
            'speaker',
            speaker,
            )
        if not os.path.exists( speaker_path ):
            raise FileNotFoundError(
                f"""
                The specified speaker file path: '{speaker}'
                does not exist.
                """
                )
    _initialize( speaker_path )
    return

def speakers() -> List[ str ]:
    speaker_dir = os.path.join(
        os.path.dirname( __file__ ),
        'speaker',
        )
    speakers = [
        os.path.basename(f)
        for f in os.listdir( speaker_dir )
        if f.endswith( '.speaker' )
        ]
    return speakers

def gesture_to_audio(
        x: Union[ Iterable[ str ], str ],
        audio_files: Optional[ Union[ Iterable[ str ], str ] ],
        normalize_audio: int = -1,
        sr: int = None,
        return_data: bool = False,
        workers: int = None,
        verbose: bool = True,
        ) -> None:

    gesture_files = make_iterable( x )
    audio_files = make_iterable( audio_files )
    if len( gesture_files ) != len( audio_files ):
        raise ValueError(
            f"""
            The number of gesture file paths: {len(gesture_files)}
            does not match the number of audio file paths: {len(audio_files)}.
            """
            )
    
    args = [
        dict(
            gesture_data = gf,
            audio_file_path = af,
            verbose_api = False,
            normalize_audio = normalize_audio,
            sr = sr,
            )
        for gf, af in zip(
            gesture_files,
            audio_files,
            )
        ]
    audio_data = multiprocess(
        _gesture_to_audio,
        args = args,
        return_data = return_data,
        workers = workers,
        verbose = verbose,
        #mp_threshold = 4,
        )
    return audio_data

def _gesture_to_audio(
        gesture_data,
        audio_file_path,
        verbose_api,
        normalize_audio,
        sr,
        ) -> np.ndarray:
    if isinstance( gesture_data, str ):
        #gesture_file = gesture_data.to_gesture_file( file_path = None )
        gesture_file = gesture_data
    else:
        raise TypeError(
            f"""
            The specified gesture data type: '{type(gesture_data)}'
            is not supported. Type must be one of the following:
            - str
            """
            )
    audio = gesture_file_to_audio(
        ges_file_path = gesture_file,
        audio_file_path = None,
        verbose_api = verbose_api,
    )
    
    audio = postprocess(
        x = audio,
        sr_out = sr,
        dBFS = normalize_audio,
        file_path = audio_file_path,
        to_numpy = True,
        )
    
    return audio

def gesture_to_motor(
        gesture_files: Union[ Iterable[ str ], str ],
        motor_files: Optional[ Union[ Iterable[ str ], str ] ],
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
            gesture_file = gf,
            motor_file = mf,
            )
        for gf, mf in zip(
            gesture_files,
            motor_files,
            )
        ]
    multiprocess(
        gesture_file_to_motor_file,
        args = args,
        return_data = False,
        workers = workers,
        verbose = verbose,
        #mp_threshold = 4,
        )
    return

def motor_to_audio(
        motor_data: Union[ MotorSequence, MotorSeries, str ],
        audio_files: Optional[ Union[ Iterable[str], str ] ] = None,
        normalize_audio: int = -1,
        sr: int = None,
        return_data: bool = False,
        workers: int = None,
        verbose: bool = True,
        ) -> np.ndarray:
    """
    Convert motor data into audio signals.

    Parameters
    ----------
    motor_data : Union[MotorScore, MotorSeries, str]
        Input data representing motor scores or series.
        Can be a MotorScore object, MotorSeries object, or a path to a file.

    audio_files : Optional[Union[Iterable[str], str]], optional
        Path or list of paths to store the generated audio files.
        If None, audio files will not be saved. Default is None.

    normalize_audio : int, optional
        Amplitude normalization factor.
        -1 indicates no normalization. Default is -1.

    sr : int, optional
        Sampling rate of the output audio.
        If None, defaults to the system's default audio sampling rate.

    return_data : bool, optional
        Flag indicating whether to return the generated audio data.
        Default is False.

    workers : int, optional
        Number of worker processes for parallel processing.
        If None, uses the system's default number of CPU cores.
        Default is None.

    verbose : bool, optional
        Verbosity mode. If True, displays progress information.
        Default is True.

    Returns
    -------
    np.ndarray
        If 'return_data' is True, returns a NumPy array of the generated audio data.

    Raises
    ------
    ValueError
        If the number of audio file paths doesn't match the number of motor data.

    FileNotFoundError
        If the specified motor file path does not exist.

    TypeError
        If the specified motor data type is not supported.
        Supported types include str, MotorScore, and MotorSeries.

    Notes
    -----
    This function converts motor data into audio signals using the VocalTractLab synthesizer.
    It processes the motor parameters and generates corresponding audio signals.

    Examples
    --------
    # Example 1: Convert MotorScore object to audio without saving files
    >>> motor_data = MotorScore(...)  # Replace '...' with actual MotorScore data
    >>> audio_data = motor_to_audio(motor_data)

    # Example 2: Convert MotorSeries object to audio and save the files
    >>> motor_series = MotorSeries(...)  # Replace '...' with actual MotorSeries data
    >>> audio_files = ['audio1.wav', 'audio2.wav']  # List of paths to save audio files
    >>> motor_to_audio(motor_series, audio_files=audio_files, return_data=False)

    # Example 3: Convert motor data from a file to audio with normalization
    >>> motor_file_path = 'path/to/motor_data.csv'  # Replace with the actual file path
    >>> audio_data = motor_to_audio(motor_file_path, normalize_audio=0.5, return_data=True)
    """

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
    """
    Generate audio from motor data.

    Parameters
    ----------
    motor_data : Union[MotorScore, MotorSeries, str]
        Input data representing motor scores or series.
        Can be a MotorScore object, MotorSeries object, or a path to a file.

    audio_file_path : Optional[str]
        Path to store the generated audio file. If None, audio will not be saved.

    normalize_audio : int
        Amplitude normalization factor. Use -1 for no normalization.

    sr : int
        Sampling rate of the output audio.

    state_samples : int, optional
        Number of samples for state duration.
        If None, defaults to a predefined constant value.

    Returns
    -------
    torch.Tensor
        A tensor representing the generated audio.

    Raises
    ------
    FileNotFoundError
        If the specified motor file path does not exist.

    TypeError
        If the specified motor data type is not supported.
        Supported types include str, MotorScore, and MotorSeries.

    Notes
    -----
    This function generates audio signals from motor data using the VocalTractLab synthesizer.
    It processes the motor parameters and synthesizes corresponding audio signals.

    Examples
    --------
    # Example 1: Generate audio from MotorScore object without saving the file
    >>> motor_data = MotorScore(...)  # Replace '...' with actual MotorScore data
    >>> audio_tensor = _motor_to_audio(motor_data, audio_file_path=None, normalize_audio=0, sr=44100)

    # Example 2: Generate audio from MotorSeries object and save the audio file
    >>> motor_series = MotorSeries(...)  # Replace '...' with actual MotorSeries data
    >>> audio_path = 'output_audio.wav'  # Path to save the audio file
    >>> _motor_to_audio(motor_series, audio_file_path=audio_path, normalize_audio=-1, sr=22050)

    # Example 3: Generate audio from a file containing motor data with custom state samples
    >>> motor_file_path = 'path/to/motor_data.csv'  # Replace with the actual file path
    >>> audio_tensor = _motor_to_audio(motor_file_path, audio_file_path=None, normalize_audio=0.8, sr=44100, state_samples=120)
    """

    if isinstance( motor_data, str ):
        if not os.path.exists( motor_data ):
            raise FileNotFoundError( 
                f"""
                The specified motor file path: '{motor_data}'
                does not exist.
                """
            )
        motor_series = MotorSeries.load(
            motor_data,
            sr = 441,
            )
    elif isinstance( motor_data, MotorSequence ):
        motor_series = motor_data.to_series(
            sr = 441,
        )
    elif isinstance( motor_data, MotorSeries ):
        motor_series = motor_data
    else:
        raise TypeError(
            f"""
            The specified motor data type: '{type(motor_data)}'
            is not supported. Type must be one of the following:
            - str
            - MotorSequence
            - MotorSeries
            """
            )
    if motor_series.sr is None:
        raise ValueError(
            f"""
            The specified motor series has no asociated sampling
            rate and thus, cannot be used for audio generation.
            Please ensure that the sampling rate is set before
            generating audio.
            """
            )
    vtl_constants = get_constants()
    if state_samples is None:
        #state_samples = vtl_constants[ 'n_samples_per_state' ]
        state_samples = int(
            vtl_constants[ 'sr_audio' ] / motor_series.sr
        )
        

    #print( motor_series.to_numpy( part='tract' ) )

    tract_params = motor_series.tract().to_numpy( transpose = False )
    glottal_params = motor_series.glottis().to_numpy( transpose = False )
    #print( tract_params.shape )
    #print( glottal_params.shape )
    #print( state_samples )

    
    audio = synth_block(
        tract_parameters = tract_params,
        glottis_parameters = glottal_params,
        state_samples = state_samples,
        verbose_api = False,
        )
    
    audio = postprocess(
        x = audio,
        sr_out = sr,
        dBFS = normalize_audio,
        file_path = audio_file_path,
        to_numpy = True,
        )
    
    return audio

def motor_to_transfer_function(
        x: Union[
            MotorSequence,
            MotorSeries,
            SupraGlottalSequence,
            SupraGlottalSeries,
            str,
            ],
        n_spectrum_samples: int = 8192,
        save_magnitude_spectrum: bool = True,
        save_phase_spectrum: bool = True,
        workers: int = None,
        verbose: bool = True,
        ):
    if isinstance( x, MotorSequence ):
        ms = x.to_series()
        sgs = ms.glottis()
    elif isinstance( x, MotorSeries ):
        sgs = x.glottis()
    elif isinstance( x, SupraGlottalSequence ):
        sgs = x.to_series()
    elif isinstance( x, str ):
        sgs = SupraGlottalSeries.load( x )
    elif isinstance( x, SupraGlottalSeries ):
        sgs = x
    else:
        raise TypeError(
            f"""
            The specified data type: '{type(x)}'
            is not supported. Type must be one of the following:
            - MotorSequence
            - MotorSeries
            - SupraGlottalSequence
            - SupraGlottalSeries
            - str
            """
            )
    args = [
        dict(
            tract_state = ts,
            n_spectrum_samples = n_spectrum_samples,
            save_magnitude_spectrum = save_magnitude_spectrum,
            save_phase_spectrum = save_phase_spectrum,
            )
        for ts in sgs.to_numpy( transpose = False )
        ]
    
    trf_data = multiprocess(
        _motor_to_transfer_function,
        args = args,
        return_data = True,
        workers = workers,
        verbose = verbose,
        #mp_threshold = 4, # TODO: not implemented yet
        )
    
    return trf_data

def _motor_to_transfer_function( **kwargs ):
    x = tract_state_to_transfer_function( **kwargs )
    x[ 'tract_state' ] = kwargs[ 'tract_state' ]
    return TransferFunction.from_dict( x )

def motor_to_tube(
        x: Union[
            MotorSequence,
            MotorSeries,
            SupraGlottalSequence,
            SupraGlottalSeries,
            str,
            ],
	    save_tube_length: bool = True,
	    save_tube_area: bool = True,
	    save_tube_articulator: bool = True,
	    save_incisor_position: bool = True,
	    save_tongue_tip_side_elevation: bool = True,
	    save_velum_opening: bool = True,
	    fast_calculation = True,
	    workers: int = None,
        verbose: bool = True,
        ) -> np.ndarray:
    
    if isinstance( x, MotorSequence ):
        ms = x.to_series()
        sgs = ms.glottis()
    elif isinstance( x, MotorSeries ):
        sgs = x.glottis()
    elif isinstance( x, SupraGlottalSequence ):
        sgs = x.to_series()
    elif isinstance( x, str ):
        sgs = SupraGlottalSeries.load( x )
    elif isinstance( x, SupraGlottalSeries ):
        sgs = x
    else:
        raise TypeError(
            f"""
            The specified data type: '{type(x)}'
            is not supported. Type must be one of the following:
            - MotorSequence
            - MotorSeries
            - SupraGlottalSequence
            - SupraGlottalSeries
            - str
            """
            )
    
    args = [
        dict(
            tract_state = ts,
            fast_calculation = fast_calculation,
            save_tube_length = save_tube_length,
            save_tube_area = save_tube_area,
            save_tube_articulator = save_tube_articulator,
            save_incisor_position = save_incisor_position,
            save_tongue_tip_side_elevation = save_tongue_tip_side_elevation,
            save_velum_opening = save_velum_opening,
            )
        for ts in sgs.to_numpy( transpose = False )
        ]
    
    tube_data = multiprocess(
        _motor_to_tube,
        args = args,
        return_data = True,
        workers = workers,
        verbose = verbose,
        #mp_threshold = 4, # TODO: not implemented yet
        )
    
    return tube_data

def _motor_to_tube( **kwargs ):
    x = tract_state_to_tube_state( **kwargs )
    x[ 'tract_state' ] = kwargs[ 'tract_state' ]
    return TubeState.from_dict( x )

def phoneme_to_audio(
        x: List[ str ],
        gesture_files: List[ str ],
        motor_files: List[ str ],
        f0_files: Optional[ List[ str ] ] = None,
        motor_f0_files: Optional[ List[ str ] ] = None,
        audio_files: Optional[ List[ str ] ] = None,
        normalize_audio = -1,
        sr = None,
        return_data = False,
        workers: int = None,
        verbose: bool = True,
        ):
    
    phoneme_to_motor(
        x = x,
        gesture_files = gesture_files,
        motor_files = motor_files,
        workers = workers,
        verbose = verbose,
        )
    
    if f0_files is not None:
        if motor_f0_files is None:
            ms_data = augment_motor_f0(
                motor_files = motor_files,
                f0_files = f0_files,
                out_files = motor_f0_files,
                return_data = True,
                workers = workers,
                verbose = verbose,
                )
        else:
            augment_motor_f0(
                motor_files = motor_files,
                f0_files = f0_files,
                out_files = motor_f0_files,
                return_data = False,
                workers = workers,
                verbose = verbose,
                )
            ms_data = motor_f0_files

    else:
        ms_data = motor_files

    audio_data = motor_to_audio(
        motor_data = ms_data,
        audio_files = audio_files,
        normalize_audio = normalize_audio,
        sr = sr,
        return_data = return_data,
        workers = workers,
        verbose = verbose,
        )
    
    return audio_data

def phoneme_to_gesture(
        x: List[ str ],
        gesture_files: List[ str ],
        workers: int = None,
        verbose: bool = True,
        ) -> np.ndarray:
    phoneme_files = make_iterable( x )
    # TODO: implement phn sequence to phn file
    gesture_files = make_iterable( gesture_files )
    if len( phoneme_files ) != len( gesture_files ):
        raise ValueError(
            f"""
            The number of phoneme file paths: {len(phoneme_files)}
            does not match the number of gesture file paths: {len(gesture_files)}.
            """
            )
    
    args = [
        dict(
            phoneme_file = pf,
            gesture_file = gf,
            verbose_api = False,
            )
        for pf, gf in zip(
            phoneme_files,
            gesture_files,
            )
        ]
    multiprocess(
        phoneme_file_to_gesture_file,
        args = args,
        return_data = False,
        workers = workers,
        verbose = verbose,
        #mp_threshold = 4,
        )
    return

def phoneme_to_motor(
        x: List[ str ],
        gesture_files: List[ str ],
        motor_files: List[ str ],
        workers: int = None,
        verbose: bool = True,
        ):
    
    phoneme_to_gesture(
        x = x,
        gesture_files = gesture_files,
        workers = workers,
        verbose = verbose,
        )
    
    gesture_to_motor(
        gesture_files = gesture_files,
        motor_files = motor_files,
        workers = workers,
        verbose = verbose,
        )
    
    return

def augment_motor_f0(
        motor_files: Union[ Iterable[ str ], str ],
        f0_files: Union[ Iterable[ str ], str ],
        out_files: Optional[ Union[ Iterable[ str ], str ] ] = None,
        target_sr: int = 441,
        return_data: bool = False,
        workers: int = None,
        verbose: bool = True,
        **kwargs,
        ):
    motor_files = make_iterable( motor_files )
    f0_files = make_iterable( f0_files )
    if len( motor_files ) != len( f0_files ):
        raise ValueError(
            f"""
            The number of motor file paths: {len(motor_files)}
            does not match the number of f0 file paths: {len(f0_files)}.
            """
            )
    if out_files is not None:
        out_files = make_iterable( out_files )
        if len( motor_files ) != len( out_files ):
            raise ValueError(
                f"""
                The number of motor file paths: {len(motor_files)}
                does not match the number of output file paths: {len(out_files)}.
                """
                )
    
    args = [
        dict(
            motor_file = mf,
            f0_file = ff,
            out_file = of,
            target_sr = target_sr,
            **kwargs,
            )
        for mf, ff, of in zip(
            motor_files,
            f0_files,
            out_files,
            )
        ]
    
    ms_data = multiprocess(
        _augment_motor_f0,
        args = args,
        return_data = return_data,
        workers = workers,
        verbose = verbose,
        #mp_threshold = 4,
        )
    return ms_data

def _augment_motor_f0(
        motor_file,
        f0_file,
        out_file,
        target_sr,
        **kwargs,
        ):
    ms = MotorSeries.load( motor_file )
    ms.resample( target_sr = target_sr )

    _, feature = audio_to_f0( f0_file )
    f0 = feature[ :, 0 ]
    tgss = TargetSeries(
        series = f0,
        sr = 100,
        tiers = [ 'F0' ],
        )
    tgss.resample( target_sr = target_sr )
    ms = ms & tgss

    if out_file is not None:
        ms.save( out_file, **kwargs )
    return ms


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