


from vocaltractlab_cython import get_constants
#import target_approximation.utils as PT
from target_approximation.utils import finalize_plot
from target_approximation.utils import get_plot
from target_approximation.utils import get_plot_limits
#import librosa
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks

from vocaltractlab.utils import multiple_formatter
from vocaltractlab.audioprocessing import amplitude_to_db



class TransferFunction():
    def __init__(
            self, 
            tract_state: np.ndarray,
            magnitude_spectrum: np.ndarray,
            phase_spectrum: np.ndarray,
            n_spectrum_samples: int,
            #name: str = 'transfer_function'
            ):
        if not isinstance( n_spectrum_samples, int ):
            raise ValueError(
                f"""
                Argument n_spectrum_samples must be an integer
                and should be a power of 2, but you passed:
                {n_spectrum_samples}
                """
                )
        self.constants = get_constants()
        self.tract_state = tract_state
        self.delta_frequency = self.constants[ 'sr_audio' ] / n_spectrum_samples
        max_bin = round( n_spectrum_samples / self.delta_frequency )
        self.n_spectrum_samples = n_spectrum_samples
        if isinstance( magnitude_spectrum, np.ndarray ):
            self.magnitude_spectrum = magnitude_spectrum[ : max_bin ]
        else:
            self.magnitude_spectrum = None
        if isinstance( phase_spectrum, np.ndarray ):
            self.phase_spectrum = phase_spectrum[ : max_bin ]
        else:
            self.phase_spectrum = None
        self.data = dict(
            frequency = self.magnitude_spectrum,
            phase = self.phase_spectrum,
            )
        self.formants = self.get_formants()
        self.f1, self.f2, self.f3, self.f4 = self.formants
        return
    
    @classmethod
    def from_dict(
            cls,
            x,
            ):
        return cls(
            tract_state = x[ 'tract_state' ],
            magnitude_spectrum = x[ 'magnitude_spectrum' ],
            phase_spectrum = x[ 'phase_spectrum' ],
            n_spectrum_samples = x[ 'n_spectrum_samples' ],
            )

    def get_formants(
            self,
            peak_distance = 1,
            # = 44100,
            ):
        sr = self.constants[ 'sr_audio' ]
        peaks, _ = find_peaks(
            self.magnitude_spectrum,
            distance = peak_distance,
            )
        peaks = [
            peak * sr/self.n_spectrum_samples
            for peak in peaks
            ]
        while peaks[ 0 ] < 100:
            del peaks[ 0 ]
        if len( peaks ) < 4:
            peaks.extend( [ 
                None for _ in range( 0, 4 - len( peaks ) )
                ] )
        elif len( peaks ) > 4:
            peaks = peaks[ : 4 ]
        return peaks

    def plot( self,
              parameters = [ 'frequency', 'phase' ],
              plot_formants = True,
              axs: list = None,
              plot_kwargs: list = [ dict( color = 'navy' ), dict( color = 'darkorange' ) ],
              **kwargs,
              ): #, scale = 'dB' ):
        figure, axs = get_plot( n_rows = len( parameters ), axs = axs )
        for index, parameter in enumerate( parameters ):
            if parameter == 'frequency':
                y = amplitude_to_db( self.data[ parameter ] )
                continuities = [ slice( 0, len(y) ) ]
                y_title = 'Intensity [dB]'
                #_min = np.min( y )
                #_max = np.max( y )
                #axs[ index ].set( ylim = [ _min - 0.1 * np.abs( _max - _min ), _max + 0.1 * np.abs( _max - _min ) ] )
                axs[ index ].set( ylim = get_plot_limits( y ) )
                axs[ index ].locator_params( axis = 'y', nbins = 4 )
            elif parameter == 'phase':
                continuities = []
                y = self.data[ parameter ]
                tmp_idx = 0
                for idx in range( 0, len(y) - 1 ):
                    if np.abs( y[idx] - y[idx+1] ) > 1.5:
                        continuities.append( slice( tmp_idx, idx+1 ) )
                        tmp_idx = idx + 1
                if tmp_idx != len( y ):
                    continuities.append( slice( tmp_idx, len( y ) ) )

                y = self.data[ parameter ]
                y_title = 'Phase'
                axs[ index ].yaxis.set_major_locator(plt.MultipleLocator(np.pi / 2))
                #axs[ index ].yaxis.set_minor_locator(plt.MultipleLocator(np.pi / 12))
                axs[ index ].yaxis.set_major_formatter(plt.FuncFormatter(multiple_formatter()))
                axs[ index ].set( ylim = [ -3.76, 3.76 ] )
            else:
                raise ValueError( 'parameters must be frequency and/or phase! Passed values are: {}'.format( parameters ) )
            x = np.arange( 0, self.n_spectrum_samples, self.delta_frequency )
            for _slice in continuities:
                axs[ index ].plot( x[ _slice ], y[ _slice ], **plot_kwargs[ index ] )
            axs[ index ].set( ylabel = y_title )
        plt.xlabel( 'Frequency [Hz]' )
        if plot_formants:
            for formant in self.formants:
                for ax in axs:
                    ax.axvline( formant, color = 'gray', ls = '--' )
        for ax in axs:
            ax.label_outer()
        finalize_plot( figure, axs, **kwargs )
        return axs