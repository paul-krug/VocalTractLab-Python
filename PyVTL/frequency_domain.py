import PyVTL.VocalTractLabApi as vtl
import librosa
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks


def multiple_formatter( denominator=2, number=np.pi, latex='\\pi' ):
	def gcd(a, b):
		while b:
			a, b = b, a%b
		return a
	def _multiple_formatter(x, pos):
		den = denominator
		num = np.int(np.rint(den*x/number))
		com = gcd(num,den)
		(num,den) = (int(num/com),int(den/com))
		if den==1:
			if num==0:
				return r'$0$'
			if num==1:
				 return r'$%s$'%latex
			elif num==-1:
				return r'$-%s$'%latex
			else:
				return r'$%s%s$'%(num,latex)
		else:
			if num==1:
				return r'$\frac{%s}{%s}$'%(latex,den)
			elif num==-1:
				return r'$\frac{-%s}{%s}$'%(latex,den)
			else:
				return r'$\frac{%s%s}{%s}$'%(num,latex,den)
	return _multiple_formatter
class Multiple:
	def __init__(self, denominator=2, number=np.pi, latex='\\pi'):
		self.denominator = denominator
		self.number = number
		self.latex = latex
	def locator(self):
		return plt.MultipleLocator(self.number / self.denominator)
	def formatter(self):
		return plt.FuncFormatter(multiple_formatter(self.denominator, self.number, self.latex))




#####################################################################################################################################################
class Transfer_Function():
	def __init__( self, 
		          magnitude_spectrum: np.ndarray,
		          phase_spectrum: np.ndarray,
		          n_spectrum_samples: int,
		          name: str = 'transfer_function'
		          ):
		if not isinstance( n_spectrum_samples, int ):
			raise ValueError( 'n_spectrum_samples must be an integer and should be a power of 2! Passed type is: {}'.format( type( n_spectrum_samples ) ) )
		self.constants = vtl.get_constants()
		self.delta_frequency = self.constants[ 'samplerate' ] / n_spectrum_samples
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
		self.data = dict( frequency = self.magnitude_spectrum, phase = self.phase_spectrum )
		self.formants = self.get_formants()
		self.f1, self.f2, self.f3, self.f4 = self.formants
		return
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	def get_formants( self, peak_distance = 1, sr = 44100 ):
		sr = self.constants[ 'samplerate' ]
		peaks, _ = find_peaks( self.magnitude_spectrum, distance = peak_distance )
		peaks = [ peak * sr/self.n_spectrum_samples for peak in peaks ]
		while peaks[ 0 ] < 100:
			del peaks[ 0 ]
		if len( peaks ) < 4:
			peaks.extend( [ None for _ in range( 0, 4 - len( peaks ) ) ] )
		elif len( peaks ) > 4:
			peaks = peaks[ : 4 ]
		return peaks
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	def plot( self, parameters = [ 'frequency', 'phase' ], plot_formants = True ): #, scale = 'dB' ):

		figure, axs = plt.subplots( len(parameters), figsize = (8, 4/3 *len(parameters) ), sharex = True, gridspec_kw = {'hspace': 0} )
		#figure.suptitle( 'Sharing both axes' )
		#parameters = self.tract.columns

		for index, parameter in enumerate( parameters ):
			if parameter == 'frequency':
				y = librosa.amplitude_to_db( self.data[ parameter ] )
				y_title = 'Intensity [dB]'
			elif parameter == 'phase':
				y = self.data[ parameter ]
				y_title = 'Phase'
				axs[ index ].yaxis.set_major_locator(plt.MultipleLocator(np.pi / 2))
				#axs[ index ].yaxis.set_minor_locator(plt.MultipleLocator(np.pi / 12))
				axs[ index ].yaxis.set_major_formatter(plt.FuncFormatter(multiple_formatter()))
			else:
				raise ValueError( 'parameters must be frequency and/or phase! Passed values are: {}'.format( parameters ) )

			axs[ index ].plot( np.arange( 0, self.n_spectrum_samples, self.delta_frequency ), y )
			axs[ index ].set( ylabel = y_title )
		plt.xlabel( 'Frequency [Hz]' )
		if plot_formants:
			for formant in self.formants:
				for ax in axs:
					ax.axvline( formant, color = 'gray', ls = '--' )
		for ax in axs:
		    ax.label_outer()
		figure.align_ylabels( axs[:] )
		plt.show()
		return