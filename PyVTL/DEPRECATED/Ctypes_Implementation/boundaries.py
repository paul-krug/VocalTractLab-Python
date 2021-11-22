import numpy as np
import pandas as pd
import librosa

from PyVTL import voice_activity_detection as vad


#####################################################################################################################################################
#---------------------------------------------------------------------------------------------------------------------------------------------------#
class Boundaries():
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	"""PyVTL phoneme boundaries""" 
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	def __init__( self, boundary_times ):
		self.initial_times = boundary_times
		self.initial_intervals = self.get_intervals( self.initial_times )
		self.times, self.degenerated_intervals = self.get_regularized_times( self.initial_times )
		self.intervals = self.get_intervals( self.times )
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	@classmethod
	def from_durations( cls, boundary_durations ): #TODO
		return cls( boundary_times )
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	@classmethod
	def uniform_from_audio_file( cls, audio_file_path: str, number_phonemes: int ):
		data, samplerate = librosa.load( audio_file_path, 44100 )
		onset, offset = vad.detect_onset_and_offset( audio_file_path )
		#duration = len( data )
		start = onset
		end =  offset
		print( 'num phon: {}'.format(number_phonemes) )
		print( 'onset: {}'.format( onset ) )
		print( 'offset: {}'.format( offset ) )
		boundary_times = [ x for x in np.arange( start, end, (end-start)/(number_phonemes) ) ] #+1*round( (end-start)/number_phonemes)
		boundary_times.append( offset )
		print( 'returning {} boundaries'.format(len(boundary_times)) )
		return cls( boundary_times )
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	def __str__( self ):
		time = self.times
		interval = self.intervals
		data = [ [ time[index], time[index+1], interval[index] ] for index, _ in enumerate( self.intervals ) ]
		return str( pd.DataFrame( data, columns=[ 'onset', 'offset', 'duration'] ) )
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	def get_intervals( self, boundary_times ):
		return [ boundary_times[ index + 1 ] - boundary_times[ index ] for index in range( 0, len( boundary_times ) - 1 ) ]
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	def get_regularized_times( self, boundary_times, regularization = 0.02 ):
		intervals = self.get_intervals( boundary_times )
		degenerated_intervals = [ index for index, interval in enumerate( intervals ) if interval < regularization ]
		times = []
		for index, time in enumerate( boundary_times ):
			if index == 0:
				times.append( boundary_times[ index ] )
			else:
				if ( boundary_times[ index ] - times[ index - 1 ] ) > regularization:
					times.append( boundary_times[ index ] )
				else:
					times.append( times[ index - 1 ] + regularization )
		return [ times, degenerated_intervals ]
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	def time_to_samples( self, n_audio_samples, audio_samplerate = 44100, hop_length = 256 ):
		n_time_samples = round( n_audio_samples / hop_length )
		boundary_function = [ 0 for x in range(0, n_time_samples + 1 ) ]
		for time in self.times:
			#print( boundary )
			#print( round( boundary / ( mfcc_settings['hop_length'] / mfcc_settings['sr'] ) ) )
			time_sample_index = round( time * audio_samplerate / hop_length )
			boundary_function[ int( time_sample_index ) ] = 1
			boundary_function = np.array( boundary_function )
		print( 'Boundary function has shape: {}'.format( boundary_function.shape ) )
		return boundary_function
#####################################################################################################################################################