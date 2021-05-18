#####################################################################################################################################################
#---------------------------------------------------------------------------------------------------------------------------------------------------#
#	- This file is a part of the Python module PyVTL
#---------------------------------------------------------------------------------------------------------------------------------------------------#
#
#	- Copyright (C) 2021, Paul Konstantin Krug, Dresden, Germany
#	- https://github.com/TUD-STKS/PyVTL
#	- Author: Paul Konstantin Krug, TU Dresden
#
#	- License:
#
#		This program is free software: you can redistribute it and/or modify
#		it under the terms of the GNU General Public License as published by
#		the Free Software Foundation, either version 3 of the License, or
#		(at your option) any later version.
#		
#		This program is distributed in the hope that it will be useful,
#		but WITHOUT ANY WARRANTY; without even the implied warranty of
#		MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#		GNU General Public License for more details.
#		
#		You should have received a copy of the GNU General Public License
#		along with this program. If not, see <http://www.gnu.org/licenses/>.
#
#---------------------------------------------------------------------------------------------------------------------------------------------------#
#####################################################################################################################################################
#---------------------------------------------------------------------------------------------------------------------------------------------------#
# pip install praat-parselmouth
# pip install dtw-python
# Requirements:
#	- python 3 (tested with version 3.7)
#	- numpy    (tested with version 1.19.5)
#	- pandas   (tested with version 1.2.1)
#	- scipy    (tested with version 1.6.0)
#
# Optional, used for visualization:
#	- matplotlib        (tested with version 3.3.3)
#	- praat-parselmouth (tested with version 0.3.3)
#---------------------------------------------------------------------------------------------------------------------------------------------------#
#####################################################################################################################################################
#---------------------------------------------------------------------------------------------------------------------------------------------------#


import os
import numpy as np
import pandas as pd
from boundaries import Boundaries
from phonemes import Phoneme

#####################################################################################################################################################
#---------------------------------------------------------------------------------------------------------------------------------------------------#
class Segment_Sequence():
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	"""PyVTL segment sequences""" 
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	def __init__( self, phonemes: list, boundaries: Boundaries ):
		self._boundaries = boundaries
		self._phonemes = [ Phoneme( phone ) for phone in phonemes ]
		self.silence_onset = self._boundaries.times[0]
		self.silence_offset = 0.2
		self.effect = [ None for _ in self._phonemes ]

		#self.data = pd.DataFrame( columns = [ 'onset', 'offset', 'duration', 'index', 'phoneme', 'effect', 'degenerated' ] )
		return
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	@classmethod
	def from_audio_file( self, audio_file_path, phonemes = None, text = None, optimize = False, language = 'de' ): #'TODO'
		if phonemes == None:
			if text == None:
				text = ASR( audio_file_path )
			phonemes = G2P.text_to_sampa( text, language )
		boundary_times = self._get_boundary_times( file_path = audio_file_path, n_phonemes = len( phonemes ) )
		boundaries = Boundaries( boundary_times )
		return cls( phonemes, boundaries )
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	@classmethod
	def from_seg_file( self, seg_file_path ): #'TODO'
		return cls( phonemes, boundaries )
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	def __del__( self ):
		print( 'Segment sequence destroyed.' )
		return
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	def __str__( self ):
		#print( 'now in str function' )
		time = self._boundaries.times
		interval = self._boundaries.intervals
		phonemes = self._phonemes
		data = [ [ time[index], time[index+1], interval[index], phonemes[index].name, self.effect[index] ] for index, _ in enumerate( self._phonemes ) ]
		return str( pd.DataFrame( data, columns=[ 'onset', 'offset', 'duration', 'phoneme', 'effect'] ) )
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	def _get_phoneme_boundary_times( self, file_path: str, n_phonemes: int ):
		data, samplerate = librosa.load( file_path, 44100 )
		onset, offset = detect_onset_and_offset( file_path )
		duration = len( data )
		start = onset
		end =  offset
		print( 'num phon: {}'.format(n_phonemes) )
		print( 'onset: {}'.format( onset ) )
		print( 'offset: {}'.format( offset ) )
		preds = [ x for x in np.arange( start, end, (end-start)/(n_phonemes) ) ] #+1*round( (end-start)/number_phonemes)
		preds.append( offset )
		print( 'returning {} boundaries'.format(len(preds)) )
		return preds
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	def assimilation( self, language = 'de' ): # Assimilation rules for German 
		phonemes = self._phonemes
		for index, phone in enumerate( phonemes ):
			if index >= 2 and phone.name == 'n' and ( phonemes[ index - 1 ].name in ['p','b'] or ( phonemes[ index - 1 ].name == '@' and phonemes[ index - 2 ].name in ['p','b'] ) ):
				self.effect[ index ] = 'm'
			if index >= 2 and phone.name == 'n' and ( phonemes[ index - 1 ].name in ['g','k'] or ( phonemes[ index - 1 ].name == '@' and phonemes[ index - 2 ].name in ['g','k'] ) ):
				self.effect[ index ] = 'N'
		return
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	def elision( self, language = 'de' ): # Elision rules for German
		phonemes = self._phonemes
		for index, phone in enumerate( phonemes ):
			if phone.name == '@' and phonemes[ index + 1 ].is_sonorant():
				self.effect[ index ] = 'elision'
		return
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	def export_seg( self, file_path, account_for_effects = False ):
		assert len( self._phonemes ) == len( self._boundaries.times ) - 1, 'Lengths do not match'
		out_file = open( file_path, 'w+' )
		out_file.write( 'name = {}; duration_s = {};\n'.format( '', self.silence_onset ) )
		for index, phoneme in enumerate( self._phonemes ):
			if account_for_effects == False or ( account_for_effects == True and self.effect[ index ] != 'elision' ):
				if self.effect[ index ] != None and self.effect[ index ] != 'elision':
					output_phone = self.effect[ index ]
				else:
					output_phone = phoneme.name
				out_file.write( 'name = {}; duration_s = {};\n'.format( output_phone, self._boundaries.intervals[ index ] ) )
		out_file.write( 'name = {}; duration_s = {};\n'.format( '', self.silence_offset ) )
		return
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	def import_seg( self, file_path ):
		return
#---------------------------------------------------------------------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------------------------------------------------------------------#
#####################################################################################################################################################