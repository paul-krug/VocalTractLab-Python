#####################################################################################################################################################
#---------------------------------------------------------------------------------------------------------------------------------------------------#
#	- This file is a part of the VocalTractLab Python module PyVTL, see https://github.com/paul-krug/VocalTractLab
#---------------------------------------------------------------------------------------------------------------------------------------------------#
#
#	- Copyright (C) 2021, Paul Konstantin Krug, Dresden, Germany
#	- https://github.com/paul-krug/VocalTractLab
#	- Author: Paul Konstantin Krug, TU Dresden
#
#	- License info:
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

#---------------------------------------------------------------------------------------------------------------------------------------------------#
#####################################################################################################################################################
#---------------------------------------------------------------------------------------------------------------------------------------------------#
# Load essential packages:
import os
import warnings
import pandas as pd
import numpy as np
import VocalTractLab.VocalTractLabApi as vtl
import VocalTractLab.function_tools as FT
import VocalTractLab.plotting_tools as PT
from VocalTractLab.plotting_tools import finalize_plot
from VocalTractLab.plotting_tools import get_plot
from VocalTractLab.plotting_tools import get_plot_limits
import matplotlib.pyplot as plt
from  itertools import chain
#import math
#---------------------------------------------------------------------------------------------------------------------------------------------------#
#####################################################################################################################################################


#####################################################################################################################################################
#---------------------------------------------------------------------------------------------------------------------------------------------------#
class Segment_Sequence():
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	"""PyVTL segment sequences""" 
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	def __init__( self, phonemes: list, durations: list, onset_duration: float = 0, offset_duration: float = None ):
		self.durations = FT.check_if_list_is_valid( durations, (int, float, np.float64) )
		self.phonemes = FT.check_if_list_is_valid( phonemes, (str) )
		self.onset_duration = onset_duration
		self.offset_duration = offset_duration
		self.effects = [ None for _ in self.phonemes ]
		self.length = self.onset_duration + np.sum( self.durations ) # + offset_uration if defined

		#self.data = pd.DataFrame( columns = [ 'onset', 'offset', 'duration', 'index', 'phoneme', 'effect', 'degenerated' ] )
		return
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	@classmethod
	def from_audio_file( cls, audio_file_path, phonemes = None, text = None, optimize = False, language = 'de' ): #'TODO'
		if phonemes == None:
			if text == None:
				text = ASR( audio_file_path )
			phonemes = G2P.text_to_sampa( text, language )
		boundaries = Boundaries.uniform_from_audio_file( audio_file_path = audio_file_path, n_phonemes = len( phonemes ) )
		return cls( phonemes, boundaries )
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	@classmethod
	def from_seg_file( cls, seg_file_path ): #'TODO'
		phonemes = []
		durations = []
		#boundary_times = []
		with open( seg_file_path ) as file:
			for line in file:
				line = list( filter( None, line.split(' ') ) )
				print( line )
				stop
				if len(line.strip()) != 0 :
					items = line.strip().split(';')
					for item in [x for x in items if x]:
						label = item.split('=')[0].strip()
						value = item.split('=')[1].strip()
						#print('Label: {}'.format(label))
						#print('Value: {}'.format(value))
						if label =='name':# and (value not in [None," ", ""]):
							phonemes.append( value )
						if label == 'duration_s':
							durations.append( float( value ) )
		return cls( phonemes, durations )
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	def __del__( self ):
		#print( 'Segment sequence destroyed.' )
		return
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	def __str__( self ):
		return str( self._get_data() )
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	def _get_data( self ):
		boundaries = [ self.onset_duration ]
		for index, duration in enumerate( self.durations ):
			boundaries.append( boundaries[-1] + duration )
		data = np.array( [ boundaries[ :-1 ], boundaries[ 1: ], self.durations, self.phonemes, self.effects ] ).T
		#print( data.shape )
		return pd.DataFrame( data, columns = [ 'onset', 'offset', 'duration', 'phoneme', 'effect' ] )
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	def _get_phoneme_boundary_times( self, file_path: str, n_phonemes: int ):
		data, samplerate = librosa.load( file_path, 44100 )
		onset, offset = detect_onset_and_offset( file_path )
		duration = len( data )
		start = onset
		end =  offset
		#print( 'num phon: {}'.format(n_phonemes) )
		#print( 'onset: {}'.format( onset ) )
		#print( 'offset: {}'.format( offset ) )
		preds = [ x for x in np.arange( start, end, (end-start)/(n_phonemes) ) ] #+1*round( (end-start)/number_phonemes)
		preds.append( offset )
		#print( 'returning {} boundaries'.format(len(preds)) )
		return preds
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	def _get_uniform_durations_from_audio_file( self, ):
		return
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
#	def export_seg( self, seg_file_path, account_for_effects = False ): Deprecated use to_seg_file instead
#		assert len( self.phonemes ) == len( self.durations ), 'Lengths do not match'
#		out_file = open( seg_file_path, 'w+' )
#		out_file.write( 'name = {}; duration_s = {};\n'.format( '', self.onset_duration ) )
#		for index, phoneme in enumerate( self.phonemes ):
#			if account_for_effects == False or ( account_for_effects == True and self.effect[ index ] != 'elision' ):
#				if self.effect[ index ] != None and self.effect[ index ] != 'elision':
#					output_phone = self.effect[ index ]
#				else:
#					output_phone = phoneme
#				out_file.write( 'name = {}; duration_s = {};\n'.format( output_phone, self.durations[ index ] ) )
#		out_file.write( 'name = {}; duration_s = {};\n'.format( '', self.offset_duration ) )
#		return
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	def get_variants( self, account_for_effects = 'on' ):
		assert len( self._phonemes ) == len( self._boundaries.times ) - 1, 'Lengths do not match'

		phonemes = []
		#boundary_times = []
		for index, phoneme in enumerate( self._phonemes ):
			if self.effect[ index ] != 'elision':
				if self.effect[ index ] != None and self.effect[ index ] != 'elision':
					phonemes.append( self.effect[ index ] )
				else:
					phonemes.append( phoneme.name )
			start = self._boundaries.times[ 0 ]
			end = self._boundaries.times[ -1 ]
			boundary_times = [ x for x in np.arange( start, end, (end-start)/( len(phonemes) ) ) ]
			boundary_times.append( end )

		return [ self, Segment_Sequence( phonemes, Boundaries( boundary_times ) ) ]
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	def import_seg( self, file_path ):
		return
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	def get_phoneme_names( self, ):
		return [ phoneme.name for phoneme in self._phonemes ]
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	def to_seg_file( self, seg_file_path, account_for_effects = False ):
		assert len( self.phonemes ) == len( self.durations), 'Lengths do not match: {}, {}'.format( self.phonemes, self.durations )
		out_file = open( seg_file_path, 'w+' )
		out_file.write( 'name = {}; duration_s = {};\n'.format( '', self.onset_duration ) )
		for index, phoneme in enumerate( self.phonemes ):
			if account_for_effects == False or ( account_for_effects == True and self.effects[ index ] != 'elision' ):
				if self.effects[ index ] != None and self.effects[ index ] != 'elision':
					output_phone = self.effects[ index ]
				else:
					output_phone = phoneme
				out_file.write( 'name = {}; duration_s = {};\n'.format( output_phone, self.durations[ index ] ) )
		if self.offset_duration != None:
			out_file.write( 'name = {}; duration_s = {};'.format( '', self.offset_duration ) )
		return seg_file_path
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	def plot( self, axs = None, **kwargs ):
			figure, axs = get_plot( n_rows = 1, axs = axs )
			data = self._get_data()
			for index, (onset, duration, phoneme) in enumerate( zip( data.loc[ :, 'onset' ], data.loc[ :, 'duration' ], data.loc[ :, 'phoneme' ] ) ):
				axs[0].axvline( onset, **PT.segment_plot_kwargs.get( 'boundaries' ) )
				axs[0].set( ylim = [ 0, 1 ] )
				axs[0].text( onset + 0.5 * duration, 0.5, phoneme, **PT.segment_plot_kwargs.get( 'phonemes' ) )
			axs[0].axvline( data[ 'offset' ].iloc[ -1 ], **PT.segment_plot_kwargs.get( 'boundaries' ) )
			axs[0].tick_params(
				axis = 'y',          # changes apply to the x-axis
				which = 'both',      # both major and minor ticks are affected
				left = False,        # ticks along the left edge are off
				right = False,       # ticks along the right edge are off
				labelleft = False)   # labels along the left edge are off
			plt.xlabel( 'Time [s]' )
			finalize_plot( figure, axs, **kwargs )
			return axs
#####################################################################################################################################################