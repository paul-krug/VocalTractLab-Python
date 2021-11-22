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

import warnings
import numpy as np
import pandas as pd
from scipy.special import binom
from scipy.special import factorial
import matplotlib.pyplot as plt
from itertools import zip_longest
from itertools import chain
from PyVTL import plotting_tools as PT
from PyVTL.plotting_tools import finalize_plot
from PyVTL.plotting_tools import get_plot
from PyVTL.plotting_tools import get_plot_limits
from PyVTL.plotting_tools import get_valid_tiers
from PyVTL import function_tools as FT
from PyVTL import tract_sequence as TS
from PyVTL.tract_sequence import Sub_Glottal_Sequence, Supra_Glottal_Sequence, Tract_Sequence
from PyVTL.audio_tools import get_f0
import PyVTL.VocalTractLabApi as vtl
from collections import Counter
from PyVTL.target_estimation import fit

def get_name_of( variable ):
	return [ k for k, v in locals().items() if v is variable ][ 0 ]


#####################################################################################################################################################
#---------------------------------------------------------------------------------------------------------------------------------------------------#
class Target():
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	"""PyVTL articulatory target""" 
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	def __init__( self,
		          onset_time: float = 0.0, 
		          duration: float = 1.0, 
		          slope: float = 0.0, 
		          offset: float = 1.0, 
		          time_constant = 0.015,
		          ):
		self.onset_time = onset_time
		self.offset_time = onset_time + duration
		self.duration = duration
		self.slope = slope
		self.offset = offset
		self.time_constant = time_constant
		self.m = slope
		self.b = offset
		self.tau = time_constant
		return
	# def __init__( self, onset_time, offset_time, slope, offset, time_constant ):
	# 	self.onset_time = onset_time
	# 	self.offset_time = offset_time
	# 	self.duration = offset_time - onset_time
	# 	self.slope = slope
	# 	self.offset = offset
	# 	self.time_constant = time_constant
	# 	self.tau = time_constant
	# 	return
#---------------------------------------------------------------------------------------------------------------------------------------------------#
#####################################################################################################################################################


'''
#####################################################################################################################################################
#---------------------------------------------------------------------------------------------------------------------------------------------------#
class Synchronized_Target():
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	"""PyVTL articulatory target""" 
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	def __init__( self, onset_time, offset_time, slope, offsets, time_constant ):
		self.targets = []
		for offset in offsets:
			self.targets.append( Target( onset_time, offset_time, slope, offset, time_constant ) )
		return
#---------------------------------------------------------------------------------------------------------------------------------------------------#
#####################################################################################################################################################
'''


#####################################################################################################################################################
#---------------------------------------------------------------------------------------------------------------------------------------------------#
class Target_Sequence():
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	"""PyVTL articulatory target""" 
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	def __init__( self,
		          onset_time = 0.0,
		          durations: list = [],
		          slopes: list = [],
		          offsets: list = [],
		          time_constants: list = [], 
		          targets: list = None,
		          name = '',
		          ):
		if targets != None:
			self.targets = targets
		else:
			self.targets = []
			function_argument_names = [ 'duration', 'slope', 'offset', 'time_constant' ]
			for idx_target, args in enumerate( zip_longest( durations, slopes, offsets, time_constants ) ):
				kwargs = {}
				kwargs[ 'onset_time' ] = onset_time
				for idx_member, arg in enumerate( args ):
					if arg != None:
						kwargs[ function_argument_names[ idx_member ] ] = arg
				self.targets.append( Target( **kwargs ) )
				onset_time += self.targets[-1].duration
		self.name = name
		return
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	@classmethod
	def from_audio_file( cls, audio_file_path, **kwargs ):
		data = get_f0( audio_file_path )
		fit_result = fit( data[ 'time' ].to_numpy(), data[ 'f0' ].to_numpy(), **kwargs )
		return cls( targets = fit_result.out_targets, name = 'f0' )
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	@classmethod
	def from_data( cls, data, name, **kwargs ):
		fit_result = fit( data[ :, 0 ], data[ :, 1 ], **kwargs )
		return cls( targets = fit_result.out_targets, name = name )
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	def __str__( self, ):
		columns = [ 'onset_time', 'duration', 'slope', 'offset', 'time_constant' ]
		return str( pd.DataFrame( [ [ tar.onset_time, tar.duration, tar.m, tar.b, tar.tau ] for tar in self.targets ], columns = columns ) )
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	def plot( self, plot_contour = True, plot_targets = True, ax = None, plot_kwargs = PT.state_plot_kwargs, **kwargs ): #, time = 'seconds'
		figure, ax = get_plot( n_rows = 1, axs = ax )
		if plot_contour:
			tam = Target_Approximation_Model()
			contour = tam.response( self.targets )
			try:
				contour_kwargs = plot_kwargs.get( self.name )
			except Exception:
				warnings.warn( 'The parameter: {} does not exist in the plot_kwargs dict, doing a standard plot now.'.format( parameter ) )
				contour_kwargs = dict( color = 'navy' )
			x = contour[ :, 0 ]
			#if time == 'samples':
			#	constants = vtl.get_constants()
			#	x *= constants[ 'samplerate_internal' ]
			ax[ 0 ].plot( x, contour[ :, 1 ], **contour_kwargs )
			ax[ 0 ].set( ylim = get_plot_limits( contour[ :, 1 ], 0.3 ) )
		if plot_targets:
			ax[ 0 ].axvline( self.targets[0].onset_time, color = 'black' )
			y_data = []
			for tar in self.targets:
				#x = tar.offset_time
				#if time == 'samples':
				#	constants = vtl.get_constants()
				#	x *= constants[ 'samplerate_internal' ]
				ax[ 0 ].axvline( tar.offset_time, color = 'black' )
				x = [ tar.onset_time, tar.offset_time ]
				y = [ tar.slope * (tar.onset_time-tar.onset_time) + tar.offset, tar.slope * (tar.offset_time-tar.onset_time) + tar.offset ]
				ax[ 0 ].plot( x, y, color = 'black', linestyle='--' )
				y_data.append( y )
			if not plot_contour:
				ax[ 0 ].set( ylim = get_plot_limits( y_data, 0.3 ) )
		ax[ 0 ].set( xlabel = 'Time [s]', ylabel = self.name )
		#ax[ 0 ].label_outer()
		finalize_plot( figure, ax, **kwargs )
		return ax
#---------------------------------------------------------------------------------------------------------------------------------------------------#
#####################################################################################################################################################


#####################################################################################################################################################
#---------------------------------------------------------------------------------------------------------------------------------------------------#
class Target_Score():
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	"""PyVTL articulatory target""" 
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	def plot( self, parameters = None, plot_contour = True, plot_targets = True, axs = None, **kwargs ):
		parameters = get_valid_tiers( parameters, self.names )
		figure, axs = get_plot( n_rows = len( parameters ), axs = axs )
		index = 0
		for target_sequence in self.target_sequences:
			if target_sequence.name in parameters:
				target_sequence.plot( plot_contour = plot_contour, plot_targets= plot_targets, ax = axs[ index ], show=False )
				index += 1
		finalize_plot( figure, axs, **kwargs )
		return
#---------------------------------------------------------------------------------------------------------------------------------------------------#
#####################################################################################################################################################





#####################################################################################################################################################
#---------------------------------------------------------------------------------------------------------------------------------------------------#
class Synchronous_Target_Score( Target_Score ):
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	"""PyVTL articulatory target""" 
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	def __init__( self,
		          durations: list,
		          names: list,
		          onset_time: float = 0.0,
		          slope_score: list = [],
		          offset_score: list = [],
		          time_constant_score: list = [],
		          ):
		onset_time = onset_time
		self.target_sequences = []
		for args in zip_longest( names, slope_score, offset_score, time_constant_score ):
			#print( 'args: {}'.format(args) )
			args_corrected = []
			for x in args:
				try:
					if x == None:
						args_corrected.append( [] )
					else:
						args_corrected.append( x )
				except Exception:
					args_corrected.append( x )
			args = args_corrected
			#print( args )
			#args = [ [] if (not isinstance( x, list) ) and (x == None) else x for x in args]
			name, slopes, offsets, time_constants = args
			#print( 'na,e:{}, slopes:{}, off:{}, tico:{}'.format( name, slopes, offsets, time_constants) )
			self.target_sequences.append( Target_Sequence( onset_time, durations, slopes, offsets, time_constants, name = name ) )
		self.names = [ target_sequence.name for target_sequence in self.target_sequences ]
		return
#---------------------------------------------------------------------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------------------------------------------------------------------#
#####################################################################################################################################################




#####################################################################################################################################################
#---------------------------------------------------------------------------------------------------------------------------------------------------#
class Supra_Glottal_Motor_Score( Target_Score ):
#---------------------------------------------------------------------------------------------------------------------------------------------------#
#
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	def __init__( self, target_scores
		          ):
		self.target_scores = target_scores
		self.target_sequences = [ target_sequence for target_score in self.target_scores for target_sequence in target_score.target_sequences ]
		self.tiers = [ target_sequence.name for target_sequence in self.target_sequences ]
		self.param_info = vtl.get_param_info( 'tract' )
		#print( self.param_info.index )
		if Counter( self.tiers ) != Counter( self.param_info.index ):
			raise ValueError( 'Tiers in Motor Score are {}, but should be {}.'.format( self.tiers, self.param_info.index ) )
		self.target_sequences.sort( key = lambda i: list( self.param_info.index ).index( i.name ) )
		self.names = [ target_sequence.name for target_sequence in self.target_sequences ]
		return
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	@classmethod
	def from_supra_glottal_sequence( cls,
		                             supra_glottal_sequence,
		                             synchronous = [ [ 'VO' ], [ 'other' ] ],
		                             durations = [ [0.5,0.2], [0.3,0.4] ],
		                             slopes = None,
		                             time_constants = None,
		                             onset_time = 0.0 ):
		if not isinstance( supra_glottal_sequence, ( TS.Supra_Glottal_Sequence, TS.Tract_Sequence ) ):
			raise ValueError( 'Passed argument supra_glottal_sequence is of type {}, but should be {} or {}.'.format( type(supra_glottal_sequence),
			TS.Supra_Glottal_Sequence, TS.Tract_Sequence ) )
		target_scores = []
		for synchronous_tiers, durations in zip( synchronous, durations ):
			if 'other' in synchronous_tiers:
				synchronous_tiers = [ x for x in supra_glottal_sequence.tract.columns if x not in sum( synchronous, [] ) ]
			#print( synchronous_tiers )
			offsets = supra_glottal_sequence.tract[ synchronous_tiers ].to_numpy().T
			#print(offsets)
			target_scores.append( Synchronous_Target_Score( durations, synchronous_tiers, onset_time, offset_score = offsets ) )
		return cls( target_scores )
#---------------------------------------------------------------------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------------------------------------------------------------------#
#####################################################################################################################################################



#####################################################################################################################################################
#---------------------------------------------------------------------------------------------------------------------------------------------------#
class Sub_Glottal_Motor_Score( Target_Score ):
#---------------------------------------------------------------------------------------------------------------------------------------------------#
#
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	def __init__( self, target_scores
		          ):
		self.target_scores = target_scores
		self.target_sequences = [ target_sequence for target_score in self.target_scores for target_sequence in target_score.target_sequences ]
		self.tiers = [ target_sequence.name for target_sequence in self.target_sequences ]
		self.param_info = vtl.get_param_info( 'glottis' )
		#print( self.param_info.index )
		if Counter( self.tiers ) != Counter( self.param_info.index ):
			raise ValueError( 'Tiers in Motor Score are {}, but should be {}.'.format( self.tiers, self.param_info.index ) )
		self.target_sequences.sort( key = lambda i: list( self.param_info.index ).index( i.name ) )
		self.names = [ target_sequence.name for target_sequence in self.target_sequences ]
		return
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	@classmethod
	def from_sub_glottal_sequence( cls,
		                           sub_glottal_sequence,
		                           synchronous = [ [ 'F0' ], [ 'PR' ], [ 'other' ] ],
		                           durations = [ [0.5,0.2], [0.1,0.7], [0.3,0.4] ],
		                           slopes = None,
		                           time_constants = None,
		                           onset_time = 0.0 ):
		if not isinstance( sub_glottal_sequence, ( TS.Sub_Glottal_Sequence, TS.Tract_Sequence ) ):
			raise ValueError( 'Passed argument sub_glottal_sequence is of type {}, but should be {} or {}.'.format( type(sub_glottal_sequence),
			TS.Sub_Glottal_Sequence, TS.Tract_Sequence ) )
		target_scores = []
		for synchronous_tiers, durations in zip( synchronous, durations ):
			if 'other' in synchronous_tiers:
				synchronous_tiers = [ x for x in sub_glottal_sequence.glottis.columns if x not in sum( synchronous, [] ) ]
			#print( synchronous_tiers )
			#print( sub_glottal_sequence )
			offsets = sub_glottal_sequence.glottis[ synchronous_tiers ].to_numpy().T
			#print(offsets)
			target_scores.append( Synchronous_Target_Score( durations, synchronous_tiers, onset_time, offset_score = offsets ) )
		return cls( target_scores )
#---------------------------------------------------------------------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------------------------------------------------------------------#
#####################################################################################################################################################



#####################################################################################################################################################
#---------------------------------------------------------------------------------------------------------------------------------------------------#
class Motor_Score( Target_Score ):
#---------------------------------------------------------------------------------------------------------------------------------------------------#
#
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	def __init__( self,
		          supra_glottal_motor_score,
		          sub_glottal_motor_score,
		          ):
		self.supra_glottal_target_sequences = supra_glottal_motor_score.target_sequences
		self.sub_glottal_target_sequences = sub_glottal_motor_score.target_sequences
		self.target_sequences = [ target_sequence for target_sequence in chain( self.supra_glottal_target_sequences, self.sub_glottal_target_sequences ) ]
		self.names = [ target_sequence.name for target_sequence in self.target_sequences ]
		return
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	@classmethod
	def from_tract_sequence( cls, tract_sequence ):
		supra_glottal_motor_score = Supra_Glottal_Motor_Score.from_supra_glottal_sequence( tract_sequence.to_supra_glottal_sequence() )
		sub_glottal_motor_score = Sub_Glottal_Motor_Score.from_sub_glottal_sequence( tract_sequence.to_sub_glottal_sequence() )
		return cls( supra_glottal_motor_score, sub_glottal_motor_score )
#---------------------------------------------------------------------------------------------------------------------------------------------------#

#---------------------------------------------------------------------------------------------------------------------------------------------------#

#---------------------------------------------------------------------------------------------------------------------------------------------------#
#####################################################################################################################################################



#####################################################################################################################################################
#---------------------------------------------------------------------------------------------------------------------------------------------------#
class Target_Approximation_Model():
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	"""PyVTL articulatory target""" 
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	def __init__( self, order = 5 ):
		self.FILTERORDER = order
		return
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	def response( self, target_sequence, discretization_rate = 44100 / 110, onset_state = None  ) -> np.array:
		trajectory = []
		start = target_sequence[ 0 ].onset_time
		end   = target_sequence[ -1 ].offset_time
		duration = end - start
		n_samples = duration * discretization_rate
		sample_times = np.arange( start, end, duration / n_samples )
		if onset_state == None:
			onset_state = target_sequence[0].offset
		current_state = [ onset_state ]
		for _ in range( 1, self.FILTERORDER ):
			current_state.append( 0.0 )

		b_begin = target_sequence[ 0 ].onset_time
		b_end = b_begin

		sample_index = 0

		for target in target_sequence:
			b_begin = b_end
			b_end = b_begin + target.duration
			c = self.calculate_coefficients( target, current_state )
			while( sample_times[ sample_index ] <= b_end ):
				constant = 0.0
				t = sample_times[ sample_index ] - b_begin
				for n in range( 0, self.FILTERORDER ):
					constant += c[ n ] * ( t**n )
				time = sample_times[ sample_index ]
				value= constant * np.exp( - (1/target.time_constant) * t ) + target.slope * t + target.offset
				trajectory.append( [ time, value ] )
				sample_index += 1
				if sample_index >= len( sample_times ):
					break
			current_state = self.calculate_state( current_state, b_end, b_begin, target );

		return np.array( trajectory )
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	def calculate_coefficients( self, target, current_state ):
		coefficients = [ 0 for _ in current_state ]
		assert len( coefficients ) == self.FILTERORDER, 'Sometimes size does matter bro...'
		coefficients[ 0 ] = current_state[ 0 ] - target.offset
		for n in range( 1, self.FILTERORDER ):
			acc = 0
			for i in range( 0, n ):
				acc += ( coefficients[ i ] * ( (-1 / target.time_constant)**(n - i) ) * binom( n, i ) * factorial( i ) )
			if n == 1:
				acc += target.slope # adaption for linear targets; minus changes in following term!
			coefficients[ n ] = ( current_state[ n ] - acc ) / factorial( n )
		return coefficients
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	def calculate_state( self, state, time, start_time, target ):
		t = time - start_time
		state_update = [ 0 for _ in range( 0, self.FILTERORDER ) ]
		c = self.calculate_coefficients( target, state)
		for n in range( 0, self.FILTERORDER ):
			acc = 0
			for i in range( 0, self.FILTERORDER ):
				q = 0
				for k in range( 0, np.min( [ self.FILTERORDER - i, n + 1 ] ) ):
					q += ( ( (-1 / target.time_constant)**(n - k) ) * binom(n, k) * c[i + k] * factorial(k + i) / factorial(i) )
				acc += ( (t**i) * q );
			state_update[ n ] = acc * np.exp( -( 1 / target.time_constant) * t)
		# correction for linear targets
		if (self.FILTERORDER > 1):
			state_update[ 0 ] += (target.offset + target.slope * t)
		if (self.FILTERORDER > 2):
			state_update[ 1 ] += target.slope
		return state_update
#---------------------------------------------------------------------------------------------------------------------------------------------------#
#	def binomial( self, n, k ):
#		result = 1
#		tmp = k
#		if (tmp > n - tmp)
#			tmp = n - tmp;
#		for i in range( 0, tmp ):
#			result *= (n - i)
#			result /= (i + 1)
#	return result
##---------------------------------------------------------------------------------------------------------------------------------------------------#
#	def factorial( self, ):

#####################################################################################################################################################



#target_1 = Target( onset_time = 0.1, offset_time= 0.3, slope=0.0, offset=0, time_constant=0.005  )
#target_2 = Target( onset_time = 0.3, offset_time= 0.7, slope=0.0, offset=3000, time_constant=0.005  )
#target_3 = Target( onset_time = 0.7, offset_time= 0.9, slope=0.0, offset=0, time_constant=0.015  )
#target_sequence = [target_1, target_2, target_3]
#tam = Target_Approximation_Model()
#trajectory = tam.response( target_sequence )
#plt.plot( trajectory[:,0], trajectory[:,1] )
#plt.show()