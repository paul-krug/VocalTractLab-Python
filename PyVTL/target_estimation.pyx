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
#cimport cVocalTractLabApi

import numpy as np
cimport numpy as np
import pandas as pd
#import ctypes
import os
import warnings
import time

from libcpp cimport bool
from libcpp.vector cimport vector  

import logging

logging.basicConfig()
log = logging.getLogger(__name__)
log.setLevel(logging.WARNING)

#import atexit

from cpython.pycapsule cimport *
#from dataclasses import dataclass
import utils_method as um
#from libcpp.vector cimport vector  
#from cpython cimport array
#from libc.stdlib cimport malloc, free
import matplotlib.pyplot as plt
from itertools import cycle

import PyVTL.targets as tg










#---------------------------------------------------------------------------------------------------------------------------------------------------#
#####################################################################################################################################################
#
#	C++ API functions:
#
#####################################################################################################################################################
#---------------------------------------------------------------------------------------------------------------------------------------------------#
cdef extern from "Data.h":
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	cdef struct Sample:
		double time, value
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	cdef struct PitchTarget:
		double slope, offset, tau, duration
#---------------------------------------------------------------------------------------------------------------------------------------------------#
cdef extern from "TargetOptimizerApi.h":
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	cdef struct FitData:
		vector[ PitchTarget ] res_targets,
		vector[ double ] res_boundaries,
		vector[ Sample ] res_trajectory,
		vector[ double ] res_ftmp,
		double res_fmin,
		double res_rmse,
		double res_corr,
		double res_time,
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	FitData estimate_targets(	
							#double *arr_times,
							#double *arr_values,
							#double *arr_boundaries,
							vector[double] input_times,
							vector[double] input_values,
							vector[double] boundaries,
							#int init_bounds,

							double weight_slope,
							double weight_offset,
							double weight_tau,
							double weight_lambda,

							double delta_slope,
							double delta_offset,
							double delta_tau,
							double delta_boundary,

							double mean_slope,
							#double mean_offset,
							double mean_tau,

							int max_iterations,
							int max_cost_evaluations,
							double rho_end,

							bool use_early_stopping,
							double epsilon,
							int patience,
						);
#---------------------------------------------------------------------------------------------------------------------------------------------------#
#####################################################################################################################################################






#---------------------------------------------------------------------------------------------------------------------------------------------------#
class Time_Signal():
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	def __init__( self, input_times, input_values ):
		self.times = input_times
		self.values = input_values
		return
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	@classmethod
	def from_pitch_tier_file( cls, pitch_tier_file_path ):
		df_pitch_tier = pd.read_csv( pitch_tier_file_path , sep = '\t', skiprows= 3 , header = None )
		df_pitch_tier.columns = [ 'time', 'value' ]
		print( df_pitch_tier )
		return cls( df_pitch_tier[ 'time' ].to_numpy(), df_pitch_tier[ 'value' ].to_numpy() )
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	@classmethod
	def from_csv( cls, csv_file_path, **kwargs ):
		df = pd.read_csv( csv_file_path , **kwargs )
		#df.columns = [ 'time', 'value' ]
		return cls( df[ 'time' ], df[ 'value' ] )
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	def get_samples( self, ):
		cdef vector[ Sample ] samples
		cdef Sample sample
		for time, value in zip( self.times, self.values ):
			sample.time = time
			sample.value = value
			samples.push_back( sample )
		return samples
#---------------------------------------------------------------------------------------------------------------------------------------------------#



#---------------------------------------------------------------------------------------------------------------------------------------------------#
class Fit_Result():
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	@um.autoargs()
	def __init__(
		self,
		in_times,
		in_values,
		in_boundaries,
		par_init_bounds,
		par_weight_slope,
		par_weight_offset,
		par_weight_tau,
		par_weight_lambda,
		par_delta_slope,
		par_delta_offset,
		par_delta_tau,
		par_delta_boundary,
		par_mean_slope,
		#par_mean_offset: float
		par_mean_tau,
		par_max_iterations,
		par_max_cost_evaluations,
		par_rho_end,
		par_use_early_stopping,
		par_epsilon,
		par_patience,
		out_targets,
		out_boundaries,
		out_trajectory,
		out_ftmp,
		out_fmin,
		out_rmse,
		out_corr,
		out_time,
		):
		pass
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	def plot( self, ):
		min_in_value = np.min( self.in_values )
		max_in_value = np.max( self.in_values ) 
		plt.scatter( self.in_times, self.in_values, c = 'darkorange' )
		plt.plot( self.out_trajectory[ :, 0 ], self.out_trajectory[ :, 1 ], color = 'navy' )
		for in_boundary, out_boundary in zip( self.in_boundaries, self.out_boundaries):
			plt.axvline( in_boundary, ymin = 0.95, color = 'black' )
			plt.axvline( out_boundary, ymax = 0.95, color = 'lightgray', ls = '--' )
		plt.ylim( [ min_in_value - 0.2 * np.abs( max_in_value - min_in_value ), max_in_value + 0.2 * np.abs( max_in_value - min_in_value ) ] )
		plt.show()
		return
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	def to_csv( self, csv_file_path ):
		return
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	def to_gestural_score( self, ges_file_path, ges_type = 'f0-gestures', unit = 'st' ): #Attention: neutral is always 0
		ges_file = open( 'ges_file_path', 'w' )
		ges_file.write( '<gestural_score>\n' )
		ges_file.write( '\t<gesture_sequence type="{}" unit="{}">\n'.format( ges_type, unit) )
		for target in self.out_targets:
			ges_file.write( '\t\t<gesture value="{}" slope="{}" duration_s="{}" time_constant_s="{}" neutral="0" />\n'.format(
				target.offset, target.slope, target.duration, target.tau ) )
		ges_file.write( '\t</gesture_sequence>\n' )
		ges_file.write( '</gestural_score>\n' )
		return
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	def to_pitch_tier( self, pitch_tier_file_path ):
		pitch_tier_file = open( 'pitch_tier_file_path', 'w' )
		pitch_tier_file.write( 'ooTextFile\n' )
		pitch_tier_file.write( 'PitchTier\n' )
		pitch_tier_file.write( '{} {} {}\n'.format( self.file_onset, self.file_offset, self.out_trajectory.shape[0] ) )
		for sample in self.out_trajectory:
			pitch_tier_file.write( '{}\t{}\n'.format( sample[0], sample[1] ) )
		return
#---------------------------------------------------------------------------------------------------------------------------------------------------#


#---------------------------------------------------------------------------------------------------------------------------------------------------#
class Fit_Window():
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	def __init__( self, times, values, boundaries, queue_position, window_length = 3, hop_length = 1, **kwargs ):
		self.times = []
		self.values = []
		self.boundaries = boundaries[ queue_position : queue_position + window_length + 1 ] # Todo Update: hop_length * queue_pos....
		for time, value in zip( times, values ):
			if time >= self.boundaries[ 0 ] and time <= self.boundaries[ -1 ]:
				self.times.append( time )
				self.values.append( value )
		self.queue_position = queue_position
		self.target_positions = range( queue_position + hop_length, queue_position + hop_length + window_length )
		self.window_length = window_length
		self.hop_lenth = hop_length
		self.fit_result = None
		self.kwargs = kwargs
		return
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	def fit( self, **kwargs ):
		self.fit_result = fit( self.times, self.values, self.boundaries, **self.kwargs )
		return
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	def plot():
		return
#---------------------------------------------------------------------------------------------------------------------------------------------------#

#---------------------------------------------------------------------------------------------------------------------------------------------------#


#---------------------------------------------------------------------------------------------------------------------------------------------------#
class Sequential_Fit():
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	def __init__( self, windows ):
		self.windows = windows
		for window in self.windows:
			window.fit()
		return
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	def plot( self, ):
		figure, axs = plt.subplots( len( self.windows ), figsize = (8, 4/3 *len( self.windows ) ), sharex = True, gridspec_kw = {'hspace': 0} )
		for window in self.windows:
			for boundary in window.boundaries:
				axs[ window.queue_position ].axvline( boundary, color = 'navy' )
			axs[ window.queue_position ].axvspan( 0, window.boundaries[0], alpha = 0.5, color = 'lightgray' )
			axs[ window.queue_position ].axvspan( window.boundaries[ -1 ], self.windows[ -1 ].boundaries[ -1 ], alpha = 0.5, color = 'lightgray' )
			axs[ window.queue_position ].scatter( window.fit_result.in_times, window.fit_result.in_values, color = 'black' )
			if window.queue_position == 0:
				axs[ window.queue_position ].axvspan( window.boundaries[0], window.boundaries[-1], alpha = 0.5, color = 'green' )
			else:
				axs[ window.queue_position ].axvspan( window.boundaries[0], window.boundaries[1], alpha = 0.5, color = 'red' )
				axs[ window.queue_position ].axvspan( window.boundaries[1], window.boundaries[-1], alpha = 0.5, color = 'green' )



		#plt.xlabel( 'Tract state' )
		for ax in axs:
			ax.label_outer()
		figure.align_ylabels( axs[:] )
		plt.show()
		return
#---------------------------------------------------------------------------------------------------------------------------------------------------#

#---------------------------------------------------------------------------------------------------------------------------------------------------#




#---------------------------------------------------------------------------------------------------------------------------------------------------#
#####################################################################################################################################################
#
#	User functions:
#
#####################################################################################################################################################
#---------------------------------------------------------------------------------------------------------------------------------------------------#

def convert_hz_to_st( values ):
	return 12 * ( np.log( values ) / np.log( 2 ) )


#---------------------------------------------------------------------------------------------------------------------------------------------------#
# 		Single core functions
#---------------------------------------------------------------------------------------------------------------------------------------------------#
def fit(	
			times,
			values,
			boundaries = [],
			init_bounds = 4,
			weight_slope = 1.0, 
			weight_offset = 0.5,
			weight_tau = 0.1, 
			weight_lambda = 0.01,
			delta_slope = 0.5,
			delta_offset = 1.0,
			delta_tau = 5,
			delta_boundary = 100.0,
			mean_slope = 0.0,
			#mean_offset = None,
			mean_tau = 20.0,
			max_iterations = 20,
			max_cost_evaluations = 1e6,
			rho_end = 1e-6,
			use_early_stopping = False,
			epsilon = 0.01,
			patience = None,
			):
	if patience == None:
		patience = -1
	if boundaries == [] and init_bounds != 0:
		step =  ( times[ -1 ] - times[ 0 ] ) / ( init_bounds - 1 )
		boundaries = np.array( [ times[ 0 ] + x * step for x in range( 0 , init_bounds ) ] )
		#print( boundaries )
		#print( len(boundaries ) )
	norm_factor_a = np.min( values )
	norm_factor_b = np.max( values ) - np.min( values ) 
	values = ( values - norm_factor_a ) / norm_factor_b
 

	#cdef np.ndarray[ np.float64_t, ndim=1 ] c_input_times = times
	#cdef np.ndarray[ np.float64_t, ndim=1 ] c_input_values = values
	#cdef np.ndarray[ np.float64_t, ndim=1 ] c_boundaries = boundaries
	cdef vector[double] c_input_times = times
	cdef vector[double] c_input_values = values
	cdef vector[double] c_boundaries = boundaries
	#c_boundaries.push_back( 5 ) 
	#print( c_boundaries )
	#print( c_boundaries.size() )
	cdef int c_init_bounds = init_bounds
	cdef double c_weight_slope = weight_slope
	cdef double c_weight_offset = weight_offset
	cdef double c_weight_tau = weight_tau
	cdef double c_weight_lambda = weight_lambda
	cdef double c_delta_slope = delta_slope
	cdef double c_delta_offset = delta_offset
	cdef double c_delta_tau = delta_tau
	cdef double c_delta_boundary = delta_boundary
	cdef double c_mean_slope = mean_slope
	#cdef double c_mean_offset
	cdef double c_mean_tau = mean_tau
	cdef int c_max_iterations = max_iterations
	cdef int c_max_cost_evaluations = max_cost_evaluations
	cdef double c_rho_end = rho_end
	cdef bool c_use_early_stopping = use_early_stopping
	cdef double c_epsilon = epsilon
	cdef int c_patience = patience
	cdef FitData fit_results
	fit_results = estimate_targets(
						#&c_input_times[0],
						#&c_input_values[0],
						#&c_boundaries[0],
						c_input_times,
						c_input_values,
						c_boundaries,
						#c_init_bounds,
						c_weight_slope,
						c_weight_offset,
						c_weight_tau,
						c_weight_lambda,
						c_delta_slope,
						c_delta_offset,
						c_delta_tau,
						c_delta_boundary,
						c_mean_slope,
						#c_mean_offset,
						c_mean_tau,
						c_max_iterations,
						c_max_cost_evaluations,
						c_rho_end,
						c_use_early_stopping,
						c_epsilon,
						c_patience,
					)
	print( 'res targets: ------------------------')
	print( fit_results.res_targets )
	fit_info = dict(
		in_times = np.array( times ),
		in_values = values * norm_factor_b + norm_factor_a,
		in_boundaries = boundaries,
		par_init_bounds = init_bounds,
		par_weight_slope = weight_slope,
		par_weight_offset = weight_offset,
		par_weight_tau = weight_tau,
		par_weight_lambda = weight_lambda,
		par_delta_slope = delta_slope,
		par_delta_offset = delta_offset,
		par_delta_tau = delta_tau,
		par_delta_boundary = delta_boundary,
		par_mean_slope = mean_slope,
		par_mean_tau = mean_tau,
		par_max_iterations = max_iterations,
		par_max_cost_evaluations = max_cost_evaluations,
		par_rho_end = rho_end,
		par_use_early_stopping = use_early_stopping,
		par_epsilon = epsilon,
		par_patience = patience,
		out_targets = [ tg.Target(	fit_results.res_boundaries.at( i ),
									fit_results.res_boundaries.at( i + 1 ),
									target['slope'] * norm_factor_b,
									target['offset'] * norm_factor_b + norm_factor_a, 
									target['tau'] / 1000, # normalization because tau is in [ms] in c++ code 
								) for i, target in enumerate( fit_results.res_targets ) ] ,
		#out_targets = np.array( [ [ target.slope * norm_factor_b + norm_factor_a, 
		#							target.offset * norm_factor_b + norm_factor_a, 
		#							target.tau, 
		#							target.duration 
		#						] for target in fit_results.res_targets ] ),
		out_boundaries = np.array( fit_results.res_boundaries ),
		out_trajectory = np.array( [ [ sample.time, sample.value * norm_factor_b + norm_factor_a ] for sample in fit_results.res_trajectory ] ),
		#out_trajectory = np.array( [ [ sample.time, sample.value ] for sample in fit_results.res_trajectory ] ),
		out_ftmp = np.array( fit_results.res_ftmp ),
		out_fmin = fit_results.res_fmin,
		out_rmse = fit_results.res_rmse,
		out_corr = fit_results.res_corr,
		out_time = fit_results.res_time,
	)
	return Fit_Result( **fit_info )
#---------------------------------------------------------------------------------------------------------------------------------------------------#
def fit_sequentially( times, 
	                  values, 
	                  boundaries = [], 
	                  init_bounds = 9, 
	                  window_length = 3, 
	                  hop_length= 1, 
	                  n_passes = 1, 
	                  **kwargs 
	                  ):
	if boundaries == [] and init_bounds != 0:
		step =  ( times[ -1 ] - times[ 0 ] ) / ( init_bounds - 1 )
		boundaries = np.array( [ times[ 0 ] + x * step for x in range( 0 , init_bounds ) ] )
	if window_length >= len( boundaries ):
		log.warning( 'window_length must be smaller than number of boundaries for sequential fit! Proceed with global fit now.' )
		return fit( times, values, boundaries, **kwargs )
	n_windows = int( ( ( len( boundaries ) - 1 - window_length ) / hop_length ) + 1 )
	#print( 'number of windows: ', n_windows )
	windows = [ Fit_Window( times, values, boundaries, queue_position = x, 
		window_length= window_length, hop_length=hop_length, **kwargs ) for x in range( 0, n_windows ) ]
	#for window in windows:
	#	print( 'window queue pos: {}'.format( window.queue_position ) )
	#	print( 'boundaries: {}'.format( window.boundaries ) )
	#raise ValueError( 'dail ')
	seq_fit = Sequential_Fit( windows )
	seq_fit.plot()
	#for window in windows:
	#	window.fit()
		#print()
		#for target in window.fit_result.out_targets:
		#	print( 'tg')
		#	print( 'onset: {}, offset: {}, slope: {}, offset: {}, tau: {}'.format( target.onset_time, target.offset_time, target.slope, target.offset, target.tau ) )
		#raise ValueError( 'dail ')
	if 'delta_boundary' in kwargs:
		if kwargs[ 'delta_boundary' ] > 0:
			boundaries = get_optimized_boundaries( windows )
	else:
		target_sequence = []
		for window in windows:
			if window.queue_position == 0:
				target_sequence.extend( window.fit_result.out_targets[ : 1 + hop_length ] )
			elif window.queue_position < n_windows - 1:
				target_sequence.append( window.fit_result.out_targets[ hop_length ] )
			elif window.queue_position == n_windows - 1:
				target_sequence.extend( window.fit_result.out_targets[ 1 : ] )
	return
#---------------------------------------------------------------------------------------------------------------------------------------------------#
def get_optimized_boundaries( windows, window_length, hop_length, o ):
	boundary_sequence = []
	n_targets = windows[-1].target_positions[-1] + 1
	target_index = cycle( [ x for x in range(0, window_length) ] )
	for index in range( 0, hop_length +1 ):
		boundary_sequence.append( windows[0].fit_result.out_targets[ target_index[ index ] ].onset_time )
	for index in range( hop_length + 1, n_targets - 1 ):
		tmp_targets = []
		for window in windows:
			if ( index in window.target_positions ) and ( index != window.target_positions[0] ):
				tmp_targets.append( window.fit_result.out_targets[ target_index[ index ] ] )
		onset_time_mean = np.mean( [ x.onset_time for x in tmp_targets ] )
		boundary_sequence.append( onset_time_mean )
	for index in range( n_targets - 1, n_targets ):
		boundary_sequence.append( windows[0].fit_result.out_targets[ target_index[ index ] ].onset_time )
	boundary_sequence.append( windows[0].fit_result.out_targets[ target_index[ index ] ].offset_time )
	return boundary_sequence
#---------------------------------------------------------------------------------------------------------------------------------------------------#
def get_optimized_targets( windows, hop_length ):
	n_windows = len( windows )
	target_sequence = []
	for window in windows:
		if window.queue_position == 0:
			target_sequence.extend( window.fit_result.out_targets[ : 1 + hop_length ] )
		elif window.queue_position < n_windows - 1:
			target_sequence.append( window.fit_result.out_targets[ hop_length ] )
		elif window.queue_position == n_windows - 1:
			target_sequence.extend( window.fit_result.out_targets[ 1 : ] )
	return target_sequence
#---------------------------------------------------------------------------------------------------------------------------------------------------#


	#for window in windows:
	#	for target in window.fit_result.out_targets:
	#		print( 'tg')
	#		print( 'onset: {}, offset: {}, slope: {}, offset: {}, tau: {}'.format( target.onset_time, target.offset_time, target.slope, target.offset, target.tau ) )




	'''

	return




	partial_signals = []
	partial_boundaries = []
	i = 0
	while i + sequence_length < len( boundaries ):
		partial_signal = []
		#for sample in time_signal:
		for time, value in zip( times, values ):
			if time >= boundaries[ i ] and time <= boundaries[ i + sequence_length ]:
				partial_signal.append( [ time, value ] )
		partial_signals.append( partial_signal )
		partial_boundaries.append( [ boundaries[ x ] for x in range( i, i + sequence_length +1 ) ] )
		i += 1
	partial_fits = []
	for signal, boundaries_par in zip( partial_signals, partial_boundaries ):
		signal = np.array( signal )
		partial_fits.append( fit( signal[:,0], signal[:,1], boundaries_par, **kwargs ) )


	tmp_avg = []
	average_factors = [1,2,3,3,3,3,2,1]

	for index, partial_fit in enumerate( partial_fits ):
		print( '-----------fit:-----')
		tgs= []
		#for tgt in range( 0, len(boundaries)-sequence_length ):
		for _ in range( 0, index ):
			tgs.append( tg.Target( 0, 0, 0, 0, 0 ) )
		for target in partial_fit.out_targets:
			print( 'onset: {}, offset: {}, slope: {}, offset: {}, tau: {}'.format( 
				target.onset_time, target.offset_time, target.slope, target.offset, target.tau ) )
			tgs.append( target )
		for _ in range( index + sequence_length, len(boundaries) ):
			tgs.append( tg.Target( 0, 0, 0, 0, 0 ) )
		tmp_avg.append( tgs )
	

	avg_targets = []

	for index in range( 0, len(boundaries)-1 ):
		if index < len(boundaries)-2: 
			avg_target = tg.Target( np.sum( [ x[ index ].onset_time for x in tmp_avg ] ) / average_factors[ index ],
									np.sum( [ x[ index + 1 ].onset_time for x in tmp_avg ] ) / average_factors[ index+1 ],
									np.sum( [ x[ index ].slope for x in tmp_avg ] ) / average_factors[ index ],
									np.sum( [ x[ index ].offset for x in tmp_avg] ) / average_factors[ index ],
									np.sum( [ x[ index ].tau for x in tmp_avg] ) / average_factors[ index ], 
								)
		else:
			avg_target = tg.Target( np.sum( [ x[ index-1 ].offset_time for x in tmp_avg ] ) / average_factors[ index-1 ],
									np.sum( [ x[ index ].offset_time for x in tmp_avg ] ) / average_factors[ index ],
									np.sum( [ x[ index ].slope for x in tmp_avg ] ) / average_factors[ index ],
									np.sum( [ x[ index ].offset for x in tmp_avg] ) / average_factors[ index ],
									np.sum( [ x[ index ].tau for x in tmp_avg] ) / average_factors[ index ], 
								)
		avg_targets.append( avg_target )



	##average_targets = []
	##average_factors = [1,2,3,3,3,3,2,1]
	#avg_targets = []
	#for index in range(0, len( boundaries )-1 ):
	#	tmp_avg = []
	#	for j, partial_fit in enumerate( partial_fits ):
	#		for target in partial_fit.out_targets:
	#			if np.abs( target.onset_time - boundaries[ index ] ) < 0.0001:
	#				print( 'boundary {}: {}'.format( index, boundaries[index] ) )
	#				print( 'fir nr: ', j)
	#				tmp_avg.append( target )
	#	#avg_target = tg.Target( boundaries[ index ], boundaries[ index + 1 ],
	#	#						np.mean( [ x.slope for x in tmp_avg ] ),
	#	#						np.mean( [ x.offset for x in tmp_avg] ),
	#	#						np.mean( [ x.tau for x in tmp_avg] ), 
	#	#					)
	#	#x = tmp_avg[ int( len(tmp_avg) / 2 ) ]
	#	#avg_target = tg.Target( boundaries[ index ], boundaries[ index + 1 ],
	#	#						x.slope,
	#	#						x.offset,
	#	#						x.tau, 
	#	#					)
	#	avg_target = tg.Target( np.mean( [ x.onset_time for x in tmp_avg ] ),
	#							np.mean( [ x.offset_time for x in tmp_avg ] ),
	#							np.mean( [ x.slope for x in tmp_avg ] ),
	#							np.mean( [ x.offset for x in tmp_avg] ),
	#							np.mean( [ x.tau for x in tmp_avg] ), 
	#						)
	#	avg_targets.append( avg_target )

	'''

	avg_targets = target_sequence
	print( avg_targets )
	partial_fits = [ window.fit_result for window in windows ]

	print( avg_targets )
	for target in avg_targets:
		print( 'onset: {}, offset: {}, slope: {}, offset: {}, tau: {}'.format( target.onset_time, target.offset_time, target.slope, target.offset, target.tau ) )
	tam = tg.Target_Approximation_Model()
	cnt = tam.response( avg_targets,  onset_state = values[0] )


	#for partial_fit in partial_fits:
		#plt.plot( partial_fit.out_trajectory[:,0], partial_fit.out_trajectory[:,1] )
	plt.plot( cnt[:,0], cnt[:,1] )
	for tgt in avg_targets:
		plt.axvline( tgt.onset_time )
	plt.axvline( avg_targets[-1].offset_time )
	plt.scatter( times, values, c = 'black' )
	plt.plot( partial_fits[0].out_trajectory[:,0], partial_fits[0].out_trajectory[:,1] )
	plt.plot( partial_fits[-1].out_trajectory[:,0], partial_fits[-1].out_trajectory[:,1] )
	plt.show()

	#return
#---------------------------------------------------------------------------------------------------------------------------------------------------#
#def get_partial_fits()
#####################################################################################################################################################