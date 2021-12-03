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
import VocalTractLab.VocalTractLabApi as vtl
from VocalTractLab.tract_sequence import Supra_Glottal_Sequence
from VocalTractLab.targets import Target
from VocalTractLab.targets import Target_Sequence
from VocalTractLab.multiprocessing_tools import _run_multiprocessing
#import matplotlib.pyplot as plt
import random
from random import getrandbits
from random import uniform
from random import randint
#from random import normal
import numpy as np
import bisect
#---------------------------------------------------------------------------------------------------------------------------------------------------#
#####################################################################################################################################################



#---------------------------------------------------------------------------------------------------------------------------------------------------#
#####################################################################################################################################################
#
#	Generation of random target sequences:
#
#####################################################################################################################################################
#---------------------------------------------------------------------------------------------------------------------------------------------------#
def generate_target_sequences(
	n_samples,
	onset_time_range = [ 0.0, 1.0 ],
	onset_state_range = [ 0, 1 ],
	duration_range = [ 0.1, 0.4 ],
	slope_range = [ -0.75, 0.75 ],
	offset_range = [ 0, 1 ],
	tau_range = [ 0.005, 0.025 ],
	n_targets_range = [1, 20],
	balance_slope = True,
	workers = None,
	):
	args = [ [ onset_time_range, onset_state_range, duration_range, slope_range, offset_range, tau_range, str(x), n_targets_range, balance_slope ]
		for x in range( 0, n_samples) ]
	target_sequence_list = _run_multiprocessing( _generate_target_sequence, args, True, workers )
	return target_sequence_list
#---------------------------------------------------------------------------------------------------------------------------------------------------#
def generate_target_contours(
	target_sequence_list,
	sr: float = 100,
	augmentation_kwargs = dict(
		resample_range = 100, # can be tuple or single float
		masking_range = [ 0.4, 0.6 ],
		random_shift_x = False,
		random_shift_y = True,
		random_masking = True,
		fixed_length = False,
		),
	workers: int = None,
	):
	args = [ [ target_sequence, sr, augmentation_kwargs ] for target_sequence in target_sequence_list ]
	contour_list = _run_multiprocessing( _generate_target_contours, args, True, workers )
	return contour_list
#---------------------------------------------------------------------------------------------------------------------------------------------------#
def generate_supra_glottal_states(
	n_samples,
	parameter_kwargs = None,
	):
	tract_param_info = vtl.get_param_info( 'tract' )
	for parameter in parameter_kwargs:
		if parameter not in tract_param_info.index:
			raise ValueError( 'Specified parameter: {} is not a supra-glottal parameter.'.format( parameter ) )
	supra_glottal_states = []
	for parameter in tract_param_info.index:
		try:
			_min = parameter_kwargs[ parameter ][0]
			_max = parameter_kwargs[ parameter ][1]
			print('success in dict try')
		except Exception:
			_min = float( tract_param_info.loc[ parameter, 'min' ] )
			_max = float( tract_param_info.loc[ parameter, 'max' ] )
		supra_glottal_states.append( np.random.uniform( _min, _max, n_samples ) )
	supra_glottal_sequence = Supra_Glottal_Sequence( np.array( supra_glottal_states ).T )
	try:
		for parameter, value in parameter_kwargs.items():
			supra_glottal_sequence.tract[ parameter ] = value
	except Exception:
		pass
	return supra_glottal_sequence
#---------------------------------------------------------------------------------------------------------------------------------------------------#



#---------------------------------------------------------------------------------------------------------------------------------------------------#
# 		single core functions
#---------------------------------------------------------------------------------------------------------------------------------------------------#
def _generate_target_sequence( args ):
	onset_time_range, onset_state_range, duration_range, slope_range, offset_range, tau_range, name, n_targets_range, balance_slope = args
	n_targets = random.randint( n_targets_range[0], n_targets_range[1] )
	onset_time = random.uniform( onset_time_range[0], onset_time_range[1] )
	if onset_state_range != None:
		onset_state = random.uniform( onset_state_range[0], onset_state_range[1] )
	else:
		onset_state = None
	#print( onset_state )
	durations = np.random.uniform( duration_range[0], duration_range[1], n_targets )
	slopes = [ 0 if ( balance_slope and getrandbits(1) ) else np.random.uniform( slope_range[0], slope_range[1] ) for _ in range(0, n_targets) ]
	offsets = np.random.uniform( offset_range[0], offset_range[1], n_targets )
	time_constants = np.random.uniform( tau_range[0], tau_range[1], n_targets )
	target_sequence = Target_Sequence( onset_time = onset_time,
	                                   durations = durations,
	                                   slopes = slopes,
	                                   offsets = offsets,
	                                   time_constants = time_constants,
	                                   name = name,
	                                   onset_state = onset_state,
	                                   )
	return target_sequence
#---------------------------------------------------------------------------------------------------------------------------------------------------#
def _generate_target_contours( args ):
	target_sequence, sr, augmentation_kwargs = args
	if augmentation_kwargs == None:
		return target_sequence.get_contour( sr )
	else:
		return _generate_augmented_target_contours( target_sequence, **augmentation_kwargs )
#---------------------------------------------------------------------------------------------------------------------------------------------------#
def _generate_augmented_target_contours(
	target_sequence,
	resample_range = 100, # can be tuple or single float
	masking_range = [ 0.4, 0.6 ],
	random_shift_x = False,
	random_shift_y = True,
	random_masking = True,
	fixed_length = False,
	):
	try:
		sr = random.uniform( resample_range[0], resample_range[1] )
	except Exception:
		sr = resample_range
	contour = target_sequence.get_contour( sr )
	if random_masking:
		for target in target_sequence.targets:
			if getrandbits(1):
				mask_percent = uniform( masking_range[0], masking_range[1] )
				mask_duration = target.duration * mask_percent
				start = uniform( target.onset_time, target.onset_time + (target.duration - mask_duration) )
				mask_start_idx = bisect.bisect_left( contour[:,0], start )
				mask_end_idx = bisect.bisect_right( contour[:,0], start + mask_duration )
				contour_x = list( contour[:,0] )
				contour_y = list( contour[:,1] )
				if fixed_length:
					contour_y[ mask_start_idx : mask_end_idx ] = [ -2 for _ in range(mask_start_idx , mask_end_idx) ]
				else:
					del contour_x[ mask_start_idx : mask_end_idx ]
					del contour_y[ mask_start_idx : mask_end_idx ]
				contour = np.array( [contour_x, contour_y] ).T
	if random_shift_x or random_shift_y:
		for index, (x, y) in enumerate( zip( contour[ :, 0 ], contour[ :, 1 ] ) ):
			if (not fixed_length) and (random_shift_x):
				contour[ index, 0 ] = x + np.random.normal( 0, 0.005 )
			if random_shift_y:
				contour[ index, 1 ] = y + np.random.normal( 0, 0.01 )
		contour[contour[:, 0].argsort()]
	return contour
#---------------------------------------------------------------------------------------------------------------------------------------------------#