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
import soundfile
import numpy as np
import pandas as pd
import parselmouth

from scipy.spatial.distance import canberra
from scipy.spatial.distance import cityblock
from scipy.spatial.distance import correlation
from scipy.spatial.distance import euclidean
from scipy.spatial.distance import sqeuclidean
#---------------------------------------------------------------------------------------------------------------------------------------------------#
#####################################################################################################################################################


#---------------------------------------------------------------------------------------------------------------------------------------------------#
def calculate_spectral_distances( reference_list, 
	                              query_list,
	                              spectral_feature = 'mfcc_13',
	                              distance_metric = correlation,
								  batch_size = None,
	                              dtw_correction = False,
								  dfw_correction = False,
	                              workers = None,
	                              ):
	reference_list, query_list = FT.check_if_input_lists_are_valid( [ ges_file_path_list, audio_file_path_list ], [ str, ( str, type(None) ) ] )
	args =  [ [reference, query, spectral_feature, distance_metric ]
		for reference, query in itertools.zip_longest( reference_list, query_list ) ]
	spectral_ditstance_list = _run_multiprocessing( _calculate_spectral_distance, args, True, workers )
	return spectral_ditstance_list
#---------------------------------------------------------------------------------------------------------------------------------------------------#



#---------------------------------------------------------------------------------------------------------------------------------------------------#
def normalize( audio, normalization = -1 ): #normalisation in dB
	norm_factor = 10**( -1 * normalization * 0.05 ) -1
	norm_max = np.max( np.abs( audio ),axis=0)
	audio /= ( norm_max + ( norm_max * norm_factor ) )
	return audio
#---------------------------------------------------------------------------------------------------------------------------------------------------#
def write( audio, audio_file_path, sr ):
	soundfile.write( audio_file_path, audio, sr )
	return
#---------------------------------------------------------------------------------------------------------------------------------------------------#

#---------------------------------------------------------------------------------------------------------------------------------------------------#
def _calculate_spectral_distance( args ):
	return spectral_distance
#---------------------------------------------------------------------------------------------------------------------------------------------------#
def get_f0( audio_file_path, lower_limit = 50, upper_limit = 250 ):
	pitch = parselmouth.Sound( audio_file_path ).to_pitch()
	pitch_times  = pitch.xs()
	pitch_values = pitch.selected_array[ 'frequency' ]
	pitch_values[ pitch_values == 0 ] = np.nan
	pitch_values[ pitch_values >= upper_limit ] = np.nan
	pitch_values[ pitch_values <= lower_limit ] = np.nan
	#data = np.array( [ pitch_times, pitch_values ] ).T
	return pd.DataFrame( np.array( [ pitch_times, pitch_values ] ).T, columns = [ 'time', 'f0' ] ).dropna()
#---------------------------------------------------------------------------------------------------------------------------------------------------#