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
import soundfile
import librosa
import numpy as np
import pandas as pd
import itertools

from scipy.spatial.distance import canberra
from scipy.spatial.distance import cityblock
from scipy.spatial.distance import correlation
from scipy.spatial.distance import euclidean
from scipy.spatial.distance import sqeuclidean

from VocalTractLab.multiprocessing_tools import _run_multiprocessing
from VocalTractLab.function_tools import check_if_input_lists_are_valid
from VocalTractLab.function_tools import make_output_path
from VocalTractLab.function_tools import save
#---------------------------------------------------------------------------------------------------------------------------------------------------#
#####################################################################################################################################################

standard_44kHz_spectrogram_kwargs = dict(
	n_fft = 512, # corresponds to 256/11025 s = 23.2199... ms @ sr = 44100 Hz
	hop_length = 110.25, # corresponds to exactly 2.5 ms @ sr = 44100 Hz
	win_length = 512,
	)

standard_22kHz_spectrogram_kwargs = dict(
	n_fft = 512, # corresponds to 256/11025 s = 23.2199... ms @ sr = 22050 Hz
	hop_length = 55.125, # corresponds to exactly 2.5 ms @ sr = 22050 Hz
	win_length = 512,
	)

standard_16kHz_spectrogram_kwargs = dict(
	n_fft = 256, # corresponds to 2/125 s = 16 ms @ sr = 16000 Hz
	hop_length = 40, # corresponds to exactly 2.5 ms @ sr = 16000 Hz
	win_length = 256,
	)

standard_16kHz_melspectrogram_40_kwargs = dict(
	n_fft = 256, # corresponds to 2/125 s = 16 ms @ sr = 16000 Hz
	hop_length = 40, # corresponds to exactly 2.5 ms @ sr = 16000 Hz
	win_length = 256,
	sr = 16000,
	n_mels = 40,
	fmin = 50,
	fmax = 8000,
	)

standard_16kHz_melspectrogram_80_kwargs = dict(
	n_fft = 256, # corresponds to 2/125 s = 16 ms @ sr = 16000 Hz
	hop_length = 40, # corresponds to exactly 2.5 ms @ sr = 16000 Hz
	win_length = 256,
	sr = 16000,
	n_mels = 80,
	fmin = 50,
	fmax = 8000,
	)

standard_16kHz_mfcc_13_kwargs = dict(
	n_fft = 256, # corresponds to 2/125 s = 16 ms @ sr = 16000 Hz
	hop_length = 40, # corresponds to exactly 2.5 ms @ sr = 16000 Hz
	win_length = 256,
	sr = 16000,
	n_mels = 80,
	fmin = 50,
	fmax = 8000,
	n_mfcc = 13,
	)

standard_16kHz_mfcc_26_kwargs = dict(
	n_fft = 256, # corresponds to 2/125 s = 16 ms @ sr = 16000 Hz
	hop_length = 40, # corresponds to exactly 2.5 ms @ sr = 16000 Hz
	win_length = 256,
	sr = 16000,
	n_mels = 80,
	fmin = 50,
	fmax = 8000,
	n_mfcc = 26,
	)

#---------------------------------------------------------------------------------------------------------------------------------------------------#
def calculate_spectral_distances( reference_list, 
	                              query_list,
	                              #spectral_feature = 'mfcc_13',
	                              distance_metric = sqeuclidean,
								  batch_size = None,
	                              dtw_correction = False,
								  dfw_correction = False,
	                              workers = None,
	                              ):
	reference_list, query_list = check_if_input_lists_are_valid( [ reference_list, query_list ], [ np.ndarray, np.ndarray ] )
	args_list = [ [ [
		reference, query, distance_metric ]
		for query in query_list
		] for reference in reference_list
	]
	spectral_ditstance_list = [ _run_multiprocessing( _calculate_spectral_distance, args, True, workers ) for args in args_list ]
	return spectral_ditstance_list
#---------------------------------------------------------------------------------------------------------------------------------------------------#
def calculate_spectral_features(
	audio_file_path_list,
	spectrogram_file_path_list = None,
	audio_kwargs = dict(
		sr = 16000,
	),
	spectrogram_kwargs = dict(
		n_fft = 256,
		hop_length = 64,
		win_length = 256,
	),
	mel_kwargs = dict(
		sr = 16000,
		n_mels = 40,
		n_fft = 256,
		hop_length = 64,
		win_length = 256,
		fmin = 50,
		fmax = 8000,
	),
	workers = None,
):
	audio_file_path_list, spectrogram_file_path_list = check_if_input_lists_are_valid( 
		[ audio_file_path_list, spectrogram_file_path_list ], [ str, ( str, type(None) ) ] 
	)
	args = [ [ audio_file_path, spectrogram_file_path, audio_kwargs, spectrogram_kwargs, mel_kwargs ]
		for audio_file_path, spectrogram_file_path in itertools.zip_longest( audio_file_path_list, spectrogram_file_path_list ) ]
	_run_multiprocessing( _calculate_spectral_features, args, False, workers )
	return
#---------------------------------------------------------------------------------------------------------------------------------------------------#
def get_f0( audio_file_path_list,
	        csv_file_path_list = None,
	        lower_f0_limit = 50,
	        upper_f0_limit = 250,
	        save_file: bool = True,
	        return_data: bool = False,
	        workers = None,
	        ):
	audio_file_path_list, csv_file_path_list = check_if_input_lists_are_valid( [ audio_file_path_list, csv_file_path_list ], 
	                                                                              [ str, ( str, type(None) ) ],
	                                                                            )
	args = [ [audio_file_path, csv_file_path, lower_f0_limit, upper_f0_limit, save_file ]
		for audio_file_path, csv_file_path in itertools.zip_longest( audio_file_path_list, csv_file_path_list ) ]
	f0_data_list = _run_multiprocessing( _get_f0, args, return_data, workers )
	return f0_data_list
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
	reference, query, distance = args
	tmp = []
	for column in range( 0, query.shape[1] ):
		vec_ref = reference[ :, column ]
		vec_qry = query[ :, column ]
		tmp.append( distance( vec_ref, vec_qry ) )
	return np.mean( tmp )
#---------------------------------------------------------------------------------------------------------------------------------------------------#
def _calculate_spectral_features( args ):
	audio_file_path, spectrogram_file_path, audio_kwargs, spectrogram_kwargs, mel_kwargs  = args
	audio, sr = librosa.load( audio_file_path, sr = audio_kwargs[ 'sr' ] )
	#start_sample = df[ 'frame_start' ] * int( sr / video_options[ 'frame_rate' ] )
	#end_sample = df[ 'frame_end' ] * int( sr / video_options[ 'frame_rate' ] )
	#roi = audio[ start_sample : end_sample ]
	#if len( roi ) > 0:
	#norm_max = np.max( np.abs( roi ) )
	#roi /= ( norm_max + ( norm_max * 0.122018 ) ) # Normalize audio to -1 dB
	spectrogram = np.abs( librosa.stft( y = audio, **spectrogram_kwargs ) )**2
	mel = librosa.feature.melspectrogram( S = spectrogram, **mel_kwargs )
	#if 'mfcc' in features:
	#mfcc = librosa.feature.mfcc( S = librosa.power_to_db( mel ), **mel_settings )
	save( spectrogram, spectrogram_file_path + '.pkl.gzip' )
	save( mel, spectrogram_file_path + '_mel_{}.pkl.gzip'.format( mel_kwargs[ 'n_mels' ] ) )
	return
#---------------------------------------------------------------------------------------------------------------------------------------------------#
def _get_f0( args ):
	try:
		import parselmouth
	except Exception:
		raise ModuleNotFoundError( 'Can not load praat-parselmouth module, which is needed for f0 estimation. Is it installed?' )
	audio_file_path, csv_file_path, lower_f0_limit, upper_f0_limit, save_file = args
	if not os.path.exists( audio_file_path ):
		warnings.warn( 'the specified audio file path does not exist: {}. API call will be skipped.'.format( audio_file_path ) )
		return
	pitch = parselmouth.Sound( audio_file_path ).to_pitch()
	pitch_times  = pitch.xs()
	pitch_values = pitch.selected_array[ 'frequency' ]
	pitch_values[ pitch_values == 0 ] = np.nan
	pitch_values[ pitch_values >= upper_f0_limit ] = np.nan
	pitch_values[ pitch_values <= lower_f0_limit ] = np.nan
	#data = np.array( [ pitch_times, pitch_values ] ).T
	df = pd.DataFrame( np.array( [ pitch_times, pitch_values ] ).T, columns = [ 'time', 'f0' ] ).dropna()
	if save_file:
		csv_file_path = make_output_path( csv_file_path, audio_file_path.rsplit( '.' )[0] + '.csv' )
		df.to_csv( csv_file_path )
	return df
#---------------------------------------------------------------------------------------------------------------------------------------------------#