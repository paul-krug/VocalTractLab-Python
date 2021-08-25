#cimport cVocalTractLabApi

import numpy as np
cimport numpy as np
import pandas as pd
import ctypes
import os
import warnings
import time

from libcpp cimport bool

import logging

logging.basicConfig()
log = logging.getLogger(__name__)
log.setLevel(logging.WARNING)

import atexit

from cpython.pycapsule cimport *
#from cpython cimport array
#from libc.stdlib cimport malloc, free

import PyVTL.tract_sequence as ts
import PyVTL.audio_tools as AT
import PyVTL.function_tools as FT
import librosa
import multiprocessing as mp
import tqdm
import itertools


#---------------------------------------------------------------------------------------------------------------------------------------------------#
#####################################################################################################################################################
#
#	C++ API functions:
#
#####################################################################################################################################################
#---------------------------------------------------------------------------------------------------------------------------------------------------#

#---------------------------------------------------------------------------------------------------------------------------------------------------#
cdef extern from "VocalTractLabApi.h":
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	int vtlCalcTongueRootAutomatically( bool automaticCalculation );
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	int vtlClose();
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	int vtlExportTractSvg(	double *tractParams,
							const char *fileName
							);
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	int vtlGesturalScoreToAudio(	const char* gesFileName,
									const char* wavFileName,
									double* audio,
									int* numSamples,
									int enableConsoleOutput
									);
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	int vtlGesturalScoreToTractSequence(	const char* gesFileName, 
											const char* tractSequenceFileName
											);
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	int vtlGetConstants(	int *audioSamplingRate,
							int *numTubeSections,
							int *numVocalTractParams,
							int *numGlottisParams
							);
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	int vtlGetGesturalScoreDuration(	const char* gesFileName,
										int* numAudioSamples,
										int* numGestureSamples
										);
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	int vtlGetGlottisParamInfo(	char *names,
								double *paramMin,
								double *paramMax,
								double *paramNeutral
								);
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	int vtlGetTractParamInfo(	char *names,
								double *paramMin,
								double *paramMax,
								double *paramNeutral
								);
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	int vtlGetTractParams(	const char *shapeName,
							double *param
							);
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	void vtlGetVersion( char *version );
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	int vtlInitialize( const char *speakerFileName );
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	int vtlInputTractToLimitedTract(	double* inTractParams,
										double* outTractParams
										);
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	int vtlSegmentSequenceToGesturalScore(	const char *segFileName,
											const char *gesFileName
											);
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	int vtlSynthBlock(	double *tractParams,
						double *glottisParams,
						int numFrames,
						int frameStep_samples,
						double *audio,
						int enableConsoleOutput
						);
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	int vtlTractSequenceToAudio(	const char* tractSequenceFileName,
									const char* wavFileName,
									double* audio,
									int* numSamples,
									);
#---------------------------------------------------------------------------------------------------------------------------------------------------#
#####################################################################################################################################################





#---------------------------------------------------------------------------------------------------------------------------------------------------#
#####################################################################################################################################################
#
#	User functions:
#
#####################################################################################################################################################
#---------------------------------------------------------------------------------------------------------------------------------------------------#


#---------------------------------------------------------------------------------------------------------------------------------------------------#
# 		Single core functions
#---------------------------------------------------------------------------------------------------------------------------------------------------#
def automatic_calculation_of_TRX_and_TRY( bool automatic_calculation = True ):
	cdef bool automaticCalculation = automatic_calculation
	value = vtlCalcTongueRootAutomatically( automaticCalculation )
	print( value )
	if value != 0:
		raise ValueError('VTL API function vtlCalcTongueRootAutomatically returned the Errorcode: {}  (See API doc for info.)' )
	print( 'aut calc was set to {}'.format(automatic_calculation) )
	return
#---------------------------------------------------------------------------------------------------------------------------------------------------#
def get_version():
	cdef char version[32]
	vtlGetVersion( version )
	print( 'Compile date of the library: "%s"' % version.decode() )
	#if self.params.verbose == True:
	#	log.info( 'Compile date of the library: "%s"' % version.decode() )
	return version.decode()
#---------------------------------------------------------------------------------------------------------------------------------------------------#
def get_constants():
	cdef int audioSamplingRate = 0
	cdef int numTubeSections = 0
	cdef int numVocalTractParams = 0
	cdef int numGlottisParams = 0
	value = vtlGetConstants( &audioSamplingRate, &numTubeSections, &numVocalTractParams, &numGlottisParams )
	if value != 0:
		raise ValueError('VTL API function vtlGetConstants returned the Errorcode: {}  (See API doc for info.)'.format( value ) )
	constants = {
		'samplerate': int( audioSamplingRate ),
		'n_tube_sections': int( numTubeSections ),
		'n_tract_params': int( numVocalTractParams ),
		'n_glottis_params': int( numGlottisParams ),
		}
	return constants
#---------------------------------------------------------------------------------------------------------------------------------------------------#
def get_param_info( str params ):
	if params not in [ 'tract', 'glottis' ]:
		log.warning( 'Unknown key in "get_param_info". Key must be "tract" or "glottis". Returning "tract" infos now.' )
		params = 'tract'
	if params == 'tract':
		key = 'n_tract_params'
	elif params == 'glottis':
		key = 'n_glottis_params'
	constants = get_constants()
	names = ( ' ' * 10 * constants[ key ] ).encode()
	cdef np.ndarray[ np.float64_t, ndim=1 ] paramMin = np.empty( constants[key], dtype='float64' )
	cdef np.ndarray[ np.float64_t, ndim=1 ] paramMax = np.empty( constants[key], dtype='float64' )
	cdef np.ndarray[ np.float64_t, ndim=1 ] paramNeutral = np.empty( constants[key], dtype='float64' )
	if params == 'tract':
		value = vtlGetTractParamInfo( names, &paramMin[0], &paramMax[0], &paramNeutral[0] )
	elif params == 'glottis':
		value = vtlGetGlottisParamInfo( names, &paramMin[0], &paramMax[0], &paramNeutral[0] )
	if value != 0:
		raise ValueError('VTL API function vtlGetTractParamInfo or vtlGetGlottisParamInfo returned the Errorcode: {}  (See API doc for info.)'.format( value ) )
	df = pd.DataFrame( np.array( [ paramMin, paramMax, paramNeutral ] ).T, columns = [ 'min', 'max', 'neutral' ] )
	df.index = names.decode().replace('\x00','').strip( ' ' ).split( ' ' )
	return df
#---------------------------------------------------------------------------------------------------------------------------------------------------#
def get_tract_params_from_shape( str shape ): #Todo: glottis shapes as well
	shapeName = shape.encode()
	constants = get_constants()
	cdef np.ndarray[ np.float64_t, ndim=2 ] param = np.empty( shape = (1, constants[ 'n_tract_params' ] ),  dtype='float64' )
	value = vtlGetTractParams( shapeName, &param[0, 0] )
	if value != 0:
		raise ValueError('VTL API function vtlGetTractParams returned the Errorcode: {}  (See API doc for info.)'.format( value ) )
	return ts.Supra_Glottal_Sequence( param )
#---------------------------------------------------------------------------------------------------------------------------------------------------#
def load_speaker_file( str speaker_file_path ):
	_close()
	_initialize( speaker_file_path )
	log.info( 'Loaded new speakerfile: {}. Overwriting existing settings with the values from the new speaker file.'.format( speaker_file_path ) )
	return
#---------------------------------------------------------------------------------------------------------------------------------------------------#

#---------------------------------------------------------------------------------------------------------------------------------------------------#


#---------------------------------------------------------------------------------------------------------------------------------------------------#
def get_gestural_score_audio_duration( str ges_file_path, return_samples = True ):
	gesFileName = ges_file_path.encode()
	cdef int numAudioSamples = 0
	numGestureSamples = NULL
	vtlGetGesturalScoreDuration( gesFileName, &numAudioSamples, numGestureSamples )
	n_samples = numAudioSamples
	if return_samples: # returning number of audio samples
		return n_samples
	#else: # returning time in seconds
	#	return n_samples / self.params.samplerate_audio
#---------------------------------------------------------------------------------------------------------------------------------------------------#


#---------------------------------------------------------------------------------------------------------------------------------------------------#
# 		User mp enabled functions
#---------------------------------------------------------------------------------------------------------------------------------------------------#
def gestural_score_to_audio(	ges_file_path_list,
								audio_file_path_list = None,
								save_file: bool = True,
								normalize_audio: int = None,
								sr: int = None,
								return_data: bool = False,
								workers: int = None,
								verbose: bool = False,
							):
	ges_file_path_list, audio_file_path_list = FT.check_if_input_lists_are_valid( [ ges_file_path_list, audio_file_path_list ], 
                                                                                  [ str, ( str, type(None) ) ],
	                                                                            )
	args =  [ [ges_file_path, audio_file_path, save_file, normalize_audio, sr, verbose]
		for ges_file_path, audio_file_path in itertools.zip_longest( ges_file_path_list, audio_file_path_list ) ]
	audio_data_list = _run_multiprocessing( _gestural_score_to_audio, args, return_data, workers )
	return audio_data_list
#---------------------------------------------------------------------------------------------------------------------------------------------------#
def gestural_score_to_tract_sequence(	ges_file_path_list,  
										tract_file_path_list = None, 
										return_data: bool = False,
										workers: int = None,
									):
	ges_file_path_list, tract_file_path_list = FT.check_if_input_lists_are_valid( [ ges_file_path_list, tract_file_path_list ], 
	                                                                              [ str, ( str, type(None) ) ],
	                                                                            )
	args = [ [ges_file_path, tract_file_path, return_data ]
		for ges_file_path, tract_file_path in itertools.zip_longest( ges_file_path_list, tract_file_path_list ) ]
	tract_sequence_list = _run_multiprocessing( _gestural_score_to_tract_sequence, args, return_data, workers )
	return tract_sequence_list
#---------------------------------------------------------------------------------------------------------------------------------------------------#
def segment_sequence_to_gestural_score(	seg_file_path_list,
										ges_file_path_list = None,
										workers: int = None,
									):
	seg_file_path_list, ges_file_path_list = FT.check_if_input_lists_are_valid( [ seg_file_path_list, ges_file_path_list ], 
	                                                                            [ str, ( str, type(None) ) ] 
	                                                                          )
	args = [ [ seg_file_path, ges_file_path ] 
		for seg_file_path, ges_file_path in itertools.zip_longest( seg_file_path_list, ges_file_path_list ) ]
	_run_multiprocessing( _segment_sequence_to_gestural_score, args, False, workers )
	return
#---------------------------------------------------------------------------------------------------------------------------------------------------#
def tract_sequence_to_audio(	tract_sequence_list,
								audio_file_path_list = None,
								save_file: bool = True,
								normalize_audio: int = None,
								sr: int = None,
								return_data: bool = False,
								workers: int = None,
								verbose: bool = False,
							):
	tract_sequence_list, audio_file_path_list = FT.check_if_input_lists_are_valid( [ tract_sequence_list, audio_file_path_list ], 
																		           [ ( str, ts.Tract_Sequence, ts.Target_Sequence ),
	                                                                                 ( str, type(None) ),
	                                                                               ]
	                                                                             )
	args = [ [tract_sequence, audio_file_path, save_file, normalize_audio, sr, verbose]
		for tract_sequence, audio_file_path in itertools.zip_longest( tract_sequence_list, audio_file_path_list ) ]
	audio_data_list = _run_multiprocessing( _tract_sequence_to_audio, args, return_data, workers )
	return audio_data_list
#---------------------------------------------------------------------------------------------------------------------------------------------------#
def tract_sequence_to_limited_tract_sequence(	tract_sequence: ts.Tract_Sequence,
												workers: int = None, 
											):
	tract_param_data = []
	args = [ [state] for state in tract_sequence.tract.to_numpy() ]
	tract_param_data = _run_multiprocessing( _tract_state_to_limited_tract_state, args, True, workers )
	return ts.Supra_Glottal_Sequence( tract_param_data )
#---------------------------------------------------------------------------------------------------------------------------------------------------#
def tract_sequence_to_svg(	tract_sequence_list,
							svg_dir_list = None,
							fps: int = 60,
							save_video = False,
							workers: int = None,
						):
	tract_sequence_list, svg_dir_list = FT.check_if_input_lists_are_valid( [ tract_sequence_list, svg_dir_list ],
	                                                                       [ ( str, ts.Tract_Sequence, ts.Target_Sequence ),
	                                                                         ( str, type(None) ),
	                                                                       ]
	                                                                     )
	args = [ [tract_sequence, svg_dir, fps ]
	for tract_sequence, svg_dir in itertools.zip_longest( tract_sequence_list, svg_dir_list ) ]
	_run_multiprocessing( _tract_sequence_to_svg, args, False, workers )
	return
#---------------------------------------------------------------------------------------------------------------------------------------------------#
#####################################################################################################################################################





#---------------------------------------------------------------------------------------------------------------------------------------------------#
#####################################################################################################################################################
#
#	Internal functions:
#
#####################################################################################################################################################
#---------------------------------------------------------------------------------------------------------------------------------------------------#


#---------------------------------------------------------------------------------------------------------------------------------------------------#
# 		constructor / destructor functions
#---------------------------------------------------------------------------------------------------------------------------------------------------#
def _initialize( str speaker_file_path ):
	speakerFileName = speaker_file_path.encode()
	value = vtlInitialize( speakerFileName )
	if value != 0:
		raise ValueError('VTL API function vtlInitialize returned the Errorcode: {}  (See API doc for info.)'.format( value ) )
	#log.info( 'VTL API initialized.' )
	return
#---------------------------------------------------------------------------------------------------------------------------------------------------#
def _close():
	value = vtlClose()
	if value != 0:
		raise ValueError('VTL API function vtlClose returned the Errorcode: {}  (See API doc for info.)'.format( value ) )
	#log.info( 'VTL API closed.' )
	return
#---------------------------------------------------------------------------------------------------------------------------------------------------#


#---------------------------------------------------------------------------------------------------------------------------------------------------#
# 		multiprocessing worker functions
#---------------------------------------------------------------------------------------------------------------------------------------------------#
#def _export_tract_svg( args ):
#	tract_state, out_file_path = args
#	cdef np.ndarray[ np.float64_t, ndim=1 ] tractParams = tract_state
#	fileName = out_file_path.encode()
#	constants = get_constants()
#	value = vtlExportTractSvg( &tractParams[0], fileName )
#	if value != 0:
#		raise ValueError('VTL API function vtlExportTractSvg returned the Errorcode: {}  (See API doc for info.)'.format( value ) )
#	return
#---------------------------------------------------------------------------------------------------------------------------------------------------#
def _gestural_score_to_audio( args ):
	# Note that returning the number of samples via numSamples is deprecated, use getGesturalScoreAudioDuration instead!
	time_start = time.time()
	ges_file_path, audio_file_path, save_file, normalize_audio, sr, verbose = args
	constants = get_constants()
	if sr == None:
		sr = constants[ 'samplerate' ]
	if not os.path.exists( ges_file_path ):
		warnings.warn( 'the specified gestural score file path does not exist: {}. API call will be skipped.'.format( ges_file_path ) )
		return
	if save_file:
		audio_file_path = FT.make_output_path( audio_file_path, ges_file_path.rsplit( '.' )[0] + '.wav' )
	if ( save_file and normalize_audio != None ) or ( save_file and sr != constants[ 'samplerate' ] ):
		save_file = False
		return_audio = True
	if save_file == False:
		wavFileName = ''.encode()
	else:
		wavFileName = audio_file_path.encode()
	gesFileName = ges_file_path.encode()
	cdef np.ndarray[ np.float64_t, ndim=1 ] audio
	cdef int enableConsoleOutput = 1 if verbose == True else 0
	audio = np.zeros( get_gestural_score_audio_duration( ges_file_path, return_samples = True ), dtype='float64' )
	time_synth_start = time.time()
	cdef int numS = 0
	value = vtlGesturalScoreToAudio( gesFileName, wavFileName, &audio[0], &numS, enableConsoleOutput )
	time_synth_end = time.time()
	print( 'elapsed synthesis time {}'.format( time_synth_end-time_synth_start ) )
	if value != 0:
		raise ValueError('VTL API function vtlGesturalScoreToAudio returned the Errorcode: {}  (See API doc for info.) \
			while processing gestural score file (input): {}, audio file (output): {}'.format(value, ges_file_path, audio_file_path) )
	if sr != constants[ 'samplerate' ]:
		audio = librosa.resample( audio, constants[ 'samplerate' ], sr )
	if normalize_audio != None:
		audio = AT.normalize( audio, normalize_audio )
	if verbose:
		log.info( 'Audio generated from gestural score file: {}'.format( ges_file_path ) )
	if save_file == False and audio_file_path not in ( None, '' ):
		AT.write( audio, audio_file_path, sr )
	time_end = time.time()
	print( 'elapsed total time {}'.format( time_end-time_start ) )
	return audio
#---------------------------------------------------------------------------------------------------------------------------------------------------#
def _gestural_score_to_tract_sequence( args ):
	ges_file_path, tract_file_path, return_sequence = args
	if not os.path.exists( ges_file_path ):
		warnings.warn( 'the specified gestural score file path does not exist: {}. API call will be skipped.'.format( ges_file_path ) )
		return
	gesFileName = ges_file_path.encode()
	tract_file_path = FT.make_output_path( tract_file_path, ges_file_path.rsplit('.')[0] + '.tract' )
	tractSequenceFileName = tract_file_path.encode()
	value = vtlGesturalScoreToTractSequence( gesFileName, tractSequenceFileName )
	if value != 0:
		raise ValueError('VTL API function vtlGesturalScoreToTractSequence returned the Errorcode: {}  (See API doc for info.) \
			while processing gestural score file (input): {}, tract sequence file (output): {}'.format(value, ges_file_path, tract_file_path) )
	log.info( 'Created tractsequence file {} from gestural score file: {}'.format( tract_file_path, ges_file_path ) )
	if return_sequence:
		return ts.Tract_Sequence.from_tract_file( tract_file_path )
	return
#---------------------------------------------------------------------------------------------------------------------------------------------------#
def _segment_sequence_to_gestural_score( args ):
	seg_file_path, ges_file_path = args
	if not os.path.exists( seg_file_path ):
		warnings.warn( 'the specified segment sequence file path does not exist: {}. API call will be skipped.'.format( seg_file_path ) )
		return
	ges_file_path = FT.make_output_path( ges_file_path, seg_file_path.rsplit('.')[0] + '.ges' )
	segFileName = seg_file_path.encode()
	gesFileName = ges_file_path.encode()
	value = vtlSegmentSequenceToGesturalScore( segFileName, gesFileName )
	if value != 0:
		raise ValueError('VTL API function vtlSegmentSequenceToGesturalScore returned the Errorcode: {}  (See API doc for info.) \
			while processing segment sequence file (input): {}, gestural score file (output): {}'.format(value, seg_file_path, ges_file_path) )
	log.info( 'Created gestural score from segment sequence file: {}'.format( seg_file_path ) )
	return
#---------------------------------------------------------------------------------------------------------------------------------------------------#
def _synth_block( args ):
	tract_sequence, state_samples, verbose = args
	constants = get_constants()
	if state_samples == None:
		state_samples = 110
	cdef int numFrames = tract_sequence.length
	cdef np.ndarray[ np.float64_t, ndim=1 ] tractParams = tract_sequence.tract.to_numpy().ravel()
	cdef np.ndarray[ np.float64_t, ndim=1 ] glottisParams = tract_sequence.glottis.to_numpy().ravel()
	cdef int frameStep_samples = state_samples
	cdef np.ndarray[ np.float64_t, ndim=1 ] audio = np.zeros( tract_sequence.length * state_samples, dtype='float64' )
	cdef int enableConsoleOutput = 1 if verbose == True else 0
	value = vtlSynthBlock( &tractParams[0], &glottisParams[0], numFrames, frameStep_samples, &audio[0], enableConsoleOutput )
	if value != 0:
		raise ValueError( 'VTL API function vtlSynthBlock returned the Errorcode: {}  (See API doc for info.)'.format( value ) )
	return audio
#---------------------------------------------------------------------------------------------------------------------------------------------------#
def _tract_sequence_to_audio( args ):
	# Note that returning the number of samples via numSamples is deprecated, use getGesturalScoreAudioDuration instead!
	tract_sequence_data, audio_file_path, save_file, normalize_audio, sr, verbose = args
	if isinstance( tract_sequence_data, str ):
		tract_file_path = tract_sequence_data
		if not os.path.exists( tract_file_path ):
			warnings.warn( 'the specified tract sequence file path does not exist: {}. API call will be skipped.'.format( tract_file_path ) )
			return
		tract_sequence = ts.Tract_Sequence.from_file( tract_file_path )
	elif isinstance( tract_sequence_data, ts.Target_Sequence ):
		target_sequence = tract_sequence_data
		tract_sequence = target_sequence.to_tract_sequence()
	else:
		tract_sequence = tract_sequence_data
	audio = _synth_block( ( tract_sequence, None, verbose ) )
	constants = get_constants()
	if sr == None:
		sr = constants[ 'samplerate' ]
	if save_file:
		audio_file_path = FT.make_output_path( audio_file_path, tract_sequence.name.rsplit( '.' )[0] + '.wav' )
	if sr != constants[ 'samplerate' ]:
		audio = librosa.resample( audio, constants[ 'samplerate' ], sr )
	if normalize_audio != None:
		audio = AT.normalize( audio, normalize_audio )
	if save_file:
		AT.write( audio, audio_file_path, sr )
	log.info( 'Audio generated from tract_sequence: {}'.format( tract_sequence.name ) )
	return audio
#---------------------------------------------------------------------------------------------------------------------------------------------------#
def _tract_state_to_limited_tract_state( args ):
	tract_state = args
	constants = get_constants()
	cdef np.ndarray[ np.float64_t, ndim=1 ] inTractParams = tract_state.ravel()
	cdef np.ndarray[ np.float64_t, ndim=1 ] outTractParams = np.zeros( constants[ 'n_tract_params' ], dtype='float64' )
	vtlInputTractToLimitedTract( &inTractParams[0], &outTractParams[0] )
	return outTractParams
#---------------------------------------------------------------------------------------------------------------------------------------------------#
def _tract_sequence_to_svg( args ):
	tract_sequence, svg_dir, fps = args
	if isinstance( tract_sequence, str ) :
		tract_file_path = tract_sequence
		if not os.path.exists( tract_file_path ) :
			warnings.warn( 'the specified tract sequence file path does not exist: {}. API call will be skipped.'.format( tract_file_path ) )
			return
		tract_sequence = ts.Tract_Sequence.from_tract_file( tract_file_path )
	#elif isinstance( tract_sequence, ts.Target_Sequence ) :
	#	target_sequence = tract_sequence
	#	tract_sequence = target_sequence.to_tract_sequence()
	#else:
	#	#pass
	#	#tract_sequence = tract_sequence_data
	svg_dir = FT.make_output_dir( svg_dir, tract_sequence.name.rsplit('.')[0] + '_svg' )
	constants = get_constants()
	cdef np.ndarray[np.float64_t, ndim = 1] tractParams = np.zeros( constants['n_tract_params'], dtype = 'float64' )
	resampled_index = [ round(index * (44100 / 110) / fps) for index in range( 0, tract_sequence.length ) ]
	resampled_tract_states = [ [ index, tract_state ] for index, tract_state in enumerate( tract_sequence.tract.to_numpy() ) if index in resampled_index ]
	for pair in resampled_tract_states:
		index, tract_state = pair
		tractParams = tract_state.ravel()
		fileName = ( svg_dir + '/tract_{}.svg'.format( index ) ).encode()
		vtlExportTractSvg( &tractParams[0], fileName )
	return







#---------------------------------------------------------------------------------------------------------------------------------------------------#
def _run_multiprocessing( function, args, return_data, workers ):
	if workers == None:
		workers = mp.cpu_count()
	pool = mp.Pool( workers )
	tasks = ( ( function, x ) for x in args)
	data = None
	if return_data:
		data = []
		for x in tqdm.tqdm( pool.imap( _worker, tasks ), total=len( args ) ):
			data.append( x )
	else:
		for x in tqdm.tqdm( pool.imap( _worker, tasks ), total=len( args ) ):
			pass
	pool.close()
	pool.join()
	return data
#---------------------------------------------------------------------------------------------------------------------------------------------------#
def _worker( args ):
	function, arg = args
	return function( arg )
#---------------------------------------------------------------------------------------------------------------------------------------------------#
#####################################################################################################################################################





#---------------------------------------------------------------------------------------------------------------------------------------------------#
#####################################################################################################################################################
#
#	Code to execute at module import:
#
#####################################################################################################################################################
#---------------------------------------------------------------------------------------------------------------------------------------------------#


#---------------------------------------------------------------------------------------------------------------------------------------------------#
atexit.register( _close )
_initialize( os.path.join( os.path.dirname(__file__), 'Speaker/JD2.speaker' ) )
#---------------------------------------------------------------------------------------------------------------------------------------------------#


#####################################################################################################################################################