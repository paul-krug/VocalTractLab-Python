#####################################################################################################################################################
#---------------------------------------------------------------------------------------------------------------------------------------------------#
#	- This file is a part of the Python module PyVTL, see https://github.com/TUD-STKS/PyVTL
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
# Requirements:
#	- python 3 (tested with version 3.7)
#	- numpy    (tested with version 1.19.5)
#	- pandas   (tested with version 1.2.1)
#	- scipy    (tested with version 1.6.0)
#---------------------------------------------------------------------------------------------------------------------------------------------------#
#####################################################################################################################################################
#---------------------------------------------------------------------------------------------------------------------------------------------------#
#__all__ = ['a', 'b', 'c']
#__version__ = '0.1'
#__author__ = 'Paul Konstantin Krug'
#---------------------------------------------------------------------------------------------------------------------------------------------------#
#####################################################################################################################################################
#---------------------------------------------------------------------------------------------------------------------------------------------------#
# Load essential packages:
import os, sys, ctypes
import multiprocessing as mp
import warnings
import itertools #import groupby
import PyVTL.core
import PyVTL.function_tools as ft
import tqdm

import pandas as pd
import numpy as np
from scipy.io import wavfile
import matplotlib.pyplot as plt
#---------------------------------------------------------------------------------------------------------------------------------------------------#
#####################################################################################################################################################


'''
#####################################################################################################################################################
#---------------------------------------------------------------------------------------------------------------------------------------------------#
# Define the relative path to the API file:
rel_path_to_vtl = './PyVTL/API/VocalTractLabApi'
# Define the relative path to the speaker file:
rel_path_to_spk = './PyVTL/Speaker/'
# Load the VocalTractLab binary 'VocalTractLabApi.dll' if you use Windows or 'VocalTractLabApi.so' if you use Linux:
try:
	if sys.platform == 'win32':
		rel_path_to_vtl += '.dll'
	else:
		rel_path_to_vtl += '.so'
	VTL = ctypes.CDLL( rel_path_to_vtl )
except:
	print( 'Could not load the VocalTractLab API, does the path "{}" exist?'.format( rel_path_to_vtl ) )
#---------------------------------------------------------------------------------------------------------------------------------------------------#
#####################################################################################################################################################
'''
api = PyVTL.core.VTL_API()
#print( 'sheesh' )
#---------------------------------------------------------------------------------------------------------------------------------------------------#
class VTL():
	def __init__( self, ):
		pass
		#global api
		#self.api = PyVTL.core.VTL_API()
#---------------------------------------------------------------------------------------------------------------------------------------------------#
#####################################################################################################################################################
#
#	User functions:
#
#####################################################################################################################################################
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	def gestural_score_to_audio(	self, 
									ges_file_path_list: list,
									audio_file_path: list = '',
									return_audio: bool = True,
									return_n_samples: bool = False,
									verbose: bool = False
								):
		ges_file_path_list, audio_file_path_list = ft.check_if_input_lists_are_valid( [ges_file_path_list, audio_file_path_list], [str, str] )
		args =  [ [ges_file_path, audio_file_path, return_audio, return_n_samples, verbose]
			for ges_file_path, audio_file_path in itertools.zip_longest( ges_file_path_list, audio_file_path_list ) ]
		audio_data_list = self._run_multiprocessing( '_gestural_score_to_audio', args )
		return audio_data_list
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	def gestural_score_to_tract_sequence(	self,
											ges_file_path_list: list,  
											tract_file_path_list: list = '', 
											return_sequence: bool = False 
										):
		ges_file_path_list, tract_file_path_list = ft.check_if_input_lists_are_valid( [ges_file_path_list, tract_file_path_list], [str, str] )
		args = [ [ges_file_path, tract_file_path, return_sequence]
			for ges_file_path, tract_file_path in itertools.zip_longest( ges_file_path_list, tract_file_path_list ) ]
		df_tract_sequence_list = self._run_multiprocessing( '_gestural_score_to_tract_sequence', args )
		return df_tract_sequence_list
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	def segment_sequence_to_gestural_score(	self,
											seg_file_path_list,
											ges_file_path_list,
										):
	seg_file_path_list, ges_file_path_list = ft.check_if_input_lists_are_valid( [seg_file_path_list, ges_file_path_list], [str, str] )
	args = [ [] 
		for seg_file_path, ges_file_path in itertools.zip_longest( seg_file_path_list, ges_file_path_list) ]
	self._run_multiprocessing( '_segment_sequence_to_gestural_score', args )
	return
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	def tract_sequence_to_audio(	self,
									tract_file_path_list: list,
									audio_file_path_list: list = '',
									return_audio = True,
									return_n_samples = False,
									verbose = False
								):
		tract_file_path_list, audio_file_path_list = ft.check_if_input_lists_are_valid( [tract_file_path_list, audio_file_path_list], [str, str] )
		args = [ [tract_file_path, audio_file_path, return_audio, return_n_samples] 
			for tract_file_path, audio_file_path in itertools.zip_longest( tract_file_path_list, audio_file_path_list ) ]
		audio_data_list = self._run_multiprocessing( '_tract_sequence_to_audio', args )
		return audio_data_list
#---------------------------------------------------------------------------------------------------------------------------------------------------#




#---------------------------------------------------------------------------------------------------------------------------------------------------#
#####################################################################################################################################################
#
#	Internal functions:
#
#####################################################################################################################################################
#---------------------------------------------------------------------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------------------------------------------------------------------#
# 		input argument related functions
#---------------------------------------------------------------------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------------------------------------------------------------------#
'''
	def _check_if_input_lists_are_valid( self, input_lists, instances_list ):
		valid_lists = []
		for input_list, instance in zip( input_lists, instances_list ):
			valid_lists.append( self._check_if_list_is_valid( input_list, instance ) )
		return valid_lists
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	def _check_lengths_of_input_lists( self, input_lists ):
		list_lengths = [ len( input_list ) for input_list in input_lists ]
		if not self._check_if_all_elements_are_equal( list_lengths ):
			warnings.warn( 'input list do not have the same lengths, shorter lists will be padded with "None".' )
			max_length = max( list_lengths )
			#print( max_length )
			for input_list in input_lists:
				while len( input_list ) < max_length:
					input_list.append( None )
			#print( input_list )
		return input_lists
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	def _check_if_object_is_iterable( self, query ):
		try:
			iter( query )
		except TypeError:
			return False
		return True
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	def _check_if_list_is_valid( self, input_list, instance ):
		is_iterable = self._check_if_object_is_iterable( input_list )
		if isinstance( input_list, str ) or is_iterable == False:
			warnings.warn( 'input is either not iterable or a single string. The input gets turned into a list now.' )
			input_list = [ input_list ]
		if input_list and all( isinstance( x, instance ) for x in input_list ):
			return input_list
		else:
			raise TypeError( 'a list containing a non-{} type object was passed, but list of {} was expected.'.format( instance, instance ) )
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	def _check_if_all_elements_are_equal( self, iterable ):
		g = itertools.groupby(iterable)
		return next(g, True) and not next(g, False)
'''
#---------------------------------------------------------------------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------------------------------------------------------------------#
# 		tract sequence related functions
#---------------------------------------------------------------------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	def _df_to_tract_seq( self, tract_file_path, df_GLP, df_VTP, glottis_model = 'Geometric glottis' ):
		f= open( tract_file_path, "w+" )
		f.write("""# The first two lines (below the comment lines) indicate the name of the vocal fold model and the number of states.\n""")
		f.write("""# The following lines contain the control parameters of the vocal folds and the vocal tract (states)\n""")
		f.write("""# in steps of 110 audio samples (corresponding to about 2.5 ms for the sampling rate of 44100 Hz).\n""")
		f.write("""# For every step, there is one line with the vocal fold parameters followed by\n""")
		f.write("""# one line with the vocal tract parameters.\n""")
		f.write("""# \n""")
		f.write("""{}\n""".format( glottis_model ) )
		f.write("""{}\n""".format(len(df_GLP)))
		for index, row in enumerate(df_GLP.index):
			for index2, column in enumerate(df_GLP.columns):
				f.write('{} '.format(df_GLP.iloc[index,index2]))
			f.write('\n')
			for index2, column in enumerate(df_VTP.columns):
				f.write('{} '.format(df_VTP.iloc[index,index2]))
			f.write('\n')
		f.close()
		if self.params.verbose:
			print('Tract Sequence saved as: "{}"'.format( tract_file_path ))
		return
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	def _get_tract_seq_len( self, tract_seq_path ):
		with open( tract_seq_path ) as file:
			for index, line in enumerate( file ):
				if index == 7:
					tract_seq_len = int( line.strip() )
					break
		if self.params.verbose:
			print( 'Entries in Tract Sequence file: {}'.format( tract_seq_len ) )
		return tract_seq_len
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	def _read_tract_seq_GLP(self, index):
		if (index > 7) and (index % 2 == 0):
			return False
		else:
			return True
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	def _read_tract_seq_VTP(self, index):
		if (index > 7) and ((index-1) % 2 == 0):
			return False
		else:
			return True
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	def _tract_seq_to_df( self, tractFilePath ):
		# Skip rows from based on condition
		df_GLP = pd.read_csv( tractFilePath, delim_whitespace = True, skiprows= lambda x: self._read_tract_seq_GLP(x) , header = None )
		df_VTP = pd.read_csv( tractFilePath, delim_whitespace = True, skiprows= lambda x: self._read_tract_seq_VTP(x) , header = None )
		#if self.params.verbose:
		#	print( 'Tract sequence opened: {}'.format( tractFilePath ) )
		return df_GLP, df_VTP
#---------------------------------------------------------------------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------------------------------------------------------------------#
# 		tube sequence related functions
#---------------------------------------------------------------------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	def _df_to_tube_seq( self, tube_file_path, df_GLP, df_param, df_area, df_length, df_art, glottis_model = 'Geometric glottis' ):
		f= open( tube_file_path, "w+" )
		f.write("""# The first two lines (below the comment lines) indicate the name of the vocal fold model and the number of states.\n""")
		f.write("""# The following lines contain a sequence of states of the vocal folds and the tube geometry\n""")
		f.write("""# in steps of 110 audio samples (corresponding to about 2.5 ms for the sampling rate of 44100 Hz).\n""")
		f.write("""# Each state is represented in terms of five lines:\n""")
		f.write("""# Line 1: glottis_param_0 glottis_param_1 ...\n""")
		f.write("""# Line 2: incisor_position_in_cm, velo_pharyngeal_opening_in_cm^2, tongue_tip_side_elevation[-1...1]\n""")
		f.write("""# Line 3: area0 area1 area2 area3 ... (Areas of the tube sections in cm^2 from glottis to mouth)\n""")
		f.write("""# Line 4: length0 length1 length2 length3 ... (Lengths of the tube sections in cm from glottis to mouth)\n""")
		f.write("""# Line 5: artic0 artic1 artic2 artic3 ... (Articulators of the tube sections between glottis and lips : 
			1 = tongue; 2 = lower incisors; 3 = lower lip; 4 = other)\n""")
		f.write("""# \n""")
		f.write("""{}\n""".format( glottis_model ) )
		f.write("""{}\n""".format( len(df_GLP) ) )
		for index, row in enumerate(df_GLP.index):
			for index2, column in enumerate(df_GLP.columns):
				f.write('{} '.format(df_GLP.iloc[index,index2]))
			f.write('\n')
			for index2, column in enumerate(df_param.columns):
				f.write('{} '.format(df_param.iloc[index,index2]))
			f.write('\n')
			for index2, column in enumerate(df_area.columns):
				f.write('{} '.format(df_area.iloc[index,index2]))
			f.write('\n')
			for index2, column in enumerate(df_length.columns):
				f.write('{} '.format(df_length.iloc[index,index2]))
			f.write('\n')
			for index2, column in enumerate(df_art.columns):
				f.write('{} '.format(df_art.iloc[index,index2]))
			f.write('\n')
		f.close()
		if self.params.verbose:
			log.info( 'Tube sequence file was saved at: "{}"'.format( tube_file_path ))
		return
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	def _get_tube_seq_len( self, tube_seq_path ):
		with open( tube_seq_path ) as file:
			for index, line in enumerate( file ):
				if index == 11:
					tube_seq_len = int( line.strip() )
					break
		if self.params.verbose:
			print( 'Entries in Tract Sequence file: {}'.format( tube_seq_len ) )
		return tube_seq_len
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	def _tube_seq_to_df( self, tube_file_path ):
		# Skip rows from based on condition
		df_GLP = pd.read_csv( tube_file_path, delim_whitespace = True, skiprows= lambda x: self._read_tube_seq_GLP(x) , header = None )
		df_tube_param  = pd.read_csv( tube_file_path, delim_whitespace = True, skiprows= lambda x: self._read_tube_seq_param(x) , header = None )
		df_tube_area   = pd.read_csv( tube_file_path, delim_whitespace = True, skiprows= lambda x: self._read_tube_seq_area(x) , header = None )
		df_tube_length = pd.read_csv( tube_file_path, delim_whitespace = True, skiprows= lambda x: self._read_tube_seq_length(x) , header = None )
		df_tube_art    = pd.read_csv( tube_file_path, delim_whitespace = True, skiprows= lambda x: self._read_tube_seq_art(x) , header = None )
		if self.params.verbose:
			print( 'Tube sequence opened: {}'.format( tube_file_path ) )
		return df_GLP, df_tube_param, df_tube_area, df_tube_length, df_tube_art
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	def _read_tube_seq_length(self, index):
		if (index > 11) and (index % 5 == 0):
			return False
		else:
			return True
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	def _read_tube_seq_art(self, index):
		if (index > 11) and ((index-1) % 5 == 0):
			return False
		else:
			return True
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	def _read_tube_seq_GLP(self, index):
		if (index > 11) and ((index-2) % 5 == 0):
			return False
		else:
			return True
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	def _read_tube_seq_param(self, index):
		if (index > 11) and ((index-3) % 5 == 0):
			return False
		else:
			return True
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	def _read_tube_seq_area(self, index):
		if (index > 11) and ((index-4) % 5 == 0):
			return False
		else:
			return True
#---------------------------------------------------------------------------------------------------------------------------------------------------#


#---------------------------------------------------------------------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	def _run_multiprocessing( self, function, args ):
		pool = mp.Pool( api.params.workers )
		tasks = ( ( function, x ) for x in args)
		#pbar = tqdm(total=len(tasks))
		data = []
		for x in tqdm.tqdm( pool.imap( worker, tasks ), total=len( args ) ):
			data.append( x )
		#data = pool.map( worker, ( ( function, x ) for x in args) )
		pool.close()
		pool.join()
		return data
	def prepare_call(self, name, args):  # creates a 'remote call' package for each argument
		for arg in args:
			yield [self.__class__.__name__, self.__dict__, name, arg]
#---------------------------------------------------------------------------------------------------------------------------------------------------#
#	def worker( self, args ):
#		self.api._gestural_score_to_tract_sequence( args )
#		return
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	def _check_if_input_lists_are_valid( self, input_lists, instances_list ):
		valid_lists = []
		for input_list, instance in zip( input_lists, instances_list ):
			valid_lists.append( self._check_if_list_is_valid( input_list, instance ) )
		return valid_lists
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	def _check_lengths_of_input_lists( self, input_lists ):
		list_lengths = [ len( input_list ) for input_list in input_lists ]
		if not self._check_if_all_elements_are_equal( list_lengths ):
			warnings.warn( 'input list do not have the same lengths, shorter lists will be padded with "None".' )
			max_length = max( list_lengths )
			#print( max_length )
			for input_list in input_lists:
				while len( input_list ) < max_length:
					input_list.append( None )
			#print( input_list )
		return input_lists
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	def _check_if_object_is_iterable( self, query ):
		try:
			iter( query )
		except TypeError:
			return False
		return True
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	def _check_if_list_is_valid( self, input_list, instance ):
		is_iterable = self._check_if_object_is_iterable( input_list )
		if isinstance( input_list, str ) or is_iterable == False:
			warnings.warn( 'input is either not iterable or a single string. The input gets turned into a list now.' )
			input_list = [ input_list ]
		if input_list and all( isinstance( x, instance ) for x in input_list ):
			return input_list
		else:
			raise TypeError( 'a list containing a non-{} type object was passed, but list of {} was expected.'.format( instance, instance ) )
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	def _check_if_all_elements_are_equal( self, iterable ):
		g = itertools.groupby(iterable)
		return next(g, True) and not next(g, False)
#---------------------------------------------------------------------------------------------------------------------------------------------------#
def worker( args ):
	global api
	function, arg = args
	return getattr( api, function )( arg )
	#function( args )
	#function, arg= args
	#return getattr( api, function ) 


def parallel_call(params):  # a helper for calling 'remote' instances
	cls = getattr(sys.modules[__name__], params[0])  # get our class type
	instance = cls.__new__(cls)  # create a new instance without invoking __init__
	instance.__dict__ = params[1]  # apply the passed state to the new instance
	method = getattr(instance, params[2])  # get the requested method
	args = params[3] if isinstance(params[3], (list, tuple)) else [params[3]]
	return method(*args)  # expand arguments, call our method and return the result


'''
#####################################################################################################################################################
#---------------------------------------------------------------------------------------------------------------------------------------------------#
class VTL():
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	"""A Python wrapper for VocalTractLab""" 
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	def __init__( self, *args ):
		self.params = vtl_params( *args )
		self.API = self.load_API()
		self.get_version()
		self.initialize()
		return
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	def __del__( self ):
		self.close()
		return
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	def load_API( self ):
		rel_path_to_vtl = os.path.join( os.path.dirname(__file__), self.params.API_dir + self.params.API_name )
		# Load the VocalTractLab binary 'VocalTractLabApi.dll' if you use Windows or 'VocalTractLabApi.so' if you use Linux:
		try:
			if sys.platform == 'win32':
				rel_path_to_vtl += '.dll'
			else:
				rel_path_to_vtl += '.so'
			API = ctypes.CDLL( rel_path_to_vtl )
		except:
			print( 'Could not load the VocalTractLab API, does the path "{}" exist?'.format( rel_path_to_vtl ) )
		return API
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	def load_speaker_file( self, speaker_file_name ):
		self.close()
		self.params.set_speaker_file( speaker_file_name )
		self.initialize()
		return
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	def automatic_calculation_of_TRX_and_TRY( self, automatic_calculation = True ):
		automaticCalculation = ctypes.c_bool( automatic_calculation )
		self.API.vtlCalcTongueRootAutomatically( automatic_calculation )
		return
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	def initialize( self ):
		speaker_file_path = ctypes.c_char_p( self.params.speaker_file_name.encode() )
		print( self.params.speaker_file_name )
		failure = self.API.vtlInitialize( speaker_file_path )
		if failure != 0:
			raise ValueError('Error in vtlInitialize! Errorcode: %i' % failure)
		else:
			if self.params.verbose:
				print('VTL successfully initialized.')
		return
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	def close( self ):
		self.API.vtlClose()
		if self.params.verbose:
			print('VTL successfully closed.')
		return
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	def get_version( self ):
		version = ctypes.c_char_p( ( ' ' * 32 ).encode() )
		self.API.vtlGetVersion( version )
		if self.params.verbose == True:
			print('Compile date of the library: "%s"' % version.value.decode() )
		return
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	def get_constants( self ):
		audioSamplingRate = ctypes.c_int(0)
		numTubeSections = ctypes.c_int(0)
		numVocalTractParams = ctypes.c_int(0)
		numGlottisParams = ctypes.c_int(0)
		self.API.vtlGetConstants( ctypes.byref( audioSamplingRate ), 
			                      ctypes.byref( numTubeSections ),
			                      ctypes.byref( numVocalTractParams ),
			                      ctypes.byref( numGlottisParams ) 
			                     )
		constants = {
		'samplerate': audioSamplingRate.value,
		'n_tube_sections': numTubeSections.value,
		'n_tract_params': numVocalTractParams.value,
		'n_glottis_params': numGlottisParams.value,
		}
		return constants
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	def get_param_info( self, params: str ):
		if params not in [ 'tract', 'glottis' ]:
			print( 'Unknown key in "get_param_info". Key must be "tract" or "glottis". Returning "tract" infos now.' )
			params = 'tract'
		if params == 'tract':
			key = 'n_tract_params'
		elif params == 'glottis':
			key = 'n_glottis_params'
		constants = self.get_constants()
		names = ctypes.c_char_p( ( ' ' * 10 * constants[ key ] ).encode() )
		paramMin        = ( ctypes.c_double * constants[ key ] )()
		paramMax        = ( ctypes.c_double * constants[ key ] )()
		paramNeutral    = ( ctypes.c_double * constants[ key ] )()
		if params == 'tract':
			self.API.vtlGetTractParamInfo( names, ctypes.byref( paramMin ), ctypes.byref( paramMax ), ctypes.byref( paramNeutral ) )
		elif params == 'glottis':
			self.API.vtlGetGlottisParamInfo( names, ctypes.byref( paramMin ), ctypes.byref( paramMax ), ctypes.byref( paramNeutral ) )
		df = pd.DataFrame( np.array( [ paramMin, paramMax, paramNeutral ] ).T, columns = [ 'min', 'max', 'neutral' ] )
		df.index = names.value.decode().strip( ' ' ).split( ' ' )
		#print( df )
		return df
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	def get_tract_params_from_shape( self, shape: str ): #Todo: glottis shapes as well
		shapeName = ctypes.c_char_p( shape.encode() )
		constants = self.get_constants()
		param = ( ctypes.c_double * constants[ 'n_tract_params' ] )()
		self.API.vtlGetTractParams( shapeName, ctypes.byref( param ) )
		return np.array( param )
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	def synth_block( self, tract_params, glottis_params, verbose = False, state_samples = None ):
		constants = self.get_constants()
		#print(constants[ 'n_tract_params' ])
		#print(constants[ 'n_glottis_params' ])
		if state_samples == None:
			state_samples = self.params.state_samples
		if len( tract_params ) != len( glottis_params ):
			print( 'TODO: Warning: Length of tract_params and glottis_params do not match. Will modify glottis_params to match.')
			# Todo: Match length
		numFrames = ctypes.c_int( len( tract_params ) )
		#print( numFrames.value )
		tractParams = (ctypes.c_double * ( numFrames.value * constants[ 'n_tract_params' ] ))()
		#print( len(tractParams) )
		tractParams[:] = tract_params.T.ravel('F')
		glottisParams = (ctypes.c_double * ( numFrames.value * constants[ 'n_glottis_params' ] ))()
		glottisParams[:] = glottis_params.T.ravel('F')
		#print( 'shape: {}'.format( np.array(tractParams).shape ) )
		#stop
		frameStep_samples = ctypes.c_int( state_samples )
		#print( frameStep_samples.value )
		audio = (ctypes.c_double * int( len( tract_params ) / self.params.samplerate_internal * self.params.samplerate_audio ) )()
		enableConsoleOutput = ctypes.c_int(1) if verbose == True else ctypes.c_int(0)
		return_val = self.API.vtlSynthBlock( tractParams, glottisParams, numFrames, frameStep_samples, audio, enableConsoleOutput )
		#print( return_val )
		return np.array( audio )
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	def export_tract_svg( self, tract_params, file_path ):
		constants = self.get_constants()
		for index, state in enumerate( tract_params ):
			tractParams = ( ctypes.c_double * ( constants[ 'n_tract_params' ] ) )()
			tractParams[:] = state.T.ravel('F')
			fileName = ctypes.c_char_p( ( file_path + '_{}.svg'.format( index ) ).encode() )
			self.API.vtlExportTractSvg( ctypes.byref( tractParams ), fileName )
		return
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	def tract_params_to_tube_data( self, tract_params ):
		constants = self.get_constants()
		tube_data = []
		for state in tract_params:
			tractParams = ( ctypes.c_double * ( constants[ 'n_tract_params' ] ) )()
			tractParams[:] = state.T.ravel('F')
			tubeLength_cm = ( ctypes.c_double * ( constants[ 'n_tube_sections' ] ) )()
			tubeArea_cm2  = ( ctypes.c_double * ( constants[ 'n_tube_sections' ] ) )()
			tubeArticulator  = ( ctypes.c_int * ( constants[ 'n_tube_sections' ] ) )()
			incisorPos_cm = ctypes.c_double( 0 )
			tongueTipSideElevation = ctypes.c_double( 0 )
			velumOpening_cm2 = ctypes.c_double( 0 )
			self.API.vtlTractToTube( ctypes.byref( tractParams ), 
			                                   ctypes.byref( tubeLength_cm ), 
											   ctypes.byref( tubeArea_cm2 ),
											   ctypes.byref( tubeArticulator ), 
											   ctypes.byref( incisorPos_cm ),
											   ctypes.byref( tongueTipSideElevation),
											   ctypes.byref( velumOpening_cm2 )
											 )
			tube_data.append( [ np.array( tubeLength_cm ), 
			np.array( tubeArea_cm2 ), 
			np.array( tubeArticulator ), 
			incisorPos_cm.value, 
			tongueTipSideElevation.value, 
			velumOpening_cm2.value ] )
		df = pd.DataFrame( tube_data, columns = [ 'tube_length_cm', 'tube_area_cm2', 'tube_articulator', 'incisor_pos_cm', 'tongue_tip_side_elevation', 'velum_opening_cm2' ])
		return df
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	def get_transfer_function( self, tract_params, n_spectrum_samples: int = 2048): # TODO, raw magnitude spec to freq spec, + formant extraction
		class TransferFunctionOptions( ctypes.Structure ):
			_fields_ = [
			("spectrumType", ctypes.c_int),
			("radiationType", ctypes.c_int),
		    ("boundaryLayer", ctypes.c_bool),
		    ("heatConduction", ctypes.c_bool),
		    ("softWalls", ctypes.c_bool),
		    ("hagenResistance", ctypes.c_bool),
		    ("innerLengthCorrections", ctypes.c_bool),
		    ("lumpedElements", ctypes.c_bool),
		    ("paranasalSinuses", ctypes.c_bool),
		    ("piriformFossa", ctypes.c_bool),
			("staticPressureDrops", ctypes.c_bool)
			]
		opts = TransferFunctionOptions()
		self.API.vtlGetDefaultTransferFunctionOptions( ctypes.byref( opts ) )
		constants = self.get_constants()
		numSpectrumSamples = ctypes.c_int( n_spectrum_samples )
		transfer_function_data = []
		for state in tract_params:
			#print( state )
			#print( state.shape )
			tractParams = ( ctypes.c_double * ( constants[ 'n_tract_params' ] ) )()
			tractParams[:] = state.T.ravel('F')
			magnitude = ( ctypes.c_double * ( n_spectrum_samples ) )()
			phase_rad = ( ctypes.c_double * ( n_spectrum_samples ) )()
			self.API.vtlGetTransferFunction( ctypes.byref( tractParams ), numSpectrumSamples, ctypes.byref( opts ), ctypes.byref( magnitude ), ctypes.byref( phase_rad ) )
			transfer_function_data.append( [ np.array( magnitude ), np.array( phase_rad ) ] )
		df = pd.DataFrame( transfer_function_data, columns = [ 'magnitude', 'phase_rad' ] )
		return df
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	def get_transfer_function_( self, tract_params, n_spectrum_samples: int = 2048): # TODO, raw magnitude spec to freq spec, + formant extraction
		constants = self.get_constants()
		numSpectrumSamples = ctypes.c_int( n_spectrum_samples )
		transfer_function_data = []
		for state in tract_params:
			#print( state )
			#print( state.shape )
			tractParams = ( ctypes.c_double * ( constants[ 'n_tract_params' ] ) )()
			tractParams[:] = state.T.ravel('F')
			magnitude = ( ctypes.c_double * ( n_spectrum_samples ) )()
			phase_rad = ( ctypes.c_double * ( n_spectrum_samples ) )()
			self.API.vtlGetTransferFunction( ctypes.byref( tractParams ), numSpectrumSamples, ctypes.byref( magnitude ), ctypes.byref( phase_rad ) )
			transfer_function_data.append( [ np.array( magnitude ), np.array( phase_rad ) ] )
		df = pd.DataFrame( transfer_function_data, columns = [ 'magnitude', 'phase_rad' ] )
		return df
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	def synthesis_reset( self ):
		self.API.vtlSynthesisReset()
		return
#---------------------------------------------------------------------------------------------------------------------------------------------------#
#	def vtlSynthesisAddTube(int numNewSamples, double *audio,
#  double *tubeLength_cm, double *tubeArea_cm2, int *tubeArticulator,
#  double incisorPos_cm, double velumOpening_cm2, double tongueTipSideElevation,
#  double *newGlottisParams):
#		return
#	def synthesis_add_tube_state( self, tube_length_state, tube_area_state, tube_articulator_state, incisor_position,
#	                              velum_opening, tongue_tip_side_elevation, glottis_state, n_new_samples ): 
#		constants = self.get_constants()
#		numNewSamples = ctypes.c_int( n_new_samples )
#		audio = (ctypes.c_double * int( n_new_samples ) )()
#		tractParams = (ctypes.c_double * ( constants[ 'n_tract_params' ] ))()
#		tractParams[:] = tract_state.T.ravel('F')
#		glottisParams = (ctypes.c_double * ( constants[ 'n_glottis_params' ] ))()
#		glottisParams[:] = glottis_state.T.ravel('F')
#		self.API.vtlSynthesisAddTract( numNewSamples, ctypes.byref( audio ), ctypes.byref( tractParams ), ctypes.byref( glottisParams ) )
#		if n_new_samples > 0:
#			return np.array( audio )
#
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	def synthesis_add_tract_state( self, tract_state, glottis_state, n_new_samples ): # Inefficient in Python, for larger synthesis use synth_block
		constants = self.get_constants()
		numNewSamples = ctypes.c_int( n_new_samples )
		audio = (ctypes.c_double * int( n_new_samples ) )()
		tractParams = (ctypes.c_double * ( constants[ 'n_tract_params' ] ))()
		tractParams[:] = tract_state.T.ravel('F')
		glottisParams = (ctypes.c_double * ( constants[ 'n_glottis_params' ] ))()
		glottisParams[:] = glottis_state.T.ravel('F')
		self.API.vtlSynthesisAddTract( numNewSamples, ctypes.byref( audio ), ctypes.byref( tractParams ), ctypes.byref( glottisParams ) )
		if n_new_samples > 0:
			return np.array( audio )
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	def test( self ): # Run without calling initialize. Completely useless function basically
		speakerFileName = ctypes.c_char_p( self.params.speaker_file_name.encode() )
		audio = (ctypes.c_double * int( self.params.samplerate_audio ) )()
		numSamples = ctypes.c_int(0)
		self.API.vtlApiTest( speakerFileName, ctypes.byref( audio ), ctypes.byref( numSamples ) )
		return np.array( audio )
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	def segment_sequence_to_gestural_score( self, segFilePath, gesFilePath ):
		segFileName = ctypes.c_char_p( segFilePath.encode() )
		gesFileName = ctypes.c_char_p( gesFilePath.encode() )
		self.API.vtlSegmentSequenceToGesturalScore( segFileName, gesFileName )
		if self.params.verbose:
			print('Created gestural score from file: {}'.format( segFilePath ))
		return
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	def gestural_score_to_tract_sequence( self, ges_file_path_list: list,  tract_file_path_list: list = '', return_sequence: bool = False ):
		ges_file_path_list, trac_file_path_list = self._check_if_input_lists_are_valid( [ges_file_path_list, tract_file_path_list], [str, str] )
		args = ( (ges_file_path, tract_file_path, return_sequence) 
			for ges_file_path, tract_file_path in itertools.zip_longest(ges_file_path_list, tract_file_path_list) )
		df_tract_sequence_list = self._run_multiprocessing( self._gestural_score_to_tract_sequence, args )
		return df_tract_sequence_list
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	def _gestural_score_to_tract_sequence( self, args ):
		ges_file_path, tract_file_path, return_sequence = args
		if not os.path.exists( ges_file_path ):
			warnings.warn( 'the specified gestural score file path does not exist: {}. API call will be skipped.'.format( ges_file_path ) )
			return
		gesFileName = ctypes.c_char_p( ges_file_path.encode() )
		if tract_file_path in (None, ''):
			tract_file_path = ges_file_path.split('.')[0] + '_tractSeq.txt'
			logging.info( 'No output file path for tract sequence was specified, saving file to {}'.format( tract_file_path ) )
		tractSequenceFileName = ctypes.c_char_p( tract_file_path.encode() )
		value = self.API.vtlGesturalScoreToTractSequence( gesFileName, tractSequenceFileName )
		if failure != 0:
			raise ValueError('VTL API function vtlGesturalScoreToTractSequence returned the Errorcode: {}. See API doc for more info.'.format(value) )
		if self.params.verbose:
			logging.info( 'Created tractsequence file {} from gestural score file: {}'.format( tract_file_path, ges_file_path ) )
		if return_sequence:
			return self.tract_seq_to_df( tract_file_path )
		return
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	def gestural_score_to_audio( self, ges_file_path: str,  audio_file_path: str = '', return_audio = True, return_n_samples = False, verbose = False ):
		if audio_file_path == '' and return_audio == False:
			print( 'Warning! Function can not return anything.' )
		wavFileName = ctypes.c_char_p( audio_file_path.encode() )
		gesFileName = ctypes.c_char_p( ges_file_path.encode() )
		if return_audio:
			audio = (ctypes.c_double * int( self.get_gestural_score_audio_duration( ges_file_path, return_samples = True ) ))()
		else:
			audio = ctypes.POINTER(ctypes.c_int)() # Null pointer
		numSamples = ctypes.POINTER(ctypes.c_int)() # TODO: add if return_n_samples
		enableConsoleOutput = ctypes.c_int(1) if verbose == True else ctypes.c_int(0)
		self.API.vtlGesturalScoreToAudio( gesFileName, wavFileName, ctypes.byref(audio), ctypes.byref(numSamples), enableConsoleOutput )
		if self.params.verbose:
			print( 'Audio generated from gestural score file: {}'.format( ges_file_path ) )
		#if audio_file_path != '':
		#	self.Write_Wav( self.Normalise_Wav(audio), audio_file_path )
		if return_audio:
			return np.array(audio)
		else:
			return None
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	def gestural_score_to_glottis_signals( self, ges_file_path: str,  sig_file_path: str ):
		gesFileName = ctypes.c_char_p( ges_file_path.encode() )
		glottisSignalsFileName = ctypes.c_char_p( sig_file_path.encode() )
		#print( 'run glottis function' )
		self.API.vtlGesturalScoreToGlottisSignals( gesFileName, glottisSignalsFileName )
		if self.params.verbose:
			print( 'Glottis signals file generated from gestural score file: {}'.format( ges_file_path ) )
		return None
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	def tract_sequence_to_audio( self, tract_seq_path: str, audio_file_path: str = '', return_audio = True, return_n_samples = False, verbose = False ):
		if audio_file_path == '' and return_audio == False:
			print( 'Warning! Function can not return anything.' )
		wavFileName = ctypes.c_char_p( audio_file_path.encode() )
		tractSequenceFileName = ctypes.c_char_p( tract_seq_path.encode() )
		if return_audio:
			audio = (ctypes.c_double * int( self.get_tract_seq_len( tract_seq_path ) * self.params.state_duration * self.params.samplerate_audio ))()
		else:
			audio = ctypes.POINTER(ctypes.c_int)()
		numSamples = ctypes.POINTER(ctypes.c_int)() # TODO: add if return_n_samples
		self.API.vtlTractSequenceToAudio( tractSequenceFileName, wavFileName, ctypes.byref(audio), ctypes.byref(numSamples) )
		if self.params.verbose:
			print('Audio generated: {}'.format(TractSeq_Filename))
		if return_audio:
			return np.array(audio)
		else:
			return None
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	def tract_state_to_limited_tract_state( self, tract_params ):
		constants = self.get_constants()
		tract_param_data = []
		for state in tract_params:
			#print( state )
			#print( state.shape )
			inTractParams = ( ctypes.c_double * ( constants[ 'n_tract_params' ] ) )()
			inTractParams[:] = state.T.ravel('F')
			outTractParams = ( ctypes.c_double * ( constants[ 'n_tract_params' ] ) )()

			self.API.vtlInputTractToLimitedTract( ctypes.byref( inTractParams ), ctypes.byref( outTractParams ) )
			tract_param_data.append( [ np.array( outTractParams ) ] )
		df = pd.DataFrame( tract_param_data, columns = [ 'tract_state' ] )
		return df
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	def get_gestural_score_audio_duration( self, ges_file_path: str, return_samples = True):
		gesFileName = ctypes.c_char_p( ges_file_path.encode() )
		numAudioSamples = ctypes.c_int(0)
		numGestureSamples = ctypes.POINTER(ctypes.c_int)()
		self.API.vtlGetGesturalScoreDuration( gesFileName, ctypes.byref( numAudioSamples ), ctypes.byref( numGestureSamples ) )
		n_samples = numAudioSamples.value
		if return_samples: # returning number of audio samples
			return n_samples
		else: # returning time in seconds
			return n_samples / self.params.samplerate_audio
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	def get_gestural_score_duration( self, ges_file_path: str, return_samples = True):
		gesFileName = ctypes.c_char_p( ges_file_path.encode() )
		numAudioSamples = ctypes.POINTER(ctypes.c_int)()
		numGestureSamples = ctypes.c_int(0)
		self.API.vtlGetGesturalScoreDuration( gesFileName, ctypes.byref( numAudioSamples ), ctypes.byref( numGestureSamples ) )
		n_samples = numGestureSamples.value
		if return_samples: # returning number of audio samples
			return n_samples
		else: # returning time in seconds
			return n_samples / self.params.samplerate_audio
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	def get_tract_seq_len( self, tract_seq_path ):
		with open( tract_seq_path ) as file:
			for index, line in enumerate( file ):
				if index == 7:
					tract_seq_len = int( line.strip() )
					break
		if self.params.verbose:
			print( 'Entries in Tract Sequence file: {}'.format( tract_seq_len ) )
		return tract_seq_len
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	def get_tube_seq_len( self, tube_seq_path ):
		with open( tube_seq_path ) as file:
			for index, line in enumerate( file ):
				if index == 11:
					tube_seq_len = int( line.strip() )
					break
		if self.params.verbose:
			print( 'Entries in Tract Sequence file: {}'.format( tube_seq_len ) )
		return tube_seq_len
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	def tract_seq_to_df( self, tractFilePath ):
		# Skip rows from based on condition
		df_GLP = pd.read_csv( tractFilePath, delim_whitespace = True, skiprows= lambda x: self.read_tract_seq_GLP(x) , header = None )
		df_VTP = pd.read_csv( tractFilePath, delim_whitespace = True, skiprows= lambda x: self.read_tract_seq_VTP(x) , header = None )
		if self.params.verbose:
			print( 'Tract sequence opened: {}'.format( tractFilePath ) )
		return df_GLP, df_VTP
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	def tube_seq_to_df( self, tube_file_path ):
		# Skip rows from based on condition
		df_GLP = pd.read_csv( tube_file_path, delim_whitespace = True, skiprows= lambda x: self.read_tube_seq_GLP(x) , header = None )
		df_tube_param  = pd.read_csv( tube_file_path, delim_whitespace = True, skiprows= lambda x: self.read_tube_seq_param(x) , header = None )
		df_tube_area   = pd.read_csv( tube_file_path, delim_whitespace = True, skiprows= lambda x: self.read_tube_seq_area(x) , header = None )
		df_tube_length = pd.read_csv( tube_file_path, delim_whitespace = True, skiprows= lambda x: self.read_tube_seq_length(x) , header = None )
		df_tube_art    = pd.read_csv( tube_file_path, delim_whitespace = True, skiprows= lambda x: self.read_tube_seq_art(x) , header = None )
		if self.params.verbose:
			print( 'Tube sequence opened: {}'.format( tube_file_path ) )
		return df_GLP, df_tube_param, df_tube_area, df_tube_length, df_tube_art
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	def read_tract_seq_GLP(self, index):
		if (index > 7) and (index % 2 == 0):
			return False
		else:
			return True
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	def read_tract_seq_VTP(self, index):
		if (index > 7) and ((index-1) % 2 == 0):
			return False
		else:
			return True
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	def read_tube_seq_length(self, index):
		if (index > 11) and (index % 5 == 0):
			return False
		else:
			return True
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	def read_tube_seq_art(self, index):
		if (index > 11) and ((index-1) % 5 == 0):
			return False
		else:
			return True
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	def read_tube_seq_GLP(self, index):
		if (index > 11) and ((index-2) % 5 == 0):
			return False
		else:
			return True
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	def read_tube_seq_param(self, index):
		if (index > 11) and ((index-3) % 5 == 0):
			return False
		else:
			return True
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	def read_tube_seq_area(self, index):
		if (index > 11) and ((index-4) % 5 == 0):
			return False
		else:
			return True
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	def df_to_tract_seq( self, tract_file_path, df_GLP, df_VTP, glottis_model = 'Geometric glottis' ):
		f= open( tract_file_path, "w+" )
		f.write("""# The first two lines (below the comment lines) indicate the name of the vocal fold model and the number of states.\n""")
		f.write("""# The following lines contain the control parameters of the vocal folds and the vocal tract (states)\n""")
		f.write("""# in steps of 110 audio samples (corresponding to about 2.5 ms for the sampling rate of 44100 Hz).\n""")
		f.write("""# For every step, there is one line with the vocal fold parameters followed by\n""")
		f.write("""# one line with the vocal tract parameters.\n""")
		f.write("""# \n""")
		f.write("""{}\n""".format( glottis_model ) )
		f.write("""{}\n""".format(len(df_GLP)))
		for index, row in enumerate(df_GLP.index):
			for index2, column in enumerate(df_GLP.columns):
				f.write('{} '.format(df_GLP.iloc[index,index2]))
			f.write('\n')
			for index2, column in enumerate(df_VTP.columns):
				f.write('{} '.format(df_VTP.iloc[index,index2]))
			f.write('\n')
		f.close()
		if self.params.verbose:
			print('Tract Sequence saved as: "{}"'.format( tract_file_path ))
		return
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	def df_to_tube_seq( self, tube_file_path, df_GLP, df_param, df_area, df_length, df_art, glottis_model = 'Geometric glottis' ):
		f= open( tube_file_path, "w+" )
		f.write("""# The first two lines (below the comment lines) indicate the name of the vocal fold model and the number of states.\n""")
		f.write("""# The following lines contain a sequence of states of the vocal folds and the tube geometry\n""")
		f.write("""# in steps of 110 audio samples (corresponding to about 2.5 ms for the sampling rate of 44100 Hz).\n""")
		f.write("""# Each state is represented in terms of five lines:\n""")
		f.write("""# Line 1: glottis_param_0 glottis_param_1 ...\n""")
		f.write("""# Line 2: incisor_position_in_cm, velo_pharyngeal_opening_in_cm^2, tongue_tip_side_elevation[-1...1]\n""")
		f.write("""# Line 3: area0 area1 area2 area3 ... (Areas of the tube sections in cm^2 from glottis to mouth)\n""")
		f.write("""# Line 4: length0 length1 length2 length3 ... (Lengths of the tube sections in cm from glottis to mouth)\n""")
		f.write("""# Line 5: artic0 artic1 artic2 artic3 ... (Articulators of the tube sections between glottis and lips : 1 = tongue; 2 = lower incisors; 3 = lower lip; 4 = other)\n""")
		f.write("""# \n""")
		f.write("""{}\n""".format( glottis_model ) )
		f.write("""{}\n""".format( len(df_GLP) ) )
		for index, row in enumerate(df_GLP.index):
			for index2, column in enumerate(df_GLP.columns):
				f.write('{} '.format(df_GLP.iloc[index,index2]))
			f.write('\n')
			for index2, column in enumerate(df_param.columns):
				f.write('{} '.format(df_param.iloc[index,index2]))
			f.write('\n')
			for index2, column in enumerate(df_area.columns):
				f.write('{} '.format(df_area.iloc[index,index2]))
			f.write('\n')
			for index2, column in enumerate(df_length.columns):
				f.write('{} '.format(df_length.iloc[index,index2]))
			f.write('\n')
			for index2, column in enumerate(df_art.columns):
				f.write('{} '.format(df_art.iloc[index,index2]))
			f.write('\n')
		f.close()
		if self.params.verbose:
			print('Tract Sequence saved as: "{}"'.format( tube_file_path ))
		return
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	def normalize_wav(self, Input_Wav, normalisation = -1 ): #normalisation in dB
		norm_factor = 10**( -1 * normalisation * 0.05 ) -1
		#print(np.max(np.abs(Input_Wav),axis=0))
		norm_max = np.max(np.abs(Input_Wav),axis=0)
		Input_Wav /= ( norm_max + ( norm_max * norm_factor ) )
		#print('Wav file normalised.')
		return Input_Wav
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	def write_wav(self, Input_Wav, Output_Path, Samplerate = None ):
		if Samplerate == None:
			Samplerate = self.params.samplerate_audio
		wav_int = np.int16(Input_Wav * (2**15 - 1))
		wavfile.write(Output_Path, Samplerate, wav_int)
		#print('Wav file saved.')
		return
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	def _run_multiprocessing( self, function, args ):
		pool = mp.Pool( self.params.workers )
		data = pool.map( function, ( (x) for x in args) )
		pool.close()
		pool.join()
		return data
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	def _check_if_input_lists_are_valid( self, input_lists, instances_list ):
		valid_lists = []
		for input_list, instance in zip( input_lists, instances_list ):
			valid_lists.append( self._check_if_list_is_valid( input_list, instance ) )
		#return self._check_lengths_of_input_lists( valid_lists )
		return valid_lists
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	def _check_lengths_of_input_lists( self, input_lists ):
		list_lengths = [ len( input_list ) for input_list in input_lists ]
		if not self._check_if_all_elements_are_equal( list_lengths ):
			warnings.warn( 'input list do not have the same lengths, shorter lists will be padded with "None".' )
			max_length = max( list_lengths )
			#print( max_length )
			for input_list in input_lists:
				while len( input_list ) < max_length:
					input_list.append( None )
			#print( input_list )
		return input_lists
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	def _check_if_object_is_iterable( self, query ):
		try:
			iter( query )
		except TypeError:
			return False
		return True
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	def _check_if_list_is_valid( self, input_list, instance ):
		is_iterable = self._check_if_object_is_iterable( input_list )
		if isinstance( input_list, str ) or is_iterable == False:
			warnings.warn( 'input is either not iterable or a single string. The input gets turned into a list now.' )
			input_list = [ input_list ]
		if input_list and all( isinstance( x, instance ) for x in input_list ):
			return input_list
		else:
			raise TypeError( 'a list containing a non-{} type object was passed, but list of {} was expected.'.format( instance, instance ) )
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	def _check_if_all_elements_are_equal( self, iterable ):
		g = itertools.groupby(iterable)
		return next(g, True) and not next(g, False)

#####################################################################################################################################################


	def tract_state_is_voiced( self, state ):
		return
#---------------------------------------------------------------------------------------------------------------------------------------------------#
#if __name__ == '__main__':
#	freeze_support()
'''