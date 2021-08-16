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
import tqdm

import pandas as pd
import numpy as np
from scipy.io import wavfile
import matplotlib.pyplot as plt
#---------------------------------------------------------------------------------------------------------------------------------------------------#
#####################################################################################################################################################



#####################################################################################################################################################
class Tract_Sequence():
	def __init__( self, tract_states: np.array, glottis_states: np.array, api: PyVTL.core.VTL_API = None ):
		if api == None:
			warnings.warn( """No api was passed. Using the default settings now to create tract sequence object.
								This object may be incompatible with your current settings""" )
			api = PyVTL.core.VTL_API()
		self.constants = api.get_constants()

		if len( tract_states.shape ) != 2 or len( glottis_state ):

		if not isinstance( tract_state, Tract_State ):
			raise TypeError( '{} type object was passed, but {} was expected.'.format( tract_state, Tract_State ) )
		if not isinstance( glottis_state, Glottis_State ):
			raise TypeError( '{} type object was passed, but {} was expected.'.format( glottis_state, Glottis_State ) )
		self.tract = tract_state
		self.glottis = glottis_state
#####################################################################################################################################################






#####################################################################################################################################################
class Supra_Glottal_Sequence():
	def __init__( self, states: np.array, api: PyVTL.core.VTL_API = None ):
		if api == None:
			warnings.warn( """No api was passed. Using the default settings now to create tract sequence object.
								This object may be incompatible with your current settings""" )
			api = PyVTL.core.VTL_API()
		self.constants = api.get_constants()

		if len( states.shape ) != 2:
			raise ValueError( "Shape of passed state is not two-dimensional. The shape should be (x,y), x: no. states, y: no. features" )



		if not isinstance( tract_state, Tract_State ):
			raise TypeError( '{} type object was passed, but {} was expected.'.format( tract_state, Tract_State ) )
		if not isinstance( glottis_state, Glottis_State ):
			raise TypeError( '{} type object was passed, but {} was expected.'.format( glottis_state, Glottis_State ) )
		self.tract = tract_state
		self.glottis = glottis_state
#####################################################################################################################################################
















#####################################################################################################################################################
class Synthesis_State():
	def __init__( self, tract_state, glottis_state ):
		if not isinstance( tract_state, Tract_State ):
			raise TypeError( '{} type object was passed, but {} was expected.'.format( tract_state, Tract_State ) )
		if not isinstance( glottis_state, Glottis_State ):
			raise TypeError( '{} type object was passed, but {} was expected.'.format( glottis_state, Glottis_State ) )
		self.tract = tract_state
		self.glottis = glottis_state
#####################################################################################################################################################


#####################################################################################################################################################
#---------------------------------------------------------------------------------------------------------------------------------------------------#
class Tract_State():
#---------------------------------------------------------------------------------------------------------------------------------------------------#
#####################################################################################################################################################
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	def __init__( self, api, state_params ):
		state_params = _check_if_list_is_valid( state_params, ( int, float ) )
		self.n_tract_params = api._get_constants()[ 'n_tract_params' ]
		param_info = api._get_params( 'tract' )
		if len( state_params ) != self.n_tract_params:
			raise ValueError( 'The number of passed tract parameters does not equal the number of tract parameters defined in the speaker file.' )
		param_info[ 'value' ] = state_params
		self.params = param_info
		return
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	def __str__( self, ):
		return str( self.params )

#####################################################################################################################################################

class Glottis_State():
	def __init__( self, api, state_params ):
		state_params = _check_if_list_is_valid( state_params, ( int, float ) )

		constants = api._get_constants()
		self.n_tract_params = constants[ 'n_tract_params' ]
		self.n_glottis_params = constants[ 'n_glottis_params' ]

		param_info = api._get_params( 'tract' )

		if len( state_params ) != n_tract_params:
			raise ValueError( 'The number of passed tract parameters does not equal the number of tract parameters defined in the speaker file.' )
		param_info[ 'value' ] = state_params
		self.params = param_info

	@classmethod
	def from_tract_file( cls, tract_file_path ):


#---------------------------------------------------------------------------------------------------------------------------------------------------#
class State_Sequence():
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
	def gestural_score_to_tract_sequence( self, ges_file_path_list: list,  tract_file_path_list: list = '', return_sequence: bool = False ):
		ges_file_path_list, trac_file_path_list = self._check_if_input_lists_are_valid( [ges_file_path_list, tract_file_path_list], [str, str] )
		args = [ [ges_file_path, tract_file_path, return_sequence]
			for ges_file_path, tract_file_path in itertools.zip_longest( ges_file_path_list, tract_file_path_list ) ]
		df_tract_sequence_list = self._run_multiprocessing( '_gestural_score_to_tract_sequence', args )
		return df_tract_sequence_list
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