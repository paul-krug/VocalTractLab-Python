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
import os
#import multiprocessing as mp
import warnings
#import PyVTL
import PyVTL.core
import pandas as pd
import numpy as np
#---------------------------------------------------------------------------------------------------------------------------------------------------#
#####################################################################################################################################################







#####################################################################################################################################################
class Sub_Glottal_Sequence():
	def __init__( self, states: np.array, api: PyVTL.core.VTL_API = None ):
		if api == None:
			warnings.warn( "No API was passed. Using default settings now to create tract sequence object, might be incompatible with your current settings" )
			api = PyVTL.core.VTL_API()
		self.constants = api.get_constants()
		if len( states.shape ) != 2:
			raise ValueError( "Shape of passed state is not two-dimensional. The shape should be (x,y), x: no. states, y: no. features" )
		if states.shape[ 1 ] != self.constants[ 'n_glottis_params' ]:
			raise ValueError( "Dimension of features is {}, but should be {}.".format( states.shape[ 1 ], self.constants[ 'n_glottis_params' ] ) )
		self.param_info = api.get_param_info( 'glottis' )
		self.glottis = pd.DataFrame( states, columns = self.param_info.index )
		api.close()
		return
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	@classmethod
	def from_tract_file( cls, tract_file_path ):
		return
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	def __str__( self, ):
		return str( self.glottis )
#####################################################################################################################################################



#####################################################################################################################################################
class Supra_Glottal_Sequence():
	def __init__( self, states: np.array, api: PyVTL.core.VTL_API = None ):
		if api == None:
			warnings.warn( "No API was passed. Using default settings now to create tract sequence object, might be incompatible with your current settings" )
			api = PyVTL.core.VTL_API()
		self.constants = api.get_constants()
		if len( states.shape ) != 2:
			raise ValueError( "Shape of passed state is not two-dimensional. The shape should be (x,y), x: no. states, y: no. features" )
		if states.shape[ 1 ] != self.constants[ 'n_tract_params' ]:
			raise ValueError( "Dimension of features is {}, but should be {}.".format( states.shape[ 1 ], self.constants[ 'n_tract_params' ] ) )
		self.param_info = api.get_param_info( 'tract' )
		self.tract = pd.DataFrame( states, columns = self.param_info.index )
		return
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	@classmethod
	def from_tract_file( cls, tract_file_path ):
		return
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	def __str__( self, ):
		return str( self.tract )
#####################################################################################################################################################



#####################################################################################################################################################
class Tract_Sequence():
	def __init__( self, tract_states: Supra_Glottal_Sequence, Sub_Glottal_Sequence ):
		if not isinstance( tract_states, Supra_Glottal_Sequence ):
			raise TypeError( '{} type object was passed, but {} was expected.'.format( tract_states, Supra_Glottal_Sequence ) )
		if not isinstance( glottis_states, Sub_Glottal_Sequence ):
			raise TypeError( '{} type object was passed, but {} was expected.'.format( glottis_states, Sub_Glottal_Sequence ) )
		for key in tract_states.constants.keys:
			if tract_states.constants[ key ] != glottis_states.constants[ key ]:
				raise ValueError( 'API constant {} is different for supra and sub glottal state sequence.'.format( key ) )

		self.param_info = { 'tract': tract_states.param_info, 'glottis': glottis_states.param_info }
		self.tract = tract_states.tract
		self.glottis = glottis_states.glottis
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	@classmethod
	def from_tract_file( cls, tract_file_path ):
		return
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	def __str__( self, ):
		return str( pd.concatenate( [ self.tract, self.glottis ] ) )
#####################################################################################################################################################