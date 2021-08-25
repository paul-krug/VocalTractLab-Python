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
import warnings
import pandas as pd
import numpy as np
import PyVTL.VocalTractLabApi as vtl
import matplotlib.pyplot as plt
#---------------------------------------------------------------------------------------------------------------------------------------------------#
#####################################################################################################################################################







#####################################################################################################################################################
class Sub_Glottal_Sequence():
	def __init__( self, states: np.array, name: str = 'sequence.sub_glottal' ):
		self.constants = vtl.get_constants()
		if len( states.shape ) != 2:
			raise ValueError( "Shape of passed state is not two-dimensional. The shape should be (x,y), x: no. states, y: no. features" )
		if states.shape[ 1 ] != self.constants[ 'n_glottis_params' ]:
			raise ValueError( "Dimension of features is {}, but should be {}.".format( states.shape[ 1 ], self.constants[ 'n_glottis_params' ] ) )
		self.param_info = vtl.get_param_info( 'glottis' )
		self.name = name
		self.glottis = pd.DataFrame( states, columns = self.param_info.index )
		self.length = len( self.glottis.index )
		return
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	@classmethod
	def from_tract_file( cls, tract_file_path ):
		df_GLP = pd.read_csv( tract_file_path, delim_whitespace = True, skiprows= lambda x: read_tract_seq_GLP(x) , header = None )
		return cls( df_GLP.to_numpy(), tract_file_path )
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	def __str__( self, ):
		return str( self.glottis )
#####################################################################################################################################################



#####################################################################################################################################################
class Supra_Glottal_Sequence():
	def __init__( self, states: np.array, name: str = 'sequence.supra_glottal' ):
		self.constants = vtl.get_constants()
		if len( states.shape ) != 2:
			raise ValueError( "Shape of passed state is not two-dimensional. The shape should be (x,y), x: no. states, y: no. features" )
		if states.shape[ 1 ] != self.constants[ 'n_tract_params' ]:
			raise ValueError( "Dimension of features is {}, but should be {}.".format( states.shape[ 1 ], self.constants[ 'n_tract_params' ] ) )
		self.param_info = vtl.get_param_info( 'tract' )
		self.name = name
		self.tract = pd.DataFrame( states, columns = self.param_info.index )
		self.length = len( self.tract.index )
		return
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	@classmethod
	def from_tract_file( cls, tract_file_path ):
		df_VTP = pd.read_csv( tract_file_path, delim_whitespace = True, skiprows= lambda x: read_tract_seq_VTP(x) , header = None )
		return cls( df_VTP.to_numpy(), tract_file_path )
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	def __str__( self, ):
		return str( self.tract )
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	def apply_biomechanical_constraints():
		self.tract = vtl.tract_sequence_to_limited_tract_sequence( self.tract ).tract
#####################################################################################################################################################



#####################################################################################################################################################
class Tract_Sequence():
	def __init__( self, tract_states: Supra_Glottal_Sequence, glottis_states: Sub_Glottal_Sequence, name: str = 'sequence.tract' ):
		if not isinstance( tract_states, Supra_Glottal_Sequence ):
			raise TypeError( '{} type object was passed, but {} was expected.'.format( tract_states, Supra_Glottal_Sequence ) )
		if not isinstance( glottis_states, Sub_Glottal_Sequence ):
			raise TypeError( '{} type object was passed, but {} was expected.'.format( glottis_states, Sub_Glottal_Sequence ) )
		for key in tract_states.constants:
			if tract_states.constants[ key ] != glottis_states.constants[ key ]:
				raise ValueError( 'API constant {} is different for supra and sub glottal state sequence.'.format( key ) )
		self.param_info = { 'tract': tract_states.param_info, 'glottis': glottis_states.param_info }
		self.name = name
		self.tract = tract_states.tract
		self.glottis = glottis_states.glottis
		lengths_difference = np.abs( tract_states.length - glottis_states.length )
		if tract_states.length > glottis_states.length:
			warnings.warn( 'lengths of supra glottal sequence is longer than sub glottal sequence. Will pad the sub glottal sequence now.' )
			self.glottis = pd.concat( 
			                             [ self.glottis, 
										   pd.DataFrame( [ self.glottis.iloc[ -1, : ] for _ in range(0, lengths_difference ) ] )
										 ], 
										 ignore_index = True 
										 )
		elif tract_states.length < glottis_states.length:
			warnings.warn( 'lengths of supra glottal sequence is shorter than sub glottal sequence. Will pad the supra glottal sequence now.' )
			self.tract = pd.concat( 
			                             [ self.tract, 
										   pd.DataFrame( [ self.tract.iloc[ -1, : ] for _ in range(0, lengths_difference ) ] )
										 ], 
										 ignore_index = True 
										 )
		if len( self.tract.index ) != len( self.glottis.index ):
			print( 'ultra fail' )
		self.length = len( self.tract.index )
		return
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	@classmethod
	def from_tract_file( cls, tract_file_path ):
		df_GLP = pd.read_csv( tract_file_path, delim_whitespace = True, skiprows= lambda x: read_tract_seq_GLP(x) , header = None )
		df_VTP = pd.read_csv( tract_file_path, delim_whitespace = True, skiprows= lambda x: read_tract_seq_VTP(x) , header = None )
		return cls( Supra_Glottal_Sequence( df_VTP.to_numpy() ), Sub_Glottal_Sequence( df_GLP.to_numpy() ), tract_file_path )
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	def __str__( self, ):
		return str( pd.concat( [ self.tract, self.glottis ], axis = 1 ) )
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	def apply_biomechanical_constraints():
		self.tract = vtl.tract_sequence_to_limited_tract_sequence( self.tract ).tract
#####################################################################################################################################################





#---------------------------------------------------------------------------------------------------------------------------------------------------#
def read_tract_seq_GLP( index ):
	if (index > 7) and (index % 2 == 0):
		return False
	else:
		return True
#---------------------------------------------------------------------------------------------------------------------------------------------------#
def read_tract_seq_VTP( index ):
	if (index > 7) and ((index-1) % 2 == 0):
		return False
	else:
		return True
#---------------------------------------------------------------------------------------------------------------------------------------------------#



#####################################################################################################################################################
class Target_Sequence():
	def __init__( self, states: np.array, name: str = 'sequence.supra_glottal' ):
		self.constants = vtl.get_constants()
		if len( states.shape ) != 2:
			raise ValueError( "Shape of passed state is not two-dimensional. The shape should be (x,y), x: no. states, y: no. features" )
		if states.shape[ 1 ] != self.constants[ 'n_tract_params' ]:
			raise ValueError( "Dimension of features is {}, but should be {}.".format( states.shape[ 1 ], self.constants[ 'n_tract_params' ] ) )
		self.param_info = vtl.get_param_info( 'tract' )
		self.name = name
		self.tract = pd.DataFrame( states, columns = self.param_info.index )
		self.length = len( self.tract.index )
		return
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	@classmethod
	def from_tract_file( cls, tract_file_path ):
		df_VTP = pd.read_csv( tract_file_path, delim_whitespace = True, skiprows= lambda x: read_tract_seq_VTP(x) , header = None )
		return cls( df_VTP.to_numpy(), tract_file_path )
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	def __str__( self, ):
		return str( self.tract )
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	def apply_biomechanical_constraints():
		self.tract = vtl.tract_sequence_to_limited_tract_sequence( self.tract ).tract
#####################################################################################################################################################