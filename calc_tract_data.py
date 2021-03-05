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
import PyVTL
from PyVTL.core import VTL
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.colors import LogNorm



vtl = VTL()
df = pd.read_pickle( 'data/random_tract_states_uniform_1M.pickle.bz2', compression = 'bz2' )
df[ 'VO' ]  = -0.1
df[ 'TRX' ] = -1000
df[ 'TRY' ] = -1000

#####################################################################################################################################################
#---------------------------------------------------------------------------------------------------------------------------------------------------#
def run( index, _slice ):
	print( 'Run nr: {}, slice [ {}:{} ]'.format( index, _slice[0], _slice[1] ) )
	states = df.loc[ _slice[0]:_slice[1] ].to_numpy()
	print( 'Calculating transfer functions.')
	df_transfer = vtl.get_transfer_function( states, n_spectrum_samples = 8192 )
	print( '  -> Done.')
	print( 'Calculating tube data.')
	df_tube_data= vtl.tract_params_to_tube_data( states )
	print( '  -> Done.')
	print( 'Concatenate the data frames.')
	df_transfer[ 'tract_state' ] = list( states )
	df_data = pd.concat( [df_transfer, df_tube_data], axis = 1 )
	del states
	del df_transfer
	del df_tube_data
	print( df_data.tract_state[0][7] )
	print( df_data.tract_state[0][14] )
	print( df_data.tract_state[0][15] )
	print( df_data.columns )
	df_data.to_pickle( 'data/data_random_tract_states_uniform_1M_{}.pickle.bz2'.format( index ), compression = 'bz2' )
	print( '  -> Done.')
	data_reload = pd.read_pickle( 'data/data_random_tract_states_uniform_1M_{}.pickle.bz2'.format( index ), compression = 'bz2' )
	print( data_reload )
	return
#---------------------------------------------------------------------------------------------------------------------------------------------------#
#####################################################################################################################################################



slices = [
           [0, 50000], 
           [50001, 100000],
           #[100001, 150000],
		   #[150001, 200000],
		   #[200001, 250000],
		   #[250001, 300000],
		   #[300001, 350000],
		   #[350001, 400000],
		   #[400001, 450000],
		   #[450001, 500000],
		 ]
for index, _slice in enumerate( slices ):
	run( index, _slice )