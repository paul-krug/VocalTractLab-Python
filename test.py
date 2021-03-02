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
from ast import literal_eval


vtl = VTL()
#vtl.get_tract_params_from_shape( shape = 'a' )
df_tract = vtl.get_param_info( params = 'tract' )
df_glottis = vtl.get_param_info( params = 'glottis' )

#tract_states = np.array( [ df_tract[ 'neutral' ] for x in range( 0, 500 ) ] )
#glottis_states = np.array( [ df_glottis[ 'neutral' ] for x in range( 0, 500 ) ] )
#print( tract_states.shape )
#print( glottis_states.shape )

#stop

#audio_synth_block = vtl.synth_block( tract_states, glottis_states, verbose=True )
#plt.plot( audio )
#plt.show()

#vtl.export_tract_svg( tract_states, 'test')
#df = vtl.tract_params_to_tube_data( tract_states )
#print( df )

#df = vtl.get_transfer_function( tract_states )
#print( df )
#plt.plot( df['phase_rad'][0] )
#plt.show()


#vtl.synthesis_reset()
#audio_1 = vtl.synthesis_add_state( tract_states[0], glottis_states[0], n_new_samples = 0 )
#audio_2 = vtl.synthesis_add_state( tract_states[0], glottis_states[0], n_new_samples = 110 )
#print( audio_1 )
#plt.plot( audio_2 )
#plt.plot( audio_synth_block )
#plt.xlim( 0, 250 )
#plt.show()



#####################################################################################################################################################
#---------------------------------------------------------------------------------------------------------------------------------------------------#
def formant_plot():
	vowels = [
	{ 'name': 'a', 'f1': 668, 'f2': 1191 },
	{ 'name': 'e', 'f1': 327, 'f2': 2157 },
	{ 'name': 'i', 'f1': 208, 'f2': 2152 },
	{ 'name': 'o', 'f1': 335, 'f2': 628 },
	{ 'name': 'u', 'f1': 254, 'f2': 796 },
	]
	dfs = [ pd.read_csv( 'formants_random_uniform_run_vo_{}_100000.txt'.format( x ) ) for x in ['q'] ]#, 'r', 's', 't'] ]
	df_formants = pd.concat( dfs, ignore_index = True )

	dfs = [ pd.read_csv( 'random_states_uniform_run_vo_{}_100000.txt'.format( x ) ) for x in ['q'] ]#, 'r', 's', 't'] ]
	df_states = pd.concat( dfs, ignore_index = True )
	df_states[ 'VO' ] = -0.1
	#print( df_states )

	#print( df_states.to_numpy().shape )
	#stop
	#df_tube_data = vtl.tract_params_to_tube_data( tract_params = df_states.to_numpy() )
	#df_tube_data.to_pickle( 'tube_data_run_vo_q_100000.pickle' )
	df_tube_data = pd.read_pickle( 'tube_data_run_vo_q_100000.pickle' )
	print( df_tube_data.shape )
	print( df_formants.shape )
	print( df_states.shape )
	#stop
	#df_tube_data = df_tube_data.applymap(lambda x: x.replace( '\n', '') if isinstance(x, str) else x)
	#df_tube_data = df_tube_data.applymap(lambda x: x.replace( '"', '') if isinstance(x, str) else x)
	#df_tube_data = df_tube_data.applymap(lambda x: x.replace( '[', '') if isinstance(x, str) else x)
	#df_tube_data = df_tube_data.applymap(lambda x: x.replace( ']', '') if isinstance(x, str) else x)
	#df_tube_data = df_tube_data.astype( np.float64 )
	#print( df_tube_data )
	#stop
	#print( df_tube_data )
	voiced_index = []
	for index in df_tube_data.index:
		#print( np.min( df_tube_data.loc[ index, 'tube_area_cm2' ] ) )
		#print( df_tube_data.loc[ index, 'tube_articulator' ] ) 
		#stop
		#print( index )
		#print( df_tube_data.loc[ index, 'tube_area_cm2' ].replace('\n','').replace(']','').replace('[','').replace('"','').split(" ").remove( '' ) )
		#tube_area = df_tube_data.loc[ index, 'tube_area_cm2' ].replace('\n','').replace(']','').replace('[','').replace('"','').split(" ")
		#tube_area = np.array( [ x for x in tube_area if x != '' ] ).astype(np.float)
		#print( np.array( df_tube_data.loc[ index, 'tube_area_cm2' ].split(' ') ).astype(np.float) )
		if np.min( df_tube_data.loc[ index, 'tube_area_cm2' ] ) > 0.1:
			if index in df_formants.index:
				voiced_index.append( index )
	print( len( voiced_index ) )
	df_voiced_formants = df_formants.iloc[ voiced_index ]
	print( df_voiced_formants )
	#stop
	n_spectrum_samples = 8192
	for index in df_voiced_formants.index:
		if df_voiced_formants.loc[ index, 'f1' ] > 1000:
			print( index )
			print( df_voiced_formants.loc[ index, 'f1' ] )
			#print( df_states.loc[ index ].values )
			transfer = vtl.get_transfer_function( np.array( [ df_states.loc[ index ].values ] ), n_spectrum_samples = 8192 )['magnitude'][0][:500]
			#print( transfer)
			plt.plot( [ x * 44100/n_spectrum_samples for x in range( 0, len(transfer) ) ], transfer )
			plt.show()
			plt.plot( df_tube_data.loc[ index, 'tube_area_cm2' ] )
			plt.show()
	stop
	#print( df_voiced_formants )
	df_voiced_formants.dropna( inplace = True)
	_,_,_,hist2d = plt.hist2d(df_voiced_formants[ 'f1' ], df_voiced_formants[ 'f2' ], bins=100, range=[[0, 2700], [0, 3000]], cmap='magma_r', norm=LogNorm(), density=True )
	plt.colorbar( hist2d )
	#plt.scatter( df_formants[ 'f1' ], df_formants[ 'f2' ],  )
	plt.xlabel( 'f1' )
	plt.ylabel( 'f2' )
	for vowel in vowels:
		plt.annotate( vowel['name'], ( vowel['f1'], vowel['f2'] ), horizontalalignment='right', verticalalignment='top' )
	plt.tight_layout()
	plt.savefig( 'formant_voicing_run_vo_qrst_g0.3.png' )
	plt.show()
	stop


	return
#---------------------------------------------------------------------------------------------------------------------------------------------------#
#####################################################################################################################################################

formant_plot()