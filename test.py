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


vtl = VTL()
#vtl.get_tract_params_from_shape( shape = 'a' )
df_tract = vtl.get_param_info( params = 'tract' )
df_glottis = vtl.get_param_info( params = 'glottis' )

tract_states = np.array( [ df_tract[ 'neutral' ] for x in range( 0, 500 ) ] )
glottis_states = np.array( [ df_glottis[ 'neutral' ] for x in range( 0, 500 ) ] )

print( tract_states.shape )
print( glottis_states.shape )

stop

audio_synth_block = vtl.synth_block( tract_states, glottis_states, verbose=True )
#plt.plot( audio )
#plt.show()

#vtl.export_tract_svg( tract_states, 'test')
#df = vtl.tract_params_to_tube_data( tract_states )
#print( df )

#df = vtl.get_transfer_function( tract_states )
#print( df )
#plt.plot( df['phase_rad'][0] )
#plt.show()


vtl.synthesis_reset()
audio_1 = vtl.synthesis_add_state( tract_states[0], glottis_states[0], n_new_samples = 0 )
audio_2 = vtl.synthesis_add_state( tract_states[0], glottis_states[0], n_new_samples = 110 )
print( audio_1 )
plt.plot( audio_2 )
plt.plot( audio_synth_block )
plt.xlim( 0, 250 )
plt.show()