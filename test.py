#import PyVTL
from PyVTL.PyVTL import VTL
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