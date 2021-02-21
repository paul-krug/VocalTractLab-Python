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

audio = vtl.synth_block( tract_states, glottis_states, verbose=True )
plt.plot( audio )
plt.show()