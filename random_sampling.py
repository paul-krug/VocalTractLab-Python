import numpy as np
import pandas as pd
from PyVTL.core import VTL
import random
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
import time


df_formants = pd.read_csv( 'formants_random_uniform_10000.txt' )
df_formants.dropna( inplace = True)
_,_,_,hist2d = plt.hist2d( df_formants[ 'f2' ], df_formants[ 'f1' ], bins = [ 30, 30 ], range=[[100, 800], [400, 1600]], cmap='magma_r' )
plt.colorbar( hist2d )
#plt.scatter( df_formants[ 'f2' ], df_formants[ 'f1' ] )
plt.xlabel( 'f1' )
plt.ylabel( 'f2' )
plt.show()
stop


vtl = VTL()

df_tract_limits = vtl.get_param_info( params = 'tract' )
print( df_tract_limits )

time_1 = time.perf_counter()

n_samples = 10000
n_spectrum_samples = 8192
random_states = []
for param in df_tract_limits.index:
	random_params = []
	for n in range( 0, n_samples ):
		sample = random.uniform( df_tract_limits.loc[ param, 'min' ], df_tract_limits.loc[ param, 'max' ] )
		random_params.append( sample )
	random_states.append( random_params )

#plt.hist( random_states[0][:] )
#plt.show()
#stop

df_random_states = pd.DataFrame( np.array( random_states ).T, columns = df_tract_limits.index )
df_random_states.to_csv( 'random_states_uniform_{}.txt'.format(n_samples), index = False )


df_transfer = vtl.get_transfer_function( np.array( random_states ).T, n_spectrum_samples )
df_transfer.to_json( 'transfer_random_uniform_{}.json'.format(n_samples) )

formants = []
for index in df_transfer.index:
	peaks, _ = find_peaks( df_transfer.loc[ index, 'magnitude' ][ : 500], distance = 5 )
	if len(peaks) >= 2:
		formants.append( [ peaks[0] * 44100/2048, peaks[1] * 44100/n_spectrum_samples ] )
	else:
		formants.append( [np.nan, np.nan])

df_formants = pd.DataFrame( formants, columns = [ 'f1', 'f2' ] ).dropna()
df_formants.to_csv( 'formants_random_uniform_{}.txt'.format(n_samples), index = False )
print( df_formants )

time_2 = time.perf_counter()
print( 'Elapsed time: {}'.format(time_2-time_1) )

#plt.scatter( df_formants[ 'f1' ], df_formants[ 'f2' ] )
plt.hist2d( df_formants[ 'f1' ], df_formants[ 'f2' ] )
plt.show()




'''
peaks, _ = find_peaks( df_transfer['magnitude'][0][ : 500], distance = 5 )
plt.plot( df_transfer['magnitude'][0] )
plt.scatter( peaks, df_transfer['magnitude'][0][peaks] )
plt.show()
print( 'f1: {} Hz'.format( peaks[0] * 44100/2048 ) )
print( 'f2: {} Hz'.format( peaks[1] * 44100/2048 ) )
'''



