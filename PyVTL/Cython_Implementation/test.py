import VocalTractLabApi as vtl
import tract_sequence as ts
import numpy as np
import matplotlib.pyplot as plt
import time
import soundfile

class wrapper():
	def __init__( self, name ):
		self.name = name
		vtl.initialize( 'JD2.speaker' )
		print( '{} initialized'.format( self.name ) )
	def __del__( self ):
		vtl.close()
		print( '{} closed'.format( self.name ) )


if __name__ == '__main__':
	args = ('1-Aber sehen will sie ihn doch.ges', None, True, -1, 16000, False )
	val = vtl.gestural_score_to_audio( '1-Aber sehen will sie ihn doch.ges', '', workers=4, normalize_audio = -1, sr = 16000, return_data = True )
	plt.plot(val[0])
	plt.show()
	print(val)
'''
stop

#vtl.get_version()
#vtl.automatic_calculation_of_TRX_and_TRY( True )
#vtl.automatic_calculation_of_TRX_and_TRY( False )
#print( vtl.get_param_info( 'tract' ) )
#print( vtl.get_param_info( 'glottis' ) )
#gp = ts.Sub_Glottal_Sequence( np.array( [ vtl.get_param_info( 'glottis' )[ 'neutral' ] ] )  )
#print(gp)
#tp1 = vtl.get_tract_params_from_shape( 'e' )
#tp2 = vtl.get_tract_params_from_shape( 'i' )
#tp = ts.Supra_Glottal_Sequence( np.concatenate( [tp1.tract.to_numpy(), tp2.tract.to_numpy() ] ) )
#print(tp)
##stop
#seq = ts.Tract_Sequence( tp, gp )


gp = ts.Sub_Glottal_Sequence(  np.array( [ vtl.get_param_info( 'glottis' )[ 'neutral' ] for _ in range(0, round(44100/110) ) ] ) )
tp = ts.Supra_Glottal_Sequence(  np.array( [ vtl.get_param_info( 'tract' )[ 'neutral' ] for _ in range(0, round(44100/110) ) ] ) )


seq = ts.Tract_Sequence( tp, gp )
print( seq )

#seq = ts.Tract_Sequence( ts.Supra_Glottal_Sequence( np.ones( shape= (round(44100/110*10), 19 ) ) )

#stop
time_start = time.time()
#for _ in range( 0, 10 ):
audio = vtl._synth_block( (seq, None, False) )
time_end = time.time()
elapsed_time = time_end-time_start

print( 'Elapsed time: {}, time per synthesis: {}'.format( elapsed_time, elapsed_time/10 ) )

soundfile.write( 'cython_test.wav', audio, 44100 )

plt.plot( audio )
plt.show()
'''



'''

stop
vtl.initialize( 'JD2.speaker' )

print( vtl.get_constants() )

vtl.close()

vtl.initialize( 'JD2.speaker' )

vtl.close()
'''