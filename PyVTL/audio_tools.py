import soundfile
import numpy as np


def normalize( audio, normalization = -1 ): #normalisation in dB
	norm_factor = 10**( -1 * normalization * 0.05 ) -1
	norm_max = np.max( np.abs( audio ),axis=0)
	audio /= ( norm_max + ( norm_max * norm_factor ) )
	return audio

def write( audio, audio_file_path, sr ):
	soundfile.write( audio_file_path, audio, sr )
	return