import numpy as np
import librosa







#---------------------------------------------------------------------------------------------------------------------------------------------------#
def detect_onset_and_offset( audio_file_path, threshold = 0.01, sr = 44100 ):
	data, samplerate = librosa.load( audio_file_path, sr )
	over_threshold = []
	for index, sample in enumerate( np.abs( data ) ):
		if sample > threshold:
			over_threshold.append( index )
	return [ over_threshold[0] / samplerate - 0.005, over_threshold[-1] / samplerate + 0.005]
#---------------------------------------------------------------------------------------------------------------------------------------------------#