#####################################################################################################################################################
#---------------------------------------------------------------------------------------------------------------------------------------------------#
#	- This file is a part of the VocalTractLab Python module PyVTL, see https://github.com/paul-krug/VocalTractLab
#---------------------------------------------------------------------------------------------------------------------------------------------------#
#
#	- Copyright (C) 2021, Paul Konstantin Krug, Dresden, Germany
#	- https://github.com/paul-krug/VocalTractLab
#	- Author: Paul Konstantin Krug, TU Dresden
#
#	- License info:
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

#---------------------------------------------------------------------------------------------------------------------------------------------------#
#####################################################################################################################################################
#---------------------------------------------------------------------------------------------------------------------------------------------------#
# Load essential packages:
from VocalTractLab.g2p import text_to_phonemes
from VocalTractLab.function_tools import check_if_list_is_valid
from VocalTractLab.segment_sequence import Segment_Sequence
from VocalTractLab.VocalTractLabApi import gestural_score_to_audio
from VocalTractLab.VocalTractLabApi import segment_sequence_to_gestural_score



def text_to_speech( sentence_list, language = 'en', **kwargs ):
	sentence_list = check_if_list_is_valid( sentence_list, str )
	phonemes_list = text_to_phonemes( sentence_list, phoneme_style = 'vtl_sampa', language = language )
	segment_sequences = get_phone_durations( phonemes_list, method = 'naive', language = language )
	for index, seg in enumerate( segment_sequences ):
		seg.to_seg_file( 'tts_tmp_{}.seg'.format( index ) )
	segment_sequence_to_gestural_score( [ 'tts_tmp_{}.seg'.format( index ) for index, _ in enumerate( segment_sequences ) ] )
	audio_data = gestural_score_to_audio( [ 'tts_tmp_{}.ges'.format( index ) for index, _ in enumerate( segment_sequences ) ],
		                                  save_file = False,
		                                  normalize_audio = -1,
		                                  sr = 16000,
		                                  return_data = True,
		                                  )
	return audio_data
def tts( text, language = 'de', **kwargs ):
	return text_to_speech( text, language = language, **kwargs )

def get_phone_durations( phonemes_list, method = 'naive', language = 'en' ):
	segment_sequences = []
	for phonemes in phonemes_list:
		phonemes_flat = []
		for _list in phonemes:
			phonemes_flat += _list
		segment_sequence = Segment_Sequence( phonemes_flat, [ 0.2 for _ in phonemes_flat ] )
		segment_sequences.append( segment_sequence )
	return segment_sequences
