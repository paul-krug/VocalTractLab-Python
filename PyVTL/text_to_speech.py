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
from PyVTL.g2p import text_to_phonemes



def text_to_speech( sentence_list, language = 'de', **kwargs ):
	sentence_list = FT.check_if_input_list_is_valid( sentence_list, str )
	phonemes_list = text_to_phonemes( sentence_list, phoneme_style = 'vtl_sampa', language = language )
	durations_list = get_phone_durations( phonemes_list, method = 'naive', language = language )
	segment_sequences = Segment_Sequence( phonemes_list, durations )
	audio = segment_sequence.to_audio_file( **kwargs )
	return audio
def tts( text, language = 'de', **kwargs ):
	return text_to_speech( text, language = language, **kwargs )
