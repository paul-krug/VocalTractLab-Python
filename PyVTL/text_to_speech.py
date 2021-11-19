
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
