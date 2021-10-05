import os
import sys
import pandas as pd
WORKING_PATH = os.getcwd()
sys.path.append( os.path.join( os.path.dirname(__file__), 'data', 'dictionaries', 'phonecodes' ) )
import phonecodes as pc

vtl_phonemes = pd.read_csv( os.path.join( os.path.dirname(__file__), 'data', 'vtl_phonemes.txt' ), sep = ',', skiprows = 4 )
#print( vtl_phonemes )

def text_to_phonemes( sentences, phoneme_style = 'vtl', language = 'en' ):
	available_styles = [ 'ipa', 'arpabet', 'xsampa', 'vtl', 'vtlsampa', 'vtl_sampa', 'disc', 'callhome' ]
	if phoneme_style not in available_styles:
		raise ValueError( 'The specified phoneme style must be one of the following: {}, not {}!'.format( available_styles, phoneme_style ) )
	elif phoneme_style in [ 'vtl', 'vtlsampa', 'vtl_sampa' ]:
		alphabet = 'xsampa'
	else:
		alphabet = phoneme_style
	if language == 'en':
		out = _en_to_phonemes( sentences, alphabet )
	if phoneme_style in [ 'vtl', 'vtlsampa', 'vtl_sampa' ]:
		out = [ get_vtl_phonemes( x ) for x in out ]
	return out

def _en_to_phonemes( sentences, phoneme_style ):
	from g2p_en import G2p
	g2p = G2p()
	phonemes_list = [ g2p( sentence ) for sentence in sentences ]
	if phoneme_style == 'arpabet':
		return phonemes_list
	else:
		phonemes_list = [ [ pc.convert( x, 'arpabet', 'ipa', language='en' ) for x in phonemes ] for phonemes in phonemes_list ]
	if phonemes_list != 'ipa':
		phonemes_list = [ [ pc.convert( x, 'ipa', phoneme_style, language='en' ) for x in phonemes ] for phonemes in phonemes_list ]
	return phonemes_list



def get_vtl_phonemes( sampa, drop = [ ',','^','_','.','/','\\','"','%','[',']','=',"'",'!' ] ):
	phones = []
	phoneme_list=[]
	#print( sampa)
	for entry in sampa:
		for char in drop:
			entry = entry.replace(char, '')
		phones.append( entry )
	#print('phones: {}'.format( phones ) )
	return phones
	#phonemes = ''.join( phones )
	found= False
	for _phonemes in phones:
		tmp_list=[]
		while len(_phonemes) > 0:
			found= False
			for entry in phonemes.vtl_sampa:
				if found == False and _phonemes.startswith( entry ):
					tmp_list.append( entry )
					_phonemes = _phonemes[ len( entry ): ]
					found = True
					#print(phonemes)
		phoneme_list.append( tmp_list )
	print(phoneme_list)

	index = 0
	while index < len( phoneme_list )-1:
		if ( phoneme_list[ index ][ -1 ] in phonemes.vtl_sampa_dict[ 'vowels' ] ) and ( phoneme_list[ index + 1 ][ 0 ] in phonemes.vtl_sampa_dict[ 'vowels' ] ):
			phoneme_list.insert( index + 1, '?' )
		index += 1
	phoneme_list_flat = list( itertools.chain( *phoneme_list ) )

	if ( phoneme_list_flat[0] in phonemes.vtl_sampa_dict['vowels'] ) or ( phoneme_list_flat[0] in phonemes.vtl_sampa_dict['diphthongs'] ):
		phoneme_list_flat.insert(0, '?')
	print( 'List of phonemes: {}'.format( phoneme_list_flat ) )
	return phoneme_list_flat