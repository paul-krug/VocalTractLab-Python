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
import os
import sys
import pandas as pd
from VocalTractLab.function_tools import check_if_list_is_valid
WORKING_PATH = os.getcwd()
sys.path.append( os.path.join( os.path.dirname(__file__), 'data', 'dictionaries', 'phonecodes' ) )
import phonecodes as pc
#---------------------------------------------------------------------------------------------------------------------------------------------------#
#####################################################################################################################################################


vtl_phonemes = pd.read_csv( os.path.join( os.path.dirname(__file__), 'data', 'vtl_phonemes.txt' ), sep = ',', skiprows = 4 )
#print( vtl_phonemes )

def text_to_phonemes( sentences, phoneme_style = 'vtl', language = 'en' ):
	sentences = check_if_list_is_valid( sentences, str )
	available_styles = [ 'ipa', 'arpabet', 'xsampa', 'vtl', 'vtlsampa', 'vtl_sampa', 'disc', 'callhome' ]
	if phoneme_style not in available_styles:
		raise ValueError( 'The specified phoneme style must be one of the following: {}, not {}!'.format( available_styles, phoneme_style ) )
	elif phoneme_style in [ 'vtl', 'vtlsampa', 'vtl_sampa' ]:
		alphabet = 'xsampa'
	else:
		alphabet = phoneme_style
	if language == 'de':
		out = _de_to_phonemes( sentences, alphabet )
	elif language == 'en':
		direct_out = _en_to_phonemes( sentences, alphabet )
		out = []
		out_phones = []
		tmp = []
		for sentence in direct_out:
			out_phones = []
			tmp = []
			for x in sentence:
				if x not in [ ' ', '.' ]:
					tmp.append( x )
				else:
					vtl_tmp = get_vtl_phonemes( tmp )
					if vtl_tmp != []:
						out_phones.append( vtl_tmp )
					tmp = []
			out.append( out_phones )
		#print( out )
	#if phoneme_style in [ 'vtl', 'vtlsampa', 'vtl_sampa' ]:
	#	out = [ get_vtl_phonemes( x ) for sentence in out for x in sentence ]
	return out

def _de_to_phonemes( sentences, phoneme_style ):
	phonemes_list = [ dictionary_lookup( sentence, language = 'de' ) for sentence in sentences ]
	if phoneme_style == 'arpabet':
		return phonemes_list
	else:
		phonemes_list = [ [ pc.convert( x, 'arpabet', 'ipa', language='en' ) for x in phonemes ] for phonemes in phonemes_list ]
	if phonemes_list != 'ipa':
		phonemes_list = [ [ pc.convert( x, 'ipa', phoneme_style, language='en' ) for x in phonemes ] for phonemes in phonemes_list ]
	return phonemes_list

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
	#print(phoneme_list)

	index = 0
	while index < len( phoneme_list )-1:
		if ( phoneme_list[ index ][ -1 ] in phonemes.vtl_sampa_dict[ 'vowels' ] ) and ( phoneme_list[ index + 1 ][ 0 ] in phonemes.vtl_sampa_dict[ 'vowels' ] ):
			phoneme_list.insert( index + 1, '?' )
		index += 1
	phoneme_list_flat = list( itertools.chain( *phoneme_list ) )

	if ( phoneme_list_flat[0] in phonemes.vtl_sampa_dict['vowels'] ) or ( phoneme_list_flat[0] in phonemes.vtl_sampa_dict['diphthongs'] ):
		phoneme_list_flat.insert(0, '?')
	#print( 'List of phonemes: {}'.format( phoneme_list_flat ) )
	return phoneme_list_flat


def lookup( text, remove_punctuation = True, language = 'de' ):
	if remove_punctuation:
		for char in [ ',','^','_','.','/','\\','"','%','[',']','=',"'",'!','?' ]:
			text = text.replace(char, '')
	dictionary = pd.read_csv( 'Dictionaries/de-ipa.txt', sep='\t' )
	dictionary['de'] = dictionary['de'].apply( lambda x: x.lower() )
	#print( dictionary)
	input_words = text.split(' ')
	#input_words = [ x.lower() for x in input_words ]
	input_words = [ x for x in input_words if x not in [ None, '', ' ' ] ]
	input_words = [ x.lower() for x in input_words ]
	print( input_words )
	ipa_words = []
	#ipa_sub_words = []
	#valid_sub_words = []
	for word in input_words:
		occurences = dictionary.loc[ ( dictionary['de'] == word ) ]
		if len( occurences ) > 0:
			ipa_words.append( occurences[0] )
		else:
			ipa_words.append( None )
	return ipa_words