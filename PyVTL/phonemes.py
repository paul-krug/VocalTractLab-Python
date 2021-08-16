#####################################################################################################################################################
#---------------------------------------------------------------------------------------------------------------------------------------------------#
#	- This file is a part of the Python module PyVTL
#---------------------------------------------------------------------------------------------------------------------------------------------------#
#
#	- Copyright (C) 2021, Paul Konstantin Krug, Dresden, Germany
#	- https://github.com/TUD-STKS/PyVTL
#	- Author: Paul Konstantin Krug, TU Dresden
#
#	- License:
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
# pip install praat-parselmouth
# pip install dtw-python
# Requirements:
#	- python 3 (tested with version 3.7)
#	- numpy    (tested with version 1.19.5)
#	- pandas   (tested with version 1.2.1)
#	- scipy    (tested with version 1.6.0)
#
# Optional, used for visualization:
#	- matplotlib        (tested with version 3.3.3)
#	- praat-parselmouth (tested with version 0.3.3)
#---------------------------------------------------------------------------------------------------------------------------------------------------#
#####################################################################################################################################################
#---------------------------------------------------------------------------------------------------------------------------------------------------#


import os
import numpy as np
import pandas as pd










#vtl_sampa = [
## Vowels
#"a:", "e:", "i:", "o:", "u:", "E:", "2:", "y:",
#"a", "e", "i", "o", "u", "E", "2", "y",
#"I", "O", "U", "9", "Y", "@", "6",
## Diphthongs
#"aI", "aU", "OY",
#"i:6", "i6", "I6", "y:6", "y6", "Y6", "e:6", "e6", "E6", "E:6", 
#"2:6", "26", "96", "a:6", "a6", "u:6", "u6", "U6", "o:6", "o6", "O6",
## Plosives
#"?", "p", "b", "t", "d", "k", "g",
## Fricatives
#"f", "v", "T", "D", "s", "z", "S", "Z", "C", "j", "x", "r", "R", "h",
## Affricates
#"pf", "ts", "tS", "dZ",
## Nasals and lateral
#"m", "n", "N", "l"
#]


vtl_sampa = [
'i:6', 'y:6', 'e:6', 'E:6', '2:6', 'a:6', 'u:6', 'o:6',
'a:', 'e:', 'i:', 'o:', 'u:', 'E:', '2:', 'y:', 'aI', 'aU', 'OY', 'i6', 'I6', 'y6', 'Y6', 'e6', 'E6', '26', '96', 'a6', 'u6', 'U6', 'o6', 'O6', 'pf', 'ts', 'tS', 'dZ',
'a', 'e', 'i', 'o', 'u', 'E', '2', 'y', 'I', 'O', 'U', '9', 'Y', '@', '6', '?', 'p', 'b', 't', 'd', 'k', 'g', 'f', 'v', 'T', 'D', 's', 'z', 'S', 'Z', 'C', 'j', 'x', 'r', 'R', 'h', 'm', 'n', 'N', 'l'
]

vtl_sampa_dict = {
	'diphthongs': [ 'i:6', 'y:6', 'e:6', 'E:6', '2:6', 'a:6', 'u:6', 'o:6', 'aI', 'aU', 'OY', 'i6', 'I6', 'y6', 'Y6', 'e6', 'E6', '26', '96', 'a6', 'u6', 'U6', 'o6', 'O6', ],
	'vowels':     [ 'a:', 'e:', 'i:', 'o:', 'u:', 'E:', '2:', 'y:', 'a', 'e', 'i', 'o', 'u', 'E', '2', 'y', 'I', 'O', 'U', '9', 'Y', '@', '6', ],
	'consonants': [ 'pf', 'ts', 'tS', 'dZ', '?', 'p', 'b', 't', 'd', 'k', 'g', 'f', 'v', 'T', 'D', 's', 'z', 'S', 'Z', 'C', 'j', 'x', 'r', 'R', 'h', 'm', 'n', 'N', 'l' ],
}

graphem_to_phoneme = {
	'ch': { 'phon': 'x', 'pos': [] }
}


vtl_sampa = [
'i:6', 'y:6', 'e:6', 'E:6', '2:6', 'a:6', 'u:6', 'o:6',
'a:', 'e:', 'i:', 'o:', 'u:', 'E:', '2:', 'y:', 'aI', 'aU', 'OY', 'i6', 'I6', 'y6', 'Y6', 'e6', 'E6', '26', '96', 'a6', 'u6', 'U6', 'o6', 'O6', 'pf', 'ts', 'tS', 'dZ',
'a', 'e', 'i', 'o', 'u', 'E', '2', 'y', 'I', 'O', 'U', '9', 'Y', '@', '6', '?', 'p', 'b', 't', 'd', 'k', 'g', 'f', 'v', 'T', 'D', 's', 'z', 'S', 'Z', 'C', 'j', 'x', 'r', 'R', 'h', 'm', 'n', 'N', 'l'
]

vtl_phonemes = [
{ 'name': 'i:6', 'type': 'vowel'    , 'class': 'diphthong'   },
{ 'name': 'y:6', 'type': 'vowel'    , 'class': 'diphthong'   },
{ 'name': 'e:6', 'type': 'vowel'    , 'class': 'diphthong'   },
{ 'name': 'E:6', 'type': 'vowel'    , 'class': 'diphthong'   },
{ 'name': '2:6', 'type': 'vowel'    , 'class': 'diphthong'   },
{ 'name': 'a:6', 'type': 'vowel'    , 'class': 'diphthong'   },
{ 'name': 'u:6', 'type': 'vowel'    , 'class': 'diphthong'   },
{ 'name': 'o:6', 'type': 'vowel'    , 'class': 'diphthong'   },
{ 'name': 'a:' , 'type': 'vowel'    , 'class': 'monophthong' },
{ 'name': 'e:' , 'type': 'vowel'    , 'class': 'monophthong' },
{ 'name': 'i:' , 'type': 'vowel'    , 'class': 'monophthong' },
{ 'name': 'o:' , 'type': 'vowel'    , 'class': 'monophthong' },
{ 'name': 'u:' , 'type': 'vowel'    , 'class': 'monophthong' },
{ 'name': 'E:' , 'type': 'vowel'    , 'class': 'monophthong' },
{ 'name': '2:' , 'type': 'vowel'    , 'class': 'monophthong' },
{ 'name': 'y:' , 'type': 'vowel'    , 'class': 'monophthong' },
{ 'name': 'aI' , 'type': 'vowel'    , 'class': 'diphthong'   },
{ 'name': 'aU' , 'type': 'vowel'    , 'class': 'diphthong'   },
{ 'name': 'OY' , 'type': 'vowel'    , 'class': 'diphthong'   },
{ 'name': 'i6' , 'type': 'vowel'    , 'class': 'diphthong'   },
{ 'name': 'I6' , 'type': 'vowel'    , 'class': 'diphthong'   },
{ 'name': 'y6' , 'type': 'vowel'    , 'class': 'diphthong'   },
{ 'name': 'Y6' , 'type': 'vowel'    , 'class': 'diphthong'   },
{ 'name': 'e6' , 'type': 'vowel'    , 'class': 'diphthong'   },
{ 'name': 'E6' , 'type': 'vowel'    , 'class': 'diphthong'   },
{ 'name': '26' , 'type': 'vowel'    , 'class': 'diphthong'   },
{ 'name': '96' , 'type': 'vowel'    , 'class': 'diphthong'   },
{ 'name': 'a6' , 'type': 'vowel'    , 'class': 'diphthong'   },
{ 'name': 'u6' , 'type': 'vowel'    , 'class': 'diphthong'   },
{ 'name': 'U6' , 'type': 'vowel'    , 'class': 'diphthong'   },
{ 'name': 'o6' , 'type': 'vowel'    , 'class': 'diphthong'   },
{ 'name': 'O6' , 'type': 'vowel'    , 'class': 'diphthong'   },
{ 'name': 'pf' , 'type': 'consonant', 'class': 'affricate'   },
{ 'name': 'ts' , 'type': 'consonant', 'class': 'affricate'   },
{ 'name': 'tS' , 'type': 'consonant', 'class': 'affricate'   },
{ 'name': 'dZ' , 'type': 'consonant', 'class': 'affricate'   },
{ 'name': 'a'  , 'type': 'vowel'    , 'class': 'monophthong' },
{ 'name': 'e'  , 'type': 'vowel'    , 'class': 'monophthong' },
{ 'name': 'i'  , 'type': 'vowel'    , 'class': 'monophthong' },
{ 'name': 'o'  , 'type': 'vowel'    , 'class': 'monophthong' },
{ 'name': 'u'  , 'type': 'vowel'    , 'class': 'monophthong' },
{ 'name': 'E'  , 'type': 'vowel'    , 'class': 'monophthong' },
{ 'name': '2'  , 'type': 'vowel'    , 'class': 'monophthong' },
{ 'name': 'y'  , 'type': 'vowel'    , 'class': 'monophthong' },
{ 'name': 'I'  , 'type': 'vowel'    , 'class': 'monophthong' },
{ 'name': 'O'  , 'type': 'vowel'    , 'class': 'monophthong' },
{ 'name': 'U'  , 'type': 'vowel'    , 'class': 'monophthong' },
{ 'name': '9'  , 'type': 'vowel'    , 'class': 'monophthong' },
{ 'name': 'Y'  , 'type': 'vowel'    , 'class': 'monophthong' },
{ 'name': '@'  , 'type': 'vowel'    , 'class': 'monophthong' },
{ 'name': '6'  , 'type': 'vowel'    , 'class': 'monophthong' },
{ 'name': '?'  , 'type': 'consonant', 'class': 'plosive'     },
{ 'name': 'p'  , 'type': 'consonant', 'class': 'plosive'     },
{ 'name': 'b'  , 'type': 'consonant', 'class': 'plosive'     },
{ 'name': 't'  , 'type': 'consonant', 'class': 'plosive'     },
{ 'name': 'd'  , 'type': 'consonant', 'class': 'plosive'     },
{ 'name': 'k'  , 'type': 'consonant', 'class': 'plosive'     },
{ 'name': 'g'  , 'type': 'consonant', 'class': 'plosive'     },
{ 'name': 'f'  , 'type': 'consonant', 'class': 'fricative' },
{ 'name': 'v'  , 'type': 'consonant', 'class': 'fricative' },
{ 'name': 'T'  , 'type': 'consonant', 'class': 'fricative' },
{ 'name': 'D'  , 'type': 'consonant', 'class': 'fricative' },
{ 'name': 's'  , 'type': 'consonant', 'class': 'fricative' },
{ 'name': 'z'  , 'type': 'consonant', 'class': 'fricative' },
{ 'name': 'S'  , 'type': 'consonant', 'class': 'fricative' },
{ 'name': 'Z'  , 'type': 'consonant', 'class': 'fricative' },
{ 'name': 'C'  , 'type': 'consonant', 'class': 'fricative' },
{ 'name': 'j'  , 'type': 'consonant', 'class': 'fricative' },
{ 'name': 'x'  , 'type': 'consonant', 'class': 'fricative' },
{ 'name': 'r'  , 'type': 'consonant', 'class': 'fricative' },
{ 'name': 'R'  , 'type': 'consonant', 'class': 'fricative' },
{ 'name': 'h'  , 'type': 'consonant', 'class': 'fricative' },
{ 'name': 'm'  , 'type': 'consonant', 'class': 'nasal' },
{ 'name': 'n'  , 'type': 'consonant', 'class': 'nasal' },
{ 'name': 'N'  , 'type': 'consonant', 'class': 'nasal' },
{ 'name': 'l'  , 'type': 'consonant', 'class': 'lateral' },
]
#df_phonemes = pd.DataFrame( vtl_phonemes, columns=['name','type','class'] )
df_phonemes = pd.read_csv( os.path.join( os.path.dirname(__file__), 'Data/vtl_phonemes.txt' ), sep=',', skiprows=4 )
#print( df_phonemes )
#stop
def get():
	#_phonemes = df_phonemes.map(str.strip)
	return df_phonemes

#---------------------------------------------------------------------------------------------------------------------------------------------------#

#####################################################################################################################################################
#---------------------------------------------------------------------------------------------------------------------------------------------------#
class Phoneme():
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	"""PyVTL phonemes""" 
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	def __init__( self, name: str ):
		exists = name in df_phonemes.name.values
		#print( df_phonemes.name )
		#print( exists )
		if exists:
			phon = df_phonemes[ df_phonemes.name == name ]
			self.name = name
			self._type = phon[ 'type' ].values[0]
			self._class = phon[ 'class' ].values[0]
			self._voiced = 'bruuh'
		else:
			print( 'Cannot create the phoneme: {}. It is not supported by VTL!'.format( name ) )
			self.__del__()
		return
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	def __del__( self ):
		#print( 'Phoneme destroyed.' )
		return
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	def __str__( self ):
		return 'name: {}, type: {}, class: {}, voiced: {}'.format( self.name, self._type, self._class, self._voiced )
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	def is_obstruent( self ):
		return self._type == 'consonant' and ( self._class == 'affricate' or self._class == 'fricative' or self._class == 'plosive' )
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	def is_sonorant( self ):
		return self._type == 'consonant' and ( self._class == 'lateral' or self._class == 'nasal' )
#---------------------------------------------------------------------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------------------------------------------------------------------#
#####################################################################################################################################################



#####################################################################################################################################################
#---------------------------------------------------------------------------------------------------------------------------------------------------#
class Consonant( Phoneme ):
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	"""PyVTL consonants""" 
#---------------------------------------------------------------------------------------------------------------------------------------------------#
#####################################################################################################################################################



#####################################################################################################################################################
#---------------------------------------------------------------------------------------------------------------------------------------------------#
class Vowel( Phoneme ):
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	"""PyVTL vowels""" 
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	def __init__( self ):
		self._name
		self._type
		self._class
		self._voiced = True
#---------------------------------------------------------------------------------------------------------------------------------------------------#
#####################################################################################################################################################


#a = Phoneme( 'i:6' )
#print(a)