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
import warnings
import numpy as np
import pandas as pd
from itertools import chain


#from VocalTractLab import plotting_tools as PT
#from VocalTractLab.plotting_tools import finalize_plot
#from VocalTractLab.plotting_tools import get_plot
#from VocalTractLab.plotting_tools import get_plot_limits
#from VocalTractLab.plotting_tools import get_valid_tiers
#from VocalTractLab import function_tools as FT
#from VocalTractLab.function_tools import is_iterable
##from VocalTractLab import tract_sequence as TS
#from VocalTractLab.tract_sequence import Sub_Glottal_Sequence, Supra_Glottal_Sequence, Motor_Sequence
#from VocalTractLab.audio_tools import get_f0
#import VocalTractLab.VocalTractLabApi as vtl
#from VocalTractLab.target_estimation import fit

from VocalTractLab.vtl_phonemes import phoneme_data
from VocalTractLab.function_tools import check_if_list_is_valid





#---------------------------------------------------------------------------------------------------------------------------------------------------#
def get_phoneme_list(
	data,
	unique_phonemes = True,
	flatten = True,
	):
	if unique_phonemes:
		return [ x[ 0 ] for x in data ]
	else:
		if flatten:
			return list( chain( *data ) )
		else:
			return [ x for x in data ]
#---------------------------------------------------------------------------------------------------------------------------------------------------#
def consonants(
	phoneme_classes = [ 'plosive', 'fricative', 'nasal', 'lateral', 'affricate' ],
	**kwargs,
	):
	phoneme_classes = check_if_list_is_valid( phoneme_classes, str )
	consonants = phoneme_data.loc[
		( phoneme_data[ 'phoneme_type' ] == 'consonant' ) &
		( phoneme_data[ 'phoneme_class' ].isin( phoneme_classes ) )
		].sampa.to_numpy()
	return get_phoneme_list(
		data = consonants,
		**kwargs,
		)
#---------------------------------------------------------------------------------------------------------------------------------------------------#
def single_consonants( **kwargs ):
	return consonants( phoneme_classes = [ 'plosive', 'fricative', 'nasal', 'lateral' ], **kwargs )
#---------------------------------------------------------------------------------------------------------------------------------------------------#
def plosives( **kwargs ):
	return consonants( phoneme_classes = [ 'plosive' ], **kwargs )
#---------------------------------------------------------------------------------------------------------------------------------------------------#
def stops( **kwargs ):
	return consonants( phoneme_classes = [ 'plosive', 'affricate' ], **kwargs )
#---------------------------------------------------------------------------------------------------------------------------------------------------#
def affricates( **kwargs ):
	return consonants( phoneme_classes = [ 'affricate' ], **kwargs )
#---------------------------------------------------------------------------------------------------------------------------------------------------#
def fricatives( **kwargs ):
	return consonants( phoneme_classes = [ 'fricative' ], **kwargs )
#---------------------------------------------------------------------------------------------------------------------------------------------------#
def obstruents( **kwargs ):
	return consonants( phoneme_classes = [ 'plosive', 'fricative', 'affricate' ], **kwargs )
#---------------------------------------------------------------------------------------------------------------------------------------------------#
def sonorants( **kwargs ):
	return consonants( phoneme_classes = [ 'nasal', 'lateral', 'central' ], **kwargs )
#---------------------------------------------------------------------------------------------------------------------------------------------------#
def approximants( **kwargs ):
	return consonants( phoneme_classes = [ 'lateral', 'central' ], **kwargs )
#---------------------------------------------------------------------------------------------------------------------------------------------------#
def nasals( **kwargs ):
	return consonants( phoneme_classes = [ 'nasal' ], **kwargs )
#---------------------------------------------------------------------------------------------------------------------------------------------------#
def laterals( **kwargs ):
	return consonants( phoneme_classes = [ 'lateral' ], **kwargs )
#---------------------------------------------------------------------------------------------------------------------------------------------------#
def diphthongs( **kwargs ):
	return vowels( phoneme_classes = [ 'diphthong' ], **kwargs )
#---------------------------------------------------------------------------------------------------------------------------------------------------#
def monophthongs( **kwargs ):
	return vowels( phoneme_classes = [ 'monophthong' ], **kwargs )
#---------------------------------------------------------------------------------------------------------------------------------------------------#
def vowels(
	phoneme_classes = [ 'monophthong', 'diphthong' ],
	**kwargs,
	):
	phoneme_classes = check_if_list_is_valid( phoneme_classes, str )
	vowels = phoneme_data.loc[
		( phoneme_data[ 'phoneme_type' ] == 'vowel' ) &
		( phoneme_data[ 'phoneme_class' ].isin( phoneme_classes ) )
		].sampa.to_numpy()
	return get_phoneme_list(
		data = vowels,
		**kwargs,
		)
#---------------------------------------------------------------------------------------------------------------------------------------------------#
def phonemes(
	**kwargs,
	):
	vowels = vowels( **kwargs )
	consonants = consonants( **kwargs )
	return chain( vowels, consonants )
#---------------------------------------------------------------------------------------------------------------------------------------------------#
def sampa():
	return phonemes( unique_phonemes = False )
#---------------------------------------------------------------------------------------------------------------------------------------------------#
#print( vowels( phoneme_classes = 'monophthong' ) )
#print( single_consonants() )

#vowels = vowels( phoneme_class = 'monophthong' )
#print( vowels )