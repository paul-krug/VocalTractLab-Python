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

import warnings
import numpy as np
import pandas as pd
from scipy.special import binom
from scipy.special import factorial
import matplotlib.pyplot as plt
from itertools import zip_longest
from itertools import chain
from collections import Counter

from VocalTractLab import plotting_tools as PT
from VocalTractLab.plotting_tools import finalize_plot
from VocalTractLab.plotting_tools import get_plot
from VocalTractLab.plotting_tools import get_plot_limits
from VocalTractLab.plotting_tools import get_valid_tiers
from VocalTractLab import function_tools as FT
from VocalTractLab.function_tools import is_iterable
#from VocalTractLab import tract_sequence as TS
from VocalTractLab.tract_sequence import Sub_Glottal_Sequence, Supra_Glottal_Sequence, Motor_Sequence
from VocalTractLab.audio_tools import get_f0
import VocalTractLab.VocalTractLabApi as vtl
from VocalTractLab.target_estimation import fit



def vowels(
	type,
	):
	phonemes = pd.read_csv()
	return