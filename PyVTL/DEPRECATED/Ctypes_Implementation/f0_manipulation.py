#####################################################################################################################################################
#---------------------------------------------------------------------------------------------------------------------------------------------------#
#	- This file is a part of the Python module PyVTL, see https://github.com/TUD-STKS/PyVTL
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
# Requirements:
#	- python 3 (tested with version 3.7)
#	- numpy    (tested with version 1.19.5)
#	- pandas   (tested with version 1.2.1)
#	- scipy    (tested with version 1.6.0)
#---------------------------------------------------------------------------------------------------------------------------------------------------#
#####################################################################################################################################################
#---------------------------------------------------------------------------------------------------------------------------------------------------#
# Load essential packages:
import pandas as pd
import numpy as np
from PyVTL import PyVTL as pv
#---------------------------------------------------------------------------------------------------------------------------------------------------#
#####################################################################################################################################################



#####################################################################################################################################################
#---------------------------------------------------------------------------------------------------------------------------------------------------#
# Define some important constants:
params = pv.vtl_params()
Frame_Duration = params.state_duration
Samplerate_Audio = params.samplerate_audio
#---------------------------------------------------------------------------------------------------------------------------------------------------#
#####################################################################################################################################################



#####################################################################################################################################################
#---------------------------------------------------------------------------------------------------------------------------------------------------#
class F0_Manipulation():
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	"""A class for the microprosody experiment""" 
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	def __init__(self, Input_Filename_Obstruents = None, verbose = False):
		self.CF0_Duration = 0.08 # 80ms
		return
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	def __del__(self):
		print('F0-Manipulation closed.')
		return
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	def Calc_CF0(self, index, start_index):
		return 11.1 *np.exp(-95.9*(index-start_index))
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	def Calc_CF0_N(self, index, start_index):
		return 11.1 *np.exp(-0.23975*(index-start_index))
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	def Calc_IF0_N(self, index, start_index, end_index):
		return -10.75*np.exp(-((index-start_index)/(np.abs(start_index-end_index)) -0.49)**2/ (0.26**2))
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	def get_obstruents( self, segFilePath ):
		obstruent_list = []
		value_list =[]
		current_time_list =[]
		current_time = 0
		previous_time = current_time
		previous_time_2 = previous_time
		with open( segFilePath ) as file:
			for line in file:
				if len(line.strip()) != 0 :
					items = line.strip().split(';')
					#print(items)
					for item in [x for x in items if x]:
						label = item.split('=')[0].strip()
						value = item.split('=')[1].strip()
						#print('Label: {}'.format(label))
						#print('Value: {}'.format(value))
						if label =='name' and (value not in [None," ", ""]):
							value_list.append(value)
						if label == 'duration_s':
							current_duration = float(value)
							previous_time_2 = previous_time
							previous_time = current_time
							current_time += current_duration
							#print('c time {}'.format(current_time))
							current_time_list.append(current_time)
						if label == 'prosody_type':
							if value == 'IF0':
								obstruent_list.append([previous_time, current_time, value])
							elif value == 'CF0':
								obstruent_list.append([previous_time_2, previous_time, value])
							else:
								raise ValueError('Error in F0_Manipulation::get_obstruents: \
								                  Wrong prosody_type label: {}, but should be IF0 or CF0.'.format(value))
		df = pd.DataFrame(obstruent_list, columns= ['Start', 'End', 'Type'])
		#print(df)
		return df, value_list, current_time_list
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	def manipulate_F0(self, df_signal, df_timestamps, amplitude = 1.0 ):
		df_signal = df_signal.copy()
		for index, row in enumerate(df_timestamps.index):
			if df_timestamps.iloc[index,:]['Type'] == 'IF0':
				d_1 = float(df_timestamps.iloc[index,:]['Start']) / Frame_Duration
				d_2 = float(df_timestamps.iloc[index,:]['End']) / Frame_Duration
				Start = round( d_1 - 0.4 * (d_2 - d_1) )
				End   = round( d_2 + 0.2 * (d_2 - d_1) )
				for i in range(Start, End):
					df_signal.iloc[i] += amplitude * self.Calc_IF0_N(index = i, start_index= Start, end_index=End)
			elif df_timestamps.iloc[index,:]['Type'] == 'CF0':
				Start_Consonant = round( float( df_timestamps.iloc[index,:]['Start'] ) / Frame_Duration )
				Start_CF0       = round( float( df_timestamps.iloc[index,:]['End'] ) / Frame_Duration )
				for i in range( Start_Consonant, Start_CF0 ):
					df_signal.iloc[i] += 11.1 / 2 * amplitude * ( 1-np.cos( 1/(Start_CF0-Start_Consonant) *np.pi * (i - Start_Consonant) ) )
				for i in range( Start_CF0, round( Start_CF0 + self.CF0_Duration/Frame_Duration ) ):
					df_signal.iloc[i] = df_signal.iloc[i] + amplitude * self.Calc_CF0_N(index = i, start_index= Start_CF0)
			else:
				print('Error, effect not specified.')
		print('Tract Sequence manipulated.')
		return df_signal
#---------------------------------------------------------------------------------------------------------------------------------------------------#
#####################################################################################################################################################