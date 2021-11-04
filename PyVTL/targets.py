
import numpy as np
from scipy.special import binom
from scipy.special import factorial
import matplotlib.pyplot as plt


#####################################################################################################################################################
#---------------------------------------------------------------------------------------------------------------------------------------------------#
class Target():
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	"""PyVTL articulatory target""" 
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	def __init__( self, onset_time, duration, slope, offset, time_constant ):
		self.onset_time = onset_time
		self.offset_time = onset_time + duration
		self.duration = duration
		self.slope = slope
		self.offset = offset
		self.time_constant = time_constant
		self.tau = time_constant
		return
	# def __init__( self, onset_time, offset_time, slope, offset, time_constant ):
	# 	self.onset_time = onset_time
	# 	self.offset_time = offset_time
	# 	self.duration = offset_time - onset_time
	# 	self.slope = slope
	# 	self.offset = offset
	# 	self.time_constant = time_constant
	# 	self.tau = time_constant
	# 	return
#---------------------------------------------------------------------------------------------------------------------------------------------------#
#####################################################################################################################################################


'''
#####################################################################################################################################################
#---------------------------------------------------------------------------------------------------------------------------------------------------#
class Synchronized_Target():
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	"""PyVTL articulatory target""" 
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	def __init__( self, onset_time, offset_time, slope, offsets, time_constant ):
		self.targets = []
		for offset in offsets:
			self.targets.append( Target( onset_time, offset_time, slope, offset, time_constant ) )
		return
#---------------------------------------------------------------------------------------------------------------------------------------------------#
#####################################################################################################################################################
'''


#####################################################################################################################################################
#---------------------------------------------------------------------------------------------------------------------------------------------------#
class Target_Sequence():
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	"""PyVTL articulatory target""" 
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	def __init__( self, targets ):
		self.targets = targets
		return
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	def plot( self, ax = None, plot_contour = True, plot_targets = True, show = True ):
		if ax == None:
			figure, ax = plt.subplots( 1, figsize = (8, 4/3) )#, sharex = True, gridspec_kw = {'hspace': 0} )
		if plot_contour:
			tam = tg.Target_Approximation_Model()
			contour = tam.response( self.targets )
			ax.plot( contour[ :, 0 ], contour[ :, 1 ], color = 'navy' )
		if plot_targets:
			ax.axvline( self.targets[0].onset_time, color = 'black' )
			for tar in self.targets:
				ax.axvline( tar.offset_time, color = 'black' )
				ax.plot( [ tar.onset_time, tar.offset_time ], 
				         [ tar.slope * (tar.onset_time-tar.onset_time) + tar.offset, tar.slope * (tar.offset_time-tar.onset_time) + tar.offset ],
				         color = 'black', linestyle='--' )
		ax.set( xlabel = 'Time [s]', ylabel = self.name )
		if show:
			plt.tight_layout()
			plt.show()
		return
#####################################################################################################################################################



#####################################################################################################################################################
#---------------------------------------------------------------------------------------------------------------------------------------------------#
class Synchronous_Target_Sequence():
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	"""PyVTL articulatory target""" 
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	def __init__( self, 
		          targets = None,
		          onset_duration = 0,
		          duration = 0.5,
		          slope = 0,
		          offset = 0,
		          time_constant = 0.015,
		          name = '',
		          ):
		self.targets = []
		for offset in offsets:
			self.targets.append( Target( onset_time, offset_time, slope, offset, time_constant ) )
		return
#---------------------------------------------------------------------------------------------------------------------------------------------------#
#####################################################################################################################################################


#####################################################################################################################################################
#---------------------------------------------------------------------------------------------------------------------------------------------------#
class Motor_Score():
#---------------------------------------------------------------------------------------------------------------------------------------------------#
#
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	def __init__( self,
		          supra_glottal_target_sequence,
		          sub_glottal_target_sequence,
		          ):
		self.tract_targets = supra_glottal_target_sequence.tract_targets
		self.glottis_targets = sub_glottal_target_sequence.glottis_targets
		return
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	def plot( self, parameters = ['LP','JA','LD','HX','HY'], n_params = 19 ):
		figure, axs = plt.subplots( len(parameters), figsize = (8, 4/3 *len(parameters) ), sharex = True, gridspec_kw = {'hspace': 0} )
		for index, parameter in enumerate( parameters ):
			axs[ index ].plot( self.tract_targets.loc[ :, parameter ] )
			axs[ index ].set( ylabel = parameter )
		plt.xlabel( 'Tract state' )
		for ax in axs:
		    ax.label_outer()
		figure.align_ylabels( axs[:] )
		plt.show()
		return
#---------------------------------------------------------------------------------------------------------------------------------------------------#

#---------------------------------------------------------------------------------------------------------------------------------------------------#
#####################################################################################################################################################



#####################################################################################################################################################
#---------------------------------------------------------------------------------------------------------------------------------------------------#
class Target_Approximation_Model():
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	"""PyVTL articulatory target""" 
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	def __init__( self, order = 5 ):
		self.FILTERORDER = order
		return
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	def response( self, target_sequence, discretization_rate = 44100 / 110, onset_state = None  ) -> np.array:
		trajectory = []
		start = target_sequence[ 0 ].onset_time
		end   = target_sequence[ -1 ].offset_time
		duration = end - start
		n_samples = duration * discretization_rate
		sample_times = np.arange( start, end, duration / n_samples )
		if onset_state == None:
			onset_state = target_sequence[0].offset
		current_state = [ onset_state ]
		for _ in range( 1, self.FILTERORDER ):
			current_state.append( 0.0 )

		b_begin = target_sequence[ 0 ].onset_time
		b_end = b_begin

		sample_index = 0

		for target in target_sequence:
			b_begin = b_end
			b_end = b_begin + target.duration
			c = self.calculate_coefficients( target, current_state )
			while( sample_times[ sample_index ] <= b_end ):
				constant = 0.0
				t = sample_times[ sample_index ] - b_begin
				for n in range( 0, self.FILTERORDER ):
					constant += c[ n ] * ( t**n )
				time = sample_times[ sample_index ]
				value= constant * np.exp( - (1/target.time_constant) * t ) + target.slope * t + target.offset
				trajectory.append( [ time, value ] )
				sample_index += 1
				if sample_index >= len( sample_times ):
					break
			current_state = self.calculate_state( current_state, b_end, b_begin, target );

		return np.array( trajectory )
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	def calculate_coefficients( self, target, current_state ):
		coefficients = [ 0 for _ in current_state ]
		assert len( coefficients ) == self.FILTERORDER, 'Sometimes size does matter bro...'
		coefficients[ 0 ] = current_state[ 0 ] - target.offset
		for n in range( 1, self.FILTERORDER ):
			acc = 0
			for i in range( 0, n ):
				acc += ( coefficients[ i ] * ( (-1 / target.time_constant)**(n - i) ) * binom( n, i ) * factorial( i ) )
			if n == 1:
				acc += target.slope # adaption for linear targets; minus changes in following term!
			coefficients[ n ] = ( current_state[ n ] - acc ) / factorial( n )
		return coefficients
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	def calculate_state( self, state, time, start_time, target ):
		t = time - start_time
		state_update = [ 0 for _ in range( 0, self.FILTERORDER ) ]
		c = self.calculate_coefficients( target, state)
		for n in range( 0, self.FILTERORDER ):
			acc = 0
			for i in range( 0, self.FILTERORDER ):
				q = 0
				for k in range( 0, np.min( [ self.FILTERORDER - i, n + 1 ] ) ):
					q += ( ( (-1 / target.time_constant)**(n - k) ) * binom(n, k) * c[i + k] * factorial(k + i) / factorial(i) )
				acc += ( (t**i) * q );
			state_update[ n ] = acc * np.exp( -( 1 / target.time_constant) * t)
		# correction for linear targets
		if (self.FILTERORDER > 1):
			state_update[ 0 ] += (target.offset + target.slope * t)
		if (self.FILTERORDER > 2):
			state_update[ 1 ] += target.slope
		return state_update
#---------------------------------------------------------------------------------------------------------------------------------------------------#
#	def binomial( self, n, k ):
#		result = 1
#		tmp = k
#		if (tmp > n - tmp)
#			tmp = n - tmp;
#		for i in range( 0, tmp ):
#			result *= (n - i)
#			result /= (i + 1)
#	return result
##---------------------------------------------------------------------------------------------------------------------------------------------------#
#	def factorial( self, ):

#####################################################################################################################################################



#target_1 = Target( onset_time = 0.1, offset_time= 0.3, slope=0.0, offset=0, time_constant=0.005  )
#target_2 = Target( onset_time = 0.3, offset_time= 0.7, slope=0.0, offset=3000, time_constant=0.005  )
#target_3 = Target( onset_time = 0.7, offset_time= 0.9, slope=0.0, offset=0, time_constant=0.015  )
#target_sequence = [target_1, target_2, target_3]
#tam = Target_Approximation_Model()
#trajectory = tam.response( target_sequence )
#plt.plot( trajectory[:,0], trajectory[:,1] )
#plt.show()