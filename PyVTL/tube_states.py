import PyVTL.VocalTractLabApi as vtl
import PyVTL.plotting_tools as PT
from PyVTL.plotting_tools import finalize_plot
from PyVTL.plotting_tools import get_plot
from PyVTL.plotting_tools import get_plot_limits
import numpy as np
import matplotlib.pyplot as plt






#####################################################################################################################################################
class Tube_State():
	def __init__( self, 
		          tube_length,
		          tube_area,
		          tube_articulator,
		          incisor_position,
		          tongue_tip_side_elevation,
		          velum_opening,
		          ):
		self.tube_length = tube_length
		self.tube_area = tube_area
		self.tube_articulator = tube_articulator
		self.incisor_position = incisor_position
		self.tongue_tip_side_elevation = tongue_tip_side_elevation
		self.velum_opening = velum_opening
		return
#---------------------------------------------------------------------------------------------------------------------------------------------------#
	def plot( self, 
	          ax = None, 
			  **kwargs,
	          ):
		figure, ax = get_plot( 1, ax )
		tube_x = [ self.tube_length[ 0 ] ]
		for length in self.tube_length[ 1: ]:
			tube_x.append( tube_x[ -1 ] + length )
		x = np.arange( 0, np.sum( self.tube_length ), 0.01 )
		y = []
		tmp_length = 0
		for index, _ in enumerate( self.tube_length ): 
			for val in x:
				if val >= tmp_length:
					if val <= tube_x[ index ]:
						y.append( self.tube_area[ index ] )
					else:
						tmp_length = tube_x[ index ]
						break
		ax.set( xlabel = 'Tube Length [cm]', ylabel = r'Cross-sectional Area [cm$^2$]' )
		#y = [ val for val in x  ]
		#x = [ self.tube_length[ 0 ] ]
		#for length in self.tube_length[ 1: ]:
		#	x.append( x[ -1 ] + length )
		ax.plot( x, y )
		finalize_plot( figure, ax, **kwargs )
		#ax.set( xlabel = 'Tube Length [cm]', ylabel = r'Cross-sectional Area [cm$^2$]' )
		return ax
#---------------------------------------------------------------------------------------------------------------------------------------------------#