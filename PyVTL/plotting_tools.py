import numpy as np
import matplotlib.pyplot as plt
from PyVTL.function_tools import make_output_path
from PyVTL.function_tools import is_iterable



state_plot_kwargs = dict( HX = dict( color = '#02010e' ),
                          HY = dict( color = '#02010e' ),
                          JX = dict( color = 'navy' ),
                          JA = dict( color = 'navy' ),
                          LP = dict( color = 'rebeccapurple' ),
                          LD = dict( color = 'rebeccapurple' ),
                          VS = dict( color = 'darkmagenta' ),
                          VO = dict( color = 'darkmagenta' ),
                          TCX = dict(  color = 'darkorange' ),
                          TCY = dict(  color = 'darkorange' ),
                          TTX = dict(  color = 'darkorange' ),
                          TTY = dict(  color = 'darkorange' ),
                          TBX = dict(  color = 'darkorange' ),
                          TBY = dict(  color = 'darkorange' ),
                          TRX = dict(  color = 'darkorange' ),
                          TRY = dict(  color = 'darkorange' ),
                          TS1 = dict(  color = 'darkorange' ),
                          TS2 = dict(  color = 'darkorange' ),
                          TS3 = dict(  color = 'darkorange' ),
                          F0 = dict( color = 'teal' ),
                          PR = dict( color = 'teal' ),
                          XB = dict( color = 'turquoise' ),
                          XT = dict( color = 'turquoise' ),
                          CA = dict( color = 'turquoise' ),
                          PL = dict( color = 'turquoise' ),
                          RA = dict( color = 'turquoise' ),
                          DP = dict( color = 'turquoise' ),
                          PS = dict( color = 'turquoise' ),
                          FL = dict( color = 'turquoise' ),
                          AS = dict( color = 'turquoise' ),
                                )



def get_plot( n_subplots, axs ):
	if axs == None:
		return plt.subplots( n_subplots, figsize = (8, 4/3 * n_subplots ), sharex = True, gridspec_kw = {'hspace': 0} )
	else:
		try:
			if len( axs ) != n_subplots:
				raise ValueError( 'Length of passed matplotlib.pyplot.axes ({}) list does not equal the number of parameters to plot ({}).'.format( len( axs ), n_subplots ) )
		except Exception:
			pass
		return [ None, axs ]


def get_plot_limits( data, offset = 0.1 ):
	_min = np.min( data )
	_max = np.max( data )
	if _min == _max:
		return [ _min - 0.75, _max + 0.75 ]
	else:
		return [ _min - offset * np.abs( _max - _min ), _max + offset * np.abs( _max - _min ) ]

def finalize_plot( figure = None, axs = None, out_file_path = None, overwrite = False, show = True ):
	if ( out_file_path != None ) or show:
		if ( figure != None ) and is_iterable( axs ):
			figure.align_ylabels( axs[:] )
		plt.tight_layout()
		if out_file_path != None:
			out_file_path = make_output_path( None, out_file_path, overwrite )
			plt.savefig( out_file_path )
		if show:
			plt.show()
	return