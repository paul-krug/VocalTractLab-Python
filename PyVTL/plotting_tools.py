import numpy as np
import matplotlib.pyplot as plt






def get_plot( n_subplots ):
	return plt.subplots( n_subplots, figsize = (8, 4/3 * n_subplots ), sharex = True, gridspec_kw = {'hspace': 0} )


def get_plot_limits( data, offset = 0.1 ):
	_min = np.min( data )
	_max = np.max( data )
	if _min == _max:
		return [ _min - offset * _min, _max + offset * _max ]
	else:
			return [ _min - offset * np.abs( _max - _min ), _max + offset * np.abs( _max - _min ) ]