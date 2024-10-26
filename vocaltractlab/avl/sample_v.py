

from copy import deepcopy
import numpy as np
from scipy.spatial import ConvexHull
from target_approximation.vocaltractlab import SupraGlottalSeries

import vocaltractlab as vtl
from ..vtl_core import ParameterSpace

def sample_v(
    parameter_space: ParameterSpace,
    n_samples = 100,
    verbose = True,
    ):

    # Handle parameter search space
    speaker = parameter_space.speaker
    vtl.load_speaker( speaker )
    min_state, max_state = parameter_space.get_range()

    # Determine the most efficient batch size
    p_open = parameter_space.occurences[ '0' ]
    # assume bad case scenario: 3 std deviation
    p_open_worst = p_open[ 'mean' ] - 3 * p_open[ 'std' ]
    sampling_batch_size = int( n_samples / p_open_worst )


    states = []
    constriction_occurences = {
        '0': 0,
        '1': 0,
        '2': 0,
        '3': 0,
        '4': 0,
        '5': 0,
        }
    while len( states ) < n_samples:
        candidates = np.random.uniform(
            low = min_state,
            high = max_state,
            size = ( sampling_batch_size, min_state.shape[0] ),
            )
        sgs = SupraGlottalSeries( candidates )
        sgs[ 'VO' ] = -0.1
        #axs = sgs.plot_distributions()

        tbs = vtl.motor_to_tube( sgs , verbose=verbose )

        for _, x in enumerate( tbs ):
            if x.constriction == 0:
                states.append( x.tract_state )
            constriction_occurences[ str( x.constriction ) ] += 1
        
    states = np.array( states )[ : n_samples ]
    return states, constriction_occurences

def get_hull(
        states,
        n_hull_layers = 5,
        ):

    sgs = SupraGlottalSeries( states )

    trf = vtl.motor_to_transfer_function( sgs )
    formant_data = np.array( [ [ x.f1 for x in trf ], [x.f2 for x in trf ] ] ).T

    fmt = deepcopy( formant_data )
    sts = deepcopy( states )
    #print( formant_data.shape )
    hull_layers = []
    hull_states = []
    hull_list = []
    for layer_index in range(0, n_hull_layers):
        hull = ConvexHull( fmt )
        hull_indices = sorted( list( set( hull.simplices.flatten() ) ) )
        hull_layers.append( fmt[ hull_indices ] )
        hull_states.append( sts[ hull_indices ] )
        fmt = np.delete( fmt, hull_indices, axis = 0 )
        sts = np.delete( sts, hull_indices, axis = 0 )
        hull_list.append( hull )

    hull_states = np.concatenate( hull_states )
    return hull_list, hull_states, hull_layers, formant_data
    
def get_formant_data( states ):
    sgs = SupraGlottalSeries( states )
    trf = vtl.motor_to_transfer_function( sgs )
    formant_data = np.array( [
        [ x.f1 for x in trf ],
        [x.f2 for x in trf ],
        ] ).T
    return formant_data