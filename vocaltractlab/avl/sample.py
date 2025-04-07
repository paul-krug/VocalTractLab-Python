

from copy import deepcopy
import numpy as np
from scipy.spatial import ConvexHull
from target_approximation.vocaltractlab import SupraGlottalSeries
#from pyMetaheuristic.algorithm import whale_optimization_algorithm

from .woa import whale_optimization_algorithm
from .sensory_optimization import Agent
from .vtl_core import ParameterSpace
from .. import core as vtl

from typing import Optional
from typing import Union
from typing import List


def get_formant_data( states ):
    sgs = SupraGlottalSeries( states )
    trf = vtl.motor_to_transfer_function( sgs )
    formant_data = np.array( [
        [ x.f1 for x in trf ],
        [x.f2 for x in trf ],
        ] ).T
    return formant_data

def get_hull(
        states,
        n_hull_layers = 5,
        ):

    #sgs = SupraGlottalSeries( states )
    #trf = vtl.motor_to_transfer_function( sgs )
    #formant_data = np.array( [ [ x.f1 for x in trf ], [x.f2 for x in trf ] ] ).T
    formant_data = get_formant_data( states )

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

def sample_c(
    shape: Union[ np.ndarray, str ],
    poa: List[ str],
    parameter_space: ParameterSpace,
    parameter_blacklist: List[ str ],
    steps: int = 5000,
    tight_constriction_length: int = 3,
    close_constriction_length: int = 3,
    extend_space = 0.1,
    ):

    # Handle parameter search space
    vtl.load_speaker(
        speaker = parameter_space.speaker,
        auto_tongue_root = parameter_space.auto_tongue_root,
        )

    optimization_parameters = [
        x for x in parameter_space.get_parameter_names()
        if x not in parameter_blacklist
        ]
    min_state, max_state = parameter_space.get_extended_range(
        extend_below = extend_space,
        extend_above = extend_space,
        parameters = optimization_parameters,
        )

    # Handle constriction constraints
    if poa == ['LF']:
        constriction = 3
        poa = [ 'L' ]
    else:
        constriction = 2

    constriction_information = dict(
        constriction = constriction,
        n_constrictions = 1,
        tight_constriction_length = tight_constriction_length,
        close_constriction_length = close_constriction_length,
        place_of_articulation = poa,
        )
    
    # Handle optimization arguments
    optimization_policy = whale_optimization_algorithm
    optimization_policy_kw = dict(
        hunting_party = 200,
        iterations = 500000,
        spiral_param = 0.5,
        verbose = False,
        min_values = min_state,
        max_values = max_state,
        )

    # Get the reference state
    if isinstance( shape, str ):
        try:
            reference_state = vtl.get_shape( shape, 'tract' )
        except ValueError:
            raise ValueError(
                f'Specified shape: {shape} not available in speaker file.'
                )
    else:
        reference_state = np.array( shape ).reshape( ( 19 ) )

    # Define optimization agent and run
    agent = Agent(
        reference_state = reference_state,
        optimization_parameters = optimization_parameters,
        optimization_policy = optimization_policy,
        optimization_policy_kw = optimization_policy_kw,
        log_dir = None,
        max_steps = steps,
        verbose = False,
        **constriction_information,
        )
    
    result = agent.run()
    return result

def sample_v(
    parameter_space: ParameterSpace,
    n_samples: int = 100,
    verbose: bool = True,
    ):

    # Handle parameter search space
    vtl.load_speaker(
        speaker = parameter_space.speaker,
        auto_tongue_root = parameter_space.auto_tongue_root,
        )
    min_state, max_state = parameter_space.get_range()

    # Get the probability of open state (0):
    p_open = parameter_space.occurences[ '0' ]
    # Determine the most efficient batch size
    # assuming bad case scenario: 3 std deviation
    p_open_worst = p_open[ 'mean' ] - 3 * p_open[ 'std' ]
    sampling_batch_size = int( n_samples / p_open_worst )


    # Sample states in batches, also count constriction occurences
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

        tbs = vtl.motor_to_tube( sgs , verbose=verbose )

        for _, x in enumerate( tbs ):
            if x.constriction == 0:
                states.append( x.tract_state )
            constriction_occurences[ str( x.constriction ) ] += 1
        
    states = np.array( states )[ : n_samples ]
    return states, constriction_occurences