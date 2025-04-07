

from copy import deepcopy
import matplotlib.pyplot as plt
import numpy as np
from tools_mp import process
#from target_approximation.vocaltractlab import SupraGlottalSeries as SGS
import uuid
import yaml


from .sample import get_hull, get_formant_data
from .sample import sample_c
from .sample import sample_v
from .vtl_core import ParameterSpace
from .vtl_core import get_valid_parameter_space
from .. import core as vtl


from typing import Optional
from typing import Union
from typing import List
from typing import Dict


class State():
    def __init__(
            self,
            x: np.ndarray,
            poa: str,
            parent: str,
            id: str,
            group: str,
            ):
        self.values = np.array( x ).squeeze()
        self.poa = poa
        self.parent = parent
        self.id = id
        self.group = group
        return
    
    def __str__(self) -> str:
        s = f'{self.poa} {self.group}'
        return s
    
    @classmethod
    def from_dict(
            cls,
            x: dict,
            ):
        s = cls(
            x = x[ 'x' ],
            poa = x[ 'poa' ],
            parent = x[ 'parent' ],
            group = x[ 'group' ],
            )
        return s
    
    def to_dict(self):
        x = dict(
            x = self.values.tolist(),
            poa = self.poa,
            parent = self.parent,
            #id = self.id,
            group = self.group,
            )
        return x

class PhonemeInventory():
    def __init__(
            self,
            parameter_space: ParameterSpace,
            states: dict = {},
            ):
        self.parameter_space = parameter_space
        self.states = states
        return
    
    def __str__(self) -> str:
        return str(self.items)
    
    @classmethod
    def create(
            cls,
            parameter_space: Optional[ ParameterSpace ] = None,
            speaker: Optional[str] = None,
            v_kwargs: Optional[ dict ] = None,
            c_kwargs: Optional[ dict ] = None,
            ):

        # Handle parameter search space
        if (parameter_space is None and speaker is None) or \
           (parameter_space is not None and speaker is not None):
            raise ValueError(
                f"""
                You must pass either a speaker name XOR a
                parameter_space object.
                """
                )
        elif parameter_space is None and speaker is not None:
            parameter_space = get_valid_parameter_space(
                speaker = speaker,
            )

        v_states = cls._sample_vowels(
            parameter_space,
            group = 'initial',
            **v_kwargs,
            )
        c_states = cls._sample_consonants(
            parameter_space,
            v_states,
            group = 'initial',
            **c_kwargs,
            )

        states = cls._merge( v_states, c_states )
        
        return cls(
            parameter_space = parameter_space,
            states = states,
            )
    
    @classmethod
    def from_yaml(
            cls,
            file_path: str,
            ):
        with open( file_path, 'r' ) as f:
            x = yaml.load( f, Loader = yaml.FullLoader )
        ps = ParameterSpace.from_dict( x[ 'parameter_space' ] )
        s = { k: State( **v ) for k, v in x[ 'states' ].items() }
        return cls(
            parameter_space = ps,
            states = s,
            )
    
    def _sample_consonants(
            self,
            parameter_space: ParameterSpace,
            v_states: Dict[str, State],
            grid_poa = [
                [ 'T2', 'T3' ],
                [ 'T3', 'T4' ],
                [ 'T4', 'T5' ],
                [ 'T5', 'T6' ],
                [ 'T6', 'T7' ],
                [ 'T3' ],
                [ 'T4' ],
                [ 'T5' ],
                [ 'T6' ],
                [ 'T7' ],
                [ 'LF' ],
                [ 'L' ],
            ],
            group: Optional[str] = None,
            **kwargs,
            ):
        
        args = []
        indices = []
        poa_list = []
        for idx, v in v_states.items():
            for poa in grid_poa:
                x = dict(
                    shape = v.values,
                    poa = poa,
                    parameter_space = parameter_space,
                    **kwargs,
                    )
                args.append(x)
                indices.append(idx)
                poa_list.append(poa)
                    
        results = process(
            function = sample_c,
            args = args,
            return_data = True,
            verbose = True,
            workers = None,
            mp_threshold = 4,
            # Not necessary because sample_c already loads the speaker
            initializer = vtl.load_speaker,
            initargs = ( vtl.active_speaker(), ),
            )
        print( 'Length of indices:', len( indices ) )
        print( 'Length of results:', len( results ) )
        
        inventory = {}
        for i, poa, result in zip(indices, poa_list, results):
            #poa = '_'.join( result[ 'place_of_articulation' ] ) # does not consider LF
            poa = '_'.join( poa )
            result_state = np.array(
                [ v for _, v in result[ 'result_state' ].items() ]
                ).squeeze()

            k = str(uuid.uuid4())
            s = State(
                x = result_state,
                poa = poa,
                parent = i,
                id = k,
                group = group,
                )
            inventory[k] = s

        return inventory

    def _sample_vowels(
            self,
            parameter_space: ParameterSpace,
            group: Optional[str] = None,
            **kwargs,
            ):
        # if n_hull_layers is in kwargs, get it and remove it
        n_hull_layers = kwargs.pop( 'n_hull_layers', None )

        x, _ = sample_v( parameter_space, **kwargs )
        
        if n_hull_layers is not None:
            _, hull_states, _, _ = get_hull(
                states = x,
                n_hull_layers = n_hull_layers,
                )
            x = hull_states

        #states = [
        #    State(
        #        x = v,
        #        poa = 'vowel',
        #        parent = None,
        #        group = group,
        #        )
        #    for v in x
        #    ]

        inventory = {}
        for v in x:
            k = str(uuid.uuid4())
            s = State(
                x = v,
                poa = 'vowel',
                parent = None,
                id = k,
                group = group,
                )
            inventory[k] = s

        return inventory
    
    def _merge(
            self,
            *args,
            ):
        x = deepcopy( args[0] )
        for arg in args[1:]:
            # if no keys overlap, update the inventory
            if len( set( x.keys() ).intersection( set( arg.keys() ) ) ) == 0:
                x.update( arg )
            else:
                raise ValueError(
                    f"""
                    Overlapping keys found
                    {set( x.keys() ).intersection( set( arg.keys() ) )}
                    """
                    )
        return x

    def _update(
            self,
            *args,
            ):
        self.states = self._merge( self.states, *args )
        return
        
    def add_vowels(
            self,
            group: Optional[str] = None,
            **kwargs,
            ):
        v_states = self._sample_vowels(
            parameter_space = self.parameter_space,
            group = group,
            **kwargs,
            )
        self._update( v_states )
        return
    
    def add_consonants(
            self,
            group: Optional[str] = None,
            **kwargs,
            ):
        c_states = self._sample_consonants(
            parameter_space = self.parameter_space,
            v_states = self.vowels(),
            group = group,
            **kwargs,
            )
        self._update( c_states )
        return
    
    def get_children(
            self,
            x: Union[str, State],
            ):
        if isinstance( x, str ):
            idx = x
        elif isinstance( x, State ):
            idx = x.id
        children = {
            k: v
            for k, v in self.states.items()
            if v.parent == idx
            }
        return children
    
    def get_parent(
            self,
            x: Union[str, State],
            ):
        if isinstance( x, str ):
            idx = self.states[ x ].parent
        elif isinstance( x, State ):
            idx = x.parent
        s = self.states[ idx ]
        return s

    #def get_vowel(
    #        self,
    #        x: Union[str, State],
    #        ):
    #    if isinstance( x, str ):
    #        s = self.states[ x ]
    #    elif isinstance( x, State ):
    #        s = x
    #    if s.poa != 'vowel':
    #        parent = self.get_parent( s )
    #        if parent.poa == 'vowel':
    #            s = parent
    #        else:
    #            raise ValueError(
    #                f"""
    #                Neither the input state nor its parent are vowels.
    #                """
    #                )
    #    return s
    
    def get_vcv(
            self,
            x: Union[str, State],
            poa: str,
            ):
        if isinstance( x, str ):
            s = self.states[ x ]
        elif isinstance( x, State ):
            s = x
        if s.poa == 'vowel':
            v = s
            consonants = [
                child
                for child in self.get_children( v )
                if child.poa != 'vowel'
                ]
            if len( consonants ) == 0:
                raise ValueError(
                    f"""
                    No consonants found for vowel {v}.
                    """
                    )
            else:
                # draw a random consonant
                c = np.random.choice( consonants )
        else:
            c = s
            v = self.get_parent( c )
        return v, c

        
    
    def get_rnd(
            self,
            n_samples: int,
            sequence: Optional[ List[ str ] ] = None,
            ):
        return

    def get_speaker(self):
        return self.parameter_space.speaker
    
    def plot_vowels(self, ax = None):
        vowels = self.vowels( return_type = 'numpy' )
        fmt_vowels = get_formant_data( vowels )
        if ax is None:
            fig, ax = plt.subplots()
        ax.scatter( fmt_vowels[:,0], fmt_vowels[:,1] )
        return ax
    
    def phonemes(self, poa, return_type = 'states'):
        x = {
            k: v
            for k, v in self.states.items()
            if v.poa ==  poa
            }
        if return_type == 'numpy':
            x = np.array( [
                v.values
                for _, v in x.items()
                ]
            )
        return x
        
    def vowels(self, **kwargs):
        x = self.phonemes( 'vowel', **kwargs )
        return x
    
    def to_dict(self):
        x = dict(
            parameter_space = self.parameter_space.to_dict(),
            states = { k: v.to_dict() for k, v in self.states.items() },
            )
        return x
    
    def to_yaml(
            self,
            file_path: str,
            ):
        x = self.to_dict()
        with open( file_path, 'w' ) as f:
            yaml.dump( x, f, sort_keys=False )
        return
        

    

        

        
    


