
import numpy as np
import pandas as pd
import vocaltractlab as vtl
from target_approximation.vocaltractlab import SupraGlottalSeries
import yaml

from typing import Optional
from typing import List


def get_occurences(
        min_state: np.ndarray,
        max_state: np.ndarray,
        n_samples: int = 1000,
        n_repeats: int = 10,
        verbose: bool = True,
        ):
    print( 'Calculating occurences' )

    candidates = np.random.uniform(
        low = min_state,
        high = max_state,
        size = (
            int( n_repeats*n_samples ),
            min_state.shape[0],
            ),
        )
    sgs = SupraGlottalSeries( candidates )
    sgs[ 'VO' ] = -0.1
    tbs = vtl.motor_to_tube( sgs , verbose=verbose )

    # chunk tbs into n_repeats * n_samples
    runs = []
    for i in range( n_repeats ):
        occ = np.zeros( 6 )
        chunk = tbs[ i*n_samples : (i+1)*n_samples ]
        for t in chunk:
            occ[ t.constriction ] += 1
        runs.append( occ )

    runs = np.array( runs )
    #print( runs.shape )
        
    constriction_occurences = {
        str(k): {
            'mean': float( runs[ :, int( k ) ].mean()/n_samples ),
            'std': float( runs[ :, int( k ) ].std()/n_samples ),
            #'n_samples': runs[ :, int( k ) ].shape[0],
            }
        for k in range( runs.shape[1] )
    }
    return constriction_occurences

def get_parameter_space(
        speaker,
        auto_tongue_root: bool,
        ):

    vtl.load_speaker( speaker )
    vtl.calculate_tongueroot_automatically( auto_tongue_root )

    info = vtl.get_param_info( 'tract' )
    info = pd.DataFrame( info )
    #print(info)
    min_state = info.loc[ :, 'min' ].to_numpy( dtype = float )
    max_state = info.loc[ :, 'max' ].to_numpy( dtype = float )

    min_state = info[ 'min' ]
    max_state = info[ 'max' ]

    parameter_names = info.loc[ :, 'name' ].to_list()

    search_space = pd.DataFrame(
        dict(
            #name = parameter_names,
            min = min_state,
            max = max_state,
            ),
        )
    search_space.index = parameter_names

    occurences = get_occurences(
        min_state = min_state,
        max_state = max_state,
        n_samples = 1000,
        n_repeats = 10,
        verbose = True,
        )

    ps = ParameterSpace(
        speaker = speaker,
        space = search_space,
        occurences = occurences,
        auto_tongue_root = auto_tongue_root,
        )
    return ps

def get_valid_parameter_space(
        speaker,
        auto_tongue_root: bool,
        sample_size = 100000,
        extend_below = 2.0,
        extend_above = 0.0,
        ):
    # speaker is set in get_parameter_space
    ps = get_parameter_space(
        speaker = speaker,
        auto_tongue_root = auto_tongue_root,
        )

    ext_range_min, ext_range_max = ps.get_extended_range(
        extend_below = extend_below,
        extend_above = extend_above,
        )

    candidates = np.random.uniform(
        low  = ext_range_min,
        high = ext_range_max,
        size = ( sample_size, ext_range_min.shape[0] ),
        )
    sgs = SupraGlottalSeries( candidates )
    #axs = sgs.plot_distributions( show = False )
    sgs_lim = vtl.limit( sgs )
    #sgs_lim.plot_distributions( axs = axs )
    #plt.show()
    #stop

    # TODO: clip "LD" to [ -0.5, 2.0 ] if desired

    min_state = np.min( sgs_lim.to_numpy( transpose = False ), axis = 0 )
    max_state = np.max( sgs_lim.to_numpy( transpose = False ), axis = 0 )

    print( 'Calculating occurences' )
    occurences = get_occurences(
        min_state = min_state,
        max_state = max_state,
        n_samples = 1000,
        n_repeats = 10,
        verbose = True,
        )

    search_space = dict(
        min = min_state,
        max = max_state,
    )
    ps = ParameterSpace(
        speaker = speaker,
        space = pd.DataFrame(
            search_space,
            index = ps.get_parameter_names(),
            ),
        occurences = occurences,
        auto_tongue_root=auto_tongue_root,
        )
    return ps

class ParameterSpace():
    def __init__(
            self,
            speaker: str,
            space: pd.DataFrame,
            occurences: dict,
            auto_tongue_root: Optional[ bool ] = None,
            ):
        self.speaker = speaker
        self.space = space
        self.occurences = occurences
        self.auto_tongue_root = auto_tongue_root
        return
    
    def __str__(self) -> str:
        _str_0 = f'Speaker: {self.speaker }\n'
        _str_1 = f'\n{self.space}'
        return _str_0 + _str_1
    
    @classmethod
    def from_dict(
            cls,
            x: dict,
            ):
        speaker = x[ 'speaker' ]
        space = pd.DataFrame( x[ 'space' ] )
        occurences = x[ 'occurences' ]
        auto_tongue_root = x[ 'auto_tongue_root' ]
        return cls(
            speaker = speaker,
            space = space,
            occurences = occurences,
            auto_tongue_root = auto_tongue_root,
            )
    
    @classmethod
    def from_yaml(
            cls,
            file_path: str,
            ):
        with open( file_path, 'r' ) as f:
            x = yaml.load( f, Loader = yaml.FullLoader )
        return cls.from_dict( x )
    
    def get_parameter_names(
            self,
            ):
        return self.space.index.to_list()
    
    def get_range(
            self,
            parameters: Optional[ List[ str ] ] = None,
            ):
        if parameters is None:
            parameters = self.get_parameter_names()
        min_state = self.space.loc[ parameters, 'min' ].to_numpy( dtype = float )
        max_state = self.space.loc[ parameters, 'max' ].to_numpy( dtype = float )
        return min_state, max_state
    
    def get_extended_range(
            self,
            extend_below = 2.0,
            extend_above = 0.0,
            parameters: Optional[ List[ str ] ] = None,
            ):
        range_min, range_max = self.get_range( parameters = parameters )
        range_nrm = range_max - range_min
        ext_range_min = range_min - extend_below * range_nrm
        ext_range_max = range_max + extend_above * range_nrm
        return ext_range_min, ext_range_max
    
    def to_dict(
            self,
            ):
        x = dict(
            speaker = self.speaker,
            space = self.space.to_dict(),
            occurences = self.occurences,
            auto_tongue_root = self.auto_tongue_root,
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
    

class VTL():
    def __init__(
            self,
            speaker: str,
            auto_tongue_root: bool,
            ):
        self.speaker = speaker
        self.auto_tongue_root = auto_tongue_root
        return