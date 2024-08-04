


from target_approximation.utils import finalize_plot
from target_approximation.utils import get_plot
from target_approximation.utils import get_plot_limits
import numpy as np
import matplotlib.pyplot as plt



class TubeState():
    def __init__(
            self,
            tract_state,
            tube_length,
            tube_area,
            tube_articulator,
            incisor_position,
            tongue_tip_side_elevation,
            velum_opening,
            ):
        self.tract_state = tract_state
        self.tube_length = tube_length
        self.tube_area = tube_area
        self.tube_articulator = tube_articulator
        self.incisor_position = incisor_position
        self.tongue_tip_side_elevation = tongue_tip_side_elevation
        self.velum_opening = velum_opening
        self.open_limit = 0.3  # 0.3 cm^2 for open tracts
        self.tight_limit = 0.001 # above 0.001 tight, below or equal closed # actual value is 0.0001 however
        self.constriction = self.get_constriction_class(
            tube_area = np.min( tube_area ),
            )
        self.constriction_data = self.get_constriction_threshold_crossings()
        return
    
    @classmethod
    def from_dict(
            cls,
            x,
            ):
        return cls(
            tract_state = x[ 'tract_state' ],
            tube_length = x[ 'tube_length' ],
            tube_area = x[ 'tube_area' ],
            tube_articulator = x[ 'tube_articulator' ],
            incisor_position = x[ 'incisor_position' ],
            tongue_tip_side_elevation = x[ 'tongue_tip_side_elevation' ],
            velum_opening = x[ 'velum_opening' ],
            )

    #def get_constriction( self, return_str = False ):
    #    constriction_strings = [ 'open', 'tight', 'closed' ]
    #    min_area = np.min( self.tube_area )
    #    constriction = None
    #    if min_area >= self.open_limit:
    #        constriction = 0
    #    elif np.isclose( min_area, 0.15 ):
    #        constriction = 3
    #    elif np.isclose( min_area, 0.25 ):
    #        constriction = 4
    #    elif min_area > self.tight_limit:
    #        constriction = 1
    #    elif np.isclose( min_area, 0.0001 ):
    #        constriction = 2
    #    if not return_str:
    #        return constriction
    #    else:
    #        return constriction_strings[ constriction ]

    def get_constriction_class(
        self,
        tube_area,
        ):
        constriction = None
        if tube_area >= self.open_limit:
            constriction = 0
        elif np.isclose( tube_area, 0.15 ):
            constriction = 3
        elif np.isclose( tube_area, 0.25 ):
            constriction = 4
        elif tube_area > self.tight_limit:
            constriction = 1
        elif np.isclose( tube_area, 0.0001 ):
            constriction = 2
        elif tube_area <= self.tight_limit:
            constriction = 5
        return constriction

    def get_tube_area_function( self ):
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
        return np.array( [ x, y ] ).T

    def get_constriction_threshold_crossings(
        self,
        n_tongue_sections = 8,
        ):
        tight_crossings = []
        close_crossings = []
        tight_crossed = False
        close_crossed = False
        tight_articulators = []
        close_articulators = []
        #for x in tube_area_function:
        #    y = x[1]
        #    if tight_crossed == False and y < self.open_limit:
        #        tight_crossings.append( x )
        #        tight_crossed = True
        #    if tight_crossed == True and y >= self.open_limit:
        #        tight_crossings.append( x )
        #        tight_crossed = False
        #    if close_crossed == False and y < self.tight_limit:
        #        close_crossings.append( x )
        #        close_crossed = True
        #    if close_crossed == True and y >= self.tight_limit:
        #        close_crossings.append( x )
        #        close_crossed = False
        articulator_token = {
            '1': 'T',# tongue;
            '2': 'I',#= lower incisors;
            '3': 'L',
            '4': 'O',# = lower lip; 4 = other
        }
        n_tongue = len( [ ar for ar in self.tube_articulator if ar == 1 ] )
        tube_articulator_tokens = []
        tongue_counter = 0
        tongue_section = -1
        for ar in self.tube_articulator:
            #print(ar)
            if ar == 1:
                if ( tongue_counter % round(n_tongue/n_tongue_sections) == 0 ) and ( tongue_section < (n_tongue_sections-1) ):
                    tongue_section += 1
                tongue_counter += 1
                tube_articulator_tokens.append(
                    'T{}'.format( tongue_section )
                    )
            else:
                tube_articulator_tokens.append(
                    articulator_token[ str( ar ) ]
                    )
        self.tube_articulator_tokens = tube_articulator_tokens

        assert len( self.tube_area ) == len( self.tube_length ), 'Not the same length, ta: {}, tl: {}'.format(
            len( self.tube_area ),
            len( self.tube_length ),
            )
        assert len( self.tube_area ) == len( self.tube_articulator_tokens ), 'Not the same length, ta: {}, ar: {}'.format(
            len( self.tube_area ),
            len( self.tube_articulator_tokens ),
            )

        x = 0
        for ta, tl, ar in zip( self.tube_area, self.tube_length, self.tube_articulator_tokens ):
            if tight_crossed == False and ta < self.open_limit:
                tight_crossings.append( x )
                tight_tb_articulators = []
                tight_tb_articulators.append(
                    dict(
                        start = x,
                        place_of_articulation = ar,
                        tube_area = ta,
                        constriction_class = self.get_constriction_class( ta ),
                        )
                    )
                #tight_articulators.append( [ x, ar ] )
                tight_crossed = True
            elif tight_crossed == True and ta < self.open_limit:
                #tight_articulators.append( [ x, ar ] )
                tight_tb_articulators.append(
                    dict(
                        start = x,
                        place_of_articulation = ar,
                        tube_area = ta,
                        constriction_class = self.get_constriction_class( ta ),
                        )
                    )
            elif tight_crossed == True and ta >= self.open_limit:
                tight_crossings.append( x )
                #tight_articulators.append( ar )
                tight_articulators.append( tight_tb_articulators )
                tight_crossed = False
            if close_crossed == False and ta < self.tight_limit:
                close_crossings.append( x )
                #close_articulators.append( [ x, ar ] )
                close_tb_articulators = []
                close_tb_articulators.append(
                    dict(
                        start = x,
                        place_of_articulation = ar,
                        tube_area = ta,
                        constriction_class = self.get_constriction_class( ta ),
                        )
                    )
                close_crossed = True
            elif close_crossed == True and ta < self.tight_limit:
                close_tb_articulators.append(
                    dict(
                        start = x,
                        place_of_articulation = ar,
                        tube_area = ta,
                        constriction_class = self.get_constriction_class( ta ),
                        )
                    )
            elif close_crossed == True and ta >= self.tight_limit:
                close_articulators.append( close_tb_articulators )
                close_crossings.append( x )
                close_crossed = False
            x += tl
        if tight_crossed == True:
            tight_crossings.append( x )
            tight_articulators.append( tight_tb_articulators )
        if close_crossed == True:
            close_crossings.append( x )
            close_articulators.append( close_tb_articulators )

        tight_constrictions = self.get_constriction_info( tight_crossings, tight_articulators )
        close_constrictions = self.get_constriction_info( close_crossings, close_articulators )
        constriction_data = dict(
            n_constrictions = len( tight_constrictions ) + len( close_constrictions ) - (1 if len(close_constrictions) > 0 else 0),
            tight_constrictions = tight_constrictions,
            close_constrictions = close_constrictions,
            )
        return constriction_data

    def get_constriction_info( self, threshold_crossings, articulators ):
        constrictions = []
        index = 0
        ar_id = 0
        while index < len( threshold_crossings ) - 1:
            #print( articulators )
            constrictions.append( 
                dict(
                    start = threshold_crossings[ index ],
                    end = threshold_crossings[ index + 1 ],
                    length = threshold_crossings[ index + 1 ] - threshold_crossings[ index ],
                    articulators = articulators[ ar_id ],
                    #area = None,
                    )
                )
            index += 2
            ar_id += 1
        #if len( threshold_crossings ) % 2 != 0:
        #    start = threshold_crossings[ len( threshold_crossings ) - 1 ]
        #    end = np.sum( self.tube_length )
        #    print( articulators )
        #    constrictions.append( 
        #        dict(
        #            start = start,
        #            end = end,
        #            length = end-start,
        #            articulators = articulators[ -1 ],
        #            #area = 0,
        #            )
        #        )
        return constrictions

    def has_precise_constriction( self ):
        tube_area_function = self.get_tube_area_function()
        tight_crossings, close_crossings = self.get_constriction_threshold_crossings( tube_area_function )
        threshold_crossings = [ len( tight_crossings ), len( close_crossings ) ]
        if self.constriction == 2:
            if not threshold_crossings in [ [1,1], [2,2] ]:
                return False
        elif self.constriction == 1:
            if not threshold_crossings in [ [2,0], [1,0] ]:
                return False
        elif self.constriction == 0:
            return False
        #if self.constriction_has_local_minimum( x, y, ):
        #    return False
        if threshold_crossings == [ 1, 0 ]:
            if tight_crossings[0][0] <= ( 1 - 0.125 ) * np.max( tube_area_function[ :, 0 ] ):
                return False
        elif threshold_crossings == [ 1, 1 ]:
            if np.abs( (np.max( tube_area_function[ :, 0 ] ) - close_crossings[0][0]) - (close_crossings[0][0] - tight_crossings[0][0]) ) >= 1:
                return False
        elif threshold_crossings == [ 2, 2 ]:
            minimum = 0.5 * ( close_crossings[0][0] + close_crossings[1][0] )
            if np.abs(  np.abs( minimum - tight_crossings[0][0] ) - np.abs( tight_crossings[1][0] - minimum ) ) >= 1:
                return False


        return True

    def plot( self, 
              axs = None, 
              **kwargs,
              ):
        #articulators = {
        #    '1': 'T',# tongue;
        #    '2': 'I',#= lower incisors;
        #    '3': 'L',
        #    '4': 'O',# = lower lip; 4 = other
        #}
        figure, axs = get_plot( n_rows = 1, axs = axs )
        tube_area_function = self.get_tube_area_function()
        axs[0].set( xlabel = 'Tube Length [cm]', ylabel = r'Cross-sectional Area [cm$^2$]' )
        #y = [ val for val in x  ]
        #x = [ self.tube_length[ 0 ] ]
        #for length in self.tube_length[ 1: ]:
        #    x.append( x[ -1 ] + length )
        axs[0].plot( tube_area_function[ :, 0 ], tube_area_function[ :, 1 ] )
        constriction_data = self.get_constriction_threshold_crossings()
        tight_constrictions = constriction_data[ 'tight_constrictions' ]
        close_constrictions = constriction_data[ 'close_constrictions' ]
        for tight_constriction in tight_constrictions:
            #axs[0].scatter( tight_crossing[0], tight_crossing[1], color = 'red', marker = 'x' )
            axs[0].scatter( tight_constriction[ 'start' ], 0.3, color = 'red', marker = 'x' )
            axs[0].scatter( tight_constriction[ 'end' ], 0.3, color = 'red', marker = 'x' )
            axs[0].plot( [tight_constriction[ 'start' ], tight_constriction[ 'start' ] + tight_constriction[ 'length' ] ], [ 0.9 , 0.9 ] )
            for element in tight_constriction[ 'articulators' ]:
                axs[0].text( s=element[ 'place_of_articulation' ], x=element[ 'start' ], y = 0.35 )
        for close_constriction in close_constrictions:
            #axs[0].scatter( close_crossing[0], close_crossing[1], color = 'green', marker = 'o' )
            #axs[0].scatter( close_crossing[0], 0.001, color = 'green', marker = 'o' )
            axs[0].scatter( close_constriction[ 'start' ], 0.001, color = 'green', marker = 'o' )
            axs[0].scatter( close_constriction[ 'end' ], 0.001, color = 'green', marker = 'o' )
            axs[0].plot( [close_constriction[ 'start' ], close_constriction[ 'start' ] + close_constriction[ 'length' ] ], [ 0.01 , 0.01 ] )
        axs[0].axhline( 0.3, color = 'gray', ls = '--' )
        axs[0].axhline( 0.001, color = 'gray', ls = '-.' )
        axs[0].set( yscale = 'log' )
        finalize_plot( figure, axs, **kwargs )
        #ax.set( xlabel = 'Tube Length [cm]', ylabel = r'Cross-sectional Area [cm$^2$]' )
        return axs