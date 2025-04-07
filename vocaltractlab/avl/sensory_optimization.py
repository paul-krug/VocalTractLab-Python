



#import vocaltractlab as vtl
from target_approximation.vocaltractlab import SupraGlottalSeries as SGS
import numpy as np
import os


#from sklearn.metrics import mean_squared_error
from scipy.spatial.distance import cosine
import time
#import matplotlib.pyplot as plt
#import argparse
from copy import deepcopy
#from itertools import chain
#from tools_io import save, load


#from random import getrandbits
#from tools_mp import multiprocess

from .. import core as vtl

from typing import Optional
from typing import List

def get_mask(
    query: SGS,
    opt_params: List[ str ],
    ) -> np.ndarray:
    mask = SGS(
        np.ones_like(
            query.to_numpy(transpose=False),
            dtype = bool,
            ),
        )
    for p in opt_params:
        mask[ p ] = False
    mask = mask.to_numpy()#.squeeze()
    return mask
    
def state_to_dict(
        state: SGS,
        ) -> dict:
    x = {
        k: v[0] for k,v in state.to_dict()[ 'series' ].items()
        }
    return x



class EarlyStoppingException( Exception ): pass

class Agent():
    def __init__(
        self,
        reference_state: np.ndarray,
        optimization_parameters: List[ str ],
        optimization_policy,
        optimization_policy_kw,
        constriction,
        n_constrictions,
        tight_constriction_length,
        close_constriction_length,
        place_of_articulation,
        log_dir: str,
        max_steps,
        verbose = True,
        ):

        # Reference state
        self.reference_state = reference_state
        self.sgs_r = SGS( reference_state )

        # Parameter space
        self.optimization_parameters = optimization_parameters
        n_opt_params = len( optimization_parameters )
        n_min_state = len( optimization_policy_kw[ 'min_values' ] )
        n_max_state = len( optimization_policy_kw[ 'max_values' ] )
        # if those lenths not all the same throw error
        if not all( [ n_opt_params == n_min_state, n_opt_params == n_max_state ] ):
            raise ValueError(
                f"""
                The number of optimization parameters ({n_opt_params}),
                the number of minimum values ({n_min_state}),
                and the number of maximum values ({n_max_state}) must be equal.
                """
                )

        # Optimizatiob constraints
        self.constriction = constriction
        self.n_constrictions = n_constrictions
        self.tight_constriction_length = tight_constriction_length
        self.close_constriction_length = close_constriction_length
        self.place_of_articulation = place_of_articulation

        # Metaheuristic algorithm and its parameters
        self.optimization_policy = optimization_policy
        self.optimization_policy_kw = deepcopy( optimization_policy_kw )

        # Results and itermediate results
        self.log_category_data = []
        self.log_loss_data = [
            dict(
                step = -1,
                total_loss = np.inf,
                constriction_loss = np.inf,
                similarity_loss = np.inf,
                #reference_state = None,
                )
            ]
        
        # Current step and maximum number of steps
        self.step = 0
        self.max_steps = max_steps

        # Output directory
        self.log_dir = log_dir
        if (log_dir is not None) and (not os.path.exists( self.log_dir )):
            os.makedirs( self.log_dir )

        # Be quiet or verbose
        self.verbose = verbose

        return

#---------------------------------------------------------------------------------------------------------------------------------------------------#
#    def get_valid_search_space_range(
#        self,
#        ):
#        supra_glottal_parameter_info = vtl.get_param_info( 'tract' )
#        supra_glottal_parameter_info.loc[ 'LD', 'min' ] = -0.5
#        supra_glottal_parameter_info.loc[ 'LD', 'max' ] = 2.0
#
#        supra_glottal_min_state = supra_glottal_parameter_info.loc[ self.optimization_parameters, 'min' ].to_numpy( dtype = float )
#        supra_glottal_max_state = supra_glottal_parameter_info.loc[ self.optimization_parameters, 'max' ].to_numpy( dtype = float )
#
#        #glottal_min_state = np.array( [ [ -0.05, -1.0 ] for x in self.optimization_set.articulatory_states if x.optimize_glottis ] ).flatten()
#        #glottal_max_state = np.array( [ [ 0.1, 1.0 ] for x in self.optimization_set.articulatory_states if x.optimize_glottis ] ).flatten()
#
#        min_state = supra_glottal_min_state #np.concatenate( [ supra_glottal_min_state, glottal_min_state ] )
#        max_state = supra_glottal_max_state #np.concatenate( [ supra_glottal_max_state, glottal_max_state ] )
#        #print(min_state)
#        #print(max_state)
#        return min_state, max_state
#---------------------------------------------------------------------------------------------------------------------------------------------------#
    def _initialize_run( self ):
        log = []
        log.append( 'Running VTL vocal tract shape optimization:' )
        log.append( 'Optimizing the following {} parameters:'.format( len( self.optimization_parameters ) ) )
        log.append( '    {}'.format( self.optimization_parameters ) )
        log.append( 'Using following optimization algorithm: {}'.format( self.optimization_policy.__name__ ) )
        log.append( 'With the following optimization parameters:' )
        for key, val in self.optimization_policy_kw.items():
            log.append( '    {}: {}'.format( key, val ) )
        log.append( 'Saving the results to: {}'.format( self.log_dir ) )
        if self.verbose:
            print( '\n'.join( log ) )
        return
#---------------------------------------------------------------------------------------------------------------------------------------------------#
    def _finalize_run( self, elapsed_time ):
        log = []

        if len( self.log_loss_data ) > 1:
            result_state = self.log_loss_data[ -1 ][ 'supra_glottal_state' ]
        else:
            result_state = None

        simulaton_result = dict(
            loss_data = self.log_loss_data,
            reference_state = state_to_dict( self.sgs_r ),
            result_state = result_state,
            computation_time = elapsed_time,
            constriction = self.constriction,
            n_constrictions = self.n_constrictions,
            tight_constriction_length = self.tight_constriction_length,
            close_constriction_length = self.close_constriction_length,
            place_of_articulation = self.place_of_articulation,
            )

        log.append( '\nSimulation finished.' )
        if self.verbose:
            print( '\n'.join( log ) )

        return simulaton_result
#---------------------------------------------------------------------------------------------------------------------------------------------------#
    def run( self ):
        self._initialize_run()
        start = time.time()
        try:
            self.optimization_policy(
                target_function = self.objective_function,
                **self.optimization_policy_kw,
                )
        except EarlyStoppingException:
            if self.verbose:
                print( 'Stopped optimization at step: {}'.format( self.step ) )
        end = time.time()
        elapsed_time = end - start
        result_params = self._finalize_run( elapsed_time )
        return result_params
#---------------------------------------------------------------------------------------------------------------------------------------------------#
    def articulator_exclusion_loss(
        self,
        constriction_articulators,
        ):
        true_positive = 0
        false_positive = 0
        for articulator in constriction_articulators:
            #print(articulator)
            #print(self.place_of_articulation)
            if articulator in self.place_of_articulation:
                #constriction_loss += 1
                true_positive +=1
            else:
                false_positive += 1
        precision = true_positive / ( true_positive + false_positive )
        recall = true_positive / len( constriction_articulators )
        try:
            f1 = 2 * ( precision * recall ) / ( precision + recall )
        except ZeroDivisionError:
            f1 = 0
        return ( 1 - f1 )
#---------------------------------------------------------------------------------------------------------------------------------------------------#
    def articulator_inclusion_loss(
        self,
        constriction_articulators,
        ):
        true_positive = 0
        false_positive = 0
        for articulator in self.place_of_articulation:
            #print(articulator)
            #print(self.place_of_articulation)
            if articulator in constriction_articulators:
                #constriction_loss += 1
                true_positive +=1
            else:
                false_positive += 1
        precision = true_positive / ( true_positive + false_positive )
        recall = true_positive / len( self.place_of_articulation )
        try:
            f1 = 2 * ( precision * recall ) / ( precision + recall )
        except ZeroDivisionError:
            f1 = 0
        return ( 1 - f1 )
#---------------------------------------------------------------------------------------------------------------------------------------------------#
    def _closure_loss(
            self,
            tube_state,
            ):
        loss = 0
        tc = tube_state.constriction_data[ 'tight_constrictions' ]
        if tc:
            tc_articulators = tc[0][ 'articulators' ]
            constriction_articulators = [
                art[ 'place_of_articulation' ] for art in tc_articulators
                ]
            if self.tight_constriction_length is not None:
                if len( constriction_articulators ) < self.tight_constriction_length:
                    loss += 25

            cc = tube_state.constriction_data['close_constrictions']
            if cc:
                cc_articulators = cc[0][ 'articulators' ]
                close_constriction_articulators = [
                    art[ 'place_of_articulation' ] for art in cc_articulators
                    if art[ 'constriction_class' ] == 2
                    ]
                if self.close_constriction_length is not None:
                    if len( close_constriction_articulators ) < self.close_constriction_length:
                        loss += 25

            loss += self.articulator_exclusion_loss( constriction_articulators )
            loss += self.articulator_inclusion_loss( constriction_articulators )
        else:
            loss += 50
        return loss
    
    def _constriction_loss(
            self,
            tube_state,
            ):
        loss = 0
        tc = tube_state.constriction_data[ 'tight_constrictions' ]
        if tc:
            tc_articulators = tc[0][ 'articulators' ]
            constriction_articulators = [
                art[ 'place_of_articulation' ] for art in tc_articulators
                ]
            tight_constriction_articulators = [
                art[ 'place_of_articulation' ] for art in tc_articulators
                if art[ 'constriction_class' ] in [ 3, 4 ]
                ]
            if self.tight_constriction_length is not None:
                if len( constriction_articulators ) < self.tight_constriction_length:
                    loss += 25
            if self.close_constriction_length is not None:
                if len( tight_constriction_articulators ) < self.close_constriction_length:
                    loss += 25    
            loss += self.articulator_exclusion_loss( constriction_articulators )
            loss += self.articulator_inclusion_loss( constriction_articulators )
        else:
            loss += 50
        return loss
    
    def calculate_sensory_loss(
        self,
        query: SGS,
        ):
        loss = 0
        #tb = vtl.motor_to_tube(
        #    query,
        #    fast_calculation = True,
        #    )[0]
        tb = vtl._motor_to_tube(
            tract_state=query.to_numpy( transpose = False ).squeeze(),
            fast_calculation = True,
            )

        # Constriction type
        if tb.constriction != self.constriction:
            loss += 100

        # Number of constrictions
        if tb.constriction_data[ 'n_constrictions' ] != self.n_constrictions:
            loss += 100

        # Lengths of constrictions and inclusion/exclusion losses
        if self.constriction == 2:
            loss += self._closure_loss( tb )
        elif self.constriction in [ 3, 4 ]:
            loss += self._constriction_loss( tb )
        else:
            raise ValueError(
                f"""
                The sensory loss only supports constriction classes 2, 3, and 4,
                but received constriction class was {self.constriction}.
                """
                )
        return loss

    def calculate_similarity_loss(
        self,
        reference: SGS,
        query: SGS,
        mask: bool = True,
        ):
        r = reference.to_numpy().squeeze()
        q = query.to_numpy().squeeze()

        if mask:
            mask = get_mask( query, self.optimization_parameters )
            r = r[ mask ]
            q = q[ mask ]

        loss = cosine( r, q )
        return loss

    def get_query_state( self, parameters_queries ):
        sgs = SGS( self.reference_state )
        for opt_p, p in zip( parameters_queries, self.optimization_parameters ):
            sgs[ p ] = opt_p
        return sgs
    
    def objective_function( self, parameters_queries ):
        # Check if stop condition is met
        if self.step >= self.max_steps:
            raise EarlyStoppingException

        # Get the supra-glottal state
        sgs_q = self.get_query_state( parameters_queries )

        # Compute losses
        constriction_loss = self.calculate_sensory_loss(
            query = sgs_q,
            )
        similarity_loss = self.calculate_similarity_loss(
            reference = self.sgs_r,
            query = sgs_q,
            #mask = True,
            mask = False,
            )
        total_loss = constriction_loss + similarity_loss

        # Handle losses
        if total_loss < self.log_loss_data[ -1 ][ 'total_loss' ]:
            self.log_loss_data.append(
                dict(
                    step = self.step,
                    total_loss = total_loss,
                    constriction_loss = constriction_loss,
                    similarity_loss = similarity_loss,
                    supra_glottal_state = deepcopy( state_to_dict( sgs_q ) ),
                    )
                )
            log = []
            log.append( 'Step: {: <{}}'.format( self.step, len( str( self.max_steps ) ) + 1 ) )
            log.append(
                'constriction loss: {: <{}}'.format( f'{constriction_loss:.5f}', 10 )
                )
            log.append( 
                'similarity loss: {}'.format( similarity_loss )
                )
            if self.verbose:
                print( ' — '.join( log ) )

        self.step += 1
        return total_loss