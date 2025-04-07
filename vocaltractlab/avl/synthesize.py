


import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

#import VocalTractLab as vtl
#import pandas as pd
import numpy as np
import os
import yaml
import tqdm

import matplotlib.pyplot as plt

import argparse
from copy import deepcopy
from itertools import chain

from tools_io import save, load
from tools_mp import process

import torch
import torchaudio

import random
from random import getrandbits
from random import uniform
from random import randint

from target_approximation import TargetSequence
from target_approximation.vocaltractlab import SupraGlottalSeries as SGS
from target_approximation.vocaltractlab import GlottalSeries
from target_approximation.vocaltractlab import SupraGlottalSequence
from target_approximation.vocaltractlab import GlottalSequence
from target_approximation.vocaltractlab import MotorSequence

#import vocaltractlab_cython as vtl

from .core import PhonemeInventory

from typing import Union

def hz_to_st(
    frequency_hz: float,
    reference: float = 1.0,
    ):
    return 12.0*np.log( frequency_hz / reference ) / np.log(2.0)

def st_to_hz(
    frequency_st: float,
    reference: float = 1.0,
    ):
  return reference*( 2**(frequency_st / 12.0) )



class MonteCarloGenerator():
    def __init__(
        self,
        phoneme_inventory: Union[ str, PhonemeInventory ],

        poa_ts2 = [ 
            'T2_T3',
            'T3_T4',
            'T4_T5',
            'T3',
            'T4',
            ],
        poa_ts3 = [ 
            #'T2_T3',
            #'T3_T4',
            #'T4_T5',
            'T5_T6',
            'T6_T7',
            #'T3',
            #'T4',
            'T5',
            'T6',
            'T7',
            ],
        ):

        # TODO: use this for the loading of pre-defined inventories
        #inventory_path = os.path.join(
        #    os.path.dirname( __file__ ),
        #    inventory_path,
        #    )
        if isinstance( phoneme_inventory, str ):
            phoneme_inventory = PhonemeInventory.from_yaml( phoneme_inventory )
        elif not isinstance( phoneme_inventory, PhonemeInventory ):
            raise ValueError(
                'phoneme_inventory must be either str or PhonemeInventory'
                )
        self.inventory = phoneme_inventory


        self.c_types = [ 
            'closure',
            'fricative',
            'lateral',
            'nasal',
            ]
        
        self.voiced_types = [ 'vowel', 'lateral', 'nasal' ]
        self.mixedvoiced_types = [ 'closure', 'fricative' ]
        self.unvoiced_types = []

        self.poa_types = [
            'T2_T3',
            'T3_T4',
            'T4_T5',
            'T5_T6',
            'T6_T7',
            'T3',
            'T4',
            'T5',
            'T6',
            'T7',
            'LF',
            'L',
            ]
        
        self.poa_ts2 = poa_ts2
        self.poa_ts3 = poa_ts3
        return

    def generate_phonemes(
        self,
        cv_ratio: float = 0.5,
        n_targets_range: list = [ 1, 50 ],
        ):
        """
        Generate a list of phonemes based on the number of targets and the ratio of consonants to vowels.
        A cv_ratio of 0.5 will result in a 50/50 distribution of vowels and consonants.
        A cv_ratio of 0.0 will result in only vowels, a cv_ratio of 1.0 will result in only consonants.
        """

        n_targets = random.randint(
            n_targets_range[0],
            n_targets_range[1],
            )
        
        if n_targets == 1:
            phoneme_categories = [ 'v' ]
        else:            
            # old code hardcoded a 50/50 distribution of vowels and consonants
            # now we use the cv_ratio to determine the distribution
            phoneme_categories = [
                'v' if random.random() < cv_ratio else 'c' for _ in range(
                    0,
                    n_targets,
                    )
                ]
            
        phonemes = []
        for category in phoneme_categories:
            if category == 'v':
                phonemes.append( 'vowel' )
            if category == 'c':
                phonemes.append(
                    random.choice( self.c_types )
                    )
                
        return phonemes
    
    def generate_voicings(
            self,
            phonemes,
            ):
        
        voicings = []

        for phoneme in phonemes:
            if phoneme in self.voiced_types:
                voicings.append( 'voiced' )
            elif phoneme in self.mixedvoiced_types:
                if getrandbits(1):
                    voicings.append( 'voiced' )
                else:
                    voicings.append( 'unvoiced' )
            elif phoneme in self.unvoiced_types:
                voicings.append( 'unvoiced' )
            else:
                raise ValueError(
                    f"""
                    Phoneme category must not be: {phoneme}
                    """
                    )
            
        return voicings

    def generate_durations(
        self,
        phonemes,
        duration_range,
        ):
        
        durations = [
            random.uniform(
                duration_range[0],
                duration_range[1],
                )
            for _ in phonemes
            ]
        
        return durations


    def generate_states(
            self,
            phonemes,
            ):

        sgs_states = []

        for phoneme in phonemes:

            if phoneme == 'vowel':
                state = random.choice(
                    self.inventory.vowels(return_type = 'numpy')
                    )
                sgs_states.append( state )
            else:
                poa = random.choice( self.poa_types )
                state = random.choice(
                    self.inventory.phonemes(
                        poa=poa,
                        return_type = 'numpy',
                        )
                    )
                sgs = SGS( state )
                #print( sgs )
                #sto
                if poa in self.poa_ts3:
                    #print( sgs )
                    if phoneme == 'fricative':
                        #sgs.states.loc[ :, 'TS3' ] = 1.0
                        sgs[ 'TS3' ] = 1.0
                    elif phoneme == 'lateral':
                        sgs[ 'TS3' ] = -1.0
                    elif phoneme == 'nasal':
                        sgs[ 'VO' ] = 0.5
                    #print(sgs)
                    #stop

                if poa in self.poa_ts2:
                    if phoneme == 'fricative':
                        sgs[ 'TS2' ] = 1.0
                    elif phoneme == 'nasal':
                        sgs[ 'VO' ] = 0.5

                if poa in [ 'L' ]:
                    if phoneme == 'nasal':
                        sgs[ 'VO' ] = 0.5

                #print( sgs )
                #stop

                sgs_states.append( 
                    sgs.to_numpy().squeeze()
                    )
                
        sgs = SGS( np.array( sgs_states ) )

        return sgs
    
    def generate_glottal_durations_old(
            self,
            phonemes,
            voicings,
            supra_glottal_durations,
            ):
        glottal_durations = []
        netto_duration_change = 0.0
        for index, phoneme in enumerate( phonemes ):
            if index < len( phonemes ) - 1:
                if voicings[ index + 1 ] == 'unvoiced':
                    if phonemes[ index + 1 ] == 'closure':
                        glottal_durations.append(
                            supra_glottal_durations[ index ] + 0.05
                            )
                        netto_duration_change += 0.05
                    elif phonemes[ index + 1 ] == 'fricative':
                        glottal_durations.append(
                            supra_glottal_durations[ index ] - 0.03
                            )
                        netto_duration_change -= 0.03
                    else:
                        glottal_durations.append(
                            supra_glottal_durations[ index ]
                            )
                else:
                    glottal_durations.append(
                            supra_glottal_durations[ index ]
                            )
            if np.sum( glottal_durations ) >= np.sum( supra_glottal_durations ):
                break
        glottal_durations.append( 
            np.sum( supra_glottal_durations ) - np.sum( glottal_durations )
            )
        #print( supra_glottal_durations )
        #print( glottal_durations )
        #stop
        return glottal_durations
    
    def generate_glottal_durations(
            self,
            phonemes,
            voicings,
            supra_glottal_durations,
            ):
        glottal_durations = []
        #netto_duration_change = 0.0
        sg_boundaries = [ 0.0 ]
        boundaries = [ 0.0 ]
        for index, phoneme in enumerate( phonemes ):
            if index < len( phonemes ) - 1:
                sg_boundary = sg_boundaries[ -1 ] + supra_glottal_durations[ index ]
                boundary = sg_boundaries[ -1 ] + supra_glottal_durations[ index ]
                #if index > 0:
                #    prev_unvoiced = voicings[ index - 1 ] == 'unvoiced'
                #else:
                #    prev_unvoiced = False
                current_unvoiced = voicings[ index ] == 'unvoiced'
                next_unvoiced = voicings[ index + 1 ] == 'unvoiced'
                current_closure = phonemes[ index ] == 'closure'
                next_closure = phonemes[ index + 1 ] == 'closure'
                current_fricative = phonemes[ index ] == 'fricative'
                next_fricative = phonemes[ index + 1 ] == 'fricative'
                if not (current_unvoiced and next_unvoiced):
                    if (current_unvoiced and current_closure) or (next_unvoiced and next_closure):
                        boundary += 0.05
                    elif (current_unvoiced and current_fricative) or (next_unvoiced and next_fricative):
                        boundary -= 0.03
                    #else:
                    #    glottal_durations.append(
                    #        supra_glottal_durations[ index ]
                    #        )
                #else:
                #    glottal_durations.append(
                #            supra_glottal_durations[ index ]
                #            )
                sg_boundaries.append( sg_boundary )
                boundaries.append( boundary )
        sg_boundaries.append( np.sum( supra_glottal_durations ) )
        boundaries.append( np.sum( supra_glottal_durations ) )
            #if np.sum( glottal_durations ) >= np.sum( supra_glottal_durations ):
            #    break

        # reordering of boundaries should not be necessary anymore
        final_boundaries = sorted( boundaries )
        #final_boundaries = boundaries
        glottal_durations = [
            final_boundaries[ idx + 1 ] - final_boundaries[ idx ]
            for idx in range( 0, len( final_boundaries ) - 1 )
            ]

        
        #glottal_durations.append( 
        #    np.sum( supra_glottal_durations ) - np.sum( glottal_durations )
        #    )
        #print( supra_glottal_durations )
        #print( glottal_durations )
        #stop
        return glottal_durations
    
    def generate_pitch_targets(
            self,
            total_duration,
            n_target_range,
            pitch_mean_range = [ 100, 200 ],
            pitch_range_st = [ -6, 6 ],
            pitch_time_constant = 0.02,
            ):
        pitch_mean_range = np.array( pitch_mean_range ).reshape( -1 )

        n_pitch_targets = random.randint( n_target_range[0], n_target_range[1] )

        if len( pitch_mean_range ) == 1:
            mean_pitch = pitch_mean_range[0]
        elif len( pitch_mean_range ) == 2:
            mean_pitch = random.uniform( pitch_mean_range[0], pitch_mean_range[1] )
        else:
            raise ValueError( 'pitch_mean_range must be int or tuple: ( lim_low, lim_high )' )

        mean_pitch_st = hz_to_st( mean_pitch )
        pitch_offsets = [
            st_to_hz( mean_pitch_st + random.uniform( pitch_range_st[0], pitch_range_st[1] ) )
            for _ in range( 0, n_pitch_targets )
            ]
        durations = [-1]
        # while durations contain values <= 0.0, calc drurations from dirichlet distribution
        # Note while dirichlet should not give negative values, it can happen in edge cases
        while any( [ x <= 0.0 for x in durations ] ):
            durations = np.random.dirichlet( [ 1 for _ in pitch_offsets ] ) * total_duration

        f0 = TargetSequence.from_offsets(
            b = np.array( pitch_offsets ),
            tau = pitch_time_constant,
            duration = durations,
            tiers = [ 'F0' ],
            )

        return f0
    
    def generate_glottal_states(
            self,
            voicings,
            ):
        
        modal_shape = [ 1.200000e+02,  8.000000e+03,  1.020000e-02,  2.035000e-02,
        5.000000e-02,  1.222044e+00,  1.000000e+00,  5.000000e-02,
        0.000000e+00,  2.500000e+01, -1.000000e+01]
        g_series = GlottalSeries( np.array( [ modal_shape for _ in voicings ] ) )

        for voicing_idx, voicing in enumerate( voicings ):
            if voicing == 'unvoiced':
                XB = 0.1
                RA = 0.0
                #glottal_sequence.states.loc[ 0, 'F0' ] = 180
                g_series[ 'XB' ][ voicing_idx ] = XB
                g_series[ 'XT' ][ voicing_idx ] = XB
                g_series[ 'CA' ][ voicing_idx ] = XB
                g_series[ 'RA' ][ voicing_idx ] = RA

        return g_series

    def generate_pressure_targets(
            self,
            durations,
            ):
        
        d_0 = 0.001
        d_2 = 0.05
        d_1 = np.sum( durations ) - (d_0 + d_2)

        pressure = TargetSequence.from_offsets(
            b = [ 0, 8000, 0 ],
            tau = 0.005,
            duration = [ d_0, d_1, d_2 ],
            #time_constants = [ 0.005, 0.005 ],
            #onset_state = 0,
            tiers = [ 'PR' ],
        )
        
        return pressure
    
    def _generate_motor_targets(
            self,
            duration_range,
            n_targets_range,
            cv_ratio,
            sg_tau,
            g_tau,
            ):
        
        phonemes = self.generate_phonemes( cv_ratio, n_targets_range )
        voicings = self.generate_voicings( phonemes )
        sg_durations = self.generate_durations( phonemes, duration_range )
        #print(np.array(sg_durations).shape)
        sg_states = self.generate_states( phonemes )
        #print( sg_states.to_numpy().shape )
        g_states = self.generate_glottal_states( voicings )

        g_durations = self.generate_glottal_durations(
            phonemes = phonemes,
            voicings = voicings,
            supra_glottal_durations = sg_durations,
            )
        f0_tgs = self.generate_pitch_targets(
            total_duration = np.sum( sg_durations ),
            n_target_range = [
                len( phonemes ), #n_targets,
                n_targets_range[1],
                ],
            )
        pr_tgs = self.generate_pressure_targets(
            durations = sg_durations,
            )

        sgs = SupraGlottalSequence.from_offsets(
            b = sg_states.to_numpy(),
            tau = sg_tau,
            duration = sg_durations,
        )

        # For debugging and sanity check (should never happen)
        # if any negative values g_durations, plot sgs and print phonemes
        if any( [ x <= 0.0 for x in g_durations ] ):
            print( phonemes )
            print( voicings)
            print( sg_durations )
            print( g_durations )
            sgs.plot()

        gs = GlottalSequence.from_offsets(
            b = g_states.to_numpy(),
            tau = g_tau,
            duration = g_durations,
        )

        ms = MotorSequence( sgs & gs & f0_tgs & pr_tgs )
        
        return ms
    
    def generate_motor_targets(
            self,
            duration_range = [ 0.055, 0.2 ],
            n_targets_range = [ 1, 50 ],
            cv_ratio = 0.5,
            n_group_range = [1,3],
            onset_range = [ 0.0, 2.5 ],
            offset_range = [ 0.0, 2.5 ],
            suspend_range = [ 0.4, 1.0 ],
            sg_tau = 0.01,
            g_tau = 0.01,
            file_path = None,
            ):
        
        if isinstance(n_group_range, int):
            n_groups = n_group_range
        elif isinstance(n_group_range, list):
            n_groups = random.randint(n_group_range[0], n_group_range[1])
        else:
            raise ValueError(
                'n_group_range must be int or list: [ lim_low, lim_high ]'
                )
        #print( n_groups )

        #ms_list = []
        onset_time = random.uniform( onset_range[0], onset_range[1] )

        # Generate the first sequence
        ms = self._generate_motor_targets(
            duration_range = duration_range,
            n_targets_range = n_targets_range,
            cv_ratio = cv_ratio,
            sg_tau = sg_tau,
            g_tau = g_tau,
            )
        # Add the other sequences
        for _ in range(n_groups-1):
            _ms = self._generate_motor_targets(
                duration_range = duration_range,
                n_targets_range = n_targets_range,
                cv_ratio = cv_ratio,
                sg_tau = sg_tau,
                g_tau = g_tau,
                )
            suspend_time = random.uniform( suspend_range[0], suspend_range[1] )
            ms.extend( suspend_time/2, padding = 'right' )
            _ms.extend( suspend_time/2, padding = 'left' )
            ms = ms + _ms
            
        onset_time = random.uniform( onset_range[0], onset_range[1] )
        offset_time = random.uniform( offset_range[0], offset_range[1] )
        ms.extend( onset_time, padding = 'left' )
        ms.extend( offset_time, padding = 'right' )

        if file_path is not None:
            ms.save( file_path )

        return ms
    
    # TODO: move this unit test thing somewhere else
    '''
    def _test_motor_len(self):
        sg_states = vtl.get_shape( 'a', 'tract' )
        g_states = vtl.get_shape( 'modal', 'glottis' )
        g_durations = [10.0]
        sg_durations = [10.0]
        sg_tau = 0.01
        g_tau = 0.01
        sgs = SupraGlottalSequence.from_offsets(
            b = sg_states.reshape( (19,1) ),
            tau = sg_tau,
            duration = sg_durations,
        )
        gs = GlottalSequence.from_offsets(
            b = g_states.reshape( (11,1) ),
            tau = g_tau,
            duration = g_durations,
        )

        ms = MotorSequence( sgs & gs )
        #m_series = ms.to_numpy( sr = 441 ).T
        m_series = ms.to_numpy( sr = float(44100/110) ).T
        sg_series = m_series[ :, :19 ]
        g_series = m_series[ :, 19: ]

        #audio = vtl.synth_block( sg_series, g_series, state_samples= 100 )
        audio = vtl.synth_block( sg_series, g_series, state_samples= 110 )
        audio /= np.max( np.abs( audio ) )
        print( audio.shape )
        # save audio with torchaudio
        torchaudio.save(
            'test_1s.wav',
            torch.tensor(
                audio,
                dtype= torch.float32,
                ).unsqueeze(0),
                44100,
                )
        plt.plot( audio )
        plt.show()

        return
    '''

    def sample(
        self,
        n_samples,
        duration_range = [ 0.055, 0.2 ],
        n_targets_range = [ 1, 50 ],

        n_group_range = [1,3],
        onset_range = [ 0.0, 2.5 ],
        offset_range = [ 0.0, 2.5 ],
        suspend_range = [ 0.4, 1.0 ],
        cv_ratio = 0.5,
        
        out_dir = None,
        verbose = True,
        workers = None,
        save_as = 'yaml.gz',
        ):
        if out_dir is None:
            args = [ {} for _ in range( 0, n_samples ) ]
            return_data = True
        else:
            os.makedirs(
                out_dir,
                exist_ok = True,
                )
            args = [
                dict(
                    file_path = os.path.join(
                        out_dir,
                        f'mc_{x}.{save_as}',
                        ),
                    )
                for x in range( 0, n_samples )
            ]
            return_data = False
        x = process(
            function = self.generate_motor_targets,
            args = args,
            return_data = return_data,
            verbose = verbose,
            workers = workers,
            mp_threshold = 4,
            )
        return x