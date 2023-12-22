import logging
logger = logging.getLogger()
logger.setLevel(logging.INFO)

import os
import unittest
from vocaltractlab_cython import gesture_file_to_audio, VtlApiError
import numpy as np

class TestGestureFileToAudio(unittest.TestCase):

    def __init__(self, methodName: str = "runTest") -> None:
        super().__init__(methodName)
        self.valid_gesture_file = os.path.join(
            os.path.dirname(__file__),
            'resources',
            'valid_gestural_score.txt',
        )
        self.invalid_gesture_file = os.path.join(
            os.path.dirname(__file__),
            'resources',
            'invalid_gestural_score.txt',
        )
        self.audio_file = os.path.join(
            os.path.dirname(__file__),
            'test_output',
            'gesture_file_to_audio.wav',
        )
        os.makedirs(
            os.path.join(
                os.path.dirname(__file__),
                'test_output',
            ),
            exist_ok=True,
            )

    def test_generate_audio_from_gesture_file(self):
        print( 'test_generate_audio_from_gesture_file' )
        # Test generating audio from a valid gesture file
        audio = gesture_file_to_audio(self.valid_gesture_file)
        self.assertIsInstance(audio, np.ndarray)  # Check if audio is a NumPy array

    def test_save_generated_audio(self):
        print( 'test_save_generated_audio' )
        # Test generating audio and saving it to a WAV file
        audio = gesture_file_to_audio(
            self.valid_gesture_file,
            self.audio_file,
            )
        self.assertIsInstance(audio, np.ndarray)  # Check if audio is a NumPy array

    def test_generate_audio_with_verbose_output(self):
        print( 'test_generate_audio_with_verbose_output' )
        # Test generating audio with verbose API output
        audio = gesture_file_to_audio(
            self.valid_gesture_file,
            verbose_api=True,
            )
        self.assertIsInstance(audio, np.ndarray)  # Check if audio is a NumPy array

    def test_invalid_gesture_file(self):
        print( 'test_invalid_gesture_file' )
        # Test generating audio from an invalid gesture file (should raise an exception)
        with self.assertRaises(VtlApiError):
            gesture_file_to_audio( self.invalid_gesture_file )