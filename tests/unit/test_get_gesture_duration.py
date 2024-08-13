import logging
logger = logging.getLogger()
logger.setLevel(logging.INFO)

import os
import unittest
from vocaltractlab_cython import get_gesture_duration, VtlApiError, get_constants

class TestGetGestureDuration(unittest.TestCase):

    def test_retrieve_duration_information(self):
        # Test retrieving duration information
        try:
            valid_gesture_file = os.path.join(
                os.path.dirname(__file__),
                'resources',
                'valid_gestural_score.txt',
            )
            duration_info = get_gesture_duration(valid_gesture_file)
            self.assertIsInstance(duration_info, dict)  # Check if duration_info is a dictionary
            # Check if the expected keys are present in the dictionary
            self.assertTrue("n_audio_samples" in duration_info)
            self.assertTrue("n_gesture_samples" in duration_info)
            self.assertTrue("duration" in duration_info)
            # Check if the values are non-negative
            self.assertTrue(duration_info["n_audio_samples"] >= 0)
            self.assertTrue(duration_info["n_gesture_samples"] >= 0)
            self.assertTrue(duration_info["duration"] >= 0.0)
        except VtlApiError as e:
            self.fail(f"Failed to retrieve duration information: {e}")