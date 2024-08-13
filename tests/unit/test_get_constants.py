import logging
logger = logging.getLogger()
logger.setLevel(logging.INFO)

import unittest
from vocaltractlab_cython import get_constants, VtlApiError

class TestGetConstants(unittest.TestCase):

    def test_retrieve_constants(self):
        # Test retrieving constants and parameters
        try:
            constants = get_constants()
            self.assertIsInstance(constants, dict)  # Check if constants is a dictionary
            # Check if the expected keys are present in the dictionary
            self.assertTrue("sr_audio" in constants)
            self.assertTrue("sr_internal" in constants)
            self.assertTrue("n_tube_sections" in constants)
            self.assertTrue("n_tract_params" in constants)
            self.assertTrue("n_glottis_params" in constants)
            self.assertTrue("n_samples_per_state" in constants)
        except VtlApiError as e:
            self.fail(f"Failed to retrieve constants: {e}")
        except ValueError as ve:
            self.fail(f"Invalid values retrieved: {ve}")