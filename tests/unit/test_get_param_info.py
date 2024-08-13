import logging
logger = logging.getLogger()
logger.setLevel(logging.INFO)

import unittest
from vocaltractlab_cython import get_param_info, VtlApiError, get_constants

class TestGetParamInfo(unittest.TestCase):

    def test_retrieve_tract_params(self):
        # Test retrieving vocal tract parameters
        try:
            vocal_tract_params = get_param_info('tract')
            self.assertIsInstance(vocal_tract_params, list)  # Check if vocal_tract_params is a list
            # Check if the list contains dictionaries with the expected keys
            for param in vocal_tract_params:
                self.assertTrue("name" in param)
                self.assertTrue("description" in param)
                self.assertTrue("unit" in param)
                self.assertTrue("min" in param)
                self.assertTrue("max" in param)
                self.assertTrue("standard" in param)
        except VtlApiError as e:
            self.fail(f"Failed to retrieve vocal tract parameters: {e}")

    def test_retrieve_glottis_params(self):
        # Test retrieving glottis parameters
        try:
            glottis_params = get_param_info('glottis')
            self.assertIsInstance(glottis_params, list)  # Check if glottis_params is a list
            # Check if the list contains dictionaries with the expected keys
            for param in glottis_params:
                self.assertTrue("name" in param)
                self.assertTrue("description" in param)
                self.assertTrue("unit" in param)
                self.assertTrue("min" in param)
                self.assertTrue("max" in param)
                self.assertTrue("standard" in param)
        except VtlApiError as e:
            self.fail(f"Failed to retrieve glottis parameters: {e}")

    def test_invalid_params_argument(self):
        # Test providing an invalid 'params' argument (should raise a ValueError)
        with self.assertRaises(ValueError):
            get_param_info('invalid_params')