import logging
logger = logging.getLogger()
logger.setLevel(logging.INFO)

import unittest
import numpy as np
from vocaltractlab_cython import get_shape, VtlApiError, get_constants

class TestGetShape(unittest.TestCase):

    def test_retrieve_vocal_tract_shape(self):
        # Test retrieving vocal tract shape parameters
        try:
            valid_shape_name = "a"
            vocal_tract_shape = get_shape(valid_shape_name, 'tract')
            self.assertIsInstance(vocal_tract_shape, np.ndarray)  # Check if vocal_tract_shape is a NumPy array
            # Check if the shape array has the expected shape (size)
            vtl_constants = get_constants()
            self.assertEqual(vocal_tract_shape.size, vtl_constants["n_tract_params"])
        except VtlApiError as e:
            self.fail(f"Failed to retrieve vocal tract shape parameters: {e}")

    def test_retrieve_glottis_shape(self):
        # Test retrieving glottis shape parameters
        try:
            valid_shape_name = "modal"
            glottis_shape = get_shape(valid_shape_name, 'glottis')
            self.assertIsInstance(glottis_shape, np.ndarray)  # Check if glottis_shape is a NumPy array
            # Check if the shape array has the expected shape (size)
            vtl_constants = get_constants()
            self.assertEqual(glottis_shape.size, vtl_constants["n_glottis_params"])
        except VtlApiError as e:
            self.fail(f"Failed to retrieve glottis shape parameters: {e}")

    def test_invalid_params_argument(self):
        # Test providing an invalid 'params' argument (should raise a ValueError)
        with self.assertRaises(ValueError):
            shape_name = "valid_shape"
            get_shape(shape_name, 'invalid_params')