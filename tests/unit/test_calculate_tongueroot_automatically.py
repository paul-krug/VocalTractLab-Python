import logging
logger = logging.getLogger()
logger.setLevel(logging.INFO)

import unittest
from vocaltractlab_cython import calculate_tongueroot_automatically, VtlApiError

class TestCalculateTonguerootAutomatically(unittest.TestCase):

    def test_enable_automatic_calculation(self):
        # Test enabling automatic calculation
        try:
            calculate_tongueroot_automatically(True)
            self.assertTrue(True)  # No exception should be raised
        except VtlApiError as e:
            self.fail(f"Failed to enable automatic calculation: {e}")

    def test_disable_automatic_calculation(self):
        # Test disabling automatic calculation
        try:
            calculate_tongueroot_automatically(False)
            self.assertTrue(True)  # No exception should be raised
        except VtlApiError as e:
            self.fail(f"Failed to disable automatic calculation: {e}")

    def test_invalid_input_type(self):
        # Test for invalid input type (non-boolean)
        with self.assertRaises(TypeError):
            calculate_tongueroot_automatically("InvalidInput")