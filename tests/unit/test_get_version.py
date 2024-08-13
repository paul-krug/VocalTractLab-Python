import logging
logger = logging.getLogger()
logger.setLevel(logging.INFO)

import unittest
from vocaltractlab_cython import get_version

class TestGetVersion(unittest.TestCase):

    def test_retrieve_version_information(self):
        # Test retrieving version information
        version = get_version()
        self.assertIsInstance(version, str)  # Check if version is a string
        self.assertGreater(len(version), 0)  # Check if the string is not empty