import unittest
import numpy as np
from vocaltractlab_cython import tract_state_to_transfer_function, get_shape

class TestTractStateToTransferFunction(unittest.TestCase):
    def test_valid_tract_state(self):
        # Create a valid vocal tract state (matching the number of tract parameters)
        vocal_tract_state = get_shape( 'a', params = 'tract' )
        
        # Compute the transfer function
        transfer_function = tract_state_to_transfer_function(vocal_tract_state)

        # Check if the computed transfer function contains magnitude and phase spectra
        with self.subTest(test_case="Valid Vocal Tract State"):
            self.assertTrue("magnitude_spectrum" in transfer_function)
            self.assertTrue("phase_spectrum" in transfer_function)
            
            # Check if the number of spectrum samples matches the default value
            self.assertEqual(transfer_function["n_spectrum_samples"], 8192)

            # Check if magnitude spectrum has the expected shape
            self.assertEqual(transfer_function["magnitude_spectrum"].shape[0], 8192)
            self.assertGreater(np.max(transfer_function["magnitude_spectrum"]), 0.0)

    def test_invalid_tract_state_length(self):
        # Create an invalid vocal tract state (wrong length)
        vocal_tract_state = np.array([0.1, 0.2])  # Incorrect length
        
        # Ensure a ValueError is raised
        with self.subTest(test_case="Invalid Vocal Tract State - Wrong Length"):
            with self.assertRaises(ValueError):
                tract_state_to_transfer_function(vocal_tract_state)

    def test_invalid_tract_state_dim(self):
        # Create an invalid vocal tract state (not a 1D array)
        vocal_tract_state = np.array([[0.1, 0.2, 0.3]])  # 2D array
        
        # Ensure a ValueError is raised
        with self.subTest(test_case="Invalid Vocal Tract State - Not 1D Array"):
            with self.assertRaises(ValueError):
                tract_state_to_transfer_function(vocal_tract_state)

if __name__ == '__main__':
    unittest.main()
