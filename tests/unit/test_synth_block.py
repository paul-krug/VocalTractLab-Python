import unittest
import numpy as np
from vocaltractlab_cython import synth_block, get_shape

class TestSynthBlock(unittest.TestCase):

    def test_valid_input(self):
        # Create valid input arrays
        tract_params = get_shape( 'a', params='tract' ).reshape(1, -1)
        glottis_params = get_shape( 'modal', params='glottis' ).reshape(1, -1)

        # Synthesize audio
        audio = synth_block(tract_params, glottis_params, state_samples=48000)

        # Check if the audio array is not empty
        self.assertTrue(len(audio) == 48000)

    def test_invalid_tract_params_shape(self):
        # Create an invalid input array with incorrect tract parameters shape
        tract_params = np.array([[0.1, 0.2, 0.3, 0.4], [0.4, 0.5, 0.6, 0.7]])

        # Ensure a ValueError is raised
        with self.assertRaises(ValueError):
            synth_block(tract_params, np.array([[0.7, 0.8, 0.9], [1.0, 1.1, 1.2]]))

    def test_invalid_glottis_params_shape(self):
        # Create an invalid input array with incorrect glottis parameters shape
        glottis_params = np.array([[0.7, 0.8], [1.0, 1.1]])

        # Ensure a ValueError is raised
        with self.assertRaises(ValueError):
            synth_block(np.array([[0.1, 0.2, 0.3], [0.4, 0.5, 0.6]]), glottis_params)

    def test_mismatched_number_of_time_steps(self):
        # Create input arrays with different numbers of time steps
        tract_params = np.array([[0.1, 0.2, 0.3], [0.4, 0.5, 0.6]])
        glottis_params = np.array([[0.7, 0.8, 0.9]])

        # Ensure a ValueError is raised
        with self.assertRaises(ValueError):
            synth_block(tract_params, glottis_params)

if __name__ == '__main__':
    unittest.main()
