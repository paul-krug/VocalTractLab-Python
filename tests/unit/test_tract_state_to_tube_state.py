import unittest
import numpy as np
from vocaltractlab_cython import tract_state_to_tube_state, get_shape

class TestTractStateToTubeState(unittest.TestCase):
    def test_valid_tract_state(self):
        # Create a valid vocal tract state (matching the number of tract parameters)
        vocal_tract_state = get_shape( 'a', params = 'tract' )

        # Compute the tube state
        tube_state = tract_state_to_tube_state(vocal_tract_state)

        # Check if the computed tube state contains the expected information
        with self.subTest(test_case="Valid Vocal Tract State"):
            self.assertTrue("tube_length" in tube_state)
            self.assertTrue("tube_area" in tube_state)
            self.assertTrue("tube_articulator" in tube_state)
            self.assertTrue("incisor_position" in tube_state)
            self.assertTrue("tongue_tip_side_elevation" in tube_state)
            self.assertTrue("velum_opening" in tube_state)

            # Check if the tube length and area arrays have the expected shape
            self.assertEqual(tube_state["tube_length"].shape[0], 40)  # Shpuld be 40
            self.assertEqual(tube_state["tube_area"].shape[0], 40)  # Should be 40

            # Check if other values are not None
            self.assertIsNotNone(tube_state["tube_articulator"])
            self.assertIsNotNone(tube_state["incisor_position"])
            self.assertIsNotNone(tube_state["tongue_tip_side_elevation"])
            self.assertIsNotNone(tube_state["velum_opening"])

    def test_invalid_tract_state_length(self):
        # Create an invalid vocal tract state (wrong length)
        vocal_tract_state = np.array([0.1, 0.2])  # Incorrect length
        
        # Ensure a ValueError is raised
        with self.subTest(test_case="Invalid Vocal Tract State - Wrong Length"):
            with self.assertRaises(ValueError):
                tract_state_to_tube_state(vocal_tract_state)

    def test_invalid_tract_state_dim(self):
        # Create an invalid vocal tract state (not a 1D array)
        vocal_tract_state = np.array([[0.1, 0.2, 0.3]])  # 2D array
        
        # Ensure a ValueError is raised
        with self.subTest(test_case="Invalid Vocal Tract State - Not 1D Array"):
            with self.assertRaises(ValueError):
                tract_state_to_tube_state(vocal_tract_state)

if __name__ == '__main__':
    unittest.main()
