import unittest
import numpy as np
import os
from vocaltractlab_cython import tract_state_to_svg, get_shape

class TestTractStateToSvg(unittest.TestCase):

    def __init__(self, methodName: str = "runTest") -> None:
        super().__init__(methodName)
        self.valid_svg_file = os.path.join(
            os.path.dirname(__file__),
            'test_output',
            'valid_shape.svg',
        )
        self.invalid_svg_file = os.path.join(
            os.path.dirname(__file__),
            'test_output',
            'invalid_shape.svg',
        )
        os.makedirs(
            os.path.join(
                os.path.dirname(__file__),
                'test_output',
            ),
            exist_ok=True,
            )
        
    def test_valid_tract_state(self):
        # Create a valid vocal tract state (matching the number of tract parameters)
        valid_vocal_tract_state = get_shape( 'a', params = 'tract' )

        # Ensure no exception is raised when exporting the SVG
        with self.subTest(test_case="Valid Vocal Tract State"):
            tract_state_to_svg(
                valid_vocal_tract_state,
                self.valid_svg_file,
                )

        # Check if the SVG file was created
        self.assertTrue(os.path.exists(self.valid_svg_file), "SVG file was not created")

    def test_invalid_tract_state(self):
        # Create an invalid vocal tract state (wrong length)
        vocal_tract_state = np.array([0.1, 0.2, 0.3])  # Incorrect length

        # Ensure a ValueError is raised
        with self.subTest(test_case="Invalid Vocal Tract State - Wrong Length"):
            with self.assertRaises(ValueError):
                tract_state_to_svg(vocal_tract_state, self.invalid_svg_file)

    def test_invalid_tract_state_dim(self):
        # Create an invalid vocal tract state (not a 1D array)
        vocal_tract_state = np.array([[0.1, 0.2, 0.3, 0.4, 0.5]])  # 2D array

        # Ensure a ValueError is raised
        with self.subTest(test_case="Invalid Vocal Tract State - Not 1D Array"):
            with self.assertRaises(ValueError):
                tract_state_to_svg(vocal_tract_state, self.invalid_svg_file)

if __name__ == '__main__':
    unittest.main()
