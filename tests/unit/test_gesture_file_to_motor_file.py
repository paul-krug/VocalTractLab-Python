import logging
logger = logging.getLogger()
logger.setLevel(logging.INFO)

import os
import unittest
from vocaltractlab_cython import gesture_file_to_motor_file, VtlApiError

class TestGestureFileToMotorFile(unittest.TestCase):

    def __init__(self, methodName: str = "runTest") -> None:
        super().__init__(methodName)
        self.valid_gesture_file = os.path.join(
            os.path.dirname(__file__),
            'resources',
            'valid_gestural_score.txt',
        )
        self.invalid_gesture_file = os.path.join(
            os.path.dirname(__file__),
            'resources',
            'invalid_gestural_score.txt',
        )
        self.valid_motor_file = os.path.join(
            os.path.dirname(__file__),
            'test_output',
            'valid_gesture_file_to_motor_file.txt',
        )
        self.invalid_motor_file = os.path.join(
            os.path.dirname(__file__),
            'test_output',
            'invalid_gesture_file_to_motor_file.txt',
        )
        os.makedirs(
            os.path.join(
                os.path.dirname(__file__),
                'test_output',
            ),
            exist_ok=True,
            )

    def test_generate_motor_file_from_gesture_file(self):
        # Test generating a motor file from a valid gesture file
        gesture_file_to_motor_file(
            self.valid_gesture_file,
            self.valid_motor_file,
            )
        self.assertTrue(os.path.exists(self.valid_motor_file))  # Check if the motor file was generated

    def test_generate_motor_file_invalid_gesture_file(self):
        # Test generating a motor file from an invalid gesture file (should raise an exception)
        with self.assertRaises(VtlApiError):
            gesture_file_to_motor_file(
                self.invalid_gesture_file,
                self.invalid_motor_file,
                )