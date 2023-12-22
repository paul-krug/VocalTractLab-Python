import logging
logger = logging.getLogger()
logger.setLevel(logging.INFO)

import os
import unittest
from vocaltractlab_cython import phoneme_file_to_gesture_file, VtlApiError

class TestPhonemeFileToGestureFile(unittest.TestCase):

    def __init__(self, methodName: str = "runTest") -> None:
        super().__init__(methodName)
        self.valid_phoneme_file = os.path.join(
            os.path.dirname(__file__),
            'resources',
            'valid_phoneme_sequence.txt',
        )
        self.invalid_phoneme_file = os.path.join(
            os.path.dirname(__file__),
            'this_file_does_not_exist',
            'invalid_phoneme_sequence.txt',
        )
        self.valid_gesture_file = os.path.join(
            os.path.dirname(__file__),
            'test_output',
            'valid_phoneme_file_to_gesture_file.txt',
        )
        self.invalid_gesture_file = os.path.join(
            os.path.dirname(__file__),
            'test_output',
            'invalid_phoneme_file_to_gesture_file.txt',
        )
        os.makedirs(
            os.path.join(
                os.path.dirname(__file__),
                'test_output',
            ),
            exist_ok=True,
            )

    def test_generate_gestural_score(self):
        # Test generating a gestural score file from a phoneme sequence file
        try:
            phoneme_file_to_gesture_file(
                self.valid_phoneme_file,
                self.valid_gesture_file,
                verbose_api=True,
                )
            self.assertTrue(os.path.exists(self.valid_gesture_file))  # Check if the gesture file was created
        except VtlApiError as e:
            self.fail(f"Failed to generate gestural score file: {e}")

    def test_invalid_phoneme_file(self):
        # Test providing a nonexistent phoneme sequence file (should raise a FileNotFoundError)
        with self.assertRaises(FileNotFoundError):
            phoneme_file_to_gesture_file(
                self.invalid_phoneme_file,
                self.invalid_gesture_file,
                )

if __name__ == '__main__':
    unittest.main()
