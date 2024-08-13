import unittest
import os
from vocaltractlab_cython import motor_file_to_audio_file, VtlApiError

class TestMotorFileToAudioFile(unittest.TestCase):

    def __init__(self, methodName: str = "runTest") -> None:
        super().__init__(methodName)
        self.valid_motor_file = os.path.join(
            os.path.dirname(__file__),
            'resources',
            'valid_motor_series.txt',
        )
        self.invalid_motor_file = os.path.join(
            os.path.dirname(__file__),
            'this_file_does_not_exist',
            'invalid_motor_series.txt',
        )
        self.valid_audio_file = os.path.join(
            os.path.dirname(__file__),
            'test_output',
            'valid_motor_file_to_audio_file.wav',
        )
        self.invalid_audio_file = os.path.join(
            os.path.dirname(__file__),
            'test_output',
            'invalid_motor_file_to_audio_file.wav',
        )
        os.makedirs(
            os.path.join(
                os.path.dirname(__file__),
                'test_output',
            ),
            exist_ok=True,
            )

    def test_valid_conversion(self):
        try:
            # Convert the motor file to an audio file
            motor_file_to_audio_file(
                self.valid_motor_file,
                self.valid_audio_file,
                )
            # Assert that the audio file was generated
            self.assertTrue(os.path.exists(self.valid_audio_file))
        except VtlApiError as e:
            self.fail(f"Conversion failed with error: {e}")

    def test_invalid_motor_file(self):
        # Attempt to convert a non-existent motor file
        with self.assertRaises(FileNotFoundError):
            motor_file_to_audio_file(
                self.invalid_motor_file,
                self.invalid_audio_file,
                )

if __name__ == '__main__':
    unittest.main()
