// ****************************************************************************
// This file is part of VocalTractLab.
// Copyright (C) 2020, Peter Birkholz, Dresden, Germany
// www.vocaltractlab.de
// author: Peter Birkholz
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program. If not, see <http://www.gnu.org/licenses/>.
//
// ****************************************************************************

// ****************************************************************************
// This file defines the entry point for the DLL application, and the functions
// defined here are C-compatible so that they can be used with the 
// MATLAB shared library interface.
// ****************************************************************************

// Make an extern "C" section so that the functions can be accessed from Matlab

#ifdef __cplusplus
extern "C"{ /* start extern "C" */
#endif

// Definition for function export, if the file is compiled as part of a dll.

#ifdef WIN32
  #ifdef _USRDLL
    #define C_EXPORT __declspec(dllexport)
  #else
    #define C_EXPORT
  #endif        // DLL
#else
  #define C_EXPORT
#endif  // WIN32

// ****************************************************************************
// The exported C-compatible functions.
// IMPORTANT: 
// All the functions defined below must be named in the VocalTractLabApi.def 
// file in the project folder, so that they are usable from MATLAB !!!
// ****************************************************************************

// ****************************************************************************
// Init. the synthesis with the given speaker file name, e.g. "JD2.speaker".
// This function should be called before any other function of this API.
// Return values:
// 0: success.
// 1: Loading the speaker file failed.
// ****************************************************************************

C_EXPORT int vtlInitialize(const char *speakerFileName);


// ****************************************************************************
// Clean up the memory and shut down the synthesizer.
// Return values:
// 0: success.
// 1: The API was not initialized.
// ****************************************************************************

C_EXPORT int vtlClose();


// ****************************************************************************
// Returns the version of this API as a string that contains the compile data.
// Reserve at least 32 chars for the string.
// ****************************************************************************

C_EXPORT void vtlGetVersion(char *version);


// ****************************************************************************
// Returns a couple of constants:
// o The audio sampling rate of the synthesized signal.
// o The number of supraglottal tube sections.
// o The number of vocal tract model parameters.
// o The number of glottis model parameters.
//
// Function return value:
// 0: success.
// 1: The API has not been initialized.
// ****************************************************************************

C_EXPORT int vtlGetConstants(int *audioSamplingRate, int *numTubeSections,
  int *numVocalTractParams, int *numGlottisParams);


// ****************************************************************************
// Returns for each vocal tract parameter the minimum value, the maximum value,
// and the neutral value. Each vector passed to this function must have at 
// least as many elements as the number of vocal tract model parameters.
// The "names" string receives the abbreviated names of the parameters separated
// by spaces. This string should have at least 10*numParams elements.
//
// Function return value:
// 0: success.
// 1: The API has not been initialized.
// ****************************************************************************

C_EXPORT int vtlGetTractParamInfo(char *names, double *paramMin, double *paramMax, double *paramNeutral);


// ****************************************************************************
// Returns for each glottis model parameter the minimum value, the maximum value,
// and the neutral value. Each vector passed to this function must have at 
// least as many elements as the number of glottis model parameters.
// The "names" string receives the abbreviated names of the parameters separated
// by spaces. This string should have at least 10*numParams elements.
//
// Function return value:
// 0: success.
// 1: The API has not been initialized.
// ****************************************************************************

C_EXPORT int vtlGetGlottisParamInfo(char *names, double *paramMin, double *paramMax, double *paramNeutral);


// ****************************************************************************
// Returns the vocal tract parameters for the given shape as defined in the
// speaker file.
// The vector passed to this function must have at least as many elements as 
// the number of vocal tract model parameters.
//
// Function return value:
// 0: success.
// 1: The API has not been initialized.
// 2: A shape with the given name does not exist.
// ****************************************************************************

C_EXPORT int vtlGetTractParams(const char *shapeName, double *param);


// ****************************************************************************
// Exports the vocal tract contours for the given vector of vocal tract
// parameters as a SVG file (scalable vector graphics).
//
// Function return value:
// 0: success.
// 1: The API has not been initialized.
// 2: Writing the SVG file failed.
// ****************************************************************************

C_EXPORT int vtlExportTractSvg(double *tractParams, const char *fileName);


// ****************************************************************************
// Provides the tube data (especially the area function) for the given vector
// of tractParams. The vectors tubeLength_cm, tubeArea_cm2, and tubeArticulator, 
// must each have as many elements as tube sections.
// The values incisorPos_cm, tongueTipSideElevation, and velumOpening_cm2 are 
// one double value each.
//
// Function return value:
// 0: success.
// 1: The API has not been initialized.
// ****************************************************************************

C_EXPORT int vtlTractToTube(double* tractParams,
  double* tubeLength_cm, double* tubeArea_cm2, int* tubeArticulator,
  double* incisorPos_cm, double* tongueTipSideElevation, double* velumOpening_cm2);


// ****************************************************************************
// Calculates the volume velocity transfer function of the vocal tract between 
// the glottis and the lips for the given vector of vocal tract parameters and
// returns the spectrum in terms of magnitude and phase.
//
// Parameters in:
// o tractParams: Is a vector of vocal tract parameters with 
//     numVocalTractParams elements.
// o numSpectrumSamples: The number of samples (points) in the requested 
//     spectrum. This number of samples includes the negative frequencies and
//     also determines the frequency spacing of the returned magnitude and
//     phase vectors. The frequency spacing is 
//     deltaFreq = SAMPLING_RATE / numSpectrumSamples.
//     For example, with the sampling rate of 44100 Hz and 
//     numSpectrumSamples = 512, the returned magnitude and phase values are 
//     at the frequencies 0.0, 86.13, 172.3, ... Hz.
//     The value of numSpectrumSamples should not be greater than 16384,
//     otherwise the returned spectrum will be bandlimited to below 10 kHz.
//
// Parameters out:
// o magnitude: Vector of spectral magnitudes at equally spaced discrete 
//     frequencies. This vector mus have at least numSpectrumSamples elements.
// o phase_rad: Vector of the spectral phase in radians at equally 
//     spaced discrete frequencies. This vector must have at least 
//     numSpectrumSamples elements.
//
// Function return value:
// 0: success.
// 1: The API has not been initialized.
// ****************************************************************************

C_EXPORT int vtlGetTransferFunction(double *tractParams, int numSpectrumSamples,
  double *magnitude, double *phase_rad);


// ****************************************************************************
// Resets the time-domain synthesis of continuous speech (using the functions
// vtlSynthesisAddTube() or vtlSynthesisAddTract()). This function must be 
// called every time you start a new synthesis.
//
// Function return value:
// 0: success.
// 1: The API has not been initialized.
// ****************************************************************************

C_EXPORT int vtlSynthesisReset();


// ****************************************************************************
// Synthesize a part of a speech signal with numNewSamples samples, during 
// which the vocal tract tube changes linearly from the tube shape passed to
// the previous call of this function to the tube shape passed to this call.
// To synthesize parts of 5 ms duration, call this function with 
// numNewSamples = 220. The synthesized signal part is written to the array 
// audio (the caller must allocate the memory for the array).
// During the *first* call of this function after vtlSynthesisReset(), no audio
// is synthesized, and numNewSamples should be 0. During the first call, only 
// the initial tube state is set.
//
// The new tube state is given in terms of the following parameters:
// o tubeLength_cm: Vector of tube sections lengths from the glottis (index 0)
//     to the mouth (index numTubeSections; see vtlGetConstants()).
// o tubeArea_cm2: According vector of tube section areas in cm^2.
// o tubeArticulator: Vector of characters (letters) that denote the articulator 
//     that confines the vocal tract at the position of the tube. We discriminate
//     1 (tongue), 2 (lower incisors), 3 (lower lip), 4 (other articulator).
// o incisorPos_cm: Position of the incisors from the glottis.
// o velumOpening_cm2: Opening of the velo-pharyngeal port in cm^2.
// o tongueTipSideElevation: Corresponds to the TS3 parameter of the vocal tract.
// o newGlottisParams: vector with parameters of the glottis model.
//
// Function return value:
// 0: success.
// 1: The API has not been initialized.
// 2: Number of generated audio samples is wrong (may happen when 
//    numNewSamples != 0 during the first call of this function after reset).
// ****************************************************************************

C_EXPORT int vtlSynthesisAddTube(int numNewSamples, double *audio,
  double *tubeLength_cm, double *tubeArea_cm2, int *tubeArticulator,
  double incisorPos_cm, double velumOpening_cm2, double tongueTipSideElevation,
  double *newGlottisParams);


// ****************************************************************************
// Synthesize a part of a speech signal with numNewSamples samples, during 
// which the vocal tract changes linearly from the tract shape passed to
// the previous call of this function to the tract shape passed to this call.
// To synthesize parts of 5 ms duration, call this function with 
// numNewSamples = 220. The synthesized signal part is written to the array 
// audio (the caller must allocate the memory for the array).
// During the *first* call of this function after vtlSynthesisReset(), no audio
// is synthesized, and numNewSamples should be 0. During the first call, only 
// the initial tube state is set.
//
// The new vocal tract state is given in terms of the following parameters:
// o tractParams: Vector of vocal tract parameters.
// o glottisParams: Vector of vocal fold model parameters.
//
// Function return value:
// 0: success.
// 1: The API has not been initialized.
// 2: Number of generated audio samples is wrong (may happen when 
//    numNewSamples != 0 during the first call of this function after reset).
// ****************************************************************************

C_EXPORT int vtlSynthesisAddTract(int numNewSamples, double *audio,
  double *tractParams, double *glottisParams);


// ****************************************************************************
// Synthesize speech with a given sequence of vocal tract model states and 
// glottis model states, and return the corresponding audio signal.
// This function makes successive calls to the function vtlSynthesisAddTract().
//
// Parameters (in/out):
// o tractParams (in): Is a concatenation of vocal tract parameter vectors
//     with the total length of (numVocalTractParams*numFrames) elements.
// o glottisParams (in): Is a concatenation of glottis parameter vectors
//     with the total length of (numGlottisParams*numFrames) elements.
// o numFrames (in): Number of successive states of the glottis and vocal tract
//     that are going to be concatenated.
// o frameStep_samples (in): The number of audio samples between adjacent 
//     frames (states). A typical value is 220, which corresponds to 5 ms.
// o audio (out): The resulting audio signal with sample values in the range 
//     [-1, +1] and with the sampling rate audioSamplingRate. The signal
//     will have (numFrames-1) * frameStep_samples samples, so the array must
//     be at least of this size.
// o enableConsoleOutput (in): Set to 1, if you want to allow output about the
//   synthesis progress in the console window. Otherwise, set it to 0.
//
// Function return value:
// 0: success.
// 1: The API has not been initialized.
// ****************************************************************************

C_EXPORT int vtlSynthBlock(double *tractParams, double *glottisParams,
  int numFrames, int frameStep_samples, double *audio, int enableConsoleOutput);


// ****************************************************************************
// Test function for this API.
// Audio should contain at least 44100 double values.
// Run this WITHOUT calling vtlInitialize() !
// ****************************************************************************

C_EXPORT int vtlApiTest(const char *speakerFileName, double *audio, int *numSamples);


// ****************************************************************************
// This function converts a segment sequence file (a TXT file containing the 
// sequence of speech segments in SAMPA and the associated durations) with the 
// name segFileName into a gestural score file (gesFileName).
// The f0 tier in the gestural score is set to a "standard" f0.
//
// Function return value:
// 0: success.
// 1: The API was not initialized.
// 2: Loading the segment sequence file failed.
// 3: Saving the gestural score file failed.
// ****************************************************************************

C_EXPORT int vtlSegmentSequenceToGesturalScore(const char *segFileName, const char *gesFileName);


// ****************************************************************************
// This function directly converts a gestural score to an audio signal or file.
// Parameters:
// o gesFileName (in): Name of the gestural score file to synthesize.
// o wavFileName (in): Name of the audio file with the resulting speech signal.
//     This can be the empty string "" if you do not want to save a WAV file.
// o audio (out): The resulting audio signal with sample values in the range 
//     [-1, +1] and with the sampling rate audioSamplingRate. Make sure that
//     this buffer is big enough for the synthesized signal. If you are not 
//     interested in the audio signal, set this pointer to NULL.
// o numSamples (out): The number of audio samples in the synthesized signal.
//     If you are not interested in this value, set this pointer to NULL.
// o enableConsoleOutput (in): Set to 1, if you want to allow output about the
//   synthesis progress in the console window. Otherwise, set it to 0.
//
// Function return value:
// 0: success.
// 1: The API was not initialized.
// 2: Loading the gestural score file failed.
// 3: Values in the gestural score file are out of range.
// 4: The WAV file could not be saved.
// ****************************************************************************

C_EXPORT int vtlGesturalScoreToAudio(const char* gesFileName, const char* wavFileName,
  double* audio, int* numSamples, int enableConsoleOutput);


// ****************************************************************************
// This function directly converts a gestural score to a tract sequence file.
// The latter is a text file containing the parameters of the vocal fold and 
// vocal tract models in steps of 5 ms.

// Parameters:
// o gesFileName (in): Name of the gestural score file to convert.
// o tractSequenceFileName (in): Name of the tract sequence file.
//
// Function return value:
// 0: success.
// 1: The API was not initialized.
// 2: Loading the gestural score file failed.
// 3: Values in the gestural score file are out of range.
// 4: The tract sequence file could not be saved.
// ****************************************************************************

C_EXPORT int vtlGesturalScoreToTractSequence(const char* gesFileName, 
  const char* tractSequenceFileName);


// ****************************************************************************
// This function converts a tract sequence file into an audio signal or file.
// Parameters:
// o tractSequenceFileName (in): Name of the tract sequence file to synthesize.
// o wavFileName (in): Name of the audio file with the resulting speech signal.
//     This can be the empty string "" if you do not want to save a WAV file.
// o audio (out): The resulting audio signal with sample values in the range 
//     [-1, +1] and with the sampling rate audioSamplingRate. Make sure that
//     this buffer is big enough for the synthesized signal. If you are not 
//     interested in the audio signal, set this pointer to NULL.
// o numSamples (out): The number of audio samples in the synthesized signal.
//     If you are not interested in this value, set this pointer to NULL.
//
// Function return value:
// 0: success.
// 1: The API was not initialized.
// 2: Synthesis of the tract sequence file failed.
// 3: The WAV file could not be saved.
// ****************************************************************************

C_EXPORT int vtlTractSequenceToAudio(const char* tractSequenceFileName,
  const char* wavFileName, double* audio, int* numSamples);


// ****************************************************************************

#ifdef __cplusplus
} /* end extern "C" */
#endif

