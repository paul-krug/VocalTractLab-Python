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

#include "VocalTractLabApi.h"
#include "Dsp.h"
#include "SoundLib.h"
#include "AudioFile.h"
#include "Synthesizer.h"
#include "SegmentSequence.h"

#include "GeometricGlottis.h"
#include "TwoMassModel.h"
#include "TriangularGlottis.h"

#include "VocalTract.h"
#include "TdsModel.h"
#include "GesturalScore.h"
#include "XmlHelper.h"
#include "XmlNode.h"
#include "TlModel.h"

#include <iostream>
#include <fstream>

enum GlottisModel
{
  GEOMETRIC_GLOTTIS,
  TWO_MASS_MODEL,
  TRIANGULAR_GLOTTIS,
  NUM_GLOTTIS_MODELS
};

static Glottis *glottis[NUM_GLOTTIS_MODELS];
static int selectedGlottis;

static VocalTract *vocalTract = NULL;
static TdsModel *tdsModel = NULL;
static Synthesizer *synthesizer = NULL;
static Tube *tube = NULL;

static bool vtlApiInitialized = false;


#if defined(WIN32) && defined(_USRDLL) 

// ****************************************************************************
/// Windows entry point for the DLL.
// ****************************************************************************

// Windows Header Files
#include <windows.h>

BOOL APIENTRY DllMain( HANDLE hModule, DWORD  ul_reason_for_call, LPVOID lpReserved)
{
  switch (ul_reason_for_call)
	{
		case DLL_PROCESS_ATTACH:
		case DLL_THREAD_ATTACH:
		case DLL_THREAD_DETACH:
		case DLL_PROCESS_DETACH:
	  break;
  }
  return TRUE;
}

#endif  // WIN32 && _USRDLL


// ****************************************************************************
// Loads the VT anatomy and the configurations for the different glottis 
// models from a speaker file.
// This function is not visible in the interface.
// ****************************************************************************

bool vtlLoadSpeaker(const char *speakerFileName, VocalTract *vocalTract, 
  Glottis *glottis[], int &selectedGlottis)
{

  // ****************************************************************
  // Load the XML data from the speaker file.
  // ****************************************************************

  vector<XmlError> xmlErrors;
  XmlNode *rootNode = xmlParseFile(string(speakerFileName), "speaker", &xmlErrors);
  if (rootNode == NULL)
  {
    xmlPrintErrors(xmlErrors);
    return false;
  }

  // ****************************************************************
  // Load the data for the glottis models.
  // ****************************************************************

  // This may be overwritten later.
  selectedGlottis = GEOMETRIC_GLOTTIS;

  XmlNode *glottisModelsNode = rootNode->getChildElement("glottis_models");
  if (glottisModelsNode != NULL)
  {
    int i;
    XmlNode *glottisNode;

    for (i=0; (i < (int)glottisModelsNode->childElement.size()) && (i < NUM_GLOTTIS_MODELS); i++)
    {
      glottisNode = glottisModelsNode->childElement[i];
      if (glottisNode->getAttributeString("type") == glottis[i]->getName())
      {
        if (glottisNode->getAttributeInt("selected") == 1)
        {
          selectedGlottis = i;
        }
        if (glottis[i]->readFromXml(*glottisNode) == false)
        {
          printf("Error: Failed to read glottis data for glottis model %d!\n", i);
          delete rootNode;
          return false;
        }
      }
      else
      {
        printf("Error: The type of the glottis model %d in the speaker file is '%s' "
          "but should be '%s'!\n", i, 
          glottisNode->getAttributeString("type").c_str(), 
          glottis[i]->getName().c_str());

        delete rootNode;
        return false;
      }
    }
  }
  else
  {
    printf("Warning: No glottis model data found in the speaker file %s!\n", speakerFileName);
  }

  // Free the memory of the XML tree !
  delete rootNode;

  // ****************************************************************
  // Load the vocal tract anatomy and vocal tract shapes.
  // ****************************************************************

  try
  {
    vocalTract->readFromXml(string(speakerFileName));
    vocalTract->calculateAll();
  }
  catch (std::string st)
  {
    printf("%s\n", st.c_str());
    printf("Error reading the anatomy data from %s.\n", speakerFileName);
    return false;
  }

  return true;
}


// ****************************************************************************
// Init. the synthesis with the given speaker file name, e.g. "JD2.speaker".
// This function should be called before any other function of this API.
// Return values:
// 0: success.
// 1: Loading the speaker file failed.
// ****************************************************************************

int vtlInitialize(const char *speakerFileName)
{
  if (vtlApiInitialized)
  {
    vtlClose();
  }

  // ****************************************************************
  // Init the vocal tract.
  // ****************************************************************

  vocalTract = new VocalTract();
  vocalTract->calculateAll();

  // ****************************************************************
  // Init the list with glottis models
  // ****************************************************************

  glottis[GEOMETRIC_GLOTTIS] = new GeometricGlottis();
  glottis[TWO_MASS_MODEL] = new TwoMassModel();
  glottis[TRIANGULAR_GLOTTIS] = new TriangularGlottis();
  
  selectedGlottis = GEOMETRIC_GLOTTIS;

  bool ok = vtlLoadSpeaker(speakerFileName, vocalTract, glottis, selectedGlottis);

  if (ok == false)
  {
    int i;
    for (i = 0; i < NUM_GLOTTIS_MODELS; i++)
    {
      delete glottis[i];
    }
    delete vocalTract;

    printf("Error in vtlInitialize(): vtlLoadSpeaker() failed.\n");
    return 1;
  }

  // ****************************************************************
  // Init the object for the time domain simulation.
  // ****************************************************************

  tdsModel = new TdsModel();

  // ****************************************************************
  // Init the Synthesizer object.
  // ****************************************************************

  synthesizer = new Synthesizer();
  synthesizer->init(glottis[selectedGlottis], vocalTract, tdsModel);

  tube = new Tube();

  // We are now initialized!
  vtlApiInitialized = true;

  return 0;
}


// ****************************************************************************
// Clean up the memory and shut down the synthesizer.
// Return values:
// 0: success.
// 1: The API was not initialized.
// ****************************************************************************

int vtlClose()
{
  if (!vtlApiInitialized)
  {
    printf("Error: The API was not initialized.\n");
    return 1;
  }

  delete synthesizer;
  delete tdsModel;

  int i;
  for (i = 0; i < NUM_GLOTTIS_MODELS; i++)
  {
    delete glottis[i];
  }

  delete vocalTract;

  vtlApiInitialized = false;
  
  return 0;
}


// ****************************************************************************
// Switch to turn off/on the automatic calculation of the tongue root 
// parameters TRX and TRY.
//
// Return values:
// 0: success.
// 1: The API was not initialized.
// ****************************************************************************

int vtlCalcTongueRootAutomatically(bool automaticCalculation)
{
    if (!vtlApiInitialized)
    {
        printf("Error: The API was not initialized.\n");
        return 1;
    }

    vocalTract->anatomy.automaticTongueRootCalc = automaticCalculation;
    vocalTract->calculateAll();

    return 0;
}



// ****************************************************************************
// Returns the version of this API as a string that contains the compile data.
// Reserve at least 32 chars for the string.
// ****************************************************************************

void vtlGetVersion(char *version)
{
  strcpy(version, __DATE__);
}


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

int vtlGetConstants(int *audioSamplingRate, int *numTubeSections,
                    int *numVocalTractParams, int *numGlottisParams,
                    int *numAudioSamplesPerTractState, double *internalSamplingRate)
{
  if (!vtlApiInitialized)
  {
    printf("Error: The API has not been initialized.\n");
    return 1;
  }

  *audioSamplingRate = SAMPLING_RATE;
  *numTubeSections = Tube::NUM_PHARYNX_MOUTH_SECTIONS;
  *numVocalTractParams = VocalTract::NUM_PARAMS;
  *numGlottisParams = (int)glottis[selectedGlottis]->controlParam.size();
  *numAudioSamplesPerTractState = Synthesizer::NUM_CHUNCK_SAMPLES;
  *internalSamplingRate = (double)SAMPLING_RATE / (double)Synthesizer::NUM_CHUNCK_SAMPLES;


  return 0;
}


// ****************************************************************************
// Returns for each supra glottal parameter the minimum value, the maximum value,
// and the standard (default) value. Each array passed to this function must have at 
// least as many elements as the number of supra glottal parameters.
// The "names" string receives the names of the parameters separated
// by tabs. This string should have at least 10*numParams elements.
// The "descriptions" string receives the descriptions of the parameters separated
// by tabs. This string should have at least 100*numParams elements.
// The "units" string receives the names of the parameter units separated
// by tabs. This string should have at least 10*numParams elements.
//
// Function return value:
// 0: success.
// 1: The API has not been initialized.
// ****************************************************************************

int vtlGetTractParamInfo(char* names, char* descriptions, char* units,
    double* paramMin, double* paramMax, double* paramStandard)
{
    if (!vtlApiInitialized)
    {
        printf("Error: The API has not been initialized.\n");
        return 1;
    }

    int i;

    strcpy(names, "");
    strcpy(descriptions, "");
    strcpy(units, "");

    for (i = 0; i < VocalTract::NUM_PARAMS; i++)
    {
        strcat(names, vocalTract->param[i].name.c_str());
        strcat(descriptions, vocalTract->param[i].description.c_str());
        strcat(units, vocalTract->param[i].unit.c_str());
        if (i != VocalTract::NUM_PARAMS - 1)
        {
            strcat(names, "\t");
            strcat(descriptions, "\t");
            strcat(units, "\t");
        }

        paramMin[i] = vocalTract->param[i].min;
        paramMax[i] = vocalTract->param[i].max;
        paramStandard[i] = vocalTract->param[i].neutral;
    }

    return 0;
}


// ****************************************************************************
// Returns for each glottis model parameter the minimum value, the maximum value,
// and the standard (default) value. Each array passed to this function must have at 
// least as many elements as the number of glottis model parameters.
// The "names" string receives the names of the parameters separated
// by tabs. This string should have at least 10*numParams elements.
// The "descriptions" string receives the descriptions of the parameters separated
// by tabs. This string should have at least 100*numParams elements.
// The "units" string receives the names of the parameter units separated
// by tabs. This string should have at least 10*numParams elements.
//
// Function return value:
// 0: success.
// 1: The API has not been initialized.
// ****************************************************************************

int vtlGetGlottisParamInfo(char* names, char* descriptions, char* units,
    double* paramMin, double* paramMax, double* paramStandard)
{
    if (!vtlApiInitialized)
    {
        printf("Error: The API has not been initialized.\n");
        return 1;
    }

    int i;
    int numGlottisParams = (int)glottis[selectedGlottis]->controlParam.size();

    strcpy(names, "");
    strcpy(descriptions, "");
    strcpy(units, "");

    for (i = 0; i < numGlottisParams; i++)
    {
        strcat(names, glottis[selectedGlottis]->controlParam[i].name.c_str());
        strcat(descriptions, glottis[selectedGlottis]->controlParam[i].description.c_str());
        strcat(units, glottis[selectedGlottis]->controlParam[i].cgsUnit.c_str());
        if (i != numGlottisParams - 1)
        {
            strcat(names, "\t");
            strcat(descriptions, "\t");
            strcat(units, "\t");
        }

        paramMin[i] = glottis[selectedGlottis]->controlParam[i].min;
        paramMax[i] = glottis[selectedGlottis]->controlParam[i].max;
        paramStandard[i] = glottis[selectedGlottis]->controlParam[i].neutral;
    }

    return 0;
}


// ****************************************************************************
// Returns the sub-glottal parameters for the given shape as defined in the
// speaker file.
// The array passed to this function must have at least as many elements as 
// the number of glottis model parameters.
//
// Function return value:
// 0: success.
// 1: The API has not been initialized.
// 2: A shape with the given name does not exist.
// ****************************************************************************

int vtlGetGlottisParams(const char *shapeName, double *glottisParams)
{
    if (!vtlApiInitialized)
    {
        printf("Error: The API has not been initialized.\n");
        return 1;
    }

    int index = glottis[selectedGlottis]->getShapeIndex(string(shapeName));
    if (index == -1)
    {
        return 2;
    }

    int numGlottisParams = (int)glottis[selectedGlottis]->controlParam.size();
    int i;
    for (i = 0; i < numGlottisParams; i++)
    {
        glottisParams[i] = glottis[selectedGlottis]->shape[index].controlParam[i];
    }

    return 0;
}



// ****************************************************************************
// Returns the supra-glottal parameters for the given shape as defined in the
// speaker file.
// The array passed to this function must have at least as many elements as 
// the number of vocal tract model parameters.
//
// Function return value:
// 0: success.
// 1: The API has not been initialized.
// 2: A shape with the given name does not exist.
// ****************************************************************************

int vtlGetTractParams(const char *shapeName, double *tractParams)
{
    if (!vtlApiInitialized)
    {
        printf("Error: The API has not been initialized.\n");
        return 1;
    }

    int index = vocalTract->getShapeIndex(string(shapeName));
    if (index == -1)
    {
        return 2;
    }

    int i;
    for (i = 0; i < VocalTract::NUM_PARAMS; i++)
    {
        tractParams[i] = vocalTract->shapes[index].param[i];
    }

    return 0;
}


// ****************************************************************************
// Exports the vocal tract contours for the given vector of vocal tract
// parameters as a SVG file (scalable vector graphics).
//
// Function return value:
// 0: success.
// 1: The API has not been initialized.
// 2: Writing the SVG file failed.
// ****************************************************************************

int vtlExportTractSvg(double *tractParams, const char *fileName)
{
  if (!vtlApiInitialized)
  {
    printf("Error: The API has not been initialized.\n");
    return 1;
  }

  // Store the current control parameter values.
  vocalTract->storeControlParams();

  // Set the given vocal tract parameters.
  int i;
  for (i = 0; i < VocalTract::NUM_PARAMS; i++)
  {
    vocalTract->param[i].x = tractParams[i];
  }
  vocalTract->calculateAll();
  
  // Save the contour as SVG file.
  bool ok = vocalTract->exportTractContourSvg(string(fileName), false, false);
  
  // Restore the previous control parameter values and 
  // recalculate the vocal tract shape.

  vocalTract->restoreControlParams();
  vocalTract->calculateAll();

  if (ok)
  {
    return 0;
  }
  else
  {
    return 2;
  }
}


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

int vtlTractToTube(double *tractParams,
  double *tubeLength_cm, double *tubeArea_cm2, int *tubeArticulator,
  double *incisorPos_cm, double *tongueTipSideElevation, double *velumOpening_cm2)
{
  if (!vtlApiInitialized)
  {
    printf("Error: The API has not been initialized.\n");
    return 1;
  }

  // ****************************************************************
  // Store the current control parameter values.
  // ****************************************************************

  vocalTract->storeControlParams();

  // ****************************************************************
  // Set the given vocal tract parameters.
  // ****************************************************************

  int i;
  for (i = 0; i < VocalTract::NUM_PARAMS; i++)
  {
    vocalTract->param[i].x = tractParams[i];
  }

  // ****************************************************************
  // Get the tube for the new vocal tract shape.
  // ****************************************************************

  Tube tube;
  vocalTract->calculateAll();
  vocalTract->getTube(&tube);

  // ****************************************************************
  // Copy the tube parameters to the user arrays.
  // ****************************************************************

  Tube::Section *ts = NULL;
  for (i = 0; i < Tube::NUM_PHARYNX_MOUTH_SECTIONS; i++)
  {
    ts = &tube.pharynxMouthSection[i];

    tubeLength_cm[i] = ts->length_cm;
    tubeArea_cm2[i] = ts->area_cm2;
    tubeArticulator[i] = ts->articulator;
  }

  *incisorPos_cm = tube.teethPosition_cm;
  *tongueTipSideElevation = tube.tongueTipSideElevation;
  *velumOpening_cm2 = tube.getVelumOpening_cm2();

  // ****************************************************************
  // Restore the previous control parameter values and 
  // recalculate the vocal tract shape.
  // ****************************************************************

  vocalTract->restoreControlParams();
  vocalTract->calculateAll();

  return 0;
}


// ****************************************************************************
// Returns the default options for the transfer function calculation. 
// 
// Parameters out:
// o opts: A struct containing the default values for the options available for
// the transfer function calculation.
//
// Function return value:
// 0: success.
// 1: The API has not been initialized.
// ****************************************************************************

int vtlGetDefaultTransferFunctionOptions(TransferFunctionOptions* opts)
{
    opts->radiationType = PARALLEL_RADIATION;
    opts->boundaryLayer = true;
    opts->heatConduction = false;
    opts->softWalls = true;
    opts->hagenResistance = false;
    opts->lumpedElements = true;
    opts->innerLengthCorrections = false;
    opts->paranasalSinuses = true;
    opts->piriformFossa = true;
    opts->staticPressureDrops = true;
    opts->spectrumType = SPECTRUM_UU;
    return 0;
}

// ****************************************************************************
// Calculates the transfer function of the vocal tract between 
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
// o opts: The options to use for the transfer function calculation. If NULL 
//     is passed, the default options will be used (see 
//     vtlGetDefaultTransferFunctionOptions()).
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

int vtlGetTransferFunction(double* tractParams, int numSpectrumSamples, TransferFunctionOptions* opts, double* magnitude, double* phase_rad)
{
    if (!vtlApiInitialized)
    {
        printf("Error: The API has not been initialized.\n");
        return 1;
    }

    int i;
    ComplexSignal s;

    if (numSpectrumSamples < 16)
    {
        numSpectrumSamples = 16;
    }

    // Calculate the vocal tract shape from the vocal tract parameters.

    for (i = 0; i < VocalTract::NUM_PARAMS; i++)
    {
        vocalTract->param[i].x = tractParams[i];
    }
    vocalTract->calculateAll();

    // Calculate the transfer function.

    TlModel* tlModel = new TlModel();

    // Set the options
    TlModel::Options tlOpts;
    if (opts == NULL)
    {
      TransferFunctionOptions tfOpts;
      vtlGetDefaultTransferFunctionOptions(&tfOpts);
      opts = &tfOpts;
    }
    tlOpts.boundaryLayer = opts->boundaryLayer;
    tlOpts.hagenResistance = opts->hagenResistance;
    tlOpts.heatConduction = opts->heatConduction;
    tlOpts.innerLengthCorrections = opts->innerLengthCorrections;
    tlOpts.lumpedElements = opts->lumpedElements;
    tlOpts.paranasalSinuses = opts->paranasalSinuses;
    tlOpts.piriformFossa = opts->piriformFossa;
    tlOpts.radiation = (TlModel::RadiationType)opts->radiationType;
    tlOpts.softWalls = opts->softWalls;

    tlModel->options = tlOpts;
    vocalTract->getTube(&tlModel->tube);
    tlModel->tube.setGlottisArea(0.0);

    tlModel->getSpectrum(TlModel::FLOW_SOURCE_TF, &s, numSpectrumSamples, Tube::FIRST_PHARYNX_SECTION);

    if (opts->spectrumType == SPECTRUM_PU)
    {
        ComplexSignal radiationSpectrum(0);
        tlModel->getSpectrum(TlModel::RADIATION, &radiationSpectrum, numSpectrumSamples, 0);
        s *= radiationSpectrum;
    }

    // Separate the transfer function into magnitude and phase.
    for (i = 0; i < numSpectrumSamples; i++)
    {
        magnitude[i] = s.getMagnitude(i);
        phase_rad[i] = s.getPhase(i);
    }

    delete tlModel;

    return 0;
}


// ****************************************************************************
// Calculates the real limited tract params (the ones that are actually used
// in the synthesis) from a given arbitrary set of tract parameters
//
// Parameters:
// o inTractParams (in): Is a vector of vocal tract parameters with 
//     numVocalTractParams elements.
// o outTractParams (out): Is a vector of vocal tract parameters with 
//     numVocalTractParams elements.
//
// Function return value:
// 0: success.
// 1: The API has not been initialized.
// ****************************************************************************

int vtlInputTractToLimitedTract(double* inTractParams, double* outTractParams)
{
    if (!vtlApiInitialized)
    {
        printf("Error: The API has not been initialized.\n");
        return 1;
    }

    // Calculate the vocal tract shape from the vocal tract parameters.
    int i;
    for (i = 0; i < VocalTract::NUM_PARAMS; i++)
    {
        vocalTract->param[i].x = inTractParams[i];
    }
    vocalTract->calculateAll();

    for (i = 0; i < VocalTract::NUM_PARAMS; i++)
    {
        outTractParams[i] = vocalTract->param[i].limitedX;
    }

    return 0;
}


// ****************************************************************************
// Resets the time-domain synthesis of continuous speech (using the functions
// vtlSynthesisAddTube() or vtlSynthesisAddTract()). This function must be 
// called every time you start a new synthesis.
//
// Function return value:
// 0: success.
// 1: The API has not been initialized.
// ****************************************************************************

int vtlSynthesisReset()
{
  if (!vtlApiInitialized)
  {
    printf("Error: The API has not been initialized.\n");
    return 1;
  }

  synthesizer->reset();
  tube->resetDynamicPart();

  return 0;
}


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

int vtlSynthesisAddTube(int numNewSamples, double* audio,
  double* tubeLength_cm, double* tubeArea_cm2, int* tubeArticulator,
  double incisorPos_cm, double velumOpening_cm2, double tongueTipSideElevation,
  double* newGlottisParams)
{
  if (!vtlApiInitialized)
  {
    printf("Error: The API has not been initialized.\n");
    return 1;
  }

  Tube::Articulator articulator[Tube::NUM_PHARYNX_MOUTH_SECTIONS];
  int i;
  for (i = 0; i < Tube::NUM_PHARYNX_MOUTH_SECTIONS; i++)
  {
    articulator[i] = (Tube::Articulator)tubeArticulator[i];
  }

  // Set the properties of the target tube.

  tube->setPharynxMouthGeometry(tubeLength_cm, tubeArea_cm2, articulator, 
    incisorPos_cm, tongueTipSideElevation);
  tube->setVelumOpening(velumOpening_cm2);
  // The aspiration strength will be set based on the glottis parameters
  // in synthesizer->add(...) below.
  tube->setAspirationStrength(0.0);

  // Synthesize the speech signal part.

  vector<double> audioVector;
  synthesizer->add(newGlottisParams, tube, numNewSamples, audioVector);

  if ((int)audioVector.size() != numNewSamples)
  {
    printf("Error in vtlSynthesisAddTube(): Number of audio samples is wrong.\n");
    return 2;
  }

  // Copy the audio samples in the given buffer.

  for (i = 0; i < numNewSamples; i++)
  {
    audio[i] = audioVector[i];
  }

  return 0;
}


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

int vtlSynthesisAddTract(int numNewSamples, double *audio,
  double *tractParams, double *glottisParams)
{
  if (!vtlApiInitialized)
  {
    printf("Error: The API has not been initialized.\n");
    return 1;
  }

  vector<double> audioVector;
  synthesizer->add(glottisParams, tractParams, numNewSamples, audioVector);

  if ((int)audioVector.size() != numNewSamples)
  {
    printf("Error in vtlSynthesisAddTube(): Number of audio samples is wrong.\n");
    return 2;
  }

  // Copy the audio samples in the given buffer.

  int i;
  for (i = 0; i < numNewSamples; i++)
  {
    audio[i] = audioVector[i];
  }

  return 0;
}


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

int vtlSynthBlock(double *tractParams, double *glottisParams,
  int numFrames, int frameStep_samples, double *audio, bool enableConsoleOutput)
{
  if (!vtlApiInitialized)
  {
    printf("Error: The API has not been initialized.\n");
    return 1;
  }

  int i;
  int samplePos = 0;
  int numGlottisParams = (int)glottis[selectedGlottis]->controlParam.size();

  if (enableConsoleOutput)
  {
    printf("Block synthesis in progress ...");
  }

  vtlSynthesisReset();

  for (i = 0; i < numFrames; i++)
  {
    if (i == 0)
    {
      // Only set the initial state of the vocal tract and glottis without generating audio.
      vtlSynthesisAddTract(0, &audio[0],
        &tractParams[i*VocalTract::NUM_PARAMS], &glottisParams[i*numGlottisParams]);
    }
    else
    {
      vtlSynthesisAddTract(frameStep_samples, &audio[samplePos],
        &tractParams[i*VocalTract::NUM_PARAMS], &glottisParams[i*numGlottisParams]);
      samplePos += frameStep_samples;
    }

    if ((enableConsoleOutput != 0) && ((i % 20) == 0))
    {
      printf(".");
    }
  }

  if (enableConsoleOutput != 0)
  {
    printf(" finished\n");
  }

  return 0;
}


// ****************************************************************************
// Test function for this API.
// Audio should contain at least 44100 double values.
// Run this WITHOUT calling vtlInitialize() !
// ****************************************************************************

int vtlApiTest(const char *speakerFileName, double *audio, int *numSamples)
{
  int failed = vtlInitialize(speakerFileName);
  if (failed != 0)
  {
    printf("Error in  in vtlApiTest(): vtlInitialize() failed.\n");
    return 1;
  }

  char version[100];
  vtlGetVersion(version);
  printf("Compile date of the library: %s\n", version);

  int audioSamplingRate = -1;
  int numTubeSections = -1;
  int numVocalTractParams = -1;
  int numGlottisParams = -1;
  int numAudioSamplesPerTractState = -1;
  double internalSamplingRate = -1.0;

  vtlGetConstants(&audioSamplingRate, &numTubeSections, &numVocalTractParams, &numGlottisParams, &numAudioSamplesPerTractState, &internalSamplingRate);

  printf("Audio sampling rate = %d\n", audioSamplingRate);
  printf("Num. of tube sections = %d\n", numTubeSections);
  printf("Num. of vocal tract parameters = %d\n", numVocalTractParams);
  printf("Num. of glottis parameters = %d\n", numGlottisParams);

  char tractParamNames[50 * 32];
  char tractParamDescriptions[500 * 32];
  char tractParamUnits[50 * 32];
  double tractParamMin[50];
  double tractParamMax[50];
  double tractParamStandard[50];

  vtlGetTractParamInfo(tractParamNames, tractParamDescriptions, tractParamUnits, tractParamMin, tractParamMax, tractParamStandard);

  char glottisParamNames[50 * 32];
  char glottisParamDescriptions[500 * 32];
  char glottisParamUnits[50 * 32];
  double glottisParamMin[50];
  double glottisParamMax[50];
  double glottisParamStandard[50];

  vtlGetGlottisParamInfo(glottisParamNames, glottisParamDescriptions, glottisParamUnits, glottisParamMin, glottisParamMax, glottisParamStandard);

  // ****************************************************************
  // Define two target tube shapes: one for /a/ and one for /i/.
  // ****************************************************************

  const int MAX_TUBES = 100;

  // These parameters are the same for both /i/ and /a/:
  double incisorPos_cm = 15.0;
  double velumOpening_cm2 = 0.0;
  double tongueTipSideElevation = 0.0;

  int i;

  // ****************************************************************
  // Define the tube for /i/.
  // ****************************************************************

  double tubeLength_cm_i[MAX_TUBES];
  double tubeArea_cm2_i[MAX_TUBES];
  int tubeArticulator_i[MAX_TUBES];

  for (i = 0; i < numTubeSections; i++)
  {
    // Full tube length is 16 cm.
    tubeLength_cm_i[i] = 16.0 / (double)numTubeSections;
    
    // Articulator is always the tongue (although not fully correct here)
    tubeArticulator_i[i] = 1;   // = tongue
    
    // Narrow mouth sections and wide pharynx sections
    if (i < numTubeSections / 2)
    {
      tubeArea_cm2_i[i] = 8.0;
    }
    else
    {
      tubeArea_cm2_i[i] = 2.0;
    }
  }

  // ****************************************************************
  // Define the tube for /a/.
  // ****************************************************************

  double tubeLength_cm_a[MAX_TUBES];
  double tubeArea_cm2_a[MAX_TUBES];
  int tubeArticulator_a[MAX_TUBES];

  for (i = 0; i < numTubeSections; i++)
  {
    // Full tube length is 16 cm.
    tubeLength_cm_a[i] = 16.0 / (double)numTubeSections;

    // Articulator is always the tongue (although not fully correct here)
    tubeArticulator_a[i] = 1;   // = tongue

    // Narrow mouth sections and wide pharynx sections
    if (i < numTubeSections / 2)
    {
      tubeArea_cm2_a[i] = 0.3;
    }
    else
    {
      tubeArea_cm2_a[i] = 8.0;
    }
  }

  // ****************************************************************
  // Set glottis parameters to default (neutral) values, which are
  // suitable for phonation.
  // ****************************************************************

  double glottisParams[Glottis::MAX_CONTROL_PARAMS];

  for (i = 0; i < numGlottisParams; i++)
  {
    glottisParams[i] = glottisParamStandard[i];
  }

  // **************************************************************************
  // Synthesize a transition from /a/ to /i/ to /a/.
  // **************************************************************************

  int numTotalSamples = 0;
  int numNewSamples = 0;

  vtlSynthesisReset();

  // Initialize with /a/ at 120 Hz.

  glottisParams[0] = 120.0;   // 120 Hz F0
  glottisParams[1] = 0.0;     // P_sub = 0 dPa.
  vtlSynthesisAddTube(0, audio, tubeLength_cm_a, tubeArea_cm2_a, tubeArticulator_a,
    incisorPos_cm, velumOpening_cm2, tongueTipSideElevation, 
    glottisParams);

  // Make 0.2 s transition to /i/ at 100 Hz.

  glottisParams[0] = 100.0;   // 100 Hz F0
  glottisParams[1] = 8000.0;  // P_sub = 8000 dPa.
  numNewSamples = (int)(0.2*audioSamplingRate);
  printf("Adding %d samples...\n", numNewSamples);

  vtlSynthesisAddTube(numNewSamples, &audio[numTotalSamples], tubeLength_cm_i, tubeArea_cm2_i, tubeArticulator_i,
    incisorPos_cm, velumOpening_cm2, tongueTipSideElevation,
    glottisParams);
  numTotalSamples += numNewSamples;

  // Make 0.2 s transition to /a/ at 80 Hz.

  glottisParams[0] = 80.0;   // 80 Hz F0
  numNewSamples = (int)(0.2*audioSamplingRate);
  printf("Adding %d samples...\n", numNewSamples);

  vtlSynthesisAddTube(numNewSamples, &audio[numTotalSamples], tubeLength_cm_a, tubeArea_cm2_a, tubeArticulator_a,
    incisorPos_cm, velumOpening_cm2, tongueTipSideElevation, glottisParams);
  numTotalSamples += numNewSamples;

  printf("Done.\n");

  *numSamples = numTotalSamples;

  // **************************************************************************
  // Clean up and close the VTL synthesis.
  // **************************************************************************

  vtlClose();

  return 0;
}


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

int vtlSegmentSequenceToGesturalScore(const char *segFileName, const char *gesFileName, bool enableConsoleOutput)
{
  if (!vtlApiInitialized)
  {
    printf("Error: The API has not been initialized.\n");
    return 1;
  }

  // Create and load the segment sequence file.
  
  SegmentSequence *segmentSequence = new SegmentSequence();
  if (segmentSequence->readFromFile(string(segFileName)) == false)
  {
    delete segmentSequence;
    printf("Error in vtlSegmentSequenceToGesturalScore(): Segment sequence file could not be loaded.\n");
    return 2;
  }

  // Create and save the gestural score.

  GesturalScore *gesturalScore = new GesturalScore(vocalTract, glottis[selectedGlottis]);
  gesturalScore->createFromSegmentSequence(segmentSequence, enableConsoleOutput);
  if (gesturalScore->saveGesturesXml(string(gesFileName)) == false)
  {
    delete segmentSequence;
    delete gesturalScore;
    printf("Error in vtlSegmentSequenceToGesturalScore(): Gestural score file could not be saved.\n");
    return 3;
  }

  delete segmentSequence;
  delete gesturalScore;

  return 0;
}


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

int vtlGesturalScoreToAudio(const char *gesFileName, const char *wavFileName,
  double *audio, int *numSamples, bool enableConsoleOutput)
{
  if (!vtlApiInitialized)
  {
    printf("Error: The API has not been initialized.\n");
    return 1;
  }

  int i;

  // ****************************************************************
  // Init and load the gestural score.
  // ****************************************************************

  GesturalScore *gesturalScore = new GesturalScore(vocalTract, glottis[selectedGlottis]);

  bool allValuesInRange = true;
  if (gesturalScore->loadGesturesXml(string(gesFileName), allValuesInRange) == false)
  {
    printf("Error in vtlGesturalScoreToAudio(): Loading the gestural score file failed!\n");
    delete gesturalScore;
    return 2;
  }

  if (allValuesInRange == false)
  {
    printf("Error in vtlGesturalScoreToAudio(): Some values in the gestural score are out of range!\n");
    delete gesturalScore;
    return 3;
  }

  // Important !!!
  gesturalScore->calcCurves();

  // ****************************************************************
  // Do the actual synthesis.
  // ****************************************************************

  vector<double> audioVector;
  Synthesizer::synthesizeGesturalScore(gesturalScore, tdsModel, audioVector, enableConsoleOutput);
  int numVectorSamples = (int)audioVector.size();

  // ****************************************************************
  // Copy the number of audio samples to the return value numSamples.
  // ****************************************************************

  if (numSamples != NULL)
  {
    *numSamples = numVectorSamples;
  }

  // ****************************************************************
  // Copy the synthesized signal into the return buffer audio.
  // ****************************************************************

  if (audio != NULL)
  {
    for (i = 0; i < numVectorSamples; i++)
    {
      audio[i] = audioVector[i];
    }
  }
   
  // ****************************************************************
  // Save the result as WAV file (if the name is not an empty string).
  // ****************************************************************

  if (wavFileName[0] != '\0')
  {
    AudioFile<double> audioFile;
    audioFile.setAudioBufferSize(1, numVectorSamples);
    audioFile.setBitDepth(16);
    audioFile.setSampleRate(SAMPLING_RATE);

    for (i = 0; i < numVectorSamples; i++)
    {
      audioFile.samples[0][i] = audioVector[i];
    }

    if (audioFile.save(string(wavFileName)) == false)
    {
      printf("Error in vtlGesturalScoreToAudio(): The WAV file could not be saved!\n");
      delete gesturalScore;
      return 4;
    }
  }

  // ****************************************************************
  // Free the memory and return.
  // ****************************************************************

  delete gesturalScore;
  return 0;
}


// ****************************************************************************
// This function directly converts a gestural score to a tract sequence file.
// The latter is a text file containing the parameters of the vocal fold and 
// vocal tract models in steps of about 2.5 ms.

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

int vtlGesturalScoreToTractSequence(const char* gesFileName, const char* tractSequenceFileName)
{
  if (!vtlApiInitialized)
  {
    printf("Error: The API has not been initialized.\n");
    return 1;
  }

  // ****************************************************************
  // Init and load the gestural score.
  // ****************************************************************

  GesturalScore* gesturalScore = new GesturalScore(vocalTract, glottis[selectedGlottis]);

  bool allValuesInRange = true;
  if (gesturalScore->loadGesturesXml(string(gesFileName), allValuesInRange) == false)
  {
    printf("Error in vtlGesturalScoreToTractSequence(): Loading the gestural score file failed!\n");
    delete gesturalScore;
    return 2;
  }

  if (allValuesInRange == false)
  {
    printf("Error in vtlGesturalScoreToTractSequence(): Some values in the gestural score are out of range!\n");
    delete gesturalScore;
    return 3;
  }

  // Important !!!
  gesturalScore->calcCurves();

  // ****************************************************************
  // Do the actual conversion.
  // ****************************************************************

  bool ok = Synthesizer::gesturalScoreToTractSequenceFile(gesturalScore, string(tractSequenceFileName));

  if (ok == false)
  {
    printf("Error in vtlGesturalScoreToTractSequence(): Saving the tract sequence file failed!\n");
    delete gesturalScore;
    return 4;
  }
   
  // ****************************************************************
  // Free the memory and return.
  // ****************************************************************

  delete gesturalScore;
  return 0;
}



// ****************************************************************************
// This function gets the duration from a gestural score.
// Parameters:
// o gesFileName (in): Name of the gestural score file.
// o audioFileDuration (out): The number of audio samples, the audio file would
//   have, if the gestural score was synthesized. This number can be slightly 
//   larger than the length of the gestural score because the audio is 
//   synthesized in chunks of a constant size. If not wanted, set to NULL.
// o gesFileDuration (out): The duration of the gestural score (in samples).
//   If not wanted, set to NULL.
//
// Function return value:
// 0: success.
// 1: The API was not initialized.
// 2: Loading the gestural score file failed.
// 3: Values in the gestural score file are out of range.
// ****************************************************************************

int vtlGetGesturalScoreDuration(const char* gesFileName, int* numAudioSamples, int* numGestureSamples)
{
    if (!vtlApiInitialized)
    {
        printf("Error: The API has not been initialized.\n");
        return 1;
    }



    // ****************************************************************
    // Init and load the gestural score.
    // ****************************************************************

    GesturalScore* gesturalScore = new GesturalScore(vocalTract, glottis[selectedGlottis]);
    static const int NUM_CHUNK_SAMPLES = 110;

    bool allValuesInRange = true;
    if (gesturalScore->loadGesturesXml(string(gesFileName), allValuesInRange) == false)
    {
        printf("Error in vtlGesturalScoreToGlottisSignals: Loading the gestural score file failed!\n");
        delete gesturalScore;
        return 2;
    }

    if (allValuesInRange == false)
    {
        printf("Error in vtlGesturalScoreToGlottisSignals: Some values in the gestural score are out of range!\n");
        delete gesturalScore;
        return 3;
    }

    // Important !!!
    gesturalScore->calcCurves();

    if (numGestureSamples != NULL)
    {
        *numGestureSamples = gesturalScore->getDuration_pt();
    }

    if (numAudioSamples != NULL)
    {
        *numAudioSamples = ( (int)( ( gesturalScore->getDuration_pt() ) / NUM_CHUNK_SAMPLES ) + 1 )  * NUM_CHUNK_SAMPLES;
    }


    return 0;
}


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

int vtlTractSequenceToAudio(const char* tractSequenceFileName, const char* wavFileName,
  double* audio, int* numSamples)
{
  if (!vtlApiInitialized)
  {
    printf("Error: The API has not been initialized.\n");
    return 1;
  }

  int i;
  vector<double> audioVector;

  bool ok = Synthesizer::synthesizeTractSequence(string(tractSequenceFileName),
    glottis[selectedGlottis], vocalTract, tdsModel, audioVector);

  if (ok == false)
  {
    printf("Error in vtlTractSequenceToAudio(): Synthesis of the tract sequence file failed.\n");
    return 2;
  }

  int numVectorSamples = (int)audioVector.size();

  // ****************************************************************
  // Copy the number of audio samples to the return value numSamples.
  // ****************************************************************

  if (numSamples != NULL)
  {
    *numSamples = numVectorSamples;
  }

  // ****************************************************************
  // Copy the synthesized signal into the return buffer audio.
  // ****************************************************************

  if (audio != NULL)
  {
    for (i = 0; i < numVectorSamples; i++)
    {
      audio[i] = audioVector[i];
    }
  }

  // ****************************************************************
  // Save the result as WAV file (if the file name is not empty).
  // ****************************************************************

  if (wavFileName[0] != '\0')
  {
    AudioFile<double> audioFile;
    audioFile.setAudioBufferSize(1, numVectorSamples);
    audioFile.setBitDepth(16);
    audioFile.setSampleRate(SAMPLING_RATE);

    for (i = 0; i < numVectorSamples; i++)
    {
      audioFile.samples[0][i] = audioVector[i];
    }

    if (audioFile.save(string(wavFileName)) == false)
    {
      printf("Error in vtlTractSequenceToAudio(): The WAV file could not be saved!\n");
      return 3;
    }
  }

  // ****************************************************************

  return 0;
}

// ****************************************************************************
