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

#include "Synthesizer.h"


// ****************************************************************************
/// Constructor.
// ****************************************************************************

Synthesizer::Synthesizer()
{
  int i;

  glottis = NULL;
  vocalTract = NULL;
  tdsModel = NULL;

  outputFlow = new double[TDS_BUFFER_LENGTH];
  outputPressure = new double[TDS_BUFFER_LENGTH];

  outputPressureFilter.createChebyshev((double)SYNTHETIC_SPEECH_BANDWIDTH_HZ / (double)SAMPLING_RATE, false, 8);

  initialShapesSet = false;

  for (i = 0; i < Glottis::MAX_CONTROL_PARAMS; i++)
  {
    prevGlottisParams[i] = 0.0;
  }
}


// ****************************************************************************
/// Destructor.
// ****************************************************************************

Synthesizer::~Synthesizer()
{
  delete [] outputFlow;
  delete [] outputPressure;
}


// ****************************************************************************
/// Initializes the synthesizer with the given objects.
// ****************************************************************************

void Synthesizer::init(Glottis *glottis, VocalTract *vocalTract, TdsModel *tdsModel)
{
  this->glottis = glottis;
  this->vocalTract = vocalTract;
  this->tdsModel = tdsModel;

  // Reset the dynamic state of all models and clear the buffers.
  reset();
}


// ****************************************************************************
/// Reset all models for a new synthesis.
// ****************************************************************************

void Synthesizer::reset()
{
  glottis->resetMotion();
  tdsModel->resetMotion();

  outputPressureFilter.resetBuffers();

  initialShapesSet = false;

  int i;
  for (i = 0; i < TDS_BUFFER_LENGTH; i++)
  {
    outputFlow[i] = 0.0;
    outputPressure[i] = 0.0;
  }
}


// ****************************************************************************
/// Generate an incremental part of the signal with a duration of numSamples
/// during which the vocal tract and glottis shapes are interpolated between 
/// the previous shapes (passed to the previous call of this function) and the 
/// given new shapes (i.e., parameters).
/// In the first call of this function (after reset()) the states of the 
/// glottis and vocal tract are initialized, and no audio is generated.
/// The value range of the generated audio samples is [-1; +1].
// ****************************************************************************

void Synthesizer::add(double *newGlottisParams, double *newTractParams, 
  int numSamples, vector<double> &audio)
{
  if (vocalTract == NULL)
  {
    return;
  }

  static Tube tube;
  int i;

  // Calculate the vocal tract with the new parameters.

  for (i = 0; i < VocalTract::NUM_PARAMS; i++)
  {
    vocalTract->param[i].x = newTractParams[i];
  }
  vocalTract->calculateAll();

  // Transform the vocal tract model into a tube.
  vocalTract->getTube(&tube);

  // Synthesize the new audio samples based on the tube model.
  add(newGlottisParams, &tube, numSamples, audio);
}


// ****************************************************************************
/// Generate an incremental part of the signal with a duration of numSamples
/// during which the tube and glottis shapes are interpolated between 
/// the previous shapes (passed to the previous call of this function) and the 
/// given new shapes.
/// In the first call of this function (after reset()) the states of the 
/// glottis and tube are initialized, and no audio is generated.
/// The value range of the generated audio samples is [-1; +1].
// ****************************************************************************

void Synthesizer::add(double *newGlottisParams, Tube *newTube, 
  int numSamples, vector<double> &audio)
{
  int i, k;

  int numGlottisParams = (int)glottis->controlParam.size();

  if (initialShapesSet == false)
  {
    prevTube = *newTube;
    for (i = 0; i < numGlottisParams; i++)
    {
      prevGlottisParams[i] = newGlottisParams[i];
    }

    initialShapesSet = true;
    return;
  }

  // ****************************************************************

  if (numSamples < 1)
  {
    return;
  }

  audio.resize(numSamples);

  double ratio, ratio1;
  // Lengths and areas of the glottis sections.
  double length_cm[Tube::NUM_GLOTTIS_SECTIONS];
  double area_cm2[Tube::NUM_GLOTTIS_SECTIONS];
  bool filtering;
  double pressure_dPa[4];
  
  double totalFlow_cm3_s;
  double mouthFlow_cm3_s;
  double nostrilFlow_cm3_s;
  double skinFlow_cm3_s;

  // ****************************************************************
  // Run through all new samples.
  // ****************************************************************

  for (i = 0; i < numSamples; i++)
  {
    ratio = (double)i / (double)numSamples;
    ratio1 = 1.0 - ratio;

    // ****************************************************************
    // Interpolate the tube.
    // ****************************************************************

    tube.interpolate(&prevTube, newTube, ratio);

    // ****************************************************************
    // Interpolate the glottis geometry.
    // ****************************************************************

    for (k = 0; k < numGlottisParams; k++)
    {
      glottis->controlParam[k].x = ratio1 * prevGlottisParams[k] + ratio * newGlottisParams[k];
    }

    glottis->calcGeometry();
    glottis->getTubeData(length_cm, area_cm2);

    tube.setGlottisGeometry(length_cm, area_cm2);
    tube.setAspirationStrength(glottis->getAspirationStrength_dB());

    // ****************************************************************
    // Do the acoustic simulation.
    // ****************************************************************

    if (tdsModel->getSampleIndex() == 0)
    {
      filtering = false;
    }
    else
    {
      filtering = true;
    }

    tdsModel->setTube(&tube, filtering);
    tdsModel->setFlowSource(0.0, -1);
    tdsModel->setPressureSource(glottis->controlParam[Glottis::PRESSURE].x, Tube::FIRST_TRACHEA_SECTION);

    // Get the four relevant pressure values for the glottis model:
    // subglottal, lower glottis, upper glottis, supraglottal.

    pressure_dPa[0] = tdsModel->getSectionPressure(Tube::LAST_TRACHEA_SECTION);
    pressure_dPa[1] = tdsModel->getSectionPressure(Tube::LOWER_GLOTTIS_SECTION);
    pressure_dPa[2] = tdsModel->getSectionPressure(Tube::UPPER_GLOTTIS_SECTION);
    pressure_dPa[3] = tdsModel->getSectionPressure(Tube::FIRST_PHARYNX_SECTION);

    // Increment the time/sample number
    glottis->incTime(1.0 / (double)SAMPLING_RATE, pressure_dPa);

    totalFlow_cm3_s = tdsModel->proceedTimeStep(mouthFlow_cm3_s, nostrilFlow_cm3_s, skinFlow_cm3_s);

    int pos = tdsModel->getSampleIndex();
    k = pos & TDS_BUFFER_MASK;
    outputFlow[k] = totalFlow_cm3_s;
    outputPressure[k] = (outputFlow[k] - outputFlow[(k - 1) & TDS_BUFFER_MASK]) / tdsModel->timeStep;
    // Scale the output to the range [-1, +1].
    audio[i] = outputPressureFilter.getOutputSample(outputPressure[k]) * 1e-7;
  }

  // ****************************************************************

  prevTube = *newTube;
  for (i = 0; i < numGlottisParams; i++)
  {
    prevGlottisParams[i] = newGlottisParams[i];
  }
}


// ****************************************************************************
/// Static function that takes the samples of type double of the source signal
/// and copies them to the position startPosInTarget in the target signal with
/// samples of type signed short.
/// The value range [-1, +1] from the source signal is linearly scaled to the
/// range [-32768, 32767] of the target signal.
// ****************************************************************************

void Synthesizer::copySignal(vector<double> &sourceSignal, Signal16 &targetSignal,
  int startPosInTarget)
{
  int i;
  signed short value = 0;
  int length = (int)sourceSignal.size();

  for (i = 0; i < length; i++)
  {
    value = (signed short)(sourceSignal[i] * 32767.0);
    targetSignal.setValue(startPosInTarget + i, value);
  }
}


// ****************************************************************************
/// Synthesis of a complete gestural score (blocking synthesis).
// ****************************************************************************

void Synthesizer::synthesizeGesturalScore(GesturalScore *gesturalScore,
  TdsModel *tdsModel, vector<double> &audio, bool enableConsoleOutput)
{
  int i;
  vector<double> signalPart;
  Glottis *glottis = gesturalScore->glottis;
  VocalTract *vocalTract = gesturalScore->vocalTract;
  double tractParams[VocalTract::NUM_PARAMS];
  double glottisParams[Glottis::MAX_CONTROL_PARAMS];
  int scoreLength_pt = gesturalScore->getDuration_pt();
  int numChunks = (int)(scoreLength_pt / NUM_CHUNCK_SAMPLES) + 1;
  double pos_s = 0.0;

  Synthesizer *synth = new Synthesizer();

  // ****************************************************************
  // Save the current state of the glottis and the vocal tract.
  // ****************************************************************

  glottis->storeControlParams();
  vocalTract->storeControlParams();

  // Important: Calc. the parameter curves from the gestural score.
  gesturalScore->calcCurves();

  // ****************************************************************
  // Generate the audio signal in small sections of about 2.5 ms 
  // length each (= 110 samples).
  // ****************************************************************

  synth->init(glottis, vocalTract, tdsModel);
  audio.resize(0);

  if (enableConsoleOutput)
  {
    printf("Synthesis of gestural score startet ");
  }

  // Get the parameters right at the beginning.
  gesturalScore->getParams(0.0, tractParams, glottisParams);
  synth->add(glottisParams, tractParams, 0, signalPart);

  for (i = 1; i <= numChunks; i++)
  {
    if (((i & 63) == 0) && (enableConsoleOutput))
    {
      printf(".");
    }

    pos_s = (double)i * NUM_CHUNCK_SAMPLES / SAMPLING_RATE;
    gesturalScore->getParams(pos_s, tractParams, glottisParams);
    synth->add(glottisParams, tractParams, NUM_CHUNCK_SAMPLES, signalPart);
    audio.insert(audio.end(), signalPart.begin(), signalPart.end());
  }

  if (enableConsoleOutput)
  {
    printf(" finished.\n");
  }

  // ****************************************************************
  // Restore the current state of the glottis and the vocal tract.
  // ****************************************************************

  glottis->restoreControlParams();
  vocalTract->restoreControlParams();

  // Free the memory.

  delete synth;
}


// ****************************************************************************
/// Synthesis (blocking) of a tube sequence from the data in a TXT file.
// ****************************************************************************

bool Synthesizer::synthesizeTubeSequence(string fileName,
  Glottis *glottis, TdsModel *tdsModel, vector<double> &audio)
{
  // ****************************************************************
  // Open the file.
  // ****************************************************************

  ifstream file(fileName);

  if (file.is_open() == false)
  {
    printf("Error in synthesizeTubeSequence(): File could not be opened.\n");
    return false;
  }

  // ****************************************************************
  // Save the current state of the glottis.
  // ****************************************************************

  glottis->storeControlParams();

  // ****************************************************************
  // Read the 11 comment lines, the glottis model type, and the 
  // number of states.
  // ****************************************************************

  int i, k;
  string line;

  for (i = 0; i < 10; i++)
  {
    getline(file, line);
  }

  // Read the type of glottis model.
  
  getline(file, line);

  if (line != glottis->getName())
  {
    printf("Error in synthesizeTubeSequence(): The selected glottis model does not correspond to the one used in the file.\n");
    return false;
  }

  // Read the number of states.
  
  double temp = 0.0;
  getline(file, line);
  if (parseTextLine(line, 1, &temp) == false)
  {
    printf("Error in synthesizeTubeSequence(): Invalid number of states.\n");
    return false;
  }

  int numStates = (int)temp;

  // ****************************************************************
  // Generate the audio signal.
  // ****************************************************************

  vector<double> signalPart;
  int numGlottisParams = (int)glottis->controlParam.size();
  Tube tube;
  double incisorPos = 0.0;
  double velumOpening = 0.0;
  double tongueTipSideElevation = 0.0;

  double glottisParams[Glottis::MAX_CONTROL_PARAMS];
  double miscData[3];    // For incisor position, nasal port area, and tongue tip side elevation.
  double tubeArea[Tube::NUM_PHARYNX_MOUTH_SECTIONS];
  double tubeLength[Tube::NUM_PHARYNX_MOUTH_SECTIONS];
  double tubeArticulatorDouble[Tube::NUM_PHARYNX_MOUTH_SECTIONS];
  Tube::Articulator tubeArticulator[Tube::NUM_PHARYNX_MOUTH_SECTIONS];
  
  bool glottisParamsOk = true;
  bool miscDataOk = true;
  bool tubeAreaOk = true;
  bool tubeLengthOk = true;
  bool tubeArticulatorOk = true;
  bool stateOk = true;

  Synthesizer *synth = new Synthesizer();
  synth->init(glottis, NULL, tdsModel);     // Vocal tract model is not needed here (= NULL).
  audio.resize(0);

  for (i = 0; (i < numStates) && (stateOk); i++)
  {
    getline(file, line);
    glottisParamsOk = parseTextLine(line, numGlottisParams, glottisParams);

    getline(file, line);
    miscDataOk = parseTextLine(line, 3, miscData);

    getline(file, line);
    tubeAreaOk = parseTextLine(line, Tube::NUM_PHARYNX_MOUTH_SECTIONS, tubeArea);

    getline(file, line);
    tubeLengthOk = parseTextLine(line, Tube::NUM_PHARYNX_MOUTH_SECTIONS, tubeLength);

    getline(file, line);
    tubeArticulatorOk = parseTextLine(line, Tube::NUM_PHARYNX_MOUTH_SECTIONS, tubeArticulatorDouble);

    // **************************************************************

    if ((glottisParamsOk) && (miscDataOk) && (tubeAreaOk) && (tubeLengthOk) && (tubeArticulatorOk))
    {
      stateOk = true;

      // Convert the type of the tube articulators.
      for (k = 0; k < Tube::NUM_PHARYNX_MOUTH_SECTIONS; k++)
      {
        tubeArticulator[k] = (Tube::Articulator)((int)tubeArticulatorDouble[k]);
      }
      incisorPos = miscData[0];
      velumOpening = miscData[1];
      tongueTipSideElevation = miscData[2];

      tube.setPharynxMouthGeometry(tubeLength, tubeArea, tubeArticulator, incisorPos, tongueTipSideElevation);
      tube.setVelumOpening(velumOpening);

      synth->add(glottisParams, &tube, NUM_CHUNCK_SAMPLES, signalPart);
      audio.insert(audio.end(), signalPart.begin(), signalPart.end());
    }
    else
    {
      stateOk = false;
    }
  }

  // ****************************************************************
  // Restore the previous state of the glottis.
  // ****************************************************************

  glottis->restoreControlParams();

  // Close the file.

  file.close();

  // Delete the Synthesizer object.

  delete synth;

  if (stateOk)
  {
    printf("The tube sequence was synthesized with %d states.\n", numStates);
    return true;
  }
  else
  {
    printf("Error: The tube sequence file was corrupted.\n");
    return false;
  }
}


// ****************************************************************************
/// Synthesis (blocking) of a sequence of glottis and vocal tract parameters
/// from the data in a TXT file.
// ****************************************************************************

bool Synthesizer::synthesizeTractSequence(string fileName,
  Glottis *glottis, VocalTract *vocalTract, TdsModel *tdsModel, vector<double> &audio)
{
  // ****************************************************************
  // Open the file.
  // ****************************************************************

  ifstream file(fileName);

  if (file.is_open() == false)
  {
    printf("Error in synthesizeTractSequence(): File could not be opened.\n");
    return false;
  }

  // ****************************************************************
  // Save the current state of the glottis and the vocal tract.
  // ****************************************************************

  glottis->storeControlParams();
  vocalTract->storeControlParams();

  // ****************************************************************
  // Read the six comment lines, the type of the glottis model, and
  // the number of states.
  // ****************************************************************

  string line;

  getline(file, line);
  getline(file, line);
  getline(file, line);
  getline(file, line);
  getline(file, line);
  getline(file, line);

  // Read the type of glottis model.
  getline(file, line);

  if (line != glottis->getName())
  {
    printf("Error in synthesizeTractSequence(): The selected glottis model does not correspond to the one used in the file.\n");
    return false;
  }

  // Read the number of states.

  double temp = 0.0;
  getline(file, line);
  if (parseTextLine(line, 1, &temp) == false)
  {
    printf("Error in synthesizeTractSequence(): Invalid number of states.\n");
    return false;
  }

  int numStates = (int)temp;

  // ****************************************************************
  // Generate the audio signal.
  // ****************************************************************

  int i;
  vector<double> signalPart;
  double tractParams[VocalTract::NUM_PARAMS];
  double glottisParams[Glottis::MAX_CONTROL_PARAMS];
  int numGlottisParams = (int)glottis->controlParam.size();
  bool glottisParamsOk = false;
  bool tractParamsOk = false;
  bool stateOk = true;

  Synthesizer *synth = new Synthesizer();
  synth->init(glottis, vocalTract, tdsModel);
  audio.resize(0);

  for (i = 0; (i < numStates) && (stateOk); i++)
  {
    getline(file, line);
    glottisParamsOk = parseTextLine(line, numGlottisParams, glottisParams);

    getline(file, line);
    tractParamsOk = parseTextLine(line, VocalTract::NUM_PARAMS, tractParams);

    if ((glottisParamsOk) && (tractParamsOk))
    {
      stateOk = true;
      synth->add(glottisParams, tractParams, NUM_CHUNCK_SAMPLES, signalPart);
      audio.insert(audio.end(), signalPart.begin(), signalPart.end());
    }
    else
    {
      stateOk = false;
    }
  }

  // ****************************************************************
  // Restore the previous state of the vocal tract and glottis.
  // ****************************************************************

  glottis->restoreControlParams();
  vocalTract->restoreControlParams();

  // Close the file.

  file.close();

  // Delete the Synthesizer object.

  delete synth;

  if (stateOk)
  {
    printf("The tract sequence was synthesized with %d states.\n", numStates);
    return true;
  }
  else
  {
    printf("Error: The tract sequence file was corrupted.\n");
    return false;
  }
}


// ****************************************************************************
/// Synthesize (blocking) a static sound with the given glottis and vocal tract.
/// If (useConstantF0 == true), the synthesis will use the (constant) f0 that is
/// set in the glottis model. Otherwise, a new f0 trajectory is imposed.
// ****************************************************************************

void Synthesizer::synthesizeStaticPhoneme(Glottis *glottis, VocalTract *vocalTract,
  TdsModel *tdsModel, bool shortLength, bool useConstantF0, vector<double> &audio)
{
  int i;
  Synthesizer *synth = new Synthesizer();
  vector<double> signalPart;

  // ****************************************************************
  // Obtain the arrays with tract and glottis parameters.
  // ****************************************************************

  double tractParams[VocalTract::NUM_PARAMS];
  double glottisParams[Glottis::MAX_CONTROL_PARAMS];

  int numGlottisParams = (int)glottis->controlParam.size();
  for (i = 0; i < numGlottisParams; i++)
  {
    glottisParams[i] = glottis->controlParam[i].x;
  }
  
  for (i = 0; i < VocalTract::NUM_PARAMS; i++)
  {
    tractParams[i] = vocalTract->param[i].x;
  }

  // ****************************************************************
  // Save the current state of the glottis and the vocal tract.
  // ****************************************************************

  glottis->storeControlParams();
  vocalTract->storeControlParams();

  // ****************************************************************
  // Generate the audio signal in three sections.
  // ****************************************************************

  synth->init(glottis, vocalTract, tdsModel);
  audio.resize(0);

  // Pressure starts at 0 dPa.

  glottisParams[Glottis::PRESSURE] = 0.0;
  if (useConstantF0 == false)
  {
    glottisParams[Glottis::FREQUENCY] = 110;
  }
  synth->add(glottisParams, tractParams, 0, signalPart);
  audio.insert(audio.end(), signalPart.begin(), signalPart.end());

  // Pressure rises up to 800 Pa;

  for (i = 0; i < 10; i++)
  {
    double factor = 0.5 * (-cos(M_PI * i / 10.0) + 1.0);
    glottisParams[Glottis::PRESSURE] = 8000.0 * factor;    // in dPa
    synth->add(glottisParams, tractParams, (int)(0.005 * SAMPLING_RATE), signalPart);
    audio.insert(audio.end(), signalPart.begin(), signalPart.end());
  }

  // Stationary part.

  glottisParams[Glottis::PRESSURE] = 8000.0;    // in dPa
  if (useConstantF0 == false)
  {
    glottisParams[Glottis::FREQUENCY] = 100;
  }

  int numSamples = 0;
  if (shortLength)
  {
    numSamples = (int)(0.200 * SAMPLING_RATE);    // 300 ms
  }
  else
  {
    numSamples = (int)(0.400 * SAMPLING_RATE);    // 600 ms
  }
  synth->add(glottisParams, tractParams, numSamples, signalPart);
  audio.insert(audio.end(), signalPart.begin(), signalPart.end());

  // Pressure is falling back to zero.

  for (i = 0; i < 10; i++)
  {
    double factor = 0.5 * (cos(M_PI * (i + 1) / 10.0) + 1.0);
    glottisParams[Glottis::PRESSURE] = 8000.0 * factor;    // in dPa
    synth->add(glottisParams, tractParams, (int)(0.005 * SAMPLING_RATE), signalPart);
    audio.insert(audio.end(), signalPart.begin(), signalPart.end());
  }

  // Pressure stays zero until the impulse response of the vocal tract
  // completely decayed.

  glottisParams[Glottis::PRESSURE] = 0.0;    // in dPa
  synth->add(glottisParams, tractParams, (int)(0.030 * SAMPLING_RATE), signalPart);
  audio.insert(audio.end(), signalPart.begin(), signalPart.end());

  // ****************************************************************
  // Restore the previous state of the vocal tract and glottis.
  // ****************************************************************

  glottis->restoreControlParams();
  vocalTract->restoreControlParams();

  delete synth;
}


// ****************************************************************************
/// Generate a text file that contains a sequence of vocal fold and vocal
/// tract control parameters (in steps of 110 samples or about 2.5 ms)
/// from the given gestural score.
// ****************************************************************************

bool Synthesizer::gesturalScoreToTractSequenceFile(GesturalScore *gesturalScore, string fileName)
{
  int i, k;
  double tractParams[VocalTract::NUM_PARAMS];
  double glottisParams[Glottis::MAX_CONTROL_PARAMS];
  int numGlottisParams = (int)gesturalScore->glottis->controlParam.size();
  int scoreLength_pt = gesturalScore->getDuration_pt();
  int numStates = (int)(scoreLength_pt / NUM_CHUNCK_SAMPLES) + 2;
  double pos_s = 0.0;

  ofstream file;
  file.open(fileName);
  if (file.is_open() == false)
  {
    printf("Error in gesturalScoreToTractSequenceFile(): The file could not be opened!\n");
    return false;
  }

  // Write some header data into the file.

  file << "# The first two lines (below the comment lines) indicate the name of the vocal fold model and the number of states." << endl;
  file << "# The following lines contain the control parameters of the vocal folds and the vocal tract (states)" << endl;
  file << "# in steps of 110 audio samples (corresponding to about 2.5 ms for the sampling rate of 44100 Hz)." << endl;
  file << "# For every step, there is one line with the vocal fold parameters followed by" << endl;
  file << "# one line with the vocal tract parameters." << endl;
  file << "#" << endl;

  // Write the name of the glottis model.
  file << gesturalScore->glottis->getName() << endl;

  // Write the number of states.
  file << numStates << endl;

  // Important: Calc. the parameter curves from the gestural score.
  gesturalScore->calcCurves();

  // ****************************************************************
  // Write the vocal tract and glottis parameters to the file every
  // 110 samples (about 2.5 ms).
  // ****************************************************************

  for (i = 0; i < numStates; i++)
  {
    pos_s = (double)i * NUM_CHUNCK_SAMPLES / SAMPLING_RATE;
    gesturalScore->getParams(pos_s, tractParams, glottisParams);
    
    // Write the new parameters to the file
    
    for (k = 0; k < numGlottisParams; k++)
    {
      file << glottisParams[k] << " ";
    }
    file << endl;

    for (k = 0; k < VocalTract::NUM_PARAMS; k++)
    {
      file << tractParams[k] << " ";
    }
    file << endl;
  }

  // Close the file.
  file.close();

  return true;
}


// ****************************************************************************
/// Generate a text file that contains a sequence of vocal fold control 
/// parameters and enhanced area functions (in steps of 110 samples or about 
/// 2.5 ms) from the given gestural score.
// ****************************************************************************

bool Synthesizer::gesturalScoreToTubeSequenceFile(GesturalScore *gesturalScore, string fileName)
{
  int i, k;
  double tractParams[VocalTract::NUM_PARAMS];
  double glottisParams[Glottis::MAX_CONTROL_PARAMS];
  int numGlottisParams = (int)gesturalScore->glottis->controlParam.size();
  int scoreLength_pt = gesturalScore->getDuration_pt();
  int numStates = (int)(scoreLength_pt / NUM_CHUNCK_SAMPLES) + 2;
  double pos_s = 0.0;
  Tube tube;

  ofstream file;
  file.open(fileName);
  if (file.is_open() == false)
  {
    printf("Error in gesturalScoreToTubeSequenceFile(): The file could not be opened!\n");
    return false;
  }

  // ****************************************************************
  // Save the current state of the vocal tract.
  // ****************************************************************

  gesturalScore->vocalTract->storeControlParams();

  // ****************************************************************
  // Write some header data into the file.
  // ****************************************************************

  file << "# The first two lines (below the comment lines) indicate the name of the vocal fold model and the number of states." << endl;
  file << "# The following lines contain a sequence of states of the vocal folds and the tube geometry" << endl;
  file << "# in steps of 110 audio samples (corresponding to about 2.5 ms for the sampling rate of 44100 Hz)." << endl;
  file << "# Each state is represented in terms of five lines:" << endl;
  file << "# Line 1: glottis_param_0 glottis_param_1 ..." << endl;
  file << "# Line 2: incisor_position_in_cm, velo_pharyngeal_opening_in_cm^2, tongue_tip_side_elevation[-1...1]" << endl;
  file << "# Line 3: area0 area1 area2 area3 ... (Areas of the tube sections in cm^2 from glottis to mouth)" << endl;
  file << "# Line 4: length0 length1 length2 length3 ... (Lengths of the tube sections in cm from glottis to mouth)" << endl;
  file << "# Line 5: artic0 artic1 artic2 artic3 ... (Articulators of the tube sections between glottis and lips : 1 = tongue; 2 = lower incisors; 3 = lower lip; 4 = other)" << endl;
  file << "#" << endl;

  // Write the name of the glottis model.
  file << gesturalScore->glottis->getName() << endl;

  // Write the number of states.
  file << numStates << endl;

  // Important: Calc. the parameter curves from the gestural score.
  gesturalScore->calcCurves();

  // ****************************************************************
  // Write the vocal tract and glottis parameters to the file every
  // 110 samples (about 2.5 ms).
  // ****************************************************************

  printf("Writing the tube sequence file started ...");

  for (i = 0; i < numStates; i++)
  {
    if ((i & 63) == 0)
    {
      printf(".");
    }

    pos_s = (double)i * NUM_CHUNCK_SAMPLES / SAMPLING_RATE;
    gesturalScore->getParams(pos_s, tractParams, glottisParams);

    // Calculate the vocal tract with the new parameters.

    for (k = 0; k < VocalTract::NUM_PARAMS; k++)
    {
      gesturalScore->vocalTract->param[k].x = tractParams[k];
    }
    gesturalScore->vocalTract->calculateAll();

    // Transform the vocal tract model into a tube.
    gesturalScore->vocalTract->getTube(&tube);

    // Write the new parameters to the file

    for (k = 0; k < numGlottisParams; k++)
    {
      file << glottisParams[k] << " ";
    }
    file << endl;

    file << tube.teethPosition_cm << " " << tube.getVelumOpening_cm2() << " "
      << tube.tongueTipSideElevation << endl;

    for (k = 0; k < Tube::NUM_PHARYNX_MOUTH_SECTIONS; k++)
    {
      file << tube.pharynxMouthSection[k].area_cm2 << " ";
    }
    file << endl;

    for (k = 0; k < Tube::NUM_PHARYNX_MOUTH_SECTIONS; k++)
    {
      file << tube.pharynxMouthSection[k].length_cm << " ";
    }
    file << endl;

    for (k = 0; k < Tube::NUM_PHARYNX_MOUTH_SECTIONS; k++)
    {
      file << tube.pharynxMouthSection[k].articulator << " ";
    }
    file << endl;
  }

  // ****************************************************************
  // Restore the current state of the vocal tract.
  // ****************************************************************

  gesturalScore->vocalTract->restoreControlParams();

  // Close the file.
  file.close();

  printf("finished.\n");

  return true;
}


// ****************************************************************************
/// Parse a text line (string) and obtain the space-separated numeric values.
// ****************************************************************************

bool Synthesizer::parseTextLine(string line, int numValues, double *values)
{
  int i;
  bool ok = true;
  istringstream iss(line);

  for (i = 0; (i < numValues) && (ok); i++)
  {
    if (!(iss >> values[i]))
    {
      ok = false;
    }
  }

  return ok;
}

// ****************************************************************************
