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

#ifndef __SYNTHESIZER_H__
#define __SYNTHESIZER_H__

#include "TdsModel.h"
#include "Tube.h"
#include "Glottis.h"
#include "VocalTract.h"
#include "GesturalScore.h"
#include "Dsp.h"
#include "IirFilter.h"
#include <vector>

using namespace std;

// ****************************************************************************
/// With this class, the user can incrementally synthesize a speech signal on 
/// the basis of a sequence of states of the vocal fold model and the vocal 
/// tract model or the tube shape.
/// In the time intervals between the provided states, the tube shapes are
/// linearly interpolated for the synthesis.
// ****************************************************************************

class Synthesizer
{
  // **************************************************************************
  // Public data.
  // **************************************************************************

public:
  // This is the default step size for the incremental synthesis 
  // corresponding to about 2.5 ms at our sampling rate of 44100 Hz.
  static const int NUM_CHUNCK_SAMPLES = 110;

  // **************************************************************************
  // Public functions.
  // **************************************************************************

public:
  Synthesizer();
  ~Synthesizer();

  void init(Glottis *glottis, VocalTract *vocalTract, TdsModel *tdsModel);
  void reset();
  void add(double *newGlottisParams, double *newTractParams, int numSamples, vector<double> &audio);
  void add(double *newGlottisParams, Tube *newTube, int numSamples, vector<double> &audio);

  // **************************************************************************

  static void copySignal(vector<double> &sourceSignal, Signal16 &targetSignal, 
    int startPosInTarget);

  static void synthesizeGesturalScore(GesturalScore *gesturalScore, 
    TdsModel *tdsModel, vector<double> &audio, bool enableConsoleOutput = true);

  static bool synthesizeTubeSequence(string fileName,
    Glottis *glottis, TdsModel *tdsModel, vector<double> &audio);

  static bool synthesizeTractSequence(string fileName,
    Glottis *glottis, VocalTract *vocalTract, TdsModel *tdsModel, vector<double> &audio);

  static void synthesizeStaticPhoneme(Glottis *glottis, VocalTract *vocalTract,
    TdsModel *tdsModel, bool shortLength, bool useConstantF0, vector<double> &audio);

  static bool gesturalScoreToTractSequenceFile(GesturalScore *gesturalScore, string fileName);
  static bool gesturalScoreToTubeSequenceFile(GesturalScore *gesturalScore, string fileName);


  // **************************************************************************
  // Private data.
  // **************************************************************************

private:
  // The three models are *not* owned by this class, but just used.
  Glottis *glottis;
  VocalTract *vocalTract;
  TdsModel *tdsModel;

  Tube prevTube;
  Tube tube;
  double prevGlottisParams[Glottis::MAX_CONTROL_PARAMS];

  static const int TDS_BUFFER_LENGTH = 256;
  static const int TDS_BUFFER_MASK = 255;

  double *outputFlow;
  double *outputPressure;
  IirFilter outputPressureFilter;
  bool initialShapesSet;

  // **************************************************************************
  // Private functions.
  // **************************************************************************

private:
  static bool parseTextLine(string line, int numValues, double *values);

};

#endif

// ****************************************************************************
