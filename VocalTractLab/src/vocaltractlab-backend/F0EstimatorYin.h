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

#ifndef __F0_ESTIMATOR_YIN_H__
#define __F0_ESTIMATOR_YIN_H__

#include <vector>
#include "Dsp.h"
#include "Signal.h"
#include "IirFilter.h"
#include "Constants.h"

using namespace std;


// ****************************************************************************
// This class estimates the fundamental frequency from a given signal.
// ****************************************************************************

class F0EstimatorYin
{
  // **************************************************************************
  // Public data.
  // **************************************************************************

public:
  static const int INTEGRATION_LENGTH = SAMPLING_RATE / 60;    // for min. F0 = 60 Hz
  static const int FRAME_LENGTH = 2*INTEGRATION_LENGTH - 1;
  // For maximal F0 of 800 Hz with a 25 ms integration window
  static const int MAX_PITCH_CANDIDATES = 32;   

  double differenceFunctionThreshold;
  double timeStep_s;

  struct FrameData
  {
    int numPitchCandidates;
    double pitchCandidateT0[MAX_PITCH_CANDIDATES];
    double pitchCandidateY[MAX_PITCH_CANDIDATES];
    /// Cummulated lowest costs up to this frame and candidate
    double lowestPathCost[MAX_PITCH_CANDIDATES];
    /// Back-pointer to the best candidate in the best path in the previous frame
    int bestPrevCandidate[MAX_PITCH_CANDIDATES];

    double rmsAmplitude;
    int zeroCrossings;
    int initialCandidate;     ///< Index in pitch candidates
    int finalCandidate;       ///< Index in pitch candidates
  };

  vector<FrameData> frames;

  // **************************************************************************
  // Public functions.
  // **************************************************************************

public:
  F0EstimatorYin();
  
  void init(Signal16 *signal, int firstRoiSample, int numRoiSamples);
  bool processChunk(int numChunkSamples);
  vector<double> finish();

  void getFrameSignal(Signal16 *signal, int centerPos, double *frame);
  void getFrameSignal(Signal *signal, int centerPos, double *frame);
  void calcNdf(double *frame, double *df, double *ndf);
  void getFrameData(double *frameSignal, double *df, double *ndf, FrameData &fd);
  double getBestLocalT0Estimate(double t_s);

  void findBestPitchPath();
  double getTransitionCost(int prevFrame, int prevCandidate, int currFrame, int currCandidate);
  double getLocalCost(int frameIndex, int candidateIndex);
  double getFinalF0(double t_s);

  void fitParabola(double *f, int rawTau, double &accurateTau, double &accurateY);
  void filterSignal(double *inputSignal, double *outputSignal, int N);

  // **************************************************************************
  // Private data.
  // **************************************************************************

private:
  // Internally we use a smaller time step to handle phase variations.
  static const double INTERNAL_TIME_STEP_S;

  Signal origSignal;
  Signal filteredSignal;
  IirFilter *filter;
  double hannWindow[FRAME_LENGTH];
  int firstRoiSample;
  int numRoiSamples;
  int firstChunkSample;


  // **************************************************************************
  // Private functions.
  // **************************************************************************

private:

};

#endif

// ****************************************************************************
