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

#ifndef __VOICE_QUALITY_ESTIMATOR_H__
#define __VOICE_QUALITY_ESTIMATOR_H__

#include <vector>
#include "Dsp.h"
#include "Signal.h"

using namespace std;


// ****************************************************************************
// This class performs a voice quality estimation along the continuum from
// pressed over normal to breathy voice quality.
// Based on the paper: Kane, Gobl (2011): "Identifying regions of non-modal
// phonation using features of the wavelet transform", Interspeech 2011.
// The voice quality measure in this paper is called "PeakSlope".
// ****************************************************************************

class VoiceQualityEstimator
{
  // **************************************************************************
  // Public data.
  // **************************************************************************

public:
  struct Slice
  {
    // The peak amplitude in each frequency band within the slice.
    double peak250;
    double peak500;
    double peak1000;
    double peak2000;
    double peak4000;
    double peak8000;
  };

  static const double SLICE_STEP_S;
  static const double MIN_PEAK_SLOPE;
  static const double MAX_PEAK_SLOPE;

  double timeStep_s;
  Signal wavelet250;
  Signal wavelet500;
  Signal wavelet1000;
  Signal wavelet2000;
  Signal wavelet4000;
  Signal wavelet8000;

  // A vector of 10 ms slices.
  vector<Slice> slices;

  // **************************************************************************
  // Public functions.
  // **************************************************************************

public:
  VoiceQualityEstimator();

  void init(Signal16 *signal, int firstRoiSample, int numRoiSamples);
  bool processChunk(int numChunkSamples);
  vector<double> finish();

  void printData(int pos_pt);
  double calcPeakSlope(double centerTime_s, bool debug = false);
  void calcSlicePeaks(int sliceIndex);
  double getFilteredSample(int pos_pt, Signal *wavelet);

  // **************************************************************************
  // Private data.
  // **************************************************************************

private:
  Signal origSignal;
  int firstRoiSlice;
  int numRoiSlices;
  int nextSlice;

  // **************************************************************************
  // Private functions.
  // **************************************************************************

private:
  void calcWavelet(Signal &wavelet, int lengthFactor);

};

#endif

// ****************************************************************************
