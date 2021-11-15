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

#ifndef _POLE_ZERO_PLAN_H_
#define _POLE_ZERO_PLAN_H_

#include "Signal.h"
#include <vector>

using namespace std;

// ****************************************************************************
// A native class containing the data and functions for a pole-zero plan.
// ****************************************************************************

class PoleZeroPlan
{
public:
  struct Location
  {
    double freq_Hz;  // Frequency
    double bw_Hz;    // Bandwidth
  };

  vector<Location> poles;
  vector<Location> zeros;
  bool higherPoleCorrection;
  int selectedPole;
  int selectedZero;

  // Functions ******************************************************

  PoleZeroPlan();
  void createExample();
  void sortLocations(vector<Location> &origList, vector<Location> &sortedList);
  void getPoleZeroSpectrum(ComplexSignal *spectrum, int spectrumLength, double upperFrequencyLimit);
  void getHigherPoleCorrection(ComplexSignal *spectrum, int spectrumLength, double effectiveLength_cm);
};

// ****************************************************************************

#endif
