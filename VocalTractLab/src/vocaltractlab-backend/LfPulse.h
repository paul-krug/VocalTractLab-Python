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

#ifndef __LF_PULSE_H__
#define __LF_PULSE_H__

#include "Signal.h"

// ****************************************************************************
/// This class represents the model of glottal flow introduced by Liljencrants
/// and Fant.
// ****************************************************************************

class LfPulse
{
  // **************************************************************************
  // Pulse parameters.
  // **************************************************************************

public:
  double F0;      // in Hz
  double AMP;     // in cm^3/s
  double OQ;      // [0, 1] Open quotient
  double SQ;      // [1, 2] Speed quotient
  double TL;      // [0, 0.2] Spectral tilt
  double SNR;	  // [0.0, 50.0] Signal to noise ratio

  // **************************************************************************
  // Public functions.
  // **************************************************************************

public:
  LfPulse();
  void resetParams();
  void getPulse(Signal& s, int numSamples, bool getDerivative);

  // **************************************************************************
  // Private functions.
  // **************************************************************************

private:
  double getEpsilon(double ta, double te);
  double getAlpha(double tp, double te, double ta, double epsilon);
  double getB(double AMP, double tp, double alpha);
};

#endif
