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

#ifndef __VOWEL_LF_H__
#define __VOWEL_LF_H__

#include "TdsModel.h"
#include "TubeSequence.h"
#include "Tube.h"
#include "LfPulse.h"
#include "TimeFunction.h"
#include "Signal.h"

// ****************************************************************************
/// This class provides the tube sequence for the time-domain synthesis of
/// a static vowel with the Liljencrant-Fant model for glottal excitation.
// ****************************************************************************

class VowelLf : public TubeSequence
{
  // **************************************************************************
  // Public functions.
  // **************************************************************************

public:
  VowelLf();
  ~VowelLf();
  void setup(const Tube &tube, const LfPulse &lfPulse, const int duration_samples);

  // Overwritten functions of the interface class

  void getTube(Tube &tube);
  void getFlowSource(double &flow_cm3_s, int &section);
  void getPressureSource(double &pressure_dPa, int &section);

  void resetSequence();
  void incPos(const double pressure_dPa[]);
  int getDuration_pt();
  int getPos_pt();

  // **************************************************************************
  // Private data.
  // **************************************************************************

private:
  int duration_samples;
  TimeFunction f0TimeFunction;
  TimeFunction ampTimeFunction;
  int pos;
  Tube *tube;
  LfPulse lfPulse;
  Signal pulseForm;
  int pulseStartPos;
};

#endif
