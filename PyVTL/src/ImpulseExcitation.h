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

#ifndef __IMPULSE_EXCITATION_H__
#define __IMPULSE_EXCITATION_H__

#include "TdsModel.h"
#include "TubeSequence.h"
#include "Tube.h"

// ****************************************************************************
/// This class provides the tube sequence for a flow impulse excitation of the
/// tube.
// ****************************************************************************

class ImpulseExcitation : public TubeSequence
{
  // **************************************************************************
  // Public functions.
  // **************************************************************************

public:
  ImpulseExcitation();
  ~ImpulseExcitation();
  void setup(const Tube &tube, int impulseSection, double impulseAmp_cm3_s);

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
  int pos;
  Tube *tube;
  int impulseSection;
  double impulseAmp_cm3_s;
};

#endif
