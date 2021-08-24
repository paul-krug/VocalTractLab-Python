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

#ifndef __TUBE_SEQUENCE_H__
#define __TUBE_SEQUENCE_H__

#include "Tube.h"

// ****************************************************************************
/// This virtual base class provides an INTERFACE to classes that generate
/// a complete area function (tube) including the glottis for each sample
/// of a speech signal.
/// This interface is implemented, for example, by classes that generate tube
/// sequences for gestural scores, for songs, for static vowels, for static 
/// fricatives, etc.
// ****************************************************************************

class TubeSequence
{
public:
  virtual ~TubeSequence() {}
  virtual void getTube(Tube &tube) = 0;
  virtual void getFlowSource(double &flow_cm3_s, int &section) = 0;
  virtual void getPressureSource(double &pressure_dPa, int &section) = 0;

  // These function set or get positions or durations in terms
  // of sample numbers at a sampling rate of 44 kHz.
  virtual void resetSequence() = 0;
  /// Requires four pressure values: subglottal, lower glottis, upper glottis, supraglottal
  virtual void incPos(const double pressure_dPa[]) = 0;
  virtual int getDuration_pt() = 0;
  virtual int getPos_pt() = 0;
};


#endif
