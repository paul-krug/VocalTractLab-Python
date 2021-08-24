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

#ifndef __GEOMETRIC_GLOTTIS_H__
#define __GEOMETRIC_GLOTTIS_H__

#include "Glottis.h"
#include "IirFilter.h"


// ****************************************************************************
// This class defines a kinematic/geometric glottis model similar to that
/// by Titze, but extended and improved by Birkholz.
// ****************************************************************************

class GeometricGlottis : public Glottis
{
  // **************************************************************************
  // Public data.
  // **************************************************************************

public:
  enum ControlParamIndex 
  { 
    // Frequency and lung pressure always must be the first two parameters
    FREQUENCY, 
    PRESSURE, 
    LOWER_END_X,
    UPPER_END_X,
    CHINK_AREA, 
    PHASE_LAG,
    RELATIVE_AMPLITUDE,
    DOUBLE_PULSING,
    PULSE_SKEWNESS,
    FLUTTER,
    ASPIRATION_STRENGTH,
    NUM_CONTROL_PARAMS  
  }; 

  enum StaticParamIndex
  {
    REST_THICKNESS,
    REST_LENGTH,
    REST_F0,
    CHINK_LENGTH,
    NUM_STATIC_PARAMS
  };

  enum DerivedParamIndex
  {
    LENGTH,
    THICKNESS,
    AMPLITUDE,
    LOWER_CORD_X,
    UPPER_CORD_X,
    LOWER_AREA,
    UPPER_AREA,
    CHINK_WIDTH,
    NUM_DERIVED_PARAMS
  };

  // **************************************************************************
  // Public functions.
  // **************************************************************************

public:
  GeometricGlottis();

  // Functions that overwrite the virtual functions in the base class.

  string getName();
  void resetMotion();
  void incTime(const double timeIncrement_s, const double pressure_dPa[]);
  void calcGeometry();
  void getTubeData(double *length_cm, double *area_cm2);
  int getApertureParamIndex();

  virtual double getAspirationStrength_dB();

  // **************************************************************************
  // Private data.
  // **************************************************************************

private:
  double phase;
  double time_s;      // absolute time
  double supraglottalPressure_dPa;

  IirFilter supraglottalPressureFilter;
};

#endif
