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

#include "StaticPhone.h"

// ****************************************************************************
/// Constructor.
// ****************************************************************************

StaticPhone::StaticPhone()
{
  useConstantF0 = true;

  tube = new Tube();
  defaultGlottis = new GeometricGlottis();
  setup(*tube, defaultGlottis, (int)(0.6*SAMPLING_RATE));
}


// ****************************************************************************
/// Destructor.
// ****************************************************************************

StaticPhone::~StaticPhone()
{
  delete tube;
}


// ****************************************************************************
/// Target values for F0 and lung pressure are given by the glottis parameters.
/// \param tube The tube geometry. The class makes a copy.
/// \param glottis Pointer to an existing glottis model.
/// \param duration_samples Duration of the phone in samples (at 44100 Hz).
// ****************************************************************************

void StaticPhone::setup(const Tube &tube, Glottis *glottis, const int duration_samples)
{
  *this->tube = tube;
  this->glottis = glottis;
  this->duration_samples = duration_samples;

  if (this->duration_samples < 0.4*SAMPLING_RATE)
  {
    this->duration_samples = (int)(0.4*SAMPLING_RATE);
  }

  double duration_s = this->duration_samples / (double)SAMPLING_RATE;

  // The time function for pulmonary pressure

  const int NUM_PRESSURE_NODES = 4;
  TimeFunction::Node p[NUM_PRESSURE_NODES] =
  {
    {0,              0.0},
    {0.04,           glottis->controlParam[Glottis::PRESSURE].x},
    {duration_s-0.2, glottis->controlParam[Glottis::PRESSURE].x},
    {duration_s,     0.0}
  };

  pressureTimeFunction.setNodes(p, NUM_PRESSURE_NODES);

  // The time function for F0; should be very close to the given
  // target F0 in the middle of the phone!!

  const int NUM_F0_NODES = 4;
  TimeFunction::Node f0[NUM_F0_NODES] =
  {
    {0.0,             0.9*glottis->controlParam[Glottis::FREQUENCY].x},
    {0.5*duration_s,  1.00*glottis->controlParam[Glottis::FREQUENCY].x},
    {0.75*duration_s, 0.8*glottis->controlParam[Glottis::FREQUENCY].x},
    {1.0*duration_s,  0.7*glottis->controlParam[Glottis::FREQUENCY].x}
  };

  f0TimeFunction.setNodes(f0, NUM_F0_NODES);

  // In the case the user wants a constant F0.
  constantF0 = glottis->controlParam[Glottis::FREQUENCY].x;

  // Reset the sequence state.
  resetSequence();
}


// ****************************************************************************
// ****************************************************************************

void StaticPhone::getTube(Tube &tube)
{
  tube = *this->tube;

  // ****************************************************************
  // Update the glottis sections of the tube.
  // ****************************************************************

  // Calc. the varying values for F0 und lung pressure.

  double t_s = (double)pos / (double)SAMPLING_RATE;
  glottis->controlParam[Glottis::PRESSURE].x  = pressureTimeFunction.getValue(t_s);
  if (useConstantF0)
  {
    glottis->controlParam[Glottis::FREQUENCY].x = constantF0;
  }
  else
  {
    glottis->controlParam[Glottis::FREQUENCY].x = f0TimeFunction.getValue(t_s);
  }

  // Set the glottis geometry.

  double length_cm[Tube::NUM_GLOTTIS_SECTIONS];
  double area_cm2[Tube::NUM_GLOTTIS_SECTIONS];

  glottis->calcGeometry();
  glottis->getTubeData(length_cm, area_cm2);

  tube.setGlottisGeometry(length_cm, area_cm2);
  tube.setAspirationStrength( glottis->getAspirationStrength_dB() );
}


// ****************************************************************************
/// There is no flow source -> return false.
// ****************************************************************************

void StaticPhone::getFlowSource(double &flow_cm3_s, int &section)
{
  section = -1;     // -1 means that there is no flow source
  flow_cm3_s = 0.0;
}


// ****************************************************************************
/// Returns the time-varying pulmonary pressure source parameters.
// ****************************************************************************

void StaticPhone::getPressureSource(double &pressure_dPa, int &section)
{
  section = Tube::FIRST_TRACHEA_SECTION;
  pressure_dPa = pressureTimeFunction.getValue( pos/(double)SAMPLING_RATE );
}


// ****************************************************************************
// ****************************************************************************

void StaticPhone::resetSequence()
{
  pos = 0;
  glottis->resetMotion();
}


// ****************************************************************************
// ****************************************************************************

void StaticPhone::incPos(const double pressure_dPa[])
{
  // Increment the time/sample number
  glottis->incTime(1.0/(double)SAMPLING_RATE, pressure_dPa);
  pos++;
}


// ****************************************************************************
// ****************************************************************************

int StaticPhone::getDuration_pt()
{
  return duration_samples;
}


// ****************************************************************************
// ****************************************************************************

int StaticPhone::getPos_pt()
{
  return pos;
}

// ****************************************************************************
