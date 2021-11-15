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

#include "VowelLf.h"

// ****************************************************************************
/// Constructor.
// ****************************************************************************

VowelLf::VowelLf()
{
  tube = new Tube();
  setup(*tube, lfPulse, (int)(0.6*SAMPLING_RATE));
}


// ****************************************************************************
/// Destructor.
// ****************************************************************************

VowelLf::~VowelLf()
{
  delete tube;
}


// ****************************************************************************
/// Target values for F0 an lung pressure are given by the glottis parameters.
// ****************************************************************************

void VowelLf::setup(const Tube &tube, const LfPulse &lfPulse, const int duration_samples)
{
  *this->tube = tube;

  // Close the glottis of the tube
  this->tube->setGlottisArea(0.0);

  this->lfPulse = lfPulse;
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
    {0.04,           lfPulse.AMP},
    {duration_s-0.2, lfPulse.AMP},
    {duration_s,     0.0}
  };

  ampTimeFunction.setNodes(p, NUM_PRESSURE_NODES);

  // The time function for F0

  const int NUM_F0_NODES = 3;
  TimeFunction::Node f0[NUM_F0_NODES] =
  {
    {0.0,   0.83*lfPulse.F0},
    {0.13*duration_s,  1.01*lfPulse.F0},
    {duration_s, 0.99*lfPulse.F0}
  };

  f0TimeFunction.setNodes(f0, NUM_F0_NODES);

  // Reset the sequence state.
  resetSequence();
}


// ****************************************************************************
// ****************************************************************************

void VowelLf::getTube(Tube &tube)
{
  tube = *this->tube;
}


// ****************************************************************************
/// There is no flow source -> return false.
// ****************************************************************************

void VowelLf::getFlowSource(double &flow_cm3_s, int &section)
{
  section = Tube::FIRST_PHARYNX_SECTION;
  flow_cm3_s = pulseForm.getValue(pos - pulseStartPos);
}


// ****************************************************************************
/// Returns the time-varying pulmonary pressure source parameters.
// ****************************************************************************

void VowelLf::getPressureSource(double &pressure_dPa, int &section)
{
  // Disable the pressure source.
  section = -1;
  pressure_dPa = 0.0;
}


// ****************************************************************************
// ****************************************************************************

void VowelLf::resetSequence()
{
  pos = 0;

  double t_s = 0.0;
  lfPulse.F0  = f0TimeFunction.getValue(t_s);
  lfPulse.AMP = ampTimeFunction.getValue(t_s);
  int numSamples = (int)(SAMPLING_RATE / lfPulse.F0);
  lfPulse.getPulse(pulseForm, numSamples, false);
  pulseStartPos = 0;
}


// ****************************************************************************
// ****************************************************************************

void VowelLf::incPos(const double pressure_dPa[])
{
  pos++;

  // Get a new glottal pulse form ?
  if (pos >= pulseStartPos + pulseForm.N)
  {
    pulseStartPos = pos;

    double t_s = (double)pos / (double)SAMPLING_RATE;
    lfPulse.F0  = f0TimeFunction.getValue(t_s);
    lfPulse.AMP = ampTimeFunction.getValue(t_s);
    int numSamples = (int)(SAMPLING_RATE / lfPulse.F0);
    lfPulse.getPulse(pulseForm, numSamples, false);
  }
}


// ****************************************************************************
// ****************************************************************************

int VowelLf::getDuration_pt()
{
  return duration_samples;
}


// ****************************************************************************
// ****************************************************************************

int VowelLf::getPos_pt()
{
  return pos;
}

// ****************************************************************************
