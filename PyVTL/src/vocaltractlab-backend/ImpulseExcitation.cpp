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

#include "ImpulseExcitation.h"

static const int DURATION_SAMPLES = SAMPLING_RATE / 2;    // Half a second


// ****************************************************************************
/// Constructor.
// ****************************************************************************

ImpulseExcitation::ImpulseExcitation()
{
  tube = new Tube();
  impulseSection = Tube::FIRST_PHARYNX_SECTION;
}


// ****************************************************************************
/// Destructor.
// ****************************************************************************

ImpulseExcitation::~ImpulseExcitation()
{
  delete tube;
}


// ****************************************************************************
/// Set the static tube shape and the section that gets the flow impulse at t=0.
// ****************************************************************************

void ImpulseExcitation::setup(const Tube &tube, int impulseSection, double impulseAmp_cm3_s)
{
  *this->tube = tube;
  this->tube->setGlottisArea(0.0);    // Close the glottis
  this->tube->setAspirationStrength(Tube::DEFAULT_ASPIRATION_STRENGTH_DB);

  this->impulseSection = impulseSection;
  this->impulseAmp_cm3_s = impulseAmp_cm3_s;

  resetSequence();
}


// ****************************************************************************
/// Returns the tube at the current position.
// ****************************************************************************

void ImpulseExcitation::getTube(Tube &tube)
{
  tube = *this->tube;
}


// ****************************************************************************
/// Returns a flow impulse once at position zero (t=0) and no flow source
/// otherwise.
// ****************************************************************************

void ImpulseExcitation::getFlowSource(double &flow_cm3_s, int &section)
{
  if (pos == 0)
  {
    flow_cm3_s = impulseAmp_cm3_s;
    section = impulseSection;
  }
  else
  {
    flow_cm3_s = 0.0;
    section = -1;
  }
}


// ****************************************************************************
/// There exists no pressure source.
// ****************************************************************************

void ImpulseExcitation::getPressureSource(double &pressure_dPa, int &section)
{
  pressure_dPa = 0.0;
  section = -1;
}


// ****************************************************************************
/// Reset the position.
// ****************************************************************************

void ImpulseExcitation::resetSequence()
{
  pos = 0;
}


// ****************************************************************************
/// Increments the position.
// ****************************************************************************

void ImpulseExcitation::incPos(const double pressure_dPa[])
{
  pos++;
}


// ****************************************************************************
/// Returns the fixed duration in samples.
// ****************************************************************************

int ImpulseExcitation::getDuration_pt()
{
  return DURATION_SAMPLES;
}


// ****************************************************************************
/// Returns the current position.
// ****************************************************************************

int ImpulseExcitation::getPos_pt()
{
  return pos;
}

// ****************************************************************************
