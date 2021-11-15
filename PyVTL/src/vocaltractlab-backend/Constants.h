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

#ifndef __PHYSICAL_CONSTANTS_H__
#define __PHYSICAL_CONSTANTS_H__

// ****************************************************************************
// Physical constants according to Flanagan (1965).
// ****************************************************************************

const double STATIC_PRESSURE_CGS = 1.013e6;   // deci-Pa = ubar
const double AMBIENT_DENSITY_CGS = 1.14e-3;   // g/cm^3
const double ADIABATIC_CONSTANT = 1.4;      
const double SOUND_VELOCITY_CGS = 3.5e4;      // cm/s
const double AIR_VISCOSITY_CGS = 1.86e-4;     // dyne-s/cm^2
const double SPECIFIC_HEAT_CGS = 0.24;        // cal/g-K
const double HEAT_CONDUCTION_CGS = 0.055e-3;  // cal/cm-s-K

const double CRITICAL_REYNOLDS_NUMBER = 1800.0;

// ****************************************************************************
// Constants for the synthesis.
// ****************************************************************************

const int SAMPLING_RATE = 44100;
const double SYNTHETIC_SPEECH_BANDWIDTH_HZ = 12000.0;

#endif

