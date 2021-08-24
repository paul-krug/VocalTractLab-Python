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

#include "PoleZeroPlan.h"
#include "Dsp.h"
#include "TlModel.h"
#include "Constants.h"

// ****************************************************************************
// Constructor.
// ****************************************************************************

PoleZeroPlan::PoleZeroPlan()
{
  // Init the variables.
  higherPoleCorrection = true;
  poles.clear();
  zeros.clear();

  selectedPole = -1;
  selectedZero = -1;

  // Testing
  createExample();
}

// ****************************************************************************
/// Creates a simple example with three formants.
// ****************************************************************************

void PoleZeroPlan::createExample()
{
  Location p, z;

  poles.clear();
  zeros.clear();

  p.freq_Hz = 500.0;
  p.bw_Hz = 50.0;
  poles.push_back(p);

  p.freq_Hz = 1500.0;
  p.bw_Hz = 70.0;
  poles.push_back(p);

  p.freq_Hz = 2500.0;
  p.bw_Hz = 90.0;
  poles.push_back(p);

  p.freq_Hz = 3500.0;
  p.bw_Hz = 110.0;
  poles.push_back(p);

  p.freq_Hz = 4500.0;
  p.bw_Hz = 130.0;
  poles.push_back(p);

    p.freq_Hz = 5500.0;
  p.bw_Hz = 150.0;
  poles.push_back(p);

  z.freq_Hz = 4000.0;
  z.bw_Hz = 200.0;
  zeros.push_back(z);
}


// ****************************************************************************
/// Sorts the poles or zeros in the given list of locations with respect to
/// frequency and puts the result in sortedList.
// ****************************************************************************

void PoleZeroPlan::sortLocations(vector<Location> &origList, vector<Location> &sortedList)
{
  int i, k, winner;
  double temp;

  // Make a simple Bubble-sort
  sortedList = origList;
  int N = (int)sortedList.size();

  for (i=0; i < N-1; i++)
  {
    // The elements from index i to index sortedList.size()-1 are still unsorted
    // -> find the smallest of thes elements and put it to index i.
    winner = i;
    for (k=i+1; k < N; k++)
    {
      if (sortedList[k].freq_Hz < sortedList[winner].freq_Hz)
      {
        winner = k;
      }
    }

    // Exchange elements at the indices i and winner.
    if (i != k)
    {
      temp = sortedList[i].freq_Hz;
      sortedList[i].freq_Hz = sortedList[winner].freq_Hz;
      sortedList[winner].freq_Hz = temp;

      temp = sortedList[i].bw_Hz;
      sortedList[i].bw_Hz = sortedList[winner].bw_Hz;
      sortedList[winner].bw_Hz = temp;
    }

  }

}


// ****************************************************************************
// Create a spectrum from the PZ-plan.
// ****************************************************************************

void PoleZeroPlan::getPoleZeroSpectrum(ComplexSignal *spectrum, int spectrumLength, double upperFrequencyLimit)
{
  ComplexValue s;
  ComplexValue pole, conjugatePole;
  ComplexValue zero, conjugateZero;
  ComplexValue numerator;
  ComplexValue denominator;
  int i, k;
  double F0 = (double)SAMPLING_RATE / (double)spectrumLength;
  int numHarmonics = (int)(upperFrequencyLimit / F0);

  if (numHarmonics >= spectrumLength/2) 
  { 
    numHarmonics = spectrumLength/2 - 1; 
  }
  spectrum->reset(spectrumLength);

  // ****************************************************************
  // Run through all discrete frequencies.
  // ****************************************************************

  for (k=0; k < numHarmonics; k++)
  {
    s = ComplexValue(0, 2.0*M_PI*k*F0);
    numerator = ComplexValue(1.0, 0.0);
    denominator = ComplexValue(1.0, 0.0);

    for (i=0; i < (int)poles.size(); i++)
    {
      pole          = ComplexValue(-poles[i].bw_Hz*M_PI,  2.0*M_PI*poles[i].freq_Hz);
      conjugatePole = ComplexValue(-poles[i].bw_Hz*M_PI, -2.0*M_PI*poles[i].freq_Hz);
      numerator     = numerator*pole*conjugatePole;
      denominator   = denominator*(s-pole)*(s-conjugatePole);
    }

    for (i=0; i < (int)zeros.size(); i++)
    {
      zero          = ComplexValue(-zeros[i].bw_Hz*M_PI,  2.0*M_PI*zeros[i].freq_Hz);
      conjugateZero = ComplexValue(-zeros[i].bw_Hz*M_PI, -2.0*M_PI*zeros[i].freq_Hz);
      numerator     = numerator*(s-zero)*(s-conjugateZero);
      denominator   = denominator*zero*conjugateZero;
    }

    spectrum->setValue(k, numerator/denominator);
  }

  // Generate the negative frequencies in the spectrum.
  generateNegativeFrequencies(spectrum);
}


// ****************************************************************************
// Get the function for the higher pole correction.
// ****************************************************************************

void PoleZeroPlan::getHigherPoleCorrection(ComplexSignal *spectrum, int spectrumLength, double effectiveLength_cm)
{
  int i;
  double x;
  double value;
  double F1 = SOUND_VELOCITY_CGS / (4.0*effectiveLength_cm);
  double a, b;
  double f;
  double F0 = (double)SAMPLING_RATE / (double)spectrumLength;
  int numHarmonics = spectrumLength/2;

  spectrum->reset(spectrumLength);
  
  int numPoles = (int)poles.size();

  // ****************************************************************
  // Calc. the coefficients a and b.
  // ****************************************************************

  a = (M_PI*M_PI) / 8.0;
  b = (M_PI*M_PI*M_PI*M_PI) / 96.0;

  for (i=1; i <= numPoles; i++)
  {
    f = 2.0*(double)i - 1.0;
    a-= 1.0 / (f*f);
    b-= 1.0 / (f*f*f*f);
  }
  b*= 0.5;

  // ****************************************************************
  // Run through all discrete frequencies.
  // ****************************************************************

  if (numPoles == 0)
  {
    for (i=0; i < numHarmonics; i++) 
    { 
      spectrum->setValue(i, 1.0, 0.0); 
    }
  }
  else
  {
    for (i=0; i < numHarmonics; i++)
    {
      x = (F0*i) / F1;
      value = exp(a*x*x + b*x*x*x*x);
      spectrum->setValue(i, value, 0);
    }
  }

  // Generate the negative frequencies in the spectrum.
  generateNegativeFrequencies(spectrum);
}

// ****************************************************************************
