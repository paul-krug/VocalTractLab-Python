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

#ifndef __TL_MODEL_H__
#define __TL_MODEL_H__

#include "Matrix2x2.h"
#include "Signal.h"
#include "Dsp.h"
#include "Tube.h"
#include "Constants.h"


// ****************************************************************************
/// This class encapsulates the simulation of vocal tract acoustics in the 
/// frequency domain on the basis of the transmission line analogy with lumped
/// elements.
// ****************************************************************************

class TlModel
{
  // ************************************************************************
  // Public data.
  // ************************************************************************

public:
  enum RadiationType 
  { 
    NO_RADIATION, 
    PISTONINSPHERE_RADIATION, 
    PISTONINWALL_RADIATION, 
    PARALLEL_RADIATION,
    NUM_RADIATION_OPTIONS
  };

  enum SpectrumType  
  { 
    INPUT_IMPEDANCE, OUTPUT_IMPEDANCE, FLOW_SOURCE_TF, PRESSURE_SOURCE_TF, RADIATION
  };

  /// Options for the acoustic simulation.

  struct Options
  {
    RadiationType radiation;
    bool boundaryLayer;
    bool heatConduction;
    bool softWalls;
    bool hagenResistance;
    bool innerLengthCorrections;
    bool lumpedElements;
    bool paranasalSinuses;
    bool piriformFossa;
    bool staticPressureDrops;
  };

  Options options;
  Tube tube;

  // ************************************************************************
  // Public functions.
  // ************************************************************************

public:
  TlModel();
  
  void getImpulseResponseWindow(Signal *window, int length);
  void getImpulseResponse(Signal *impulseResponse, int lengthExponent);
  void getSpectrum(SpectrumType type, ComplexSignal *spectrum, int spectrumLength, int section);

  int getMostConstrictedSection();
  double getMeanFlow(double lungPressure_dPa);
  void setLungPressure(double lungPressure_dPa);
  void getFormants(double *formantFreq, double *formantBW, int &numFormants, 
    const int MAX_FORMANTS, bool &frictionNoise, bool &isClosure, bool &isNasal);

  static double getCircumference(double area);

  // ************************************************************************
  // Private data.
  // ************************************************************************

private:
  /// Maximal number of sampling points in the spectrum
  static const int    MAX_NUM_FREQ = 4096;
  static const double MIN_AREA_CM2;
  static const double MIN_FREQ_RAD;

  // Options and tube geometry used in the previous calculation
  Options prevOptions;
  Tube prevTube;

  /// The product of the tube section matrices within a branch.
  Matrix2x2 matrixProduct[Tube::NUM_SECTIONS][MAX_NUM_FREQ];

  bool resetCalculations;   ///< Must the calculations be reset, because some parameter has changed
  double f0;                ///< Current frequency resolution
  int numFreq;
  double lungPressure_dPa;   ///< Currently set lung pressure in Pa

  double discreteOmega[MAX_NUM_FREQ];
  ComplexValue mouthRadiationImpedance[MAX_NUM_FREQ];
  ComplexValue noseRadiationImpedance[MAX_NUM_FREQ];
  ComplexValue lungTerminationImpedance[MAX_NUM_FREQ];
  ComplexValue radiationCharacteristic[MAX_NUM_FREQ];


  // ************************************************************************
  // Private functions.
  // ************************************************************************

private:
  void prepareCalculations();

  ComplexValue getRadiationCharacteristic(double omega);
  ComplexValue getRadiationImpedance(double omega, double radiationArea_cm2);
  void getLumpedSectionImpedances(double omega, Tube::Section *ts, ComplexValue &Za, ComplexValue &Zb);
  Matrix2x2 getSectionMatrix(double omega, int section);
  ComplexValue getJunctionImpedance(double omega, double A1_cm2, double A2_cm2);

  ComplexValue getInputImpedance(int freqIndex, int section);
  ComplexValue getOutputImpedance(int freqIndex, int section);
  ComplexValue getPressureSourceTF(int freqIndex, int section);
  ComplexValue getFlowSourceTF(int freqIndex, int section);
};

#endif
