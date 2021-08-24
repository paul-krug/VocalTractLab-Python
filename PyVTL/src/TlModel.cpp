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

#include "TlModel.h"
#include <cmath>

const double TlModel::MIN_AREA_CM2  = 0.01e-2; // = 0.01 mm^2
const double TlModel::MIN_FREQ_RAD = 0.0001;

static const double RS_FACTOR = sqrt(AMBIENT_DENSITY_CGS*AIR_VISCOSITY_CGS*0.5);
static const double LS_FACTOR = AMBIENT_DENSITY_CGS;
static const double CP_FACTOR = 1.0 / (AMBIENT_DENSITY_CGS*SOUND_VELOCITY_CGS*SOUND_VELOCITY_CGS);
static const double GP_FACTOR = ((ADIABATIC_CONSTANT - 1.0)*
  sqrt((0.5*HEAT_CONDUCTION_CGS) / (SPECIFIC_HEAT_CGS*AMBIENT_DENSITY_CGS))) /
  (AMBIENT_DENSITY_CGS*SOUND_VELOCITY_CGS*SOUND_VELOCITY_CGS);


// ****************************************************************************
/// Constructor.
// ****************************************************************************

TlModel::TlModel()
{
  // Initial options.

  options.radiation = PARALLEL_RADIATION;
  options.boundaryLayer       = true;
  options.heatConduction      = false;
  options.softWalls           = true;
  options.hagenResistance     = false;
  options.lumpedElements      = true;
  options.innerLengthCorrections = false;
  options.paranasalSinuses    = true;
  options.piriformFossa       = true;
  options.staticPressureDrops = true;

  prevOptions = options;
  prevTube = tube;

  resetCalculations = true;
  f0 = 20.0;
  numFreq = (int)(8000.0 / f0);
  lungPressure_dPa = 0.0;

  // Init the matrices.

  int i, k;

  for (i=0; i < Tube::NUM_SECTIONS; i++)
  {
    for (k=0; k < MAX_NUM_FREQ; k++) 
    { 
      matrixProduct[i][k].unitMatrix(); 
    }
  }

}


// ****************************************************************************
/// Calculates a window with the form of a half cosine-cycle between 0 and Pi/2.
/// \param window Pointer to the resulting signal
/// \param length Desired window length
// ****************************************************************************

void TlModel::getImpulseResponseWindow(Signal *window, int length)
{
  int i;

  window->reset(length);
  
  for (i=0; i < length; i++) 
  { 
    window->x[i] = 0.5 - 0.5*cos(M_PI + (M_PI*(double)i) / (double)(length - 1));
  }
}


// ****************************************************************************
/// Returns the impulse response of the vocal tract with the length 
/// N = 2^lengthExponent.
/// \param impulseResponse Pointer to the impulse response signal
/// \param lengthExponent The number of samples is 2^lengthExponent
// ****************************************************************************

void TlModel::getImpulseResponse(Signal *impulseResponse, int lengthExponent)
{
  int signalLength = 1 << lengthExponent;

  ComplexSignal flowSourceTF(signalLength);
  ComplexSignal radiationSpectrum(signalLength);
  Signal window(signalLength);
  
  impulseResponse->reset(signalLength);

  getSpectrum(FLOW_SOURCE_TF, &flowSourceTF, signalLength, Tube::FIRST_PHARYNX_SECTION);
  getSpectrum(RADIATION, &radiationSpectrum, signalLength, 0);

  flowSourceTF*= radiationSpectrum;
  complexIFFT(flowSourceTF, lengthExponent, true);
  getImpulseResponseWindow(&window, signalLength);

  for (int i=0; i < signalLength; i++)
  {
    impulseResponse->x[i] = flowSourceTF.re[i]*window.x[i];
  }
}


// ****************************************************************************
/// Calculates a spectrum of a specific type and with a given length (number
/// of sampling points).
/// \param type The type of spectrum (Transfer function, imput impedance, ...)
/// \param Pointer to the resulting complex spectrum
/// \param Number of spectrum points (including the negative frequencies)
/// \param section Section index, depending on the type of transfer function
// ****************************************************************************

void TlModel::getSpectrum(SpectrumType type, ComplexSignal *spectrum, int spectrumLength, int section)
{
  double TL_CUTOFF_FREQ = 22050;

  if(this->options.lumpedElements == true)
  {
    // For the cutoff frequency, we *cannot* take the frequency
    // SYNTHETIC_SPEECH_BANDWIDTH_HZ (=12000), because the calculations
    // with *lumped* network elements start to generate very noisy
    // spectra from about 11300 Hz on. Maybe the wavelength (= 3 cm)
    // is getting too short at that frequency, but I am not really sure.

    TL_CUTOFF_FREQ = 10000;
  }
  
  
  int i;
  ComplexValue v;
 
  double newF0 = (double)SAMPLING_RATE / (double)spectrumLength;
  int newNumFreq = (int)(TL_CUTOFF_FREQ / newF0);
  
  if (newNumFreq > MAX_NUM_FREQ) 
  { 
    newNumFreq = MAX_NUM_FREQ; 
  }
  if (newNumFreq >= spectrumLength / 2) 
  { 
    newNumFreq = spectrumLength / 2 - 1; 
  }

  if ((newF0 != f0) || (newNumFreq != numFreq)) 
  { 
    resetCalculations = true; 
  }
  
  f0 = newF0;
  numFreq = newNumFreq;

  // Did any of the options change ?

  if ((options.boundaryLayer        != prevOptions.boundaryLayer)  ||
      (options.hagenResistance      != prevOptions.hagenResistance) ||
      (options.heatConduction       != prevOptions.heatConduction) ||
      (options.innerLengthCorrections  != prevOptions.innerLengthCorrections) ||
      (options.paranasalSinuses     != prevOptions.paranasalSinuses) ||
      (options.piriformFossa        != prevOptions.piriformFossa) ||
      (options.lumpedElements       != prevOptions.lumpedElements) ||
      (options.radiation            != prevOptions.radiation) ||
      (options.softWalls            != prevOptions.softWalls) ||
      (options.staticPressureDrops  != prevOptions.staticPressureDrops))
  { 
    resetCalculations = true; 
  }

  // Did the tube geometry change ?

  if (tube != prevTube)
  {
    resetCalculations = true; 
  }

  // ****************************************************************

  if (resetCalculations) 
  { 
    prepareCalculations(); 
  }

  // ****************************************************************

  spectrum->reset(spectrumLength);

  for (i=0; i < numFreq; i++)
  {
    v = 0.0;
    if (type == INPUT_IMPEDANCE)    { v = getInputImpedance(i, section); } else
    if (type == OUTPUT_IMPEDANCE)   { v = getOutputImpedance(i, section); } else
    if (type == PRESSURE_SOURCE_TF) { v = getPressureSourceTF(i, section); } else
    if (type == FLOW_SOURCE_TF)     { v = getFlowSourceTF(i, section); } else
    if (type == RADIATION)          { v = radiationCharacteristic[i]; }

    spectrum->setValue(i, v);
  }

  // Set the values above the cutoff frequency to a very small value.
  const double SMALL_VALUE = 0.000000001;
  for (i = numFreq; i <= spectrumLength / 2; i++)
  {
    spectrum->setValue(i, SMALL_VALUE);
  }

  generateNegativeFrequencies(spectrum);    // Fill up the negative frequencies
}


// ****************************************************************************
/// Returns the index of the most constricted tube section in the vocal tract.
// ****************************************************************************

int TlModel::getMostConstrictedSection()
{
  int i;
  int k = Tube::FIRST_PHARYNX_SECTION;

  for (i = Tube::FIRST_PHARYNX_SECTION+1; i <= Tube::LAST_MOUTH_SECTION; i++)
  {
    if (tube.section[i]->area_cm2 <= tube.section[k]->area_cm2) 
    { 
      k = i; 
    }
  }
  return k;
}


// ****************************************************************************
/// Returns an approximation for the mean flow through the vocal tract
/// considering the pulmonary pressure and Bernoulli head pressure losses
/// at the glottal and a supraglottal constrictions.
/// \param lungPressure_dPa The pulmonary pressure in deci-Pa.
// ****************************************************************************

double TlModel::getMeanFlow(double lungPressure_dPa)
{
  double A_g = tube.section[ Tube::LOWER_GLOTTIS_SECTION ]->area_cm2;
  double A_c = tube.section[getMostConstrictedSection()]->area_cm2;

  if (A_g < MIN_AREA_CM2) 
  { 
    A_g = MIN_AREA_CM2; 
  }
  if (A_c < MIN_AREA_CM2) 
  { 
    A_c = MIN_AREA_CM2; 
  }

  setLungPressure(lungPressure_dPa);

  double flow = sqrt(2.0*lungPressure_dPa / (AMBIENT_DENSITY_CGS*(1.0/(A_g*A_g) + 1.0/(A_c*A_c))));
  return flow;
}


// ****************************************************************************
/// Sets a new value for the pulmonary pressure.
/// \param lungPressure_dPa The pulmonary pressure in deci-Pa.
// ****************************************************************************

void TlModel::setLungPressure(double lungPressure_dPa)
{
  if (lungPressure_dPa != this->lungPressure_dPa) 
  { 
    resetCalculations = true; 
  }
  this->lungPressure_dPa = lungPressure_dPa;
}


// ****************************************************************************
/// Returns information about the formants in the vocal tract transfer function
/// (frequency and bandwidth) together with additional information about
/// wheather the spectrum shows signs of nasality, a closure or a critical
/// constriction.
/// \param formantFreq An array for the formant frequencies for at least 
/// MAX_FORMANTS entries
/// \param formantBW An array for the formant bandwidths for at least 
/// MAX_FORMANTS entries
/// \param numFormants The actual number of extracted formants
/// \param MAX_FORMANTS Maximal number of formants to be extracted
/// \param frictionNoise Shows the spectrum signs of a critical constriction?
/// \param isClosure Has the vocal tract a closure?
/// \param isNasal Are there clear anti-formants in the spectrum indicating 
/// nasality?
// ****************************************************************************

void TlModel::getFormants(double *formantFreq, double *formantBW, int &numFormants, 
  const int MAX_FORMANTS, bool &frictionNoise, bool &isClosure, bool &isNasal)
{
  // Extract the transfer function **********************************
  
  const int SPECTRUM_LENGTH_EXPONENT = 11;
  const int SPECTRUM_LENGTH = 1 << SPECTRUM_LENGTH_EXPONENT;
  ComplexSignal s;

  getSpectrum(FLOW_SOURCE_TF, &s, SPECTRUM_LENGTH, Tube::FIRST_PHARYNX_SECTION);

  // ****************************************************************
  // Determine the formants.
  // ****************************************************************

  numFormants = 0;

  int i, k;
  double f0 = (double)SAMPLING_RATE / (double)SPECTRUM_LENGTH;
  int firstSample = (int)(150.0 / f0);    // Minimum frequency = 150 Hz
  int lastSample = (int)(7000.0 / f0);    // Maximum frequency = 7000 Hz
  double a0, a1, a2;
  double formantAmp[32];    // Expect not more than 32 formant peaks
  double currHeight;
  double threshold;
  bool leftOK, rightOK;
  const double ABS_MIN_THRESHOLD = 0.316;  // Absolute threshold of -10 dB for all peaks
  const double EPSILON = 0.000001;
  double leftBwFreq = 0.0;
  double rightBwFreq = 0.0;
  double t;
  double den;

  for (i=firstSample; i < lastSample; i++)
  {
    // Is there a local maximum ? ***********************************
    a0 = s.getMagnitude(i-1);
    a1 = s.getMagnitude(i);
    a2 = s.getMagnitude(i+1);

    if ((a1 >= a0) && (a1 > a2) && (a1 >= ABS_MIN_THRESHOLD) && (numFormants < MAX_FORMANTS))
    {
      // Parabolic interpolation ************************************

      formantFreq[numFormants] = f0*(i + 0.5*(a2-a0) / (2.0*a1 - a0 - a2));
      formantAmp[numFormants]  = a1 + (a2-a0)*(a2-a0) / (8.0*(2.0*a1 - a0 - a2));
      formantBW[numFormants] = 0.0;

      // ************************************************************
      // Can we find samples to the left and the right of the current 
      // position, where the magnitude dropped 1 dB below the current
      // magnitude, without higher magnitudes inbetween ?
      // ************************************************************

      leftOK = false;
      rightOK = false;
      currHeight = s.getMagnitude(i);
      threshold = currHeight*0.891;     // Factor for -1 dB

      k = i-1;
      while ((k > firstSample) && (s.getMagnitude(k) <= currHeight) && (s.getMagnitude(k) > threshold)) { k--; }
      if (s.getMagnitude(k) <= threshold) { leftOK = true; }

      k = i+1;
      while ((k < lastSample) && (s.getMagnitude(k) <= currHeight) && (s.getMagnitude(k) > threshold)) { k++; }
      if (s.getMagnitude(k) <= threshold) { rightOK = true; }

      // If so, the regard the peak as a new formant ****************

      if ((leftOK) && (rightOK))
      {
        // Determine the bandwidth **********************************

        leftOK = false;
        rightOK = false;
        threshold = currHeight*0.708;     // Factor for -3 dB

        k = i-1;
        while ((k > firstSample) && (s.getMagnitude(k) <= currHeight) && (s.getMagnitude(k) > threshold)) { k--; }
        if (s.getMagnitude(k) <= threshold) 
        { 
          leftOK = true; 
          den = s.getMagnitude(k+1) - s.getMagnitude(k);
          if (den < EPSILON) { den = EPSILON; }
          t = (threshold - s.getMagnitude(k)) / den;      // 0 <= t <= 1
          leftBwFreq = ((double)k + t)*f0;
        }

        k = i+1;
        while ((k < lastSample) && (s.getMagnitude(k) <= currHeight) && (s.getMagnitude(k) > threshold)) { k++; }
        if (s.getMagnitude(k) <= threshold) 
        { 
          rightOK = true; 
          den = s.getMagnitude(k) - s.getMagnitude(k-1);    // den < 0
          if (den > -EPSILON) { den = -EPSILON; }
          t = (threshold - s.getMagnitude(k-1)) / den;      // 0 <= t <= 1
          rightBwFreq = ((double)(k-1) + t)*f0;
        }

        // Set the new bandwidth value ******************************

        if ((rightOK) && (leftOK))
        {
          formantBW[numFormants] = rightBwFreq - leftBwFreq;
        }
        else
        if (leftOK)
        {
          formantBW[numFormants] = 2.0*(formantFreq[numFormants] - leftBwFreq);
        }
        else
        if (rightOK)
        {
          formantBW[numFormants] = 2.0*(rightBwFreq - formantFreq[numFormants]);
        }
        else
        {
          formantBW[numFormants] = 100.0;     // Default value for the error case
        }

        numFormants++;
      }       // We had the necessary drop of -1 dB to each side

    }
  }

  // ****************************************************************
  // Is there a closure in the vocal tract tube ? Assume that, when
  // only one or no formant peaks are above 0 dB below 4 kHz!
  // ****************************************************************

  k = 0;
  for (i=0; (i < numFormants) && (i < 3); i++)
  {
    if ((formantAmp[i] >= 1.0) && (formantFreq[i] < 4000.0)) { k++; }
  }
  if (k < 2) { isClosure = true; } else { isClosure = false; }

  // ****************************************************************
  // Is friction noise being generated due to a critical constriction?
  // ****************************************************************

  double oldGlottalArea0 = tube.section[Tube::LOWER_GLOTTIS_SECTION]->area_cm2;
  double oldGlottalArea1 = tube.section[Tube::LOWER_GLOTTIS_SECTION + 1]->area_cm2;

  tube.section[Tube::LOWER_GLOTTIS_SECTION]->area_cm2 = 0.3;   // = 0.3 cm^2 (Glottal opening for fricatives)
  tube.section[Tube::LOWER_GLOTTIS_SECTION + 1]->area_cm2 = 0.3;

  double meanFlow = getMeanFlow(8000.0);    // Mean volume velocity with 8000 dPa subglottal pressure
  
  double A_c = tube.section[getMostConstrictedSection()]->area_cm2;
  if (A_c < MIN_AREA_CM2) 
  { 
    A_c = MIN_AREA_CM2; 
  }
  
  double diameter = 2.0*sqrt(A_c/M_PI);
  double velocity = meanFlow / A_c;
  
  double Re = velocity*diameter*AMBIENT_DENSITY_CGS/AIR_VISCOSITY_CGS;  // Reynolds-number
  const double Re_crit = 1800.0;      // A little bit higher than normal
  const double Re_threshold = 70000000.0;    //70000000.0; <- original

  if (Re*Re - Re_crit*Re_crit > Re_threshold) { frictionNoise = true; } else { frictionNoise = false; }

  tube.section[Tube::LOWER_GLOTTIS_SECTION]->area_cm2 = oldGlottalArea0;
  tube.section[Tube::UPPER_GLOTTIS_SECTION]->area_cm2 = oldGlottalArea1;

  // ****************************************************************
  // Is the nasal port open ?
  // ****************************************************************

  const double NASALITY_THRESHOLD_CM2 = 0.01;   // = 1 mm^2
  if (tube.section[Tube::FIRST_NOSE_SECTION]->area_cm2 > NASALITY_THRESHOLD_CM2) 
  { 
    isNasal = true; 
  } 
  else 
  { 
    isNasal = false; 
  }
}


// ****************************************************************************
/// Returns the circumference of a circle with the given area.
/// \param area The area of the circle.
// ****************************************************************************

double TlModel::getCircumference(double area)
{
  return 2.0*sqrt(area*M_PI);
}


// ****************************************************************************
/// Prepares the determination of spectral data by the pre-calculation of all
/// necessary matrices in the given frequency raster.
/// The frequency raster must have been set beforehand by f0 and numFreq.
/// For each tube section and all discrete frequency values, the matrix product
/// will be calculated from the first tube section in the branch to the 
/// current tube section.
// ****************************************************************************

void TlModel::prepareCalculations()
{
  int i, k, m;
  double omega;
  Matrix2x2 M, K;
  Tube::Section *ts = NULL;
  ComplexValue fossaInputImpedance;
  ComplexValue inputImpedance;

  // ****************************************************************
  // Keep in mind the current options and tube geometry
  // ****************************************************************

  prevOptions = options;
  prevTube = tube;
  
  // ****************************************************************
  // Pre-calculate some frequency data.
  // ****************************************************************

  resetCalculations = false;

  if (numFreq > MAX_NUM_FREQ) 
  { 
    numFreq = MAX_NUM_FREQ; 
  }

  for (i=0; i < numFreq; i++)
  {
    omega = 2.0*M_PI*f0*(double)i;
    
    discreteOmega[i] = omega;
    mouthRadiationImpedance[i]  = getRadiationImpedance(omega, tube.section[Tube::LAST_MOUTH_SECTION]->area_cm2);
    noseRadiationImpedance[i] = getRadiationImpedance(omega, tube.section[Tube::LAST_NOSE_SECTION]->area_cm2);
    lungTerminationImpedance[i] = 0.0;
    radiationCharacteristic[i]  = getRadiationCharacteristic(omega);
  }

  // ****************************************************************
  // Calculate the mean flow through the vocal tract.
  // ****************************************************************

  double A_g = tube.section[Tube::LOWER_GLOTTIS_SECTION]->area_cm2;    // glottal area
  int minAreaSection = getMostConstrictedSection();
  double A_c = tube.section[minAreaSection]->area_cm2;

  if (A_g < MIN_AREA_CM2) { A_g = MIN_AREA_CM2; }
  if (A_c < MIN_AREA_CM2) { A_c = MIN_AREA_CM2; }

  // The mean flow through the vocal tract
  double meanFlow = sqrt(2.0*lungPressure_dPa / (AMBIENT_DENSITY_CGS*(1.0/(A_g*A_g) + 1.0/(A_c*A_c))));

  // When the nasal tract is coupled -> meanFlow = 0.0
  if (tube.section[Tube::FIRST_NOSE_SECTION]->area_cm2 > MIN_AREA_CM2) 
  { 
    meanFlow = 0.0; 
  }

  // ****************************************************************
  // Calculate the product matrices for the tube sections.
  // ****************************************************************

  for (i=0; i < numFreq; i++)
  {
    omega = discreteOmega[i];

    // **************************************************************
    // Input impedance of the piriform fossa.
    // **************************************************************

    K.unitMatrix();
    for (k = Tube::FIRST_FOSSA_SECTION; k <= Tube::LAST_FOSSA_SECTION; k++)
    {
      ts = tube.section[k];
      K*= getSectionMatrix(omega, k);
      matrixProduct[k][i] = K;
    }
    fossaInputImpedance = K.A / K.C;

    // **************************************************************
    // The matrix products of the subglottal system and pharynx.
    // **************************************************************

    K.unitMatrix();

    for (k=Tube::FIRST_TRACHEA_SECTION; k <= Tube::LAST_PHARYNX_SECTION; k++)   
    {
      ts = tube.section[k];

      // Add an "inner length correction" (additional inductivity)
      // between the previous and the current tube section as
      // described in Sondhi (1983).
      if ((k > Tube::FIRST_PHARYNX_SECTION) && (k <= Tube::LAST_MOUTH_SECTION) && (options.innerLengthCorrections))
      {
        M.unitMatrix();
        M.B = getJunctionImpedance(omega, tube.section[k-1]->area_cm2, tube.section[k]->area_cm2);
        K*= M;
      }

      K*= getSectionMatrix(omega, k);
      matrixProduct[k][i] = K;

      // Add the differential small-signal resistance at the glottis

      if ((k == Tube::LAST_TRACHEA_SECTION) && (options.staticPressureDrops))
      {
        M.unitMatrix();
        M.B = AMBIENT_DENSITY_CGS*meanFlow / (A_g*A_g);
        K*= M;
      }

      // Put a small-signal flow resistance at the entrance of the supraglottal constriction

      if ((k == minAreaSection-1) && (options.staticPressureDrops))
      {
        M.unitMatrix();
        M.B = AMBIENT_DENSITY_CGS*meanFlow / (A_c*A_c);
        K*= M;
      }

      // Consider the piriform fossa as a side branch.
      
      if ((k == Tube::FIRST_PHARYNX_SECTION + Tube::FOSSA_COUPLING_SECTION) && (options.piriformFossa))
      {
        M.unitMatrix();
        M.C = 1.0 / fossaInputImpedance;
        K*= M;
      }
    }       // Loop for the sections of the trachea + glottis + pharynx

    // **************************************************************
    // The matrix products of the mouth cavity.
    // **************************************************************

    K.unitMatrix();
    for (k = Tube::FIRST_MOUTH_SECTION; k <= Tube::LAST_MOUTH_SECTION; k++)
    {
      ts = tube.section[k];

      // Add an "inner length correction" (additional inductivity)
      // between the previous and the current tube section as
      // described in Sondhi (1983).
      if ((k > Tube::FIRST_PHARYNX_SECTION) && (k <= Tube::LAST_MOUTH_SECTION) && (options.innerLengthCorrections))
      {
        M.unitMatrix();
        M.B = getJunctionImpedance(omega, tube.section[k-1]->area_cm2, tube.section[k]->area_cm2);
        K*= M;
      }

      K*= getSectionMatrix(omega, k);
      matrixProduct[k][i] = K;

      // Put a small-signal flow resistance at the entrance of the supraglottal constriction

      if ((k == minAreaSection-1) && (options.staticPressureDrops))
      {
        M.unitMatrix();
        M.B = AMBIENT_DENSITY_CGS*meanFlow / (A_c*A_c);
        K*= M;
      }
    }     // Loop for the mouth sections

    // **************************************************************
    // The matrix products of the nasal cavity.
    // **************************************************************

    K.unitMatrix();
    for (k = Tube::FIRST_NOSE_SECTION; k <= Tube::LAST_NOSE_SECTION; k++)   
    {
      ts = tube.section[k];
      K*= getSectionMatrix(omega, k);
      matrixProduct[k][i] = K;

      // Coupling of the paranasal sinuses ? ************************

      if (options.paranasalSinuses)
      {
        for (m=0; m < Tube::NUM_SINUS_SECTIONS; m++)
        {
          if (k == Tube::FIRST_NOSE_SECTION + Tube::SINUS_COUPLING_SECTION[m])
          {
            M = getSectionMatrix(omega, Tube::FIRST_SINUS_SECTION + m);
            inputImpedance = M.A / M.C;
            M.unitMatrix();
            M.C = 1.0 / inputImpedance;
            K*= M;
          }
        }
      }
    }     // Loop for the nose sections

  }   // Loop for the frequencies
}


// ****************************************************************************
/// Returns the spectral value of the radiation characteristic at the angular
/// frequency omega.
/// \param omega Angular frequency
// ****************************************************************************

ComplexValue TlModel::getRadiationCharacteristic(double omega)
{
  const double MICROPHONE_DISTANCE_CM = 25.0;    // 25 cm
  ComplexValue H;
  
  H = ComplexValue(0.0, (AMBIENT_DENSITY_CGS*omega) / (4.0*M_PI*MICROPHONE_DISTANCE_CM));
  
  double angle = (omega*MICROPHONE_DISTANCE_CM) / SOUND_VELOCITY_CGS;
  double re = cos(angle);
  double im = -sin(angle);
  H*= ComplexValue(re, im);
  return H;
}


// ****************************************************************************
/// Returns the value of the radiation impedance for the given angular
/// frequency and radiation area.
/// \param omega Angular frequency
/// \param radiationArea_cm2 Radiation area in cm^2
// ****************************************************************************

ComplexValue TlModel::getRadiationImpedance(double omega, double radiationArea_cm2)
{
  // Init. the radiation impedance with 0 (no radiation)
  ComplexValue radiationImpedance(0.0, 0.0);

  if (omega < MIN_FREQ_RAD) 
  { 
    omega = MIN_FREQ_RAD; 
  }
  if (radiationArea_cm2 < MIN_AREA_CM2) 
  { 
    radiationArea_cm2 = MIN_AREA_CM2; 
  }

  // ****************************************************************
  // Approximation by a piston in an infinit wall.
  // ****************************************************************

  if (options.radiation == PISTONINWALL_RADIATION)
  {
    double a = sqrt(radiationArea_cm2 / M_PI);   // Radius of the mouth opening in cm
    double k = omega / SOUND_VELOCITY_CGS;      // Wave number
    double ka = k*a;

    double analogImpedance = (AMBIENT_DENSITY_CGS * SOUND_VELOCITY_CGS) / radiationArea_cm2;
    radiationImpedance = analogImpedance * ComplexValue(0.5*ka*ka, (8.0*ka) / (3.0*M_PI));
  }
  else

    // ****************************************************************
    // Approximation by a piston in a sphere (Wakita, Fant).
    // ****************************************************************

  if (options.radiation == PISTONINSPHERE_RADIATION)
  {
    double K;
    double freq = omega / (2.0*M_PI);
    if (freq < 1600) 
    { 
      K = (0.6*freq) / 1600.0 + 1.0; 
    } 
    else 
    { 
      K = 1.6; 
    }
    double R0 = (AMBIENT_DENSITY_CGS*omega*omega*K) / (4.0*M_PI*SOUND_VELOCITY_CGS);
    double L0 = (8.0*AMBIENT_DENSITY_CGS) / (3.0*M_PI*sqrt(M_PI*radiationArea_cm2));

    radiationImpedance = ComplexValue(R0, omega*L0);
  }

  else

  // Approximation by a parallel impedance of a resistor and a coil
  
  if (options.radiation == PARALLEL_RADIATION)
  {
    ComplexValue Z1 = (128.0*AMBIENT_DENSITY_CGS*SOUND_VELOCITY_CGS) / (9.0*M_PI*M_PI*radiationArea_cm2);
    ComplexValue Z2 = ComplexValue(0.0, (omega*8.0*AMBIENT_DENSITY_CGS) / (3.0*M_PI*sqrt(M_PI*radiationArea_cm2)));
 
    radiationImpedance = (Z1*Z2) / (Z1 + Z2);
  }

  return radiationImpedance;
}


// ****************************************************************************
/// Returns the series and parallel impedance of a two-port network of a tube 
/// section with lumped elements.
/// \param omega Angular frequency
/// \param ts Pointer to the tube section
/// \param Za Return value for the series impedance
/// \param Zb Return value for the parallel impedance
// ****************************************************************************

void TlModel::getLumpedSectionImpedances(double omega, Tube::Section *ts, ComplexValue &Za, ComplexValue &Zb)
{
  if (ts == NULL)
  {
    Za = 0.0;
    Zb = 0.0;
    return;
  }
  if (omega < MIN_FREQ_RAD) 
  { 
    omega = MIN_FREQ_RAD; 
  }

  double area = ts->area_cm2;
  if (area < MIN_AREA_CM2)  
  { 
    area = MIN_AREA_CM2; 
  }
  
  double circ = getCircumference(area);
  double length = ts->length_cm;

  double Rs = RS_FACTOR*(0.5*length*circ*sqrt(omega)) / (area*area);
  double Ls = LS_FACTOR*(0.5*length) / (area);
  double Cp = CP_FACTOR*length*area;
  double Gp = GP_FACTOR*length*circ*sqrt(omega);
  double Rf = (4.0*AIR_VISCOSITY_CGS*length*M_PI) / (area*area);  // flow resistance for a circular tube

  // Series impedance.

  Za = ComplexValue(0, omega*Ls);
  if (options.boundaryLayer)   { Za+= Rs; }
  if (options.hagenResistance) { Za+= Rf; }

  // Parallel impedance.

  Zb = ComplexValue(0, omega*Cp);
  if (options.heatConduction) { Zb+= Gp; }
  if (options.softWalls)
  { 
    Zb+= 1.0 / (ComplexValue(ts->wallResistance_cgs, omega*ts->wallMass_cgs - ts->wallStiffness_cgs/omega) / (length*circ));
  }
  Zb = 1.0 / Zb;        // From conductivity to impedance
}


// ****************************************************************************
/// Returns the 2x2 matrix of a tube section in the form:
/// (P_in, U_in) = (A B C D)*(P_out, U_out).
/// \param omega Angular frequency
/// \param section Section index
// ****************************************************************************

Matrix2x2 TlModel::getSectionMatrix(double omega, int section)
{
  Tube::Section *ts = tube.section[section];

  Matrix2x2 M;

  if (omega < MIN_FREQ_RAD) 
  { 
    omega = MIN_FREQ_RAD; 
  }

  double area = ts->area_cm2;
  if (area < MIN_AREA_CM2)  
  { 
    area = MIN_AREA_CM2; 
  }

  double circ = getCircumference(area);
  double length = ts->length_cm;

  // ****************************************************************
  // Special treatment of Helmholtz resonators.
  // ****************************************************************

  if ((section >= Tube::FIRST_SINUS_SECTION) && (section <= Tube::LAST_SINUS_SECTION))
  {
    double innerSurface = 4.0*M_PI*pow((3.0*ts->volume_cm3) / (4.0*M_PI), 2.0/3.0);

    double Rs = RS_FACTOR*length*circ*sqrt(omega) / (area*area);
    double Ls = LS_FACTOR*length / area;
    double Cp = CP_FACTOR*ts->volume_cm3;
    double Gp = innerSurface / 65000.0;

    ComplexValue Za = ComplexValue(0, omega*Ls);
    if (options.boundaryLayer) { Za+= Rs; }

    ComplexValue Zb = 1.0 / ComplexValue(Gp, omega*Cp);   // 1.0 / Sum of the conductivities

    M.D = -1.0;
    M.B = -Zb;
    M.C = 1.0/Zb;
    M.A = 1.0 + Za/Zb;
  }
  else

  // ****************************************************************
  // Approximation with lumped elements (quasi-static sound fields)
  // ****************************************************************

  if (options.lumpedElements)
  {
    ComplexValue Za, Zb;
    getLumpedSectionImpedances(omega, ts, Za, Zb);

    // The matrix elements
    M.A = M.D = 1.0 + (Za/Zb);
    M.B = (Za*Za)/Zb + 2.0*Za;
    M.C = 1.0 / Zb;
  }
  else

  // ****************************************************************
  // Real wave propagation.
  // ****************************************************************
  {
    double Rs = RS_FACTOR*circ*sqrt(omega) / (area*area);
    double Ls = LS_FACTOR / area;
    double Cp = CP_FACTOR*area;
    double Gp = GP_FACTOR*circ*sqrt(omega);
    double Rf = 8.0*AIR_VISCOSITY_CGS*M_PI / (area*area);

    // The entire impedance in series (per length-unit)
    ComplexValue z(0, omega*Ls);
    if (options.boundaryLayer)   { z+= Rs; }
    if (options.hagenResistance) { z+= Rf; }

    // The entire conductivity in parallel (per length-unit)
    ComplexValue y(0, omega*Cp);
    if (options.heatConduction) { y+= Gp; }
    if (options.softWalls)
    { 
      y += 1.0 / (ComplexValue(ts->wallResistance_cgs, omega*ts->wallMass_cgs - ts->wallStiffness_cgs / omega) / circ);
    }

    ComplexValue Z0 = std::sqrt(z/y);
    ComplexValue gamma = std::sqrt(z*y);

    // The matrix elements.

    M.A = M.D = std::cosh(gamma*length);
    ComplexValue SINH = std::sinh(gamma*length);
    M.B = SINH * Z0;
    M.C = SINH / Z0;
  }

  return M;
}


// ****************************************************************************
/// Calculates the junction impedance between two adjacent tube sections with
/// the given areas according to SONDHI (1983).
/// This additional inductivity is a kind of "inner length correction" that
/// affects the formant frequencies when there are abrupt changes in the area
/// function.
// ****************************************************************************

ComplexValue TlModel::getJunctionImpedance(double omega, double A1_cm2, double A2_cm2)
{
  ComplexValue Z;
  double a, b;    // Radii of the bigger and smaller tube section

  // Safety checks.

  if (A1_cm2 < MIN_AREA_CM2)
  {
    A1_cm2 = MIN_AREA_CM2;
  }
  
  if (A2_cm2 < MIN_AREA_CM2)
  {
    A2_cm2 = MIN_AREA_CM2;
  }

  // ****************************************************************

  if (A1_cm2 > A2_cm2)
  {
    a = sqrt(A1_cm2 / M_PI);
    b = sqrt(A2_cm2 / M_PI);
  }
  else
  {
    a = sqrt(A2_cm2 / M_PI);
    b = sqrt(A1_cm2 / M_PI);
  }

  double H = 1 - b/a;
  Z = ComplexValue(0, omega*8.0*AMBIENT_DENSITY_CGS*H / (3.0*M_PI*M_PI*b));

  return Z;
}


// ****************************************************************************
/// Calculates the input impedance looking into the given tube section to the
/// right.
/// \param freqIndex Index of the spectral line
/// \param section Index of the tube section
// ****************************************************************************

ComplexValue TlModel::getInputImpedance(int freqIndex, int section)
{
  Matrix2x2 M;
  ComplexValue loadImpedance;

  if (freqIndex == 0) { freqIndex = 1; }

  // ****************************************************************
  // Is the tube section in the branch lung/pharynx ?
  // ****************************************************************

  if ((section >= Tube::FIRST_TRACHEA_SECTION) && (section <= Tube::LAST_PHARYNX_SECTION))
  {
    ComplexValue noseInputImpedance  = getInputImpedance(freqIndex, Tube::FIRST_NOSE_SECTION);
    ComplexValue mouthInputImpedance = getInputImpedance(freqIndex, Tube::FIRST_MOUTH_SECTION);
    loadImpedance = (noseInputImpedance*mouthInputImpedance) / (noseInputImpedance+mouthInputImpedance);

    M.unitMatrix();
    if (section > Tube::FIRST_TRACHEA_SECTION)
    {
      M = matrixProduct[section-1][freqIndex];
      M.invert();
    }
    M*= matrixProduct[Tube::LAST_PHARYNX_SECTION][freqIndex];
  }
  else

  // ****************************************************************
  // Is the tube section in the mouth ?
  // ****************************************************************

  if ((section >= Tube::FIRST_MOUTH_SECTION) && (section <= Tube::LAST_MOUTH_SECTION))
  {
    loadImpedance = mouthRadiationImpedance[freqIndex];

    M.unitMatrix();
    if (section > Tube::FIRST_MOUTH_SECTION)
    {
      M = matrixProduct[section-1][freqIndex];
      M.invert();
    }
    M*= matrixProduct[Tube::LAST_MOUTH_SECTION][freqIndex];
  }
  else

  // ****************************************************************
  // Is the tube section in the nasal cavity ?
  // ****************************************************************

  if ((section >= Tube::FIRST_NOSE_SECTION) && (section <= Tube::LAST_NOSE_SECTION))
  {
    loadImpedance = noseRadiationImpedance[freqIndex];

    M.unitMatrix();
    if (section > Tube::FIRST_NOSE_SECTION)
    {
      M = matrixProduct[section-1][freqIndex];
      M.invert();
    }
    M*= matrixProduct[Tube::LAST_NOSE_SECTION][freqIndex];
  }

  return (M.A*loadImpedance + M.B) / (M.C*loadImpedance + M.D);
}


// ****************************************************************************
/// Calculates the output impedance of a tube section, when you look into it
/// to the left.
/// \param freqIndex Index of the spectral line
/// \param section Index of the tube section
// ****************************************************************************

ComplexValue TlModel::getOutputImpedance(int freqIndex, int section)
{
  Matrix2x2 M, K;
  ComplexValue lungImpedance = lungTerminationImpedance[freqIndex];

  if (freqIndex == 0) { freqIndex = 1; }

  // ****************************************************************
  // Is the tube section in the branch lung/pharynx ?
  // ****************************************************************

  if ((section >= Tube::FIRST_TRACHEA_SECTION) && (section <= Tube::LAST_PHARYNX_SECTION))
  {
    M = matrixProduct[section][freqIndex];
  }
  else

  // ****************************************************************
  // Is the tube section in the mouth ?
  // ****************************************************************

  if ((section >= Tube::FIRST_MOUTH_SECTION) && (section <= Tube::LAST_MOUTH_SECTION))
  {
    M = matrixProduct[Tube::LAST_PHARYNX_SECTION][freqIndex];
    K.unitMatrix();
    K.C = 1.0 / getInputImpedance(freqIndex, Tube::FIRST_NOSE_SECTION);
    M*= K;          // Coupling matrix for the nasal cavity
    M*= matrixProduct[section][freqIndex];
  }
  else

  // ****************************************************************
  // Is the tube section in the nasal cavity ?
  // ****************************************************************

  if ((section >= Tube::FIRST_NOSE_SECTION) && (section <= Tube::LAST_NOSE_SECTION))
  {
    M = matrixProduct[Tube::LAST_PHARYNX_SECTION][freqIndex];
    K.unitMatrix();
    K.C = 1.0 / getInputImpedance(freqIndex, Tube::FIRST_MOUTH_SECTION);
    M*= K;        // Coupling matrix for the mouth cavity
    M*= matrixProduct[section][freqIndex];
  }

  return (lungImpedance*M.D + M.B) / (M.A + lungImpedance*M.C);
}


// ****************************************************************************
/// Returns the value of the transfer function from a sound pressure source 
/// right before the given section to the flow at the mouth-/nose-opening 
/// (U_out/P_in).
/// When section == -1, the pressure source position is at the right end of the
/// last oral tube section.
/// \param freqIndex Index of the spectral line
/// \param section Index of the tube section
// ****************************************************************************

ComplexValue TlModel::getPressureSourceTF(int freqIndex, int section)
{
  Matrix2x2 M, K;
  ComplexValue result;
  ComplexValue factor;
  ComplexValue inputImpedance;
  ComplexValue outputImpedance;
  ComplexValue Za, Zb;
  
  if (freqIndex == 0) { freqIndex = 1; }

  // ****************************************************************
  // Is the source in the branch lung/pharynx ?
  // ****************************************************************

  if ((section >= Tube::FIRST_TRACHEA_SECTION) && (section <= Tube::LAST_PHARYNX_SECTION))
  {
    result = 0.0;
    Matrix2x2 pharynxMatrix;

    // The input impedance looking to the right.

    inputImpedance = getInputImpedance(freqIndex, section);

    // The output impedance looking to the left.

    if (section == Tube::FIRST_TRACHEA_SECTION)
    {
      outputImpedance = lungTerminationImpedance[freqIndex];
    }
    else
    {
      outputImpedance = getOutputImpedance(freqIndex, section-1);
    }

    factor = 1.0 / (outputImpedance + inputImpedance);

    // The matrix from section -> Tube::LAST_PHARYNX_SECTION

    M.unitMatrix();
    if (section > Tube::FIRST_TRACHEA_SECTION) 
    { 
      M = matrixProduct[section-1][freqIndex]; 
    }
    M.invert();
    M*= matrixProduct[Tube::LAST_PHARYNX_SECTION][freqIndex];
    pharynxMatrix = M;

    // The TF from section -> Tube::LAST_MOUTH_SECTION
    
    K = pharynxMatrix;

    M.unitMatrix();
    M.C = 1.0 / getInputImpedance(freqIndex, Tube::FIRST_NOSE_SECTION);
    K*= M;      // Kopplungsmatrix zum Nasenraum
    K*= matrixProduct[Tube::LAST_MOUTH_SECTION][freqIndex];
    result+= factor / (K.C*mouthRadiationImpedance[freqIndex] + K.D);

    // The TF from section -> Tube::LAST_NOSE_SECTION
    
    K = pharynxMatrix;

    M.unitMatrix();
    M.C = 1.0 / getInputImpedance(freqIndex, Tube::FIRST_MOUTH_SECTION);
    K*= M;      // Kopplungsmatrix zum Mundraum
    K*= matrixProduct[Tube::LAST_NOSE_SECTION][freqIndex];
    result+= factor / (K.C*noseRadiationImpedance[freqIndex] + K.D);
  }
  else

  // ****************************************************************
  // Is the source in the mouth cavity ?
  // ****************************************************************

  if (((section >= Tube::FIRST_MOUTH_SECTION) && (section <= Tube::LAST_MOUTH_SECTION)) || 
      (section == -1))
  {
    result = 0.0;

    // The input impedance looking to the right.

    if (section == -1)
    {
      inputImpedance = mouthRadiationImpedance[freqIndex];
    }
    else
    {
      inputImpedance = getInputImpedance(freqIndex, section);
    }

    // The output impedance looking to the left.

    if (section == -1)
    {
      outputImpedance = getOutputImpedance(freqIndex, Tube::LAST_MOUTH_SECTION);
    }
    else
    if (section == Tube::FIRST_MOUTH_SECTION)
    {
      Za = getOutputImpedance(freqIndex, Tube::LAST_PHARYNX_SECTION);
      Zb = getInputImpedance(freqIndex, Tube::FIRST_NOSE_SECTION);
      outputImpedance = (Za*Zb) / (Za+Zb);
    }
    else
    {
      outputImpedance = getOutputImpedance(freqIndex, section-1);
    }

    factor = 1.0 / (outputImpedance + inputImpedance);

    // The matrix from section -> Tube::LAST_MOUTH_SECTION

    if (section == -1)
    {
      M.unitMatrix();
    }
    else
    if (section > Tube::FIRST_MOUTH_SECTION)
    {
      M = matrixProduct[section-1][freqIndex];
      M.invert();
      M*= matrixProduct[Tube::LAST_MOUTH_SECTION][freqIndex];
    }
    else
    {
      M = matrixProduct[Tube::LAST_MOUTH_SECTION][freqIndex];
    }

    result+= factor / (M.C*mouthRadiationImpedance[freqIndex] + M.D);

    // The matrix from section -> Tube::LAST_NOSE_SECTION

    M.unitMatrix();
    if (section == -1)
    {
      M = matrixProduct[Tube::LAST_MOUTH_SECTION][freqIndex];
    }
    else
    if (section > Tube::FIRST_MOUTH_SECTION)
    {
      M = matrixProduct[section-1][freqIndex];
    }
    M.invert();

    K.unitMatrix();
    K.C = 1.0 / getOutputImpedance(freqIndex, Tube::LAST_PHARYNX_SECTION);
    M*= K;        // Coupling matrix for the sub-velar system

    M*= matrixProduct[Tube::LAST_NOSE_SECTION][freqIndex];

    result+= factor / (M.C*noseRadiationImpedance[freqIndex] + M.D);
  }

  return result;
}


// ****************************************************************************
/// Calculates the volume velocity transfer function from a flow source in the
/// center of the given tube section to the flow at the mouth-/nose opening
/// (U_out/U_in).
/// \param freqIndex Index of the spectral line
/// \param section Index of the tube section
// ****************************************************************************

ComplexValue TlModel::getFlowSourceTF(int freqIndex, int section)
{
  Matrix2x2 M, K;
  ComplexValue result;
  ComplexValue factor;
  ComplexValue inputImpedance;
  ComplexValue outputImpedance;
  ComplexValue sourceImpedance;
  ComplexValue Za, Zb;

  if (freqIndex == 0) { freqIndex = 1; }

  // ****************************************************************
  // If the section is the first section of the nasal cavity, treat 
  // the nasal cavity as if it was closed at the posterior end.
  // ****************************************************************

  if (section == Tube::FIRST_NOSE_SECTION)
  {
    K = matrixProduct[Tube::LAST_NOSE_SECTION][freqIndex];
    result = 1.0 / (K.C*noseRadiationImpedance[freqIndex] + K.D);
    // Nothing more to do -> Return right here!
    return result;
  }
  

  // ****************************************************************
  // Is the source in the branch lung/pharynx ?
  // ****************************************************************

  if ((section >= Tube::FIRST_TRACHEA_SECTION) && (section <= Tube::LAST_PHARYNX_SECTION))
  {
    result = 0.0;

    // The input impedance looking to the right

    if (section == Tube::LAST_PHARYNX_SECTION)
    {
      Za = getInputImpedance(freqIndex, Tube::FIRST_MOUTH_SECTION);
      Zb = getInputImpedance(freqIndex, Tube::FIRST_NOSE_SECTION);
      inputImpedance = (Za*Zb) / (Za+Zb);
    }
    else
    { 
      inputImpedance = getInputImpedance(freqIndex, section+1);
    }

    // The output impedance looking to the left

    if (section == Tube::FIRST_TRACHEA_SECTION)
    {
      outputImpedance = lungTerminationImpedance[freqIndex];
    }
    else
    {
      outputImpedance = getOutputImpedance(freqIndex, section-1);
    }

    // The parallel source impedance of the tube section

    getLumpedSectionImpedances(discreteOmega[freqIndex], tube.section[section], Za, Zb);
    inputImpedance+= Za;
    outputImpedance+= Za;
    sourceImpedance = Zb;

    // The factor for the TF

    factor = (sourceImpedance*outputImpedance) / (sourceImpedance + outputImpedance);
    factor = factor / (inputImpedance + factor);

    // **************************************************************
    // The matrix from section -> Tube::LAST_PHARYNX_SECTION
    // **************************************************************

    M = matrixProduct[section][freqIndex];
    M.invert();
    M*= matrixProduct[Tube::LAST_PHARYNX_SECTION][freqIndex];

    Matrix2x2 pharynxMatrix;
    pharynxMatrix.unitMatrix();
    pharynxMatrix.B = Za;
    pharynxMatrix*= M;

    // **************************************************************
    // The TF from section -> Tube::LAST_MOUTH_SECTION
    // **************************************************************

    K = pharynxMatrix;

    M.unitMatrix();
    M.C = 1.0 / getInputImpedance(freqIndex, Tube::FIRST_NOSE_SECTION);
    K*= M;      // Coupling matrix to the nose cavity

    K*= matrixProduct[Tube::LAST_MOUTH_SECTION][freqIndex];
    result+= factor / (K.C*mouthRadiationImpedance[freqIndex] + K.D);

    // **************************************************************
    // The TF from section -> Tube::LAST_NOSE_SECTION
    // **************************************************************

    K = pharynxMatrix;

    M.unitMatrix();
    M.C = 1.0 / getInputImpedance(freqIndex, Tube::FIRST_MOUTH_SECTION);
    K*= M;      // Coupling matrix to the mouth cavity

    K*= matrixProduct[Tube::LAST_NOSE_SECTION][freqIndex];
    result+= factor / (K.C*noseRadiationImpedance[freqIndex] + K.D);
  }
  else

  // ****************************************************************
  // Is the source in the mouth cavity ?
  // ****************************************************************

  if ((section >= Tube::FIRST_MOUTH_SECTION) && (section <= Tube::LAST_MOUTH_SECTION))
  {
    result = 0.0;

    // The input impedance looking to the right

    if (section == Tube::LAST_MOUTH_SECTION)
    {
      inputImpedance = mouthRadiationImpedance[freqIndex];
    }
    else
    { 
      inputImpedance = getInputImpedance(freqIndex, section+1);
    }

    // The output impedance looking to the left

    if (section == Tube::FIRST_MOUTH_SECTION)
    {
      Za = getOutputImpedance(freqIndex, Tube::LAST_PHARYNX_SECTION);
      Zb = getInputImpedance(freqIndex, Tube::FIRST_NOSE_SECTION);
      outputImpedance = (Za*Zb) / (Za+Zb);
    }
    else
    {
      outputImpedance = getOutputImpedance(freqIndex, section-1);
    }

    // The parallel source impedance of the tube section

    getLumpedSectionImpedances(discreteOmega[freqIndex], tube.section[section], Za, Zb);
    inputImpedance+= Za;
    outputImpedance+= Za;
    sourceImpedance = Zb;

    // The factor for the TF

    factor = (sourceImpedance*outputImpedance) / (sourceImpedance + outputImpedance);
    factor = factor / (inputImpedance + factor);

    // The matrix from section -> Tube::LAST_MOUTH_SECTION

    M = matrixProduct[section][freqIndex];
    M.invert();
    M*= matrixProduct[Tube::LAST_MOUTH_SECTION][freqIndex];

    K.unitMatrix();
    K.B = Za;
    K*= M;

    result+= factor / (K.C*mouthRadiationImpedance[freqIndex] + K.D);

    // The matrix from section -> Tube::LAST_NOSE_SECTION

    M = matrixProduct[section][freqIndex];
    M.invert();

    K.unitMatrix();
    K.C = 1.0 / getOutputImpedance(freqIndex, Tube::LAST_PHARYNX_SECTION);
    M*= K;        // Coupling matrix for the sub-velar system

    M*= matrixProduct[Tube::LAST_NOSE_SECTION][freqIndex];

    K.unitMatrix();
    K.B = Za;
    K*= M;

    result+= factor / (K.C*noseRadiationImpedance[freqIndex] + K.D);
  }

  return result;
}


// ****************************************************************************

