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

#ifndef __DSP_H__
#define __DSP_H__

#include "Signal.h"
#include <complex>

typedef std::complex<double> ComplexValue;

enum WindowType { RECTANGULAR_WINDOW, HAMMING_WINDOW, RIGHT_HALF_OF_HAMMING_WINDOW, 
                  LEFT_HALF_OF_HAMMING_WINDOW, RIGHT_HALF_OF_HANN_WINDOW, GAUSS_WINDOW };

// ****************************************************************************

int modulo(int x, int y);
double getSignalEnergy(const Signal& signal, int startPos, int numSamples);
double getSignalEnergy(const Signal16& signal, int startPos, int numSamples);
double getMeanSignalPower(Signal& signal, int startPos, int numSamples);

// Fourier Transform **********************************************************

void rectangularToPolar(ComplexSignal& s, int length);
void polarToRectangular(ComplexSignal& s, int length);
void generateNegativeFrequencies(ComplexSignal *spectrum);

void realDFT(Signal& timeSignal, ComplexSignal& freqSignal, int length, bool normalize);
void realIDFT(ComplexSignal& freqSignal, Signal& timeSignal, int length, bool normalize);

void realIFFT(ComplexSignal& s, int lengthExponent, bool normalize);
void realFFT(ComplexSignal& s, int lengthExponent, bool normalize);

void complexIFFT(ComplexSignal& s, int lengthExponent, bool normalize);
void complexFFT(ComplexSignal& s, int lengthExponent, bool normalize);

void complexIDFT(ComplexSignal& freqSignal, ComplexSignal& timeSignal, int length, bool normalize);
void complexDFT(ComplexSignal& timeSignal, ComplexSignal& freqSignal, int length, bool normalize);

int getFrameLengthExponent(int windowLength_pt);

void getWindow(Signal& window, int length, WindowType type);

// Linear Predictive Coding ***************************************************

void getLPCCoefficients(const double *signal, int numSamples, double *coeff, int N);
void getLPCResidual(const double *signal, double *residual, long l, const double *coeff, long N);
void predictSignal(double *signal, const double *residual, long l, const double *coeff, long N);
void LPCToPolynomCoefficients(double *LPCCoeff, double *polynomCoeff, long N);

// Calculation of the zeros of a polynom. *************************************

void getSquareRoots(double beta, double gamma, ComplexValue &x0, ComplexValue &x1);
ComplexValue getPolynomValue(double *a, long N, ComplexValue x);
void getPolynomRoots(double *a, int &N, ComplexValue *roots);
void getRealPolynomRoots(double *a, int &N, double *roots, int& numRealRoots);

// Transformation between different frequency scales. *************************

inline double HzToBark(double f) { return (26.81*f) / (1960.0+f) - 0.53; }
inline double BarkToHz(double z) { return (1960.0*(z+0.53)) / (26.28-z); }
inline double HzToMel(double f) { return HzToBark(f)*100; }
inline double MelToHz(double m) { return BarkToHz(m/100.0); }
double hertzToSemitones(double freq_Hz);
double semitonesToHertz(double freq_st);


#endif
