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

#ifndef __IIRFILTER_H__
#define __IIRFILTER_H__

#include <complex>
#include "Signal.h"

const int MAX_IIR_ORDER = 32;
const int IIR_BUFFER_MASK   = 63; 
const int IIR_BUFFER_LENGTH = 64;


typedef std::complex<double> ComplexValue;

// ****************************************************************************
/// This class represents a recursive infinit impulse response filter.
/// The transfer function has the form:
///
/// H(z) = [a0 + a1*z^(-1) + a2*z^(-2) + ...] / [1 - b1*z^(-1) - b2*z^(-2) - ...].
///
/// The recursion formula is thus:
///
/// y[n] = a0*x[n] + a1*x[n-1] + a2*x[n-2] + ... + b1*y[n-1] + b2*y[n-2] + ...
///
// ****************************************************************************

class IirFilter
{
  public:
    IirFilter();
    void resetBuffers(double initialValue = 0.0);

    double getOutputSample(double nextInputSample);
    ComplexValue getFrequencyResponse(double freqRatio);
    void getFrequencyResponse(ComplexSignal *spectrum, int spectrumLength);
    void getFrequencyResponse(ComplexSignal *spectrum, int spectrumLength, int SR, double F0);

    void setGain(double gain);
    void setCoefficients(const double *A, const double *B, const int newOrder);
    bool combineWithFilter(const IirFilter *f, bool cascade);

    // Creation of some simple filters.

    void createUnityFilter();
    void createSinglePoleLowpass(double cutoffFreqRatio);
    void createSinglePoleHighpass(double cutoffFreqRatio);
    void createSecondOrderLowpass(double freqRatio, double Q);
    void createChebyshev(double cutoffFreqRatio, bool isHighpass, int numPoles);

public:
    double a[MAX_IIR_ORDER+1];
    double b[MAX_IIR_ORDER+1];
    int order;

private:
    int pos;
    double inputBuffer[IIR_BUFFER_LENGTH];
    double outputBuffer[IIR_BUFFER_LENGTH];

    void clearCoefficients();
};


#endif
