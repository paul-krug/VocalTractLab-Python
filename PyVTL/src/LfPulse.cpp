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

#include "LfPulse.h"
#include <random>


// ****************************************************************************
/// Constructor. Sets the default parameter values.
// ****************************************************************************

LfPulse::LfPulse()
{
  resetParams();
}


// ****************************************************************************
// ****************************************************************************

void LfPulse::resetParams()
{
  AMP = 300.0;     // cm^3/s
  F0  = 120.0;     // Hz
  OQ  = 0.5;
  SQ  = 3.0;
  TL  = 0.02;
  SNR = 50.0;
}


// ****************************************************************************
// ****************************************************************************

void LfPulse::getPulse(Signal& s, int numSamples, bool getDerivative)
{
  const double MIN_TA = 0.01;   // Don't make is smaller!
  s.reset(numSamples);

  // ****************************************************************
  // Transform the parameters OQ, SQ and TL into the times te, tp and
  // ta.
  // ****************************************************************

  double T0 = 1.0;        // Scaling comes later
  double te = OQ;
  double tp = (te*SQ) / (1.0 + SQ);
  double ta = TL;
  if (ta < MIN_TA) { ta = MIN_TA; }
  if (ta > T0 - te) { ta = T0 - te; }

  double epsilon = getEpsilon(ta, te);
  double alpha   = getAlpha(tp, te, ta, epsilon);
  double B       = getB(AMP, tp, alpha);

  double w = 3.1415926 / tp;
  double t;

  int i;

  // ****************************************************************
  // Calculate the first derivative of the glottal flow.
  // ****************************************************************

  if (getDerivative)
  {
    double c1 = (B*exp(alpha*te)*sin(w*te)) / (epsilon*ta);
    double c2 = exp(-epsilon*(T0 - te));

    for (i=0; i < s.N; i++)
    {
      t = (double)i / (double)s.N;
      if (t <= te)
      {
        s.x[i] = B*exp(alpha*t)*sin(w*t);
      }
      else
      {
        s.x[i] = c1*(exp(-epsilon*(t-te)) - c2);
      }
    }
  }
  else

  // ****************************************************************
  // Calculate the glottal flow waveform.
  // ****************************************************************

  {
    double u1_te = (B*(exp(alpha*te)*(alpha*sin(w*te) - w*cos(w*te)) + w)) / (w*w + alpha*alpha);
    double preFactor = (B*exp(alpha*te)*sin(w*te)*exp(epsilon*te)) / (epsilon*ta);
    double F2_te = preFactor*(-exp(-epsilon*te)/epsilon - te*exp(-epsilon*T0));

	// to generate noise
	std::random_device rd{};
	std::mt19937 randomNumberGenerator{ rd() };
	// inputSample is a random number with the standard deviation
	// 1/sqrt(12) and range limited to [-1.0, 1.0]
	std::normal_distribution<double> normalDistribution(0.0, 1.0 / sqrt(12.0));
	double inputSample;

    for (i=0; i < s.N; i++)
    {
      t = (double)i / (double)s.N;
      if (t <= te)
      {
        s.x[i] = (B*(exp(alpha*t)*(alpha*sin(w*t) - w*cos(w*t)) + w)) / (w*w + alpha*alpha);
      }
      else
      {
        s.x[i] = u1_te + preFactor*(-exp(-epsilon*t)/epsilon - t*exp(-epsilon*T0)) - F2_te;
      }

	  // add noise 
	  inputSample = 0.0;
	  do
	  {
		  inputSample = normalDistribution(randomNumberGenerator);
	  } while ((inputSample < -1.0) || (inputSample > 1.0));
	  s.x[i] *= (1. + pow(10, -SNR / 20.) * inputSample);
    }
  }
}


// ****************************************************************************
/// Calculates the epsilon value in the exponent of the right pulse part.
// ****************************************************************************

double LfPulse::getEpsilon(double ta, double te)
{
  const int MAX_STEPS = 40;
  const double MIN_TA = 0.0001;
  const double MIN_C  = 0.001;
  double c = 1.0 - te;
  if (c < MIN_C) 
  { 
    c = MIN_C; 
  }
  
  if (ta < MIN_TA) 
  { 
    ta = MIN_TA; 
  }
  
  if (ta > c-0.00001) 
  { 
    ta = c-0.00001; 
  }

  // Newton-Iteration

  double epsilon = 1.0 / ta;      // 1st Approximation
  int numSteps = 0;
  double h, h2;

  do
  {
    h  = 1.0 - exp(-epsilon*c) - epsilon*ta;
    h2 = c*exp(-epsilon*c) - ta;
    epsilon = epsilon - h/h2;
  } while ((numSteps < MAX_STEPS) && (fabs(h) > 0.00001));

  return epsilon;  
}


// ****************************************************************************
/// Calculates the alpha value.
// ****************************************************************************

double LfPulse::getAlpha(double tp, double te, double ta, double epsilon)
{
  double w = M_PI / tp;
  double w2 = w*w;
  double SIN = sin(w*te);
  double COS = cos(w*te);
  double X = (SIN*exp(epsilon*te)*(-exp(-epsilon)/epsilon - exp(-epsilon))) / (epsilon*ta);
  double Y = (SIN*exp(epsilon*te)*(-exp(-epsilon*te)/epsilon - te*exp(-epsilon))) / (epsilon*ta);
  double a[2];
  double h[2];
  double newAlpha, newH;

  a[0] = 0.0;
  a[1] = 0.0;

  int numSteps = 0;
  do
  {
    numSteps++;
    a[1]+= 1.0;
    h[1] = w/(w2+a[1]*a[1]) + exp(a[1]*te)*((a[1]*SIN - w*COS)/(w2+a[1]*a[1]) + X - Y);
  } while ((numSteps < 20) && (h[1] >= 0.0));

  // h(a[0]) should now be > 0 and h(a[1]) should now be < 0.

  if (h[1] >= 0)
  {
    return 0.0;
  }

  // Use the approximation algorithm "Regula falsi"

  numSteps = 0;

  do
  {
    h[0] = w/(w2+a[0]*a[0]) + exp(a[0]*te)*((a[0]*SIN - w*COS)/(w2+a[0]*a[0]) + X - Y);
    h[1] = w/(w2+a[1]*a[1]) + exp(a[1]*te)*((a[1]*SIN - w*COS)/(w2+a[1]*a[1]) + X - Y);
    newAlpha = a[0] -  (h[0]*(a[1]-a[0])) / (h[1]-h[0]);
    newH = w/(w2+newAlpha*newAlpha) + exp(newAlpha*te)*((newAlpha*SIN - w*COS)/(w2+newAlpha*newAlpha) + X - Y);

    if (newH < 0.0)
    {
      a[1] = newAlpha;
    }
    else
    {
      a[0] = newAlpha;
    }
    numSteps++;
  } while ((numSteps < 20) && (fabs(newH) > 0.00001));

  return newAlpha;
}


// ****************************************************************************
/// Returns the amplitude of the first part of u'(t) for a given amplitude of
/// u(t).
// ****************************************************************************

double LfPulse::getB(double AMP, double tp, double alpha)
{
  double w = 3.1415926 / tp;
  double help = (exp(alpha*tp)*(alpha*sin(w*tp) - w*cos(w*tp)) + w) / (w*w + alpha*alpha);
  double B = AMP / help;
  return B;
}

// ****************************************************************************
