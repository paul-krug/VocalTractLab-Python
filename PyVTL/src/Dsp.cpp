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

#include "Dsp.h"
#include <cmath>

// Reference frequency for the conversion between Hz and st.
static const double REFERENCE_FREQUENCY = 1.0;


// ****************************************************************************
// Calculates x mod y.
// ****************************************************************************

int modulo(int x, int y)
{
  if (y < 1) { y = 1; }
  
  if (x >= 0)
  {
    return (x % y);
  }

  // x < 0:
  return (y - ((-x) % y));
}

// ****************************************************************************
// Returns the signal energy for the samples between startPos and 
// startPos + numSamples-1.
// ****************************************************************************

double getSignalEnergy(const Signal& signal, int startPos, int numSamples)
{
  int i;
  double energy = 0.0;
  double value;
  
  if (numSamples < 0) { numSamples = 0; }
  int endPos = startPos+numSamples - 1;

  for (i=startPos; i <= endPos ; i++) 
  { 
    value = signal.x[modulo(i, signal.N)];
    energy+= value*value;
  }

  return energy;
}

// ****************************************************************************
// Returns the signal energy for the samples between startPos and 
// startPos + numSamples-1.
// ****************************************************************************

double getSignalEnergy(const Signal16& signal, int startPos, int numSamples)
{
  int i;
  double energy = 0.0;
  double value;
  
  if (numSamples < 0) { numSamples = 0; }
  int endPos = startPos+numSamples - 1;

  for (i=startPos; i <= endPos ; i++) 
  { 
    value = signal.x[modulo(i, signal.N)];
    energy+= value*value;
  }

  return energy;
}

// ****************************************************************************
// Returns the mean signal intensity between startPos and startPos + numSamples-1.
// ****************************************************************************

double getMeanSignalPower(Signal& signal, int startPos, int numSamples)
{
  double energy = getSignalEnergy(signal, startPos, numSamples);
  if (numSamples < 1) { numSamples = 1; }
  double power = energy / numSamples;

  return power;
}


// ****************************************************************************
// Changes real and imaginary parts in the signal into amplitude and phase.
// ****************************************************************************

void rectangularToPolar(ComplexSignal& s, int length)
{
  s.setMinLength(length);
  double re, im;

  for (int i=0; i < length; i++)
  {
    re = s.re[i];
    im = s.im[i];
    s.re[i] = sqrt(re*re + im*im);
    s.im[i] = atan2(im, re);
  }
}

// ****************************************************************************
// Changes amplitude and phase in the signal into real and imaginary parts.
// ****************************************************************************

void polarToRectangular(ComplexSignal& s, int length)
{
  s.setMinLength(length);
  double mag, phase;

  for (int i=0; i < length; i++)
  {
    mag   = s.re[i];
    phase = s.im[i];
    s.re[i] = mag*cos(phase);
    s.im[i] = mag*sin(phase);
  }
}

// ****************************************************************************
// Generates the negative frequencies from the positive frequencies.
// ****************************************************************************

void generateNegativeFrequencies(ComplexSignal *spectrum)
{
  if (spectrum == NULL) { return; }
  int i;
  int N = spectrum->N;

  for (i = N/2+1; i < N; i++)
  {
    spectrum->re[i] = spectrum->re[N-i];
    spectrum->im[i] =-spectrum->im[N-i];
  }
}


// ****************************************************************************
// Calculates the real DFT by simple correlation.
// ****************************************************************************

void realDFT(Signal& timeSignal, ComplexSignal& freqSignal, int length, bool normalize)
{
  // allocate at least N/2+1 values.
  timeSignal.setMinLength(length);
  freqSignal.setMinLength(length/2+1);

  int i, k;
  double angle;
  double re, im;
  int l2 = length/2;

  for (k=0; k <= l2; k++)
  {
    re = 0.0;
    im = 0.0;

    for (i=0; i < length; i++)
    {
      angle = (2.0*M_PI*k*i) / (double)length;
      re+= timeSignal.x[i]*cos(angle);
      // Imaginary part with negative sign for consistency with the complex DFT.
      im-= timeSignal.x[i]*sin(angle); 
    }

    if (normalize)
    {
      im = -im;
      re/= (double)l2;
      im/= (double)l2;
      if ((k == 0) || (k == l2)) { re/= 2; }
    }

    freqSignal.re[k] = re;
    freqSignal.im[k] = im;
  }
}

// ****************************************************************************
// Calculates the real inverse DFT by correlation.
// ****************************************************************************

void realIDFT(ComplexSignal& freqSignal, Signal& timeSignal, int length, bool normalize)
{
  // allocate at least N/2+1 values.
  timeSignal.setMinLength(length);
  freqSignal.setMinLength(length/2+1);

  int i, k;
  double re, im;
  double angle;

  for (i=0; i < length; i++) { timeSignal.x[i] = 0.0; }

  for (k=0; k <= length/2; k++)
  {
    if (normalize)
    {
      im = -freqSignal.im[k] / (double)(length/2);   // Negative sign for imaginary part!
      re =  freqSignal.re[k] / (double)(length/2);
      if ((k == 0) || (k == length/2)) { re/= 2.0; }
    }
    else
    {
      re = freqSignal.re[k];
      im = freqSignal.im[k];
    }

    for (i=0; i < length; i++)
    {
      angle = (2.0*M_PI*k*i) / (double)length;
      timeSignal.x[i]+= re*cos(angle) + im*sin(angle);
    }
  }
}

// ****************************************************************************
// Calc. the complex fast FT of the signal s with the length
// N = 2^lengthExponent.
// The resulting signal is written back to s.
// ****************************************************************************

void complexFFT(ComplexSignal& s, int lengthExponent, bool normalize)
{
  int i, j, k;

  int N = 1 << lengthExponent;
  s.setMinLength(N);

  // Bit reordering.
  double tr, ti;
  int nm1 = N - 1;
  int nd2 = N / 2;

  j = nd2;
  for (i=1; i <= N-2; i++)
  {
    if (i < j)
    {
      tr = s.re[j];
      ti = s.im[j];
      s.re[j] = s.re[i];
      s.im[j] = s.im[i];
      s.re[i] = tr;
      s.im[i] = ti;
    }
    k = nd2;

    while (k <= j)
    {
      j-= k;
      k/= 2;
    }
    j+= k;
  }   // next i


  // ****************************************************************

  int l, le, le2;
  double ur, ui;
  double sr, si;
  int jm1, ip;

  for (l=1; l <= lengthExponent; l++)
  {
    le  = 1 << l;
    le2 = le / 2;
    ur  = 1.0;
    ui  = 0.0;
    sr  = cos(M_PI / (double)le2);
    si  = -sin(M_PI / (double)le2);

    for (j=1; j <= le2; j++)
    {
      jm1 = j - 1;
      for (i=jm1; i <= nm1; i+=le)
      {
        ip = i + le2;
        tr = s.re[ip]*ur - s.im[ip]*ui;
        ti = s.re[ip]*ui + s.im[ip]*ur;
        s.re[ip] = s.re[i] - tr;
        s.im[ip] = s.im[i] - ti;
        s.re[i]+= tr;
        s.im[i]+= ti;
      }   // next i
      
      tr = ur;
      ur = tr*sr - ui*si;
      ui = tr*si + ui*sr;
    }   // next j
  }   // next l

  // Normalize the results? *****************************************

  if (normalize)
  {
    for (i=0; i < N; i++) 
    { 
      s.re[i]/= (double)N;
      s.im[i]/= (double)N;
    }    
  }
}

// ****************************************************************************
// Calc. the inverse complex FFT of the complex Signal s.
// The signal length is 2^lengthExponent.
// ****************************************************************************

void complexIFFT(ComplexSignal& s, int lengthExponent, bool normalize)
{
  int i, k;

  int N = 1 << lengthExponent;
  s.setMinLength(N);

  for (k=0; k < N; k++) { s.im[k] = -s.im[k]; }
 
  complexFFT(s, lengthExponent, normalize);

  for (i=0; i < N; i++) { s.im[i] = -s.im[i]; }
}

// ****************************************************************************
// Calc. the inverse DFT by correlation (slow).
// ****************************************************************************

void complexIDFT(ComplexSignal& freqSignal, ComplexSignal& timeSignal, int length, bool normalize)
{
  int i, k;
  double angle;
  double S, C;

  freqSignal.setMinLength(length);

  timeSignal.reset(length);

  // **************************************************************

  for (i=0; i < length; i++)
  {
    for (k=0; k < length; k++)
    {
      angle = (2.0*M_PI*k*i) / (double)length;
      S = sin(angle);
      C = cos(angle);
      timeSignal.re[i]+= freqSignal.re[k]*C - freqSignal.im[k]*S;
      timeSignal.im[i]+= freqSignal.im[k]*C + freqSignal.re[k]*S;
    }

    if (normalize)
    {
      timeSignal.re[i]/= length;
      timeSignal.im[i]/= length;
    }
  }
}


// ****************************************************************************
// Calc. the DFT for a complex signal by correlation.
// ****************************************************************************


void complexDFT(ComplexSignal& timeSignal, ComplexSignal& freqSignal, int length, bool normalize)
{
  int i, k;
  double angle;
  double S, C;

  timeSignal.setMinLength(length);

  freqSignal.reset(length);

  // ****************************************************************

  for (k=0; k < length; k++)
  {
    for (i=0; i < length; i++)
    {
      angle = (2.0*M_PI*k*i) / (double)length;
      S = sin(angle);
      C = cos(angle);
      freqSignal.re[k]+= timeSignal.re[i]*C + timeSignal.im[i]*S;
      freqSignal.im[k]+= timeSignal.im[i]*C - timeSignal.re[i]*S;
    }

    if (normalize)
    {
      freqSignal.re[k]/= length;
      freqSignal.im[k]/= length;
    }
  }
}


// ****************************************************************************
/// Returns the smallest exponent e for which windowLength <= 2^e.
/// For example, for windowLength = 500, e = 9, because 2^9 = 512.
// ****************************************************************************

int getFrameLengthExponent(int windowLength_pt)
{
  int e = 1;
  while (((int)1 << e) < windowLength_pt)
  {
    e++;
  }

  return e;
}


// ****************************************************************************
// Calc. the fast FT for the real part of signal s.
// The imaginary part of the input is ignored.
// The result is also written into s and includes the negative frequencies.
// N = 2^lengthExponent is the length of the input/output signal.
// This function is about 30% faster than the FFT with the complex time signal.
// ****************************************************************************

void realFFT(ComplexSignal& s, int lengthExponent, bool normalize)
{
  int N = 1 << lengthExponent;
  s.setMinLength(N);

  int i, j, im, ip2, ipm, jm1, ip;

  // ******************************************************
  
  for (i=0; i < N/2; i++)
  {
    s.re[i] = s.re[2*i];
    s.im[i] = s.re[2*i+1];
  }

  // ******************************************************

  complexFFT(s, lengthExponent-1, false); // Do not normalize here yet.

  // ******************************************************

  int nm1 = N - 1;
  int nd2 = N / 2;
  int n4  = (N/4) - 1;

  for (i=1; i <= n4; i++)
  {
    im = nd2 - i;
    ip2 = i + nd2;
    ipm = im + nd2;
    s.re[ip2] =  (s.im[i] + s.im[im]) / 2;
    s.re[ipm] =   s.re[ip2];
    s.im[ip2] = -(s.re[i] - s.re[im]) / 2;
    s.im[ipm] =  -s.im[ip2];

    s.re[i]   = (s.re[i] + s.re[im]) / 2;
    s.re[im]  =  s.re[i];
    s.im[i]   = (s.im[i] - s.im[im]) / 2;
    s.im[im]  = -s.im[i];
  }   // next i

  s.re[(N*3)/4] = s.im[N/4];
  s.re[nd2]     = s.im[0];
  s.im[(N*3)/4] = 0.0;
  s.im[nd2]     = 0.0;
  s.im[N/4]     = 0.0;
  s.im[0]       = 0.0;

  // ******************************************************

  int l   = lengthExponent;
  int le  = 1 << l;
  int le2 = le / 2;
  double ur = 1.0;
  double ui = 0.0;
  double sr = cos(M_PI / (double)le2);
  double si = -sin(M_PI / (double)le2);
  double ti;
  double tr;

  for (j=1; j <= le2; j++)
  {
    jm1 = j - 1;
    for (i=jm1; i <= nm1; i+=le)
    {
      ip = i + le2;
      tr = s.re[ip]*ur - s.im[ip]*ui;
      ti = s.re[ip]*ui + s.im[ip]*ur;
      s.re[ip] = s.re[i] - tr;
      s.im[ip] = s.im[i] - ti;
      s.re[i]+= tr;
      s.im[i]+= ti;
    }   // next i
      
    tr = ur;
    ur = tr*sr - ui*si;
    ui = tr*si + ui*sr;
  }   // next j

  // Normalize the result?

  if (normalize)
  {
    for (i=0; i < N; i++) 
    { 
      s.re[i]/= (double)N;
      s.im[i]/= (double)N;
    }    
  }
}


// ****************************************************************************
// Calc. the fast inverse FT for a real signal.
// N = 2^lengthExponent is the length of the IDFT, and the input signal is 
// filled from index 0 to index N/2 with spectral coefficients.
// The remaining values (N/2+1 .. N-1) automatically added.
// The resulting time signal is in s.re[], and s.im[] is zero.
// ****************************************************************************

void realIFFT(ComplexSignal& s, int lengthExponent, bool normalize)
{
  int N = 1 << lengthExponent;
  s.setMinLength(N);

  int i, k;

  // ******************************************************

  for (k=N/2+1; k < N; k++)
  {
    s.re[k] = s.re[N-k];
    s.im[k] =-s.im[N-k];
  }

  for (k=0; k < N; k++) { s.re[k]+= s.im[k]; }

  realFFT(s, lengthExponent, false);

  // Postprocessing ***************************************
  
  for (i=0; i < N; i++)
  {
    s.re[i] = s.re[i] + s.im[i];
    s.im[i] = 0.0;
  }

  // Normalize the result?

  if (normalize)
  {
    for (i=0; i < N; i++) { s.re[i]/= (double)N; }    
  }
}

// ****************************************************************************
// ****************************************************************************

void getWindow(Signal& window, int length, WindowType type)
{
  int i;
  window.reset(length);

  // ****************************************************************

  if (type == RECTANGULAR_WINDOW)
  {
    for (i=0; i < length; i++) { window.x[i] = 1.0; }
  }
  
  // ****************************************************************

  else
  if (type == HAMMING_WINDOW)
  {
    for (i=0; i < length; i++) 
    { 
      window.x[i] = 0.54 - 0.46*cos((2.0*M_PI*(double)i) / (double)(length-1)); 
    }
  }

  // ****************************************************************

  else
  if (type == RIGHT_HALF_OF_HAMMING_WINDOW)
  {
    for (i=0; i < length; i++) 
    { 
      window.x[i] = 0.54 - 0.46*cos(M_PI + (M_PI*(double)i) / (double)(length-1));
    }
  }

  // ****************************************************************

  else
  if (type == LEFT_HALF_OF_HAMMING_WINDOW)
  {
    for (i=0; i < length; i++) 
    { 
      window.x[i] = 0.54 - 0.46*cos((M_PI*(double)i) / (double)(length-1));
    }
  }

  // ****************************************************************

  else
  if (type == RIGHT_HALF_OF_HANN_WINDOW)
  {
    for (i=0; i < length; i++) 
    { 
      window.x[i] = 0.5 - 0.5*cos(M_PI + (M_PI*(double)i) / (double)(length-1));
    }
  }

  // ****************************************************************

  else
  if (type == GAUSS_WINDOW)
  {
    const double yEdge = 0.01;
    double s = (double)(length*length) / (4.0*log(yEdge));
    for (i=0; i < length; i++) 
    { 
      window.x[i] = exp(((i-length/2)*(i-length/2)) / s);
    }
  }

  // ****************************************************************

  else
  {
    for (i=0; i < length; i++) { window.x[i] = 1.0; }
  }
}


// ----------------------------------------------------------------------------
// Linear Predictive Coding.
// ----------------------------------------------------------------------------


// ****************************************************************************
// Calc. LPC coefficients. coeff[0] is set to 1.
// N is the number of coefficients except coeff[0].
// ****************************************************************************

void getLPCCoefficients(const double *signal, int numSamples, double *coeff, int N)
{
  const int MAX_COEFF = 256; 
  
  int i, j, p;
  double r[MAX_COEFF];
  double alpha[MAX_COEFF];
  double beta[MAX_COEFF];
  double z[MAX_COEFF];      // Reflection coefficients
  double E, q;

  if (N > MAX_COEFF-1) { N = MAX_COEFF-1; }

  // ****************************************************************

  for (i=0; i <= N; i++)
  {
    r[i] = 0.0;
    for (j=0; j < numSamples-i; j++) { r[i]+= signal[j]*signal[j+i]; }
  }

  // Levinson-Durbin algorithm.
  
  E = r[0];
  alpha[0] = 1.0;
  z[0] = 0.0;

  for (p=1; p <= N; p++)
  {
    q = 0.0;
    for (i=0; i < p; i++) { q+= alpha[i]*r[p-i]; }
    if (E == 0.0) { E = 0.0001; }
    z[p] = -q / E;
    alpha[p] = 0.0;
    for (i=0; i <= p; i++) { beta[i] = alpha[i] + z[p]*alpha[p-i]; }
    for (i=0; i <= p; i++) { alpha[i] = beta[i]; }

    E = E*(1.0 - z[p]*z[p]);
  }

  coeff[0] = 1;
  for (i=1; i <= N; i++) { coeff[i] = -alpha[i]; }
}

// ****************************************************************************
// Calc. the LPC residual by inverse filtering.
// ****************************************************************************

void getLPCResidual(const double *signal, double *residual, long l, const double *coeff, long N)
{
  long i, j;
  for (i=0; i < l; i++)
  {
    residual[i] = signal[i];
    for (j=1; j <= N; j++)
    {
      if (i-j >= 0) { residual[i]-= signal[i-j]*coeff[j]; }
    }
  }
}


// ****************************************************************************
// Calc. the original signal from the residual and the coefficients.
// ****************************************************************************

void predictSignal(double *signal, const double *residual, long l, const double *coeff, long N)
{
  long i, j;

  for (i=0; i < l; i++)
  {
    signal[i] = residual[i];
    for (j=1; j <= N; j++)
    {
      if (i-j >= 0) { signal[i]+= coeff[j]*signal[i-j]; }
    }
  }
}

// ****************************************************************************
// Translates the coefficients of the predictor polynom 1-a1*(z^-1) - a2*(z^-2) - 
// ... -aN*(z^-N) by multiplication with z^N into coefficients of the form
// a0*z^N + a1*z^(N-1) + ... + aN
// For this we just need to negate the coefficients 1..N.
// ****************************************************************************

void LPCToPolynomCoefficients(double *LPCCoeff, double *polynomCoeff, long N)
{
  polynomCoeff[0] = LPCCoeff[0];
  long i;
  for (i=1; i <= N; i++)
  {
    polynomCoeff[i] = -LPCCoeff[i];
  }
}

// ----------------------------------------------------------------------------
// Calculation of zeros.
// ----------------------------------------------------------------------------

// ****************************************************************************
// Returns the zeros of the equation x^2+beta*x+gamma.
// ****************************************************************************

void getSquareRoots(double beta, double gamma, ComplexValue &x0, ComplexValue &x1)
{
  double re = -0.5*beta;
  double root = 0.25*beta*beta - gamma;

  // Two complex zeros.
  if (root <= 0.0)
  {
    root = -root;
    double im = sqrt(root);
    x0 = ComplexValue(re, im);
    x1 = ComplexValue(re, -im);
  }
  // Two real zeros.
  else
  {
    root = sqrt(root);
    x0 = ComplexValue(re + root, 0.0);
    x1 = ComplexValue(re - root, 0.0);
  }
}

// ****************************************************************************
// Returns the value of the polynom at the argument x.
// ****************************************************************************

ComplexValue getPolynomValue(double *a, long N, ComplexValue x)
{
  ComplexValue result(0.0, 0.0);
  ComplexValue power(1.0, 0.0);

  long i;

  for (i=N; i >= 0; i--)
  {
    result+= a[i]*power;
    power*= x;
  }

  return result;
}


// ****************************************************************************
// Calculates simultaneously the quadratic factors of the polynom of order N.
// If N is uneven, it is made even by setting a zero at x=0.
// ****************************************************************************

void getPolynomRoots(double *a, int &N, ComplexValue *roots)
{
  const long MAX_M = 128;   // max. Anzahl Quadratfaktoren

  int l, j, i, m;
  bool fertig = false;
  int E;                    // Zähler für noch ungenaue Quadratfaktoren
  double p, q;              // Koeffizienten des akt. Quadratvektors
  double c[2];              // Koeffizienten des Restglieds im Hornerschema
  double temp;
  double u, w;
  double S, T, Sm, Tm, SNew, TNew;
  double h, k;              // Die Korrekturen
  double D;                 // Die Determinante
  double startAngle[MAX_M]; // Wo liegt jeweils der 1. Näherungswert auf dem Einheitskreis ?


  // Falls der Polynomgrad N ungerade ist, dann eine weitere Nullstelle
  // durch Multiplikation des Polynoms mit x hinzufügen ***

  if ((N & 1) == 1)
  {
    N++;
    a[N] = 0.0;
  }

  int M = N / 2;       // Anzahl der Quadratfaktoren
  
  // Faktor[i] = x^2 + beta[i]*x + gamma[i]
  double beta[MAX_M];         // Hier werden die Indizes 1 .. M benutzt
  double gamma[MAX_M];

  // Die komplexen Einheitswurzeln als Startnäherung vorgeben
  for (m=1; m <= M-1; m++)
  {
    startAngle[m] = (M_PI*m) / (double)M;
    beta[m] = 2.0*cos(startAngle[m]);
    gamma[m] = 1.0;
  }
  startAngle[M] = 0.0;
  beta[M] = 0.0;
  gamma[M] = -1;

  // relative Genauigkeitsschranke ************************
  double epsilon = 0.0001;
  double epsilon1 = 2.0*epsilon;
  double epsilon2 = 2.0*epsilon;
  double E1 = 0.0;
  double E2 = 0.0;
  const long MAX = 32;    // max. Anzahl der Iterationsschritte für alle Quadratfaktoren

  fertig = false;
  l = 1;

  while ((l <= MAX) && (!fertig))
  {
    E = 0;                // Zähler für die korrigierten Quadratfaktoren
    j = 1;

    // Alle M Quadratfaktoren durchlaufen
    while (j <= M)
    {
      p = beta[j];
      q = gamma[j];

      // Doppelzeiliges Hornerschema zur Bestimmung von c[0] und c[1]
      c[0] = a[0];
      c[1] = a[1] - p*a[0];
      for (i=2; i <= N; i++)
      {
        temp = c[1];
        c[1] = a[i] - p*c[1] - q*c[0];
        c[0] = temp;
      }
      c[1] = c[1] + p*c[0];

      // Im ersten Durchlauf von l die Ungenauigkeiten aufsummieren
      if (l == 1) { E2+= fabs(c[0]) + fabs(c[1]); }

      // Weitere Korrektur ist noch nötig
      if (fabs(c[0]) + fabs(c[1]) >= epsilon2)
      {
        u = -0.5*p;
        w = u*u - q;
        S = a[0];
        T = 0;

        // Berechnung des Produktes als S + v*T
        for (m=1; m <= M; m++)
        {
          if (m != j)
          {
            Tm = beta[m] - p;
            Sm = u*Tm + gamma[m] - q;
            SNew = S*Sm + w*T*Tm;
            TNew = S*Tm + T*Sm;
            S = SNew;
            T = TNew;
          }
        }

        D = S*S - T*T*w;
        
        // Determinante ist nahe Null => Startwerte abändern und nochmal!
        if (fabs(D) < epsilon)
        {
          startAngle[j]+= 0.012345;  // Winkel auf komplexem Zahlenkreis um ein paar Grad erhöhen
          beta[j] = 2.0*cos(startAngle[j]);
          gamma[j] = 1.0;
          // j wird nicht erhöht !
        }
        else
        // Die Korrekturen für den j-ten Quadratfaktor können berechnet werden
        {
          h = (c[0]*(S - u*T) - T*c[1]) / D;
          k = (c[1]*(S + u*T) + c[0]*T*q) / D;
          beta[j]+= h;
          gamma[j]+= k;

          if (fabs(h) + fabs(k) >= epsilon1) { E++; }

          if (l == 1)
          {
            E1+= fabs(h) + fabs(k);
            epsilon1 = (epsilon*E1) / M;
            epsilon2 = (epsilon*E2) / M;
          }

          j++;    // Auf zum nächsten Quadratfaktor
        }
      }   // weitere Korrektur war nötig
      else { j++; }

    }   // Durchlauf der Quadratfaktoren

    if (E == 0) { fertig = true; }
    l++;

  }   // Iterationsschritte


  // Aus den M beta[j] und gamma[j] die Nullstellen berechnen.

  int nextRoot = 0;
  ComplexValue comp[2];
  ComplexValue re, wurzel;

  for (j=1; j <= M; j++)
  {
    getSquareRoots(beta[j], gamma[j], comp[0], comp[1]);  // Nullstellen des Quadratfaktors
    roots[nextRoot++] = comp[0];
    roots[nextRoot++] = comp[1];
  }
}

// ****************************************************************************
// Calculates simultaneously the quadratic factors of the polynom of order N.
// If N is uneven, it is made even by setting a zero at x=0.
// Returns only the real zeros of the polynom.
// ****************************************************************************

void getRealPolynomRoots(double *a, int &N, double *roots, int& numRealRoots)
{
  const int MAX_M = 128;   // max. Anzahl Quadratfaktoren

  int l, j, i, m;
  bool fertig = false;
  int E;                    // Zähler für noch ungenaue Quadratfaktoren
  double p, q;              // Koeffizienten des akt. Quadratvektors
  double c[2];              // Koeffizienten des Restglieds im Hornerschema
  double temp;
  double u, w;
  double S, T, Sm, Tm, SNew, TNew;
  double h, k;              // Die Korrekturen
  double D;                 // Die Determinante
  double startAngle[MAX_M]; // Wo liegt jeweils der 1. Näherungswert auf dem Einheitskreis ?


  // Falls der Polynomgrad N ungerade ist, dann eine weitere Nullstelle
  // durch Multiplikation des Polynoms mit x hinzufügen ***

  if ((N & 1) == 1)
  {
    N++;
    a[N] = 0.0;
  }

  int M = N / 2;       // Anzahl der Quadratfaktoren
  
  // Faktor[i] = x^2 + beta[i]*x + gamma[i]
  double beta[MAX_M];         // Hier werden die Indizes 1 .. M benutzt
  double gamma[MAX_M];

  // Die komplexen Einheitswurzeln als Startnäherung vorgeben
  for (m=1; m <= M-1; m++)
  {
    startAngle[m] = (M_PI*m) / (double)M;
    beta[m] = 2.0*cos(startAngle[m]);
    gamma[m] = 1.0;
  }
  startAngle[M] = 0.0;
  beta[M] = 0.0;
  gamma[M] = -1;

  // relative Genauigkeitsschranke.

  double epsilon = 0.0001;
  double epsilon1 = 2.0*epsilon;
  double epsilon2 = 2.0*epsilon;
  double E1 = 0.0;
  double E2 = 0.0;
  const long MAX = 32;    // max. Anzahl der Iterationsschritte für alle Quadratfaktoren

  fertig = false;
  l = 1;

  while ((l <= MAX) && (!fertig))
  {
    E = 0;                // Zähler für die korrigierten Quadratfaktoren
    j = 1;

    // Alle M Quadratfaktoren durchlaufen
    while (j <= M)
    {
      p = beta[j];
      q = gamma[j];

      // Doppelzeiliges Hornerschema zur Bestimmung von c[0] und c[1]
      c[0] = a[0];
      c[1] = a[1] - p*a[0];
      for (i=2; i <= N; i++)
      {
        temp = c[1];
        c[1] = a[i] - p*c[1] - q*c[0];
        c[0] = temp;
      }
      c[1] = c[1] + p*c[0];

      // Im ersten Durchlauf von l die Ungenauigkeiten aufsummieren
      if (l == 1) { E2+= fabs(c[0]) + fabs(c[1]); }

      // Weitere Korrektur ist noch nötig
      if (fabs(c[0]) + fabs(c[1]) >= epsilon2)
      {
        u = -0.5*p;
        w = u*u - q;
        S = a[0];
        T = 0;

        // Berechnung des Produktes als S + v*T
        for (m=1; m <= M; m++)
        {
          if (m != j)
          {
            Tm = beta[m] - p;
            Sm = u*Tm + gamma[m] - q;
            SNew = S*Sm + w*T*Tm;
            TNew = S*Tm + T*Sm;
            S = SNew;
            T = TNew;
          }
        }

        D = S*S - T*T*w;
        
        // Determinante ist nahe Null => Startwerte abändern und nochmal!
        if (fabs(D) < epsilon)
        {
          startAngle[j]+= 0.012345;  // Winkel auf komplexem Zahlenkreis um ein paar Grad erhöhen
          beta[j] = 2.0*cos(startAngle[j]);
          gamma[j] = 1.0;
          // j wird nicht erhöht !
        }
        else
        // Die Korrekturen für den j-ten Quadratfaktor können berechnet werden
        {
          h = (c[0]*(S - u*T) - T*c[1]) / D;
          k = (c[1]*(S + u*T) + c[0]*T*q) / D;
          beta[j]+= h;
          gamma[j]+= k;

          if (fabs(h) + fabs(k) >= epsilon1) { E++; }

          if (l == 1)
          {
            E1+= fabs(h) + fabs(k);
            epsilon1 = (epsilon*E1) / M;
            epsilon2 = (epsilon*E2) / M;
          }

          j++;    // Auf zum nächsten Quadratfaktor
        }
      }   // weitere Korrektur war nötig
      else { j++; }

    }   // Durchlauf der Quadratfaktoren

    if (E == 0) { fertig = true; }
    l++;

  }   // Iterationsschritte


  // Aus den M beta[j] und gamma[j] die (reellen) Nullstellen berechnen.

  numRealRoots = 0;
  double re, wurzel;

  for (j=1; j <= M; j++)
  {
    re     = -0.5*beta[j];
    wurzel = 0.25*beta[j]*beta[j] - gamma[j];

    if (wurzel >= 0.0)
    {
      wurzel = sqrt(wurzel);
      roots[numRealRoots++] = re + wurzel;
      roots[numRealRoots++] = re - wurzel;
    }
  }
}

// ****************************************************************************
// ****************************************************************************

double hertzToSemitones(double freq_Hz)
{
  if (freq_Hz < 1.0) { freq_Hz = 1.0; }
  return 12.0*log(freq_Hz/REFERENCE_FREQUENCY) / log(2.0);
}

// ****************************************************************************
// ****************************************************************************

double semitonesToHertz(double freq_st)
{
  return REFERENCE_FREQUENCY*pow(2, freq_st / 12.0);
}


// ****************************************************************************
