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

#include "VoiceQualityEstimator.h"
#include "Constants.h"
#include <cstdio>

const double VoiceQualityEstimator::SLICE_STEP_S = 0.01;    // = 10 ms
const double VoiceQualityEstimator::MIN_PEAK_SLOPE = -10.0;
const double VoiceQualityEstimator::MAX_PEAK_SLOPE = 0.0;

// ****************************************************************************
// ****************************************************************************

VoiceQualityEstimator::VoiceQualityEstimator()
{
  timeStep_s = 0.01;    // 10 ms

  firstRoiSlice = 0;
  numRoiSlices = 0;
  nextSlice = 0;

  // Init the wavelets for the six center frequencies.
  calcWavelet(wavelet8000, 1);
  calcWavelet(wavelet4000, 2);
  calcWavelet(wavelet2000, 4);
  calcWavelet(wavelet1000, 8);
  calcWavelet(wavelet500, 16);
  calcWavelet(wavelet250, 32);
}


// ****************************************************************************
/// This function initializes the estimation process and performs the
/// pre-processing operations.
// ****************************************************************************

void VoiceQualityEstimator::init(Signal16 *signal, int firstRoiSample, int numRoiSamples)
{
  int i;
  int signalLength = signal->N;

  // ****************************************************************
  // Set the size of the slice vector and init the slice vector.
  // ****************************************************************

  slices.resize((unsigned int)(((double)signalLength / (double)SAMPLING_RATE) / SLICE_STEP_S));
  for (i=0; i < (int)slices.size(); i++)
  {
    slices[i].peak250 = 0.0;
    slices[i].peak500 = 0.0;
    slices[i].peak1000 = 0.0;
    slices[i].peak2000 = 0.0;
    slices[i].peak4000 = 0.0;
    slices[i].peak8000 = 0.0;
  }

  // ****************************************************************
  // Set the region of interest.
  // ****************************************************************

  firstRoiSlice = (int)(((double)firstRoiSample / (double)SAMPLING_RATE) / SLICE_STEP_S);
  int lastRoiSlice = (int)(((double)(firstRoiSample + numRoiSamples) / (double)SAMPLING_RATE) / SLICE_STEP_S);
  
  firstRoiSlice-= 3;
  if (firstRoiSlice < 0)
  {
    firstRoiSlice = 0;
  }

  numRoiSlices = lastRoiSlice + 3 - firstRoiSlice + 1;
  if (firstRoiSlice + numRoiSlices > (int)slices.size())
  {
    numRoiSlices = (int)slices.size() - firstRoiSlice;    
  }

  // Next slice to be processed.
  nextSlice = firstRoiSlice;


  // ****************************************************************
  // Make an internal copy of the original signal.
  // ****************************************************************

  origSignal.reset(signalLength);

  for (i=0; i < signalLength; i++)
  {
    origSignal.x[i] = (double)signal->x[i];
  }

}


// ****************************************************************************
/// This function allows to process the VQ estimation in one or more chunks
/// of data such that the calling GUI can stay reactive.
/// Call this function one or several times after init(...) to perform the 
/// time-consuming operations on the frames for the next chunk of samples.
/// When all samples of the region of interest in the signal were processed
/// by one or more chunk operations, this function returns true.
///
/// \return true, if all samples in the region of interest were processed.
/// If false is returned, processChunk(...) should be called again for the
/// remaining samples, before F0 estimation is finished with finish(...).
// ****************************************************************************

bool VoiceQualityEstimator::processChunk(int numChunkSamples)
{
  const int NUM_CHUNCK_SLICES = 1;
  int i;
  bool finished = false;

  for (i=0; (i < NUM_CHUNCK_SLICES) && (finished == false); i++)
  {
    if ((nextSlice < (int)slices.size()) && (nextSlice < firstRoiSlice + numRoiSlices))
    {
      calcSlicePeaks(nextSlice);
      nextSlice++;
    }
    else
    {
      finished = true;
    }    
  }

  return finished;
}


// ****************************************************************************
/// After the whole input signal was processed (processChunk(...)==0), this
/// function returns a vector with the voice quality estimates every timeStep_s 
/// seconds from the beginning of the input signal on.
// ****************************************************************************

vector<double> VoiceQualityEstimator::finish()
{
  vector<double> vqValues;
  int i;
  double centerTime_s;

  int numVqValues = (int)(origSignal.N / ((double)SAMPLING_RATE * timeStep_s));
  vqValues.resize(numVqValues);

  for (i=0; i < numVqValues; i++)
  {
    centerTime_s = (double)i*timeStep_s;
    vqValues[i] = calcPeakSlope(centerTime_s);

//    vqValues[i] = MIN_PEAK_SLOPE + (MAX_PEAK_SLOPE - MIN_PEAK_SLOPE)*
  //    0.5*(1.0 + cos(2.0*M_PI*i/10.0 + ((double)rand()/(double)RAND_MAX)*2.0*M_PI));
  }

  return vqValues;
}


// ****************************************************************************
// ****************************************************************************

void VoiceQualityEstimator::printData(int pos_pt)
{
  int sliceIndex = (int)(((double)pos_pt / (double)SAMPLING_RATE) / SLICE_STEP_S);
  if ((sliceIndex >= 0) && (sliceIndex < (int)slices.size()))
  {
    Slice *s = &slices[sliceIndex];
    printf("slice no. = %d : %2.2f  %2.2f  %2.2f  %2.2f  %2.2f  %2.2f\n",
      sliceIndex, s->peak250, s->peak500, s->peak1000, s->peak2000, s->peak4000, s->peak8000);
  }
}


// ****************************************************************************
/// Calculates the peak slope parameter for the signal part centered around
/// the given time.
// ****************************************************************************

double VoiceQualityEstimator::calcPeakSlope(double centerTime_s, bool debug)
{
  // Take the peaks in seven consecutive slices as one frame (= 70 ms).
  int centerSlice = (int)(centerTime_s / SLICE_STEP_S);
  int firstSlice = centerSlice - 3;
  int lastSlice = centerSlice + 3;
  if (firstSlice < 0)
  {
    firstSlice = 0;
  }
  if (lastSlice >= (int)slices.size())
  {
    lastSlice = (int)slices.size() - 1;    
  }

  Slice frame;
  Slice *s;
  int i;

  frame.peak250 = 0.0;
  frame.peak500 = 0.0;
  frame.peak1000 = 0.0;
  frame.peak2000 = 0.0;
  frame.peak4000 = 0.0;
  frame.peak8000 = 0.0;

  for (i=firstSlice; i <= lastSlice; i++)
  {
    s = &slices[i];
    
    if (s->peak250 > frame.peak250)
    {
      frame.peak250 = s->peak250;
    }

    if (s->peak500 > frame.peak500)
    {
      frame.peak500 = s->peak500;
    }

    if (s->peak1000 > frame.peak1000)
    {
      frame.peak1000 = s->peak1000;
    }

    if (s->peak2000 > frame.peak2000)
    {
      frame.peak2000 = s->peak2000;
    }

    if (s->peak4000 > frame.peak4000)
    {
      frame.peak4000 = s->peak4000;
    }

    if (s->peak8000 > frame.peak8000)
    {
      frame.peak8000 = s->peak8000;
    }
  }

  // ****************************************************************
  // Make the linear regression with the six points.
  // ****************************************************************

  const double EPSILON = 0.000000001;
  const int NUM_POINTS = 5;   // 6
  double x[NUM_POINTS] = { 0.0, 1.0, 2.0, 3.0, 4.0 }; //, 5.0 };
  double y[NUM_POINTS] = { frame.peak250, frame.peak500, frame.peak1000, frame.peak2000, frame.peak4000 }; //, frame.peak8000 };

  // Calc. the log10(.) of the amplitudes to make the peak slope
  // value resistant against amplitude changes of the audio signal.

  for (i=0; i < NUM_POINTS; i++)
  {
    if (y[i] < EPSILON)
    {
      y[i] = EPSILON;
    }
    y[i] = 20.0*log10(y[i]);
  }
    
  // Calc. the mean value of the x- and y-data.

  double meanX = 0.0;
  double meanY = 0.0;
  for (i=0; i < NUM_POINTS; i++)
  {
    meanX+= x[i];
    meanY+= y[i];
  }
  meanX/= NUM_POINTS;
  meanY/= NUM_POINTS;

  double numerator = 0.0;
  double denominator = 0.0;

  for (i=0; i < NUM_POINTS; i++)
  {
    numerator+= x[i]*y[i];
    denominator+= x[i]*x[i];
  }

  numerator-= NUM_POINTS*meanX*meanY;
  denominator-= NUM_POINTS*meanX*meanX;

  double result = 0.0;    // Default value
  if (fabs(denominator) > EPSILON)
  {
    result = numerator / denominator;    
  }

  // Scale down the results according to the original paper!
//  result/= 32768.0;


  if (debug)
  {
    printf("slices %d...%d : %2.2f  %2.2f  %2.2f  %2.2f  %2.2f  reg=%2.2f\n",
    firstSlice, lastSlice, y[0], y[1], y[2], y[3], y[4], result);
  }


  return result;
}


// ****************************************************************************
/// Calculates the maximum in each of the six frequency bands within the given
/// time slice of 10 ms.
// ****************************************************************************

void VoiceQualityEstimator::calcSlicePeaks(int sliceIndex)
{
  if ((sliceIndex < 0) || (sliceIndex >= (int)slices.size()))
  {
    return;
  }

  double startTime_s = SLICE_STEP_S*(double)sliceIndex;
  double endTime_s   = SLICE_STEP_S*(double)(sliceIndex + 1);

  int firstSample = (int)(startTime_s*SAMPLING_RATE);
  int lastSample  = (int)(endTime_s*SAMPLING_RATE) - 1;
  int i;
  double d;
  Slice *s = &slices[sliceIndex];

  s->peak250 = 0.0;
  s->peak500 = 0.0;
  s->peak1000 = 0.0;
  s->peak2000 = 0.0;
  s->peak4000 = 0.0;
  s->peak8000 = 0.0;

  for (i=firstSample; i <= lastSample; i++)
  {
    d = getFilteredSample(i, &wavelet250);
    if (d > s->peak250)
    {
      s->peak250 = d;
    }

    d = getFilteredSample(i, &wavelet500);
    if (d > s->peak500)
    {
      s->peak500 = d;
    }

    d = getFilteredSample(i, &wavelet1000);
    if (d > s->peak1000)
    {
      s->peak1000 = d;
    }

    d = getFilteredSample(i, &wavelet2000);
    if (d > s->peak2000)
    {
      s->peak2000 = d;
    }

    d = getFilteredSample(i, &wavelet4000);
    if (d > s->peak4000)
    {
      s->peak4000 = d;
    }

    d = getFilteredSample(i, &wavelet8000);
    if (d > s->peak8000)
    {
      s->peak8000 = d;
    }
  }
}


// ****************************************************************************
/// Returns the sample value of the original signal filtered with the given
/// wavelet at the given position.
// ****************************************************************************

double VoiceQualityEstimator::getFilteredSample(int pos_pt, Signal *wavelet)
{
  double sum = 0.0;
  int startPos = pos_pt - wavelet->N/2;
  int N = wavelet->N;
  int i;

  // Condition at the beginning and end of the signal.
  if ((startPos < 0) || (startPos + N > origSignal.N))
  {
    return 0.0;
  }

  for (i=0; i < N; i++)
  {
    sum+= origSignal.x[startPos + i] * wavelet->x[i];
  }

  return sum;
}


// ****************************************************************************
/// Calculates a symmetrical wavelet with the given length factor (=1,2,4,8,16,
/// 32) with respect to the 8 kHz mother wavelet.
// ****************************************************************************

void VoiceQualityEstimator::calcWavelet(Signal &wavelet, int lengthFactor)
{
  // Center frequency of the mother wavelet filter.
  const double F_MAX = 8000.0;
  // Half of the periods to the left and half to the right (wavelet is symmetric)
  const int NUM_PERIODS = 6;
  double periodLength_s = lengthFactor / F_MAX;
  double waveletLength_s = NUM_PERIODS * periodLength_s;
  int waveletLength_pt = (int)((double)SAMPLING_RATE * waveletLength_s);
  double tau_s = 0.5 / F_MAX;
  int i;
  double t_s;

  wavelet.reset(waveletLength_pt);

  for (i=0; i < waveletLength_pt; i++)
  {
    t_s = (double)(-waveletLength_pt/2 + i) / (double)SAMPLING_RATE;
    t_s/= lengthFactor;

    wavelet.x[i] = -cos(2.0*M_PI*F_MAX*t_s) * exp(-t_s*t_s / (2.0*tau_s*tau_s));
  }
}

// ****************************************************************************
