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

#include "F0EstimatorYin.h"
#include <limits>
#include <cmath>
#include "Constants.h"

// Static constants.

const double F0EstimatorYin::INTERNAL_TIME_STEP_S = 0.002;    // 2 ms


// ****************************************************************************
/// Constructor.
// ****************************************************************************

F0EstimatorYin::F0EstimatorYin()
{
  int i;

  // ****************************************************************
  // Init the public variables.
  // ****************************************************************

  differenceFunctionThreshold = 0.1;
  timeStep_s = 0.01;

  // ****************************************************************
  // Init the private variables.
  // ****************************************************************

  firstRoiSample = 0;
  numRoiSamples = 0;
  firstChunkSample = 0;

  // ****************************************************************
  // Create a Hann window over the 30 ms region in the center of a
  // frame for the calculation of the rms amplitude.
  // ****************************************************************

  for (i=0; i < FRAME_LENGTH; i++)
  {
    hannWindow[i] = 0.0;
  }

  int windowLength = (int)(0.03*(double)SAMPLING_RATE);
  if (windowLength > FRAME_LENGTH)
  {
    windowLength = FRAME_LENGTH;
  }
  int windowOffset = (FRAME_LENGTH - windowLength) / 2;

  for (i=0; i < windowLength; i++)
  {
    hannWindow[windowOffset + i] = 0.5*(1.0 - cos(2.0*M_PI*i / (windowLength-1)));
  }

  // ****************************************************************
  // Create the band-pass filter for the initial filtering.
  // ****************************************************************
  
  const double UPPER_CUTOFF_FREQ = 1000.0;
  const double LOWER_CUTOFF_FREQ = 60.0;    // 40.0
  
  filter = new IirFilter();
  filter->createChebyshev(UPPER_CUTOFF_FREQ / (double)SAMPLING_RATE, false, 2);
  IirFilter highPass;
  highPass.createChebyshev(LOWER_CUTOFF_FREQ / (double)SAMPLING_RATE, true, 4);
  filter->combineWithFilter(&highPass, true);
}


// ****************************************************************************
/// This function initializes the estimation process and performs the
/// pre-processing operations (filtering).
/// firstRoiSample and numRoiSamples define the region of interest in signal,
/// for which F0 shall be estimated. The rest of the signal is considered
/// as unvoiced.
// ****************************************************************************

void F0EstimatorYin::init(Signal16 *signal, int firstRoiSample, int numRoiSamples)
{
  this->firstRoiSample = firstRoiSample;
  this->numRoiSamples = numRoiSamples;
  firstChunkSample = firstRoiSample;

  int i;
  int signalLength = signal->N;

  // ****************************************************************
  // Filter the whole input signal (independent of the ROI).
  // ****************************************************************

  origSignal.reset(signalLength);
  filteredSignal.reset(signalLength);

  for (i=0; i < signalLength; i++)
  {
    origSignal.x[i] = (double)signal->x[i];
  }
  filterSignal(origSignal.x, filteredSignal.x, signalLength);

  // ****************************************************************
  // Reset all internal frames for the whole signal to represent
  // unvoiced regions.
  // ****************************************************************

  int numFrames = (int)((double)signalLength / (double)(SAMPLING_RATE*INTERNAL_TIME_STEP_S));
  frames.resize(numFrames);
  
  for (i=0; i < numFrames; i++)
  {
    // Always put the unvoiced pitch candidate at the beginning of the
    // candidate list by default, so that it exists also for frames
    // outside the region of interest (needed by Viterbi Search).
    frames[i].numPitchCandidates = 1;
    frames[i].pitchCandidateT0[0] = 0.0;
    frames[i].pitchCandidateY[0] = 1.0;

    frames[i].rmsAmplitude = 0.0;
    frames[i].zeroCrossings = 0;
    frames[i].initialCandidate = -1;
    frames[i].finalCandidate = -1;
  }
}


// ****************************************************************************
/// This function allows to process the F0 estimation in one or more chunks
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

bool F0EstimatorYin::processChunk(int numChunkSamples)
{
  double frame[FRAME_LENGTH];
  double df[INTEGRATION_LENGTH];
  double ndf[INTEGRATION_LENGTH];
  int centerPos_pt;
  int i;

  int lastChunkSample = firstChunkSample + numChunkSamples - 1;
  if (lastChunkSample > firstRoiSample + numRoiSamples - 1)
  {
    lastChunkSample = firstRoiSample + numRoiSamples - 1;    
  }
  int firstFrame = (int)(firstChunkSample / ((double)SAMPLING_RATE*INTERNAL_TIME_STEP_S));
  int lastFrame = (int)((lastChunkSample) / ((double)SAMPLING_RATE*INTERNAL_TIME_STEP_S));

//  printf("Processing frames %d to %d\n", firstFrame, lastFrame);

  // ****************************************************************
  // Calculate one frame every ms to enable step 6 in the YIN-paper.
  // ****************************************************************

  for (i=firstFrame; i <= lastFrame; i++)
  {
    centerPos_pt = (int)(i*INTERNAL_TIME_STEP_S*(double)SAMPLING_RATE);
    getFrameSignal(&filteredSignal, centerPos_pt, frame);
    calcNdf(frame, df, ndf);

    getFrameData(frame, df, ndf, frames[i]);
  }

  // ****************************************************************
  // Increment the internal chunk start position and return true,
  // when the whole ROI was processed.
  // ****************************************************************

  firstChunkSample+= numChunkSamples;
  if (firstChunkSample >= firstRoiSample + numRoiSamples)
  {
    return true;
  }
  else
  {
    return false;
  }
}


// ****************************************************************************
/// After the whole input signal was processed (processChunk(...)==0), this
/// function returns a vector with the F0 estimates every timeStep_s seconds
/// from the beginning of the input signal on.
// ****************************************************************************

vector<double> F0EstimatorYin::finish()
{
  const double EPSILON = 0.0000001;
  int i;
  // The returned vector with F0 values.
  vector<double> f0Values;

  int numF0Values = (int)(origSignal.N / ((double)SAMPLING_RATE * timeStep_s));
  f0Values.resize(numF0Values);
  
  // ****************************************************************
  // Use the Viterbi algorithm to make a V/U-decision for each frame
  // and obtain the optimal sequence of candidates.
  // ****************************************************************

  findBestPitchPath();

  for (i=0; i < numF0Values; i++)
  {
    f0Values[i] = getFinalF0((double)i*timeStep_s);
  }
  
  return f0Values;
}


// ****************************************************************************
/// Returns the samples of the frame of length FRAME_LENGTH around the sample
/// at centerPos in the given signal.
// ****************************************************************************

void F0EstimatorYin::getFrameSignal(Signal16 *signal, int centerPos, double *frame)
{
  int i;
  int offset = INTEGRATION_LENGTH;

  for (i = 0; i < FRAME_LENGTH; i++)
  {
    frame[i] = signal->getValue(centerPos - offset + i);
  }
}


// ****************************************************************************
/// Returns the samples of the frame of length FRAME_LENGTH around the sample
/// at centerPos in the given signal.
// ****************************************************************************

void F0EstimatorYin::getFrameSignal(Signal *signal, int centerPos, double *frame)
{
  int i;
  for (i = 0; i < FRAME_LENGTH; i++)
  {
    frame[i] = signal->getValue(centerPos - INTEGRATION_LENGTH + i);
  }
}


// ****************************************************************************
/// This function calculates the NDF directly via the equation in the appendix 
/// of the YIN paper, and not via the cross-correlation function.
// ****************************************************************************

void F0EstimatorYin::calcNdf(double *frame, double *df, double *ndf)
{
  int tau, k;
  int startPos;
  double sum;
  double d;

  for (tau=0; tau < INTEGRATION_LENGTH; tau++)
  {
    sum = 0.0;
    startPos = (INTEGRATION_LENGTH - 1 - tau) / 2;
    for (k=0; k < INTEGRATION_LENGTH; k++)
    {
      d = frame[startPos + k] - frame[startPos + k + tau];
      sum+= d*d;
    }

    df[tau] = sum;
  }

  // Calculate the normalized difference function

  ndf[0] = 1.0;
  sum = 0.0;
  for (tau=1; tau < INTEGRATION_LENGTH; tau++)
  {
    sum+= df[tau];
    ndf[tau] = df[tau]*tau / sum;
  }
}


// ****************************************************************************
/// Calculates some properties of the current frame, like the pitch candidates,
/// the energy, and the zero-crossing rate.
// ****************************************************************************

void F0EstimatorYin::getFrameData(double *frameSignal, double *df, double *ndf, FrameData &fd)
{
  const double MAX_F0 = 800.0;
  const int MIN_TAU = (int)((double)SAMPLING_RATE / MAX_F0);
  const double MIN_DELTA_Y = 0.1;
  int i;

  // ****************************************************************
  // Find all distinct local minima (with y < 1) in the NDF and
  // take them as pitch candidates.
  // Apply a hysteresis to avoid minima due to noise.
  // ****************************************************************

  // The very first candidate is always at T0=0 and indicates an
  // unvoiced frame.

  fd.pitchCandidateT0[0] = 0.0;
  fd.pitchCandidateY[0] = 1.0;    // Must be 1.0 -> optimal cost for unvoiced decision.
  fd.numPitchCandidates = 1;

  // ****************************************************************

  double dummy;
  double accurateTau;
  double accurateY;
  int maxPos = 0;
  int minPos = 0;

  for (i=0; i < INTEGRATION_LENGTH; i++)
  {
    if ((i > minPos) && (minPos > maxPos) && 
        (ndf[maxPos]-ndf[minPos] > MIN_DELTA_Y) && (ndf[i]-ndf[minPos] > MIN_DELTA_Y) && 
        (fd.numPitchCandidates < MAX_PITCH_CANDIDATES) && 
        (minPos > MIN_TAU) && 
        (ndf[minPos] <= 1.0))
    {
      // Obtain the accurate y-value of the dip from the NDF and the
      // accurate lag-value from the DF by parabolic interpolation 
      // (step 5 in YIN).
      
      fitParabola(ndf, minPos, dummy, accurateY);
      fitParabola(df, minPos, accurateTau, dummy);

      fd.pitchCandidateT0[fd.numPitchCandidates] = accurateTau / (double)SAMPLING_RATE;
      fd.pitchCandidateY[fd.numPitchCandidates] = accurateY;
      fd.numPitchCandidates++;
      maxPos = i;
      minPos = i;
    }

    // minPos must always have equal or greater lag than maxPos.
    if (ndf[i] > ndf[maxPos])
    {
      maxPos = i;
      minPos = i;
    }
    if (ndf[i] < ndf[minPos])
    {
      minPos = i;
    }
  }

  // ****************************************************************
  // Determine the initial best pitch candidate with respect to
  // this single frame.
  // ****************************************************************

  // Find the candidate with the lowest T0 below the the y-threshold.
  
  fd.initialCandidate = -1;
  for (i=0; i < fd.numPitchCandidates; i++)
  {
    if (fd.pitchCandidateY[i] <= differenceFunctionThreshold)
    {
      if ((fd.initialCandidate == -1) || (fd.pitchCandidateT0[i] < fd.pitchCandidateT0[ fd.initialCandidate ]))
      {
        fd.initialCandidate = i;
      }
    }
  }

  // If no such candidate exists, take the candidate with the lowest y of all.

  if ((fd.initialCandidate == -1) && (fd.numPitchCandidates > 0))
  {
    fd.initialCandidate = 0;
    for (i=1; i < fd.numPitchCandidates; i++)
    {
      if (fd.pitchCandidateY[i] < fd.pitchCandidateY[ fd.initialCandidate ])
      {
        fd.initialCandidate = i;
      }
    }
  }

  // ****************************************************************
  // Calculate the rms amplitude and zero crossing rate.
  // ****************************************************************

  fd.rmsAmplitude = 0.0;
  fd.zeroCrossings = 0;

  double d;
  double sum = 0.0;

  for (i=0; i < FRAME_LENGTH; i++)
  {
    d = frameSignal[i]*hannWindow[i];
    sum+= d*d;
  }
  sum/= (double)FRAME_LENGTH;
  fd.rmsAmplitude = sqrt(sum);

}


// ****************************************************************************
/// Returns the best local fundamental period estimate around the time t_s
/// according to step 6 in the YIN-paper.
/// Call this function not before the initial best candidates were determined
/// for each frame in getFrameData(...).
// ****************************************************************************

double F0EstimatorYin::getBestLocalT0Estimate(double t_s)
{
  double T0 = 0.0;

  int centerFrameIndex = (int)(t_s/INTERNAL_TIME_STEP_S + 0.5);
  if (centerFrameIndex < 0)
  {
    centerFrameIndex = 0;
  }
  if (centerFrameIndex > (int)frames.size() - 1)
  {
    centerFrameIndex = (int)frames.size() - 1;
  }

  int numSurroundingFrames = (int)((double)INTEGRATION_LENGTH / ((double)SAMPLING_RATE*INTERNAL_TIME_STEP_S));
  int leftFrameIndex = centerFrameIndex - (numSurroundingFrames - 1) / 2;
  int rightFrameIndex = centerFrameIndex + (numSurroundingFrames - 1) / 2;

  if (leftFrameIndex < 0) 
  {
    leftFrameIndex = 0;
  }
  if (rightFrameIndex > (int)frames.size() - 1)
  {
    rightFrameIndex = (int)frames.size() - 1;
  }


  // ****************************************************************
  // Find the best initial pitch candidate in the immediate 
  // neighbourhood.
  // ****************************************************************

  int i;
  int initialCandidate;
  double bestT0 = 0.0;
  double bestY = 1000000.0;

  for (i=leftFrameIndex; i <= rightFrameIndex; i++)
  {
    initialCandidate = frames[i].initialCandidate;
    if (initialCandidate != -1)
    {
      if (frames[i].pitchCandidateY[initialCandidate] < bestY)
      {
        bestY  = frames[i].pitchCandidateY[initialCandidate];
        bestT0 = frames[i].pitchCandidateT0[initialCandidate];
      }
    }
  }

  // ****************************************************************
  // In the current frame, find the best pitch candidate in the 
  // near region (+/- 20%) around bestT0.
  // ****************************************************************

  double lowerT0Limit = 0.8*bestT0;
  double upperT0Limit = 1.2*bestT0;
  
  // bestT0 and bestY now for the current frame.
  bestT0 = 0.0;
  bestY = 1000000.0;

  FrameData *fd = &frames[centerFrameIndex];

  for (i=0; i < fd->numPitchCandidates; i++)
  {
    if ((fd->pitchCandidateT0[i] >= lowerT0Limit) && 
      (fd->pitchCandidateT0[i] <= upperT0Limit) &&
      (fd->pitchCandidateY[i] < bestY))
    {
      bestT0 = fd->pitchCandidateT0[i];
      bestY  = fd->pitchCandidateY[i];
    }
  }

  // If no pitch candidate was found in the search region for the
  // current frame, return the initial pitch candidate as 
  // fallback solution.

  if (fabs(bestT0) < 0.0000001)
  {
    if (fd->initialCandidate != -1)
    {
      bestT0 = fd->pitchCandidateT0[ fd->initialCandidate ];
    }
  }

  T0 = bestT0;

  return T0;
}


// ****************************************************************************
/// Perform a Viterbi search to find the best path through all pitch candidates
/// in all frames including the V/U-decision.
// ****************************************************************************

void F0EstimatorYin::findBestPitchPath()
{
  int numFrames = (int)frames.size();
  if (numFrames < 1)
  {
    return;
  }

  int i, k, m;
  double localCost;
  double transitionCost;
  double lowestPathCost;
  int bestPrevCandidate;
  int numCurrCandidates;
  int numPrevCandidates;

  // ****************************************************************
  // Initialize the search with the candidates of the first frame.
  // ****************************************************************

  for (k=0; k < frames[0].numPitchCandidates; k++)
  {
    frames[0].lowestPathCost[k] = getLocalCost(0, k);
    frames[0].bestPrevCandidate[k] = -1;    // There is no prev. frame.
  }

  // ****************************************************************
  // Continue the search for the other frames.
  // ****************************************************************

  for (i=1; i < numFrames; i++)
  {
    numCurrCandidates = frames[i].numPitchCandidates;

    for (k=0; k < numCurrCandidates; k++)
    {
      // Get the local cost for the frame i and candidate k
      localCost = getLocalCost(i, k);

      // Find the minimum of the path cost from each of the 
      // candidates of the previous frame.

      numPrevCandidates = frames[i-1].numPitchCandidates;
      lowestPathCost = numeric_limits<double>::max();
      bestPrevCandidate = -1;

      for (m=0; m < numPrevCandidates; m++)
      {
        transitionCost = getTransitionCost(i-1, m, i, k);
        if (frames[i-1].lowestPathCost[m] + transitionCost + localCost < lowestPathCost)
        {
          lowestPathCost = frames[i-1].lowestPathCost[m] + transitionCost + localCost;
          bestPrevCandidate = m;          
        }
      }

      frames[i].lowestPathCost[k] = lowestPathCost;
      frames[i].bestPrevCandidate[k] = bestPrevCandidate;
    }
  }

  // ****************************************************************
  // Propagate back from the last frame to obtain the path with the
  // lowest cost using the back-pointers.
  // ****************************************************************

  // Find the candidate with the lowest path cost in the last frame.
  
  FrameData *fd = &frames[numFrames-1];
  m = 0;   // The best candidate in the last frame.

  for (k=1; k < fd->numPitchCandidates; k++)
  {
    if (fd->lowestPathCost[k] < fd->lowestPathCost[m])
    {
      m = k;
    }
  }

  // Do the back-propagation and assign the finalCandidate member
  // for each frame.

  i = numFrames - 1;
  while (m != -1)
  {
    frames[i].finalCandidate = m;
    m = frames[i].bestPrevCandidate[m];
    i--;
  }

}


// ****************************************************************************
/// Returns the path cost for the transition from (prevFrame, prevCandidate)
/// to (currFrame, currCandidate).
// ****************************************************************************

double F0EstimatorYin::getTransitionCost(int prevFrame, int prevCandidate, int currFrame, int currCandidate)
{
  const double EPSILON = 1.0;
  const double OCTAVE_CHANGE_COST = 2.0;
  const double AMPLITUDE_TRANSITION_COST = 0.3;   //0.5;
  const double FIXED_VOICING_STATE_TRANSITION_COST = 0.2;   // 0.5
  static const int NUM_FRAMES_PER_20_MS = (int)(0.02 / INTERNAL_TIME_STEP_S);

  double cost = 0.0;

  // ****************************************************************
  // Transition between two potentially voiced frames.
  // ****************************************************************

  if ((prevCandidate > 0) && (currCandidate > 0))
  {
    // The cost is proportional to the absolute difference of the
    // pitch values in semitones.

    double prevT0 = frames[prevFrame].pitchCandidateT0[prevCandidate];
    double currT0 = frames[currFrame].pitchCandidateT0[currCandidate];

    cost = OCTAVE_CHANGE_COST*fabs(log(prevT0 / currT0) / log(2.0));
  }
  else

  // ****************************************************************
  // Transition from a voiced to an unvoiced frame or
  // from an unvoiced to a voiced frame.
  // ****************************************************************

  if (((prevCandidate > 0) && (currCandidate == 0)) || 
      ((prevCandidate == 0) && (currCandidate > 0)))
  {
    int rightFrame = currFrame + NUM_FRAMES_PER_20_MS / 2;
    int leftFrame = rightFrame - NUM_FRAMES_PER_20_MS;
    if (leftFrame < 0)
    {
      leftFrame = 0;
    }
    if (rightFrame >= (int)frames.size())
    {
      rightFrame = (int)frames.size() - 1;
    }

    double ratio = frames[rightFrame].rmsAmplitude / (frames[leftFrame].rmsAmplitude + EPSILON);

    // Transition from an unvoiced to a voiced frame.
    if (prevCandidate == 0)
    {
      cost = AMPLITUDE_TRANSITION_COST / ratio;
    }
    else
    // Transition from a voiced to an unvoiced frame.
    {
      cost = AMPLITUDE_TRANSITION_COST * ratio;
    }
    cost+= FIXED_VOICING_STATE_TRANSITION_COST;

  }
  else

  // ****************************************************************
  // Transition between two potentially unvoiced frames.
  // ****************************************************************
  {
    cost = 0;     // no cost here
  }

  return cost;
}


// ****************************************************************************
/// Returns the local cost of the given candidate in the given frame for the
/// cheapest path of pitch values.
// ****************************************************************************

double F0EstimatorYin::getLocalCost(int frameIndex, int candidateIndex)
{
  const double INFINITY_COST = 1000000.0;    // Extremely high cost!
  const double VOICE_RMS_THRESHOLD = 100.0;
  double cost = 0.0;
  FrameData *fd = &frames[frameIndex];
  int i;

  // ****************************************************************
  // The first pitch candidate is always the candidate for the 
  // unvoiced state (T0=0).
  // ****************************************************************

  if (candidateIndex == 0)
  {
    // The cost is one minus the the y-value of the lowest dip.
    // First, find the y-value of the lowest dip.
    double minY = 1.0;
    for (i=0; i < fd->numPitchCandidates; i++)
    {
      if (fd->pitchCandidateY[i] < minY)
      {
        minY = fd->pitchCandidateY[i];
      }
    }
    cost = 1.0 - minY;
  }
  else

  // ****************************************************************
  // Evaluate the cost for a potentially voiced frame.
  // ****************************************************************

  {
    // Find the candidate with the lowest T0 of those below the the y-threshold.

    int k = -1;
    for (i=0; i < fd->numPitchCandidates; i++)
    {
      if (fd->pitchCandidateY[i] <= differenceFunctionThreshold)
      {
        if ((k == -1) || (fd->pitchCandidateT0[i] < fd->pitchCandidateT0[k]))
        {
          k = i;
        }
      }
    }

    // There is no perfect pitch candidate.
    if (k == -1)
    {
      cost = fd->pitchCandidateY[candidateIndex];
    }
    else
    // The perfect pitch candidate is k.
    {
      if (candidateIndex == k)
      {
        // Prefer the perfect candidate with a reduced cost !
        cost = fd->pitchCandidateY[candidateIndex] - 0.1;
      }
      else
      {
        cost = fd->pitchCandidateY[candidateIndex];
      }
    }

    // The rms amplitude must be above VOICE_RMS_THRESHOLD for the
    // frame to be classified as voiced.
    if (fd->rmsAmplitude < VOICE_RMS_THRESHOLD)
    {
      cost = INFINITY_COST;
    }
  }



cost = 0.2*cost;    // 0.2


  return cost;
}


// ****************************************************************************
/// Returns the final F0 value at the time t_s based on the finalCandidate
/// estimates in the frames.
/// When 0 is returned, the frame is considered as unvoiced.
// ****************************************************************************

double F0EstimatorYin::getFinalF0(double t_s)
{
  if (frames.size() < 1)
  {
    return 0.0;
  }

  // Which is the frame closest to t_s ?

  int frameIndex = (int)(t_s/INTERNAL_TIME_STEP_S + 0.5);
  if (frameIndex < 0)
  {
    frameIndex = 0;
  }
  if (frameIndex > (int)frames.size() - 1)
  {
    frameIndex = (int)frames.size() - 1;
  }

  int finalCandidate = frames[frameIndex].finalCandidate;
  if (finalCandidate == -1)
  {
    return 0.0;
  }

  const double EPSILON = 0.0000001;
  double T0 = frames[frameIndex].pitchCandidateT0[finalCandidate];
  double f0 = 0.0;

  if (fabs(T0) < EPSILON)
  {
    f0 = 0.0;
  }
  else
  {
    f0 = 1.0 / T0;
  }

  return f0;
}


// ****************************************************************************
/// Fits a parabola to the local minimum at the sample index rawTau in the
/// function f, and returns the tau and y of the minimum of the parabola.
// ****************************************************************************

void F0EstimatorYin::fitParabola(double *f, int rawTau, double &accurateTau, double &accurateY)
{
  const int SEARCH_RANGE = 10;

  // Default values
  accurateTau = rawTau;
  accurateY = f[rawTau];

  // If there is no local minimum at the position rawTau we search
  // in the near surrounding.

  int counter = 0;
  int increment = 0;
  int pos = rawTau;
  bool ok = false;

  while ((!ok) && (counter < SEARCH_RANGE))
  {
    pos+= increment;
    if ((pos > 0) && (pos < INTEGRATION_LENGTH-1))
    {
      if ((f[pos] <= f[pos-1]) && (f[pos] < f[pos+1])) 
      { 
        ok = true;
        double s = f[pos+1] - f[pos-1];
        double t = f[pos-1] - 2.0*f[pos] + f[pos+1];
        accurateTau = pos - (0.5*s)/t;
        accurateY = f[pos] - (s*s)/(8.0*t);
      }
    }

    if (!ok)
    {
      counter++;
      if (increment >= 0) 
      { 
        increment = -increment-1; 
      } 
      else 
      { 
        increment = -increment+1; 
      }
    }
  }
}


// ****************************************************************************
/// Apply a band-pass filter between 40 Hz and 1 kHz to the input signal.
// ****************************************************************************

void F0EstimatorYin::filterSignal(double *inputSignal, double *outputSignal, int N)
{
  filter->resetBuffers(inputSignal[0]);
  int i;

  for (i=0; i < N; i++)
  {
    outputSignal[i] = filter->getOutputSample(inputSignal[i]);
  }
}

// ****************************************************************************

