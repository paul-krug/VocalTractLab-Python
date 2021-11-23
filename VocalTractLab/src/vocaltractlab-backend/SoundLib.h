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

#ifndef __SOUND_LIB__
#define __SOUND_LIB__

#include <string>
#include "Signal.h"
#include "Dsp.h"

using namespace std;

// ****************************************************************************

class SoundInterface
{
public:
  virtual ~SoundInterface() {}
  virtual bool init(int samplingRate) = 0;
  virtual bool startPlayingWave(signed short *data, int numSamples, bool loop = false) = 0;
  virtual bool stopPlaying() = 0;
  virtual bool startRecordingWave(signed short *data, int numSamples) = 0;
  virtual bool stopRecording() = 0;
  static SoundInterface *getInstance();
};


// ****************************************************************************
// deprecated interface (for backward compatibility reasons)
// ****************************************************************************

inline bool initSound(int samplingRate)
{
  return SoundInterface::getInstance()->init(samplingRate);
}

inline bool waveStartPlaying(signed short *data, int numSamples, bool loop)
{
  return SoundInterface::getInstance()->startPlayingWave(data, numSamples, loop);
}

inline bool waveStopPlaying()
{
  return SoundInterface::getInstance()->stopPlaying();
}

inline bool waveStartRecording(signed short *data, int numSamples)
{
  return SoundInterface::getInstance()->startRecordingWave(data, numSamples);
}

inline bool waveStopRecording()
{
  return SoundInterface::getInstance()->stopRecording();
}

// ****************************************************************************

#endif
