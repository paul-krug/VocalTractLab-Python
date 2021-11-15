#pragma once
#include <vector>

// global constants
const unsigned FILTERORDER(5);			// TAM
//const unsigned RANDOMITERATIONS(20);		// BOBYQA
const double SAMPLERATE{ 44100.0 / 110.0 }; 		// Hz
// vector of syllable bounds
using BoundaryVector = std::vector<double>;

// sample of a discrete time signal
struct Sample
{
	double time;
	double value;
};

// discrete time signal types
using SampleTimes = std::vector<double>;
using TimeSignal = std::vector<Sample>;

// pitchTarget according to the TAM
struct PitchTarget
{
	double slope;
	double offset;
	double tau; // time constant
	double duration;
};

// vector of pitch targets
using TargetVector = std::vector<PitchTarget>;

struct SignalStat
{
	double minTime;
	double maxTime;
	double minValue;
	double maxValue;
};

class Data // Singleton design pattern: Central object to store all data. No more than one instance allowed.
{
public:
	static Data& getInstance() // Call to get the only existing instance
	{
		static Data instance;
		return instance;
	}
	void reset()
	{
		initialBoundaries.clear();
		optimalBoundaries.clear();
		pitchTargets.clear();
		optimalF0.clear();
		originalF0.clear();
	}


public:
	BoundaryVector initialBoundaries;
	BoundaryVector optimalBoundaries;
	std::vector<PitchTarget> pitchTargets;
	Sample onset;
	TimeSignal optimalF0;
	TimeSignal originalF0;

private:
	Data() = default;
	~Data() = default;
	Data(const Data&) = delete;
	Data& operator=(const Data&) = delete;
};