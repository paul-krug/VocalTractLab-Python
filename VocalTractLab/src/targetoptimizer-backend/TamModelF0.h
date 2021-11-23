#pragma once
#include <vector>
#include "Data.h"


class TamModelF0 {
public:
	// constructors
	TamModelF0(const BoundaryVector& bounds, const double onsetValue);

	// public member functions
	void setPitchTargets(const TargetVector& targets);
	void setBoundaries(const BoundaryVector& boundaries);
	TimeSignal calculateF0(const double samplingPeriod) const;
	TimeSignal calculateF0(const SampleTimes& times) const;

	TargetVector getPitchTargets() const;
	BoundaryVector getBoundaries() const;
	Sample getOnset() const;
	void setOnset( const double onsetValue );

private:
	// private member functions
	TimeSignal applyFilter(const SampleTimes& times) const;

	// data members
	Sample m_onset;
	TargetVector m_targets;
	BoundaryVector m_bounds;
};
