#pragma once
#include <vector>
#include "Data.h"


// vector types for CdlpFilter
using FilterState = std::vector<double>;
using FilterCoefficients = std::vector<double>;

// Nth order critical damped low pass filter for target approximation
class CdlpFilter {
public:
	// constructors
	CdlpFilter(const unsigned order = FILTERORDER) : m_filterOrder(order) {};

	// public member functions
	TimeSignal response(const SampleTimes& sampleTimes, const TargetVector& targets, const Sample onset) const;

private:
	// private member functions
	FilterCoefficients calculateCoefficients(const PitchTarget& target, const FilterState& state) const;
	FilterState calculateState(const FilterState& state, const double time, const double startTime, const PitchTarget& target) const;
	static double binomial(const unsigned n, const unsigned k);
	static double factorial(unsigned n);

	// data members
	unsigned m_filterOrder;
};