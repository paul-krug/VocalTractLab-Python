#include "TamModelF0.h"
#include "CdlpFilter.h"


TamModelF0::TamModelF0(const BoundaryVector& bounds, const double onsetValue) : m_bounds(bounds)
{
	// determine onset
	m_onset.time = bounds[0];
	m_onset.value = onsetValue;

	// determine target durations
	for (std::vector<double>::const_iterator bd = bounds.begin() + 1; bd != bounds.end(); ++bd)
	{
		double duration = *bd - *(bd - 1);
		PitchTarget pt = { 0.0, 0.0, 0.0, duration };
		m_targets.push_back(pt);
	}
}

void TamModelF0::setPitchTargets(const TargetVector& targets)
{
	m_targets = targets;
}

void TamModelF0::setBoundaries(const BoundaryVector& boundaries)
{
	m_bounds = boundaries;
	m_onset.time = boundaries.at(0);
}

TimeSignal TamModelF0::calculateF0(const double samplingPeriod) const
{
	// get length of signal
	double start = m_onset.time;
	double end = start;
	for (std::vector<PitchTarget>::const_iterator pt = m_targets.begin(); pt != m_targets.end(); ++pt)
	{
		end += pt->duration;
	}

	// sampling
	SampleTimes times;
	for (double t = start; t <= end; t += samplingPeriod)
	{
		times.push_back(t);
	}

	// return sampled f0
	return applyFilter(times);
}

TimeSignal TamModelF0::calculateF0(const SampleTimes& times) const
{
	// return sampled f0
	return applyFilter(times);
}

TargetVector TamModelF0::getPitchTargets() const
{
	return m_targets;
}

BoundaryVector TamModelF0::getBoundaries() const
{
	return m_bounds;
}

Sample TamModelF0::getOnset() const
{
	return m_onset;
}

void TamModelF0::setOnset( const double onsetValue )
{
	m_onset.value = onsetValue;
}

TimeSignal TamModelF0::applyFilter(const SampleTimes& times) const
{
	CdlpFilter lowPass(FILTERORDER);
	return lowPass.response(times, m_targets, m_onset);
}