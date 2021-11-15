#include "CdlpFilter.h"
#include <algorithm>
#include <cmath>
#include <sstream>
#include "dlib/error.h"

TimeSignal CdlpFilter::response(const SampleTimes& sampleTimes, const TargetVector& targets, const Sample onset) const
{
	// local container
	TimeSignal f0;

	// keep state at syllable bound
	FilterState currentState;
	currentState.push_back(onset.value);
	for (unsigned i = 1; i < m_filterOrder; ++i)
	{
		currentState.push_back(0.0);
	}

	// keep index of current sample
	unsigned sampleIndex(0);

	// track syllable bounds of each target
	double bBegin = onset.time;
	double bEnd = bBegin;

	for (std::vector<PitchTarget>::const_iterator pt = targets.begin(); pt != targets.end(); ++pt)
	{
		// update bounds
		bBegin = bEnd;
		bEnd = bBegin + pt->duration;

		// filter coefficients
		FilterCoefficients c = calculateCoefficients(*pt, currentState);

		while (sampleTimes[sampleIndex] <= bEnd)
		{
			double acc(0.0);
			double t = sampleTimes[sampleIndex] - bBegin;	// current samplePoint, time shift
			for (unsigned n = 0; n < m_filterOrder; ++n)
			{
				acc += (c[n] * std::pow(t, n));
			}

			double time = sampleTimes[sampleIndex];
			double value = acc * std::exp(-(1000.0 / pt->tau) * t) + pt->slope * t + pt->offset;

			Sample s = { time,value };
			f0.push_back(s);
			sampleIndex++;

			if (sampleIndex >= sampleTimes.size())
			{
				// We have calculated the f0 for all points of interest, so return
				return f0;
			}
		}

		// update filter state
		currentState = calculateState(currentState, bEnd, bBegin, *pt);
	}

	return f0;
}

FilterCoefficients CdlpFilter::calculateCoefficients(const PitchTarget& target, const FilterState& state) const
{
	FilterCoefficients coeffs(state.size(), 0.0);

	if (coeffs.size() != m_filterOrder)
	{
		std::ostringstream msg;
		msg << "[calculateCoefficients] Wrong size of coefficient vector! " << coeffs.size() << " != " << m_filterOrder;
		throw dlib::error(msg.str());
	}

	coeffs[0] = state[0] - target.offset;	// 0th coefficient
	for (unsigned n = 1; n < m_filterOrder; ++n)	// other coefficients
	{
		double acc(0.0);
		for (unsigned i = 0; i < n; ++i)
		{
			acc += (coeffs[i] * std::pow(-(1000.0 / target.tau), n - i) * binomial(n, i) * factorial(i));
		}

		if (n == 1)
		{
			acc += target.slope; // adaption for linear targets; minus changes in following term!
		}

		coeffs[n] = (state[n] - acc) / factorial(n);
	}

	return coeffs;
}

FilterState CdlpFilter::calculateState(const FilterState& state, const double time, const double startTime, const PitchTarget& target) const
{
	// setup
	double t(time - startTime); // sample time
	const unsigned& N(m_filterOrder);
	FilterState stateUpdate(N);

	// filter coefficients
	FilterCoefficients c = calculateCoefficients(target, state);

	for (unsigned n = 0; n < N; ++n)
	{
		// calculate value of nth derivative
		double acc(0.0);
		for (unsigned i = 0; i < N; ++i)
		{
			// pre-calculate q-value
			double q(0.0);
			for (unsigned k = 0; k < std::min(N - i, n + 1); ++k)
			{
				q += (std::pow(-(1000.0 / target.tau), n - k) * binomial(n, k) * c[i + k] * factorial(k + i) / factorial(i));
			}

			acc += (std::pow(t, i) * q);
		}

		stateUpdate[n] = acc * std::exp(-(1000.0 / target.tau) * t);
	}

	// correction for linear targets
	if (N > 1)
	{
		stateUpdate[0] += (target.offset + target.slope * t);
	}
	if (N > 2)
	{
		stateUpdate[1] += target.slope;
	}

	return stateUpdate;
}

double CdlpFilter::binomial(const unsigned n, const unsigned k)
{
	double result = 1;
	unsigned int tmp = k;

	if (tmp > n - tmp)
		tmp = n - tmp;

	for (unsigned i = 0; i < tmp; ++i)
	{
		result *= (n - i);
		result /= (i + 1);
	}

	return result;
}

double CdlpFilter::factorial(unsigned n)
{
	return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}
