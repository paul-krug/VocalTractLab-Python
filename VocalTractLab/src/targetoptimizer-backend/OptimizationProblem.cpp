#include "OptimizationProblem.h"
#include <iostream>


OptimizationProblem::OptimizationProblem(
	const ParameterSet& parameters,
	const TimeSignal& originalF0,
	const BoundaryVector& bounds
	):
m_computationTime(-1.0),
m_parameters(parameters), 
m_originalF0(originalF0), 
m_bounds(bounds), 
m_modelOptimalF0(bounds, 
originalF0[0].value
) {};

void OptimizationProblem::setOptimum(
	const BoundaryVector& boundaries,
	const TargetVector& targets,
	const double onset_value,
    const double computationTime,
    const std::vector<double> optimizationSolutions )
{
	//m_modelOptimalF0.setOnset( m_originalF0[0].time * targets[0].slope + targets[0].offset ); // added
	m_modelOptimalF0.setOnset( onset_value );
	m_modelOptimalF0.setBoundaries( boundaries );
	m_modelOptimalF0.setPitchTargets(targets);
	m_computationTime = computationTime;
	m_optimizationSolutions = optimizationSolutions;
}

double OptimizationProblem::getComputationTime() const
{
	return m_computationTime;
}

ParameterSet OptimizationProblem::getParameters() const
{
	return m_parameters;
}

TimeSignal OptimizationProblem::getModelF0() const
{
	double samplingfrequency = SAMPLERATE; // Hz
	double dt = 1.0 / samplingfrequency;
	return m_modelOptimalF0.calculateF0(dt);
}

TargetVector OptimizationProblem::getPitchTargets() const
{
	return m_modelOptimalF0.getPitchTargets();
}

BoundaryVector OptimizationProblem::getBoundaries() const
{
	return m_modelOptimalF0.getBoundaries();
}

void OptimizationProblem::setBoundaries( const BoundaryVector& boundaries )
{
	m_bounds = boundaries;
}

Sample OptimizationProblem::getOnset() const
{
	return m_modelOptimalF0.getOnset();
}

TimeSignal OptimizationProblem::getOriginalF0() const
{
	return m_originalF0;
}

std::vector<double> OptimizationProblem::getOptimizationSolutions() const
{
	return m_optimizationSolutions;
}

double OptimizationProblem::getOriginalF0_Onset() const
{
	return m_originalF0.front().time;
}

double OptimizationProblem::getOriginalF0_Offset() const
{
	return m_originalF0.back().time;
}



SampleTimes OptimizationProblem::extractTimes(const TimeSignal& f0)
{
	SampleTimes times;
	for (unsigned i = 0; i < f0.size(); ++i)
	{
		times.push_back(f0[i].time);
	}

	return times;
}

DlibVector OptimizationProblem::signal2dlibVec(const TimeSignal& f0)
{
	DlibVector values;
	values.set_size(f0.size());
	for (unsigned i = 0; i < f0.size(); ++i)
	{
		values(i) = f0[i].value;
	}

	return values;
}

double OptimizationProblem::getCorrelationCoefficient() const
{
	SampleTimes times = extractTimes(m_originalF0);
	TimeSignal modelF0 = m_modelOptimalF0.calculateF0(times);

	// Dlib format
	DlibVector orig = signal2dlibVec(m_originalF0);
	DlibVector model = signal2dlibVec(modelF0);

	// return correlation between filtered and original f0
	DlibVector x = orig - dlib::mean(orig);
	DlibVector y = model - dlib::mean(model);
	return (dlib::dot(x, y)) / ((std::sqrt(dlib::sum(dlib::squared(x)))) * (std::sqrt(dlib::sum(dlib::squared(y)))));

}

double OptimizationProblem::getSquareCorrelationCoefficient( const TamModelF0& tamF0 ) const
{
	SampleTimes times = extractTimes(m_originalF0);
	TimeSignal modelF0 = tamF0.calculateF0(times);

	// Dlib format
	DlibVector orig = signal2dlibVec(m_originalF0);
	DlibVector model = signal2dlibVec(modelF0);

	// return correlation between filtered and original f0
	DlibVector x = orig - dlib::mean(orig);
	DlibVector y = model - dlib::mean(model);
	return dlib::dot(x, y) * dlib::dot(x, y) / ( dlib::sum(dlib::squared(x)) * dlib::sum(dlib::squared(y)) );

}

double OptimizationProblem::getCostFunction() const
{
	return costFunction(m_modelOptimalF0);
}

double OptimizationProblem::getRootMeanSquareError() const
{
	SampleTimes times = extractTimes(m_originalF0);
	TimeSignal modelF0 = m_modelOptimalF0.calculateF0(times);

	// Dlib format
	DlibVector orig = signal2dlibVec(m_originalF0);
	DlibVector model = signal2dlibVec(modelF0);

	// return RMSE between filtered and original f0
	return std::sqrt(dlib::mean(dlib::squared(model - orig)));
}

double OptimizationProblem::getMeanSquareError( const TamModelF0& tamF0 ) const
{
	SampleTimes times = extractTimes(m_originalF0);
	TimeSignal modelF0 = tamF0.calculateF0(times);

	// Dlib format
	DlibVector orig = signal2dlibVec(m_originalF0);
	DlibVector model = signal2dlibVec(modelF0);

	// return MSE between filtered and original f0
	return dlib::mean(dlib::squared(model - orig));
}

std::tuple<double, double> OptimizationProblem::getOptStats( const BoundaryVector& boundaries, const TargetVector& targets ) const
{
	TamModelF0 tamF0(boundaries, m_originalF0[0].value);
	//TamModelF0 tamF0(boundaries, m_originalF0[0].time * targets[0].slope + targets[0].offset);
	tamF0.setPitchTargets(targets);

	double MSE = getMeanSquareError( tamF0 );
	double SSC = getSquareCorrelationCoefficient( tamF0 );

	return std::make_tuple(MSE, SSC);
}

double OptimizationProblem::operator() (const DlibVector& arg) const
{
	// convert data
	TargetVector targets;
	BoundaryVector initialBounds = m_bounds;
	initialBounds.back() = getOriginalF0_Offset();
	BoundaryVector boundaries = initialBounds;
	bool relativeDelta = true;



	double modified_Bound = 0.0;
	int number_Targets = int( ( arg.size() - 1 ) / m_parameters.searchSpaceParameters.numberOptVar ); // - 1 because of onset optimization
	int number_optVar = int( m_parameters.searchSpaceParameters.numberOptVar );

	if ( m_parameters.searchSpaceParameters.optimizeBoundaries == true )
	{
		for (unsigned i = 0; i < number_Targets; ++i)
		{
			if ( relativeDelta ){
				if ( arg(number_optVar * i + 3 + 1) >= 0) // + 1 because of onset optimization
				{
					boundaries.at(i) += arg(number_optVar * i + 3 + 1)*( initialBounds.at(i+1)-initialBounds.at(i) )*0.01; // + 1 because of onset optimization
				}
				else
				{
					if ( i==0 )
					{
						boundaries.at(i) += arg(number_optVar * i + 3 + 1) * 0.1 * 0.01; // + 1 because of onset optimization
					}
					else
					{
						boundaries.at(i) += arg(number_optVar * i + 3 + 1)*( initialBounds.at(i)-initialBounds.at(i-1) ) *0.01; // + 1 because of onset optimization
					}
				}
			}
			else
			{
				boundaries.at(i) += arg(m_parameters.searchSpaceParameters.numberOptVar * i + 3 + 1)/1000; // + 1 because of onset optimization
			}
			if ( (i==0) && (boundaries.at(0) > m_originalF0[0].time))
			{
				boundaries.at(0) = m_originalF0[0].time;
			}
		}
		std::sort( boundaries.begin(), boundaries.end() );
		boundaries.back()  = getOriginalF0_Offset();
	}

	for (unsigned i = 0; i < number_Targets; ++i)
	{
		PitchTarget pt;
		pt.slope = arg(m_parameters.searchSpaceParameters.numberOptVar * i + 0 + 1); // + 1 because of onset optimization
		pt.offset = arg(m_parameters.searchSpaceParameters.numberOptVar * i + 1 + 1); // + 1 because of onset optimization
		pt.tau = arg(m_parameters.searchSpaceParameters.numberOptVar * i + 2 + 1); // + 1 because of onset optimization
		pt.duration = ( boundaries.at(i+1) - boundaries.at(i) );

		targets.push_back(pt);
	}

	// DEBUG #Hack
	if (targets.back().duration == 0.0)
	{
		targets.back().duration += 0.001;
		targets.end()[-2].duration -= 0.001;
		boundaries.end()[-2] -= 0.001;
	}
	// END DEBUG
	TamModelF0 tamF0(boundaries, arg( 0 ) );
	//TamModelF0 tamF0(boundaries, m_originalF0[0].value);
	//TamModelF0 tamF0(boundaries, m_originalF0[0].time * targets[0].slope + targets[0].offset);
	tamF0.setPitchTargets(targets);


	return costFunction(tamF0);
}

double OptimizationProblem::costFunction(const TamModelF0& tamF0) const // TODO: pow(x,2) is much slower than x*x (depending on compiler opt), Change it
{
	// get model f0
	SampleTimes times = extractTimes(m_originalF0);
	TimeSignal modelF0 = tamF0.calculateF0(times);

	// calculate error
	double error = 0.0;

	for (unsigned i = 0; i < modelF0.size(); ++i)
	{
		error += std::pow((m_originalF0[i].value - modelF0[i].value), 2.0);
	}

	// calculate penalty term
	double penalty = 0.0;
	TargetVector targets = tamF0.getPitchTargets();
	for (unsigned i = 0; i < targets.size(); ++i)

	{
		penalty += (m_parameters.regularizationParameters.weightSlope * std::pow(targets[i].slope - m_parameters.searchSpaceParameters.meanSlope, 2.0));
		penalty += (m_parameters.regularizationParameters.weightOffset * std::pow(targets[i].offset - m_parameters.searchSpaceParameters.meanOffset, 2.0));
		penalty += (m_parameters.regularizationParameters.weightTau * std::pow(targets[i].tau - m_parameters.searchSpaceParameters.meanTau, 2.0));
	}

	return error + m_parameters.regularizationParameters.lambda * penalty;
}

std::ostream& operator<<(std::ostream& os, OptimizationProblem& op)
{
	using namespace std;
	auto tmp = op.getParameters();
	os << tmp << endl;
	os << "Orig. F0: " << endl;
	for (auto f0val:op.getOriginalF0())
	{
		os << f0val.time << "\t" << f0val.value << endl;
	}
	os << "Boundaries: " << endl;
	for (auto bound : op.getBoundaries())
	{
		os << bound << endl;
	}
	return os;
}
