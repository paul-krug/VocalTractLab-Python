#pragma once
#include "dlib/matrix.h"
#include "Data.h"
#include "TamModelF0.h"

#include <tuple>


struct RegularizationParameters
{
	double lambda{ 0.01 };
	double weightSlope{ 1.0 };
	double weightOffset{ 0.5 };
	double weightTau{ 0.1 };
};

struct SearchSpaceParameters
{
	double deltaSlope{ 50.0 };
	double deltaOffset{ 20.0 };
	double deltaTau{ 5.0 };
	double deltaBoundary{ 100.0 };
	int    initBounds{ 0 };
	bool   optimizeBoundaries{ deltaBoundary > 0 };
	int    numberOptVar{ optimizeBoundaries ? 4 : 3 };
	double meanSlope{ 0.0 };
	double meanOffset{ 0.0 };
	double meanTau{ 20.0 };
};

// parameter set defining an optimisation problem
struct ParameterSet
{
	RegularizationParameters regularizationParameters;
	SearchSpaceParameters searchSpaceParameters;

	friend std::ostream& operator<<(std::ostream& os, ParameterSet& ps)
	{
		using namespace std;
		return		os << "Reg. Parameters:" << "\n"
			<< ps.regularizationParameters.lambda << "\n"
			<< ps.regularizationParameters.weightSlope << "\n"
			<< ps.regularizationParameters.weightOffset << "\n"
			<< ps.regularizationParameters.weightTau << "\n"
			<< "Search Space Parameters:" << "\n"
			<< ps.searchSpaceParameters.deltaSlope << "\n"
			<< ps.searchSpaceParameters.deltaOffset << "\n"
			<< ps.searchSpaceParameters.deltaTau << "\n"
			<< ps.searchSpaceParameters.deltaBoundary << "\n"
			<< ps.searchSpaceParameters.optimizeBoundaries << "\n"
			<< ps.searchSpaceParameters.numberOptVar << "\n"
			<< ps.searchSpaceParameters.meanSlope << "\n"
			<< ps.searchSpaceParameters.meanOffset << "\n"
			<< ps.searchSpaceParameters.meanTau;
	}
};

// dlib linear algebra column vector for optimization tasks
typedef dlib::matrix<double, 0, 1> DlibVector;

// optimization problem for calculating pitch targets
class OptimizationProblem {
public:
	// constructors
	OptimizationProblem(const ParameterSet& parameters, const TimeSignal& originalF0, const BoundaryVector& bounds);

	// public member functions
	void setOptimum(const BoundaryVector& boundaries, const TargetVector& targets,
	                const double computationTime, const std::vector<double> optimizationSolutions);
	void setBoundaries(const BoundaryVector& boundaries);

	ParameterSet getParameters() const;
	TimeSignal getModelF0() const;
	TargetVector getPitchTargets() const;
	BoundaryVector getBoundaries() const;
	Sample getOnset() const;
	TimeSignal getOriginalF0() const;
	double getOriginalF0_Onset() const;
	double getOriginalF0_Offset() const;
	double getComputationTime() const;
	double getCorrelationCoefficient() const;
	double getCostFunction() const;
	double getRootMeanSquareError() const;
	double getSquareCorrelationCoefficient( const TamModelF0& tamF0 ) const;
	double getMeanSquareError( const TamModelF0& tamF0 ) const;
	std::tuple<double, double> getOptStats( const BoundaryVector& boundaries, const TargetVector& targets ) const;
	std::vector<double> getOptimizationSolutions() const;

	// operator called by optimizer
	double operator() (const DlibVector& arg) const;
	friend std::ostream& operator<<(std::ostream& os, OptimizationProblem& op);


private:
	// private member functions
	double costFunction(const TamModelF0& tamF0) const;
	static SampleTimes extractTimes(const TimeSignal& f0);
	static DlibVector signal2dlibVec(const TimeSignal& f0);

	// data members
	ParameterSet m_parameters;
	TimeSignal m_originalF0;
	BoundaryVector m_bounds;
	double m_computationTime;
	std::vector<double> m_optimizationSolutions;

	// store result
	TamModelF0 m_modelOptimalF0;
};