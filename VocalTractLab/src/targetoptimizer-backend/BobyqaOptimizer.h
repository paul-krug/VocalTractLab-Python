#pragma once
#include "OptimizationProblem.h"
#include <string>



struct OptimizerOptions
{
	int maxIterations{ 1000 };
	bool useEarlyStopping{ false };
	double epsilon{ 0.1 };  // Minimum relative improvement to update the optimum if early stopping is used
	int patience{ recommendedPatience(1) };  // If after this many iterations the minimum has not improved, the optimization stops. 
	int maxCostEvaluations{ 100000 };
	double rhoEnd{ 0.001 };
	static int recommendedPatience(int numTargets) { return 15 * numTargets; }
};

// solver for an optimization problem utilizing BOBYQA algorithm
class BobyqaOptimizer {
public:
	// constructors
	BobyqaOptimizer() = default;

	// public member functions
	void optimize(OptimizationProblem& op, OptimizerOptions optOpt);


private:
	// private member functions
	static double getRandomValue(const double min, const double max);
};