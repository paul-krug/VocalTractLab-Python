#include <iostream>
#include <string>
#include <vector>

#include "BobyqaOptimizer.h"
#include "TargetOptimizerApi.h"
//#include "TextGridReader.h"

//#include <fstream>




FitData estimate_targets(	//int init_bounds,
						//double *arr_times,
						//double *arr_values,
						//double *arr_boundaries,
						std::vector<double> input_times,
						std::vector<double> input_values,
						std::vector<double> boundaries,

						//int init_bounds,

						double weight_slope,
						double weight_offset,
						double weight_tau,
						double weight_lambda,

						double delta_slope,
						double delta_offset,
						double delta_tau,
						double delta_boundary,

						double mean_slope,
						//double mean_offset,
						double mean_tau,

						int max_iterations,
						int max_cost_evaluations,
						double rho_end,

						bool use_early_stopping,
						double epsilon,
						int patience
					)
{
	//std::vector<double> input_times( arr_times, arr_times + sizeof arr_times / sizeof arr_times[0] );
	//std::vector<double> input_values( arr_values, arr_values + sizeof arr_values / sizeof arr_values[0] );
	//std::vector<double> boundaries;
	//for ( const auto & boundary : arr_boundaries )
	//{
	//	boundaries.push_back( boundary );
	//}


	//( arr_boundaries, arr_boundaries + sizeof arr_boundaries / sizeof arr_boundaries[0] );
	if ( input_times.size() != input_values.size() )
	{
		std::cout << "Input time and value vectors need to have the same size!" << std::endl;
		//return;
	}
	std::cout << "weight_slope:" << weight_slope << std::endl;
	std::cout << "weight_offset:" << weight_offset << std::endl;
	std::cout << "use_early_stopping:" << use_early_stopping << std::endl;
	std::cout << "max_iterations:" << max_iterations << std::endl;


	std::cout << "boundaries size:" << boundaries.size() << std::endl;
	std::cout << "times size:" << input_times.size() << std::endl;
	std::cout << "values size:" << input_values.size() << std::endl;

	TimeSignal f0;
	Sample sample;
	std::cout << "Debug message 1" << std::endl;
	for ( unsigned i = 0; i < input_times.size(); ++i ) {
		sample.time = input_times.at( i );
		sample.value = input_values.at( i );
		f0.push_back( sample );
	}
	std::cout << "time signal size:" << f0.size() << std::endl;
	std::cout << "Debug message 2" << std::endl;
	//calculate mean f0
	double meanF0 = 0.0;
	for ( std::vector<Sample>::const_iterator sp = f0.begin(); sp != f0.end(); ++sp )
	{
		meanF0 += sp->value;
	}
	meanF0 /= f0.size();

	ParameterSet parameters;
	std::cout << "Debug message 3" << std::endl;
	parameters.regularizationParameters.weightSlope  = weight_slope;
	parameters.regularizationParameters.weightOffset = weight_offset;
	parameters.regularizationParameters.weightTau    = weight_tau;
	parameters.regularizationParameters.lambda       = weight_lambda;

	parameters.searchSpaceParameters.deltaSlope    = delta_slope;
	parameters.searchSpaceParameters.deltaOffset   = delta_offset;
	parameters.searchSpaceParameters.deltaTau      = delta_tau;
	parameters.searchSpaceParameters.deltaBoundary = delta_boundary;

	parameters.searchSpaceParameters.meanSlope  = mean_slope;
	parameters.searchSpaceParameters.meanOffset = meanF0;
	parameters.searchSpaceParameters.meanTau    = mean_tau;
			
	parameters.searchSpaceParameters.optimizeBoundaries = ( parameters.searchSpaceParameters.deltaBoundary != 0 );
	parameters.searchSpaceParameters.numberOptVar       = ( parameters.searchSpaceParameters.optimizeBoundaries ? 4 : 3 );

	//parameters.searchSpaceParameters.initBounds = init_bounds;
	//std::cout << "Debug message 4" << std::endl;
	//if ( parameters.searchSpaceParameters.initBounds < 2 )
	//{
	//	std::cout << "You need at least two boundaries! Specify a TextGrid file or set initBounds >=2." << std::endl;
	//	return;
	//}
	//std::cout << "Debug message 5" << std::endl;
	BoundaryVector bounds = boundaries;

	OptimizerOptions optOpt;
	std::cout << "Debug message 6" << std::endl;
	optOpt.maxIterations = max_iterations;
	optOpt.maxCostEvaluations = max_cost_evaluations;
	optOpt.rhoEnd = rho_end;
	optOpt.useEarlyStopping = use_early_stopping;
	optOpt.epsilon = epsilon;
	std::cout << "Debug message 7" << std::endl;
	if ( patience < 0 )
	{
		optOpt.patience = OptimizerOptions::recommendedPatience(bounds.size() - 1);
	} else {
		optOpt.patience = patience;
	}

	//if (parameters.searchSpaceParameters.initBounds != 0)
	//{
	//	std::cout << "bounds init" << std::endl;
	//	double pitch_start = f0.at(0).time;
	//	double pitch_end = f0.back().time;
	//	double pitch_interval = pitch_end - pitch_start;
	//	double step = pitch_interval / (parameters.searchSpaceParameters.initBounds - 1);
	//	std::vector<double> initBoundaries;
	//	for (int i = 0; i < parameters.searchSpaceParameters.initBounds; ++i)
	//	{
	//		initBoundaries.push_back(pitch_start + i * step);
	//	}
	//	bounds = initBoundaries;
	//	initBoundaries.clear();
	//}

	// main functionality
	std::cout << "bounds size:" << bounds.size() << std::endl;
	std::cout << "Debug message 8" << std::endl;
	for(int i = 0; i < bounds.size(); ++i)
	{
		std::cout << "bound entry nr:" << i << " is: " << bounds.at(i) << std::endl;
	}
	OptimizationProblem problem( parameters, f0, bounds );
	BobyqaOptimizer optimizer;
	std::cout << "Debug message 9" << std::endl;
	optimizer.optimize( problem, optOpt );
	TargetVector optTargets = problem.getPitchTargets();
	BoundaryVector optBoundaries = problem.getBoundaries();
	std::cout << "Debug message 10" << std::endl;
	TimeSignal optF0 = problem.getModelF0();
	Sample onset = problem.getOnset();
	std::vector<double> ftmp_vector = problem.getOptimizationSolutions();
	std::cout << "Debug message 11" << std::endl;

	FitData fit_results;
	fit_results.res_targets = problem.getPitchTargets();
	fit_results.res_boundaries = problem.getBoundaries();
	fit_results.res_trajectory = problem.getModelF0();
	fit_results.res_ftmp = problem.getOptimizationSolutions();
	fit_results.res_fmin = problem.getCostFunction();
	fit_results.res_rmse = problem.getRootMeanSquareError();
	fit_results.res_corr = problem.getCorrelationCoefficient();
	fit_results.res_time = problem.getComputationTime();
	fit_results.res_onset= problem.getOnset().value;

	return fit_results;
}