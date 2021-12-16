#include "BobyqaOptimizer.h"

#include <omp.h>
#include "dlib/optimization.h"
#include "dlib/threads.h"
#include <iostream> //eclude again
//#include <string>
#include <fstream>
#include <chrono>
#include <ctime>
#include <numeric>



void BobyqaOptimizer::optimize(OptimizationProblem& op, OptimizerOptions optOpt)
{
	unsigned number_Targets = op.getPitchTargets().size();
	ParameterSet ps = op.getParameters();
	bool optimizeBoundaries = ps.searchSpaceParameters.optimizeBoundaries;
	BoundaryVector initialBoundaries = op.getBoundaries();
	initialBoundaries.back() = op.getOriginalF0_Offset();
	BoundaryVector tmpBoundaries = initialBoundaries;
	BoundaryVector optBoundaries;
	TargetVector   tmpTargets = op.getPitchTargets();

	TimeSignal original_contour = op.getOriginalF0();

	std::vector<double> contour_values;
	for( const auto& sample: original_contour ) {
    	contour_values.push_back( sample.value );
	}

	double contour_onset_value = contour_values.front();
	double contour_std = get_standard_deviation( contour_values );
	double onset_min_bound = contour_onset_value - contour_std;
	double onset_max_bound = contour_onset_value + contour_std;

	int number_optVar = ps.searchSpaceParameters.numberOptVar; // Optimize the 3 target parameters by default

	double tmpMSE = 1e6;
	double tmpSCC = 0;


	const unsigned RANDOMITERATIONS = optOpt.maxIterations;
	const int patience = optOpt.patience;
	double epsilon = optOpt.epsilon;
	bool useEarlyStopping = optOpt.useEarlyStopping;

	const long max_f_evals = optOpt.maxCostEvaluations;
	const double rho_end = optOpt.rhoEnd;


	bool relativeDelta = true;




	// precalculate search space bounds (ssp_bounds)
	std::vector<double> min_bounds;
	std::vector<double> max_bounds;
	min_bounds.push_back(ps.searchSpaceParameters.meanSlope - ps.searchSpaceParameters.deltaSlope); //mmin
	min_bounds.push_back(ps.searchSpaceParameters.meanOffset - ps.searchSpaceParameters.deltaOffset); //bmin
	min_bounds.push_back(ps.searchSpaceParameters.meanTau - ps.searchSpaceParameters.deltaTau); //tmin
	max_bounds.push_back(ps.searchSpaceParameters.meanSlope + ps.searchSpaceParameters.deltaSlope); // mmax
	max_bounds.push_back(ps.searchSpaceParameters.meanOffset + ps.searchSpaceParameters.deltaOffset); //bmax
	max_bounds.push_back(ps.searchSpaceParameters.meanTau + ps.searchSpaceParameters.deltaTau); //tmax

	if (optimizeBoundaries)
	{
		min_bounds.push_back(-1 * ps.searchSpaceParameters.deltaBoundary); //boundary_min
		max_bounds.push_back(ps.searchSpaceParameters.deltaBoundary); //boundary_max
	}

	DlibVector lowerBound, upperBound;
	lowerBound.set_size(number_Targets * number_optVar + 1); // + 1 because of onset optimization
	upperBound.set_size(number_Targets * number_optVar + 1); // + 1 because of onset optimization

	lowerBound( 0 ) = onset_min_bound;
	upperBound( 0 ) = onset_max_bound;

	for (unsigned i = 0; i < number_Targets; ++i)
	{
		for (unsigned ssp_bound = 0; ssp_bound < number_optVar; ++ssp_bound)
		{
			lowerBound(number_optVar * i + ssp_bound + 1) = min_bounds.at(ssp_bound); // + 1 because of onset optimization
			upperBound(number_optVar * i + ssp_bound + 1) = max_bounds.at(ssp_bound); // + 1 because of onset optimization
		}
	}


	double trustRegion = max_bounds.at(0) - min_bounds.at(0);

	for (unsigned i = 1; i < number_optVar; ++i)
	{
		trustRegion = std::min(trustRegion, max_bounds.at(i) - min_bounds.at(i));
	}

	// optmization setup
	long npt(2 * lowerBound.size() + 1);	// number of interpolation points
	const double rho_begin = (trustRegion - 0.9) / 2.0; // initial trust region radius


	// initialize
	double fmin(1e6);
	DlibVector xtmp;
	dlib::mutex mu;

	// set up OpenMP
	int numThreads = omp_get_max_threads();
	omp_set_num_threads(numThreads);
	
	bool SearchFinished = false;
	int boundaryResetCounter = 0;
	int iteration = 0;
	int convergence = 0;

	std::vector<double> ftmp_vector(RANDOMITERATIONS, -1.0);
	

	//std::ofstream LOG;


	std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

	#pragma omp parallel for schedule(dynamic)

	for (int it = 0; it < RANDOMITERATIONS; ++it)
	{
		double ftmp;
		if (!SearchFinished)
		{
			DlibVector x;
			#pragma omp critical (getRandomValues)
			{
				std::cout << "random iteration: " << it << " begins." << std::endl;
				// random initialization
				srand(it);
				x.set_size(number_Targets * number_optVar + 1); // + 1 because of onset optimization
				x( 0 ) = getRandomValue( onset_min_bound, onset_max_bound );
				for (unsigned i = 0; i < number_Targets; ++i)
				{
					for (unsigned ssp_bound = 0; ssp_bound < number_optVar; ++ssp_bound)
					{
						x(number_optVar * i + ssp_bound + 1) = getRandomValue(min_bounds.at(ssp_bound), max_bounds.at(ssp_bound)); // + 1 because of onset optimization
					}
				}
			}
			try
			{
				#pragma omp critical (printProblem)
				{
					// DEBUG message
#ifdef DEBUG_MSG
					std::cout << "Random Iteration: " << it << "\n"
						<< "Opt.Problem:\n" << op
						//<< "X:\n" << x
						<< "Length of X: " << x.size() << "\n"
						<< "NPT: " << npt << "\n"
						<< "Lower Bound:\n" << lowerBound
						<< "upperBound:\n" << upperBound
						<< "rho Begin: " << rho_begin << "\n"
						<< "rho End: " << rho_end << "\n"
						<< "Max Evals: " << max_f_evals << std::endl;
					std::cout << "FTMP-Vector:" << std::endl;
					for (auto item:ftmp_vector)
					{
						std::cout << item << ", ";
					}
					std::cout << "\n" << std::endl;
#endif
				}
				// optimization algorithm: BOBYQA
				ftmp = dlib::find_min_bobyqa(op, x, npt, lowerBound, upperBound, rho_begin, rho_end, max_f_evals);
			}
			catch (dlib::bobyqa_failure & err)
			{
				// DEBUG message
				#ifdef DEBUG_MSG
				std::cout << "\t[optimize] WARNING: no convergence during optimization in iteration: " << it << std::endl << err.info << std::endl;
				#endif
			}
			// write optimization results back
			#pragma omp critical (updateMinValue)
			{
				ftmp_vector.at(it) = ftmp;
				#ifdef DEBUG_MSG
				std::cout << "Iteration nr: " << it << " fmin: " << fmin << " ftmp: " << ftmp << std::endl;
				#endif
				if (ftmp < fmin && ftmp > 0.0)	// opt returns 0 on error
				{
					if (useEarlyStopping)
					{ 
						if (ftmp < fmin - fmin * epsilon)
						{
							convergence = 0;
							fmin = ftmp;
							xtmp = x;
						}
					}
					else
					{
						fmin = ftmp;
						xtmp = x;
					}
					//std::cout << "ftmp: " << ftmp << std::endl;
				}
				++convergence;

				if (useEarlyStopping && convergence >= patience)
				{
					SearchFinished = true;
					std::cout << "" << std::endl;
					std::cout << "Search stopped early!" << std::endl;
				}

				++iteration;
#ifdef USE_WXWIDGETS
				if (waitbar_ != nullptr)
				{
					waitbar_->Update(iteration);
				}
#endif
				auto end = std::chrono::system_clock::now();
				std::time_t end_time = std::chrono::system_clock::to_time_t(end);
				std::cout << "random iteration: " << it << " ended. System time: " << std::ctime(&end_time) << std::endl;
			}
		}
	}





	if (fmin == 1e6)
	{
		std::cout << "  dlib error!!!!!!!!!!!!! " << std::endl;
		throw dlib::error("[optimize] BOBYQA algorithms didn't converge! Try to increase number of evaluations"); //reason can also be that non-valied numbers are passed as input, e.g. nan
	}

	// convert result to TargetVector
	TargetVector optTargets;

	if (!optimizeBoundaries)
	{
		optBoundaries = initialBoundaries;
	}
	else
	{

		//tmpBoundaries.back() = op.getOriginalF0_Offset();
		for (unsigned i = 0; i < number_Targets; ++i)
		{
			if (relativeDelta) {
				if (xtmp(number_optVar * i + 3 + 1) >= 0) // + 1 because of onset optimization
				{
					tmpBoundaries.at(i) += xtmp(number_optVar * i + 3 + 1) * (initialBoundaries.at(i + 1) - initialBoundaries.at(i)) * 0.01; // + 1 because of onset optimization
				}
				else
				{
					if (i == 0)
					{
						tmpBoundaries.at(i) += xtmp(number_optVar * i + 3 + 1) * 0.1 * 0.01; // + 1 because of onset optimization
					}
					else
					{
						tmpBoundaries.at(i) += xtmp(number_optVar * i + 3 + 1) * (initialBoundaries.at(i) - initialBoundaries.at(i - 1)) * 0.01; // + 1 because of onset optimization
					}
				}
			}
			else // for absolute delta (currently disabled)
			{
				tmpBoundaries.at(i) += xtmp(number_optVar * i + 3 + 1) / 1000; //divide by 1000 because delta is ms // + 1 because of onset optimization
			}

			if ((i == 0) && (tmpBoundaries.at(0) > op.getOriginalF0_Onset()))
			{
				tmpBoundaries.at(0) = op.getOriginalF0_Onset();
			}
		}
		std::sort(tmpBoundaries.begin(), tmpBoundaries.end());
		tmpBoundaries.back() = op.getOriginalF0_Offset();
		op.setBoundaries(tmpBoundaries);
		optBoundaries = tmpBoundaries;
	}
	for (unsigned i = 0; i < number_Targets; ++i)
	{
		PitchTarget pt;
		pt.slope = xtmp(number_optVar * i + 0 + 1); // + 1 because of onset optimization
		pt.offset = xtmp(number_optVar * i + 1 + 1); // + 1 because of onset optimization
		pt.tau = xtmp(number_optVar * i + 2 + 1); // + 1 because of onset optimization
		pt.duration = optBoundaries.at(i + 1) - optBoundaries.at(i);
		optTargets.push_back(pt);
	}
	double opt_onset_value = xtmp( 0 );


	// store optimum
	std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
	double computationTime = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count();
	op.setOptimum(optBoundaries, optTargets, opt_onset_value, computationTime, ftmp_vector);
	std::cout << "Elapsed time = " << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << "[s]" << std::endl;


		// DEBUG message
//#ifdef DEBUG_MSG
	std::cout << "\t[optimize] mse = " << fmin << std::endl;
//#endif
}

/*#ifdef USE_WXWIDGETS
void BobyqaOptimizer::optimize(OptimizationProblem& op, OptimizerOptions optOpt, wxGenericProgressDialog* waitbar)
{
	waitbar_ = waitbar;
	return optimize(op, optOpt);
}
#endif*/


double BobyqaOptimizer::getRandomValue(const double min, const double max)
{
	return min + ((double)rand() / RAND_MAX) * (max - min);
}

double BobyqaOptimizer::get_standard_deviation( const std::vector<double> v )
{
	double sum = std::accumulate(std::begin(v), std::end(v), 0.0);
	double m =  sum / v.size();
	double accum = 0.0;
	std::for_each (std::begin(v), std::end(v), [&](const double d) {
    	accum += (d - m) * (d - m);
	});
	double stdev = sqrt(accum / (v.size()-1));
	return stdev;
}
