#pragma once
#include "Data.h"
 
//more about this in reference 1
#ifdef DLLDIR_EX
   #define DLLDIR __declspec(dllexport)   // export DLL information
 
#else
   #define DLLDIR  __declspec(dllimport)   // import DLL information
 
#endif 

struct FitData {
	TargetVector res_targets;
	BoundaryVector res_boundaries;
	TimeSignal res_trajectory;
	std::vector<double> res_ftmp;
	double res_fmin;
	double res_rmse;
	double res_corr;
	double res_time;
	double res_onset;

	//double res_f0_time_array;
	//double res_f0_value_array;
	//double res_boundary_array;
	//double res_t_slope_array;
	//double res_t_offset_array;
	//double res_t_tau_array;
	//double res_t_duration_arra;
};

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
						int patience);