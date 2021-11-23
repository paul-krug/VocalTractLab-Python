# CMakeList.txt : CMake project for TargetOptimizerApi, include source and define
# project specific logic here.
#
cmake_minimum_required (VERSION 3.9)

project ("TargetOptimizerApi")

find_package(OpenMP REQUIRED)

add_definitions( -D_USE_MATH_DEFINES -D_CRT_SECURE_NO_WARNINGS )
add_library( TargetOptimizerApi SHARED TargetOptimizerApi.def BobyqaOptimizer.cpp BobyqaOptimizer.h CdlpFilter.cpp CdlpFilter.h Data.h ./dlib/all/source.cpp dlib/error.h dlib/matrix.h OptimizationProblem.cpp OptimizationProblem.h TamModelF0.cpp TamModelF0.h TargetOptimizerApi.cpp TargetOptimizerApi.h )

target_link_libraries( TargetOptimizerApi PRIVATE OpenMP::OpenMP_CXX )
