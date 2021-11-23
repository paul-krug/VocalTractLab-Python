# CMakeList.txt : CMake project for VocalTractLabApi, include source and define
# project specific logic here.
#
cmake_minimum_required (VERSION 3.8)

project ("VocalTractLabApi")

add_definitions(-D_USE_MATH_DEFINES -D_CRT_SECURE_NO_WARNINGS)
add_library (VocalTractLabApi SHARED VocalTractLabApi.def AnatomyParams.cpp AnatomyParams.h AudioFile.h Constants.h Dsp.cpp Dsp.h F0EstimatorYin.cpp F0EstimatorYin.h GeometricGlottis.cpp GeometricGlottis.h Geometry.cpp Geometry.h GesturalScore.cpp GesturalScore.h Glottis.cpp Glottis.h IirFilter.cpp IirFilter.h ImpulseExcitation.cpp ImpulseExcitation.h LfPulse.cpp LfPulse.h Matrix2x2.cpp Matrix2x2.h PoleZeroPlan.cpp PoleZeroPlan.h Sampa.cpp Sampa.h SegmentSequence.cpp SegmentSequence.h Signal.cpp Signal.h   Splines.cpp Splines.h StaticPhone.cpp StaticPhone.h Surface.cpp Surface.h Synthesizer.cpp Synthesizer.h TdsModel.cpp TdsModel.h TimeFunction.cpp TimeFunction.h TlModel.cpp TlModel.h TriangularGlottis.cpp TriangularGlottis.h Tube.cpp Tube.h TubeSequence.h TwoMassModel.cpp TwoMassModel.h VocalTract.cpp VocalTract.h VocalTractLabApi.cpp VocalTractLabApi.h VoiceQualityEstimator.cpp VoiceQualityEstimator.h VowelLf.cpp VowelLf.h XmlHelper.cpp XmlHelper.h XmlNode.cpp XmlNode.h)
