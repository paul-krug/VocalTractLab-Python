cdef extern from "VocalTractLabApi.h":
	int vtlCalcTongueRootAutomatically(bool automaticCalculation);
	void vtlGetVersion(char *version);
	int vtlInitialize(const char *speakerFileName);
	int vtlClose();
	int vtlGetConstants(int *audioSamplingRate, int *numTubeSections, int *numVocalTractParams, int *numGlottisParams);