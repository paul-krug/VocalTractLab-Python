// ****************************************************************************
// This file is part of VocalTractLab.
// Copyright (C) 2020, Peter Birkholz, Dresden, Germany
// The new solver (Cholesky factorization) was added by Johann Marwitz.
// www.vocaltractlab.de
// author: Peter Birkholz
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program. If not, see <http://www.gnu.org/licenses/>.
//
// ****************************************************************************

#ifndef __TDS_MODEL__
#define __TDS_MODEL__

#include <cmath>
#include <string>
#include <random>

#include "Dsp.h"
#include "Geometry.h"
#include "IirFilter.h"
#include "GeometricGlottis.h"
#include "Tube.h"
#include "Constants.h"
#include "TimeFunction.h"


// ****************************************************************************
/// Class for the simulation of vocal tract acoustics in the time domain on the
/// basis of a branched tube model of the vocal tract.
/// All simulation is calculated in CGS-units.
// ****************************************************************************

class TdsModel
{
public:

  // ************************************************************************
  // Constants.
  // ************************************************************************

  /// Two additional output currents at the mouth and nose opening
  static const int NUM_BRANCH_CURRENTS = Tube::NUM_SECTIONS + 4;

  // Max. number of non-zero places per row in the matrix (for Gauss-Seidel)
  static const int MAX_CONCERNED_MATRIX_COLUMNS = 16;
  static const int MAX_CONSTRICTIONS = 4;
  static const int CONSTRICTION_BUFFER_SIZE = 65536;    // For more than 1 s
  static const int CONSTRICTION_BUFFER_SIZE_MASK = 65535;

  static const int NUM_NOISE_BUFFER_SAMPLES = 8;
  static const int NOISE_BUFFER_MASK = 7;

  // Max. number of non-zero places per row in the matrix when it is saved in symmetric envelope structure (for cholesky factorization)
  static const int MAX_CONCERNED_MATRIX_COLUMNS_SYMMETRIC_ENVELOPE = 57;
  // Max. number of non-zero places per	column in the matrix when it is saved in symmetric envelope structure (for cholesky factorization)
  static const int MAX_CONCERNED_MATRIX_ROWS_SYMMETRIC_ENVELOPE = 10;

  // Max. number of non-zero places per row in the symmetric saved matrix
  static const int MAX_CONCERNED_MATRIX_COLUMNS_SYMMETRIC = 3;
  // Max. number of non-zero places per	column the symmetric saved matrix
  static const int MAX_CONCERNED_MATRIX_ROWS_SYMMETRIC = 4;

  // For the temporal discretization
  static const double THETA;
  static const double THETA1;

  static const double MIN_AREA_CM2;
  static const double NOISE_CUTOFF_FREQ;


  // ************************************************************************
  /// Some options
  // ************************************************************************

  enum SolverType
  {
    SOR_GAUSS_SEIDEL,
    CHOLESKY_FACTORIZATION,
    NUM_SOLVER_TYPES
  };

  struct Options
  {
    bool turbulenceLosses;        ///< Consider fluid dynamic losses due to turbulence
    bool softWalls;               ///< Consider losses due to soft walls
    bool generateNoiseSources;    ///< Generate noise sources during the simulation
    bool radiationFromSkin;       ///< Allow sound radiation from the skin
    bool piriformFossa;           ///< Include the piriform fossa
    bool innerLengthCorrections;  ///< Additional inductivities between adjacent sections
    bool transvelarCoupling;      ///< Sound transmission through the velum tissue?
    SolverType solverType;
  };

  // ************************************************************************
  /// Information for one individual lumped noise source
  // ************************************************************************

  struct NoiseSource
  {
    bool   isFirstOrder;
    double cutoffFreq;
    double targetAmp1kHz;   // Amp. for 1 kHz bandwidth of the noise shaping filter
    double currentAmp1kHz;
    // Input and output sample buffers for the spectral shaping filter
    double inputBuffer[NUM_NOISE_BUFFER_SAMPLES];
    double outputBuffer[NUM_NOISE_BUFFER_SAMPLES];
    double sample;        ///< The current sampling point of the noise source
  };

  // ************************************************************************
  /// Information about a supraglottal constriction
  // ************************************************************************

  struct Constriction
  {
    int firstSection;
    int lastSection;
    int narrowestSection;
    double obstaclePos;
    int obstacleSection;
    double area;          // in cm^2
    double flow;          // in cm^3 / s
    double velocity;      // in cm / s
    double cutoffFreq;   // in Hz
    double gain;
    double fullAmp;
    Tube::Articulator articulator;
  };


  // ************************************************************************
  /// Structure for one individual branch current in the electrical
  /// network. The identity of a branch current is defined by the
  /// indices of the tube sections from where it comes and where
  /// it goes.
  // ************************************************************************

  struct BranchCurrent
  {
    int sourceSection;
    int targetSection;

    double magnitude;
    double magnitudeRate;
    double noiseMagnitude;    ///< The flow low-pass filtered at 1000 Hz
  };

  // ************************************************************************
  /// An individual short homogeneous tube section.
  // ************************************************************************

  struct TubeSection
  {
//    bool isDynamic;        ///< Can the network components R, C, L, ... change ?

    double pos;
    double area;
    double length;
    double volume;
    Tube::Articulator articulator;

    NoiseSource monopoleSource; ///< Is created in the center of the tube section
    NoiseSource dipoleSource;   ///< Is created at the entrance of the tube section

    double pressure;
    double pressureRate;

    /// \name Indices of the inflowing and outflowing currents
    /// @{
    int currentIn;
    int currentOut[2];
    /// @}

    /// \name Wall properties
    /// @{
    double Mw;    ///< Mass per unit-area
    double Bw;    ///< Resistance per unit-area
    double Kw;    ///< Stiffness per unit-area
    /// @}

    /// \name Wall currents
    /// @{
    double wallCurrent;       ///< Current flow "into" the walls
    double wallCurrentRate;   ///< 1st derivative of the wall-flow
    double wallCurrentRate2;  ///< 2nd derivative of the wall-flow
    /// @}

    double L;        ///< Inductivity
    double C;        ///< Capacity
    double R[2];     ///< Ohm's resistance left and right
    double S;        ///< Pressure source at the inlet of the section (A constant in the pressure-difference eq.)

    // For the wall vibration
    double alpha;
    double beta;

    // Temporary values
    double D;
    double E;
  };

  // ************************************************************************
  // Public variables.
  // ************************************************************************

  double flowSourceAmp;
  int    flowSourceSection;   ///< Where is the periodic volume velocity source ?
  double pressureSourceAmp;
  int    pressureSourceSection;  ///< Where is the static pressure source ?

  NoiseSource lipsDipoleSource;   ///< Last noise source at the mouth opening

  TubeSection tubeSection[Tube::NUM_SECTIONS];
  BranchCurrent branchCurrent[NUM_BRANCH_CURRENTS];

  // Help variables to effectively solve the system of eqs. with Gauss-Seidel
  int numFilledRowValues[NUM_BRANCH_CURRENTS];
  int filledRowIndex[NUM_BRANCH_CURRENTS][MAX_CONCERNED_MATRIX_COLUMNS];

  // Help variables to effectively solve the system of eqs. with cholesky factorization
  int numFilledRowValuesSymmetricEnvelope[NUM_BRANCH_CURRENTS];
  int filledRowIndexSymmetricEnvelope[NUM_BRANCH_CURRENTS][MAX_CONCERNED_MATRIX_COLUMNS_SYMMETRIC_ENVELOPE];
  int numFilledColumnValuesSymmetricEnvelope[NUM_BRANCH_CURRENTS];
  int filledColumnIndexSymmetricEnvelope[NUM_BRANCH_CURRENTS][MAX_CONCERNED_MATRIX_ROWS_SYMMETRIC_ENVELOPE];

  bool doNetworkInitialization;
  double timeStep;
  /// Aspiration strength from -40 dB to 0 dB.
  double aspirationStrength_dB;

  double matrix[NUM_BRANCH_CURRENTS][NUM_BRANCH_CURRENTS];
  double factorizationMatrix[NUM_BRANCH_CURRENTS][NUM_BRANCH_CURRENTS];
  double solutionVector[NUM_BRANCH_CURRENTS];
  double flowVector[NUM_BRANCH_CURRENTS];

  int SORIterations;      ///< For external evaluation of the iterations needed to solve the system of eqs.

  IirFilter glottalToneFilter;
  IirFilter transglottalPressureFilter;
  double glottalBernoulliFactor;

  IirFilter transvelarCouplingFilter1;
  IirFilter transvelarCouplingFilter2;

  Constriction constriction[MAX_CONSTRICTIONS];
  int numConstrictions;
  
  // The constriction that is nearest to this tube section will be
  // monitored (constriction data are written to the buffer).
  int constrictionMonitorTubeSection;
  Constriction *constrictionBuffer;

  Options options;


  // ************************************************************************
  // Public methods of the class
  // ************************************************************************

public:
  TdsModel();
  ~TdsModel();

  void initModel();
  void resetMotion();
  bool saveConstrictionBuffer(std::string fileName);

  // **************************************************************
  /// \name These functions should be called for each time step
  // **************************************************************
  /// @{
  void setTube(Tube *tube, bool filtering = false);
  void getTube(Tube *tube);
  void setFlowSource(double flow_cm3_s, int section);
  void setPressureSource(double pressure_dPa, int section = Tube::FIRST_TRACHEA_SECTION);
  double proceedTimeStep(double &mouthFlow_cm3_s, double &nostrilFlow_cm3_s, 
    double &skinFlow_cm3_s, const string &matrixFileName = "");
  /// @}

  // **************************************************************

  void solveEquationsSor(const string &matrixFileName = "");
  void solveEquationsCholesky();
  int getSampleIndex() { return position; }
  void getSectionFlow(int sectionIndex, double &inflow, double &outflow);
  double getSectionPressure(int sectionIndex);

  // ************************************************************************
  // Private data.
  // ************************************************************************

private:
  int position;         ///< Internal counter for the sampling position
  double teethPosition; ///< Position of the teeth (from the glottis)
  // Elevation of the tongue tip side (corresponding to TS3 of the vocal tract model)
  double tongueTipSideElevation;
  std::mt19937 randomNumberGenerator;

  // ************************************************************************
  // Private functions.
  // ************************************************************************

private:
  void prepareTimeStep();

  double getJunctionInductance(double A1_cm2, double A2_cm2);

  void calcMatrix();
  void updateVariables();
  
  void resetConstriction(Constriction *c);
  void calcNoiseSources();
  void calcNoiseSample(NoiseSource* s, double ampThreshold);

  double getCurrentIn(const int section);
  double getCurrentOut(const int section);
  double getCurrentIn(const TubeSection *ts);
  double getCurrentOut(const TubeSection *ts);
};


#endif
