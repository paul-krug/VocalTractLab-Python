// ****************************************************************************
// This file is part of VocalTractLab.
// Copyright (C) 2020, Peter Birkholz, Dresden, Germany
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

#include "TwoMassModel.h"
#include "Constants.h"
#include <cmath>

// ****************************************************************************
// Constructor.
// ****************************************************************************

TwoMassModel::TwoMassModel()
{
  int i;

  // ****************************************************************
  // Control parameters.
  // ****************************************************************

  Parameter cp[NUM_CONTROL_PARAMS] =
  {
    { "F0", "Fundamental frequency", "Hz", 1, "Hz", 40.0, 600.0, 120.0, 0.0 },
    { "PR", "Lung pressure", "dPa", 1.0, "dPa", 0.0, 20000.0, 8000.0, 0.0 },
    { "XB", "Lower rest displacement", "cm", 10.0, "mm", -0.05, 0.3, 0.01, 0.0 },
    { "XT", "Upper rest displacement", "cm", 10.0, "mm", -0.05, 0.3, 0.01, 0.0 },
    { "EAA", "Extra arytenoid area", "cm^2", 100.0, "mm^2", -0.25, 0.25, 0.0, 0.0 },
    { "DF", "Damping factor", "", 1, "", 0.3, 3.0, 1.0, 0.0 },
  };

  controlParam.clear();
  for (i=0; i < NUM_CONTROL_PARAMS; i++)
  {
    controlParam.push_back( cp[i] );
    controlParam[i].x = controlParam[i].neutral;
  }

  // ****************************************************************
  // Static parameters.
  // ****************************************************************

  Parameter sp[NUM_STATIC_PARAMS] =
  {
    { "COL", "Cord length", "cm", 10, "mm", 0.5, 2.0, 1.3, 0.0 },
    { "LRT", "Lower rest thickness", "cm", 10, "mm", 0.1, 0.5, 0.25, 0.0 },
    { "URT", "Upper rest thickness", "cm", 10, "mm", 0.01, 0.2, 0.05, 0.0 },
    { "LRM", "Lower rest mass", "g", 1, "g", 0.01, 0.2, 0.125, 0.0 },
    { "URM", "Upper rest mass", "g", 1, "g", 0.01, 0.2, 0.025, 0.0 },
    { "LDR", "Lower damping ratio", "", 1, "", 0.0, 3.0, 0.1, 0.0 },
    { "UDR", "Upper damping ratio", "", 1, "", 0.0, 3.0, 0.6, 0.0 },
    { "LSK", "Lower spring k", "dyne/cm", 0.001, "kdyne/cm", 0.0, 400000.0, 80000.0, 0.0 },
    { "USK", "Upper spring k", "dyne/cm", 0.001, "kdyne/cm", 0.0, 400000.0, 8000.0, 0.0 },
    { "LSE", "Lower spring eta", "1/cm^2", 1.0, "1/cm^2", 0.0, 1000.0, 100.0, 0.0 },
    { "USE", "Upper spring eta", "1/cm^2", 1.0, "1/cm^2", 0.0, 1000.0, 100.0, 0.0 },
    { "LCK", "Lower contact spring k", "dyne/cm", 0.001, "kdyne/cm", 0.0, 400000.0, 240000.0, 0.0 },
    { "UCK", "Upper contact spring k", "dyne/cm", 0.001, "kdyne/cm", 0.0, 400000.0, 24000.0, 0.0 },
    { "LCE", "Lower contact spring eta", "1/cm^2", 1.0, "1/cm^2", 0.0, 1000.0, 500.0, 0.0 },
    { "UCE", "Upper contact spring eta", "1/cm^2", 1.0, "1/cm^2", 0.0, 1000.0, 500.0, 0.0 },
    { "CSK", "Coupling spring k", "dyne/cm", 0.001, "kdyne/cm", 0.0, 400000.0, 25000.0, 0.0 },
    { "CRW", "Critical width", "cm", 10, "mm", 0.0, 0.2, 0.0, 0.0 },
    { "NF0", "Natural F0", "Hz", 1, "Hz", 20.0, 500.0, 158.0, 0.0 },
    { "dF0", "Derivative d(F0)/dQ", "Hz", 1, "Hz", 20.0, 500.0, 100.0, 0.0 },
    { "CL", "Chink length", "cm", 10, "mm", 0.0, 0.5, 0.2, 0.0 }
  };

  staticParam.clear();
  for (i=0; i < NUM_STATIC_PARAMS; i++)
  {
    staticParam.push_back( sp[i] );
    staticParam[i].x = staticParam[i].neutral;
  }

  // ****************************************************************
  // Derived parameters.
  // ****************************************************************

  Parameter dp[NUM_DERIVED_PARAMS] =
  {
    { "LRD", "Lower relative displacement", "cm", 10, "mm", -0.3, 0.3, 0.0, 0.0 },
    { "URD", "Upper relative displacement", "cm", 10, "mm", -0.3, 0.3, 0.0, 0.0 },
    { "LAD", "Lower absolute displacement", "cm", 10, "mm", 0.0, 1.0, 0.0, 0.0 },
    { "UAD", "Upper absolute displacement", "cm", 10, "mm", 0.0, 1.0, 0.0, 0.0 },
    { "COL", "Cord length", "cm", 10, "mm", 0.0, 3.0, 0.0, 0.0 },
    { "LT", "Lower thickness", "cm", 10, "mm", 0.0, 2.0, 0.0, 0.0 },
    { "UT", "Upper thickness", "cm", 10, "mm", 0.0, 2.0, 0.0, 0.0 },
    { "LA", "Lower area", "cm^2", 100, "mm^2", 0.0, 4.0, 0.0, 0.0 },
    { "UA", "Upper area", "cm^2", 100, "mm^2", 0.0, 4.0, 0.0, 0.0 },
    { "CW", "Chink width", "cm", 10, "mm", 0.0, 1.0, 0.0, 0.0 },
    { "TQ", "Tension Q", "", 1, "", 0.0, 6.0, 0.0, 0.0 }
  };

  derivedParam.clear();
  for (i=0; i < NUM_DERIVED_PARAMS; i++)
  {
    derivedParam.push_back( dp[i] );
    derivedParam[i].x = derivedParam[i].neutral;
  }

  // ****************************************************************
  // Create a shape for the default parameter values.
  // ****************************************************************

  Shape s;
  s.name = "default";
  s.controlParam.resize(NUM_CONTROL_PARAMS);
  for (i=0; i < NUM_CONTROL_PARAMS; i++)
  {
    s.controlParam[i] = controlParam[i].neutral;
  }

  shape.push_back(s);

  // ****************************************************************
  // Stop the motion and calculate the geometry.
  // ****************************************************************

  resetMotion();
  calcGeometry();

  // Declare the current state as saved.
  clearUnsavedChanges();
}


// ****************************************************************************
// ****************************************************************************

string TwoMassModel::getName()
{
  return string("Two-mass model");
}


// ****************************************************************************
/// Stops the motion of the vocal folds.
// ****************************************************************************

void TwoMassModel::resetMotion()
{
  int i;

  for (i=0; i < BUFFER_LENGTH; i++)
  {
    relativeDisplacementBuffer[0][i] = 0.0;
    relativeDisplacementBuffer[1][i] = 0.0;
  }

  pos = 0;

  supraglottalPressureFilter.createChebyshev(25.0/(double)SAMPLING_RATE, false, 4);
  supraglottalPressureFilter.resetBuffers();
}


// ****************************************************************************
/// Performs a time step of the digital simulation.
/// Requires four pressure values: subglottal, lower glottis, upper glottis, 
/// supraglottal.
// ****************************************************************************

void TwoMassModel::incTime(const double timeIncrement_s, const double pressure_dPa[])
{
  int i;

  double F0 = controlParam[FREQUENCY].x;
  double time_s = (double)pos*timeIncrement_s;

  // Add jitter ("flutter") according to Klatt and Klatt (1990)
  // This flutter pattern has been confirmed by measurements of Birkholz!
  F0 += 0.25*(F0 / 100.0)*(sin(2.0*M_PI*12.7*time_s) + sin(2.0*M_PI*7.1*time_s) +
    sin(2.0*M_PI*4.7*time_s));

  double Q = getTensionParameter(F0);

  double subglottalPressure_dPa   = pressure_dPa[0];
  double lowerEdgePressure_dPa    = pressure_dPa[1];
  double upperEdgePressure_dPa    = pressure_dPa[2];
  double supraglottalPressure_dPa = pressure_dPa[3];

  // criticalWidth is the width of the glottis, at which they are
  // supposed to be in mechanical contact. This width is zero in
  // the original model, but greater than zero in Pelorson (1996).

  double criticalWidth = staticParam[CRITICAL_WIDTH].x;
  double criticalX = 0.5*criticalWidth;

  // ****************************************************************
  // Displacement of the masses.
  // ****************************************************************

  double relX[2];
  double restX[2];
  double absX[2];
  double minRelX[2];
  double prevRelX[2];

  relX[0] = relativeDisplacementBuffer[0][pos & BUFFER_MASK];
  relX[1] = relativeDisplacementBuffer[1][pos & BUFFER_MASK];

  restX[0] = controlParam[REST_DISP_1].x;
  restX[1] = controlParam[REST_DISP_2].x;

  absX[0] = restX[0] + relX[0];
  absX[1] = restX[1] + relX[1];

  minRelX[0] = criticalX - restX[0];
  minRelX[1] = criticalX - restX[1];

  prevRelX[0] = relativeDisplacementBuffer[0][(pos-1) & BUFFER_MASK];
  prevRelX[1] = relativeDisplacementBuffer[1][(pos-1) & BUFFER_MASK];

  // ****************************************************************
  // Length, thickness, and mass of the two masses.
  // ****************************************************************

  double cordLength;
  double thickness[2];
  getLengthAndThickness(Q, cordLength, thickness);

  double mass[2];
  mass[0] = staticParam[MASS_1].x / Q;
  mass[1] = staticParam[MASS_2].x / Q;

  // ****************************************************************
  // Spring coefficients.
  // ****************************************************************

  double springK[2];
  double springEta[2];
  double contactSpringK[2];
  double contactSpringEta[2];

  springK[0] = staticParam[SPRING_K_1].x * Q;
  springK[1] = staticParam[SPRING_K_2].x * Q;
  springEta[0] = staticParam[SPRING_ETA_1].x;
  springEta[1] = staticParam[SPRING_ETA_2].x;

  contactSpringK[0] = staticParam[CONTACT_SPRING_K_1].x * Q;
  contactSpringK[1] = staticParam[CONTACT_SPRING_K_2].x * Q;
  contactSpringEta[0] = staticParam[CONTACT_SPRING_ETA_1].x;
  contactSpringEta[1] = staticParam[CONTACT_SPRING_ETA_2].x;

  double kc = staticParam[COUPLING_SPRING_K].x * Q*Q;

  // ****************************************************************
  // Damping resistance of the two masses.
  // ****************************************************************

  double dampingFactor = controlParam[DAMPING_FACTOR].x;
  double dampingRatio[2];

  dampingRatio[0] = staticParam[DAMPING_RATIO_1].x;
  dampingRatio[1] = staticParam[DAMPING_RATIO_2].x;    

  if (absX[0] <= criticalX)
  {
    // Add a damping ration value of 1.0 ("add critical damping")
    dampingRatio[0]+= 1.0;
  }
  if (absX[1] <= criticalX)
  {
    // Add a damping ration value of 1.0 ("add critical damping")
    dampingRatio[1]+= 1.0;
  }

  double r[2];
  // Include the squared damping factor similar to Sondhi and Schroeter (1987).
  r[0] = 2.0*dampingRatio[0]*sqrt( mass[0]*springK[0] ) * dampingFactor*dampingFactor;
  r[1] = 2.0*dampingRatio[1]*sqrt( mass[1]*springK[1] ) * dampingFactor*dampingFactor;

  // ****************************************************************
  // Driving forces for the two masses.
  // ****************************************************************

  double force[2];

  // Upper and lower slit open
  if ((absX[0] > criticalWidth) && (absX[1] > criticalWidth))
  {
    force[0] = lowerEdgePressure_dPa*cordLength*thickness[0];
    force[1] = upperEdgePressure_dPa*cordLength*thickness[1];
  }
  else

  // Lower slit closed, upper slit open
  if ((absX[0] <= criticalWidth) && (absX[1] > criticalWidth))
  {
    force[0] = subglottalPressure_dPa*cordLength*thickness[0];
    force[1] = upperEdgePressure_dPa*cordLength*thickness[1];
  }
  else

  // Lower slit open, upper slit closed
  if ((absX[0] > criticalWidth) && (absX[1] <= criticalWidth))
  {
    force[0] = lowerEdgePressure_dPa*cordLength*thickness[0];
    force[1] = lowerEdgePressure_dPa*cordLength*thickness[1];
  }
  else

  // Lower and upper slit closed
  {
    force[0] = subglottalPressure_dPa*cordLength*thickness[0];
    force[1] = supraglottalPressure_dPa*cordLength*thickness[1];
  }

  // ****************************************************************
  // Calculate the new positions.
  // ****************************************************************

  // The nonlinear restoring forces based on the existing samples
  
  double deltaX;
  double nonlinearRestoringForce[2];

  for (i=0; i < 2; i++)
  {
    // Switch off the contact spring forces?
    if (relX[i] > minRelX[i])
    {
      contactSpringK[i] = 0.0;
      contactSpringEta[i] = 0.0;      
    }

    deltaX = relX[i] - minRelX[i];
    nonlinearRestoringForce[i] = springK[i]*springEta[i]*relX[i]*relX[i]*relX[i] +
      contactSpringK[i]*contactSpringEta[i]*deltaX*deltaX*deltaX;
  }

  double T = timeIncrement_s;
  
  double A = mass[0] + r[0]*T + T*T*(springK[0] + contactSpringK[0]) + kc*T*T;
  double B = -kc*T*T;
  double C = -kc*T*T;
  double D = mass[1] + r[1]*T + T*T*(springK[1] + contactSpringK[1]) + kc*T*T;
  double E = force[0]*T*T + 2.0*mass[0]*relX[0] - mass[0]*prevRelX[0] + r[0]*T*relX[0] +
    T*T*contactSpringK[0]*minRelX[0] - nonlinearRestoringForce[0]*T*T;
  double F = force[1]*T*T + 2.0*mass[1]*relX[1] - mass[1]*prevRelX[1] + r[1]*T*relX[1] +
    T*T*contactSpringK[1]*minRelX[1] - nonlinearRestoringForce[1]*T*T;

  double det = A*D - B*C;
  const double EPSILON = 0.000000001;
  if (fabs(det) < EPSILON)
  {
    det = EPSILON;
  }

  double newX[2];
  newX[0] = (E*D - B*F) / det;
  newX[1] = (A*F - E*C) / det;

  relativeDisplacementBuffer[0][(pos+1) & BUFFER_MASK] = newX[0];
  relativeDisplacementBuffer[1][(pos+1) & BUFFER_MASK] = newX[1];


  // ****************************************************************
  // Increment the current position.
  // ****************************************************************

  pos++;
}


// ****************************************************************************
/// Calculates the current glottis geometry from the static and control
/// parameters.
// ****************************************************************************

void TwoMassModel::calcGeometry()
{
  int i;

  // ****************************************************************
  // Calculate displacement values.
  // ****************************************************************

  double restX[2];
  double relX[2];
  double absX[2];

  restX[0] = controlParam[REST_DISP_1].x;
  restX[1] = controlParam[REST_DISP_2].x;

  relX[0] = relativeDisplacementBuffer[0][pos & BUFFER_MASK];
  relX[1] = relativeDisplacementBuffer[1][pos & BUFFER_MASK];

  for (i=0; i < 2; i++)
  {
    absX[i] = restX[i] + relX[i];
    if (absX[i] < 0.0)
    {
      absX[i] = 0.0;
    }
  }
  
  // ****************************************************************
  // Area of the chink between the arytenoids.
  // ****************************************************************

  double passiveChinkWidth = 2.0*controlParam[REST_DISP_2].x;
  if (passiveChinkWidth < 0.0) 
  { 
    passiveChinkWidth = 0.0; 
  }
  double chinkArea = passiveChinkWidth * staticParam[CHINK_LENGTH].x + controlParam[ARY_AREA].x;
  if (chinkArea < 0.0)
  {
    chinkArea = 0.0;
  }

  // ****************************************************************
  // ****************************************************************

  double Q = getTensionParameter( controlParam[FREQUENCY].x );
  double cordLength;
  double thickness[2];

  getLengthAndThickness(Q, cordLength, thickness);

  double area[2];
  area[0] = 2.0*cordLength*absX[0] + chinkArea;
  area[1] = 2.0*cordLength*absX[1] + chinkArea;
  

  // ****************************************************************
  // Set the derived values.
  // ****************************************************************

  double chinkLength = staticParam[CHINK_LENGTH].x;
  if (chinkLength < 0.000001)
  {
    chinkLength = 0.000001;
  }

  derivedParam[RELATIVE_DISP_1].x = relX[0];
  derivedParam[RELATIVE_DISP_2].x = relX[1];
  derivedParam[ABSOLUTE_DISP_1].x = absX[0];
  derivedParam[ABSOLUTE_DISP_2].x = absX[1];
  derivedParam[CURRENT_LENGTH].x = cordLength;
  derivedParam[CURRENT_THICKNESS_1].x = thickness[0];
  derivedParam[CURRENT_THICKNESS_2].x = thickness[1];
  derivedParam[CURRENT_AREA_1].x = area[0];
  derivedParam[CURRENT_AREA_2].x = area[1];
  derivedParam[CHINK_WIDTH].x = chinkArea / chinkLength;
  derivedParam[CURRENT_TENSION].x = Q;
}


// ****************************************************************************
/// Returns the length and area of the two glottal tube sections.
// ****************************************************************************

void TwoMassModel::getTubeData(double *length_cm, double *area_cm2)
{
  length_cm[0] = derivedParam[CURRENT_THICKNESS_1].x;
  length_cm[1] = derivedParam[CURRENT_THICKNESS_2].x;

  area_cm2[0] = derivedParam[CURRENT_AREA_1].x;
  area_cm2[1] = derivedParam[CURRENT_AREA_2].x;
}


// ****************************************************************************
/// Returns the index of the parameter that best represents the glottal
/// aperture.
// ****************************************************************************

int TwoMassModel::getApertureParamIndex()
{
  return REST_DISP_2;
}

// ****************************************************************************
/// Calculates the tension parameter Q from the given fundamental frequency.
// ****************************************************************************

double TwoMassModel::getTensionParameter(double f0)
{
  double f0DivQ = staticParam[F0_DIV_Q].x;
  
  // Avoid division by zero.
  if (f0DivQ < 0.000001)
  {
    f0DivQ = 0.000001;
  }

  double Q = 1.0 + (f0 - staticParam[NATURAL_F0].x) / f0DivQ;
  if (Q < 0.05)
  {
    Q = 0.05;
  }

  return Q;
}

// ****************************************************************************
/// Returns the length of the vocal cords and the thickness of the lower and 
/// the upper mass for the given tension parameter Q.
// ****************************************************************************

void TwoMassModel::getLengthAndThickness(const double Q, double &length_cm, double thickness[])
{
  double factor = sqrt(Q);
  length_cm = staticParam[REST_LENGTH].x * factor;
  thickness[0] = staticParam[REST_THICKNESS_1].x / factor;
  thickness[1] = staticParam[REST_THICKNESS_2].x / factor;
}

// ****************************************************************************
