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

#include "TriangularGlottis.h"
#include "Constants.h"
#include <cmath>

// ****************************************************************************
// Constructor.
// ****************************************************************************

TriangularGlottis::TriangularGlottis()
{
  int i;

  // ****************************************************************
  // Control parameters.
  // ****************************************************************

  Parameter cp[NUM_CONTROL_PARAMS] =
  {
    { "F0", "Fundamental frequency", "Hz", 1, "Hz", 40.0, 600.0, 120.0, 0.0 },
    { "PR", "Lung pressure", "dPa", 1.0, "dPa", 0.0, 20000.0, 8000.0, 0.0 },
    { "XB", "Lower rest displacement", "cm", 10, "mm", -0.05, 0.3, 0.01, 0.0 },
    { "XT", "Upper rest displacement", "cm", 10, "mm", -0.05, 0.3, 0.01, 0.0 },
    { "AA", "Arytenoid area", "cm^2", 100, "mm^2", -0.1, 0.5, 0.0, 0.0 },
    { "AS", "Aspiration strength", "dB", 1, "dB", -40.0, 0.0, -40.0, -40.0 },
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
    { "LRT", "Lower rest thickness", "cm", 10, "mm", 0.1, 0.5, 0.24, 0.0 },
    { "URT", "Upper rest thickness", "cm", 10, "mm", 0.01, 0.2, 0.06, 0.0 },
    { "LRM", "Lower rest mass", "g", 1, "g", 0.01, 0.2, 0.12, 0.0 },
    { "URM", "Upper rest mass", "g", 1, "g", 0.01, 0.2, 0.03, 0.0 },
    { "LDR", "Lower damping ratio", "", 1, "", 0.0, 3.0, 0.1, 0.0 },
    { "UDR", "Upper damping ratio", "", 1, "", 0.0, 3.0, 0.6, 0.0 },
    { "LSK", "Lower spring k", "dyne/cm", 0.001, "kdyne/cm", 0.0, 400000.0, 80000.0, 0.0 },
    { "USK", "Upper spring k", "dyne/cm", 0.001, "kdyne/cm", 0.0, 400000.0, 8000.0, 0.0 },
    { "LCK", "Lower contact spring k", "dyne/cm", 0.001, "kdyne/cm", 0.0, 500000.0, 240000.0, 0.0 },
    { "UCK", "Upper contact spring k", "dyne/cm", 0.001, "kdyne/cm", 0.0, 500000.0, 24000.0, 0.0 },
    { "CSK", "Coupling spring k", "dyne/cm", 0.001, "kdyne/cm", 0.0, 400000.0, 25000.0, 0.0 },
    { "IL",  "Inlet length", "cm", 10, "mm", 0.0, 1.0, 0.05, 0.0 },
    { "OL",  "Outlet length", "cm", 10, "mm", 0.0, 0.5, 0.01, 0.0 },
    { "NF0", "Natural f0", "Hz", 1, "Hz", 20.0, 500.0, 129.0, 0.0 },
    { "dF0", "Derivative d(F0)/dQ", "Hz", 1, "Hz", 20.0, 500.0, 125.51, 0.0 },
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
    { "LRD", "Lower relative displacement", "cm", 10.0, "mm", -0.3, 0.3, 0.0, 0.0 },
    { "URD", "Upper relative displacement", "cm", 10.0, "mm", -0.3, 0.3, 0.0, 0.0 },
    { "LAD", "Lower absolute displacement", "cm", 10.0, "mm", 0.0, 1.0, 0.0, 0.0 },
    { "UAD", "Upper absolute displacement", "cm", 10.0, "mm", 0.0, 1.0, 0.0, 0.0 },
    { "COL", "Cord length", "cm", 10.0, "mm", 0.0, 3.0, 0.0, 0.0 },
    { "LT", "Lower thickness", "cm", 10.0, "mm", 0.0, 2.0, 0.0, 0.0 },
    { "UT", "Upper thickness", "cm", 10.0, "mm", 0.0, 2.0, 0.0, 0.0 },
    { "LA", "Lower area", "cm^2", 100, "mm^2", 0.0, 4.0, 0.0, 0.0 },
    { "UA", "Upper area", "cm^2", 100, "mm^2", 0.0, 4.0, 0.0, 0.0 },
    { "TQ", "Tension Q", "", 1, "", 0.0, 6.0, 0.0, 0.0 },
    { "CA", "Contact area", "cm^2", 100, "mm^2", 0.0, 2.0, 0.0, 0.0 }
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

string TriangularGlottis::getName()
{
  return string("Triangular glottis");
}


// ****************************************************************************
/// Stops the motion of the vocal folds.
// ****************************************************************************

void TriangularGlottis::resetMotion()
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

void TriangularGlottis::incTime(const double timeIncrement_s, const double pressure_dPa[])
{
  int i;

  double F0 = controlParam[FREQUENCY].x;
  double time_s = (double)pos*timeIncrement_s;
  
  // Add jitter ("flutter") according to Klatt and Klatt (1990)
  // This flutter pattern has been confirmed by measurements of Birkholz!
  F0+= 0.25*(F0/100.0)*(sin(2.0*M_PI*12.7*time_s) + sin(2.0*M_PI*7.1*time_s) + 
                        sin(2.0*M_PI*4.7*time_s));

  double Q = getTensionParameter(F0);

  double subglottalPressure_dPa   = pressure_dPa[0];
  double lowerEdgePressure_dPa    = pressure_dPa[1];
  double upperEdgePressure_dPa    = pressure_dPa[2];
  double supraglottalPressure_dPa = pressure_dPa[3];

  // Low-pass filter the pressure right above the glottis for the 
  // calculation of the vocal fold's vibration amplitude.
//  double supraglottalPressure_dPa = supraglottalPressureFilter.getOutputSample(pressure_dPa[3]);
  

  // ****************************************************************
  // Displacement of the masses.
  // ****************************************************************

  double relX[2];
  double restX[2];
  double backX[2];      // Abs. displacement at the vocal processes
  double prevRelX[2];

  relX[0] = relativeDisplacementBuffer[0][pos & BUFFER_MASK];
  relX[1] = relativeDisplacementBuffer[1][pos & BUFFER_MASK];

  restX[0] = controlParam[REST_DISP_1].x;
  restX[1] = controlParam[REST_DISP_2].x;

  backX[0] = restX[0] + relX[0];
  backX[1] = restX[1] + relX[1];

  prevRelX[0] = relativeDisplacementBuffer[0][(pos-1) & BUFFER_MASK];
  prevRelX[1] = relativeDisplacementBuffer[1][(pos-1) & BUFFER_MASK];

  // ****************************************************************
  // Length, thickness, and mass of the two masses and open-close
  // dimensions.
  // ****************************************************************

  double cordLength;
  double thickness[2];
  double mass[2];
  double openLength[2];
  double contactLength[2];
  double meanOpenWidth[2];
  double meanContactZ[2];

  getLengthAndThickness(Q, cordLength, thickness);
  getOpenCloseDimensions(openLength, contactLength, meanOpenWidth, meanContactZ);

  mass[0] = staticParam[MASS_1].x / Q;
  mass[1] = staticParam[MASS_2].x / Q;

  // ****************************************************************
  // Spring coefficients.
  // ****************************************************************

  double springK[2];
  double contactSpringK[2];
  // alpha is the proportion of the glottis that is in contact (0 <= alpha <= 1).
  double alpha[2];

  alpha[0] = contactLength[0] / cordLength;
  alpha[1] = contactLength[1] / cordLength;

  springK[0] = staticParam[SPRING_K_1].x * Q;
  springK[1] = staticParam[SPRING_K_2].x * Q;
  contactSpringK[0] = staticParam[CONTACT_SPRING_K_1].x * Q;
  contactSpringK[1] = staticParam[CONTACT_SPRING_K_2].x * Q;

  // In the original implementation in VTL 2.1 and also in the original 
  // two-mass model by Ishizaka and Flanagan (1972), the coupling spring
  // constant varied with Q^2. However, we found that this prevents a
  // reliable oscillation at high frequencies (> 400 Hz). In constrast,
  // a constant coupling spring (independent of Q) allows oscillation 
  // up to 600 Hz and possibly higher!

  double kc = staticParam[COUPLING_SPRING_K].x; // *Q*Q;

  // ****************************************************************
  // Damping resistance of the two masses.
  // ****************************************************************

  double dampingRatio[2];

  // Add a damping ratio value of 1.0 ("add critical damping"), when
  // the glottis is fully closed (alpha = 1)
  dampingRatio[0] = staticParam[DAMPING_RATIO_1].x + alpha[0]*1.0;
  dampingRatio[1] = staticParam[DAMPING_RATIO_2].x + alpha[1]*1.0;

  double r[2];
  r[0] = 2.0*dampingRatio[0]*sqrt( mass[0]*springK[0] );
  r[1] = 2.0*dampingRatio[1]*sqrt( mass[1]*springK[1] );

  // ****************************************************************
  // Driving forces for the two masses.
  // For the closed part of the masses, the subglottal pressure is
  // always the driving force.
  // ****************************************************************

  double force[2];

  // The pressures acting on the faces of the masses.

  force[0] = lowerEdgePressure_dPa*openLength[0]*thickness[0];
  force[1] = upperEdgePressure_dPa*openLength[1]*thickness[1];

  // The forces on the masses due to the pressures at the inlet and 
  // outlet regions.

  double inletLength = staticParam[INLET_LENGTH].x;
  double outletLength = staticParam[OUTLET_LENGTH].x;
  // One 0.5 because of torque balance, the other 0.5 because of mean value.
  force[0]+= 0.5*0.5*(subglottalPressure_dPa + lowerEdgePressure_dPa)*inletLength*cordLength;
  force[1]+= 0.5*0.5*(supraglottalPressure_dPa + upperEdgePressure_dPa)*outletLength*cordLength;


  // ****************************************************************
  // restXStar[] is the lateral rest position at the mid-point of the
  // closed vocal fold section.
  // ****************************************************************

  double restXStar[2];
  for (i=0; i < 2; i++)
  {
    // When the vocal folds are apart, the glottis is assumed to be
    // triangular.
    if (restX[i] >= 0.0)
    {
      restXStar[i] = restX[i]*(1.0 - meanContactZ[i]/cordLength);
    }
    else
    // When the vocal folds are pressed together, they are assumed
    // to be parallel.
    {
      restXStar[i] = restX[i];
    }
  }

  // ****************************************************************
  // Calculate the new positions.
  // ****************************************************************

  double T = timeIncrement_s;
  
  double A = mass[0] + r[0]*T + T*T*(springK[0] + contactSpringK[0]*alpha[0]) + kc*T*T;
  double B = -kc*T*T;
  double C = -kc*T*T;
  double D = mass[1] + r[1]*T + T*T*(springK[1] + contactSpringK[1]*alpha[1]) + kc*T*T;

  double E = force[0]*T*T + 2.0*mass[0]*relX[0] - mass[0]*prevRelX[0] + r[0]*T*relX[0] -
    T*T*contactSpringK[0]*alpha[0]*restXStar[0];
  
  double F = force[1]*T*T + 2.0*mass[1]*relX[1] - mass[1]*prevRelX[1] + r[1]*T*relX[1] -
    T*T*contactSpringK[1]*alpha[1]*restXStar[1];

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

void TriangularGlottis::calcGeometry()
{
  // ****************************************************************
  // Calculate displacement values.
  // ****************************************************************

  double relX[2];
  double backX[2];

  relX[0] = relativeDisplacementBuffer[0][pos & BUFFER_MASK];
  relX[1] = relativeDisplacementBuffer[1][pos & BUFFER_MASK];
  backX[0] = controlParam[REST_DISP_1].x + relX[0];
  backX[1] = controlParam[REST_DISP_2].x + relX[1];

  // ****************************************************************
  // Area of the chink between the arytenoids.
  // ****************************************************************

  double chinkArea = controlParam[ARY_AREA].x;
  if (chinkArea < 0.0)
  {
    chinkArea = 0.0;
  }

  // ****************************************************************
  // ****************************************************************

  double Q = getTensionParameter( controlParam[FREQUENCY].x );
  double cordLength;
  double thickness[2];
  double area[2];
  double openLength[2];
  double contactLength[2];
  double meanOpenWidth[2];
  double meanContactZ[2];
  double contactArea;

  getLengthAndThickness(Q, cordLength, thickness);
  getOpenCloseDimensions(openLength, contactLength, meanOpenWidth, meanContactZ);
  area[0] = openLength[0]*meanOpenWidth[0] + chinkArea;
  area[1] = openLength[1]*meanOpenWidth[1] + chinkArea;

  contactArea = getContactArea(backX, relX, cordLength, thickness[0] + thickness[1]);

  // ****************************************************************
  // Set the derived values.
  // ****************************************************************

  derivedParam[RELATIVE_DISP_1].x = relX[0];
  derivedParam[RELATIVE_DISP_2].x = relX[1];
  derivedParam[ABSOLUTE_DISP_1].x = backX[0];
  derivedParam[ABSOLUTE_DISP_2].x = backX[1];
  derivedParam[CURRENT_LENGTH].x = cordLength;
  derivedParam[CURRENT_THICKNESS_1].x = thickness[0];
  derivedParam[CURRENT_THICKNESS_2].x = thickness[1];
  derivedParam[CURRENT_AREA_1].x = area[0];
  derivedParam[CURRENT_AREA_2].x = area[1];
  derivedParam[CURRENT_TENSION].x = Q;
  derivedParam[CONTACT_AREA].x = contactArea;
}


// ****************************************************************************
/// Returns the length and area of the two glottal tube sections.
// ****************************************************************************

void TriangularGlottis::getTubeData(double *length_cm, double *area_cm2)
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

int TriangularGlottis::getApertureParamIndex()
{
  return REST_DISP_2;
}


// ****************************************************************************
// ****************************************************************************

double TriangularGlottis::getAspirationStrength_dB()
{
  return controlParam[ASPIRATION_STRENGTH].x;
}


// ****************************************************************************
/// Calculates the tension parameter Q from the given fundamental frequency.
// ****************************************************************************

double TriangularGlottis::getTensionParameter(double f0)
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

void TriangularGlottis::getLengthAndThickness(const double Q, double &length_cm, double thickness[])
{
  double factor = sqrt(Q);
  length_cm = staticParam[REST_LENGTH].x * factor;
  thickness[0] = staticParam[REST_THICKNESS_1].x / factor;
  thickness[1] = staticParam[REST_THICKNESS_2].x / factor;
}

// ****************************************************************************
/// Returns some measures of the triangular glottis.
// ****************************************************************************

void TriangularGlottis::getOpenCloseDimensions(double openLength[], double contactLength[], 
  double meanOpenWidth[], double meanContactZ[])
{
  int i;
  double Q = getTensionParameter( controlParam[FREQUENCY].x );
  double cordLength;
  double thickness[2];

  // Get the current length and thickness of the masses.
  getLengthAndThickness(Q, cordLength, thickness);

  // Absolute displacements the the frontal and dorsal ends.
  double restX[2];
  double relX[2];
  double backX[2];
  double frontX[2];

  restX[0] = controlParam[REST_DISP_1].x;
  restX[1] = controlParam[REST_DISP_2].x;
  relX[0]  = relativeDisplacementBuffer[0][pos & BUFFER_MASK];
  relX[1]  = relativeDisplacementBuffer[1][pos & BUFFER_MASK];
  backX[0] = restX[0] + relX[0];
  backX[1] = restX[1] + relX[1];

  // When the vocal folds are pressed together, they are assumed to 
  // be parallel, and otherwise triangular.

  if (restX[0] < 0.0)
  {
    frontX[0] = backX[0];
  }
  else
  {
    frontX[0] = relX[0];
  }
  
  if (restX[1] < 0.0)
  {
    frontX[1] = backX[1];
  }
  else
  {
    frontX[1] = relX[1];
  }

  // ****************************************************************

  for (i=0; i < 2; i++)
  {
    // Default is a completely closed glottis
    openLength[i] = 0.0;
    meanOpenWidth[i] = 0.0;
    contactLength[i] = cordLength;
    meanContactZ[i] = 0.5*cordLength;

    // Left and right edge don't touch.
    if ((backX[i] > 0.0) && (frontX[i] > 0.0))
    {
      openLength[i] = cordLength;
      meanOpenWidth[i] = backX[i] + frontX[i];
      contactLength[i] = 0.0;
      meanContactZ[i] = 0.0;
    }
    else
    // Glottis completely closed.
    if ((backX[i] <= 0.0) && (frontX[i] <= 0.0))
    {
      openLength[i] = 0.0;
      meanOpenWidth[i] = 0.0;
      contactLength[i] = cordLength;
      meanContactZ[i] = 0.5*cordLength;
    }
    else
    {
      const double EPSILON = 0.000000001;
      if (fabs(restX[i]) < EPSILON)
      {
        restX[i] = EPSILON;
      }
      double zApex = cordLength*(1.0 + relX[i] / restX[i]);
      if ((zApex >= 0.0) && (zApex <= cordLength))
      {
        // Glottis partly opened at the back
        if (backX[i] > 0.0)
        {
          openLength[i] = zApex;
          meanOpenWidth[i] = backX[i];
          contactLength[i] = cordLength - zApex;
          meanContactZ[i] = 0.5*(zApex + cordLength);
        }
        else
        // Glottis partly opened at the front
        {
          openLength[i] = cordLength - zApex;
          meanOpenWidth[i] = frontX[i];
          contactLength[i] = zApex;
          meanContactZ[i] = 0.5*zApex;
        }
      }
    }
  }

}


// ****************************************************************************
/// Returns the contact area  between the left and right vocal fold.
// ****************************************************************************

double TriangularGlottis::getContactArea(double backX[], double frontX[], double length, double thickness)
{
  const double EPSILON = 0.000000001;
  const int NUM_SLICES = 32;
  double dy = thickness / (double)NUM_SLICES;
  int i;
  double x0, x1;
  double t;
  double temp;
  double dA;
  double area = 0.0;
  double denom;

  // Run through several thin slices in vertical direction.

  for (i=0; i < NUM_SLICES; i++)
  {
    t = (double)i / (double)NUM_SLICES;
    x0 = (1.0-t)*backX[1] + t*backX[0];
    x1 = (1.0-t)*frontX[1] + t*frontX[0];

    // Let x0 be the smaller of x0 and x1
    if (x0 >= x1)
    {
      double dummy = x0;
      x0 = x1;
      x1 = dummy;
    }

    dA = 0.0;
    if ((x0 <= 0.0) && (x1 <= 0.0))
    {
      dA = length*dy;      
    }
    else
    if ((x0 <= 0.0) && (x1 > 0.0))
    {
      denom = x1 - x0;
      if (denom < EPSILON)
      {
        denom = EPSILON;
      }
      temp = -x0*length / denom;
      dA = temp*dy;
    }

    area+= dA;
  }

  return area;
}

// ****************************************************************************
