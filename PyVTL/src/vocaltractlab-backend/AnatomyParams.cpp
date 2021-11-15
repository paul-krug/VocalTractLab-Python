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

#include "AnatomyParams.h"
#include "Dsp.h"

VocalTract *AnatomyParams::referenceVocalTract = new VocalTract();


// ****************************************************************************
// Constructor.
// ****************************************************************************

AnatomyParams::AnatomyParams()
{
  Param p[NUM_ANATOMY_PARAMS] =
  {
    { "Lip width", "W0", "cm", 0.5, 1.5, 0.5 },
    { "Mandible height", "H5", "cm", 0.9, 1.8, 0.9 },
    { "Lower molars height", "H4", "cm", 0.3, 0.7, 0.3 },
    { "Upper molars height", "H3", "cm", 0.3, 0.7, 0.3 },
    { "Palate height", "H2", "cm", 0.8, 2.0, 0.8 },
    { "Palate depth", "D0", "cm", 2.6, 5.0, 2.6 }, 
    { "Hard palate length", "W1 minus incisor width", "cm", 2.2, 5.2, 2.2 }, 
    { "Soft palate length", "W2", "cm", 2.3, 3.1, 2.3 },
    { "Pharynx length", "H0", "cm", 3.5, 8.0, 3.5 }, 
    { "Larynx length", "H1", "cm", 2.0, 4.0, 2.0 }, 
    { "Larynx width", "W4", "cm", 2.0, 3.5, 2.0 },
    { "Vocal fold length", "W3", "cm", 0.5, 2.0, 0.5 },
    { "Oral-pharyngeal angle", "A0", "deg", -105.0, -90.0, -105.0 },
  };

  int i;
  for (i=0; i < NUM_ANATOMY_PARAMS; i++)
  {
    param[i] = p[i];
  }

  // Get the initial parameter values from the original (unmodified) vocal tract
  getFrom(referenceVocalTract);
}

// ****************************************************************************
// Restrict the values of all parameters.
// ****************************************************************************

void AnatomyParams::restrictParams()
{
  int i;

  for (i=0; i < NUM_ANATOMY_PARAMS; i++)
  {
    if (param[i].x < param[i].min) 
    { 
      param[i].x = param[i].min; 
    }
    
    if (param[i].x > param[i].max) 
    { 
      param[i].x = param[i].max;
    }
  }

  // Force some additional constraints.

  double d = param[PALATE_HEIGHT].x + param[UPPER_MOLARS_HEIGHT].x + 
             param[LOWER_MOLARS_HEIGHT].x + param[MANDIBLE_HEIGHT].x;
  
  if (param[PHARYNX_LENGTH].x < d) 
  { 
    param[PHARYNX_LENGTH].x = d; 
  }
}

// ****************************************************************************
// Calc. the parameter values from a given age and gender.
// ****************************************************************************

void AnatomyParams::calcFromAge(int age_month, bool isMale)
{
  double ANGFHJ, ANGPJ, SNHY, ATLPTM, SNPNS, ANSPNS, LIPS, THIHY;
  double HYGLOTT, ANGFHV, LENVC, SYMPH, PALWI, PALHIT;
  double A, B, C, D, E, F;

  // 1 year is the lowest age.

  if (age_month < 12) 
  { 
    age_month = 12; 
  }

  double t = (double)age_month / 12.0;      // Time in years

  // ****************************************************************
  // ANGFHJ:
  // Angle of Frankfort horizontal with respect to the body of the 
  // jaw.
  // ****************************************************************

  A = 26.856;
  B = -0.170;
  ANGFHJ = A + B*t;

  // ****************************************************************
  // ANGFHV:
  // Angle in degrees of posterior pharynx wall with respect to
  // a line perpendicular to the Frankfort horizontal.
  // ****************************************************************

  if (age_month < 30)
  {
    if (isMale) 
    { 
      A = 4.88; 
      B = 1.467; 
    } 
    else 
    { 
      A = 13.037; 
      B = 1.305; 
    }
  }
  else
  {
    if (isMale) 
    { 
      A = -2.086; 
      B = 0.545; 
    } 
    else 
    { 
      A = 4.383; 
      B = 0.203; 
    }
  }

  ANGFHV = A + B*t;

  // ****************************************************************
  // ANGPJ:
  // Angle of palatal line with respect to the body of the jaw.
  // ****************************************************************

  A = 27.85;
  B = -0.401;
  ANGPJ = A + B*t;

  // ****************************************************************
  // ANSPNS:
  // Distance from the anterior nasal spine to the posterior nasal
  // spine. 
  // ****************************************************************

  if (isMale)
  {
    A = 49.152;
    B = -0.889;
    C = 0.512;
    D = 6.002;
    E = 10.395;
    F = 0.920;
  }
  else
  {
    A = 44.414;
    B = -1.210;
    C = 0.699;
    D = 7.807;
    E = 5.218;
    F = 0.727;
  }

  ANSPNS = A/(1.0 + exp(B-C*C*t)) + D/(1.0+exp(E-F*F*t));
  ANSPNS/= 10.0;    // Transformation from mm in cm

  // ****************************************************************
  // ATLPTM:
  // Distance from the anterior tubercle of the atlas to the 
  // pterygomaxillary fissure measured parallel to the Frankfort 
  // horizontal. 
  // ****************************************************************

  if (isMale) 
  { 
    A = 27.043; 
    B = 0.213; 
  } 
  else 
  { 
    A = 27.727; 
    B = 0.094; 
  }

  ATLPTM = A + B*t;
  ATLPTM/= 10.0;    // from mm to cm

  // ****************************************************************
  // HYGLOTT:
  // Distance from lower border of hyoid bone to the glottis measured
  // along the dorsal surface of the epiglottis.
  // ****************************************************************

  if (isMale)
  {
    A = 8.351;
    B = -4.020;
    C = -1.316;
    D = 2.756;
    E = 0.083;
    HYGLOTT = A + B*exp(C*t) + D*exp(E*t);
  }
  else
  {
    A = 11.143;
    B = 0.298;
    C = 1.4;
    D = -1.235;
    HYGLOTT = A + B*t - exp(C + D*t);
  }
  HYGLOTT/= 10.0;   // from mm in cm

  // ****************************************************************
  // LENVC:
  // Length of the vocal cords.
  // ****************************************************************

  if (isMale)
  {
    A = 7.618;
    B = -0.311;
    C = 1.14;
    D = 11.05;
    E = 13.003;
    F = 0.939;
    LENVC = A/(1+exp(B - C*C*t)) + D/(1+exp(E-F*F*t));
  }
  else
  {
    A = 5.997;
    B = 0.39;
    C = 0.804;
    D = -2.435;
    LENVC = A + B*t - exp(C + D*t);
  }
  LENVC/= 10.0;   // from mm in cm

  // ****************************************************************
  // LIPS:
  // Horizontal thickness of soft tissue covering the anterior edge 
  // of the maxilla at the indentation below the anterior nasal spine.
  // ****************************************************************

  if (isMale)
  {
    A = 27.439;
    B = -17.612;
    C = -0.030;
  }
  else
  {
    A = 22.07;
    B = -13.41;
    C = -0.041;
  }
  LIPS = A + B*exp(C*t);
  LIPS/= 10.0;

  // ****************************************************************
  // PALHIT:
  // Height of the palate from the margin between the teeth and the 
  // gums to the highest point on the arch of the hard palate,
  // measured normal to the occlusal plane.
  // ****************************************************************

  A = 17.813;
  B = -8.859;
  C = -0.077;
  PALHIT = A + B*exp(C*t);
  PALHIT/= 10.0;

  // ****************************************************************
  // PALWI:
  // Width of the palate measured between the lingual surfaces of the
  // molars.
  // ****************************************************************

  if (isMale)
  {
    A = 35.471;
    B = -8.472;
    C = -0.098;
  }
  else
  {
    A = 32.709;
    B = -7.378;
    C = -0.174;
  }
  PALWI = A + B*exp(C*t);
  PALWI/= 10.0;

  // ****************************************************************
  // SNHY:
  // Distance from the sella-nasion line to the midpoint of the 
  // vertical dimension of the hyoid bone measured perpendicular to
  // the Frankfort horizontal.
  // ****************************************************************

  if (isMale)
  {
    A = 62.731;
    B = -0.526;
    C = 1.103;
    D = 82.559;
    E = 2.655;
    F = 0.423;
  }
  else
  {
    A = 64.002;
    B = -0.621;
    C = 1.143;
    D = 38.254;
    E = 3.06;
    F = 0.597;
  }
  SNHY = A/(1+exp(B - C*C*t)) + D/(1+exp(E-F*F*t));
  SNHY/= 10.0;

  // ****************************************************************
  // SNPNS:
  // Distance from the sella-nasion line to the posterior nasal spine
  // measured perpendicular to the Frankfort horizontal.
  // ****************************************************************

  if (isMale)
  {
    A = 42.342;
    B = 0.169;
    C = 0.372;
    D = 6.91;
    E = 1.154;
    F = 1.748;
  }
  else
  {
    A = 22.865;
    B = -0.507;
    C = 1.326;
    D = 21.141;
    E = 1.35;
    F = 0.519;
  }
  SNPNS = A/(1+exp(B - C*C*t)) + D/(1+exp(E-F*F*t));
  SNPNS/= 10.0;

  // ****************************************************************
  // SYMPH:
  // Height of the mandibular symphysis as measured from the highest
  // and most prominent point on the lower alveolar arch to the 
  // lowest point on the symphysis.
  // ****************************************************************

  A = 23.411;
  B = -1.310;
  C = -0.626;
  D = 4.481;
  E = 23.061;
  F = 1.369;
  SYMPH = A/(1+exp(B - C*C*t)) + D/(1+exp(E-F*F*t));
  SYMPH/= 10.0;

  // ****************************************************************
  // THIHY:
  // Distance from lower border of hyoid bone to the glottis measured
  // along the dorsal surface of the epiglottis.
  // ****************************************************************

  A = 5.633;
  B = 0.185;
  THIHY = A + B*t;
  THIHY/= 10.0;

  // ****************************************************************
  // Combine the above date into the anatomic variables.
  // ****************************************************************

  param[LIP_WIDTH].x = LIPS - 0.5;
  param[HARD_PALATE_LENGTH].x = ANSPNS - 1.0; // Subtract the width of the upper incisor.
  param[SOFT_PALATE_LENGTH].x = ATLPTM - 0.2;

  // The 8.0 deg result from the angle between the maxillary alveolar
  // margin and the palatal line!
  param[ORAL_PHARYNGEAL_ANGLE].x = -(ANGFHV + ANGFHJ - ANGPJ + 90.0) + 8.0;

  if (param[ORAL_PHARYNGEAL_ANGLE].x > -90.0) 
  { 
    param[ORAL_PHARYNGEAL_ANGLE].x = -90.0; 
  }

  param[PHARYNX_LENGTH].x = SNHY - SNPNS - 0.5*THIHY - 0.2;
  param[LARYNX_LENGTH].x = HYGLOTT*sin(120.0*M_PI/180.0) + THIHY + 0.1;
  param[VOCAL_FOLD_LENGTH].x = LENVC;
  param[LARYNX_WIDTH].x = param[SOFT_PALATE_LENGTH].x;
  
  param[PALATE_HEIGHT].x = PALHIT;
  
  // Teeth height according to Goldstein (1980), p. 167

  if (age_month < 7*12)
  {
    param[UPPER_MOLARS_HEIGHT].x = 0.3;
    param[LOWER_MOLARS_HEIGHT].x = param[UPPER_MOLARS_HEIGHT].x;
  }
  else
  {
    param[UPPER_MOLARS_HEIGHT].x = 0.5;
    param[LOWER_MOLARS_HEIGHT].x = param[UPPER_MOLARS_HEIGHT].x;
  }

  param[MANDIBLE_HEIGHT].x = SYMPH*sin(70.0*M_PI/180.0) - 1.0;
  param[PALATE_DEPTH].x = PALWI;
}

// ****************************************************************************
// Calc. the parameter values from a given vocal tract.
// ****************************************************************************

void AnatomyParams::getFrom(VocalTract *tract)
{
  VocalTract::Anatomy *a = &tract->anatomy;

  param[LIP_WIDTH].x = a->lipsWidth_cm;
  param[MANDIBLE_HEIGHT].x = a->jawHeight_cm[3];
  param[LOWER_MOLARS_HEIGHT].x = a->lowerTeethHeight_cm[3];
  param[UPPER_MOLARS_HEIGHT].x = a->upperTeethHeight_cm[3];

  param[PALATE_HEIGHT].x = a->palateHeight_cm[3];
  param[PALATE_DEPTH].x  = -2.0*a->palatePoints[3].z;

  param[HARD_PALATE_LENGTH].x = a->palatePoints[8].x;
  param[SOFT_PALATE_LENGTH].x = -a->pharynxFulcrum.x;

  param[PHARYNX_LENGTH].x = a->palateHeight_cm[3] -
    0.5*(tract->param[VocalTract::HY].min + tract->param[VocalTract::HY].max);

  param[LARYNX_LENGTH].x = -0.5*(a->larynxNarrowPoints[3].y + a->larynxWidePoints[3].y);

  // This is the maximum width of the larynx.
  param[LARYNX_WIDTH].x = a->larynxWidePoints[0].x;

  param[VOCAL_FOLD_LENGTH].x = 
    0.5*(a->larynxNarrowPoints[3].x - a->larynxNarrowPoints[4].x +
         a->larynxWidePoints[3].x - a->larynxWidePoints[4].x);

  param[ORAL_PHARYNGEAL_ANGLE].x = a->pharynxRotationAngle_deg;
}

// ****************************************************************************
// Apply the current parameter values for the given vocal tract.
// ****************************************************************************

void AnatomyParams::setFor(VocalTract *tract)
{
  double widthScale, heightScale;
  double modelDepthScale;
  double upperTeethHeightScale, lowerTeethHeightScale;
  double d;
  int i;
  VocalTract::Anatomy *anatomy = &tract->anatomy;
  VocalTract::Anatomy *origAnatomy = &referenceVocalTract->anatomy;
  Point3D *P = NULL;    // Array of 3D points

  AnatomyParams origAnatomyParams;
  origAnatomyParams.getFrom(referenceVocalTract);

  // Restrict the parameters first.
  restrictParams();

  // ****************************************************************
  // Scale the lips.
  // ****************************************************************

  anatomy->lipsWidth_cm = param[LIP_WIDTH].x;

  // ****************************************************************
  // Scale the palate and mandible.
  // ****************************************************************

  // Length of the "inner palate" from the incisors to behind the 5th teeth.
  double innerDentitionLengthT5 = origAnatomy->palatePoints[8].x - origAnatomy->palatePoints[3].x;
  double innerPalateLength = param[HARD_PALATE_LENGTH].x;
  double innerJawLength    = param[HARD_PALATE_LENGTH].x * origAnatomy->jawPoints[8].x / origAnatomy->palatePoints[8].x;

  upperTeethHeightScale = param[UPPER_MOLARS_HEIGHT].x / origAnatomy->upperTeethHeight_cm[3];
  lowerTeethHeightScale = param[LOWER_MOLARS_HEIGHT].x / origAnatomy->lowerTeethHeight_cm[3];

  double posteriorPalateHeight = origAnatomy->palateHeight_cm[0];

  // ****************************************************************
  // The target palate length is shorter than the target milk teeth 
  // dentition.
  // ****************************************************************

  if (innerPalateLength < innerDentitionLengthT5)
  {
// posteriorPalateHeight = ...

    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // MUST BE REWORKED LATER !!!!
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    widthScale = innerPalateLength / (origAnatomy->palatePoints[8].x - origAnatomy->palatePoints[3].x);
    modelDepthScale = -param[PALATE_DEPTH].x / (origAnatomy->palatePoints[3].z*2.0);
    // Scale the milk teeth
    for (i=3; i < 9; i++)
    {
      anatomy->palatePoints[i].x = origAnatomy->palatePoints[0].x + 
        widthScale*(origAnatomy->palatePoints[i].x - origAnatomy->palatePoints[3].x);
      anatomy->palatePoints[i].z = modelDepthScale*origAnatomy->palatePoints[i].z;
      anatomy->palatePoints[i].y = origAnatomy->palatePoints[i].y;
    }

    // Make the first 3 ribs like the 4th rib and get rid of the molars!
    for (i=0; i < 3; i++)
    {
      anatomy->palatePoints[i] = anatomy->palatePoints[3];
      anatomy->palateAngle_deg[i] = anatomy->palateAngle_deg[3];
      anatomy->palateHeight_cm[i] = anatomy->palateHeight_cm[3];

      anatomy->upperTeethHeight_cm[i] = 0.0;
    }
  }
  else

  // ****************************************************************
  // The palate length is long enough for all milk teeth, and for
  // either one, two or three of the molars.
  // ****************************************************************

  {
    int missingTeeth = 0;

    // ****************************************************************
    // Palate.
    // ****************************************************************

    P = &origAnatomy->palatePoints[0];

//    widthScale = 1.0;
    heightScale = param[PALATE_HEIGHT].x / origAnatomy->palateHeight_cm[3];

    missingTeeth = 0;
    if (innerPalateLength < P[8].x - P[2].x) { missingTeeth = 3; } else
    if (innerPalateLength < P[8].x - P[1].x) { missingTeeth = 2; } else
    if (innerPalateLength < P[8].x - P[0].x) { missingTeeth = 1; }

    double innerDentitionLength = P[8].x - P[missingTeeth].x;
    double teethOffsetX = innerPalateLength - innerDentitionLength;

    modelDepthScale = -param[PALATE_DEPTH].x / (P[3].z*2.0);

    // ****************************************************************
    // Calculate the posterior palate height.
    // ****************************************************************

    posteriorPalateHeight = origAnatomy->palateHeight_cm[0];
    double x0, x1;
    for (i=0; i < VocalTract::NUM_JAW_RIBS-1; i++)
    {
      x0 = innerPalateLength - (P[8].x - P[i].x);
      x1 = innerPalateLength - (P[8].x - P[i+1].x);

      if ((x0 <= 0.0) && (x1 >= 0.0))
      {
        double den = x1 - x0;
        if (den < 0.000001) 
        { 
          den = 0.000001; 
        }
        double t = (0.0 - x0) / den;
        posteriorPalateHeight = (1.0-t)*origAnatomy->palateHeight_cm[i] + t*origAnatomy->palateHeight_cm[i+1];
      }
    }
    posteriorPalateHeight*= heightScale;

    // ****************************************************************
    // Set all palate parameters.
    // ****************************************************************

    for (i=8; i >= 0; i--)
    {
      if (i >= missingTeeth)
      {
        anatomy->palatePoints[i].x = teethOffsetX + (P[i].x - P[missingTeeth].x);
        anatomy->palatePoints[i].z = modelDepthScale*P[i].z;
        anatomy->palatePoints[i].y = P[i].y;

        anatomy->palateAngle_deg[i] = origAnatomy->palateAngle_deg[i];
        anatomy->palateHeight_cm[i] = origAnatomy->palateHeight_cm[i]*heightScale;
        anatomy->upperTeethHeight_cm[i] = origAnatomy->upperTeethHeight_cm[i]*upperTeethHeightScale;
      }
      else

      // Put the ribs without existing teeth at the place x=0 and make
      // them look like the rib for which the first tooth exists.
      // These ribs are changed due to velum movement during articulation 
      // anyway. And set the height of the non-existing teeth to zero!
      {
        anatomy->palatePoints[i] = anatomy->palatePoints[missingTeeth];
        anatomy->palatePoints[i].x = 0.000001;

        anatomy->palateAngle_deg[i] = anatomy->palateAngle_deg[missingTeeth];
        anatomy->palateHeight_cm[i] = posteriorPalateHeight;
        anatomy->upperTeethHeight_cm[i] = 0.0;
      }

      anatomy->upperTeethWidthBottom_cm[i] = origAnatomy->upperTeethWidthBottom_cm[i]*modelDepthScale;
      anatomy->upperTeethWidthTop_cm[i]    = origAnatomy->upperTeethWidthTop_cm[i]*modelDepthScale;
    }

    // ****************************************************************
    // Mandible.
    // ****************************************************************

    P = &origAnatomy->jawPoints[0];
    
    // Let modelDepthScale be the same as before!
    // But the height scale is different:
    heightScale = param[MANDIBLE_HEIGHT].x / origAnatomy->jawHeight_cm[3];

    missingTeeth = 0;
    if (innerJawLength < P[8].x - P[2].x) { missingTeeth = 3; } else
    if (innerJawLength < P[8].x - P[1].x) { missingTeeth = 2; } else
    if (innerJawLength < P[8].x - P[0].x) { missingTeeth = 1; }

    innerDentitionLength = P[8].x - P[missingTeeth].x;
    teethOffsetX = innerPalateLength - innerDentitionLength;

    // ****************************************************************
    // Set all jaw parameters.
    // ****************************************************************

    for (i=8; i >= 0; i--)
    {
      if (i >= missingTeeth)
      {
        anatomy->jawPoints[i].x = teethOffsetX + (P[i].x - P[missingTeeth].x);
        anatomy->jawPoints[i].z = modelDepthScale*P[i].z;
        anatomy->jawPoints[i].y = P[i].y;

        anatomy->jawAngle_deg[i] = origAnatomy->jawAngle_deg[i];
        anatomy->jawHeight_cm[i] = origAnatomy->jawHeight_cm[i]*heightScale;
        anatomy->lowerTeethHeight_cm[i] = origAnatomy->lowerTeethHeight_cm[i]*lowerTeethHeightScale;
      }
      else

      // Put the ribs without existing teeth at the place x=0 and make
      // them look like the rib for which the first tooth exists.
      // These ribs are changed due to velum movement during articulation 
      // anyway. And set the height of the non-existing teeth to zero!
      {
        anatomy->jawPoints[i] = anatomy->jawPoints[missingTeeth];
        anatomy->jawPoints[i].x = 0.000001;

        anatomy->jawAngle_deg[i] = anatomy->jawAngle_deg[missingTeeth];
        anatomy->jawHeight_cm[i] = anatomy->jawHeight_cm[missingTeeth];
        anatomy->lowerTeethHeight_cm[i] = 0.0;
      }

      anatomy->lowerTeethWidthBottom_cm[i] = origAnatomy->lowerTeethWidthBottom_cm[i]*modelDepthScale;
      anatomy->lowerTeethWidthTop_cm[i]    = origAnatomy->lowerTeethWidthTop_cm[i]*modelDepthScale;
    }
  }

  // ****************************************************************

  // Length of the lower incisor roots.

  heightScale = param[MANDIBLE_HEIGHT].x / origAnatomyParams.param[MANDIBLE_HEIGHT].x;
  anatomy->toothRootLength_cm = origAnatomy->toothRootLength_cm * heightScale;

  // Move the vertical jaw rest position when the height of 
  // the teeth changed.

  double deltaUpperTeethHeight_cm = param[UPPER_MOLARS_HEIGHT].x - origAnatomy->upperTeethHeight_cm[3];
  double deltaLowerTeethHeight_cm = param[LOWER_MOLARS_HEIGHT].x - origAnatomy->lowerTeethHeight_cm[3];
  
  anatomy->jawFulcrum = origAnatomy->jawFulcrum;
  anatomy->jawRestPos = origAnatomy->jawRestPos;
  anatomy->jawRestPos.y-= deltaUpperTeethHeight_cm + deltaLowerTeethHeight_cm;


  // ****************************************************************
  // Scale the velum points in all three states and the center of
  // rotation (between mouth and pharynx). The lowest point of the
  // lowest velum state is kept fixed! Scaling is only done upwards!
  // ****************************************************************

  double origLowestVelumPointY  = origAnatomy->velumLowPoints[0].y;
  double origHighestVelumPointY = origAnatomy->velumHighPoints[VocalTract::NUM_VELUM_RIBS-2].y;
  double newHighestVelumPointY  = posteriorPalateHeight; //anatomy->palateHeight_cm[1];

  heightScale = (newHighestVelumPointY - origLowestVelumPointY) / (origHighestVelumPointY - origLowestVelumPointY);
  widthScale = param[SOFT_PALATE_LENGTH].x / (-origAnatomy->pharynxFulcrum.x);

  anatomy->pharynxFulcrum.x = origAnatomy->pharynxFulcrum.x*widthScale;
  anatomy->pharynxFulcrum.y = origLowestVelumPointY + 
    (origAnatomy->pharynxFulcrum.y - origLowestVelumPointY)*heightScale;

  for (i=0; i < VocalTract::NUM_VELUM_RIBS-1; i++)
  {
    anatomy->velumLowPoints[i].x = origAnatomy->velumLowPoints[i].x*widthScale;
    anatomy->velumMidPoints[i].x = origAnatomy->velumMidPoints[i].x*widthScale;
    anatomy->velumHighPoints[i].x = origAnatomy->velumHighPoints[i].x*widthScale;

    anatomy->velumLowPoints[i].y = origLowestVelumPointY + 
      (origAnatomy->velumLowPoints[i].y - origLowestVelumPointY)*heightScale;
    anatomy->velumMidPoints[i].y = origLowestVelumPointY + 
      (origAnatomy->velumMidPoints[i].y - origLowestVelumPointY)*heightScale;
    anatomy->velumHighPoints[i].y = origLowestVelumPointY + 
      (origAnatomy->velumHighPoints[i].y - origLowestVelumPointY)*heightScale;
  }

  // Scale the maximum nasal port area.

  d = (param[PALATE_DEPTH].x *param[SOFT_PALATE_LENGTH].x) / 
    (origAnatomyParams.param[PALATE_DEPTH].x * origAnatomyParams.param[SOFT_PALATE_LENGTH].x);

  anatomy->maxNasalPortArea_cm2 = origAnatomy->maxNasalPortArea_cm2 * d;


  // ****************************************************************
  // Scale the depth of the entire vocal tract with modelDepthScale,
  // that was used for depth scaling of the jaws.
  // ****************************************************************

  anatomy->uvulaDepth_cm = origAnatomy->uvulaDepth_cm*modelDepthScale;
  anatomy->epiglottisDepth_cm = origAnatomy->epiglottisDepth_cm*modelDepthScale;
  anatomy->pharynxUpperDepth_cm = origAnatomy->pharynxUpperDepth_cm*modelDepthScale;
  anatomy->pharynxLowerDepth_cm = origAnatomy->pharynxLowerDepth_cm*modelDepthScale;
  anatomy->larynxUpperDepth_cm = origAnatomy->larynxUpperDepth_cm*modelDepthScale;
  anatomy->larynxLowerDepth_cm = origAnatomy->larynxLowerDepth_cm*modelDepthScale;
  
  // ****************************************************************
  // Calc. the shape of the larynx.
  // ****************************************************************

  double meanHeight = -0.5*(origAnatomy->larynxNarrowPoints[3].y + origAnatomy->larynxWidePoints[3].y);
  heightScale = param[LARYNX_LENGTH].x / meanHeight;
  
  widthScale  = param[LARYNX_WIDTH].x / origAnatomy->larynxWidePoints[0].x;
  
  for (i=0; i < 8; i++)
  {
    anatomy->larynxNarrowPoints[i].x = origAnatomy->larynxNarrowPoints[i].x*widthScale;
    anatomy->larynxNarrowPoints[i].y = origAnatomy->larynxNarrowPoints[i].y*heightScale;

    anatomy->larynxWidePoints[i].x = origAnatomy->larynxWidePoints[i].x*widthScale;
    anatomy->larynxWidePoints[i].y = origAnatomy->larynxWidePoints[i].y*heightScale;
  }

  // Scale the lower part of the larynx according to the new vocal folds length

  widthScale = param[VOCAL_FOLD_LENGTH].x / (origAnatomy->larynxNarrowPoints[3].x - origAnatomy->larynxNarrowPoints[4].x);
  anatomy->larynxNarrowPoints[3].x = anatomy->larynxNarrowPoints[4].x + param[VOCAL_FOLD_LENGTH].x;
  anatomy->larynxNarrowPoints[2].x = anatomy->larynxNarrowPoints[5].x + 
    widthScale*(origAnatomy->larynxNarrowPoints[2].x - origAnatomy->larynxNarrowPoints[5].x);

  widthScale = param[VOCAL_FOLD_LENGTH].x / (origAnatomy->larynxWidePoints[3].x - origAnatomy->larynxWidePoints[4].x);
  anatomy->larynxWidePoints[3].x = anatomy->larynxWidePoints[4].x + param[VOCAL_FOLD_LENGTH].x;
  anatomy->larynxWidePoints[2].x = anatomy->larynxWidePoints[5].x + 
    widthScale*(origAnatomy->larynxWidePoints[2].x - origAnatomy->larynxWidePoints[5].x);

  // Make sure that the larynx tube is not too constricted in the
  // transformed larynx. The constriction size should not be smaller 
  // than in the reference VT.

  double minWidth = origAnatomy->larynxNarrowPoints[2].x - origAnatomy->larynxNarrowPoints[5].x;
  if (anatomy->larynxNarrowPoints[2].x - anatomy->larynxNarrowPoints[5].x < minWidth)
  {
    anatomy->larynxNarrowPoints[5].x = anatomy->larynxNarrowPoints[2].x - minWidth;
  }

  minWidth = origAnatomy->larynxWidePoints[2].x - origAnatomy->larynxWidePoints[5].x;
  if (anatomy->larynxWidePoints[2].x - anatomy->larynxWidePoints[5].x < minWidth)
  {
    anatomy->larynxWidePoints[5].x = anatomy->larynxWidePoints[2].x - minWidth;
  }

  // ****************************************************************
  // Set the angle between the mouth and the pharynx cavity.
  // ****************************************************************

  anatomy->pharynxRotationAngle_deg = param[ORAL_PHARYNGEAL_ANGLE].x;

  // ****************************************************************
  // Scale the epiglottis and uvula according to the pharynx length.
  // ****************************************************************

  heightScale = param[PHARYNX_LENGTH].x / origAnatomyParams.param[PHARYNX_LENGTH].x;

  anatomy->uvulaWidth_cm  = origAnatomy->uvulaWidth_cm*heightScale;
  anatomy->uvulaHeight_cm = origAnatomy->uvulaHeight_cm*heightScale;
  anatomy->uvulaDepth_cm  = origAnatomy->uvulaDepth_cm*modelDepthScale;

  anatomy->epiglottisWidth_cm  = origAnatomy->epiglottisWidth_cm*heightScale;
  anatomy->epiglottisHeight_cm = origAnatomy->epiglottisHeight_cm*heightScale;
  anatomy->epiglottisDepth_cm  = origAnatomy->epiglottisDepth_cm*modelDepthScale;

  // ****************************************************************
  // Scale the size of the tongue body and tongue tip circles.
  // ****************************************************************

  widthScale  = (param[HARD_PALATE_LENGTH].x + param[SOFT_PALATE_LENGTH].x) / 
    (origAnatomyParams.param[HARD_PALATE_LENGTH].x + origAnatomyParams.param[SOFT_PALATE_LENGTH].x);

  heightScale = param[PHARYNX_LENGTH].x / origAnatomyParams.param[PHARYNX_LENGTH].x;

  anatomy->tongueCenterRadiusX_cm = origAnatomy->tongueCenterRadiusX_cm * widthScale;
  anatomy->tongueCenterRadiusY_cm = origAnatomy->tongueCenterRadiusY_cm * heightScale;
  anatomy->tongueTipRadius_cm     = origAnatomy->tongueTipRadius_cm * widthScale;

  // The anatomy parameters of the tongue root are adjusted
  // PRECISELY in adjustTongueRootCalculation() further below.
  // Here, they are just set to preliminary working values.

  anatomy->automaticTongueRootCalc = origAnatomy->automaticTongueRootCalc;
  anatomy->tongueRootTrxSlope     = origAnatomy->tongueRootTrxSlope;
  anatomy->tongueRootTrxIntercept = origAnatomy->tongueRootTrxIntercept;
  anatomy->tongueRootTrySlope     = origAnatomy->tongueRootTrySlope;
  anatomy->tongueRootTryIntercept = origAnatomy->tongueRootTryIntercept;

  // ****************************************************************
  // Change the min/max-values of the vocal tract parameters.
  // ****************************************************************

  heightScale = param[PHARYNX_LENGTH].x / origAnatomyParams.param[PHARYNX_LENGTH].x;

  widthScale  = (param[HARD_PALATE_LENGTH].x + param[SOFT_PALATE_LENGTH].x) / 
    (origAnatomyParams.param[HARD_PALATE_LENGTH].x + origAnatomyParams.param[SOFT_PALATE_LENGTH].x);

  VocalTract::Param *origParam = referenceVocalTract->param;
  VocalTract::Param *newParam  = tract->param;
  
  // HX is a relative value (0.0 ... 1.0) and needs no scaling.

  // HY.

  newParam[VocalTract::HY].min = transformY(&origAnatomyParams, origParam[VocalTract::HY].min);
  newParam[VocalTract::HY].max = transformY(&origAnatomyParams, origParam[VocalTract::HY].max);
  // HY may not go above the lower end of the mandible...
  if (newParam[VocalTract::HY].min > d) { newParam[VocalTract::HY].min = d; }
  if (newParam[VocalTract::HY].max > d) { newParam[VocalTract::HY].max = d; }

  // JX is a relative displacement value.
  
  newParam[VocalTract::JX].min = origParam[VocalTract::JX].min * widthScale;
  newParam[VocalTract::JX].max = origParam[VocalTract::JX].max * widthScale;

  // JA is an angle and is invariant to the length transformations.

  // LP is a relative value and needs no transformation.

  // LD.

  newParam[VocalTract::LD].min = origParam[VocalTract::LD].min*heightScale;
  newParam[VocalTract::LD].max = origParam[VocalTract::LD].max*heightScale;

  // VS and VO are relative values and need no transformation.

  // TCX and TCY.

  newParam[VocalTract::TCX].min = transformX(&origAnatomyParams, origParam[VocalTract::TCX].min);
  newParam[VocalTract::TCX].max = transformX(&origAnatomyParams, origParam[VocalTract::TCX].max);

  newParam[VocalTract::TCY].min = transformY(&origAnatomyParams, origParam[VocalTract::TCY].min);
  newParam[VocalTract::TCY].max = transformY(&origAnatomyParams, origParam[VocalTract::TCY].max);

  // TTX and TTY.

  newParam[VocalTract::TTX].min = transformX(&origAnatomyParams, origParam[VocalTract::TTX].min);
  newParam[VocalTract::TTX].max = transformX(&origAnatomyParams, origParam[VocalTract::TTX].max);

  newParam[VocalTract::TTY].min = transformY(&origAnatomyParams, origParam[VocalTract::TTY].min);
  newParam[VocalTract::TTY].max = transformY(&origAnatomyParams, origParam[VocalTract::TTY].max);

  // TBX and TBY.

  newParam[VocalTract::TBX].min = transformX(&origAnatomyParams, origParam[VocalTract::TBX].min);
  newParam[VocalTract::TBX].max = transformX(&origAnatomyParams, origParam[VocalTract::TBX].max);

  newParam[VocalTract::TBY].min = transformY(&origAnatomyParams, origParam[VocalTract::TBY].min);
  newParam[VocalTract::TBY].max = transformY(&origAnatomyParams, origParam[VocalTract::TBY].max);

  // TRX and TRY.

  newParam[VocalTract::TRX].min = transformX(&origAnatomyParams, origParam[VocalTract::TRX].min);
  newParam[VocalTract::TRX].max = transformX(&origAnatomyParams, origParam[VocalTract::TRX].max);

  newParam[VocalTract::TRY].min = transformY(&origAnatomyParams, origParam[VocalTract::TRY].min);
  newParam[VocalTract::TRY].max = transformY(&origAnatomyParams, origParam[VocalTract::TRY].max);

  // TS1, TS2, TS3, TS4. Scale the tongue side height linearly with the VT depth.

  double depthScale = param[PALATE_DEPTH].x / origAnatomyParams.param[PALATE_DEPTH].x;
  for (i=0; i < 4; i++)
  {
    newParam[VocalTract::TS1+i].min = origParam[VocalTract::TS1+i].min*depthScale;
    newParam[VocalTract::TS1+i].max = origParam[VocalTract::TS1+i].max*depthScale;
  }

  // Leave MA1, MA2, MA3 (the minimum area ranges) as they are ...


  // ****************************************************************
  // Scale the actual vocal tract parameter values.
  // ****************************************************************

  double oldTractParams[VocalTract::NUM_PARAMS];
  double newTractParams[VocalTract::NUM_PARAMS];

  // The actual parameter value
  
  for (i=0; i < VocalTract::NUM_PARAMS; i++) 
  { 
    oldTractParams[i] = origParam[i].x; 
  }
  
  adaptArticulation(oldTractParams, newTractParams);

  for (i=0; i < VocalTract::NUM_PARAMS; i++) 
  { 
    newParam[i].x = newTractParams[i]; 
  }

  // ****************************************************************
  // Scale the neutral parameter value.
  // ****************************************************************

  for (i=0; i < VocalTract::NUM_PARAMS; i++) 
  { 
    oldTractParams[i] = origParam[i].neutral; 
  }

  adaptArticulation(oldTractParams, newTractParams);

  for (i=0; i < VocalTract::NUM_PARAMS; i++) 
  { 
    newParam[i].neutral = newTractParams[i]; 
  }

  // ****************************************************************
  // Transform all the vocal tract shapes of the reference speaker.
  // ****************************************************************

  tract->shapes.clear();
  VocalTract::Shape shape;
  int k;

  for (i=0; i < (int)referenceVocalTract->shapes.size(); i++)
  {
    adaptArticulation(referenceVocalTract->shapes[i].param, newTractParams);

    shape.name = referenceVocalTract->shapes[i].name;
    for (k=0; k < VocalTract::NUM_PARAMS; k++)
    {
      shape.param[k] = newTractParams[k];
    }

    tract->shapes.push_back(shape);
  }

  // ****************************************************************
  // Reinitialize the vocal tract with the modified anatomy.
  // ****************************************************************

  tract->initReferenceSurfaces();
  tract->calculateAll();

  adjustTongueRootCalculation(tract);
  tract->calculateAll();
}


// ****************************************************************************
/// Adjust the slopes and intercepts of the linear equations to predict the
/// tongue root parameters (TRX, TRY) from the hyoid and tongue body positions.
// ****************************************************************************

void AnatomyParams::adjustTongueRootCalculation(VocalTract *targetVocalTract)
{
  // We consider two extreme shapes (1 and 2). Shape 1 is similar
  // to /a/ with the tongue body sitting directly on the hyoid.
  // Shape 2 is similar to /i/ with the tongue body directly below
  // the middle point of the hard palate.
  // The predicted TRX, TRY values for these two shapes in the 
  // original model are transformed to corresponding TRX and 
  // TRY values in the new model.
  
  Point2D origH, origTc1, origTc2, origTr1, origTr2;
  Point2D newH, newTc1, newTc2, newTr1, newTr2;
  double distance;
  double d1, d2;
  double tx, ty;

  int i;

  // ****************************************************************
  // Save the original parameters to reset them at the end.
  // ****************************************************************

  double origSavedParams[VocalTract::NUM_PARAMS];
  double newSavedParams[VocalTract::NUM_PARAMS];

  for (i=0; i < VocalTract::NUM_PARAMS; i++)
  {
    origSavedParams[i] = referenceVocalTract->param[i].x;
    newSavedParams[i] = targetVocalTract->param[i].x;
  }

  // ****************************************************************

  VocalTract::Anatomy *newAnatomy = &targetVocalTract->anatomy;
  VocalTract::Anatomy *origAnatomy = &referenceVocalTract->anatomy;

  // Set all parameters to average values, especially HX and HY!
  
  for (i=0; i < VocalTract::NUM_PARAMS; i++)
  {
    referenceVocalTract->param[i].x = 0.5*(referenceVocalTract->param[i].min + referenceVocalTract->param[i].max);
    targetVocalTract->param[i].x = 0.5*(targetVocalTract->param[i].min + targetVocalTract->param[i].max);
  }
  referenceVocalTract->calculateAll();
  targetVocalTract->calculateAll();
  
  // Hyoid position.
  origH = referenceVocalTract->surface[VocalTract::LOWER_COVER].getVertex(VocalTract::NUM_LARYNX_RIBS-1, VocalTract::NUM_LOWER_COVER_POINTS-1).toPoint2D();
  newH = targetVocalTract->surface[VocalTract::LOWER_COVER].getVertex(VocalTract::NUM_LARYNX_RIBS-1, VocalTract::NUM_LOWER_COVER_POINTS-1).toPoint2D();

  // Tongue center position for the /a/ shape directly above the hyoid.
  origTc1.set( origH.x, origH.y + origAnatomy->tongueCenterRadiusY_cm );
  newTc1.set( newH.x, newH.y + newAnatomy->tongueCenterRadiusY_cm );

  // Tongue center position for the /a/ shape directly above the hyoid.
  origTc2.set( origAnatomy->palatePoints[3].x, origAnatomy->palateHeight_cm[3] - origAnatomy->tongueCenterRadiusY_cm );
  newTc2.set( newAnatomy->palatePoints[3].x, newAnatomy->palateHeight_cm[3] - newAnatomy->tongueCenterRadiusY_cm );

  // Calculate the tongue root positions for both shapes in the
  // original anatomy.
  
  distance = (origH - origTc1).magnitude();
  origTr1.x = origAnatomy->tongueRootTrxSlope * distance + origAnatomy->tongueRootTrxIntercept;
  origTr1.y = origAnatomy->tongueRootTrySlope * origTc1.x + origAnatomy->tongueRootTryIntercept;

  distance = (origH - origTc2).magnitude();
  origTr2.x = origAnatomy->tongueRootTrxSlope * distance + origAnatomy->tongueRootTrxIntercept;
  origTr2.y = origAnatomy->tongueRootTrySlope * origTc2.x + origAnatomy->tongueRootTryIntercept;

  // Where are the tongue root positions in the new shapes ?
  // Assume that they have the same relative positions to the
  // tongue body as in the original vocal tract.

  tx = (origTr1.x - origTc1.x) / origAnatomy->tongueCenterRadiusX_cm;
  ty = (origTr1.y - origTc1.y) / origAnatomy->tongueCenterRadiusY_cm;
  newTr1.x = newTc1.x + newAnatomy->tongueCenterRadiusX_cm * tx;
  newTr1.y = newTc1.y + newAnatomy->tongueCenterRadiusY_cm * ty;

  tx = (origTr2.x - origTc2.x) / origAnatomy->tongueCenterRadiusX_cm;
  ty = (origTr2.y - origTc2.y) / origAnatomy->tongueCenterRadiusY_cm;
  newTr2.x = newTc2.x + newAnatomy->tongueCenterRadiusX_cm * tx;
  newTr2.y = newTc2.y + newAnatomy->tongueCenterRadiusY_cm * ty;

  // ****************************************************************
  // We have 2 equations for the tongue root x-position:
  // Tr1.x = slopeX*d1 + interceptX
  // Tr2.x = slopeX*d2 + interceptX
  //
  // These can be solved for slopeX and interceptX:
  // slopeX = (Tr2.x - Tr1.x) / (d2 - d1)
  // interceptX = Tr1.x - slopeX*d1
  // ****************************************************************

  d1 = (newH - newTc1).magnitude();
  d2 = (newH - newTc2).magnitude();
  newAnatomy->tongueRootTrxSlope     = (newTr2.x - newTr1.x) / (d2 - d1);
  newAnatomy->tongueRootTrxIntercept = newTr1.x - newAnatomy->tongueRootTrxSlope * d1;

  // ****************************************************************
  // We have 2 equations for the tongue root y-position:
  // Tr1.y = slopeY*Tc1.x + interceptY
  // Tr2.y = slopeY*Tc2.x + interceptY
  //
  // These can be solved for slopeY and interceptY:
  // slopeY = (Tr2.y - Tr1.y) / (Tc2.x - Tc1.x)
  // interceptY = Tr1.y - slopeY*Tc1.x
  // ****************************************************************

  newAnatomy->tongueRootTrySlope     = (newTr2.y - newTr1.y) / (newTc2.x - newTc1.x);
  newAnatomy->tongueRootTryIntercept = newTr1.y - newAnatomy->tongueRootTrySlope * newTc1.x;

  newAnatomy->automaticTongueRootCalc = origAnatomy->automaticTongueRootCalc;


  // ****************************************************************
  // Reset the previous VT parameters.
  // ****************************************************************

  for (i=0; i < VocalTract::NUM_PARAMS; i++)
  {
    referenceVocalTract->param[i].x = origSavedParams[i];
    targetVocalTract->param[i].x = newSavedParams[i];
  }

  referenceVocalTract->calculateAll();
  targetVocalTract->calculateAll();
}


// ****************************************************************************
/// Loads the given vocal tract model to use it as reference speaker for the
/// transformations. It should be a male adult vocal tract.
// ****************************************************************************

bool AnatomyParams::loadReferenceVocalTract(const string &fileName)
{
  bool ok = true;

  try
  {
    referenceVocalTract->readFromXml(fileName);
    referenceVocalTract->calculateAll();
  }
  catch (std::string st)
  {
    printf("Error reading the anatomy data from %s.\n", fileName.c_str());
    ok = false;
  }

  return ok;
}


// ****************************************************************************
// Transforms the x-coord. oldX from the anatomy given by origParams to the 
// current anatomy by linear scalings.
// ****************************************************************************

double AnatomyParams::transformX(AnatomyParams *origAnatomy, double origX)
{
  double newX;

  if (origX < 0.0)
  {
    // Horizontal scaling according to the length ratios of the velum
    newX = origX * param[SOFT_PALATE_LENGTH].x / origAnatomy->param[SOFT_PALATE_LENGTH].x;
  }
  else
  {
    // Horizontal scaling according to the length ratios of the palate
    newX = origX * param[HARD_PALATE_LENGTH].x / origAnatomy->param[HARD_PALATE_LENGTH].x;
  }

  return newX;
}

// ****************************************************************************
// Transforms the y-coord. oldY from the anatomy given by origParams to the 
// current anatomy by linear scalings.
// ****************************************************************************

double AnatomyParams::transformY(AnatomyParams *origAnatomy, double origX, double origY)
{
  double t = (origX - (-origAnatomy->param[SOFT_PALATE_LENGTH].x)) / 
    (origAnatomy->param[SOFT_PALATE_LENGTH].x + origAnatomy->param[HARD_PALATE_LENGTH].x);
  if (t < 0.0) { t = 0.0; }
  if (t > 1.0) { t = 1.0; }

  double newLeftY = param[PALATE_HEIGHT].x + 
    (origY - origAnatomy->param[PALATE_HEIGHT].x) * param[PHARYNX_LENGTH].x / origAnatomy->param[PHARYNX_LENGTH].x;

  double frontHeight = param[PALATE_HEIGHT].x + param[UPPER_MOLARS_HEIGHT].x + 
    param[LOWER_MOLARS_HEIGHT].x + param[MANDIBLE_HEIGHT].x;
  double origFrontHeight = origAnatomy->param[PALATE_HEIGHT].x + origAnatomy->param[UPPER_MOLARS_HEIGHT].x + 
    origAnatomy->param[LOWER_MOLARS_HEIGHT].x + origAnatomy->param[MANDIBLE_HEIGHT].x;

  double newRightY = param[PALATE_HEIGHT].x + 
    (origY - origAnatomy->param[PALATE_HEIGHT].x) * frontHeight / origFrontHeight;

  return (1.0-t)*newLeftY + t*newRightY;
}

// ****************************************************************************
// Transforms the y-coord. oldY from the anatomy given by origParams to the 
// current anatomy by linear scalings.
// ****************************************************************************

double AnatomyParams::transformY(AnatomyParams *origAnatomy, double origY)
{
  double newY = param[PALATE_HEIGHT].x + 
    (origY - origAnatomy->param[PALATE_HEIGHT].x) * param[PHARYNX_LENGTH].x / origAnatomy->param[PHARYNX_LENGTH].x;
  return newY;
}

// ****************************************************************************
// Adapt the articulation of a shape from the original vocal tract to the
// current anatomy parameters.
// ****************************************************************************

void AnatomyParams::adaptArticulation(double *oldParams, double *newParams)
{
  int i;
  AnatomyParams origAnatomyParams;
  origAnatomyParams.getFrom(referenceVocalTract);
  
  double heightScale = param[PHARYNX_LENGTH].x / origAnatomyParams.param[PHARYNX_LENGTH].x;
  double widthScale = (param[SOFT_PALATE_LENGTH].x + param[HARD_PALATE_LENGTH].x) / 
    (origAnatomyParams.param[SOFT_PALATE_LENGTH].x + origAnatomyParams.param[HARD_PALATE_LENGTH].x);

  newParams[VocalTract::HX] = oldParams[VocalTract::HX];
  newParams[VocalTract::HY] = transformY(&origAnatomyParams, 0.0, oldParams[VocalTract::HY]);

  newParams[VocalTract::JX] = oldParams[VocalTract::JX] * widthScale;
  newParams[VocalTract::JA] = oldParams[VocalTract::JA];

  newParams[VocalTract::LP] = oldParams[VocalTract::LP];
  newParams[VocalTract::LD] = oldParams[VocalTract::LD] * heightScale;

  newParams[VocalTract::VS] = oldParams[VocalTract::VS];
  newParams[VocalTract::VO] = oldParams[VocalTract::VO];

  newParams[VocalTract::TCX] = transformX(&origAnatomyParams, oldParams[VocalTract::TCX]);
  newParams[VocalTract::TCY] = transformY(&origAnatomyParams, oldParams[VocalTract::TCX], oldParams[VocalTract::TCY]);
  
  newParams[VocalTract::TTX] = transformX(&origAnatomyParams, oldParams[VocalTract::TTX]);
  newParams[VocalTract::TTY] = transformY(&origAnatomyParams, oldParams[VocalTract::TTX], oldParams[VocalTract::TTY]);

  newParams[VocalTract::TBX] = transformX(&origAnatomyParams, oldParams[VocalTract::TBX]);
  newParams[VocalTract::TBY] = transformY(&origAnatomyParams, oldParams[VocalTract::TBX], oldParams[VocalTract::TBY]);

  newParams[VocalTract::TRX] = transformX(&origAnatomyParams, oldParams[VocalTract::TRX]);
  newParams[VocalTract::TRY] = transformY(&origAnatomyParams, oldParams[VocalTract::TRX], oldParams[VocalTract::TRY]);

  // Scale the tongue side height linearly with the VT depth. *******

  double depthScale = param[PALATE_DEPTH].x / origAnatomyParams.param[PALATE_DEPTH].x;
  for (i=0; i < 4; i++)
  {
    newParams[VocalTract::TS1+i] = oldParams[VocalTract::TS1+i]*depthScale;
  }
}

// ****************************************************************************
