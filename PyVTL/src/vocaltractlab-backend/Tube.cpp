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

#include "Tube.h"
#include <iostream>
#include <cstdio>
#include <cmath>


// Data for the relaxed cheeck
//const double Tube::STANDARD_WALL_MASS_CGS       = 2.1;
//const double Tube::STANDARD_WALL_RESISTANCE_CGS = 800.0;
//const double Tube::STANDARD_WALL_STIFFNESS_CGS  = 84500.0;

// Data for the tense-cheeck
//const double Tube::STANDARD_WALL_MASS_CGS = 1.5;
//const double Tube::STANDARD_WALL_RESISTANCE_CGS = 1060.0;
//const double Tube::STANDARD_WALL_STIFFNESS_CGS = 33300.0;

// Data for the neck
//const double Tube::STANDARD_WALL_MASS_CGS = 2.4;
//const double Tube::STANDARD_WALL_RESISTANCE_CGS = 2320.0;
//const double Tube::STANDARD_WALL_STIFFNESS_CGS = 491000.0;


// I gave the wall a higher resistance than proposed in the
// literature. With a lower resistance, F1 increases very 
// strongly (making it hard to get, e.g., a good /u/).
// The wall stiffness mainly controls how long a voicebar
// lasts in the closure phase of a plosive. It has little
// effect on the F1 damping.

const double Tube::STANDARD_WALL_MASS_CGS = 2.4;
const double Tube::STANDARD_WALL_RESISTANCE_CGS = 5000.0;
const double Tube::STANDARD_WALL_STIFFNESS_CGS = 100000.0;


const double Tube::MIN_AREA_CM2 = 0.01e-2;  // 0.01 mm^2

// This is the default strength of the aspiration source with respect 
// to the maximum strength (=0 dB).
const double Tube::DEFAULT_ASPIRATION_STRENGTH_DB = -40.0;


// ****************************************************************
// Subglottal sections.
// We take a tube of 23 cm for the subglottal system, which has
// a rather constant cross-section along this part according to 
// Weibels data (1963), and which gives us F1_sub = 550 Hz and
// F2_sub = 1270 Hz.
// The large-scale study of Lulich et al. (2012) found subglottal
// resonances of men at 554, 1327, and 2179 Hz.
// ****************************************************************

const double Tube::DEFAULT_SUBGLOTTAL_CAVITY_LENGTH_CM = 23.0;
const double Tube::DEFAULT_NASAL_CAVITY_LENGTH_CM = 11.4;
const double Tube::DEFAULT_PIRIFORM_FOSSA_LENGTH_CM = 3.0;
const double Tube::DEFAULT_PIRIFORM_FOSSA_VOLUME_CM3 = 2.0;

// Index of the nose sections, to which the paranasal sinuses are attached
const int Tube::SINUS_COUPLING_SECTION[Tube::NUM_SINUS_SECTIONS] =
{ 
  8, 9, 11, 12 
};

const char Tube::ARTICULATOR_LETTER[Tube::NUM_ARTICULATORS] = { 'V', 'T', 'I', 'L', 'N' };

// ****************************************************************************
/// Constructor. Initializes the area function.
// ****************************************************************************

Tube::Tube()
{
  subglottalCavityLength_cm = 0.0;
  nasalCavityLength_cm = 0.0;
  piriformFossaLength_cm = 0.0;
  piriformFossaVolume_cm3 = 0.0;

  // Init the static parts.
  staticPartsInitialized = false;

  initSubglottalCavity();
  initNasalCavity();
  initPiriformFossa();

  staticPartsInitialized = true;

  resetDynamicPart();
  createSectionList();
  
  setVelumOpening(0.0);
  setGlottisArea(0.0);
  setAspirationStrength(DEFAULT_ASPIRATION_STRENGTH_DB);

  teethPosition_cm = 15.0;
  tongueTipSideElevation = 0.0;

  calcPositions();
}


// ****************************************************************************
/// Init the subglottal system.
// ****************************************************************************

void Tube::initSubglottalCavity(const double length_cm)
{
  int i;
  Section *ts = NULL;

  const double EPSILON = 0.000000001;
  if (fabs(length_cm - subglottalCavityLength_cm) < EPSILON)
  {
    return;
  }

  subglottalCavityLength_cm = length_cm;

  // ****************************************************************

  for (i=0; i < NUM_TRACHEA_SECTIONS; i++)
  {
    ts = &tracheaSection[i];

    ts->pos_cm     = 0.0;
    ts->area_cm2   = 2.5;     // 2.5 cm^2 trachea area
    ts->length_cm  = length_cm / (double)NUM_TRACHEA_SECTIONS;
    ts->volume_cm3 = ts->area_cm2 * ts->length_cm;

    ts->wallMass_cgs       = 0.25;
    ts->wallResistance_cgs = 1000.0;
    ts->wallStiffness_cgs  = STANDARD_WALL_STIFFNESS_CGS;

    ts->articulator = OTHER_ARTICULATOR;
  }

  ts = &tracheaSection[0];
  ts->area_cm2   = 4.0;     // 4 cm^2 trachea area
  ts->length_cm  = length_cm / (double)NUM_TRACHEA_SECTIONS;
  ts->volume_cm3 = ts->area_cm2 * ts->length_cm;

  ts = &tracheaSection[1];
  ts->area_cm2   = 3.0;     // 3 cm^2 trachea area
  ts->length_cm  = length_cm / (double)NUM_TRACHEA_SECTIONS;
  ts->volume_cm3 = ts->area_cm2 * ts->length_cm;
}


// ****************************************************************************
/// Init the nasal cavity.
// ****************************************************************************

void Tube::initNasalCavity(const double length_cm)
{
  const double NOSE_WALL_MASS_CGS = STANDARD_WALL_MASS_CGS;
  const double NOSE_WALL_RESISTANCE_CGS = STANDARD_WALL_RESISTANCE_CGS;
  const double NOSE_WALL_STIFFNESS_CGS = STANDARD_WALL_STIFFNESS_CGS;

  int i;
  Section *ts = NULL;

  const double EPSILON = 0.000000001;
  if (fabs(length_cm - nasalCavityLength_cm) < EPSILON)
  {
    return;
  }

  nasalCavityLength_cm = length_cm;

  // ****************************************************************
  // Nose sections and nasal sinuses.
  // Measurements by Dang & Honda (1994), Subject 1
  // The values were taken from the center of the individual tube 
  // sections.
  // ****************************************************************

  const double NOSE_LENGTH_CM[NUM_NOSE_SECTIONS] =
    { 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6,
      0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6 };

  const double NOSE_AREA_CM2[NUM_NOSE_SECTIONS] =
    { 1.63, 2.07, 2.72, 3.59, 4.24, 3.26, 3.04, 3.04, 2.72, 2.5,
      2.39, 2.39, 1.85, 0.76, 1.41, 1.74, 1.30, 1.74, 0.76 };


  const double SINUS_VOLUME_CM3[NUM_SINUS_SECTIONS] = { 11.3, 6.8, 33.0, 6.2 };
  const double NECK_LENGTH_CM[NUM_SINUS_SECTIONS]   = { 0.3,   0.3,   0.45,  1.0 };
  const double NECK_AREA_CM2[NUM_SINUS_SECTIONS]    = { 0.185, 0.185, 0.145, 0.11 };

  double lengthFactor = length_cm / DEFAULT_NASAL_CAVITY_LENGTH_CM;

  // The nasal cavity

  for (i=0; i < NUM_NOSE_SECTIONS; i++)
  {
    ts = &noseSection[i];

    ts->pos_cm = 0.0;
    
    // The area is only initialized once, because it can also be
    // changed by setVelumOpening(...).
    
    if (staticPartsInitialized == false)
    {
      ts->area_cm2 = NOSE_AREA_CM2[i];
    }

    ts->length_cm = lengthFactor * NOSE_LENGTH_CM[i];
    ts->volume_cm3 = ts->area_cm2 * ts->length_cm;

    ts->wallMass_cgs = NOSE_WALL_MASS_CGS;
    ts->wallResistance_cgs = NOSE_WALL_RESISTANCE_CGS;
    ts->wallStiffness_cgs = NOSE_WALL_STIFFNESS_CGS;

    ts->articulator = OTHER_ARTICULATOR;
  }

  // The paranasal sinuses (they never change).

  if (staticPartsInitialized == false)
  {
    int k;

    for (i=0; i < NUM_SINUS_SECTIONS; i++)
    {
      ts = &sinusSection[i];

      k = SINUS_COUPLING_SECTION[i];
      ts->pos_cm     = 0.0;
      ts->area_cm2   = NECK_AREA_CM2[i];
      ts->length_cm  = NECK_LENGTH_CM[i];
      ts->volume_cm3 = SINUS_VOLUME_CM3[i];

  	  ts->wallMass_cgs       = 0.0;
	    ts->wallResistance_cgs = 6500.0;
	    ts->wallStiffness_cgs  = STANDARD_WALL_STIFFNESS_CGS;

      ts->articulator = OTHER_ARTICULATOR;
    }
  }
}


// ****************************************************************************
/// Init the sinus piriformis.
/// length_cm: The effective acousic length (including end correction).
/// volume_cm3: The total volume of both the left and the right sinus.
// ****************************************************************************

void Tube::initPiriformFossa(const double length_cm, const double volume_cm3)
{
  int i;
  Section *ts = NULL;

  const double EPSILON = 0.000000001;
  if ((fabs(length_cm - piriformFossaLength_cm) < EPSILON) && 
      (fabs(volume_cm3 - piriformFossaVolume_cm3) < EPSILON))
  {
    return;
  }

  // Keep in mind the set values.
  piriformFossaLength_cm = length_cm;
  piriformFossaVolume_cm3 = volume_cm3;

  double maxArea_cm2 = 2.0 * volume_cm3 / length_cm;
  double sectionLength_cm = length_cm / (double)NUM_FOSSA_SECTIONS;

  for (i = 0; i < NUM_FOSSA_SECTIONS; i++)
  {
    ts = &fossaSection[i];

    ts->pos_cm = 0.0;
    ts->area_cm2 = maxArea_cm2 * (1.0 - (i + 0.5) / (double)NUM_FOSSA_SECTIONS);
    ts->length_cm = sectionLength_cm;
    ts->volume_cm3 = ts->area_cm2 * ts->length_cm;

    ts->wallMass_cgs = STANDARD_WALL_MASS_CGS;
    ts->wallResistance_cgs = STANDARD_WALL_RESISTANCE_CGS;
    ts->wallStiffness_cgs = STANDARD_WALL_STIFFNESS_CGS;

    ts->articulator = OTHER_ARTICULATOR;
  }
}


// ****************************************************************************
/// Resets the dynamic part of the area function.
// ****************************************************************************

void Tube::resetDynamicPart()
{
  int i;
  Section *ts = NULL;

  // ****************************************************************
  // Glottal sections.
  // ****************************************************************

  for (i=0; i < NUM_GLOTTIS_SECTIONS; i++)
  {
    ts = &glottisSection[i];

    ts->pos_cm     = 0.0;
    ts->area_cm2   = 0.1;     // 0.1 cm^2 glottal area
    ts->length_cm  = 0.3;     // 0.3 cm tube section length
    ts->volume_cm3 = ts->area_cm2 * ts->length_cm;

    ts->wallMass_cgs       = STANDARD_WALL_MASS_CGS;
    ts->wallResistance_cgs = STANDARD_WALL_RESISTANCE_CGS;
    ts->wallStiffness_cgs = STANDARD_WALL_STIFFNESS_CGS;

    ts->articulator = VOCAL_FOLDS;
  }

  // ****************************************************************
  // Pharynx and mouth sections.
  // ****************************************************************

  for (i=0; i < NUM_PHARYNX_MOUTH_SECTIONS; i++)
  {
    ts = &pharynxMouthSection[i];

    ts->pos_cm     = 0.0;
    ts->area_cm2   = 4.0;    // 4 cm^2
    ts->length_cm  = 16.0 / (double)NUM_PHARYNX_MOUTH_SECTIONS;
    ts->volume_cm3 = ts->area_cm2 * ts->length_cm;

    ts->wallMass_cgs       = STANDARD_WALL_MASS_CGS;
    ts->wallResistance_cgs = STANDARD_WALL_RESISTANCE_CGS;
    ts->wallStiffness_cgs = STANDARD_WALL_STIFFNESS_CGS;

    ts->articulator = OTHER_ARTICULATOR;
  }

}


// ****************************************************************************
/// Sets the geometry of the pharynx and mouth sections.
/// \param length_cm Lengths of the NUM_PHARYNX_MOUTH_SECTIONS in cm.
/// \param area_cm2 Areas of the NUM_SUPRAGLOTTAL_SECTIONS in cm^2.
// ****************************************************************************

void Tube::setPharynxMouthGeometry(const double *length_cm, const double *area_cm2,
  const Articulator *articulator, const double incisorPos_cm, const double tongueTipSideElevation)
{
  int i;
  Section *ts = NULL;
  double pos_cm = 0.0;

  for (i=0; i < NUM_PHARYNX_MOUTH_SECTIONS; i++)
  {
    ts = &pharynxMouthSection[i];

    ts->pos_cm     = pos_cm;
    ts->length_cm  = length_cm[i];
    ts->area_cm2   = area_cm2[i];
    if (ts->area_cm2 < MIN_AREA_CM2) { ts->area_cm2 = MIN_AREA_CM2; }
    ts->volume_cm3 = ts->area_cm2 * ts->length_cm;

    ts->articulator = articulator[i];
    
    pos_cm+= ts->length_cm;
  }

  this->teethPosition_cm = incisorPos_cm;
  this->tongueTipSideElevation = tongueTipSideElevation;

  calcPositions();
}


// ****************************************************************************
/// Sets the geometry of the lower and upper glottis sections.
/// \param length_cm Lengths of the NUM_GLOTTIS_SECTIONS in cm.
/// \param area_cm2 Areas of the NUM_GLOTTIS_SECTIONS in cm^2.
// ****************************************************************************

void Tube::setGlottisGeometry(const double *length_cm, const double *area_cm2)
{
  int i;
  Section *ts = NULL;
  double pos_cm = 0.0;

  for (i=0; i < NUM_GLOTTIS_SECTIONS; i++)
  {
    ts = &glottisSection[i];

    ts->pos_cm     = pos_cm;
    ts->length_cm  = length_cm[i];
    ts->area_cm2   = area_cm2[i];
    if (ts->area_cm2 < MIN_AREA_CM2) { ts->area_cm2 = MIN_AREA_CM2; }
    ts->volume_cm3 = ts->area_cm2 * ts->length_cm;

    ts->articulator = VOCAL_FOLDS;

    pos_cm+= ts->length_cm;
  }

  calcPositions();
}


// ****************************************************************************
/// Convenience function to set both sections of the glottis to the same given 
/// area.
// ****************************************************************************

void Tube::setGlottisArea(const double area_cm2)
{
  double l_cm[NUM_GLOTTIS_SECTIONS] = { 0.3, 0.3 };
  double a_cm2[NUM_GLOTTIS_SECTIONS]  = { area_cm2, area_cm2 };

  setGlottisGeometry(l_cm, a_cm2);
}

// ****************************************************************************
/// Modifies the first few areas of the nose for a naso-pharyngeal port of the
/// given size.
// ****************************************************************************

void Tube::setVelumOpening(const double openingArea_cm2)
{
  const int N = 4;  // Only the first N-1 sections of the nose change their shape
  int i;
  Section *ts = NULL;
  double targetArea_cm2 = noseSection[N].area_cm2;

  for (i=0; i < N; i++)
  {
    ts = &noseSection[i];
    ts->area_cm2 = openingArea_cm2 + ((double)(i*i)*(targetArea_cm2 - openingArea_cm2)) / (double)(N*N);
    if (ts->area_cm2 < MIN_AREA_CM2) { ts->area_cm2 = MIN_AREA_CM2; }
    ts->volume_cm3 = ts->area_cm2 * ts->length_cm;
  }
}


// ****************************************************************************
/// Sets the strength (amplification) of the aspiration noise source.
/// Null dB corresponds to the maximum strength (aspiration for voiceless
/// plosives). Negative numbers correspond to less strength.
// ****************************************************************************

void Tube::setAspirationStrength(const double aspirationStrength_dB)
{
  this->aspirationStrength_dB = aspirationStrength_dB;
}


// ****************************************************************************
/// Interpolate linearly between leftTube and rightTube to configure this 
/// current tube.
/// The parameter ratio = [0, 1] determines the contribution of rightTube,
/// i.e., thisTube = (1.0-ratio)*leftTube + ratio*rightTube.
// ****************************************************************************

void Tube::interpolate(const Tube *leftTube, const Tube *rightTube, const double ratio)
{
  int i;
  double length_cm[NUM_PHARYNX_MOUTH_SECTIONS];
  double area_cm2[NUM_PHARYNX_MOUTH_SECTIONS];
  Articulator articulator[NUM_PHARYNX_MOUTH_SECTIONS];

  double ratio1 = 1.0 - ratio;

  for (i=0; i < NUM_PHARYNX_MOUTH_SECTIONS; i++)
  {
    area_cm2[i]  = ratio1*leftTube->pharynxMouthSection[i].area_cm2 + ratio*rightTube->pharynxMouthSection[i].area_cm2;
    length_cm[i] = ratio1*leftTube->pharynxMouthSection[i].length_cm + ratio*rightTube->pharynxMouthSection[i].length_cm;

    if (ratio < 0.5)
    {
      articulator[i] = leftTube->pharynxMouthSection[i].articulator;
    }
    else
    {
      articulator[i] = rightTube->pharynxMouthSection[i].articulator;
    }
  }

  // Interpolate the teeth position and the tongue tip side elevation.

  double incisorPos_cm = ratio1*leftTube->teethPosition_cm + ratio*rightTube->teethPosition_cm;
  double tongueTipSideElevation = ratio1 * leftTube->tongueTipSideElevation + ratio * rightTube->tongueTipSideElevation;

  // Set the dimensions of the main vocal tract.

  this->setPharynxMouthGeometry(length_cm, area_cm2, articulator, incisorPos_cm, tongueTipSideElevation);

  // Interpolate the velum area.

  double leftOpening_cm2 = leftTube->getVelumOpening_cm2();
  double rightOpening_cm2 = rightTube->getVelumOpening_cm2();
  
  this->setVelumOpening( ratio1*leftOpening_cm2 + ratio*rightOpening_cm2 );

  // Interpolate the aspiration strength.

  double leftAspirationStrength_dB = leftTube->aspirationStrength_dB;
  double rightAspirationStrength_dB = rightTube->aspirationStrength_dB;

  this->aspirationStrength_dB = ratio1*leftAspirationStrength_dB + ratio*rightAspirationStrength_dB;

  // Interpolate the length of the subglottal system, the nasal cavity,
  // and the piriform fossae.

  double leftSubglottalCavityLength_cm, rightSubglottalCavityLength_cm;
  double leftNasalCavityLength_cm, rightNasalCavityLength_cm;
  double leftPiriformFossaLength_cm, rightPiriformFossaLength_cm;
  double leftPiriformFossaVolume_cm3, rightPiriformFossaVolume_cm3;

  leftTube->getStaticTubeDimensions(
    leftSubglottalCavityLength_cm,
    leftNasalCavityLength_cm, 
    leftPiriformFossaLength_cm,
    leftPiriformFossaVolume_cm3);

  rightTube->getStaticTubeDimensions(
    rightSubglottalCavityLength_cm,
    rightNasalCavityLength_cm,
    rightPiriformFossaLength_cm,
    rightPiriformFossaVolume_cm3);

  this->initSubglottalCavity(ratio1*leftSubglottalCavityLength_cm + ratio*rightSubglottalCavityLength_cm);
  this->initNasalCavity(ratio1*leftNasalCavityLength_cm + ratio*rightNasalCavityLength_cm);
  this->initPiriformFossa(ratio1*leftPiriformFossaLength_cm + ratio*rightPiriformFossaLength_cm,
    ratio1*leftPiriformFossaVolume_cm3 + ratio*rightPiriformFossaVolume_cm3);
}


// ****************************************************************************
/// Returns the current naso-pharyngeal port area.
// ****************************************************************************

double Tube::getVelumOpening_cm2() const
{
  return noseSection[0].area_cm2;
}


// ****************************************************************************
/// Returns a couple of static tube dimensions.
// ****************************************************************************

void Tube::getStaticTubeDimensions(double &subglottalCavityLength_cm,
  double &nasalCavityLength_cm, double &piriformFossaLength_cm,
  double &piriformFossaVolume_cm3) const
{
  subglottalCavityLength_cm = this->subglottalCavityLength_cm;
  nasalCavityLength_cm = this->nasalCavityLength_cm;
  piriformFossaLength_cm = this->piriformFossaLength_cm;
  piriformFossaVolume_cm3 = this->piriformFossaVolume_cm3;
}


// ****************************************************************************
/// Creates the list with pointers to ALL tube sections of the model for
/// convenient iterations of all sections.
// ****************************************************************************

void Tube::createSectionList()
{
  int i;

  for (i=0; i < NUM_TRACHEA_SECTIONS; i++)
  {
    section[FIRST_TRACHEA_SECTION + i] = &tracheaSection[i];
  }

  for (i=0; i < NUM_GLOTTIS_SECTIONS; i++)
  {
    section[FIRST_GLOTTIS_SECTION + i] = &glottisSection[i];
  }

  for (i=0; i < NUM_PHARYNX_MOUTH_SECTIONS; i++)
  {
    section[FIRST_PHARYNX_SECTION + i] = &pharynxMouthSection[i];
  }

  for (i=0; i < NUM_NOSE_SECTIONS; i++)
  {
    section[FIRST_NOSE_SECTION + i] = &noseSection[i];
  }

  for (i=0; i < NUM_FOSSA_SECTIONS; i++)
  {
    section[FIRST_FOSSA_SECTION + i] = &fossaSection[i];
  }

  for (i=0; i < NUM_SINUS_SECTIONS; i++)
  {
    section[FIRST_SINUS_SECTION + i] = &sinusSection[i];
  }

}


// ****************************************************************************
/// Calculate the absolute positions of the tube sections.
// ****************************************************************************

void Tube::calcPositions()
{
  int i;
  double pos;
  const double REFERENCE_POS = 0.0;

  pos = REFERENCE_POS;

  // ****************************************************************
  // The two glottal tube sections.
  // ****************************************************************

  pos-= glottisSection[1].length_cm;
  glottisSection[1].pos_cm = pos;

  pos-= glottisSection[0].length_cm;
  glottisSection[0].pos_cm = pos;

  // ****************************************************************
  // Trachea sections.
  // ****************************************************************

  for (i = NUM_TRACHEA_SECTIONS-1; i >= 0; i--)
  {
    pos-= tracheaSection[i].length_cm;
    tracheaSection[i].pos_cm = pos;
  }

  // ****************************************************************
  // Pharynx and mouth segments.
  // ****************************************************************

  double fossaPos = 0.0;
  double nosePos = 0.0;
  pos = REFERENCE_POS;
  
  for (i=0; i < NUM_PHARYNX_MOUTH_SECTIONS; i++)
  {
    pharynxMouthSection[i].pos_cm = pos;
    pos+= pharynxMouthSection[i].length_cm;

    if (i == FOSSA_COUPLING_SECTION) { fossaPos = pos; }
    if (i == NUM_PHARYNX_SECTIONS - 1) { nosePos = pos; }
  }

  // ****************************************************************
  // Nose and piriform fossa segments.
  // ****************************************************************

  pos = nosePos;
  for (i=0; i < NUM_NOSE_SECTIONS; i++)
  {
    noseSection[i].pos_cm = pos;
    pos+= noseSection[i].length_cm;
  }

  pos = fossaPos;
  for (i=0; i < NUM_FOSSA_SECTIONS; i++)
  {
    fossaSection[i].pos_cm = pos;
    pos+= fossaSection[i].length_cm;
  }

  // ****************************************************************
  // Position of the paranasal sinuses.
  // ****************************************************************

  int couplingSection;

  for (i=0; i < NUM_SINUS_SECTIONS; i++)
  {
    couplingSection = SINUS_COUPLING_SECTION[i];
    sinusSection[i].pos_cm = noseSection[couplingSection].pos_cm + noseSection[couplingSection].length_cm;
  }

}

// ****************************************************************************
// ****************************************************************************

void Tube::print()
{
  int i;

  for (i=0; i < NUM_SECTIONS; i++)
  {
    if (i == FIRST_TRACHEA_SECTION)
    {
      printf("\n# Trachea sections\n");
    }
    else
    if (i == FIRST_GLOTTIS_SECTION)
    {
      printf("\n# Glottis sections\n");
    }
    else
    if (i == FIRST_PHARYNX_SECTION)
    {
      printf("\n# Pharynx-mouth sections\n");
    }
    else
    if (i == FIRST_NOSE_SECTION)
    {
      printf("\n# Nose sections\n");
    }
    else
    if (i == FIRST_FOSSA_SECTION)
    {
      printf("\n# Piriform fossa sections\n");
    }
    else
    if (i == FIRST_SINUS_SECTION)
    {
      printf("\n# Paranasal sinus sections\n");
    }

    printf("#%2d: x=%6.2f cm l=%6.2f cm A=%6.2f cm^2 V=%6.2f cm^3\n", 
      i, 
      section[i]->pos_cm,
      section[i]->length_cm,
      section[i]->area_cm2,
      section[i]->volume_cm3
    );
  }
}

// ****************************************************************************
/// Overwrite the assignment operator!
// ****************************************************************************

void Tube::operator=(const Tube &t)
{
  int i;

  for (i=0; i < NUM_SECTIONS; i++)
  {
    *this->section[i] = *t.section[i];
  }

  this->teethPosition_cm = t.teethPosition_cm;
  this->aspirationStrength_dB = t.aspirationStrength_dB;
  this->tongueTipSideElevation = t.tongueTipSideElevation;
}

// ****************************************************************************
/// Are the tubes identical ?
// ****************************************************************************

bool Tube::operator==(const Tube &t)
{
  if ((this->teethPosition_cm != t.teethPosition_cm) || 
    (this->tongueTipSideElevation != t.tongueTipSideElevation))
  {
    return false;
  }

  bool equal = true;
  int i;
  Section *a, *b;

  for (i=0; i < NUM_SECTIONS; i++)
  {
    a = this->section[i];
    b = t.section[i];

    if ((a->area_cm2           != b->area_cm2) ||
        (a->length_cm          != b->length_cm) ||
        (a->pos_cm             != b->pos_cm) ||
        (a->volume_cm3         != b->volume_cm3) ||
        (a->wallMass_cgs       != b->wallMass_cgs) ||
        (a->wallResistance_cgs != b->wallResistance_cgs) ||
        (a->wallStiffness_cgs  != b->wallStiffness_cgs) ||
        (a->articulator        != b->articulator))
    {
      equal = false;
      break;
    }
  }

  return equal;
}

// ****************************************************************************
/// Are the tubes not identical ?
// ****************************************************************************

bool Tube::operator!=(const Tube &t)
{
  return !(*this == t);
}

// ****************************************************************************
