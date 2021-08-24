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

#ifndef __TUBE_H__
#define __TUBE_H__

// ****************************************************************************
/// This class represents a discrete vocal tract area function or tube 
/// including all side branches and the glottis.
/// Only the pharynx and mouth sections and the glottal sections are meant to be 
/// changed during speaking. All other sections should be left in their 
/// original state. The first few sections of the nose cavity are changed
/// indirectly by calling setNasalPortArea().
// ****************************************************************************

class Tube
{
  // **************************************************************************
  // Public data.
  // **************************************************************************

public:
  // Each tube section is in the region of one of these articulators.
  // This information is needed for the synthesis of fricatives.
  enum Articulator
  {
    VOCAL_FOLDS,
    TONGUE,
    LOWER_INCISORS,
    LOWER_LIP,
    OTHER_ARTICULATOR,
    NUM_ARTICULATORS
  };

  struct Section
  {
    // Maybe later add a "shape factor" as defined by Dang, Honda, Suzuki (JASA, 1994)
    // to model elliptical sections (shape factor = 1 for a circle).

    double pos_cm;         ///< Relative position of the tube section
    double area_cm2;
    double length_cm;
    double volume_cm3;
    double wallMass_cgs;
    double wallStiffness_cgs;
    double wallResistance_cgs;
    Articulator articulator;
  };

  static const int NUM_TRACHEA_SECTIONS  = 23;    // 1 cm per section = 23 cm in total
  static const int NUM_GLOTTIS_SECTIONS  = 2;
  static const int NUM_PHARYNX_SECTIONS  = 16;
  static const int NUM_MOUTH_SECTIONS    = 24;
  static const int NUM_PHARYNX_MOUTH_SECTIONS = NUM_PHARYNX_SECTIONS + NUM_MOUTH_SECTIONS;
  static const int NUM_NOSE_SECTIONS     = 19;
  static const int NUM_SINUS_SECTIONS    = 4;
  static const int NUM_FOSSA_SECTIONS    = 5;
  static const int NUM_SECTIONS          = NUM_TRACHEA_SECTIONS + NUM_GLOTTIS_SECTIONS + 
    NUM_PHARYNX_MOUTH_SECTIONS + NUM_NOSE_SECTIONS + NUM_SINUS_SECTIONS + NUM_FOSSA_SECTIONS;

  static const int FIRST_TRACHEA_SECTION = 0;
  static const int FIRST_GLOTTIS_SECTION = FIRST_TRACHEA_SECTION + NUM_TRACHEA_SECTIONS;
  static const int LOWER_GLOTTIS_SECTION = FIRST_GLOTTIS_SECTION;
  static const int UPPER_GLOTTIS_SECTION = FIRST_GLOTTIS_SECTION + 1;
  static const int FIRST_PHARYNX_SECTION = FIRST_GLOTTIS_SECTION + NUM_GLOTTIS_SECTIONS;
  static const int FIRST_MOUTH_SECTION   = FIRST_PHARYNX_SECTION + NUM_PHARYNX_SECTIONS;
  static const int FIRST_NOSE_SECTION    = FIRST_PHARYNX_SECTION + NUM_PHARYNX_MOUTH_SECTIONS;
  static const int FIRST_FOSSA_SECTION   = FIRST_NOSE_SECTION + NUM_NOSE_SECTIONS;
  static const int FIRST_SINUS_SECTION   = FIRST_FOSSA_SECTION + NUM_FOSSA_SECTIONS;

  static const int LAST_TRACHEA_SECTION  = FIRST_TRACHEA_SECTION + NUM_TRACHEA_SECTIONS - 1;
  static const int LAST_PHARYNX_SECTION  = FIRST_PHARYNX_SECTION + NUM_PHARYNX_SECTIONS - 1;
  static const int LAST_MOUTH_SECTION    = FIRST_MOUTH_SECTION + NUM_MOUTH_SECTIONS - 1;
  static const int LAST_NOSE_SECTION     = FIRST_NOSE_SECTION + NUM_NOSE_SECTIONS - 1;
  static const int LAST_FOSSA_SECTION    = FIRST_FOSSA_SECTION + NUM_FOSSA_SECTIONS - 1;
  static const int LAST_SINUS_SECTION    = FIRST_SINUS_SECTION + NUM_SINUS_SECTIONS - 1;

  static const int FOSSA_COUPLING_SECTION = 3;    // Piriform fossa coupling section
  static const int SINUS_COUPLING_SECTION[NUM_SINUS_SECTIONS];

  static const double STANDARD_WALL_MASS_CGS;
  static const double STANDARD_WALL_RESISTANCE_CGS;
  static const double STANDARD_WALL_STIFFNESS_CGS;

  static const double MIN_AREA_CM2;
  static const double DEFAULT_ASPIRATION_STRENGTH_DB;

  static const double DEFAULT_SUBGLOTTAL_CAVITY_LENGTH_CM;
  static const double DEFAULT_NASAL_CAVITY_LENGTH_CM;
  static const double DEFAULT_PIRIFORM_FOSSA_LENGTH_CM;
  static const double DEFAULT_PIRIFORM_FOSSA_VOLUME_CM3;

  static const char ARTICULATOR_LETTER[NUM_ARTICULATORS];

  // ****************************************************************

  Section tracheaSection[NUM_TRACHEA_SECTIONS];
  Section glottisSection[NUM_GLOTTIS_SECTIONS];
  Section pharynxMouthSection[NUM_PHARYNX_MOUTH_SECTIONS];
  Section noseSection[NUM_NOSE_SECTIONS];
  Section sinusSection[NUM_SINUS_SECTIONS];
  Section fossaSection[NUM_FOSSA_SECTIONS];

  /// List of pointers to ALL tube sections of the model
  Section *section[NUM_SECTIONS];

  /// Position of the incisors.
  double teethPosition_cm;

  /// Aspiration strength at the glottis.
  double aspirationStrength_dB;

  /// Elevation of the tongue tip side (corresponding to the TS3 
  /// parameter of the vocal tract model).
  double tongueTipSideElevation;

  // **************************************************************************
  // Public functions.
  // **************************************************************************

public:
  Tube();

  void initSubglottalCavity(const double length_cm = DEFAULT_SUBGLOTTAL_CAVITY_LENGTH_CM);
  void initNasalCavity(const double length_cm = DEFAULT_NASAL_CAVITY_LENGTH_CM);
  void initPiriformFossa(const double length_cm = DEFAULT_PIRIFORM_FOSSA_LENGTH_CM, 
    const double volume_cm3 = DEFAULT_PIRIFORM_FOSSA_VOLUME_CM3);

  void resetDynamicPart();

  void setPharynxMouthGeometry(const double *length_cm, const double *area_cm2,
    const Articulator *articulator, const double incisorPos_cm, const double tongueTipSideElevation);
  void setGlottisGeometry(const double *length_cm, const double *area_cm2);
  void setGlottisArea(const double area_cm2);
  void setVelumOpening(const double area_cm2);
  void setAspirationStrength(const double aspirationStrength_dB);

  void interpolate(const Tube *leftTube, const Tube *rightTube, const double ratio);

  double getVelumOpening_cm2() const;

  void getStaticTubeDimensions(double &subglottalCavityLength_cm,
    double &nasalCavityLength_cm, double &piriformFossaLength_cm, 
    double &piriformFossaVolume_cm3) const;

  void calcPositions();
  void print();

  void operator=(const Tube &t);
  bool operator==(const Tube &t);
  bool operator!=(const Tube &t);

  // **************************************************************************
  // Private data.
  // **************************************************************************

private:
  bool staticPartsInitialized;
  double subglottalCavityLength_cm;
  double nasalCavityLength_cm;
  double piriformFossaLength_cm;
  double piriformFossaVolume_cm3;

  // **************************************************************************
  // Private functions.
  // **************************************************************************

private:
  void createSectionList();

};

#endif

