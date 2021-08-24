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

#ifndef __VOCALTRACT_H__
#define __VOCALTRACT_H__

#include <string>
#include "Surface.h"
#include "Splines.h"
#include "Tube.h"
#include "XmlNode.h"


// ****************************************************************************
/// This class represents the geometrical 3D vocal tract model.
// Note: All length and area variables are in cm or cm^2 !!
// ****************************************************************************

class VocalTract
{
  // **************************************************************************
  /// Public data (constants, data types, variables).
  // **************************************************************************

public:
  static const int NUM_EPIGLOTTIS_RIBS = 4;
  static const int NUM_EPIGLOTTIS_POINTS = 5;
  static const int NUM_UVULA_RIBS = 4;
  static const int NUM_UVULA_POINTS = 5;

  static const int NUM_LARYNX_RIBS  = 5;
  static const int NUM_PHARYNX_RIBS = 3;
  static const int NUM_VELUM_RIBS   = 6;
  static const int NUM_PALATE_RIBS  = 9;    ///< for both upper and lower jaw
  static const int NUM_JAW_RIBS     = NUM_PALATE_RIBS;
  static const int NUM_THROAT_RIBS  = NUM_PHARYNX_RIBS;

  static const int NUM_UPPER_COVER_RIBS = NUM_LARYNX_RIBS + NUM_PHARYNX_RIBS + NUM_VELUM_RIBS + NUM_JAW_RIBS;
  static const int NUM_UPPER_COVER_POINTS = 6;

  static const int NUM_LOWER_COVER_RIBS = NUM_LARYNX_RIBS + NUM_THROAT_RIBS + NUM_JAW_RIBS;
  static const int NUM_LOWER_COVER_POINTS = 5;

  static const int NUM_DYNAMIC_TONGUE_RIBS = 33;
  static const int NUM_STATIC_TONGUE_RIBS  = 4;
  static const int NUM_TONGUE_RIBS = NUM_DYNAMIC_TONGUE_RIBS + NUM_STATIC_TONGUE_RIBS;
  static const int NUM_TONGUE_POINTS = 11;
  static const int MAX_TONGUE_RIBS_GLOBAL = 128;

  static const int NUM_TEETH_RIBS = 3*(NUM_JAW_RIBS-1) + 1;
  static const int NUM_TEETH_POINTS = 5;
    
  static const int NUM_LIP_RIBS = NUM_JAW_RIBS;
  static const int NUM_INNER_LIP_POINTS = 5;
  static const int NUM_OUTER_LIP_POINTS = 5;
  static const int NUM_LIP_POINTS = NUM_INNER_LIP_POINTS + NUM_OUTER_LIP_POINTS;

  static const int NUM_FILL_RIBS = NUM_JAW_RIBS + 3;
  static const int NUM_FILL_POINTS = 4;

  static const int NUM_RADIATION_RIBS = NUM_LIP_RIBS + 4;
  static const int NUM_RADIATION_POINTS = 6;

  // ****************************************************************
  // Constants for the cross-sections.
  // ****************************************************************
  
  static const int NUM_PROFILE_SAMPLES = 96;  ///< Number of points in lateral direction
  static const double PROFILE_LENGTH;
  static const double PROFILE_SAMPLE_LENGTH;
  static const double INVALID_PROFILE_SAMPLE;
  static const double MIN_PROFILE_VALUE;
  static const double MAX_PROFILE_VALUE;
  static const double EXTREME_PROFILE_VALUE;

  // ****************************************************************
  // Other constants.
  // ****************************************************************

  static const int NUM_CENTERLINE_POINTS_EXPONENT = 7;
  static const int NUM_CENTERLINE_POINTS = (1 << NUM_CENTERLINE_POINTS_EXPONENT) + 1;
  static const int NUM_TUBE_SECTIONS = Tube::NUM_PHARYNX_MOUTH_SECTIONS;
  /// Coupling section for the nasal cavity.
  static const int NUM_PHARYNX_SECTIONS = Tube::NUM_PHARYNX_SECTIONS;
  /// Coupling section for the piriformis fossa.
  static const int PIRIFORM_FOSSA_SECTION = Tube::NUM_FOSSA_SECTIONS;
  
  // ****************************************************************
  // Identifier for the different surfaces.
  // ****************************************************************

  enum SurfaceIndex 
  { 
    UPPER_TEETH, LOWER_TEETH, 
    UPPER_COVER, LOWER_COVER,
    UPPER_LIP, LOWER_LIP,
    PALATE, MANDIBLE, LOWER_TEETH_ORIGINAL,
    LOW_VELUM, MID_VELUM, HIGH_VELUM,
    NARROW_LARYNX_FRONT, NARROW_LARYNX_BACK,
    WIDE_LARYNX_FRONT, WIDE_LARYNX_BACK,
    TONGUE,
    UPPER_COVER_TWOSIDE,
    LOWER_COVER_TWOSIDE,
    UPPER_TEETH_TWOSIDE, 
    LOWER_TEETH_TWOSIDE, 
    UPPER_LIP_TWOSIDE, 
    LOWER_LIP_TWOSIDE,
    LEFT_COVER,        ///< Fill surface between the upper and lower teeth 
    RIGHT_COVER,       ///< Fill surface between the upper and lower teeth 
    UVULA_ORIGINAL,
    UVULA,
    UVULA_TWOSIDE,
    EPIGLOTTIS_ORIGINAL,
    EPIGLOTTIS,
    EPIGLOTTIS_TWOSIDE,
    RADIATION,
    NUM_SURFACES
  };

  // ****************************************************************
  // Time-variant vocal tract parameters.
  // ****************************************************************

  struct Param
  {
    double x;         ///< Parameter value set by the user or defined as target
    double limitedX;  ///< Parameter value limited by biomechanical constraints
    double min;
    double max;
    double neutral;
    string abbr;      ///< Abbreviation, e.g. "TCX"
    string name;      ///< Long name, e.g. "Horizontal tongue body position"
  };

  // ****************************************************************

  enum ParamIndex
  {
    HX, HY, JX, JA,
    LP, LD, VS, VO,
    TCX, TCY, TTX, TTY,
    TBX, TBY, TRX, TRY,
    TS1, TS2, TS3,
    NUM_PARAMS
  };

  // ****************************************************************
  // Anatomical, articulation-invariant vocal tract shape parameters.
  // ****************************************************************

  struct Anatomy
  {
    Point3D palatePoints[NUM_PALATE_RIBS];
    double palateAngle_deg[NUM_PALATE_RIBS];
    double palateHeight_cm[NUM_PALATE_RIBS];
    double upperTeethHeight_cm[NUM_PALATE_RIBS];
    double upperTeethWidthTop_cm[NUM_PALATE_RIBS];
    double upperTeethWidthBottom_cm[NUM_PALATE_RIBS];
    
    Point2D jawFulcrum;
    Point2D jawRestPos;
    double toothRootLength_cm;
    Point3D jawPoints[NUM_JAW_RIBS];
    double jawAngle_deg[NUM_JAW_RIBS];
    double jawHeight_cm[NUM_JAW_RIBS];
    double lowerTeethHeight_cm[NUM_JAW_RIBS];
    double lowerTeethWidthTop_cm[NUM_JAW_RIBS];
    double lowerTeethWidthBottom_cm[NUM_JAW_RIBS];

    double tongueTipRadius_cm;
    double tongueCenterRadiusX_cm;
    double tongueCenterRadiusY_cm;
    bool automaticTongueRootCalc;
    double tongueRootTrxSlope;
    double tongueRootTrxIntercept;
    double tongueRootTrySlope;
    double tongueRootTryIntercept;

    double lipsWidth_cm;

    double uvulaWidth_cm;
    double uvulaHeight_cm;
    double uvulaDepth_cm;
    Point2D velumLowPoints[NUM_VELUM_RIBS-1];
    Point2D velumMidPoints[NUM_VELUM_RIBS-1];
    Point2D velumHighPoints[NUM_VELUM_RIBS-1];
    double maxNasalPortArea_cm2;

    Point2D pharynxFulcrum;
    double pharynxRotationAngle_deg;
    double pharynxTopRibY_cm;
    double pharynxUpperDepth_cm;
    double pharynxLowerDepth_cm;
    double pharynxBackWidth_cm;

    double epiglottisWidth_cm;
    double epiglottisHeight_cm;
    double epiglottisDepth_cm;
    double epiglottisAngle_deg;

    double larynxUpperDepth_cm;
    double larynxLowerDepth_cm;
    Point2D larynxWidePoints[8];
    Point2D larynxNarrowPoints[8];

    double piriformFossaLength_cm;
    double piriformFossaVolume_cm3;
    double subglottalCavityLength_cm;
    double nasalCavityLength_cm;

    // Factors controlling the intrinsic, direction-dependent
    // velocities of the articulators (i.e. parameters).
    double positiveVelocityFactor[NUM_PARAMS];
    double negativeVelocityFactor[NUM_PARAMS];
  } anatomy;

  // ****************************************************************
  // A vocal tract shape or configuration, for example for the 
  // vowels [i], [a], [u].
  // ****************************************************************

  struct Shape
  {
    string name;
    double param[NUM_PARAMS];    ///< Parameter value
  };

  // ****************************************************************
  // A line point on the center line between glottis and lips.
  // ****************************************************************

  struct CenterLinePoint
  {
    Point2D point;
    Point2D normal;
    double pos;
    double min, max;
    double reserved;
  };

  // ****************************************************************
  // Specification of a tongue rib.
  // ****************************************************************

  struct TongueRib
  {
    // The following variables must be set.
    Point2D point;
    double leftSideHeight;
    double rightSideHeight;

    // The following variables are for internal use only.
    Point2D left, right;
    Point2D normal;
    double minX, maxX;
    double minY, maxY;
  };

  // ****************************************************************

  struct CrossSection
  {
    double area;
    double circ;
    double pos;
    Tube::Articulator articulator;
  };

  // ****************************************************************

  struct TubeSection
  {
    double area;
    double circ;
    double pos;
    double length;
    Tube::Articulator articulator;
  };

  // ****************************************************************
  // EMA points.
  // ****************************************************************

  enum EmaSurface
  {
    EMA_SURFACE_TONGUE,
    EMA_SURFACE_UPPER_COVER,
    EMA_SURFACE_LOWER_COVER,
    EMA_SURFACE_UPPER_LIP,
    EMA_SURFACE_LOWER_LIP,
    NUM_EMA_SURFACES    
  };

  struct EmaPoint
  {
    string name;
    // EMA points are always in the midsagittal plane on the surface.
    EmaSurface emaSurface;
    int vertexIndex;
  };

  // ****************************************************************
  // Variables.
  // ****************************************************************

  Surface surface[NUM_SURFACES];
  Param param[NUM_PARAMS];
  vector<Shape> shapes;
  vector<EmaPoint> emaPoints;

  TongueRib tongueRib[MAX_TONGUE_RIBS_GLOBAL];

  // Guiding lines for the lip corners.
  LineStrip3D wideLipCornerPath;
  LineStrip3D narrowLipCornerPath;
  LineStrip3D lipCornerPath;

  // Center line, cross sections, and tube sections.
  
  double centerLineLength;
  CenterLinePoint roughCenterLine[NUM_CENTERLINE_POINTS];
  CenterLinePoint centerLine[NUM_CENTERLINE_POINTS];
  CrossSection crossSection[NUM_CENTERLINE_POINTS];
  TubeSection tubeSection[NUM_TUBE_SECTIONS];

  // Position and opening of the velo-pharyngal port
  double nasalPortPos_cm;
  double nasalPortArea_cm2;

  // Position of the incisors along the center line.
  double incisorPos_cm;
  
  // Picewise linear approximations of the outer teeth edges
  Point3D upperGumsInnerEdge[NUM_JAW_RIBS];
  Point3D upperGumsOuterEdge[NUM_JAW_RIBS];
  Point3D lowerGumsInnerEdge[NUM_JAW_RIBS];
  Point3D lowerGumsOuterEdge[NUM_JAW_RIBS];
  Point3D lowerGumsInnerEdgeOrig[NUM_JAW_RIBS];
  Point3D lowerGumsOuterEdgeOrig[NUM_JAW_RIBS];


  // **************************************************************************
  /// Public functions.
  // **************************************************************************

public:
  VocalTract();
  ~VocalTract();

  // ****************************************************************
  // Initialization.
  // ****************************************************************

  void init();    ///< Is automatically called by the constructor.
  void initSurfaceGrids();
  void initReferenceSurfaces();
  void initLarynx();
  void initJaws();
  void initVelum();
  void setDefaultEmaPoints();
  Point3D getEmaPointCoord(int index);
  void getEmaSurfaceVertexRange(int emaSurface, int *min, int *max);

  // ****************************************************************
  // Reading and writing of anatomy and shapes from/to xml-files.
  // ****************************************************************

  void readAnatomyXml(XmlNode *anatomyNode);
  void readShapesXml(XmlNode *shapeListNode);
  void readFromXml(const string &speakerFileName);
  void writeAnatomyXml(std::ostream &os, int indent);
  void writeShapesXml(std::ostream &os, int indent);
  void writeToXml(std::ostream &os, int indent);

  bool saveAsObjFile(const string &fileName, bool saveBothSides = true);

  // ****************************************************************
  // Calculate the surfaces, center line, and area functions
  // ****************************************************************

  void setParams(double *controlParams);
  void calculateAll();
  
  // ****************************************************************
  // Calculate all geometric surfaces.
  // ****************************************************************

  void calcSurfaces();
  void calcLips();
  void getImportantLipPoints(Point3D &onset, Point3D &corner, Point3D &F0, Point3D &F1, double &yClose);
  void calcRadiation(Point3D lipCorner);
  void getHyoidTongueTangent(Point2D &H, Point2D &T);
  void calcTongue();
  double tongueSideParamToElevation_cm(double paramValue);
  double tongueSideParamToMinArea_cm2(double paramValue);
  void restrictTongueParams();
  Point2D limitEllipsePos(Point2D C, double rx, double ry, LineStrip2D &border, Point2D A);
  void calcTongueRibs();
  void verifyTongueRibNormal(int rigid, int flexible);
  
  // ****************************************************************
  // Calculate the acoustic center line of the vocal tract.
  // ****************************************************************

  void calcCenterLine();
  void verifyCenterLineNormal(int left, int center, int right);

  // ****************************************************************
  // Calculate the piecewise linear area function.
  // ****************************************************************

  void calcCrossSections();
  void getCrossProfiles(Point2D P, Point2D v, double *upperProfile, double *lowerProfile, 
    bool considerTongue, Tube::Articulator &articulator, bool debug = false);
  void insertUpperProfileLine(Point2D P0, Point2D P1, int surfaceIndex, 
    double *upperProfile, int *upperProfileSurface);
  void insertLowerProfileLine(Point2D P0, Point2D P1, int surfaceIndex, 
    double *lowerProfile, int *lowerProfileSurface);
  void insertLowerCoverProfileLine(Point2D P0, Point2D P1, int surfaceIndex,
    double *upperProfile, int *upperProfileSurface, double *lowerProfile, int *lowerProfileSurface);
  void getCrossSection(double *upperProfile, double *lowerProfile, CrossSection *section);
  
  bool exportCrossSections(const string &fileName);
  bool exportTractContourSvg(const string &fileName, bool addCenterLine, bool addCutVectors);
  void addRibPointsSvg(ostream &os, Surface *s, int rib, int firstRibPoint, int lastRibPoint);
  void addRibsSvg(ostream &os, Surface *s, int firstRib, int lastRib, int ribPoint);

  // ****************************************************************
  // Calculate the piecewise constant area function.  
  // ****************************************************************

  void crossSectionsToTubeSections();

  // ****************************************************************
  // Functions mainly for calls from outside this class.
  // ****************************************************************

  int getShapeIndex(const string &name);
  bool isVowelShapeName(const string &name);
  double getPharynxBackX(double y);
  void getTube(Tube *tube);
  double getCenterLinePos(Point2D Q, int &bestIndex, double &bestT);
  void getCutVector(double pos, Point2D &P, Point2D &v);
  void restrictParam(int index);
  bool hasUnsavedChanges();
  void clearUnsavedChanges();

  void storeControlParams();
  void restoreControlParams();

  // **************************************************************************
  /// Private data.
  // **************************************************************************

private:
  bool intersectionsPrepared[NUM_SURFACES];  // For the fast intersection method
  bool hasStoredControlParams;
  double storedControlParams[NUM_PARAMS];

  LineStrip2D upperOutline;
  LineStrip2D lowerOutline;
  LineStrip2D tongueOutline;
  LineStrip2D epiglottisOutline;
};

// ****************************************************************************

#endif
