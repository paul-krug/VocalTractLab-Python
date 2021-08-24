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

#include "VocalTract.h"
#include "Dsp.h"
#include "Geometry.h"
#include "XmlHelper.h"

#include <iomanip>
#include <iostream>
#include <cmath>
#include <fstream>


using namespace std;

// ****************************************************************************
// Initialize the constants.
// ****************************************************************************

const double VocalTract::PROFILE_LENGTH = 5.5;      // cm
const double VocalTract::PROFILE_SAMPLE_LENGTH = PROFILE_LENGTH / (double)NUM_PROFILE_SAMPLES;
const double VocalTract::INVALID_PROFILE_SAMPLE = 1000000.0;

const double VocalTract::MIN_PROFILE_VALUE = -2.75;  // cm; Should not be smaller.
const double VocalTract::MAX_PROFILE_VALUE = 10.0;  // Must be so high, so that the upper cover is always intersected !!
const double VocalTract::EXTREME_PROFILE_VALUE = 1000000.0;

static bool makeFasterIntersections = true;


// ****************************************************************************
// Constructor.
// ****************************************************************************

VocalTract::VocalTract()
{
  init();
}


// ****************************************************************************
/// Destructor.
// ****************************************************************************

VocalTract::~VocalTract() 
{ 
}



// ****************************************************************************
// Initialize the vocal tract with the current anatomy data.
// ****************************************************************************

void VocalTract::init()
{
  // ****************************************************************
  // Init all sufaces.
  // ****************************************************************

  initSurfaceGrids();
  setDefaultEmaPoints();

  // ****************************************************************
  // Initialize the tongue temporarily with a dummy shape.
  // ****************************************************************

  int i, k;
  Surface *tongue = &surface[TONGUE];

  for (i=0; i < tongue->numRibs; i++)
  {
    for (k=0; k < tongue->numRibPoints; k++)
    {
      tongue->setVertex(i, k, Point3D(-0.31, -1.02, 0));
    }
  }

  // ****************************************************************
  // Create the xml-string that defines the anatomy.
  // ****************************************************************

  std::string anatomyString = 
    "<anatomy>\n"
    "  <!--****************************************************************************-->\n"
    "  <palate>\n"
    "    <p0 x=\"0.2\" z=\"-2.3\" teeth_height=\"0.5\" top_teeth_width=\"1.05\" bottom_teeth_width=\"1.05\" palate_height=\"1.3\" palate_angle_deg=\"39.5\"/>\n"
    "    <p1 x=\"0.9\" z=\"-2.2\" teeth_height=\"0.5\" top_teeth_width=\"1.05\" bottom_teeth_width=\"1.05\" palate_height=\"1.15\" palate_angle_deg=\"39.5\"/>\n"
    "    <p2 x=\"1.8\" z=\"-2.0\" teeth_height=\"0.5\" top_teeth_width=\"1.0\" bottom_teeth_width=\"1.0\" palate_height=\"1.425\" palate_angle_deg=\"60.8\"/>\n"
    "    <p3 x=\"2.8\" z=\"-1.8\" teeth_height=\"0.5\" top_teeth_width=\"1.0\" bottom_teeth_width=\"1.0\" palate_height=\"1.6\" palate_angle_deg=\"60.8\"/>\n"
    "    <p4 x=\"3.5\" z=\"-1.6\" teeth_height=\"0.6\" top_teeth_width=\"0.8\" bottom_teeth_width=\"0.8\" palate_height=\"1.4\" palate_angle_deg=\"60.8\"/>\n"
    "    <p5 x=\"4.15\" z=\"-1.4\" teeth_height=\"0.7\" top_teeth_width=\"0.7\" bottom_teeth_width=\"0.7\" palate_height=\"0.7\" palate_angle_deg=\"38.0\"/>\n"
    "    <p6 x=\"4.55\" z=\"-1.1\" teeth_height=\"0.8\" top_teeth_width=\"0.65\" bottom_teeth_width=\"0.3\" palate_height=\"0.15\" palate_angle_deg=\"23.4\"/>\n"
    "    <p7 x=\"4.7\" z=\"-0.6\" teeth_height=\"0.8\" top_teeth_width=\"0.8\" bottom_teeth_width=\"0.2\" palate_height=\"0.0\" palate_angle_deg=\"0.0\"/>\n"
    "    <p8 x=\"4.7\" z=\"0.0\" teeth_height=\"0.8\" top_teeth_width=\"0.85\" bottom_teeth_width=\"0.2\" palate_height=\"0.0\" palate_angle_deg=\"0.0\"/>\n"
    "  </palate>\n"
    "  <!--****************************************************************************-->\n"
    "  <jaw fulcrum_x=\"-6.5\" fulcrum_y=\"2.0\" rest_pos_x=\"0.0\" rest_pos_y=\"-1.2\" tooth_root_length=\"0.8\">\n"
    "    <p0 x=\"0.2\" z=\"-2.3\" teeth_height=\"0.5\"  top_teeth_width=\"1.05\" bottom_teeth_width=\"1.05\" jaw_height=\"1.5\" jaw_angle_deg=\"69.5\"/>\n"
    "    <p1 x=\"1.2\" z=\"-2.2\" teeth_height=\"0.5\"  top_teeth_width=\"1.1\" bottom_teeth_width=\"1.1\" jaw_height=\"1.5\" jaw_angle_deg=\"69.5\"/>\n"
    "    <p2 x=\"2.2\" z=\"-1.9\" teeth_height=\"0.5\"  top_teeth_width=\"1.05\" bottom_teeth_width=\"1.05\" jaw_height=\"1.5\" jaw_angle_deg=\"69.5\"/>\n"
    "    <p3 x=\"3.2\" z=\"-1.6\" teeth_height=\"0.5\"  top_teeth_width=\"0.9\" bottom_teeth_width=\"0.9\" jaw_height=\"1.5\" jaw_angle_deg=\"69.5\"/>\n"
    "    <p4 x=\"3.9\" z=\"-1.4\" teeth_height=\"0.5\"  top_teeth_width=\"0.75\" bottom_teeth_width=\"0.75\" jaw_height=\"1.0\" jaw_angle_deg=\"42.2\"/>\n"
    "    <p5 x=\"4.3\" z=\"-1.1\" teeth_height=\"0.55\"  top_teeth_width=\"0.6\" bottom_teeth_width=\"0.7\" jaw_height=\"0.4\" jaw_angle_deg=\"35.8\"/>\n"
    "    <p6 x=\"4.5\" z=\"-0.7\" teeth_height=\"0.6\"  top_teeth_width=\"0.3\" bottom_teeth_width=\"0.8\" jaw_height=\"0.13\" jaw_angle_deg=\"31.4\"/>\n"
    "    <p7 x=\"4.55\" z=\"-0.5\" teeth_height=\"0.7\"  top_teeth_width=\"0.2\" bottom_teeth_width=\"0.9\" jaw_height=\"0.0\" jaw_angle_deg=\"0.0\"/>\n"
    "    <p8 x=\"4.55\" z=\"0.0\" teeth_height=\"0.7\"  top_teeth_width=\"0.2\" bottom_teeth_width=\"0.9\" jaw_height=\"0.0\" jaw_angle_deg=\"0.0\"/>\n"
    "  </jaw>\n"
    "  <!--****************************************************************************-->\n"
    "  <tongue>\n"
	  "    <tip radius=\"0.2000\"/>\n"
	  "    <body radius_x=\"1.8000\" radius_y=\"1.8000\"/>\n"
		"    <root automatic_calc=\"1\" trx_slope=\"0.938\" trx_intercept=\"-5.11\" try_slope=\"0.831\" try_intercept=\"-3.03\"/>\n"
	  "  </tongue>\n"
    "  <!--****************************************************************************-->\n"
    "  <lips width=\"1.3\"/>\n"
    "  <!--****************************************************************************-->\n"
    "  <velum velum_angle_deg=\"50.0\" uvula_width=\"0.7\" uvula_height=\"0.9\" uvula_depth=\"0.7\" max_nasal_port_area=\"2.0\" >\n"
    "    <low points=\"-1.25 -0.7 -1 -0.3 -0.7 0 -0.3 0.35 0 0.55 \"/>\n"
    "    <mid points=\"-1.75 0 -1.55 0.3 -1.15 0.5 -0.6 0.7 0 0.87 \"/>\n"
    "    <high points=\"-1.75 0.5 -1.5 0.95 -1.1 1.2 -0.6 1.25 0 1.3 \"/>\n"
    "  </velum>\n"
    "  <!--****************************************************************************-->\n"
    "  <pharynx fulcrum_x=\"-2.372\" fulcrum_y=\"1.214\" rotation_angle_deg=\"-98.0\" top_rib_y=\"-1.4\" upper_depth=\"3.8\" lower_depth=\"3.4\" back_side_width=\"1.5\"/>\n"
    "  <!--****************************************************************************-->\n"
    "  <larynx upper_depth=\"1.0\" lower_depth=\"1.0\" epiglottis_width=\"0.5\" epiglottis_height=\"1.6\" epiglottis_depth=\"1.4\" epiglottis_angle_deg=\"100.0\">\n"
    "    <narrow points=\"1.8 0 1.05 -0.2 1.55 -1.2 2.68 -3.2 1.48 -3.2 1.1 -1.2 0 -1 0 0 \"/>\n"
    "    <wide points=\"3.3 0 2.13 -0.2 1.9 -1.2 2.68 -3.2 1.48 -3.2 1.45 -1.2 0 -1 0 0 \"/>\n"  
    "  </larynx>\n"
    "  <piriform_fossa length = \"2.5\" volume = \"1.5\"/>\n"
    "  <subglottal_cavity length=\"23.0\"/>\n"
	  "  <nasal_cavity length=\"11.4\"/>\n"
    "  <!--****************************************************************************-->\n"
    "  <param index=\"0\"  name=\"HX\"   min=\"0.0\"   max=\"1.0\"   neutral=\"1.0\"   positive_velocity_factor=\"1.0\"   negative_velocity_factor=\"1.0\"/>\n"  
    "  <param index=\"1\"  name=\"HY\"   min=\"-6.0\"  max=\"-3.5\"  neutral=\"-4.75\"   positive_velocity_factor=\"1.0\"   negative_velocity_factor=\"1.0\"/>\n" 
    "  <param index=\"2\"  name=\"JX\"   min=\"-0.5\"  max=\"0.0\"   neutral=\"0.0\"   positive_velocity_factor=\"1.0\"   negative_velocity_factor=\"1.0\"/>\n"  
    "  <param index=\"3\"  name=\"JA\"   min=\"-7.0\"  max=\"0.0\"   neutral=\"-2.0\"   positive_velocity_factor=\"1.0\"   negative_velocity_factor=\"1.0\"/>\n"  
    "  <param index=\"4\"  name=\"LP\"   min=\"-1.0\"  max=\"1.0\"   neutral=\"-0.07\"   positive_velocity_factor=\"1.0\"   negative_velocity_factor=\"1.0\"/>\n" 
    "  <param index=\"5\"  name=\"LD\"   min=\"-2.0\"  max=\"4.0\"   neutral=\"0.95\"   positive_velocity_factor=\"1.0\"   negative_velocity_factor=\"1.0\"/>\n" 
    "  <param index=\"6\"  name=\"VS\"  min=\"0.0\"   max=\"1.0\"   neutral=\"0.0\"   positive_velocity_factor=\"1.0\"   negative_velocity_factor=\"1.0\"/>\n" 
    "  <param index=\"7\"  name=\"VO\"  min=\"-0.1\"   max=\"1.0\"   neutral=\"-0.1\"   positive_velocity_factor=\"1.0\"   negative_velocity_factor=\"1.0\"/>\n" 
    "  <param index=\"8\"  name=\"TCX\"  min=\"-3.0\"  max=\"4.0\"   neutral=\"-0.4\"   positive_velocity_factor=\"1.0\"   negative_velocity_factor=\"1.0\"/>\n"  
    "  <param index=\"9\"  name=\"TCY\"  min=\"-3.0\"  max=\"1.0\"   neutral=\"-1.46\"   positive_velocity_factor=\"1.0\"   negative_velocity_factor=\"1.0\"/>\n" 
    "  <param index=\"10\" name=\"TTX\"   min=\"1.5\"   max=\"5.5\"   neutral=\"3.5\"   positive_velocity_factor=\"1.0\"   negative_velocity_factor=\"1.0\"/>\n" 
    "  <param index=\"11\" name=\"TTY\"  min=\"-3.0\"  max=\"2.5\"   neutral=\"-1.0\"   positive_velocity_factor=\"1.0\"   negative_velocity_factor=\"1.0\"/>\n" 
    "  <param index=\"12\" name=\"TBX\"  min=\"-3.0\"  max=\"4.0\"   neutral=\"2.0\"   positive_velocity_factor=\"1.0\"   negative_velocity_factor=\"1.0\"/>\n"  
    "  <param index=\"13\" name=\"TBY\"  min=\"-3.0\"  max=\"5.0\"   neutral=\"0.5\"   positive_velocity_factor=\"1.0\"   negative_velocity_factor=\"1.0\"/>\n"  
    "  <param index=\"14\" name=\"TRX\"  min=\"-4.0\"  max=\"2.0\"   neutral=\"0.0\"   positive_velocity_factor=\"1.0\"   negative_velocity_factor=\"1.0\"/>\n"  
    "  <param index=\"15\" name=\"TRY\"  min=\"-6.0\"  max=\"0.0\"   neutral=\"0.0\"   positive_velocity_factor=\"1.0\"   negative_velocity_factor=\"1.0\"/>\n"  
    "  <param index=\"16\" name=\"TS1\"  min=\"0.0\" max=\"1.0\"  neutral=\"0.0\"   positive_velocity_factor=\"1.0\"   negative_velocity_factor=\"1.0\"/>\n"  
    "  <param index=\"17\" name=\"TS2\"  min=\"0.0\" max=\"1.0\"  neutral=\"0.0\"   positive_velocity_factor=\"1.0\"   negative_velocity_factor=\"1.0\"/>\n" 
    "  <param index=\"18\" name=\"TS3\"  min=\"-1.0\" max=\"1.0\"  neutral=\"0.0\"   positive_velocity_factor=\"1.0\"   negative_velocity_factor=\"1.0\"/>\n" 
    "</anatomy>";

  // ****************************************************************
  // Init. the anatomy from the above xml-string.
  // ****************************************************************

  XmlNode *node = xmlParseString(anatomyString, "anatomy");
  if (node == NULL)
  {
    printf("Fatal error: No <anatomy> node!\n");
  }

  try
  {
    readAnatomyXml(node);
  }
  catch (std::string st)
  {
    printf("Fatal error: %s\n", st.c_str());
  }

  // Delete the XML-tree.
  delete node;

  // ****************************************************************
  // Init. the long names of the parameters.
  // ****************************************************************

  param[HX].name = "Horz. hyoid pos.";
  param[HY].name = "Vert. hyoid pos.";
  param[JX].name = "Horz. jaw pos.";
  param[JA].name = "Jaw angle (deg.)";
  param[LP].name = "Lip protrusion";
  param[LD].name = "Lip distance";
  param[VS].name = "Velum shape";
  param[VO].name = "Velic opening";
  param[TCX].name = "Tongue body X";
  param[TCY].name = "Tongue body Y";
  param[TTX].name = "Tongue tip X";
  param[TTY].name = "Tongue tip Y";
  param[TBX].name = "Tongue blade X";
  param[TBY].name = "Tongue blade Y";
  param[TRX].name = "Tongue root X";
  param[TRY].name = "Tongue root Y";
  param[TS1].name = "Tongue side elevation 1";
  param[TS2].name = "Tongue side elevation 2";
  param[TS3].name = "Tongue side elevation 3";


  // ****************************************************************
  // Calculate all (center line, area function, tube function).
  // ****************************************************************

  calculateAll();

  // ****************************************************************

  hasStoredControlParams = false;
  for (i = 0; i < NUM_PARAMS; i++)
  {
    storedControlParams[i] = param[i].neutral;
  }
}


// ****************************************************************************
// Init the surfaces of the vocal tract.
// ****************************************************************************

void VocalTract::initSurfaceGrids()
{
  // ****************************************************************
  // Allocate memory for all surfaces.
  // ****************************************************************

  surface[UPPER_COVER].init(NUM_UPPER_COVER_RIBS, NUM_UPPER_COVER_POINTS);
  surface[LOWER_COVER].init(NUM_LOWER_COVER_RIBS, NUM_LOWER_COVER_POINTS);

  surface[UPPER_TEETH].init(NUM_TEETH_RIBS, NUM_TEETH_POINTS);
  surface[LOWER_TEETH].init(NUM_TEETH_RIBS, NUM_TEETH_POINTS);

  surface[UPPER_LIP].init(NUM_LIP_RIBS, NUM_LIP_POINTS);
  surface[LOWER_LIP].init(NUM_LIP_RIBS, NUM_LIP_POINTS);

  surface[TONGUE].init(NUM_TONGUE_RIBS, NUM_TONGUE_POINTS);

  surface[LEFT_COVER].init(NUM_FILL_RIBS, NUM_FILL_POINTS);
  surface[RIGHT_COVER].init(NUM_FILL_RIBS, NUM_FILL_POINTS);

  surface[EPIGLOTTIS].init(NUM_EPIGLOTTIS_RIBS, NUM_EPIGLOTTIS_POINTS);
  surface[UVULA].init(NUM_UVULA_RIBS, NUM_UVULA_POINTS);

  surface[RADIATION].init(NUM_RADIATION_RIBS, NUM_RADIATION_POINTS);

  // Two-sided grids.

  surface[UPPER_COVER_TWOSIDE].init(NUM_UPPER_COVER_RIBS, NUM_UPPER_COVER_POINTS*2-1);
  surface[LOWER_COVER_TWOSIDE].init(NUM_LOWER_COVER_RIBS, NUM_LOWER_COVER_POINTS*2-1);

  surface[UPPER_TEETH_TWOSIDE].init(NUM_TEETH_RIBS*2-1, NUM_TEETH_POINTS);
  surface[LOWER_TEETH_TWOSIDE].init(NUM_TEETH_RIBS*2-1, NUM_TEETH_POINTS);

  surface[UPPER_LIP_TWOSIDE].init(NUM_LIP_RIBS*2-1, NUM_LIP_POINTS);
  surface[LOWER_LIP_TWOSIDE].init(NUM_LIP_RIBS*2-1, NUM_LIP_POINTS);

  surface[EPIGLOTTIS_TWOSIDE].init(NUM_EPIGLOTTIS_RIBS, NUM_EPIGLOTTIS_POINTS*2-1);
  surface[UVULA_TWOSIDE].init(NUM_UVULA_RIBS, NUM_UVULA_POINTS*2-1);

  // ****************************************************************
  // Set the crease angle (just needed for rendering).
  // ****************************************************************

  surface[UPPER_COVER].creaseAngle_deg  = 170.0;
  surface[LOWER_COVER].creaseAngle_deg  = 80.0;
  surface[UPPER_TEETH].creaseAngle_deg  = 40.0;
  surface[LOWER_TEETH].creaseAngle_deg  = 40.0;
  surface[UPPER_LIP].creaseAngle_deg    = 90.0;
  surface[LOWER_LIP].creaseAngle_deg    = 90.0;
  surface[TONGUE].creaseAngle_deg       = 90.0;
  surface[LEFT_COVER].creaseAngle_deg   = 170;
  surface[RIGHT_COVER].creaseAngle_deg  = 170;
  surface[EPIGLOTTIS].creaseAngle_deg   = 170;
  surface[UVULA].creaseAngle_deg        = 170;

  surface[UPPER_COVER_TWOSIDE].creaseAngle_deg = 170.0;
  surface[LOWER_COVER_TWOSIDE].creaseAngle_deg = 80.0;
  surface[UPPER_TEETH_TWOSIDE].creaseAngle_deg = 40.0;
  surface[LOWER_TEETH_TWOSIDE].creaseAngle_deg = 40.0;
  surface[UPPER_LIP_TWOSIDE].creaseAngle_deg   = 90.0;
  surface[LOWER_LIP_TWOSIDE].creaseAngle_deg   = 90.0;
  surface[EPIGLOTTIS_TWOSIDE].creaseAngle_deg  = 170;
  surface[UVULA_TWOSIDE].creaseAngle_deg       = 170;

  // ****************************************************************
  // Swap the triangle orientation (surface normals) for some
  // surfaces that are rendered.
  // ****************************************************************

  surface[LOWER_COVER].swapTriangleOrientation();
  surface[UPPER_TEETH].swapTriangleOrientation();
  surface[UPPER_LIP].swapTriangleOrientation();
  surface[TONGUE].swapTriangleOrientation();
  surface[LEFT_COVER].swapTriangleOrientation();
  surface[EPIGLOTTIS].swapTriangleOrientation();

  surface[LOWER_COVER_TWOSIDE].swapTriangleOrientation();
  surface[UPPER_TEETH_TWOSIDE].swapTriangleOrientation();
  surface[UPPER_LIP_TWOSIDE].swapTriangleOrientation();

  surface[EPIGLOTTIS_TWOSIDE].swapTriangleOrientation();

  // ****************************************************************
  // Sufaces that are used internally only (not used for rendering).
  // ****************************************************************

  surface[NARROW_LARYNX_FRONT].init(NUM_LARYNX_RIBS, NUM_LOWER_COVER_POINTS);
  surface[NARROW_LARYNX_BACK].init(NUM_LARYNX_RIBS, NUM_UPPER_COVER_POINTS);
  surface[WIDE_LARYNX_FRONT].init(NUM_LARYNX_RIBS, NUM_LOWER_COVER_POINTS);
  surface[WIDE_LARYNX_BACK].init(NUM_LARYNX_RIBS, NUM_UPPER_COVER_POINTS);

  surface[LOWER_TEETH_ORIGINAL].init(NUM_TEETH_RIBS, NUM_TEETH_POINTS);
  surface[LOW_VELUM].init(NUM_VELUM_RIBS, NUM_UPPER_COVER_POINTS);
  surface[MID_VELUM].init(NUM_VELUM_RIBS, NUM_UPPER_COVER_POINTS);
  surface[HIGH_VELUM].init(NUM_VELUM_RIBS, NUM_UPPER_COVER_POINTS);
  surface[PALATE].init(NUM_JAW_RIBS, NUM_UPPER_COVER_POINTS);

  surface[MANDIBLE].init(NUM_JAW_RIBS, NUM_LOWER_COVER_POINTS);

  surface[EPIGLOTTIS_ORIGINAL].init(NUM_EPIGLOTTIS_RIBS, NUM_EPIGLOTTIS_POINTS);
  surface[UVULA_ORIGINAL].init(NUM_UVULA_RIBS, NUM_UVULA_POINTS);
}


// ****************************************************************************
/// Must be called after any of the anatomy parameters was changed.
// ****************************************************************************

void VocalTract::initReferenceSurfaces()
{
  initLarynx();
  initJaws();
  initVelum();
}


// ****************************************************************************
// Initialize the larynx surfaces and the epiglottis.
// ****************************************************************************

void VocalTract::initLarynx()
{
  const double EPSILON = 0.000001;
  int i, j, k;
  Point2D *L;
  double d;

  // ****************************************************************
  // Check conditions for the larynx points.
  // ****************************************************************

  bool ok = true;

  for (k=0; k < 2; k++)
  {
    if (k == 0) { L = anatomy.larynxNarrowPoints; } else { L = anatomy.larynxWidePoints; }

    if (L[0].y != L[7].y) { L[7].y = L[0].y; ok = false; }
    if (L[2].y != L[5].y) { L[5].y = L[2].y; ok = false; }
    if (L[3].y != L[4].y) { L[4].y = L[3].y; ok = false; }

    if (L[2].y < L[3].y) { L[2].y = L[5].y = L[3].y; ok = false; }
    if (L[6].y < L[5].y) { L[6].y = L[5].y; ok = false; }
    if (L[1].y < L[6].y) { L[1].y = L[6].y; ok = false; }
    if (L[0].y < L[1].y) { L[0].y = L[7].y = L[1].y; ok = false; }

    if (L[7].x > L[0].x) { L[7].x = L[0].x; ok = false; }
    if (L[5].x > L[2].x) { L[5].x = L[2].x; ok = false; }
    if (L[4].x > L[3].x) { L[4].x = L[3].x; ok = false; }
  }
  if (ok == false)
  {
    printf("Attention: The larynx points defined in the anatomy structure do not obey "
      "the required conditions! The conditions were enforced at runtime.\n");
  }

  // ****************************************************************
  // Go through all ribs.
  // ****************************************************************

  Point3D C[3];
  double w[3] = { 1.0, 0.7, 1.0 };
  BezierCurve3D curve;
  Surface *s = NULL;
  double y[NUM_LARYNX_RIBS];      // y-coord. of the current rib
  double xl[NUM_LARYNX_RIBS];     // x-coord. of the current rib on the left (posterior) contour
  double xr[NUM_LARYNX_RIBS];     // x-coord. of the current rib on the right (anterior) contour
  double xm[NUM_LARYNX_RIBS];     // x-coord. of the border between the posterior and anterior larynx part
  double zm[NUM_LARYNX_RIBS];     // leftmost z-coord. of the current rib

  for (k=0; k < 2; k++)
  {
    if (k == 0) 
    { 
      L = anatomy.larynxNarrowPoints; 
    } 
    else 
    { 
      L = anatomy.larynxWidePoints; 
    }

    y[0]  = L[3].y; 
    xl[0] = L[4].x;
    xr[0] = L[3].x;
    xm[0] = 0.5*(xl[0] + xr[0]);
    zm[0] = -0.5*anatomy.larynxLowerDepth_cm;

    y[1]  = L[2].y; 
    xl[1] = L[5].x;
    xr[1] = L[2].x;
    xm[1] = 0.5*(xl[1] + xr[1]);
    zm[1] = -0.5*anatomy.larynxUpperDepth_cm;

    y[2]  = L[6].y; 
    xl[2] = L[6].x;
    d = L[1].y - L[2].y;
    if (d < EPSILON) { d = EPSILON; }
    xr[2] = L[2].x + (L[1].x - L[2].x)*(L[6].y - L[2].y)/d;
    xm[2] = 0.25*xl[2] + 0.75*xr[2];
    zm[2] = -0.5*anatomy.pharynxLowerDepth_cm;

    y[3]  = L[1].y; 
    d = L[7].y - L[6].y;
    if (d < EPSILON) { d = EPSILON; }
    xl[3] = L[6].x + (L[7].x - L[6].x)*(L[1].y - L[6].y)/d;
    xr[3] = L[1].x;
    xm[3] = 0.2*xl[3] + 0.8*xr[3];
    zm[3] = -0.5*anatomy.pharynxLowerDepth_cm;

    y[4]  = L[0].y; 
    xl[4] = L[7].x;
    xr[4] = L[0].x;
    xm[4] = 0.5*(xl[4] + xr[4]);
    zm[4] = -0.5*anatomy.pharynxLowerDepth_cm;

    // Go through all ribs and model the ribs with quater ellipses

    for (i=0; i < NUM_LARYNX_RIBS; i++)
    {
      // Larynx back side

      if (k == 0) { s = &surface[NARROW_LARYNX_BACK]; } else { s = &surface[WIDE_LARYNX_BACK]; }

      C[0].set(xm[i], y[i], zm[i]);
      C[1].set(xl[i], y[i], zm[i]);
      C[2].set(xl[i], y[i], 0.0);
      curve.setPoints(3, C, w);

      for (j=0; j < NUM_UPPER_COVER_POINTS; j++)
      {
        d = curve.getUniformParam((double)j/(double)(NUM_UPPER_COVER_POINTS-1));
        s->setVertex(i, j, curve.getPoint(d));
      }

      // Larynx front side

      if (k == 0) { s = &surface[NARROW_LARYNX_FRONT]; } else { s = &surface[WIDE_LARYNX_FRONT]; }
      
      C[0].set(xm[i], y[i], zm[i]);
      C[1].set(xr[i], y[i], zm[i]);
      C[2].set(xr[i], y[i], 0.0);
      curve.setPoints(3, C, w);

      for (j=0; j < NUM_LOWER_COVER_POINTS; j++)
      {
        d = curve.getUniformParam((double)j/(double)(NUM_LOWER_COVER_POINTS-1));
        s->setVertex(i, j, curve.getPoint(d));
      }
    }
  }     // for k=0..1

  // ****************************************************************
  // Epiglottis.
  // ****************************************************************

  double width  = anatomy.epiglottisWidth_cm;
  double height = anatomy.epiglottisHeight_cm;
  double depth  = anatomy.epiglottisDepth_cm;
  Point3D P[NUM_EPIGLOTTIS_POINTS];

  // Ribs 0, 1 and 2
  d = 0.0;
  P[0].set(0.0, d, 0.0);
  P[1].set(-0.25*width, d, -0.75*0.5*depth);
  P[2].set(-0.5*width, d, -0.5*depth);
  P[3].set(-width, d, -0.75*0.5*depth);
  P[4].set(-width, d, 0.0);
  for (i=0; i < NUM_EPIGLOTTIS_POINTS; i++)
  {
    P[i].y = 0.0;
    surface[EPIGLOTTIS_ORIGINAL].setVertex(0, i, P[i]);
    P[i].y = 0.25*height;
    surface[EPIGLOTTIS_ORIGINAL].setVertex(1, i, P[i]);
    P[i].y = 0.75*height;
    surface[EPIGLOTTIS_ORIGINAL].setVertex(2, i, P[i]);
  }

  // Rib 3
  d = height;
  P[0].set(-0.5*width, d, 0.0);
  P[1].set(-0.5*width, d, -0.75*0.5*depth);
  P[2].set(-0.5*width, d, -0.75*0.5*depth);
  P[3].set(-0.5*width, d, -0.75*0.5*depth);
  P[4].set(-0.5*width, d, 0.0);
  for (i=0; i < NUM_EPIGLOTTIS_POINTS; i++)
  {
    surface[EPIGLOTTIS_ORIGINAL].setVertex(3, i, P[i]);
  }
}


// ****************************************************************************
// Initialize the hard palate, the mandible and the corresponding teeth rows.
// ****************************************************************************

void VocalTract::initJaws()
{
  Surface *upperJaw   = &surface[PALATE];
  Surface *lowerJaw   = &surface[MANDIBLE];
  Surface *upperTeeth = &surface[UPPER_TEETH];
  Surface *lowerTeeth = &surface[LOWER_TEETH_ORIGINAL];

  const double MIN_ANGLE_DEG = 0.000001;
  const double MAX_ANGLE_DEG = 89.999999;
  BezierCurve3D curve;
  Point3D C[3];
  double  w[3] = { 1.0, 1.0, 1.0 };
  Point3D Q;
  double angle_deg;
  double height_cm;
  int i, j;
  double u, t;
  double th[NUM_TEETH_RIBS];    // Teeth height
  double ttw[NUM_TEETH_RIBS];   // Teeth top width
  double tbw[NUM_TEETH_RIBS];   // Teeth bottom width
  Point3D tp[NUM_TEETH_RIBS];   // Point 0 of the teeth ribs (connection to jaw)
  Point3D tn[NUM_TEETH_RIBS];   // Teeth normal in the xz-plane; directed to the outside of the mouth
  double bottomX;
  double delta;

  // ****************************************************************
  // Upper and lower jaw.
  // ****************************************************************

  for (i=0; i < NUM_JAW_RIBS; i++)
  {
    // Rib of the upper jaw.

    C[0] = anatomy.palatePoints[i];
    
    height_cm = anatomy.palateHeight_cm[i];
    angle_deg = anatomy.palateAngle_deg[i];
    if (angle_deg < MIN_ANGLE_DEG) { angle_deg = MIN_ANGLE_DEG; }
    if (angle_deg > MAX_ANGLE_DEG) { angle_deg = MAX_ANGLE_DEG; }

    C[1].x = C[0].x;
    C[1].y = height_cm;
    C[1].z = C[0].z + height_cm / tan(angle_deg*M_PI/180.0);
    if (C[1].z > 0.0) { C[1].z = 0.0; }

    C[2].x = C[0].x;
    C[2].y = height_cm;
    C[2].z = 0.0;
    
    w[0] = 1.0; w[1] = 1.0; w[2] = 1.0;
    curve.setPoints(3, C, w);

    for (j=0; j < NUM_UPPER_COVER_POINTS; j++)
    {
      u = curve.getUniformParam((double)j/(double)(NUM_UPPER_COVER_POINTS-1));
      Q = curve.getPoint(u);
      upperJaw->setVertex(i, j, Q);
    }

    // Rib of the lower jaw.

    C[0] = anatomy.jawPoints[i];
    
    height_cm = anatomy.jawHeight_cm[i];
    angle_deg = anatomy.jawAngle_deg[i];
    if (angle_deg < MIN_ANGLE_DEG) { angle_deg = MIN_ANGLE_DEG; }
    if (angle_deg > MAX_ANGLE_DEG) { angle_deg = MAX_ANGLE_DEG; }

    // Make the most posterior jaw rib(s) a bit oblique.

    delta = (anatomy.jawPoints[8].x - anatomy.jawPoints[0].x)*1.0/4.5;
    if (C[0].x - anatomy.jawPoints[0].x < delta)
    {
      bottomX = anatomy.jawPoints[0].x + delta;
    }
    else
    {
      bottomX = C[0].x;
    }

    C[1].x = bottomX;
    C[1].y = -height_cm;
    C[1].z = C[0].z + height_cm / tan(angle_deg*M_PI/180.0);
    if (C[1].z > 0.0) { C[1].z = 0.0; }

    C[2].x = bottomX;
    C[2].y = -height_cm;
    C[2].z = 0.0;
    
    w[0] = 1.0; w[1] = 1.0; w[2] = 1.0;
    curve.setPoints(3, C, w);

    for (j=0; j < NUM_LOWER_COVER_POINTS; j++)
    {
      u = curve.getUniformParam((double)j/(double)(NUM_LOWER_COVER_POINTS-1));
      Q = curve.getPoint(u);
      lowerJaw->setVertex(i, j, Q);
    }
  }

  // ****************************************************************
  // Upper teeth.
  // ****************************************************************

  static const double GROOVE_DEPTH = 0.15;     // cm

  // Teeth points and teeth normals at the teeth interspaces
  tp[0]   = anatomy.palatePoints[0];
  tn[0].x = (anatomy.palatePoints[1].z - anatomy.palatePoints[0].z);
  tn[0].y = 0.0;
  tn[0].z = -(anatomy.palatePoints[1].x - anatomy.palatePoints[0].x);
  tn[0].normalize();

  tp[NUM_TEETH_RIBS-1] = anatomy.palatePoints[NUM_JAW_RIBS-1];
  tn[NUM_TEETH_RIBS-1].set(1.0, 0.0, 0.0);

  for (i=1; i < NUM_JAW_RIBS-1; i++)
  {
    j = i*3;
    tp[j]  = anatomy.palatePoints[i];
    
    tn[j].x = (anatomy.palatePoints[i+1].z - anatomy.palatePoints[i-1].z);
    tn[j].y = 0.0;
    tn[j].z = -(anatomy.palatePoints[i+1].x - anatomy.palatePoints[i-1].x);
    tn[j].normalize();
  }

  for (i=0; i < NUM_JAW_RIBS-1; i++)
  {
    t = 0.75;

    th[i*3+1]  = anatomy.upperTeethHeight_cm[i];
    ttw[i*3+1] = t*anatomy.upperTeethWidthTop_cm[i] + (1.0-t)*anatomy.upperTeethWidthTop_cm[i+1];
    tbw[i*3+1] = t*anatomy.upperTeethWidthBottom_cm[i] + (1.0-t)*anatomy.upperTeethWidthBottom_cm[i+1];
    tp[i*3+1]  = t*tp[i*3] + (1.0-t)*tp[i*3+3];
    tn[i*3+1]  = t*tn[i*3] + (1.0-t)*tn[i*3+3];
    tn[i*3+1].normalize();

    th[i*3+2]  = anatomy.upperTeethHeight_cm[i];
    ttw[i*3+2] = (1.0-t)*anatomy.upperTeethWidthTop_cm[i] + t*anatomy.upperTeethWidthTop_cm[i+1];
    tbw[i*3+2] = (1.0-t)*anatomy.upperTeethWidthBottom_cm[i] + t*anatomy.upperTeethWidthBottom_cm[i+1];
    tp[i*3+2]  = (1.0-t)*tp[i*3] + t*tp[i*3+3];
    tn[i*3+2]  = (1.0-t)*tn[i*3] + t*tn[i*3+3];
    tn[i*3+2].normalize();

    // Make the width of teeth with height=0 unvisible.
    if (anatomy.upperTeethHeight_cm[i] == 0.0)
    {
      ttw[i*3+1] = 0.0;
      tbw[i*3+1] = 0.0;
      ttw[i*3+2] = 0.0;
      tbw[i*3+2] = 0.0;
    }
  }

  // Ribs at the teeth interspaces.
  t = 0.8;
  th[0]  = th[1] - GROOVE_DEPTH;
  ttw[0] = 0.0;
  tbw[0] = 0.0;

  th[NUM_TEETH_RIBS-1]  = t*th[NUM_TEETH_RIBS-2];
  ttw[NUM_TEETH_RIBS-1] = t*ttw[NUM_TEETH_RIBS-2];
  tbw[NUM_TEETH_RIBS-1] = t*tbw[NUM_TEETH_RIBS-2];

  for (i=1; i < NUM_JAW_RIBS-1; i++)
  {
    j = i*3;

    if (th[j-1] < th[j+1])
    {
      th[j] = th[j-1] - GROOVE_DEPTH;
    }
    else
    {
      th[j] = th[j+1] - GROOVE_DEPTH;
    }
    ttw[j] = t*0.5*(ttw[j-1] + ttw[j+1]);
    tbw[j] = t*0.5*(tbw[j-1] + tbw[j+1]);
  }

  // Create the upper teeth ribs
  for (i=0; i < NUM_TEETH_RIBS; i++)
  {
    if (ttw[i] < 0.0) { ttw[i] = 0.0; }
    if (tbw[i] < 0.0) { tbw[i] = 0.0; }
    if (th[i] < 0.0) { th[i] = 0.0; }

    Q = tp[i];
    upperTeeth->setVertex(i, 0, Q);
    upperTeeth->setVertex(i, 4, Q);
    
    Q = tp[i] + tn[i]*(ttw[i] - tbw[i]);
    Q.y = -th[i];
    upperTeeth->setVertex(i, 1, Q);

    Q = tp[i] + tn[i]*ttw[i];
    Q.y = -th[i];
    upperTeeth->setVertex(i, 2, Q);
    Q.y = 0.0;
    upperTeeth->setVertex(i, 3, Q);
  }

  // Rough approximation of the upper and lower outer edge of the upper teeth.

  for (i=0; i < NUM_JAW_RIBS; i++)
  {
    upperGumsInnerEdge[i] = anatomy.palatePoints[i];
    upperGumsOuterEdge[i] = anatomy.palatePoints[i] + anatomy.upperTeethWidthTop_cm[i]*tn[3*i];
  }

  // ****************************************************************
  // Lower teeth.
  // ****************************************************************

  // Teeth points and teeth normals at the teeth interspaces.
  tp[0]  = anatomy.jawPoints[0];
  tn[0].x = (anatomy.jawPoints[1].z - anatomy.jawPoints[0].z);
  tn[0].y = 0.0;
  tn[0].z = -(anatomy.jawPoints[1].x - anatomy.jawPoints[0].x);
  tn[0].normalize();

  tp[NUM_TEETH_RIBS-1] = anatomy.jawPoints[NUM_JAW_RIBS-1];
  tn[NUM_TEETH_RIBS-1].set(1.0, 0.0, 0.0);

  for (i=1; i < NUM_JAW_RIBS-1; i++)
  {
    j = i*3;
    tp[j]  = anatomy.jawPoints[i];
    
    tn[j].x = (anatomy.jawPoints[i+1].z - anatomy.jawPoints[i-1].z);
    tn[j].y = 0.0;
    tn[j].z = -(anatomy.jawPoints[i+1].x - anatomy.jawPoints[i-1].x);
    tn[j].normalize();
  }

  for (i=0; i < NUM_JAW_RIBS-1; i++)
  {
    t = 0.75;

    th[i*3+1]  = anatomy.lowerTeethHeight_cm[i];
    ttw[i*3+1] = t*anatomy.lowerTeethWidthTop_cm[i] + (1.0-t)*anatomy.lowerTeethWidthTop_cm[i+1];
    tbw[i*3+1] = t*anatomy.lowerTeethWidthBottom_cm[i] + (1.0-t)*anatomy.lowerTeethWidthBottom_cm[i+1];
    tp[i*3+1]  = t*tp[i*3] + (1.0-t)*tp[i*3+3];
    tn[i*3+1]  = t*tn[i*3] + (1.0-t)*tn[i*3+3];
    tn[i*3+1].normalize();

    th[i*3+2]  = anatomy.lowerTeethHeight_cm[i];
    ttw[i*3+2] = (1.0-t)*anatomy.lowerTeethWidthTop_cm[i] + t*anatomy.lowerTeethWidthTop_cm[i+1];
    tbw[i*3+2] = (1.0-t)*anatomy.lowerTeethWidthBottom_cm[i] + t*anatomy.lowerTeethWidthBottom_cm[i+1];
    tp[i*3+2]  = (1.0-t)*tp[i*3] + t*tp[i*3+3];
    tn[i*3+2]  = (1.0-t)*tn[i*3] + t*tn[i*3+3];
    tn[i*3+2].normalize();

    // Make the width of teeth with height=0 unvisible.

    if (anatomy.lowerTeethHeight_cm[i] == 0.0)
    {
      ttw[i*3+1] = 0.0;
      tbw[i*3+1] = 0.0;
      ttw[i*3+2] = 0.0;
      tbw[i*3+2] = 0.0;
    }
  }

  // Ribs at the teeth interspaces.

  t = 0.8;
  th[0]  = th[1] - GROOVE_DEPTH;
  ttw[0] = 0.0;
  tbw[0] = 0.0;

  th[NUM_TEETH_RIBS-1]  = t*th[NUM_TEETH_RIBS-2];
  ttw[NUM_TEETH_RIBS-1] = t*ttw[NUM_TEETH_RIBS-2];
  tbw[NUM_TEETH_RIBS-1] = t*tbw[NUM_TEETH_RIBS-2];

  for (i=1; i < NUM_JAW_RIBS-1; i++)
  {
    j = i*3;
    if (th[j-1] < th[j+1])
    {
      th[j] = th[j-1] - GROOVE_DEPTH;
    }
    else
    {
      th[j] = th[j+1] - GROOVE_DEPTH;
    }
    ttw[j] = t*0.5*(ttw[j-1] + ttw[j+1]);
    tbw[j] = t*0.5*(tbw[j-1] + tbw[j+1]);
  }

  // Create the lower teeth ribs.

  for (i=0; i < NUM_TEETH_RIBS; i++)
  {
    if (ttw[i] < 0.0) { ttw[i] = 0.0; }
    if (tbw[i] < 0.0) { tbw[i] = 0.0; }
    if (th[i] < 0.0) { th[i] = 0.0; }

    Q = tp[i];
    lowerTeeth->setVertex(i, 0, Q);
    lowerTeeth->setVertex(i, 4, Q);
    
    Q = tp[i] + tn[i]*(tbw[i] - ttw[i]);
    Q.y = th[i];
    lowerTeeth->setVertex(i, 1, Q);

    Q = tp[i] + tn[i]*tbw[i];
    Q.y = th[i];
    lowerTeeth->setVertex(i, 2, Q);
    Q.y = 0.0;
    lowerTeeth->setVertex(i, 3, Q);
  }

  // Rough approximation of the upper and lower outer edge of the lower teeth.

  for (i=0; i < NUM_JAW_RIBS; i++)
  {
    lowerGumsInnerEdgeOrig[i] = anatomy.jawPoints[i];
    lowerGumsOuterEdgeOrig[i] = anatomy.jawPoints[i] + anatomy.lowerTeethWidthBottom_cm[i]*tn[3*i];
  }

  // ****************************************************************
  // Determination of the line, along which the mouth corner is 
  // allowed to move.
  // ****************************************************************

  wideLipCornerPath.reset(0);
  narrowLipCornerPath.reset(0);

  // Add outer teeth edge points from the posterior side of the second
  // premolar to the gap between the outer incisor and the corner tooth

  for (i=4; i <= 6; i++)
  {
    Q = upperGumsOuterEdge[i];
    Q.y = 0.0;
    wideLipCornerPath.addPoint(Q);
    narrowLipCornerPath.addPoint(Q);
  }

  // 1 cm in front of the upper incisors
  double palateLength = anatomy.palatePoints[8].x - anatomy.palatePoints[0].x;
  double palateDepth  = -2.0*anatomy.palatePoints[0].z;
  
  double lipProtrusion = 1.0*palateLength / 4.5;    // = 1 cm for a palate length of 4.5 cm
  double lipFrontMinZ  = -1.0*palateDepth / 4.6;    // Reference palate depth = 4.6 cm
  double lipFrontMaxZ  = -0.4*palateDepth / 4.6;

  double maxLipCornerX = upperGumsOuterEdge[8].x + lipProtrusion;

  wideLipCornerPath.addPoint(Point3D(maxLipCornerX, 0.0, lipFrontMinZ));
  narrowLipCornerPath.addPoint(Point3D(maxLipCornerX, 0.0, lipFrontMaxZ));
}


// ****************************************************************************
// Initialize the three states of the velum.
// ****************************************************************************

void VocalTract::initVelum()
{
  int i, j;
  double t, u;
  Point3D C[3];       // Control point of the Bezier-curve for a rib
  double w[3] = { 1.0, 1.0, 1.0 };
  BezierCurve3D curve;
  Point3D Q;
  Point3D posteriorEndPoint, anteriorEndPoint;

  // The z-coord. of the middle point of the square splines that form 
  // the ribs is linearly interpolated between sourceZ and targetZ !

  double sourceZ = -0.5*anatomy.pharynxUpperDepth_cm;
  t = anatomy.palateAngle_deg[0];
  if (t <= 0.00001) { t = 0.00001; }
  if (t > 89.99999) { t = 89.99999; }
  double targetZ = anatomy.palatePoints[0].z + anatomy.palateHeight_cm[1]/tan(t*M_PI/180.0);

  // The posterior (where the velum and the pharynx back wall meet) and 
  // anterior (on the posterior maxillary plane; x=0) lower end points.

  posteriorEndPoint.x = getPharynxBackX(anatomy.pharynxTopRibY_cm) + anatomy.pharynxBackWidth_cm;
  posteriorEndPoint.y = anatomy.pharynxTopRibY_cm;
  posteriorEndPoint.z = -0.5*anatomy.pharynxUpperDepth_cm;

  anteriorEndPoint.x = 0.0;
  anteriorEndPoint.y = 0.0;
  anteriorEndPoint.z = anatomy.palatePoints[0].z;

  // ****************************************************************
  // Go through all velum ribs.
  // ****************************************************************

  for (i=0; i < NUM_VELUM_RIBS; i++)
  {
    if (i == 0) 
    { 
      C[0] = posteriorEndPoint; 
    }
    else
    {
      t = (double)(i-1) / (double)(NUM_VELUM_RIBS-2);
      C[0] = (1.0-t)*posteriorEndPoint + t*anteriorEndPoint;
    }

    t = (double)i / (double)(NUM_VELUM_RIBS-1);

    // **************************************************************
    // Low velum state.
    // **************************************************************

    if (i == 0)
    {
      C[2].x = getPharynxBackX(anatomy.velumLowPoints[0].y);
      C[2].y = anatomy.velumLowPoints[0].y;
    }
    else
    {
      C[2] = anatomy.velumLowPoints[i-1].toPoint3D();
    }
    C[2].z = 0.0;

    C[1].x = C[2].x;
    C[1].y = C[2].y;
    C[1].z = (1.0-t)*sourceZ + t*targetZ;

    curve.setPoints(3, C, w);

    for (j=0; j < NUM_UPPER_COVER_POINTS; j++)
    {
      u = (double)j/(double)(NUM_UPPER_COVER_POINTS-1);
      u = curve.getUniformParam(u);
      Q = curve.getPoint(u);
      surface[LOW_VELUM].setVertex(i, j, Q);
    }

    // **************************************************************
    // Mid velum state.
    // **************************************************************

    if (i == 0)
    {
      C[2].x = getPharynxBackX(anatomy.velumMidPoints[0].y);
      C[2].y = anatomy.velumMidPoints[0].y;
    }
    else
    {
      C[2] = anatomy.velumMidPoints[i-1].toPoint3D();
    }
    C[2].z = 0.0;

    C[1].x = C[2].x;
    C[1].y = C[2].y;
    C[1].z = (1.0-t)*sourceZ + t*targetZ;

    curve.setPoints(3, C, w);

    for (j=0; j < NUM_UPPER_COVER_POINTS; j++)
    {
      u = (double)j/(double)(NUM_UPPER_COVER_POINTS-1);
      u = curve.getUniformParam(u);
      Q = curve.getPoint(u);
      surface[MID_VELUM].setVertex(i, j, Q);
    }

    // **************************************************************
    // High velum state.
    // **************************************************************

    if (i == 0)
    {
      C[2].x = getPharynxBackX(anatomy.velumHighPoints[0].y);
      C[2].y = anatomy.velumHighPoints[0].y;
    }
    else
    {
      C[2] = anatomy.velumHighPoints[i-1].toPoint3D();
    }
    C[2].z = 0.0;

    C[1].x = C[2].x;
    C[1].y = C[2].y;
    C[1].z = (1.0-t)*sourceZ + t*targetZ;

    curve.setPoints(3, C, w);

    for (j=0; j < NUM_UPPER_COVER_POINTS; j++)
    {
      u = (double)j/(double)(NUM_UPPER_COVER_POINTS-1);
      u = curve.getUniformParam(u);
      Q = curve.getPoint(u);
      surface[HIGH_VELUM].setVertex(i, j, Q);
    }
  }

  // ****************************************************************
  // Uvula.
  // ****************************************************************

  double width  = anatomy.uvulaWidth_cm;
  double height = anatomy.uvulaHeight_cm;
  double depth  = anatomy.uvulaDepth_cm;
  Point3D P;

  // Ribs 0 and 1
  for (i=0; i < NUM_UVULA_POINTS; i++)
  {
    t = M_PI*(double)i / (double)(NUM_UVULA_POINTS-1);
    P.set(0.5*width*cos(t) - 0.5*width, 0.0, -0.5*depth*sin(t));
    surface[UVULA_ORIGINAL].setVertex(0, i, P);
    P.y = -0.5*height;
    surface[UVULA_ORIGINAL].setVertex(1, i, P);
  }

  // Rib 2
  for (i=0; i < NUM_UVULA_POINTS; i++)
  {
    t = M_PI*(double)i / (double)(NUM_UVULA_POINTS-1);
    P.set(0.75*0.5*width*cos(t) - 0.5*width, -0.75*height, -0.75*0.5*depth*sin(t));
    surface[UVULA_ORIGINAL].setVertex(2, i, P);
  }

  // Rib 3
  for (i=0; i < NUM_UVULA_POINTS; i++)
  {
    P.set(-0.5*width, -height, 0.0);
    surface[UVULA_ORIGINAL].setVertex(3, i, P);
  }
}


// ****************************************************************************
/// Init EMA points at default positions.
// ****************************************************************************

void VocalTract::setDefaultEmaPoints()
{
  EmaPoint p;

  emaPoints.clear();

  p.name = "TB";
  p.emaSurface = EMA_SURFACE_TONGUE;
  p.vertexIndex = 10;
  emaPoints.push_back(p);

  p.name = "TM";
  p.emaSurface = EMA_SURFACE_TONGUE;
  p.vertexIndex = 20;
  emaPoints.push_back(p);

  p.name = "TT";
  p.emaSurface = EMA_SURFACE_TONGUE;
  p.vertexIndex = 30;
  emaPoints.push_back(p);

  p.name = "UL";
  p.emaSurface = EMA_SURFACE_UPPER_LIP;
  p.vertexIndex = NUM_LIP_POINTS - 1;
  emaPoints.push_back(p);

  p.name = "LL";
  p.emaSurface = EMA_SURFACE_LOWER_LIP;
  p.vertexIndex = NUM_LIP_POINTS - 1;
  emaPoints.push_back(p);

  p.name = "JAW";
  p.emaSurface = EMA_SURFACE_LOWER_COVER;
  p.vertexIndex = NUM_LOWER_COVER_RIBS - 1;
  emaPoints.push_back(p);
}


// ****************************************************************************
/// Returns the 3D-coord. of the given EMA point.
// ****************************************************************************

Point3D VocalTract::getEmaPointCoord(int index)
{
  Point3D P(0.0, 0.0, 0.0);

  if ((index < 0) || (index >= (int)emaPoints.size()))
  {
    return P;
  }

  EmaPoint *e = &emaPoints[index];
  Surface *s = NULL;
  int i = e->vertexIndex;

  if (e->emaSurface == EMA_SURFACE_TONGUE)
  {
    s = &surface[TONGUE];
    if (i < 0) { i = 0; }
    if (i >= s->numRibs) { i = s->numRibs-1; }
    P = s->getVertex(i, s->numRibPoints/2);
  }
  else
  if (e->emaSurface == EMA_SURFACE_UPPER_COVER)
  {
    s = &surface[UPPER_COVER_TWOSIDE];
    if (i < 0) { i = 0; }
    if (i >= s->numRibs) { i = s->numRibs-1; }
    P = s->getVertex(i, s->numRibPoints/2);
  }
  else
  if (e->emaSurface == EMA_SURFACE_LOWER_COVER)
  {
    s = &surface[LOWER_COVER_TWOSIDE];
    if (i < 0) { i = 0; }
    if (i >= s->numRibs) { i = s->numRibs-1; }
    P = s->getVertex(i, s->numRibPoints/2);
  }
  else
  if (e->emaSurface == EMA_SURFACE_UPPER_LIP)
  {
    s = &surface[UPPER_LIP_TWOSIDE];
    if (i < 0) { i = 0; }
    if (i >= s->numRibPoints) { i = s->numRibPoints-1; }
    P = s->getVertex(s->numRibs/2, i);
  }
  else
  if (e->emaSurface == EMA_SURFACE_LOWER_LIP)
  {
    s = &surface[LOWER_LIP_TWOSIDE];
    if (i < 0) { i = 0; }
    if (i >= s->numRibPoints) { i = s->numRibPoints-1; }
    P = s->getVertex(s->numRibs/2, i);
  }

  return P;
}


// ****************************************************************************
/// Returns the minimum and maximum vertex index for the given EMA surface.
// ****************************************************************************

void VocalTract::getEmaSurfaceVertexRange(int emaSurface, int *min, int *max)
{
  if ((emaSurface < 0) || (emaSurface >= NUM_EMA_SURFACES) || 
    (min == NULL) || (max == NULL))
  {
    return;
  }

  // Default setting.
  *min = 0;
  *max = 0;

  Surface *s = NULL;

  if (emaSurface == EMA_SURFACE_TONGUE)
  {
    s = &surface[TONGUE];
    *min = 0;
    *max = s->numRibs-1;
  }
  else
  if (emaSurface == EMA_SURFACE_UPPER_COVER)
  {
    s = &surface[UPPER_COVER_TWOSIDE];
    *min = 0;
    *max = s->numRibs-1;
  }
  else
  if (emaSurface == EMA_SURFACE_LOWER_COVER)
  {
    s = &surface[LOWER_COVER_TWOSIDE];
    *min = 0;
    *max = s->numRibs-1;
  }
  else
  if (emaSurface == EMA_SURFACE_UPPER_LIP)
  {
    s = &surface[UPPER_LIP_TWOSIDE];
    *min = 0;
    *max = s->numRibPoints-1;
  }
  else
  if (emaSurface == EMA_SURFACE_LOWER_LIP)
  {
    s = &surface[LOWER_LIP_TWOSIDE];
    *min = 0;
    *max = s->numRibPoints-1;
  }
}


// ****************************************************************************
// Reads the anatomy data for the vocal tract from the <anatomy> element node.
// ****************************************************************************

void VocalTract::readAnatomyXml(XmlNode *anatomyNode)
{
  XmlNode *sectionNode;
  XmlNode *node;
  int i, k;
  char st[1024];
  string numberString;
  istringstream is;

  try
  {
    // **************************************************************
    // Palate.
    // **************************************************************

    sectionNode = XmlHelper::getChildNode(anatomyNode, "palate");

    for (i=0; i < NUM_PALATE_RIBS; i++)
    {
      sprintf(st, "p%d", i);
      node = XmlHelper::getChildNode(sectionNode, st);

      XmlHelper::readAttribute(node, "x", anatomy.palatePoints[i].x);
      anatomy.palatePoints[i].y = 0.0;
      XmlHelper::readAttribute(node, "z", anatomy.palatePoints[i].z);

      XmlHelper::readAttribute(node, "teeth_height", anatomy.upperTeethHeight_cm[i]);
      XmlHelper::readAttribute(node, "top_teeth_width", anatomy.upperTeethWidthTop_cm[i]);
      XmlHelper::readAttribute(node, "bottom_teeth_width", anatomy.upperTeethWidthBottom_cm[i]);
      XmlHelper::readAttribute(node, "palate_height", anatomy.palateHeight_cm[i]);
      XmlHelper::readAttribute(node, "palate_angle_deg", anatomy.palateAngle_deg[i]);
    }

    // **************************************************************
    // Jaw.
    // **************************************************************

    sectionNode = XmlHelper::getChildNode(anatomyNode, "jaw");
    
    XmlHelper::readAttribute(sectionNode, "fulcrum_x", anatomy.jawFulcrum.x);
    XmlHelper::readAttribute(sectionNode, "fulcrum_y", anatomy.jawFulcrum.y);
    XmlHelper::readAttribute(sectionNode, "rest_pos_x", anatomy.jawRestPos.x);
    XmlHelper::readAttribute(sectionNode, "rest_pos_y", anatomy.jawRestPos.y);
    XmlHelper::readAttribute(sectionNode, "tooth_root_length", anatomy.toothRootLength_cm);

    for (i=0; i < NUM_JAW_RIBS; i++)
    {
      sprintf(st, "p%d", i);
      node = XmlHelper::getChildNode(sectionNode, st);

      XmlHelper::readAttribute(node, "x", anatomy.jawPoints[i].x);
      anatomy.jawPoints[i].y = 0.0;
      XmlHelper::readAttribute(node, "z", anatomy.jawPoints[i].z);

      XmlHelper::readAttribute(node, "teeth_height", anatomy.lowerTeethHeight_cm[i]);
      XmlHelper::readAttribute(node, "top_teeth_width", anatomy.lowerTeethWidthTop_cm[i]);
      XmlHelper::readAttribute(node, "bottom_teeth_width", anatomy.lowerTeethWidthBottom_cm[i]);
      XmlHelper::readAttribute(node, "jaw_height", anatomy.jawHeight_cm[i]);
      XmlHelper::readAttribute(node, "jaw_angle_deg", anatomy.jawAngle_deg[i]);
    }

    // **************************************************************
    // Tongue.
    // **************************************************************

    sectionNode = XmlHelper::getChildNode(anatomyNode, "tongue");

    node = XmlHelper::getChildNode(sectionNode, "tip");
    XmlHelper::readAttribute(node, "radius", anatomy.tongueTipRadius_cm);

    node = XmlHelper::getChildNode(sectionNode, "body");
    XmlHelper::readAttribute(node, "radius_x", anatomy.tongueCenterRadiusX_cm);
    XmlHelper::readAttribute(node, "radius_y", anatomy.tongueCenterRadiusY_cm);

    node = XmlHelper::getChildNode(sectionNode, "root");
    XmlHelper::readAttribute(node, "automatic_calc", k);
    if (k == 0)
    {
      anatomy.automaticTongueRootCalc = false;
    }
    else
    {
      anatomy.automaticTongueRootCalc = true;
    }
    
    XmlHelper::readAttribute(node, "trx_slope", anatomy.tongueRootTrxSlope);
    XmlHelper::readAttribute(node, "trx_intercept", anatomy.tongueRootTrxIntercept);
    XmlHelper::readAttribute(node, "try_slope", anatomy.tongueRootTrySlope);
    XmlHelper::readAttribute(node, "try_intercept", anatomy.tongueRootTryIntercept);

    // **************************************************************
    // Lips.
    // **************************************************************

    sectionNode = XmlHelper::getChildNode(anatomyNode, "lips");

    XmlHelper::readAttribute(sectionNode, "width", anatomy.lipsWidth_cm);

    // **************************************************************
    // Velum.
    // **************************************************************

    sectionNode = XmlHelper::getChildNode(anatomyNode, "velum");

    XmlHelper::readAttribute(sectionNode, "uvula_width", anatomy.uvulaWidth_cm);
    XmlHelper::readAttribute(sectionNode, "uvula_height", anatomy.uvulaHeight_cm);
    XmlHelper::readAttribute(sectionNode, "uvula_depth", anatomy.uvulaDepth_cm);
    XmlHelper::readAttribute(sectionNode, "max_nasal_port_area", anatomy.maxNasalPortArea_cm2);

    // **************************************************************
    // Low velum points.
    // **************************************************************

    node = XmlHelper::getChildNode(sectionNode, "low");
    XmlHelper::readAttribute(node, "points", numberString);

    is.clear();     // clear eof-flag
    is.str(numberString);
    for (i=0; i < 5; i++)   //NUM_VELUM_RIBS-1
    {
      is >> anatomy.velumLowPoints[i].x;
      if (is.eof()) 
      { 
        throw std::string("There are not enough point coordinates in the 'points' attribute of the <") + 
          node->name.c_str() + "> tag!"; 
      }
      is >> anatomy.velumLowPoints[i].y;
    }

    // **************************************************************
    // Mid velum points.
    // **************************************************************
    
    node = XmlHelper::getChildNode(sectionNode, "mid");
    XmlHelper::readAttribute(node, "points", numberString);

    is.clear();     // clear eof-flag
    is.str(numberString);
    for (i=0; i < 5; i++)   //NUM_VELUM_RIBS-1
    {
      is >> anatomy.velumMidPoints[i].x;
      if (is.eof()) 
      { 
        throw std::string("There are not enough point coordinates in the 'points' attribute of the <") + 
          node->name.c_str() + "> tag!"; 
      }
      is >> anatomy.velumMidPoints[i].y;
    }

    // **************************************************************
    // High velum points.
    // **************************************************************
    
    node = XmlHelper::getChildNode(sectionNode, "high");
    XmlHelper::readAttribute(node, "points", numberString);

    is.clear();     // clear eof-flag
    is.str(numberString);
    for (i=0; i < 5; i++)   //NUM_VELUM_RIBS-1
    {
      is >> anatomy.velumHighPoints[i].x;
      if (is.eof()) 
      { 
        throw std::string("There are not enough point coordinates in the 'points' attribute of the <") + 
          node->name.c_str() + "> tag!"; 
      }
      is >> anatomy.velumHighPoints[i].y;
    }

    // **************************************************************
    // Pharynx.
    // **************************************************************

    sectionNode = XmlHelper::getChildNode(anatomyNode, "pharynx");

    XmlHelper::readAttribute(sectionNode, "fulcrum_x", anatomy.pharynxFulcrum.x);
    XmlHelper::readAttribute(sectionNode, "fulcrum_y", anatomy.pharynxFulcrum.y);
    XmlHelper::readAttribute(sectionNode, "rotation_angle_deg", anatomy.pharynxRotationAngle_deg);
    XmlHelper::readAttribute(sectionNode, "top_rib_y", anatomy.pharynxTopRibY_cm);
    XmlHelper::readAttribute(sectionNode, "upper_depth", anatomy.pharynxUpperDepth_cm);
    XmlHelper::readAttribute(sectionNode, "lower_depth", anatomy.pharynxLowerDepth_cm);
    XmlHelper::readAttribute(sectionNode, "back_side_width", anatomy.pharynxBackWidth_cm);

    // **************************************************************
    // Larynx.
    // **************************************************************

    sectionNode = XmlHelper::getChildNode(anatomyNode, "larynx");

    XmlHelper::readAttribute(sectionNode, "upper_depth", anatomy.larynxUpperDepth_cm);
    XmlHelper::readAttribute(sectionNode, "lower_depth", anatomy.larynxLowerDepth_cm);
    XmlHelper::readAttribute(sectionNode, "epiglottis_width", anatomy.epiglottisWidth_cm);
    XmlHelper::readAttribute(sectionNode, "epiglottis_height", anatomy.epiglottisHeight_cm);
    XmlHelper::readAttribute(sectionNode, "epiglottis_depth", anatomy.epiglottisDepth_cm);
    XmlHelper::readAttribute(sectionNode, "epiglottis_angle_deg", anatomy.epiglottisAngle_deg);

    // **************************************************************
    // Narrow larynx points.
    // **************************************************************

    node = XmlHelper::getChildNode(sectionNode, "narrow");
    XmlHelper::readAttribute(node, "points", numberString);

    is.clear();     // clear eof-flag
    is.str(numberString);
    for (i=0; i < 8; i++)
    {
      is >> anatomy.larynxNarrowPoints[i].x;
      if (is.eof()) 
      { 
        throw std::string("There are not enough point coordinates in the 'points' attribute of the <") +
          node->name.c_str() + "> tag!"; 
      }
      is >> anatomy.larynxNarrowPoints[i].y;
    }

    // **************************************************************
    // Wide larynx points.
    // **************************************************************

    node = XmlHelper::getChildNode(sectionNode, "wide");
    XmlHelper::readAttribute(node, "points", numberString);

    is.clear();     // clear eof-flag
    is.str(numberString);
    for (i=0; i < 8; i++)
    {
      is >> anatomy.larynxWidePoints[i].x;
      if (is.eof()) 
      { 
        throw std::string("There are not enough point coordinates in the 'points' attribute of the <") +
          node->name.c_str() + "> tag!"; 
      }
      is >> anatomy.larynxWidePoints[i].y;
    }

    // **************************************************************
    // Piriform fossa.
    // **************************************************************

    sectionNode = XmlHelper::getChildNode(anatomyNode, "piriform_fossa");

    XmlHelper::readAttribute(sectionNode, "length", anatomy.piriformFossaLength_cm);
    XmlHelper::readAttribute(sectionNode, "volume", anatomy.piriformFossaVolume_cm3);

    // **************************************************************
    // Subglottal cavity.
    // **************************************************************

    sectionNode = XmlHelper::getChildNode(anatomyNode, "subglottal_cavity");

    XmlHelper::readAttribute(sectionNode, "length", anatomy.subglottalCavityLength_cm);

    // **************************************************************
    // Nasal cavity.
    // **************************************************************

    sectionNode = XmlHelper::getChildNode(anatomyNode, "nasal_cavity");

    XmlHelper::readAttribute(sectionNode, "length", anatomy.nasalCavityLength_cm);

    // ****************************************************************
    // Read the parameter short names, min, max, and neutral values,
    // as well as the intrinsic velocity factors.
    // ****************************************************************

    int numParamNodes = anatomyNode->numChildElements("param");
    int m;
    bool paramRead[NUM_PARAMS];

    for (i=0; i < NUM_PARAMS; i++) { paramRead[i] = false; }

    // Go through all parameter definitions in the xml-file.

    for (i=0; i < numParamNodes; i++)
    {
      node = anatomyNode->getChildElement("param", i);
      if (node == NULL) 
      { 
        throw std::string("One of the <param> tags is empty!"); 
      }

      XmlHelper::readAttribute(node, "index", m);
      if ((m < 0) || (m >= NUM_PARAMS)) 
      { 
        throw std::string("Invalid index in the element <") + 
          node->name.c_str() + "> (out of range)"; 
      }
      paramRead[m] = true;    // Mark this parameter as read/initialized

      XmlHelper::readAttribute(node, "name", param[m].abbr);
      XmlHelper::readAttribute(node, "min", param[m].min);
      XmlHelper::readAttribute(node, "max", param[m].max);
      XmlHelper::readAttribute(node, "neutral", param[m].neutral);
      XmlHelper::readAttribute(node, "positive_velocity_factor", anatomy.positiveVelocityFactor[m]);
      XmlHelper::readAttribute(node, "negative_velocity_factor", anatomy.negativeVelocityFactor[m]);

      // Set the current parameter values to their neutral values
      
      param[m].x = param[m].neutral;
      param[m].limitedX = param[m].neutral;
    }

    // All parameters must have been read !!! 

    bool ok = true;
    for (i=0; i < NUM_PARAMS; i++) 
    { 
      if (paramRead[i] == false) { ok = false; }
    }
    if (ok == false)
    {
      throw "Fatal error: Not all of the vocal tract parameters could be initialized from the Xml file!";
    }

  }
  // ****************************************************************
  // Catch all exceptions that might occur during the reading
  // ****************************************************************

  catch (std::string st)
  {
    throw;
  }
  catch (...)
  {
    throw string("Unknown error occured in readAnatomyNode().");
  }

  // ****************************************************************
  // Pre-calculate some surfaces by means of the anatomy parameters.
  // ****************************************************************

  initReferenceSurfaces();
}


// ****************************************************************************
/// Reads the list of vocal tract shapes from the <shapes> element node.
// ****************************************************************************

void VocalTract::readShapesXml(XmlNode *shapeListNode)
{
  int i, k, m, n;
  XmlNode *shapeNode;
  XmlNode *paramNode;
  string paramName;
  Shape s;
  
  // Clear the current shape list.
  shapes.clear();

  // ****************************************************************
  // Run through all shapes.
  // ****************************************************************

  int numShapes = shapeListNode->numChildElements("shape");
  int numParams = 0;

  for (i=0; i < numShapes; i++)
  {
    // Assign neutral parameter and dominance values for initialization.

    for (k=0; k < NUM_PARAMS; k++)
    {
      s.param[k] = param[k].neutral;
    }

    // Read the xml data for the new shape.

    shapeNode = shapeListNode->getChildElement("shape", i);

    if (shapeNode->hasAttribute("name"))
    {
      s.name = shapeNode->getAttributeString("name");
    }
    else
    {
      s.name = "--";
    }

    // Run through all parameters of the current shape.

    numParams = shapeNode->numChildElements("param");
    for (k=0; k < numParams; k++)
    {
      paramNode = shapeNode->getChildElement("param", k);
      paramName = paramNode->getAttributeString("name");

      n = -1;
      for (m=0; (m < NUM_PARAMS) && (n == -1); m++)
      {
        if (paramName == param[m].abbr) 
        { 
          n = m; 
        }
      }
      if (n != -1)
      {
        s.param[n] = paramNode->getAttributeDouble("value");
      }
    }

    // Append the new shape to the shape list.
    shapes.push_back(s);
  }

}


// ****************************************************************************
/// Read the speaker anatomy and vocal tract shape list from an xml file.
// ****************************************************************************

void VocalTract::readFromXml(const string &speakerFileName)
{
  XmlNode *rootNode = xmlParseFile(speakerFileName, "speaker");
  if (rootNode == NULL)
  {
    throw std::string("Error parsing the file ") + speakerFileName + ".";
  }

  XmlNode *vocalTractNode = rootNode->getChildElement("vocal_tract_model");
  if (vocalTractNode == NULL)
  {
    throw std::string("The file ") + speakerFileName + " has no <vocal_tract_model> element!";
  }

  // ****************************************************************
  // Read the <anatomy> part of the vocal tract model.
  // ****************************************************************

  XmlNode *anatomyNode = vocalTractNode->getChildElement("anatomy");
  if (anatomyNode == NULL)
  {
    throw std::string("The file ") + speakerFileName + " has no <anatomy> element!";
  }
  try
  {
    readAnatomyXml(anatomyNode);
  }
  catch (std::string st)
  {
    throw;
  }

  // ****************************************************************
  // Read the <shapes> part of the vocal tract model.
  // ****************************************************************

  XmlNode *shapeListNode = vocalTractNode->getChildElement("shapes");
  if (shapeListNode == NULL)
  {
    throw std::string("The file ") + speakerFileName + " has no <shapes> element!";
  }
  try
  {
    readShapesXml(shapeListNode);
  }
  catch (std::string st)
  {
    throw;
  }

  // Delete the XML-tree.

  delete rootNode;
}


// ****************************************************************************
// Write the anatomy parameters in xml-format into the ouput stream os.
// ****************************************************************************

void VocalTract::writeAnatomyXml(ostream &os, int indent)
{
  int i;
  char st[1024];

  os << setprecision(4);    // floating point precision
  os << string(indent, ' ') << "<anatomy>" << endl;
  indent+= 2;

  // ****************************************************************
  // Palate data.
  // ****************************************************************
  
  os << string(indent, ' ') << "<palate>" << endl;
  indent+= 2;

  for (i=0; i < NUM_PALATE_RIBS; i++)
  {
    sprintf(st, "<p%d x=\"%2.4f\" z=\"%2.4f\" teeth_height=\"%2.4f\" top_teeth_width=\"%2.4f\" "
      "bottom_teeth_width=\"%2.4f\" palate_height=\"%2.4f\" palate_angle_deg=\"%2.4f\"/>", 
      i, anatomy.palatePoints[i].x, anatomy.palatePoints[i].z, anatomy.upperTeethHeight_cm[i],
      anatomy.upperTeethWidthTop_cm[i], anatomy.upperTeethWidthBottom_cm[i], anatomy.palateHeight_cm[i],
      anatomy.palateAngle_deg[i]);
    
    os << string(indent, ' ') << string(st) << endl;
  }

  indent-= 2;
  os << string(indent, ' ') << "</palate>" << endl;

  // ****************************************************************
  // Jaw data.
  // ****************************************************************

  sprintf(st, "<jaw fulcrum_x=\"%2.4f\" fulcrum_y=\"%2.4f\" rest_pos_x=\"%2.4f\" rest_pos_y=\"%2.4f\" tooth_root_length=\"%2.4f\">", 
    anatomy.jawFulcrum.x, anatomy.jawFulcrum.y,
    anatomy.jawRestPos.x, anatomy.jawRestPos.y,
    anatomy.toothRootLength_cm);
  os << string(indent, ' ') << string(st) << endl;
  indent+= 2;

  for (i=0; i < NUM_PALATE_RIBS; i++)
  {
    sprintf(st, "<p%d x=\"%2.4f\" z=\"%2.4f\" teeth_height=\"%2.4f\" top_teeth_width=\"%2.4f\" "
      "bottom_teeth_width=\"%2.4f\" jaw_height=\"%2.4f\" jaw_angle_deg=\"%2.4f\"/>", 
      i, anatomy.jawPoints[i].x, anatomy.jawPoints[i].z, anatomy.lowerTeethHeight_cm[i],
      anatomy.lowerTeethWidthTop_cm[i], anatomy.lowerTeethWidthBottom_cm[i], anatomy.jawHeight_cm[i],
      anatomy.jawAngle_deg[i]);
    
    os << string(indent, ' ') << string(st) << endl;
  }

  indent-= 2;
  os << string(indent, ' ') << "</jaw>" << endl;
  
  // ****************************************************************
  // Lips data.
  // ****************************************************************

  sprintf(st, "<lips width=\"%2.4f\"/>", anatomy.lipsWidth_cm);
  os << string(indent, ' ') << string(st) << endl;

  // ****************************************************************
  // Tongue data.
  // ****************************************************************

  os << string(indent, ' ') << "<tongue>" << endl;
  indent+= 2;

  sprintf(st, "<tip radius=\"%2.4f\"/>", anatomy.tongueTipRadius_cm);
  os << string(indent, ' ') << string(st) << endl;

  sprintf(st, "<body radius_x=\"%2.4f\" radius_y=\"%2.4f\"/>", 
    anatomy.tongueCenterRadiusX_cm,
    anatomy.tongueCenterRadiusY_cm);
  os << string(indent, ' ') << string(st) << endl;

  sprintf(st, "<root automatic_calc=\"%d\" trx_slope=\"%2.4f\" trx_intercept=\"%2.4f\" try_slope=\"%2.4f\" try_intercept=\"%2.4f\"/>", 
    (int)anatomy.automaticTongueRootCalc,
    anatomy.tongueRootTrxSlope,
    anatomy.tongueRootTrxIntercept,
    anatomy.tongueRootTrySlope,
    anatomy.tongueRootTryIntercept);
  os << string(indent, ' ') << string(st) << endl;

  indent-= 2;
  os << string(indent, ' ') << "</tongue>" << endl;

  // ****************************************************************
  // Velum data.
  // ****************************************************************

  sprintf(st, "<velum uvula_width=\"%2.4f\" uvula_height=\"%2.4f\" uvula_depth=\"%2.4f\" "
    "max_nasal_port_area=\"%2.4f\" >",
    anatomy.uvulaWidth_cm, anatomy.uvulaHeight_cm, anatomy.uvulaDepth_cm,
    anatomy.maxNasalPortArea_cm2);

  os << string(indent, ' ') << string(st) << endl;
  indent+= 2;

  // Low velum
  os << string(indent, ' ') << "<low points=\"";
  for (i=0; i < 5; i++) { os << anatomy.velumLowPoints[i].x << " " << anatomy.velumLowPoints[i].y << " "; }
  os << "\"/>" << endl;

  // Mid velum
  os << string(indent, ' ') << "<mid points=\"";
  for (i=0; i < 5; i++) { os << anatomy.velumMidPoints[i].x << " " << anatomy.velumMidPoints[i].y << " "; }
  os << "\"/>" << endl;

  // High velum
  os << string(indent, ' ') << "<high points=\"";
  for (i=0; i < 5; i++) { os << anatomy.velumHighPoints[i].x << " " << anatomy.velumHighPoints[i].y << " "; }
  os << "\"/>" << endl;

  indent-= 2;
  os << string(indent, ' ') << "</velum>" << endl;

  // ****************************************************************
  // Pharynx.
  // ****************************************************************

  sprintf(st, "<pharynx fulcrum_x=\"%2.4f\" fulcrum_y=\"%2.4f\" rotation_angle_deg=\"%2.4f\" "
    "top_rib_y=\"%2.4f\" upper_depth=\"%2.4f\" lower_depth=\"%2.4f\" back_side_width=\"%2.4f\"/>",
    anatomy.pharynxFulcrum.x, anatomy.pharynxFulcrum.y, anatomy.pharynxRotationAngle_deg,
    anatomy.pharynxTopRibY_cm, anatomy.pharynxUpperDepth_cm, anatomy.pharynxLowerDepth_cm,
    anatomy.pharynxBackWidth_cm);

  os << string(indent, ' ') << string(st) << endl;

  // ****************************************************************
  // Larynx.
  // ****************************************************************

  sprintf(st, "<larynx upper_depth=\"%2.4f\" lower_depth=\"%2.4f\" epiglottis_width=\"%2.4f\" "
    "epiglottis_height=\"%2.4f\" epiglottis_depth=\"%2.4f\" epiglottis_angle_deg=\"%2.4f\">",
    anatomy.larynxUpperDepth_cm, anatomy.larynxLowerDepth_cm, anatomy.epiglottisWidth_cm,
    anatomy.epiglottisHeight_cm, anatomy.epiglottisDepth_cm, anatomy.epiglottisAngle_deg);

  os << string(indent, ' ') << string(st) << endl;
  indent+= 2;

  // Narrow larynx
  os << string(indent, ' ') << "<narrow points=\"";
  for (i=0; i < 8; i++) { os << anatomy.larynxNarrowPoints[i].x << " " << anatomy.larynxNarrowPoints[i].y << " "; }
  os << "\"/>" << endl;

  // Wide larynx
  os << string(indent, ' ') << "<wide points=\"";
  for (i=0; i < 8; i++) { os << anatomy.larynxWidePoints[i].x << " " << anatomy.larynxWidePoints[i].y << " "; }
  os << "\"/>" << endl;

  indent-= 2;
  os << string(indent, ' ') << "</larynx>" << endl;

  // ****************************************************************
  // Piriform fossae.
  // ****************************************************************

  sprintf(st, "<piriform_fossa length=\"%2.4f\" volume=\"%2.4f\"/>", 
    anatomy.piriformFossaLength_cm, anatomy.piriformFossaVolume_cm3);
  os << string(indent, ' ') << string(st) << endl;

  // ****************************************************************
  // Subglottal cavity.
  // ****************************************************************

  sprintf(st, "<subglottal_cavity length=\"%2.4f\"/>", anatomy.subglottalCavityLength_cm);
  os << string(indent, ' ') << string(st) << endl;

  // ****************************************************************
  // Nasal cavity.
  // ****************************************************************

  sprintf(st, "<nasal_cavity length=\"%2.4f\"/>", anatomy.nasalCavityLength_cm);
  os << string(indent, ' ') << string(st) << endl;

  // ****************************************************************
  // Write the parameter definitions.
  // ****************************************************************

  for (i=0; i < NUM_PARAMS; i++)
  {
    sprintf(st, "<param index=\"%d\"  name=\"%s\"  min=\"%2.3f\"  max=\"%2.3f\"  neutral=\"%2.3f\"  "
      "positive_velocity_factor=\"%2.3f\"  negative_velocity_factor=\"%2.3f\" />",
      i, param[i].abbr.c_str(), param[i].min, param[i].max, param[i].neutral,
      anatomy.positiveVelocityFactor[i], anatomy.negativeVelocityFactor[i]);
    os << string(indent, ' ') << st << endl;
  }

  // End of the <anatomy> element.
  
  indent-= 2;
  os << string(indent, ' ') << "</anatomy>" << endl;
}


// ****************************************************************************
/// Write the data for the vocal tract shapes in xml-format into the ouput 
/// stream os.
// ****************************************************************************

void VocalTract::writeShapesXml(std::ostream &os, int indent)
{
  int i, k;
  char st[1024];
  Shape *s = NULL;

  os << string(indent, ' ') << "<shapes>" << endl;
  indent+= 2;

  // ****************************************************************
  // Run through all vocal tract shapes.
  // ****************************************************************

  for (i=0; i < (int)shapes.size(); i++)
  {
    s = &(shapes[i]);
    os << string(indent, ' ') << "<shape name=\"" << s->name <<"\">" << endl;
    indent+= 2;

    for (k=0; k < NUM_PARAMS; k++)
    {
      sprintf(st, "<param name=\"%s\" value=\"%2.4f\"/>", 
        param[k].abbr.c_str(), s->param[k]);
      os << string(indent, ' ') << string(st) << endl;
    }

    indent-= 2;
    os << string(indent, ' ') << "</shape>" << endl;
  }

  indent-= 2;
  os << string(indent, ' ') << "</shapes>" << endl;
}


// ****************************************************************************
/// Writes the vocal tract anatomy data and the vocal tract shapes to the
/// given xml stream.
// ****************************************************************************

void VocalTract::writeToXml(std::ostream &os, int indent)
{
  os << string(indent, ' ') << "<vocal_tract_model>" << endl;
  indent+= 2;

  writeAnatomyXml(os, indent); 
  writeShapesXml(os, indent); 

  indent-= 2;
  os << string(indent, ' ') << "</vocal_tract_model>" << endl;
}


// ****************************************************************************
// ****************************************************************************

bool VocalTract::saveAsObjFile(const string &fileName, bool saveBothSides)
{
  int i, k;
  Point3D normal;

  const int NUM_EXPORT_SURFACES = 11;
  Surface *exportSurface[NUM_EXPORT_SURFACES];
  int firstVertex[NUM_EXPORT_SURFACES];
  int numVertices[NUM_EXPORT_SURFACES];
  Surface *s = NULL;

  if (saveBothSides)
  {
    exportSurface[0]  = &surface[TONGUE];
	  exportSurface[1]  = &surface[UPPER_TEETH_TWOSIDE];
	  exportSurface[2]  = &surface[LOWER_TEETH_TWOSIDE];
    exportSurface[3]  = &surface[UPPER_LIP_TWOSIDE];
    exportSurface[4]  = &surface[LOWER_LIP_TWOSIDE];
    exportSurface[5]  = &surface[UPPER_COVER_TWOSIDE];
    exportSurface[6]  = &surface[LOWER_COVER_TWOSIDE];
    exportSurface[7]  = &surface[LEFT_COVER];
    exportSurface[8]  = &surface[RIGHT_COVER];
    exportSurface[9]  = &surface[EPIGLOTTIS_TWOSIDE];
    exportSurface[10] = &surface[UVULA_TWOSIDE];
  }
  else
  {
    exportSurface[0]  = &surface[TONGUE];
    exportSurface[1]  = &surface[VocalTract::UPPER_TEETH];
	  exportSurface[2]  = &surface[VocalTract::LOWER_TEETH];
	  exportSurface[3]  = &surface[VocalTract::UPPER_LIP];
    exportSurface[4]  = &surface[VocalTract::LOWER_LIP];
    exportSurface[5]  = &surface[VocalTract::UPPER_COVER];
    exportSurface[6]  = &surface[VocalTract::LOWER_COVER];
    exportSurface[7]  = &surface[VocalTract::LEFT_COVER];
    exportSurface[8]  = NULL;
    exportSurface[9]  = &surface[VocalTract::EPIGLOTTIS];
    exportSurface[10] = &surface[VocalTract::UVULA];
  }

  string groupName[NUM_EXPORT_SURFACES] =
  {
    "TONGUE",
    "UPPER TEETH",
    "LOWER TEETH",
    "UPPER LIP",
    "LOWER LIP",
    "UPPER COVER",
    "LOWER COVER",
    "LEFT COVER",
    "RIGHT COVER",
    "EPIGLOTTIS",
    "UVULA"
  };

  string mtlName[NUM_EXPORT_SURFACES] =
  {
    "TONGUE_MTL",
    "TEETH_MTL",
    "TEETH_MTL",
    "LIP_MTL",
    "LIP_MTL",
    "COVER_MTL",
    "COVER_MTL",
    "COVER_MTL",
    "COVER_MTL",
    "COVER_MTL",
    "COVER_MTL"
  };

  // ****************************************************************
  // Create the material file name with and without the path.
  // ****************************************************************

  string longMtlFileName  = fileName;
  string shortMtlFileName = fileName;

  size_t numChars = longMtlFileName.length();
  if (numChars > 3)
  {
    // Replace obj by mtl.
    longMtlFileName[numChars-3] = 'm';
    longMtlFileName[numChars-2] = 't';
    longMtlFileName[numChars-1] = 'l';
  }
  
  size_t l = longMtlFileName.find_last_of('/');
  if (l == string::npos)
  {
    l = longMtlFileName.find_last_of('\\');
  }

  if (l != string::npos)
  {
    shortMtlFileName = longMtlFileName.substr(l+1);
  }

  // ****************************************************************
  // Create the material file with the material definitions.
  // ****************************************************************

  ofstream mtlFile(longMtlFileName.c_str());
  if (!mtlFile) 
  { 
    return false; 
  }

  double transTeeth = 0.5;
  double transLip   = 0.6;
  double transCover = 0.4;

  mtlFile << "newmtl TONGUE_MTL" << endl;
  mtlFile << "Ka 1.0 0.5 0.5" << endl;    // ambient reflectivity
  mtlFile << "Kd 1.0 0.5 0.5" << endl;    // diffuse reflectivity
  mtlFile << "Ks 1.0 1.0 1.0" << endl;    // specular reflectivity
  mtlFile << "Ns 50" << endl;             // specular exponent
  mtlFile << "d 1.0" << endl;             // opaque
  
  mtlFile << "newmtl TEETH_MTL" << endl;
  mtlFile << "Ka 1.0 1.0 1.0" << endl;    // ambient reflectivity
  mtlFile << "Kd 0.9 0.9 0.9" << endl;    // diffuse reflectivity
  mtlFile << "Ks 1.0 1.0 1.0" << endl;    // specular reflectivity
  mtlFile << "Ns 50" << endl;             // specular exponent
  mtlFile << "Tf " << transTeeth << " " << transTeeth << " " << transTeeth << endl;

  mtlFile << "newmtl LIP_MTL" << endl;
  mtlFile << "Ka 1.0 0.6 0.6" << endl;    // ambient reflectivity
  mtlFile << "Kd 1.0 0.6 0.6" << endl;    // diffuse reflectivity
  mtlFile << "Ks 1.0 1.0 1.0" << endl;    // specular reflectivity
  mtlFile << "Ns 50" << endl;             // specular exponent
  mtlFile << "Tf " << transLip << " " << transLip << " " << transLip << endl;
  
  mtlFile << "newmtl COVER_MTL" << endl;
  mtlFile << "Ka 0.5 0.5 0.5" << endl;    // ambient reflectivity
  mtlFile << "Kd 0.8 0.8 0.8" << endl;    // diffuse reflectivity
  mtlFile << "Ks 1.0 1.0 1.0" << endl;    // specular reflectivity
  mtlFile << "Ns 50" << endl;             // specular exponent
  mtlFile << "Tf " << transCover << " " << transCover << " " << transCover << endl;

  mtlFile.close();


  // ****************************************************************
  // Calc. the first vertex and the number of vertices of all
  // surfaces, when they are lined up in a single sequence.
  // ****************************************************************

  firstVertex[0] = 0;
  numVertices[0] = exportSurface[0]->numVertices;

  for (i=1; i < NUM_EXPORT_SURFACES; i++)
  {
    firstVertex[i] = firstVertex[i-1] + numVertices[i-1];

    if (exportSurface[i] != NULL)
    {
      numVertices[i] = exportSurface[i]->numVertices;
    }
    else
    {
      numVertices[i] = 0;
    }
  }

  // ****************************************************************
  // Calculate the normal vectors on all surfaces to be saved.
  // ****************************************************************

  for (i=0; i < NUM_EXPORT_SURFACES; i++)
  {
    if (exportSurface[i] != NULL)
    {
      exportSurface[i]->calculateNormals();
    }
  }

  // ****************************************************************
  // Some plane normals must be adjusted in order to avoid sharp 
  // edges.
  // ****************************************************************

  // Normals at the interface between the upper and lower cover

  Point3D P;
  int numRibs = NUM_LARYNX_RIBS + NUM_PHARYNX_RIBS;
  Surface *upperCover = exportSurface[5];
  Surface *lowerCover = exportSurface[6];

  for (i=0; i < numRibs; i++)
  {
    P = upperCover->getNormal(i, 0) + lowerCover->getNormal(i, 0);
    P.normalize();
    upperCover->setNormal(i, 0, P);
    lowerCover->setNormal(i, 0, P);

    if (saveBothSides) 
    { 
      P = upperCover->getNormal(i, upperCover->numRibPoints-1) + 
          lowerCover->getNormal(i, lowerCover->numRibPoints-1);
      P.normalize();
      upperCover->setNormal(i, upperCover->numRibPoints-1, P); 
      lowerCover->setNormal(i, lowerCover->numRibPoints-1, P);
    }
  }

  // Normals at the edge of the filling surfaces.

  Surface *leftCover = exportSurface[7];
  Surface *rightCover = exportSurface[8];

  int firstRib = NUM_LARYNX_RIBS + NUM_PHARYNX_RIBS - 1;
  int lastRib  = NUM_LARYNX_RIBS + NUM_PHARYNX_RIBS + NUM_VELUM_RIBS - 1;

  if (leftCover != NULL)
  {
    P = upperCover->getNormal(firstRib, 0);
    for (i=firstRib; i <= lastRib; i++) { upperCover->setNormal(i, 0, P); }
    leftCover->setNormal(0, 0, P);
    leftCover->setNormal(0, 1, P);
    leftCover->setNormal(0, 2, P);
    leftCover->setNormal(0, 3, P);
    
    leftCover->setNormal(1, 0, P);
    leftCover->setNormal(1, 1, P);

    P = lowerCover->getNormal(NUM_LARYNX_RIBS + NUM_THROAT_RIBS, 0);
    leftCover->setNormal(1, 2, P);
    leftCover->setNormal(1, 3, P);
  }

  if (rightCover != NULL)
  {
    P = upperCover->getNormal(firstRib, upperCover->numRibPoints-1);
    for (i=firstRib; i <= lastRib; i++) { upperCover->setNormal(i, upperCover->numRibPoints-1, P); }
    rightCover->setNormal(0, 0, P);
    rightCover->setNormal(0, 1, P);
    rightCover->setNormal(0, 2, P);
    rightCover->setNormal(0, 3, P);
    
    rightCover->setNormal(1, 0, P);
    rightCover->setNormal(1, 1, P);

    P = lowerCover->getNormal(NUM_LARYNX_RIBS + NUM_THROAT_RIBS, lowerCover->numRibPoints-1);
    rightCover->setNormal(1, 2, P);
    rightCover->setNormal(1, 3, P);
  }
  
  // Beginning and end of the uvula ribs.

  if (saveBothSides)
  {
    Surface *s = &surface[UVULA_TWOSIDE];
    if (s != NULL)
    {
      for (i=0; i < s->numRibs; i++)
      {
        P = s->getNormal(i, 0) + s->getNormal(i, s->numRibPoints-1);
        P.normalize();
        s->setNormal(i, 0, P);
        s->setNormal(i, s->numRibPoints-1, P);
      }
    }
  }

  // Beginning and end of the epiglottis ribs.

  if (saveBothSides)
  {
    Surface *s = &surface[EPIGLOTTIS_TWOSIDE];
    if (s != NULL)
    {
      for (i=0; i < s->numRibs; i++)
      {
        P = s->getNormal(i, 0) + s->getNormal(i, s->numRibPoints-1);
        P.normalize();
        s->setNormal(i, 0, P);
        s->setNormal(i, s->numRibPoints-1, P);
      }
    }
  }

  // ****************************************************************
  // Open the file and output the name of the material file.
  // ****************************************************************

  ofstream file(fileName);
  if (!file) 
  { 
    return false; 
  }

  file << "mtllib " << shortMtlFileName << endl;

  // ****************************************************************
  // Output the vertices.
  // ****************************************************************

  for (k=0; k < NUM_EXPORT_SURFACES; k++)
  {
    s = exportSurface[k];
    if (s != NULL)
    {
      for (i=0; i < numVertices[k]; i++)
      {
        file << "v  "
          << (double)s->vertex[i].coord.x << "  " 
          << (double)s->vertex[i].coord.y << "  " 
          << (double)s->vertex[i].coord.z << endl;
      }
    }
  }

  // ****************************************************************
  // Output the vertex normals.
  // ****************************************************************

  for (k=0; k < NUM_EXPORT_SURFACES; k++)
  {
    s = exportSurface[k];
    if (s != NULL)
    {
      for (i=0; i < numVertices[k]; i++)
      {
        normal = s->getNormal(s->vertex[i].rib, s->vertex[i].ribPoint);
        file << "vn  "
          << (double)normal.x << "  " 
          << (double)normal.y << "  " 
          << (double)normal.z << endl;
      }
    }
  }

  // ****************************************************************
  // Output the faces (triangles).
  // ****************************************************************

  int index0, index1, index2;   // Indices of the three corners

  for (k=0; k < NUM_EXPORT_SURFACES; k++)
  {
    s = exportSurface[k];
    if (s != NULL)
    {
      file << "g " << groupName[k] << endl;
      file << "usemtl " << mtlName[k] << endl;

      for (i=0; i < s->numTriangles; i++)
      {
        index0 = 1 + firstVertex[k] + s->triangle[i].vertex[0];
        index1 = 1 + firstVertex[k] + s->triangle[i].vertex[1];
        index2 = 1 + firstVertex[k] + s->triangle[i].vertex[2];

        file << "f  "
          << index0 << "//" << index0 << "  " 
          << index1 << "//" << index1 << "  " 
          << index2 << "//" << index2 << endl;
      }
    }
  }

  file << endl;

  // The file is closed automatically.
  return true;
}


// ****************************************************************************
/// Convenience function to set all control parameters at once.
// ****************************************************************************

void VocalTract::setParams(double *controlParams)
{
  int i;
  for (i = 0; i < NUM_PARAMS; i++)
  {
    param[i].x = controlParams[i];
  }
}


// ****************************************************************************
/// Calculates all surfaces, the center line, and the area functions.
// ****************************************************************************

void VocalTract::calculateAll()
{
  // ****************************************************************
  // Set the limited parameter values equal to the set values.
  // ****************************************************************

  int i;
  for (i=0; i < NUM_PARAMS; i++)
  {
    restrictParam(i);
    param[i].limitedX = param[i].x;
  }

  // ****************************************************************
  // Do the calculations.
  // ****************************************************************

  calcSurfaces();
  calcCenterLine();
  calcCrossSections();
  crossSectionsToTubeSections();
}


// ****************************************************************************
/// Calculate all surfaces of the model.
// ****************************************************************************

void VocalTract::calcSurfaces()
{
  int i, k;
  int rib;
  Point3D P, Q, R, A, v;
  Point3D vertex;
  double s, t;
  double dx, dy;
  double x, y, z;
  double cosinus;
  double sinus;
  double angle_rad;
  double shearCoeff = 
    cos(anatomy.pharynxRotationAngle_deg*M_PI/180.0) / 
    sin(anatomy.pharynxRotationAngle_deg*M_PI/180.0);
  double shearX;

  // ****************************************************************
  // Because the surfaces change their shape here, all intersections
  // need to be newly prepared before the intersectioning.
  // ****************************************************************

  for (i=0; i < NUM_SURFACES; i++)
  {
    intersectionsPrepared[i] = false;
  }

  // ****************************************************************
  // Restrict the parameter values.
  // ****************************************************************

  for (i=0; i < NUM_PARAMS; i++)
  {
    restrictParam(i);
  }

  // ****************************************************************
  // Additional restrictions for parameter values:
  // o HY may not be higher than the lowest point of the most posterior
  //   jaw rib.
  // o The most anterior larynx point (hyoid) may not be right of
  //   the most posterior jaw rib -> constrain HX!
  // ****************************************************************

  angle_rad = param[JA].x*M_PI / 180.0;
  cosinus = cos(angle_rad);
  sinus = sin(angle_rad);

  // Vertex of the lowest point of the most posterior jaw rib.
  vertex = surface[MANDIBLE].getVertex(0, NUM_LOWER_COVER_POINTS - 1);

  // (dx, dy) is the vertex position relativ to the fulcrum before the rotation
  dx = vertex.x + anatomy.jawRestPos.x + param[JX].x - anatomy.jawFulcrum.x;
  dy = vertex.y + anatomy.jawRestPos.y - anatomy.jawFulcrum.y;
  vertex.x = cosinus*dx - sinus*dy + anatomy.jawFulcrum.x;
  vertex.y = sinus*dx + cosinus*dy + anatomy.jawFulcrum.y;
  vertex.z = vertex.z;

  if (param[HY].x > vertex.y)
  {
    param[HY].limitedX = vertex.y;
  }

  // ****************************************************************

  P = surface[NARROW_LARYNX_FRONT].getVertex(NUM_LARYNX_RIBS - 1, NUM_LOWER_COVER_POINTS - 1);
  Q = surface[WIDE_LARYNX_FRONT].getVertex(NUM_LARYNX_RIBS - 1, NUM_LOWER_COVER_POINTS - 1);

  // Horizontal offset for the larynx surfaces
  double xOffset = getPharynxBackX(param[HY].limitedX);

  t = Q.x - P.x;
  if (fabs(t) < 0.000001)
  {
    t = 0.000001;
  }
  double maxHX = (vertex.x - xOffset - P.x) / t;

  if (param[HX].x > maxHX)
  {
    param[HX].limitedX = maxHX;
  }

  // Also consider the normal allowed range.
  if (param[HX].limitedX < param[HX].min)
  {
    param[HX].limitedX = param[HX].min;
  }
  if (param[HX].limitedX > param[HX].max)
  {
    param[HX].limitedX = param[HX].max;
  }

  // ****************************************************************
  // The upper cover is the combination of the larynx, the pharynx
  // back side, the velum and the palate.
  // ****************************************************************
  
  rib = 0;

  // ****************************************************************
  // Larynx back side (use shearing to adapt to an oblique pharynx 
  // back wall).
  // ****************************************************************

  t = param[HX].limitedX;
  x = getPharynxBackX(param[HY].limitedX);   // Horizontal offset for the larynx surfaces

  for (i=0; i < NUM_LARYNX_RIBS; i++)
  {
    for (k=0; k < NUM_UPPER_COVER_POINTS; k++)
    {
      P = surface[NARROW_LARYNX_BACK].getVertex(i, k);
      Q = surface[WIDE_LARYNX_BACK].getVertex(i, k);
      P = P*(1.0-t) + Q*t;
      shearX = P.y*shearCoeff;
      P.x+= x + shearX;
      P.y+= param[HY].limitedX;

      surface[UPPER_COVER].setVertex(rib, k, P);
    }
    rib++;
  }

  // ****************************************************************
  // Pharynx back side.
  // ****************************************************************

  const double pharynxX_cm[NUM_UPPER_COVER_POINTS] = {  1.0, 0.69, 0.41, 0.19, 0.05, 0.0 };
  const double pharynxZ_cm[NUM_UPPER_COVER_POINTS] = { -1.0, -0.95, -0.81, -0.59, -0.31, 0.0 };
  double scaleX, scaleZ;

  for (i=0; i < NUM_PHARYNX_RIBS; i++)
  {
    t = (double)i/(double)(NUM_PHARYNX_RIBS-1);
    scaleZ = 0.5*((1.0-t)*anatomy.pharynxLowerDepth_cm + t*anatomy.pharynxUpperDepth_cm);
    scaleX = anatomy.pharynxBackWidth_cm;

    t = (double)(i+1)/(double)NUM_PHARYNX_RIBS;
    y = (1.0-t)*param[HY].limitedX + t*anatomy.pharynxTopRibY_cm;
    x = getPharynxBackX(y);

    for (k=0; k < NUM_UPPER_COVER_POINTS; k++)
    {
      P.x = x + scaleX*pharynxX_cm[k];
      P.y = y;
      P.z = scaleZ*pharynxZ_cm[k];
      surface[UPPER_COVER].setVertex(rib, k, P);
    }
    rib++;
  }

  // ****************************************************************
  // Soft palate.
  // ****************************************************************

  // Relative velum shape parameter value.
  s = (param[VS].x - param[VS].min) / (param[VS].max - param[VS].min);

  // Relative velic opening parameter value.
  t = param[VO].x;
  // Restrict t to [0, 1].
  if (t > 1.0) { t = 1.0; }   
  if (t < 0.0) { t = 0.0; }

  for (i=0; i < NUM_VELUM_RIBS; i++)
  {
    for (k=0; k < NUM_UPPER_COVER_POINTS; k++)
    {
      P = surface[HIGH_VELUM].getVertex(i, k);
      Q = surface[MID_VELUM].getVertex(i, k);
      R = surface[LOW_VELUM].getVertex(i, k);

      // The point between the high and mid closed shapes.
      A = (1.0-s)*P + s*Q;
      // The point between the interpolated closed and the open shape.
      A = (1.0-t)*A + t*R;

      surface[UPPER_COVER].setVertex(rib, k, A);
    }
    rib++;
  }


  // ****************************************************************
  // Hard palate.
  // ****************************************************************

  double factor;
  int lastVelumRib = rib-1;

  for (i=0; i < NUM_JAW_RIBS; i++)
  {
    // Interpolate the y-values of the posterior hard palate ribs
    // between their original values and the most anterior velum rib
    factor = anatomy.palatePoints[i].x / 1.2;   // Max. infuence for 1.2 cm
    if (factor < 1.0)
    {
      for (k=0; k < NUM_UPPER_COVER_POINTS; k++)
      {
        P = surface[UPPER_COVER].getVertex(lastVelumRib, k);
        Q = surface[PALATE].getVertex(i, k);
        vertex.x = Q.x;
        vertex.y = (1.0-factor)*P.y + factor*Q.y;
        if (Q.x < 0.01) { vertex.z = P.z; } else { vertex.z = Q.z; }
        surface[UPPER_COVER].setVertex(rib, k, vertex);
      }
    }
    else
    {
      for (k=0; k < NUM_UPPER_COVER_POINTS; k++)
      {
        vertex = surface[PALATE].getVertex(i, k);
        surface[UPPER_COVER].setVertex(rib, k, vertex);
      }
    }

    rib++;
  }

  // *****************************************************************
  // The lower cover is the combination of the larynx, the pharynx
  // front side and the lower jaw.
  // *****************************************************************

  rib = 0;

  // ****************************************************************
  // Larynx front side.
  // ****************************************************************

  t = param[HX].limitedX;
  x = getPharynxBackX(param[HY].limitedX);   // Horizontal offset for the larynx surfaces

  for (i=0; i < NUM_LARYNX_RIBS; i++)
  {
    for (k=0; k < NUM_LOWER_COVER_POINTS; k++)
    {
      P = surface[NARROW_LARYNX_FRONT].getVertex(i, k);
      Q = surface[WIDE_LARYNX_FRONT].getVertex(i, k);
      P = P*(1.0-t) + Q*t;
      shearX = P.y*shearCoeff;
      P.x+= x + shearX;
      P.y+= param[HY].limitedX;

      surface[LOWER_COVER].setVertex(rib, k, P);
    }
    rib++;
  }

  rib+= NUM_THROAT_RIBS;

  // ****************************************************************
  // Mandible.
  // ****************************************************************

  angle_rad = param[JA].x*M_PI / 180.0;
  cosinus = cos(angle_rad);
  sinus = sin(angle_rad);

  for (i=0; i < NUM_JAW_RIBS; i++)
  {
    for (k=0; k < NUM_LOWER_COVER_POINTS; k++)
    {
      vertex = surface[MANDIBLE].getVertex(i, k);

      // (dx, dy) is the vertex position relativ to the fulcrum before the rotation
      dx = vertex.x + anatomy.jawRestPos.x + param[JX].x - anatomy.jawFulcrum.x;
      dy = vertex.y + anatomy.jawRestPos.y - anatomy.jawFulcrum.y;
      x = cosinus*dx - sinus*dy + anatomy.jawFulcrum.x;
      y = sinus*dx + cosinus*dy + anatomy.jawFulcrum.y;
      z = vertex.z;
      
      surface[LOWER_COVER].setVertex(rib, k, x, y, z);
    }
    rib++;

    // Rotate the approximation of the lower gum edges

    vertex = lowerGumsOuterEdgeOrig[i];

    // (dx, dy) is the vertex position relativ to the fulcrum before the rotation
    dx = vertex.x + anatomy.jawRestPos.x + param[JX].x - anatomy.jawFulcrum.x;
    dy = vertex.y + anatomy.jawRestPos.y - anatomy.jawFulcrum.y;
    x = cosinus*dx - sinus*dy + anatomy.jawFulcrum.x;
    y = sinus*dx + cosinus*dy + anatomy.jawFulcrum.y;
    z = vertex.z;

    lowerGumsOuterEdge[i].set(x, y, z);

    vertex = lowerGumsInnerEdgeOrig[i];

    // (dx, dy) is the vertex position relativ to the fulcrum before the rotation
    dx = vertex.x + anatomy.jawRestPos.x + param[JX].x - anatomy.jawFulcrum.x;
    dy = vertex.y + anatomy.jawRestPos.y - anatomy.jawFulcrum.y;
    x = cosinus*dx - sinus*dy + anatomy.jawFulcrum.x;
    y = sinus*dx + cosinus*dy + anatomy.jawFulcrum.y;
    z = vertex.z;

    lowerGumsInnerEdge[i].set(x, y, z);
  }

  // ****************************************************************
  // Front side of the pharynx (throat).
  // ****************************************************************

  Point3D C[3];   // Points on the rational Bezier curve
  double weight[3] = { 1.0, 0.71, 1.0 };
  BezierCurve3D curve;

  R = surface[LOWER_COVER].getVertex(NUM_LARYNX_RIBS-1, NUM_LOWER_COVER_POINTS-1);
  C[2] = R;

  for (i=0; i < NUM_THROAT_RIBS; i++)
  {
    Q = surface[UPPER_COVER].getVertex(NUM_LARYNX_RIBS+i, 0);
    C[0] = Q;
    C[1].set(R.x, R.y, Q.z);

    curve.setPoints(3, C, weight);

    for (k=0; k < NUM_LOWER_COVER_POINTS; k++)
    {
      Q = curve.getPoint((double)k/(double)(NUM_LOWER_COVER_POINTS-1));
      surface[LOWER_COVER].setVertex(NUM_LARYNX_RIBS+i, k, Q);
    }
  }

  // ****************************************************************
  // Lower teeth.
  // ****************************************************************

  for (i=0; i < NUM_TEETH_RIBS; i++)
  {
    for (k=0; k < NUM_TEETH_POINTS; k++)
    {
      vertex = surface[LOWER_TEETH_ORIGINAL].getVertex(i, k);
      
      // (dx, dy) is the vertex position relativ to the fulcrum before the rotation
      dx = vertex.x + anatomy.jawRestPos.x + param[JX].x - anatomy.jawFulcrum.x;
      dy = vertex.y + anatomy.jawRestPos.y - anatomy.jawFulcrum.y;
      x = cosinus*dx - sinus*dy + anatomy.jawFulcrum.x;
      y = sinus*dx + cosinus*dy + anatomy.jawFulcrum.y;
      z = vertex.z;

      surface[LOWER_TEETH].setVertex(i, k, x, y, z);
    }
  }

  // ****************************************************************
  // The lips.
  // Make sure that upperGumsOuterEdge[] and lowerGumsOuterEdge[]
  // were already calculated !!
  // ****************************************************************

  calcLips();

  // ****************************************************************
  // The gap between the upper and lower teeth
  // ****************************************************************
  
  // 1st rib
  P = surface[UPPER_COVER].getVertex(NUM_LARYNX_RIBS+NUM_PHARYNX_RIBS, 0);
  surface[LEFT_COVER].setVertex(0, 0, P);
  surface[LEFT_COVER].setVertex(0, 1, P);
  P = surface[LOWER_COVER].getVertex(NUM_LARYNX_RIBS+NUM_PHARYNX_RIBS-1, 0);
  surface[LEFT_COVER].setVertex(0, 2, P);
  surface[LEFT_COVER].setVertex(0, 3, P);

  // 2nd rib
  P = surface[UPPER_COVER].getVertex(NUM_LARYNX_RIBS+NUM_PHARYNX_RIBS+NUM_VELUM_RIBS-1, 0);
  surface[LEFT_COVER].setVertex(1, 0, P);
  surface[LEFT_COVER].setVertex(1, 1, P);
  P = surface[LOWER_COVER].getVertex(NUM_LARYNX_RIBS+NUM_PHARYNX_RIBS, 0);
  surface[LEFT_COVER].setVertex(1, 2, P);
  surface[LEFT_COVER].setVertex(1, 3, P);

  // 3rd rib
  P = surface[UPPER_COVER].getVertex(NUM_LARYNX_RIBS+NUM_PHARYNX_RIBS+NUM_VELUM_RIBS, 0);
  surface[LEFT_COVER].setVertex(2, 0, P);
  surface[LEFT_COVER].setVertex(2, 1, P);
  P = surface[LOWER_COVER].getVertex(NUM_LARYNX_RIBS+NUM_PHARYNX_RIBS, 0);
  surface[LEFT_COVER].setVertex(2, 2, P);
  surface[LEFT_COVER].setVertex(2, 3, P);

  // Last teeth rib of the upper and lower edge of the fill surface

  double length;
  double zLip[2];
  Point3D innerPoint, outerPoint;
  Point3D *outerEdge;
  Point3D *innerEdge;

  zLip[0] = surface[UPPER_LIP].getVertex(0, 0).z;
  zLip[1] = surface[LOWER_LIP].getVertex(0, 0).z;

  // k=0 -> upper gums; k=1 -> lower gums ***************************

  for (k=0; k < 2; k++)
  {
    if (k == 0) 
    { 
      innerEdge = upperGumsInnerEdge; 
      outerEdge = upperGumsOuterEdge; 
    } 
    else 
    { 
      innerEdge = lowerGumsInnerEdge; 
      outerEdge = lowerGumsOuterEdge; 
    }

    for (i=0; i < NUM_JAW_RIBS-1; i++)
    {
      P = innerEdge[i]; 
      Q = innerEdge[i+1];

      length = Q.z - P.z;

      if (P.z <= zLip[k]) 
      { 
        innerPoint = innerEdge[i];
        outerPoint = outerEdge[i]; 
      }
      if ((P.z <= zLip[k]) && (Q.z >= zLip[k]) && (length > 0.0))
      {
        t = (zLip[k] - P.z) / length;
        innerPoint = innerEdge[i] + t*(innerEdge[i+1]-innerEdge[i]);
        outerPoint = outerEdge[i] + t*(outerEdge[i+1]-outerEdge[i]);
      }

      if (k == 0)
      {
        surface[LEFT_COVER].setVertex(3+i, 0, innerPoint);
        surface[LEFT_COVER].setVertex(3+i, 1, outerPoint);
      }
      else
      {
        surface[LEFT_COVER].setVertex(3+i, 2, outerPoint);
        surface[LEFT_COVER].setVertex(3+i, 3, innerPoint);
      }
    }

    // Last rib of the fill surface
    if (k == 0)
    {
      surface[LEFT_COVER].setVertex(NUM_FILL_RIBS-1, 0, innerPoint);
      surface[LEFT_COVER].setVertex(NUM_FILL_RIBS-1, 1, outerPoint);
    }
    else
    {
      surface[LEFT_COVER].setVertex(NUM_FILL_RIBS-1, 2, outerPoint);
      surface[LEFT_COVER].setVertex(NUM_FILL_RIBS-1, 3, innerPoint);
    }
  }

  // ****************************************************************
  // Calculation of the contours of the upper and lower cover in the
  // midsagittal plane.
  // ****************************************************************

  int upperCoverPoint = NUM_UPPER_COVER_POINTS-1;
  int lowerCoverPoint = NUM_LOWER_COVER_POINTS-1;
  int tonguePoint     = NUM_TONGUE_POINTS/4;    // nicht genau den Punkt in der Mitte nehmen
  int lipRib          = NUM_LIP_RIBS-1;
  int teethRib        = NUM_TEETH_RIBS-1;

  // ****************************************************************
  // The upper contour.
  // ****************************************************************

  upperOutline.reset(0);
  for (i=0; i < NUM_UPPER_COVER_RIBS; i++)
  {
    upperOutline.addPoint(surface[UPPER_COVER].getVertex(i, upperCoverPoint).toPoint2D());
  }

  for (i=0; i < NUM_TEETH_POINTS-1; i++)
  {
    upperOutline.addPoint(surface[UPPER_TEETH].getVertex(teethRib, i).toPoint2D());
  }

  for (i=0; i < NUM_INNER_LIP_POINTS+2; i++)
  {
    upperOutline.addPoint(surface[UPPER_LIP].getVertex(lipRib, i).toPoint2D());
  }

  // Delete the last (wrapped) points.
  while (upperOutline.getControlPoint(upperOutline.getNumPoints()-1).x <= 
         upperOutline.getControlPoint(upperOutline.getNumPoints()-2).x) { upperOutline.delPoint(); }

  // ****************************************************************
  // The lower contour.
  // ****************************************************************

  lowerOutline.reset(0);

  for (i=0; i < NUM_LOWER_COVER_RIBS; i++)
  {
    lowerOutline.addPoint(surface[LOWER_COVER].getVertex(i, lowerCoverPoint).toPoint2D());
  }

  for (i=0; i < NUM_TEETH_POINTS-1; i++)
  {
    lowerOutline.addPoint(surface[LOWER_TEETH].getVertex(teethRib, i).toPoint2D());
  }

  for (i=0; i < NUM_INNER_LIP_POINTS+2; i++)
  {
    lowerOutline.addPoint(surface[LOWER_LIP].getVertex(lipRib, i).toPoint2D());
  }

  while (lowerOutline.getControlPoint(lowerOutline.getNumPoints()-1).x <= 
         lowerOutline.getControlPoint(lowerOutline.getNumPoints()-2).x) { lowerOutline.delPoint(); }

  // ****************************************************************
  // Calculation of the tongue geometry.
  // ****************************************************************

  calcTongue();

  // ****************************************************************
  // Epiglottis. Must be calculated AFTER the tongue.
  // ****************************************************************

  P = surface[LOWER_COVER].getVertex(NUM_LARYNX_RIBS-2, NUM_LOWER_COVER_POINTS-1);
  Q = surface[LOWER_COVER].getVertex(NUM_LARYNX_RIBS-1, NUM_LOWER_COVER_POINTS-1);
  v = Q - P;
  v.normalize();

  A = P + anatomy.epiglottisWidth_cm*v;

  // Align the first epiglottal rib with v
  for (k=0; k < NUM_EPIGLOTTIS_POINTS; k++)
  {
    P = surface[EPIGLOTTIS_ORIGINAL].getVertex(0, k);
    x = v.x*P.x - v.y*P.y;
    y = v.y*P.x + v.x*P.y;
    z = P.z;
    surface[EPIGLOTTIS].setVertex(0, k, x + A.x, y + A.y, z + A.z);
  }

  // The 2nd, 3rd and 4th rib
  double minAngle_rad = anatomy.epiglottisAngle_deg*M_PI/180.0;

  for (i=1; i < 10; i++)        // Check the first few tongue ribs
  {
    P = surface[TONGUE].getVertex(i, NUM_TONGUE_POINTS/2) - A;
    if (P.y < anatomy.epiglottisHeight_cm)
    {
      if (P.y < 0.01) { P.y = 0.01; }
      angle_rad = atan2(P.y, P.x);
      if (angle_rad > minAngle_rad) { minAngle_rad = angle_rad; }
    }
  }

  minAngle_rad-= 0.5*M_PI;       // Subtract 90 deg
  if (minAngle_rad > 0.25*M_PI) { minAngle_rad = 0.25*M_PI; }
  sinus = sin(minAngle_rad);
  cosinus = cos(minAngle_rad);

  for (i=1; i < NUM_EPIGLOTTIS_RIBS; i++)
  {
    for (k=0; k < NUM_EPIGLOTTIS_POINTS; k++)
    {
      P = surface[EPIGLOTTIS_ORIGINAL].getVertex(i, k);
      x = cosinus*P.x - sinus*P.y;
      y = sinus*P.x + cosinus*P.y;
      z = P.z;
      surface[EPIGLOTTIS].setVertex(i, k, x + A.x, y + A.y, z + A.z);
    }
  }

  // ****************************************************************
  // Uvula. Must be calculated AFTER the tongue.
  // ****************************************************************

  A = surface[UPPER_COVER].getVertex(NUM_LARYNX_RIBS+NUM_PHARYNX_RIBS+1, NUM_UPPER_COVER_POINTS-1);

  for (i=0; i < NUM_UVULA_RIBS; i++)
  {
    for (k=0; k < NUM_UVULA_POINTS; k++)
    {
      P = surface[UVULA_ORIGINAL].getVertex(i, k);
      P+= A;
      surface[UVULA].setVertex(i, k, P);
    }
  }

  // ****************************************************************
  // The tongue contour.
  // ****************************************************************

  tongueOutline.reset(0);
  for (i=0; i < NUM_TONGUE_RIBS; i++)
  {
    tongueOutline.addPoint(surface[TONGUE].getVertex(i, tonguePoint).toPoint2D());
  }

  // ****************************************************************
  // The epiglottis contour (for center line calculations).
  // ****************************************************************

  epiglottisOutline.reset(0);
  for (i=0; i < NUM_EPIGLOTTIS_RIBS; i++)
  {
    epiglottisOutline.addPoint(surface[EPIGLOTTIS].getVertex(i, NUM_EPIGLOTTIS_POINTS-1).toPoint2D());
  }
  for (i=NUM_EPIGLOTTIS_RIBS-2; i >= 0; i--)
  {
    epiglottisOutline.addPoint(surface[EPIGLOTTIS].getVertex(i, 0).toPoint2D());
  }


  // ****************************************************************
  // Two-sided grids for the upper and lower cover.        
  // ****************************************************************

  Surface *source, *target;

  source = &surface[UPPER_COVER];
  target = &surface[UPPER_COVER_TWOSIDE];

  for (i=0; i < target->numRibs; i++)
  {
    for (k=0; k < source->numRibPoints; k++)
    {
      P = source->getVertex(i, k);
      target->setVertex(i, k, P);
      P.z = -P.z;
      target->setVertex(i, target->numRibPoints-1-k, P);
    }
  }

  source = &surface[LOWER_COVER];
  target = &surface[LOWER_COVER_TWOSIDE];

  for (i=0; i < target->numRibs; i++)
  {
    for (k=0; k < source->numRibPoints; k++)
    {
      P = source->getVertex(i, k);
      target->setVertex(i, k, P);
      P.z = -P.z;
      target->setVertex(i, target->numRibPoints-1-k, P);
    }
  }

  // ****************************************************************
  // Two-sided grids for the upper and lower lips and teeth.
  // ****************************************************************

  source = &surface[UPPER_TEETH];
  target = &surface[UPPER_TEETH_TWOSIDE];

  for (i=0; i < source->numRibs; i++)
  {
    for (k=0; k < source->numRibPoints; k++)
    {
      P = source->getVertex(i, k);
      target->setVertex(i, k, P);
      P.z = -P.z;
      target->setVertex(target->numRibs-1-i, k, P);
    }
  }

  source = &surface[LOWER_TEETH];
  target = &surface[LOWER_TEETH_TWOSIDE];

  for (i=0; i < source->numRibs; i++)
  {
    for (k=0; k < source->numRibPoints; k++)
    {
      P = source->getVertex(i, k);
      target->setVertex(i, k, P);
      P.z = -P.z;
      target->setVertex(target->numRibs-1-i, k, P);
    }
  }

  source = &surface[UPPER_LIP];
  target = &surface[UPPER_LIP_TWOSIDE];

  for (i=0; i < source->numRibs; i++)
  {
    for (k=0; k < source->numRibPoints; k++)
    {
      P = source->getVertex(i, k);
      target->setVertex(i, k, P);
      P.z = -P.z;
      target->setVertex(target->numRibs-1-i, k, P);
    }
  }

  source = &surface[LOWER_LIP];
  target = &surface[LOWER_LIP_TWOSIDE];

  for (i=0; i < source->numRibs; i++)
  {
    for (k=0; k < source->numRibPoints; k++)
    {
      P = source->getVertex(i, k);
      target->setVertex(i, k, P);
      P.z = -P.z;
      target->setVertex(target->numRibs-1-i, k, P);
    }
  }

  // ****************************************************************
  // The right filling surface between the upper and lower teeth.
  // ****************************************************************

  source = &surface[LEFT_COVER];
  target = &surface[RIGHT_COVER];

  for (i=0; i < source->numRibs; i++)
  {
    for (k=0; k < source->numRibPoints; k++)
    {
      P = source->getVertex(i, k);
      P.z = -P.z;
      target->setVertex(i, k, P);
    }
  }

  // ****************************************************************
  // The epiglottis and uvula.
  // ****************************************************************

  source = &surface[EPIGLOTTIS];
  target = &surface[EPIGLOTTIS_TWOSIDE];

  for (i=0; i < target->numRibs; i++)
  {
    for (k=0; k < source->numRibPoints; k++)
    {
      P = source->getVertex(i, k);
      target->setVertex(i, k, P);
      P.z = -P.z;
      target->setVertex(i, target->numRibPoints-1-k, P);
    }
  }

  source = &surface[UVULA];
  target = &surface[UVULA_TWOSIDE];

  for (i=0; i < target->numRibs; i++)
  {
    for (k=0; k < source->numRibPoints; k++)
    {
      P = source->getVertex(i, k);
      target->setVertex(i, k, P);
      P.z = -P.z;
      target->setVertex(i, target->numRibPoints-1-k, P);
    }
  }

}


// ****************************************************************************
// Calculation of the lip shape.
// ****************************************************************************

void VocalTract::calcLips()
{
  Surface *upperLip   = &surface[UPPER_LIP];
  Surface *lowerLip   = &surface[LOWER_LIP];

  const Point3D XY_NORMAL(0.0, 0.0, 1.0);   // ORIGINAL

  const double EPSILON = 0.000001;
  int i, k;
  int N;
  Point3D onset, corner, F0, F1, P0, P1, Q0, Q1, Q, R;
  double yClose;
  double t;
  Point3D C[3];
  double  weight[3];
  BezierCurve3D curve;
  LineStrip3D upperPath, lowerPath;
  double length;
  int numRibs;
  double radius;
  double minZ, maxX;
  Point3D innerOrigin[NUM_LIP_RIBS];    // Where along the teeth row do the individual ribs begin ?
  Point3D outerOrigin[NUM_LIP_RIBS];    // Where along the teeth row do the individual ribs begin ?
  int numFixedRibs;
  int numDivisions[NUM_LIP_RIBS];
  double maxLength;
  Point3D tempInnerOrigin[NUM_LIP_RIBS];
  Point3D tempOuterOrigin[NUM_LIP_RIBS];

  double LIP_RADIUS = anatomy.lipsWidth_cm - 0.4;
  if (LIP_RADIUS < 0.4) { LIP_RADIUS = 0.4; }

  // ****************************************************************
  // ****************************************************************

  getImportantLipPoints(onset, corner, F0, F1, yClose);

  // ****************************************************************
  // Outer edge of the upper and lower lip
  // -> upperPath and lowerPath.
  // ****************************************************************

  upperPath.reset(0);
  lowerPath.reset(0);

  // From lip onset to the lip separation point.

  N = lipCornerPath.getNumPoints();
  for (i=0; i < N-1; i++)
  {
    P0 = lipCornerPath.getControlPoint(i);
    P1 = lipCornerPath.getControlPoint(i+1);
    length = P1.x - P0.x;
    if (length < EPSILON) { length = EPSILON; }

    if ((P0.x <= onset.x) && (P1.x >= onset.x))
    {
      upperPath.addPoint(onset);
      lowerPath.addPoint(onset);
    }

    if ((P0.x > onset.x) && (P0.x < corner.x))
    {
      upperPath.addPoint(P0);
      lowerPath.addPoint(P0);
    }

    if ((P0.x <= corner.x) && (P1.x >= corner.x))
    {
      upperPath.addPoint(corner);
      lowerPath.addPoint(corner);
    }
  }

  // Upper lip from lip corner to the front.

  C[0] = corner;
  C[1] = F0;
  C[2] = Point3D(F0.x, F0.y, 0.0);

  weight[0] = 1.0; weight[1] = 1.5; weight[2] = 1.0;

  curve.setPoints(3, C, weight);
  for (i=1; i <= 7; i++) 
  { 
    upperPath.addPoint(curve.getPoint((double)i / 7.0)); 
  }

  // Lower lip from lip corner to the front.

  C[0] = corner;
  C[1] = F1;
  C[2] = Point3D(F1.x, F1.y, 0.0);

  weight[0] = 1.0; weight[1] = 1.5; weight[2] = 1.0;

  curve.setPoints(3, C, weight);
  for (i=1; i <= 7; i++) 
  { 
    lowerPath.addPoint(curve.getPoint((double)i / 7.0)); 
  }

  // ****************************************************************
  // Ribs of the upper lip.
  // ****************************************************************

  numRibs = 0;
  for (i=1; i < NUM_JAW_RIBS; i++)
  {
    P0 = upperGumsOuterEdge[i-1];
    P1 = upperGumsOuterEdge[i];
    length = P1.x - P0.x;
    if (length < EPSILON) { length = EPSILON; }

    Q0 = upperGumsInnerEdge[i-1];
    Q1 = upperGumsInnerEdge[i];

    if ((P0.x <= onset.x) && (P1.x >= onset.x))
    {
      t = (onset.x - P0.x) / length;
      outerOrigin[numRibs] = P0 + t*(P1-P0);
      innerOrigin[numRibs] = Q0 + t*(Q1-Q0);
      numRibs++;
    }

    if ((P0.x <= corner.x) && (P1.x >= corner.x))
    {
      t = (corner.x - P0.x) / length;
      outerOrigin[numRibs] = P0 + t*(P1-P0);
      innerOrigin[numRibs] = Q0 + t*(Q1-Q0);
      numRibs++;
    }

    if (P1.x > onset.x) 
    { 
      outerOrigin[numRibs] = P1; 
      innerOrigin[numRibs] = Q1; 
      numRibs++;
    }

  }


  // Insert the remaining ribs successively between the "fixed" ribs.

  numFixedRibs = numRibs;
  for (i=0; i < NUM_LIP_RIBS; i++) { numDivisions[i] = 1; }

  while (numRibs < NUM_LIP_RIBS)
  {
    // Find the segment with the greatest length
    k = -1;
    maxLength = -1000000.0;
    for (i=0; i < numFixedRibs-1; i++)
    {
      length = (outerOrigin[i+1].z - outerOrigin[i].z) / numDivisions[i];
      if (length > maxLength)
      {
        maxLength = length;
        k = i;
      }
    }
    numDivisions[k]++;
    numRibs++;
  }

  numRibs = 0;
  for (i=0; i < numFixedRibs; i++)
  {
    if (i < numFixedRibs-1)
    {
      P0 = outerOrigin[i];
      P1 = outerOrigin[i+1];

      Q0 = innerOrigin[i];
      Q1 = innerOrigin[i+1];
      for (k=0; k < numDivisions[i]; k++)
      {
        t = (double)k / (double)numDivisions[i];
        tempOuterOrigin[numRibs] = P0 + (P1-P0)*t;
        tempInnerOrigin[numRibs] = Q0 + (Q1-Q0)*t;
        numRibs++;
      }
    }
    else
    {
      tempOuterOrigin[numRibs] = outerOrigin[i];
      tempInnerOrigin[numRibs] = innerOrigin[i];
      numRibs++;
    }
  }

  for (i=0; i < numRibs; i++) 
  { 
    outerOrigin[i] = tempOuterOrigin[i]; 
    innerOrigin[i] = tempInnerOrigin[i]; 
  }


  // Calculate the points on the upper lip ribs.

  minZ = outerOrigin[0].z;
  maxX = upperPath.getPoint(1.0).x;

  for (i=0; i < numRibs; i++)
  {
    P0 = outerOrigin[i];
    t = (double)i / (double)(numRibs - 1);
    P1 = upperPath.getPoint(t);

    // Generate the inner lip points.

    C[0] = P0;
    C[2] = P1;

    t = 0.7;
    Q = P1; Q.y = P0.y;
    Q = (1.0-t)*Q + t*P1;
  
  	t = 0.0;

    R = P0; R.y = P1.y;
    R = (1.0-t)*R + t*P1;

    t = (P0.x - outerOrigin[0].x) / (outerOrigin[NUM_LIP_RIBS-1].x - outerOrigin[0].x);
    t = t*t*t;
    C[1] = (1.0-t)*Q + t*R;

    weight[0] = 1.0;
    weight[1] = 2.0;
    weight[2] = 1.0;
    curve.setPoints(3, C, weight);

    upperLip->setVertex(i, 0, innerOrigin[i]);
    for (k=0; k < NUM_INNER_LIP_POINTS-1; k++)
    {
      t = (double)k / (double)(NUM_INNER_LIP_POINTS-2);
      upperLip->setVertex(i, k+1, curve.getPoint(t));
    }

    // Generate the outer lip points.

    Q = P1;
    radius = LIP_RADIUS*(Q.x - corner.x) / (maxX - corner.x);
    if (radius < 0.0) { radius = 0.0; }
    Q.y+= radius;

    for (k = 0; k < NUM_OUTER_LIP_POINTS; k++)
    {
      t = -0.5*M_PI + 0.5*M_PI*(double)(k+1) / (double)NUM_OUTER_LIP_POINTS;
      R.x = Q.x + radius*cos(t);
      R.y = Q.y + radius*sin(t);
      R.z = Q.z;

      upperLip->setVertex(i, NUM_INNER_LIP_POINTS + k, R);
    }
  }

  // ****************************************************************
  // Ribs of the lower lip. 
  // ****************************************************************

  numRibs = 0;
  for (i=1; i < NUM_JAW_RIBS; i++)
  {
    P0 = lowerGumsOuterEdge[i-1];
    P1 = lowerGumsOuterEdge[i];
    length = P1.x - P0.x;
    if (length < EPSILON) { length = EPSILON; }

    Q0 = lowerGumsInnerEdge[i-1];
    Q1 = lowerGumsInnerEdge[i];

    if ((P0.x <= onset.x) && (P1.x >= onset.x))
    {
      t = (onset.x - P0.x) / length;
      outerOrigin[numRibs] = P0 + t*(P1-P0);
      innerOrigin[numRibs] = Q0 + t*(Q1-Q0);
      numRibs++;
    }

    if ((P0.x <= corner.x) && (P1.x >= corner.x))
    {
      t = (corner.x - P0.x) / length;
      outerOrigin[numRibs] = P0 + t*(P1-P0);
      innerOrigin[numRibs] = Q0 + t*(Q1-Q0);
      numRibs++;
    }

    if (P1.x > onset.x) 
    { 
      outerOrigin[numRibs] = P1; 
      innerOrigin[numRibs] = Q1; 
      numRibs++;
    }
  }

  // Insert the remaining ribs successively between the "fixed" ribs.

  numFixedRibs = numRibs;
  for (i=0; i < NUM_LIP_RIBS; i++) { numDivisions[i] = 1; }

  while (numRibs < NUM_LIP_RIBS)
  {
    // Find the segment with the greatest length.
    k = -1;
    maxLength = -1000000.0;
    for (i=0; i < numFixedRibs-1; i++)
    {
      length = (outerOrigin[i+1].z - outerOrigin[i].z) / numDivisions[i];
      if (length > maxLength)
      {
        maxLength = length;
        k = i;
      }
    }
    numDivisions[k]++;
    numRibs++;
  }

  numRibs = 0;
  for (i=0; i < numFixedRibs; i++)
  {
    if (i < numFixedRibs-1)
    {
      P0 = outerOrigin[i];
      P1 = outerOrigin[i+1];

      Q0 = innerOrigin[i];
      Q1 = innerOrigin[i+1];
      for (k=0; k < numDivisions[i]; k++)
      {
        t = (double)k / (double)numDivisions[i];
        tempOuterOrigin[numRibs] = P0 + (P1-P0)*t;
        tempInnerOrigin[numRibs] = Q0 + (Q1-Q0)*t;
        numRibs++;
      }
    }
    else
    {
      tempOuterOrigin[numRibs] = outerOrigin[i];
      tempInnerOrigin[numRibs] = innerOrigin[i];
      numRibs++;
    }
  }

  for (i=0; i < numRibs; i++) 
  { 
    outerOrigin[i] = tempOuterOrigin[i]; 
    innerOrigin[i] = tempInnerOrigin[i]; 
  }

  // Calculate the points on the lower lip ribs.

  minZ = outerOrigin[0].z;
  maxX = lowerPath.getPoint(1.0).x;

  for (i=0; i < numRibs; i++)
  {
    P0 = outerOrigin[i];
    t = (double)i / (double)(numRibs - 1);
    P1 = lowerPath.getPoint(t);

    // The inner lip points.

    C[0] = P0;
    C[2] = P1;

    t = 0.7;
    Q = P1; Q.y = P0.y;
    Q = (1.0-t)*Q + t*P1;
    
    t = 0.0;

    R = P0; R.y = P1.y;
    R = (1.0-t)*R + t*P1;

    t = (P0.x - outerOrigin[0].x) / (outerOrigin[NUM_LIP_RIBS-1].x - outerOrigin[0].x);
    t = t*t*t;
    C[1] = (1.0-t)*Q + t*R;

    weight[0] = 1.0;
    weight[1] = 2.0;
    weight[2] = 1.0;
    curve.setPoints(3, C, weight);

    lowerLip->setVertex(i, 0, innerOrigin[i]);
    for (k=0; k < NUM_INNER_LIP_POINTS-1; k++)
    {
      t = (double)k / (double)(NUM_INNER_LIP_POINTS-2);
      lowerLip->setVertex(i, k+1, curve.getPoint(t));
    }

    // The outer lip points.

    Q = P1;
    radius = LIP_RADIUS*(Q.x - corner.x) / (maxX - corner.x);
    if (radius < 0.0) { radius = 0.0; }
    Q.y-= radius;

    for (k = 0; k < NUM_OUTER_LIP_POINTS; k++)
    {
      t = 0.5*M_PI - 0.5*M_PI*(double)(k+1) / (double)NUM_OUTER_LIP_POINTS;
      R.x = Q.x + radius*cos(t);
      R.y = Q.y + radius*sin(t);
      R.z = Q.z;

      lowerLip->setVertex(i, NUM_INNER_LIP_POINTS + k, R);
    }
  }

  // ****************************************************************
  // Calculate the radiation surface.
  // ****************************************************************

  calcRadiation(corner);
}


// ****************************************************************************
// Calculates some important lip "landmarks", from which the whole lip surfaces
// can be constructed.
// ****************************************************************************

void VocalTract::getImportantLipPoints(Point3D &onset, Point3D &corner, 
                                       Point3D &F0, Point3D &F1, double &yClose)
{
  int i;
  int N = narrowLipCornerPath.getNumPoints();

  // ****************************************************************
  // Interpolation between the narrow and the wide lip corner path 
  // and transformation into the global coordinate system.
  // ****************************************************************

  double angle_rad = 0.5 * param[JA].x*M_PI/180.0;
  double cosinus = cos(angle_rad);
  double sinus   = sin(angle_rad);
  double t, t0; 
  double dx, dy;
  Point3D Q, R;

  lipCornerPath.reset(N);

  t = (param[LD].x - param[LD].min) / (param[LD].max - param[LD].min);
  for (i=0; i < N; i++)
  {
    Q  = (1.0-t)*narrowLipCornerPath.getControlPoint(i) + t*wideLipCornerPath.getControlPoint(i);

    // (dx, dy) is the vertex position relativ to the fulcrum before the rotation.
    dx = Q.x + anatomy.jawRestPos.x + 0.5*param[JX].x - anatomy.jawFulcrum.x;
    dy = Q.y + 0.5*anatomy.jawRestPos.y - anatomy.jawFulcrum.y;
    R.x = cosinus*dx - sinus*dy + anatomy.jawFulcrum.x;
    R.y = sinus*dx + cosinus*dy + anatomy.jawFulcrum.y;
    R.z = Q.z;

    lipCornerPath.setPoint(i, R);
  }

  // Curve parameter, at which the lips separate from the teeth.

  double lipCornerSeparationX = upperGumsOuterEdge[6].x;    
  if (lowerGumsOuterEdge[6].x < lipCornerSeparationX)
  {
    lipCornerSeparationX = lowerGumsOuterEdge[6].x;
  }

  t0 = narrowLipCornerPath.getIntersection(Point3D(lipCornerSeparationX, 0.0, 0.0), Point3D(1.0, 0.0, 0.0));

  // Curve parameter for the lip corner point.

  t = (param[LP].x - param[LP].min) / (param[LP].max - param[LP].min);

  // The lip corner point.

  corner = lipCornerPath.getPoint(t);

  // Onset of the lip surface.

  if (t > t0) 
  { 
    onset = lipCornerPath.getPoint(t0); 
  } 
  else 
  { 
    onset = corner; 
  }

  // ****************************************************************
  // The important points on the upper and lower lip.
  // ****************************************************************

  Point3D T0, T1;
  double L;

  T0 = upperGumsOuterEdge[8];    // In the midsagittal plane
  T1 = lowerGumsOuterEdge[8];

  const double X_OFFSET = 0.1;
  double Y_OFFSET = 0.35 - t*0.3;
  const double L_MAX = 1.0;
  const double L_MIN = 0.3;

  // x-coordinate.

  t = (param[LP].x - param[LP].min) / (param[LP].max - param[LP].min);
  L = (1.0-t)*L_MAX + t*L_MIN;

  F0.x = corner.x + L + 0.3;
  F1.x = corner.x + L;
  if (F0.x < T0.x + X_OFFSET) { F0.x = T0.x + X_OFFSET; }
  if (F1.x < T1.x + X_OFFSET) { F1.x = T1.x + X_OFFSET; }

  // y-coordinate.
  
  const double MIN_T = -0.05;

  t = (param[LP].x - param[LP].min) / (param[LP].max - param[LP].min);
  yClose = 0.5*(T0.y + T1.y) + (1.0 - t)*Y_OFFSET;

  t = 0.5*param[LD].x;
  // Full lip closure ?
  if (t < MIN_T) 
  { 
    t = MIN_T; 
  }
  F0.y = yClose + t;
  F1.y = yClose - t;

  if (F0.y > T0.y) { F0.y = T0.y; }
  if (F1.y > T0.y) { F1.y = T0.y; }
  if (F0.y < T1.y) { F0.y = T1.y; }
  if (F1.y < T1.y) { F1.y = T1.y; }

  // z-coordinate.

  Point3D lipTangent = upperGumsOuterEdge[5] - upperGumsOuterEdge[0];
  F0.z = corner.z + (lipTangent.z/lipTangent.x)*(F0.x - corner.x);
  F1.z = corner.z + (lipTangent.z/lipTangent.x)*(F1.x - corner.x);
}


// ****************************************************************************
// Calculate the radiation surface.
// ****************************************************************************

void VocalTract::calcRadiation(Point3D lipCorner)
{
  int i, k;
  Point3D v, n;
  Point3D U, L, M, P;
  Surface *upperLip = &surface[UPPER_LIP];
  Surface *lowerLip = &surface[LOWER_LIP];
  int lipPoint = NUM_INNER_LIP_POINTS;
  double h = 0.0;
  double x[NUM_RADIATION_POINTS];
  double y[NUM_RADIATION_POINTS];
  double angle_rad;
  double minZ = 0.0, newZ = 0.0;

  for (i=0; i < NUM_RADIATION_POINTS; i++)
  {
    angle_rad = -0.5*M_PI + M_PI*(double)i/(double)(NUM_RADIATION_POINTS-1);
    x[i] = cos(angle_rad);
    y[i] = 0.5*sin(angle_rad) + 0.5;
  }

  for (i=0; i < NUM_RADIATION_RIBS; i++)
  {
    if (i < NUM_LIP_RIBS)
    {
      n.set(0.0, 0.0, -1.0);

      U = upperLip->getVertex(i, lipPoint);
      L = lowerLip->getVertex(i, lipPoint);
      M = 0.5*(U + L);
      v = U - L;
      h = U.y - L.y;
      if (h <= 0.0) { h = 0.0; }

      if (U.x <= lipCorner.x)
      {
        h = 0;
      }
      else
      {
        newZ = M.z - h;
        if (newZ > minZ) { newZ = minZ; } else { minZ = newZ; }
        h = M.z - minZ;
      }
    }
    else
    {
      angle_rad = 0.5*M_PI - 0.5*M_PI*(double)(i - NUM_LIP_RIBS + 1) / (double)(NUM_RADIATION_RIBS - NUM_LIP_RIBS);
      n.set(cos(angle_rad), 0.0, -sin(angle_rad));
      // Keep L, v and n from the last lip rib.
    }

    for (k=0; k < NUM_RADIATION_POINTS; k++)
    {
      P = L + y[k]*v + x[k]*h*n;
      surface[RADIATION].setVertex(i, k, P);
    }
  }
}


// ****************************************************************************
/// Calculates the tanget on the posterior side of the tongue body circle 
/// that runs through the hyoid.
/// H is the position of the hyoid and T the point of contact with the 
/// tongue circle.
// ****************************************************************************

void VocalTract::getHyoidTongueTangent(Point2D &H, Point2D &T)
{
  // H is the hyoid and C and rx/ry are the center and radius
  // of the tongue body circles.

  Point2D C(param[TCX].limitedX, param[TCY].limitedX);
  double rx = anatomy.tongueCenterRadiusX_cm;
  double ry = anatomy.tongueCenterRadiusY_cm;
  double alpha;

  // ****************************************************************
  // The Bezier curve for the tongue root.
  // ****************************************************************

  H = surface[LOWER_COVER].getVertex(NUM_LARYNX_RIBS-1, NUM_LOWER_COVER_POINTS-1).toPoint2D();
  
  alpha = getEllipseTangent(H, C, rx, ry, true);

  T.x = C.x + rx*cos(alpha);
  T.y = C.y + ry*sin(alpha);
}


// ****************************************************************************
// Calculate the 3D-shape of the tongue.
// ****************************************************************************

void VocalTract::calcTongue()
{
  int i;
  TongueRib *rib = &tongueRib[0];
  const double EPSILON = 0.000001;
  LineStrip2D targetCurve;
  BezierCurve3D rootCurve, bladeCurve;
  double alpha[4];
  double delta, t;
  double angle;

  // If the parameters (TRX, TRY) for the tongue root are supposed to
  // be calculated automatically, do it here.

  if (anatomy.automaticTongueRootCalc)
  {
    // The tongue body center position.
    Point2D C(param[TCX].limitedX, param[TCY].limitedX);
    // The position of the hyoid.
    Point2D H = surface[LOWER_COVER].getVertex(NUM_LARYNX_RIBS-1, NUM_LOWER_COVER_POINTS-1).toPoint2D();
    double distance = (C - H).magnitude();

    param[TRX].x = anatomy.tongueRootTrxSlope * distance + anatomy.tongueRootTrxIntercept;
    param[TRY].x = anatomy.tongueRootTrySlope * param[TCX].x + anatomy.tongueRootTryIntercept;
  }

  // Restrict the tongue parameters so that the tongue does not extend
  // to much beyond the vocal tract covers.

  restrictTongueParams();

  // H is the hyoid and C0, C1, r0 and r1 are the centers and radii 
  // of the two tongue circles.

  Point2D C, D;
  Point2D C0(param[TCX].limitedX, param[TCY].limitedX);
  Point2D C1(param[TTX].limitedX, param[TTY].limitedX);
  double r0x = anatomy.tongueCenterRadiusX_cm;
  double r0y = anatomy.tongueCenterRadiusY_cm;
  double r1 = anatomy.tongueTipRadius_cm;
  Point3D Q[3];
  const double WEIGHT[3] = { 1.0, 2.0, 1.0 };

  const int NUM_TONGUE_SIDE_POINTS = 5;
  int    tongueSideIndex[NUM_TONGUE_SIDE_POINTS];
  double tongueSideElevation_cm[NUM_TONGUE_SIDE_POINTS];

  LineStrip2D leftSideHeight;
  LineStrip2D rightSideHeight;

  // ****************************************************************
  // The Bezier curve for the tongue root.
  // ****************************************************************

  Q[0] = surface[LOWER_COVER].getVertex(NUM_LARYNX_RIBS-1, NUM_LOWER_COVER_POINTS-1);
  
  alpha[0] = getEllipseTangent(Q[0].toPoint2D(), C0, r0x, r0y, true);

  Q[2].x = C0.x + r0x*cos(alpha[0]);
  Q[2].y = C0.y + r0y*sin(alpha[0]);
  Q[2].z = 0.0;
  
  Q[1].x = param[TRX].x;
  Q[1].y = param[TRY].x;

  alpha[0] = getEllipseTangent(Q[1].toPoint2D(), C0, r0x, r0y, true);

  Q[2].x = C0.x + r0x*cos(alpha[0]);
  Q[2].y = C0.y + r0y*sin(alpha[0]);
  Q[2].z = 0.0;

  rootCurve.setPoints(3, &Q[0], &WEIGHT[0]);

  // ****************************************************************
  // The Bezier curve for the tongue blade.
  // ****************************************************************

  Q[1].set(param[TBX].x, param[TBY].x, 0.0);

  alpha[1] = getEllipseTangent(Q[1].toPoint2D(), C0, r0x, r0y, false);

  // Make sure that alpha[1] is always smaller/equal than alpha[0]!
  if ((alpha[1] > alpha[0]) && (alpha[1] - alpha[0] < 0.25*M_PI))
  {
    alpha[1] = alpha[0];
  }
  Q[0].x = C0.x + r0x*cos(alpha[1]);
  Q[0].y = C0.y + r0y*sin(alpha[1]);
  Q[0].z = 0.0;


  alpha[2] = getCircleTangent(Q[1].toPoint2D(), C1, r1, true);
  Q[2].x = C1.x + r1*cos(alpha[2]);
  Q[2].y = C1.y + r1*sin(alpha[2]);
  Q[2].z = 0.0;

  bladeCurve.setPoints(3, &Q[0], &WEIGHT[0]);

  // ****************************************************************
  // Calculation of the entire mid-sagittal tongue contour from the 
  // four individual segments.
  // ****************************************************************

  targetCurve.reset(0);
  tongueSideIndex[0] = 0;

  // The Bezier curve for the tongue root.

  for (i=0; i < 32; i++)
  {
    targetCurve.addPoint(rootCurve.getPoint((double)i/32.0).toPoint2D());
  }
  tongueSideIndex[1] = targetCurve.getNumPoints() - 8; // 16;

  // Circle segment for the tongue back.

  if (alpha[0] < 0.0) { alpha[0]+= 2.0*M_PI; }
  delta = alpha[1] - alpha[0];
  if (delta > -EPSILON) { delta = -EPSILON; }

  for (i=0; i < 32; i++)
  {
    angle = alpha[0] + delta*(double)i / 32.0;
    C.x = C0.x + r0x*cos(angle);
    C.y = C0.y + r0y*sin(angle);
    targetCurve.addPoint(C);
  }
  tongueSideIndex[2] = targetCurve.getNumPoints() - 8; // 16;

  // The Bezier curve for the tongue blade.

  for (i=0; i < 32; i++)
  {
    targetCurve.addPoint(bladeCurve.getPoint((double)i/32.0).toPoint2D());
  }
//  tongueSideIndex[3] = targetCurve.getNumPoints() - 16;
  tongueSideIndex[3] = targetCurve.getNumPoints() - 1;

  // Circle segment for the tongue tip.

  alpha[3] = 0.0;
  if (alpha[2] < 0.0) { alpha[2]+= 2.0*M_PI; }
  delta = alpha[3] - alpha[2];
  if (delta > -EPSILON) { delta = -EPSILON; }

  for (i=0; i < 8; i++)
  {
    angle = alpha[2] + delta*(double)i / 7.0;
    C.x = C1.x + r1*cos(angle);
    C.y = C1.y + r1*sin(angle);
    targetCurve.addPoint(C);
  }
  tongueSideIndex[4] = targetCurve.getNumPoints() - 1;
  
  // ****************************************************************
  // Left tongue side spline.
  // ****************************************************************

  tongueSideElevation_cm[0] = 1.0;
  tongueSideElevation_cm[1] = tongueSideParamToElevation_cm(param[TS1].x);
  tongueSideElevation_cm[2] = tongueSideParamToElevation_cm(param[TS2].x);
  tongueSideElevation_cm[3] = tongueSideParamToElevation_cm(param[TS3].x);
  tongueSideElevation_cm[4] = tongueSideElevation_cm[3];
  if (tongueSideElevation_cm[4] > 0.0) { tongueSideElevation_cm[4] = 0.0; }

  leftSideHeight.reset(0);
  for (i=0; i < NUM_TONGUE_SIDE_POINTS; i++)
  {
    t = targetCurve.getCurveParam(tongueSideIndex[i]);
    leftSideHeight.addPoint(Point2D(t, tongueSideElevation_cm[i]));
  }

  // ****************************************************************
  // Right tongue side spline.
  // ****************************************************************

  rightSideHeight.reset(0);
  for (i=0; i < NUM_TONGUE_SIDE_POINTS; i++)
  {
    t = targetCurve.getCurveParam(tongueSideIndex[i]);
    rightSideHeight.addPoint(Point2D(t, tongueSideElevation_cm[i]));
  }

  // ****************************************************************
  // The individual points in the midsagittal plane may not exceed
  // the contour of the vocal tract hull.
  // ****************************************************************

  // Center of a circle segment.
  Point2D center(0.0, -1.5);
  // Normal vector for the translation of a tongue point.
  Point2D n;    
  Point2D intersection;
  int N = targetCurve.getNumPoints();

  for (i=0; i < N; i++)
  {
    n.set(0.0, 0.0);

    C = targetCurve.getControlPoint(i);
   
    if ((C.x <= center.x) && (C.y <= center.y)) { n.set(1.0, 0.0); } else
    if ((C.x >= center.x) && (C.y >= center.y)) { n.set(0.0, -1.0); } else
    if ((C.x <= center.x) && (C.y >= center.y)) { n = center - C; }

    if ((n.x != 0.0) || (n.y != 0.0))
    {
      n.normalize();

      if (upperOutline.getSpecialIntersection(C, n, t, intersection))
      {
        if (t > 0.0) 
        { 
          targetCurve.setPoint(i, intersection); 
        }
      }
    }
  }

  // ****************************************************************
  // Calculate the tongue ribs.
  // ****************************************************************

  for (i=0; i < NUM_DYNAMIC_TONGUE_RIBS; i++)
  {
    t = (double)i/(double)(NUM_DYNAMIC_TONGUE_RIBS-1);
    rib[i].point = targetCurve.getPoint(t);
    rib[i].leftSideHeight  = leftSideHeight.getFunctionValue(t);
    rib[i].rightSideHeight = rightSideHeight.getFunctionValue(t);
  }

  // ****************************************************************
  // Low-pass filtering of the left and right tongue sides.
  // ****************************************************************

  const double FILTER_COEFF = 0.5;

  for (i=1; i < NUM_DYNAMIC_TONGUE_RIBS; i++)
  {
    rib[i].leftSideHeight+= FILTER_COEFF*(rib[i-1].leftSideHeight - rib[i].leftSideHeight);
    rib[i].rightSideHeight+= FILTER_COEFF*(rib[i-1].rightSideHeight - rib[i].rightSideHeight);
  }

  // ****************************************************************
  // Get the normal vectors of the tongue ribs.
  // ****************************************************************

  N = 4;
  Point2D n0, n1;

  for (i=N; i < NUM_DYNAMIC_TONGUE_RIBS-N; i++)
  {
    rib[i].normal = (rib[i+1].point - rib[i-1].point).turnLeft();
    rib[i].normal.normalize();
  }

  // Adjust the normals at the tongue root.

  n0.set(-1.0, 0.0);
  n1 = rib[N].normal;
  for (i=0; i < N; i++)
  {
    t = (double)i / (double)N;
    rib[i].normal = (1.0-t)*n0 + t*n1;
  }

  // Adjust the normals at the tongue tip.

  n0 = rib[NUM_DYNAMIC_TONGUE_RIBS-N-1].normal;
  n1.set(0.0, 1.0);
  for (i=0; i < N; i++)
  {
    t = (double)(i+1) / (double)N;
    rib[NUM_DYNAMIC_TONGUE_RIBS-N+i].normal = (1.0-t)*n0 + t*n1;
  }

  // ****************************************************************
  // Run through all tongue ribs forwards and backwards and correct
  // any possible intersections between neighboring ribs.
  // ****************************************************************

  for (i=0; i < NUM_DYNAMIC_TONGUE_RIBS; i++)
  {
    rib[i].minY = param[VocalTract::TS1].min - 0.5;
    rib[i].maxY = param[VocalTract::TS1].max + 0.5;
  }

  Point2D origNormal[NUM_TONGUE_RIBS];
  Point2D newNormal[NUM_TONGUE_RIBS];

  for (i=0; i < NUM_DYNAMIC_TONGUE_RIBS; i++) 
  { 
    origNormal[i] = rib[i].normal; 
  }
  
  // Do the adjustment in a forward run.

  for (i=0; i < NUM_DYNAMIC_TONGUE_RIBS-1; i++) 
  { 
    verifyTongueRibNormal(i, i+1); 
  }
 
  for (i=0; i < NUM_DYNAMIC_TONGUE_RIBS; i++) 
  {
    newNormal[i] = rib[i].normal;
    rib[i].normal = origNormal[i];
  }

  // Do the adjustment in a reverse run.

  for (i=NUM_DYNAMIC_TONGUE_RIBS-1; i > 0; i--) 
  { 
    verifyTongueRibNormal(i, i-1); 
  }

  for (i=0; i < NUM_DYNAMIC_TONGUE_RIBS; i++) 
  {
    newNormal[i]+= rib[i].normal;
    rib[i].normal = newNormal[i].normalize();
  }


  // ****************************************************************
  // Determine the lower boundary including the mouth floor.
  // ****************************************************************

  LineStrip2D lowerBoundary;

  lowerBoundary.reset(0);
  C = surface[LOWER_COVER].getVertex(NUM_LARYNX_RIBS + NUM_THROAT_RIBS, NUM_LOWER_COVER_POINTS-1).toPoint2D();
  C.y+= 0.5;
  lowerBoundary.addPoint(C);

  if (lowerOutline.getClosestIntersection(C, Point2D(1.0, 0.0), t, intersection))
  {
    lowerBoundary.addPoint(intersection);
  }
  else
  {
    lowerBoundary.addPoint(C);              // Error case !!
  }

  double minY = lowerBoundary.getControlPoint(1).y;
  for (i = NUM_LARYNX_RIBS + NUM_THROAT_RIBS; i < NUM_LOWER_COVER_RIBS; i++)
  {
    C = surface[LOWER_COVER].getVertex(i, NUM_LOWER_COVER_POINTS-1).toPoint2D();
    if (C.y > minY) { lowerBoundary.addPoint(C); }
  }

  for (i=0; i < 3; i++)
  {
    lowerBoundary.addPoint(surface[LOWER_TEETH].getVertex(NUM_TEETH_RIBS-1, i).toPoint2D());
  }

  C = lowerBoundary.getPoint(lowerBoundary.getNumPoints()-1);
  C.x+= 10.0;
  lowerBoundary.addPoint(C);

  // ****************************************************************
  // Calculate the rib that finishes the dynamic ribs.
  // ****************************************************************

  double tongueTipRadius = anatomy.tongueTipRadius_cm;

  i = NUM_DYNAMIC_TONGUE_RIBS-1;
  if (rib[i].leftSideHeight < rib[i].rightSideHeight) 
  { 
    t = rib[i].leftSideHeight; 
  } 
  else 
  { 
    t = rib[i].rightSideHeight; 
  }
  Point2D finalPoint = rib[i].point + rib[i].normal*t;

  if (param[TTY].limitedX - tongueTipRadius < finalPoint.y) 
  { 
    finalPoint.y = param[TTY].limitedX - tongueTipRadius; 
  }

  i = NUM_DYNAMIC_TONGUE_RIBS;
  rib[i].normal.set(0.0, 1.0);
  rib[i].leftSideHeight = 0.0;
  rib[i].rightSideHeight = 0.0;
  rib[i].point = finalPoint;

  n.set(0.0, 1.0);
  for (i=NUM_TONGUE_RIBS/2; i <= NUM_DYNAMIC_TONGUE_RIBS; i++)
  {
    // The tongue body points.
    if (lowerBoundary.getClosestIntersection(rib[i].point, n, t, intersection))
    {
      if (t > 0.0) { rib[i].point = intersection; }
    }

    // The tongue sides.
    if (lowerBoundary.getClosestIntersection(rib[i].point, rib[i].normal, t, intersection))
    {
      if (rib[i].leftSideHeight < t)  { rib[i].leftSideHeight = t; }
      if (rib[i].rightSideHeight < t) { rib[i].rightSideHeight = t; }
    }

  }

  // ****************************************************************
  // Calculate the most anterior three ribs.
  // ****************************************************************

  // Die drittletzte Rippe im 45-Grad-Winkel setzen.

  C = rib[NUM_TONGUE_RIBS-4].point;
  D = C;
  D.x-= 2.0*tongueTipRadius;
  D.y-= 2.0*tongueTipRadius;

  if (lowerBoundary.getClosestIntersection(C, D-C, t, intersection))
  {
    if (t < 1.0) { D = intersection; }
  }

  i = NUM_TONGUE_RIBS-3;
  rib[i].point = D;
  rib[i].normal.set(0.0, 1.0);
  rib[i].leftSideHeight = 0.0;
  rib[i].rightSideHeight = 0.0;

  // Define the two last tongue ribs.

  C = rib[NUM_TONGUE_RIBS-3].point;
  if (lowerBoundary.getClosestIntersection(C, Point2D(0.0, -1.0), t, intersection))
  {
    rib[NUM_TONGUE_RIBS-2].point = intersection;
  }
  else
  {
    rib[NUM_TONGUE_RIBS-2].point = C;    // Error case.
  }

  if (rib[NUM_TONGUE_RIBS-2].point.x < lowerBoundary.getControlPoint(1).x)
  {
    rib[NUM_TONGUE_RIBS-1].point = lowerBoundary.getControlPoint(1);
  }
  else
  {
    rib[NUM_TONGUE_RIBS-1].point = rib[NUM_TONGUE_RIBS-2].point;
  }

  i = NUM_TONGUE_RIBS-2;
  rib[i].normal.set(0.0, 1.0);
  rib[i].leftSideHeight = 0.0;
  rib[i].rightSideHeight = 0.0;

  i = NUM_TONGUE_RIBS-1;
  rib[i].normal.set(0.0, 1.0);
  rib[i].leftSideHeight = 0.0;
  rib[i].rightSideHeight = 0.0;

  // ****************************************************************
  // Calculate the actual rib shapes from their descriptions.
  // ****************************************************************

  calcTongueRibs();

  // ****************************************************************
  // Add the bulging at the tongue root (as for our model speaker JD).
  // ****************************************************************

  const int NUM_BULGE_RIBS = 8;
  const double BULGE_HEIGHT_CM = 0.5;
  Surface *s = &surface[TONGUE];
  double b;
  int k;
  Point3D P;

  for (i=0; i < NUM_BULGE_RIBS; i++)
  {
    n = rib[i].normal;
    for (k=0; k < s->numRibPoints; k++)
    {
      b = BULGE_HEIGHT_CM*(0.5-0.5*cos(2.0*M_PI*(double)i/(double)(NUM_BULGE_RIBS-1)))*
                          (0.5-0.5*cos(2.0*M_PI*(double)k/(double)(s->numRibPoints-1)));
      P = s->getVertex(i, k);
      P.x+= n.x*b;
      P.y+= n.y*b;
      s->setVertex(i, k, P);
    }
  }
}


// ****************************************************************************
/// Returns the elevation of the tongue side with respect to the midline for
/// a tongue side parameter value. The parameter value is in [-1 ... +1].
// ****************************************************************************

double VocalTract::tongueSideParamToElevation_cm(double paramValue)
{
  double elevation_cm = 0.0;

  // Positive values raise the tongue sides.
  if (paramValue >= 0.0)
  {
    // For the parameter values between 0 and 0.3, raise the tongue
    // side from 0 to 1.0 cm. For higher values of the parameter
    // the elevation keeps constant and only the "bracing force" 
    // increases.

    const double MAX_ELEVATION_CM = 1.0;
    elevation_cm = MAX_ELEVATION_CM * (paramValue / 0.3);
    if (elevation_cm > MAX_ELEVATION_CM)
    {
      elevation_cm = MAX_ELEVATION_CM;
    }
  }
  else
  
  // Negative values lower the tongue sides (only for laterals).
  {
    // For the parameter values between 0 and -0.3, lower the tongue
    // side from 0 to -1.0 cm. For lower values of the parameter
    // the elevation keeps constant and only the "opening force" 
    // increases.

    const double MIN_ELEVATION_CM = -1.0;
    elevation_cm = MIN_ELEVATION_CM * (paramValue / -0.3);
    if (elevation_cm < MIN_ELEVATION_CM)
    {
      elevation_cm = MIN_ELEVATION_CM;
    }
  }

  return elevation_cm;
}


// ****************************************************************************
// Returns the "minimum area" value for a given tongue side parameter value.
// The minimum area ensures a certain minimum area in the area function
// to compensate for inaccuracies due to the discrete triangle surfaces of the 
// tongue and palate.
// ****************************************************************************

double VocalTract::tongueSideParamToMinArea_cm2(double paramValue)
{
  double minArea_cm2 = 0.0;

  // Parameter values between 0.2 and 0.4 increase the bracing force
  // and so the min. area up to 0.15 cm^2 (which is needed for fricatives).

  if (paramValue > 0.2)
  {
    const double MAX_AREA_CM2 = 0.15;
    minArea_cm2 = MAX_AREA_CM2 * (paramValue - 0.2) / 0.2;
    if (minArea_cm2 > MAX_AREA_CM2)
    {
      minArea_cm2 = MAX_AREA_CM2;
    }
  }
  else

  // Parameter values between -0.05 and -0.25 (for laterals) increase 
  // the bracing force and so the lateral min. area up to 0.25 cm^2 
  // (which is needed for laterals).

  if (paramValue < -0.05)
  {
    const double MAX_AREA_CM2 = 0.25;
    minArea_cm2 = MAX_AREA_CM2 * (paramValue + 0.05) / (-0.2);
    if (minArea_cm2 > MAX_AREA_CM2)
    {
      minArea_cm2 = MAX_AREA_CM2;
    }
  }

  return minArea_cm2;
}


// ****************************************************************************
/// Constraints the tongue parameters so that the tongue does not stick out
/// to much beyond the vocal tract hull.
// ****************************************************************************

void VocalTract::restrictTongueParams()
{
  const double EPSILON = 0.000001;
  const double DELTA = 0.3;   // The tongue circles can go 3 mm beyond the hull.
  Point2D C, P0, P1, P2;
  LineStrip2D upperBorder;
  LineStrip2D lowerBorderTT;    // Lower boundary for the tongue tip circle.
  LineStrip2D lowerBorderTB;   // Lower boundary for the tongue body circle.
  int i;
  int N;
  double r;
  double tongueTipRadius = anatomy.tongueTipRadius_cm;


  // Keep in mind the user defined values of some vocal tract parameters
  
  double origTTX = param[TTX].x;
  double origTTY = param[TTY].x;
  double origTCX = param[TCX].x;
  double origTCY = param[TCY].x;

  // Pre-limitation of the tongue circle parameter values.

  restrictParam(TCX);
  restrictParam(TCY);
  restrictParam(TTX);
  restrictParam(TTY);

  // ****************************************************************
  // Define the point A towards which the circles retract when they 
  // are outside the upper border.
  // ****************************************************************

  Point2D Q0 = surface[UPPER_COVER].getVertex(NUM_UPPER_COVER_RIBS-NUM_PALATE_RIBS+1, NUM_UPPER_COVER_POINTS-1).toPoint2D();
  Point2D Q1 = surface[UPPER_COVER].getVertex(NUM_UPPER_COVER_RIBS-1, NUM_UPPER_COVER_POINTS-1).toPoint2D();
  Point2D A(Q0.x, Q1.y - (Q1.x - Q0.x));

  // ****************************************************************
  // The upper border for the tongue circles is combined of the
  // rear pharyngeal wall, the palatal outline, the incisors, and
  // a vertical anterior line through the incisors.
  // ****************************************************************

  for (i=NUM_LARYNX_RIBS; i < NUM_UPPER_COVER_RIBS; i++)
  {
    upperBorder.addPoint(surface[UPPER_COVER].getVertex(i, NUM_UPPER_COVER_POINTS-1).toPoint2D());
  }
  for (i=0; i < 2; i++)
  {
    upperBorder.addPoint(surface[UPPER_TEETH].getVertex(NUM_TEETH_RIBS-1, i).toPoint2D());
  }

  P0 = upperBorder.getControlPoint( upperBorder.getNumPoints() - 1 );
  upperBorder.addPoint( Point2D(P0.x, P0.y - 10.0) );

  // Extend the upper border a little bit radially from A.
  N = upperBorder.getNumPoints();
  for (i=0; i < N; i++)
  {
    C = upperBorder.getControlPoint(i);
    C+= (C - A).normalize()*DELTA;
    if (C.x > P0.x)
    {
      C.x = P0.x;
    }
    upperBorder.setPoint(i, C);
  }

  // ****************************************************************
  // Lower border for the tongue tip.
  // ****************************************************************

  lowerBorderTT.reset(1);
  for (i=NUM_LARYNX_RIBS - 1; i < NUM_LOWER_COVER_RIBS; i++)
  {
    lowerBorderTT.addPoint(surface[LOWER_COVER].getVertex(i, NUM_LOWER_COVER_POINTS-1).toPoint2D());
  }

  for (i=0; i < 3; i++)
  {
    lowerBorderTT.addPoint(surface[LOWER_TEETH].getVertex(NUM_TEETH_RIBS-1, i).toPoint2D());
  }

  C = lowerBorderTT.getControlPoint(1);
  C.x-= 10.0;
  lowerBorderTT.setPoint(0, C);

  // ****************************************************************
  // Lower border for the tongue body.
  // ****************************************************************

  lowerBorderTB.reset(1);

  P0 = surface[LOWER_COVER].getVertex(NUM_LARYNX_RIBS - 1, NUM_LOWER_COVER_POINTS - 1).toPoint2D();
  P2 = surface[LOWER_COVER].getVertex(NUM_LOWER_COVER_RIBS - 1, NUM_LOWER_COVER_POINTS - 1).toPoint2D();
  P1 = Point2D(P2.x - 1.0, P0.y);

  if (P1.x < P0.x) { P1.x = P0.x; }

  lowerBorderTB.addPoint(P0);
  lowerBorderTB.addPoint(P1);
  lowerBorderTB.addPoint(P2);

  for (i = 0; i < 3; i++)
  {
    lowerBorderTB.addPoint(surface[LOWER_TEETH].getVertex(NUM_TEETH_RIBS - 1, i).toPoint2D());
  }

  C = lowerBorderTB.getControlPoint(1);
  C.x -= 10.0;
  lowerBorderTB.setPoint(0, C);


  // ****************************************************************
  // Restrict the tongue body and the tongue tip positions.
  // ****************************************************************

  double rx = anatomy.tongueCenterRadiusX_cm;
  double ry = anatomy.tongueCenterRadiusY_cm;

  // Tongue tip must be right of the tongue body.
  if (param[TTX].x < param[TCX].x + rx + tongueTipRadius)
  {
    param[TTX].x = param[TCX].x + rx + tongueTipRadius;
  }

  // Limit the tongue tip position with respect to the border.
  C = limitEllipsePos( Point2D(param[TTX].x, param[TTY].x), 
    tongueTipRadius, tongueTipRadius, upperBorder, A);
  param[TTX].x = C.x;
  param[TTY].x = C.y;

  // The tongue body must be left of the tongue tip.
  if (param[TCX].x > param[TTX].x - rx - tongueTipRadius)
  {
    param[TCX].x = param[TTX].x - rx - tongueTipRadius;
  }

  // Limit the tongue body position with respect to the border.
  C = limitEllipsePos( Point2D(param[TCX].x, param[TCY].x), 
    anatomy.tongueCenterRadiusX_cm, anatomy.tongueCenterRadiusY_cm, upperBorder, A);
  param[TCX].x = C.x;
  param[TCY].x = C.y;

  // Tongue tip must be right of the tongue body (once again).
  if (param[TTX].x < param[TCX].x + rx + tongueTipRadius)
  {
    param[TTX].x = param[TCX].x + rx + tongueTipRadius;
  }

  // Limit the tongue tip with respect to the lower border.
  // Therefore, set the anchor point to a very high position,
  // so that retractions happen essentially upwards.

  A.y = 1000.0;

  C = limitEllipsePos( Point2D(param[TTX].x, param[TTY].x), 
    tongueTipRadius, tongueTipRadius, lowerBorderTT, A);
  param[TTX].x = C.x;
  param[TTY].x = C.y;

  // Limit the tongue body with respect to the lower border.
  C = limitEllipsePos( Point2D(param[TCX].x, param[TCY].x), 
    anatomy.tongueCenterRadiusX_cm, anatomy.tongueCenterRadiusY_cm, lowerBorderTB, A);
  param[TCX].x = C.x;
  param[TCY].x = C.y;


  // ****************************************************************
  // Limit the middle point of the tongue root spline.
  // ****************************************************************

  restrictParam(TRX);
  restrictParam(TRY);

  r = anatomy.tongueCenterRadiusX_cm;
  if (anatomy.tongueCenterRadiusY_cm > r) 
  { 
    r = anatomy.tongueCenterRadiusY_cm; 
  }

  // C is the point in the bottom-left part of the circle, where
  // the tangent has the gradient -1 (see p. 30 of Birkholz' diss.).

  double alpha_rad = atan2(anatomy.tongueCenterRadiusY_cm, anatomy.tongueCenterRadiusX_cm);
  if ((alpha_rad <= M_PI) && (alpha_rad >= -M_PI))
  {
    alpha_rad+= M_PI;
  }
  C.x = param[TCX].x + cos(alpha_rad) * (anatomy.tongueCenterRadiusX_cm + EPSILON);
  C.y = param[TCY].x + sin(alpha_rad) * (anatomy.tongueCenterRadiusY_cm + EPSILON);

  // Hyoid

  Point2D H = surface[LOWER_COVER].getVertex(NUM_LARYNX_RIBS-1, NUM_LOWER_COVER_POINTS-1).toPoint2D();

  // Restrict the y-coordinate

  if (param[TRY].x < H.y) { param[TRY].x = H.y; }
  if (param[TRY].x > param[TCY].x) { param[TRY].x = param[TCY].x; }

  // Restrict the x-coordinate

  // Top-right slope
  if (param[TRX].x - C.x > C.y - param[TRY].x) 
  { 
    param[TRX].x = C.x + C.y - param[TRY].x; 
  }

  // Vertical left boundary.
  double minX = param[TCX].x - rx - (param[TCY].x - H.y);
  if (param[TRX].x < minX) 
  { 
    param[TRX].x = minX;
  }

  // Vertical right boundary.
  double maxX = param[TCX].x - r;
  if (H.x + 0.5 > maxX) { maxX = H.x + 0.5; }
  if (param[TRX].x > maxX) { param[TRX].x = maxX; }

  // ****************************************************************
  // Limit the middle point of the tongue blade spline.
  // ****************************************************************
  
  restrictParam(TBX);
  restrictParam(TBY);

  if (param[TBX].x < param[TCX].x) { param[TBX].x = param[TCX].x; }
  if (param[TBX].x > param[TTX].x) { param[TBX].x = param[TTX].x; }

  // Later, the ellipsoid shape should be considered here !!
  r = anatomy.tongueCenterRadiusX_cm;
  if (anatomy.tongueCenterRadiusY_cm > r) 
  { 
    r = anatomy.tongueCenterRadiusY_cm; 
  }

  r+= EPSILON;
  C.x = param[TCX].x;
  C.y = param[TCY].x + 1.415*r;
  if (param[TBY].x - C.y > param[TBX].x - C.x)  { param[TBY].x = C.y + param[TBX].x - C.x; }
  if (param[TBY].x - C.y < -(param[TBX].x-C.x)) { param[TBY].x = C.y -(param[TBX].x - C.x); }

  r = tongueTipRadius + EPSILON;
  C.x = param[TTX].x;
  C.y = param[TTY].x + 1.415*r;
  if (param[TBY].x - C.y < param[TBX].x - C.x)  { param[TBY].x = C.y + param[TBX].x - C.x; }

  // ****************************************************************
  // Limit the remaining tongue parameters.
  // ****************************************************************

  restrictParam(TS1);
  restrictParam(TS2);
  restrictParam(TS3);

  // ****************************************************************
  // Most of the original parameter values were limited above.
  // So, make these values the limited values.
  // ****************************************************************

  param[TCX].limitedX = param[TCX].x;
  param[TCY].limitedX = param[TCY].x;
  param[TTX].limitedX = param[TTX].x;
  param[TTY].limitedX = param[TTY].x;
  param[TBX].limitedX = param[TBX].x;
  param[TBY].limitedX = param[TBY].x;
  param[TRX].limitedX = param[TRX].x;
  param[TRY].limitedX = param[TRY].x;
  param[TS1].limitedX = param[TS1].x;
  param[TS2].limitedX = param[TS2].x;
  param[TS3].limitedX = param[TS3].x;

  // ****************************************************************
  // Only some parameters have different set and limited values.
  // ****************************************************************

  param[TTX].x = origTTX;
  param[TTY].x = origTTY;
  param[TCX].x = origTCX;
  param[TCY].x = origTCY;
}


// ****************************************************************************
/// Translates an ellipse with center point C and radii rx and ry on the line
/// towards the anchor point A until it is within the given border line.
/// The possibly translated center point of the ellipse is returned.
// ****************************************************************************

Point2D VocalTract::limitEllipsePos(Point2D C, double rx, double ry, LineStrip2D &border, Point2D A)
{
  const double EPSILON = 0.000001;
  int i;
  int N = border.getNumPoints();
  Point2D P0, P1, Q, w, u;
  double rx2 = rx*rx;
  double ry2 = ry*ry;
  double p, q;
  double root;
  double s, t;
  double alpha_rad;
  double denominator;

  // This is the Euclidian distance from the ellipse center to the
  // anchor point A, which will subsequently be reduced, if the
  // ellipse should be outside the border.

  double finalT = (C - A).magnitude();
  Point2D v = (C - A).normalize();

  // ****************************************************************
  // Run through all points (and line segments) of the border line).
  // ****************************************************************

  for (i=0; i < N; i++)
  {
    P0 = border.getControlPoint(i);
    Q = A - P0;

    // Check if the ellipse position should be moved towards A until
    // it just touches the border point P0.

    denominator = v.x*v.x*ry2 + v.y*v.y*rx2;
    if (fabs(denominator) > EPSILON)
    {
      p = 2.0*(Q.x*v.x*ry2 + Q.y*v.y*rx2) / denominator;
      q = (Q.x*Q.x*ry2 + Q.y*Q.y*rx2 - rx2*ry2) / denominator;
      root = 0.25*p*p - q;
      if (root >= 0.0)
      {
        root = sqrt(root);
        t = -0.5*p - root;
        if ((t >= 0.0) && (t < finalT))
        {
          finalT = t;
        }
      }
    }

    // Check if the ellipse position should be moved towards A until
    // it just touches the border line segment between P0 and P1.

    if (i < N-1)
    {
      P1 = border.getControlPoint(i + 1);
      w = P1 - P0;
      alpha_rad = atan2(-ry*w.x, rx*w.y);
      u.set( rx*cos(alpha_rad), ry*sin(alpha_rad) );

      denominator = v.x*w.y - v.y*w.x;
      if (fabs(denominator) > EPSILON)
      {
        // Assume the line segment is tangent at one side of the ellipse.
        Q = A - P0 + u;
        s = (Q.y*v.x - Q.x*v.y) / denominator;
        if ((s >= 0) && (s <= 1.0))
        {
          t = (Q.y*w.x - Q.x*w.y) / denominator;
          if ((t >= 0.0) && (t < finalT))
          {
            finalT = t;
          }
        }

        // Assume the line segment is tangent at the other side of the ellipse.
        Q = A - P0 - u;
        s = (Q.y*v.x - Q.x*v.y) / denominator;
        if ((s >= 0) && (s <= 1.0))
        {

          t = (Q.y*w.x - Q.x*w.y) / denominator;
          if ((t >= 0.0) && (t < finalT))
          {
            finalT = t;
          }
        }
      }

    }
  }

  C = A + v*finalT;

  return C;
}


// ****************************************************************************
// Calculate the shape of the individual tongue ribs.
// ****************************************************************************

void VocalTract::calcTongueRibs()
{
  const double INVALID = INVALID_PROFILE_SAMPLE;

  // 4.6 cm is the default width between the upper incisors
  double widthFactor = -2.0*anatomy.palatePoints[0].z / 4.6;    
  double HALF_WIDTH = 1.75*widthFactor;

  Tube::Articulator articulator;
  Surface *tongue = &surface[TONGUE];
  TongueRib *rib = &tongueRib[0];
  int i, k;
  Point3D P;
  double t;


  // ****************************************************************
  // Calc. the left and right target point for all ribs.
  // ****************************************************************

  for (i=0; i < NUM_TONGUE_RIBS; i++)
  {
    rib[i].left.set(-HALF_WIDTH, rib[i].leftSideHeight);
    rib[i].right.set(HALF_WIDTH, rib[i].rightSideHeight);
  }

  // ****************************************************************
  // The very first tongue rib corresponds to the hyoid rib of the
  // lower cover.
  // ****************************************************************

  int targetRib = NUM_LARYNX_RIBS-1;
  int targetPoint;
  Point3D Q;

  for (k=0; k < NUM_TONGUE_POINTS; k++)
  {
    targetPoint = 2 + k*(NUM_LOWER_COVER_POINTS*2-1-1-4) / (NUM_TONGUE_POINTS-1);
    if (targetPoint < NUM_LOWER_COVER_POINTS)
    {
      Q = surface[LOWER_COVER].getVertex(targetRib, targetPoint);
    }
    else
    {
      Q = surface[LOWER_COVER].getVertex(targetRib, 2*NUM_LOWER_COVER_POINTS-1-targetPoint);
      Q.z = -Q.z;
    }
    tongue->setVertex(0, k, Q);
  }
  rib[0].minX = tongue->getVertex(0, 0).z;
  rib[0].maxX = tongue->getVertex(0, tongue->numRibPoints-1).z;


  // ****************************************************************
  // Combine all ribs to a surface.
  // ****************************************************************

  const int N2 = NUM_PROFILE_SAMPLES/2;
  double z;
  double a2_left[NUM_TONGUE_RIBS], a3_left[NUM_TONGUE_RIBS];
  double a2_right[NUM_TONGUE_RIBS], a3_right[NUM_TONGUE_RIBS];
  Point2D leftEdge, rightEdge;
  double upperProfile[NUM_TONGUE_RIBS][NUM_PROFILE_SAMPLES];
  double lowerProfile[NUM_TONGUE_RIBS][NUM_PROFILE_SAMPLES];
  double x, y;
  bool lastInside, isInside;
  int index;
  double ratio;
  double minWidth;
  bool ok;
  
  for (i=1; i < NUM_TONGUE_RIBS; i++)
  {
    getCrossProfiles(rib[i].point, rib[i].normal, upperProfile[i], lowerProfile[i], false, articulator);

    // The left side (determine minX).
    
    leftEdge = rib[i].left;
    a2_left[i] = 3.0*leftEdge.y / (leftEdge.x*leftEdge.x);
    a3_left[i] = -2.0*leftEdge.y / (leftEdge.x*leftEdge.x*leftEdge.x);

    rib[i].minX = leftEdge.x;
    lastInside = false;
    ok = false;

    for (k=0; k < N2; k++)
    {
      x = -0.5*PROFILE_LENGTH + (double)k*PROFILE_SAMPLE_LENGTH;

      isInside = false;
      if ((upperProfile[i][k] != INVALID) && (lowerProfile[i][k] != INVALID) && (x >= leftEdge.x))
      {
        y = a2_left[i]*x*x + a3_left[i]*x*x*x;
        if ((y >= lowerProfile[i][k]) && (y <= upperProfile[i][k])) { isInside = true; }
      }

      if ((lastInside == false) && (isInside == true))
      {
        rib[i].minX = x;
        ok = true;
      }

      lastInside = isInside;
    }

    if (ok == false) { rib[i].minX = -0.1; }

    // The right side (determine maxX).

    rightEdge = rib[i].right;
    a2_right[i] = 3.0*rightEdge.y / (rightEdge.x*rightEdge.x);
    a3_right[i] = -2.0*rightEdge.y / (rightEdge.x*rightEdge.x*rightEdge.x);

    rib[i].maxX = rightEdge.x;
    lastInside = false;
    ok = false;

    for (k=NUM_PROFILE_SAMPLES; k > N2; k--)
    {
      x = (double)k*PROFILE_SAMPLE_LENGTH - 0.5*PROFILE_LENGTH;

      isInside = false;
      if ((upperProfile[i][k] != INVALID) && (lowerProfile[i][k] != INVALID) && (x <= rightEdge.x))
      {
        y = a2_right[i]*x*x + a3_right[i]*x*x*x;
        if ((y >= lowerProfile[i][k]) && (y <= upperProfile[i][k])) { isInside = true; }
      }

      if ((lastInside == false) && (isInside == true))
      {
        rib[i].maxX = x;
        ok = true;
      }

      lastInside = isInside;
    }

    if (ok == false) { rib[i].maxX = 0.1; }

    // **************************************************************
    // Every tongue rib has a minimum width. At the tongue root, it 
    // is somewhat smaller than in the rest of the vocal tract.
    // **************************************************************

    if (i < 5)
    {
      minWidth = 1.0 + 0.3*(double)i/5.0;
    }
    else
    if (i >= NUM_TONGUE_RIBS-10)
    {
      k = i - (NUM_TONGUE_RIBS-10);
      minWidth = 1.3 - 0.5*(double)k/10.0;
    }
    else
    {
      minWidth = 1.3;
    }
    minWidth*= widthFactor;     // In order to account for different VT depths

    if (rib[i].maxX <  minWidth) { rib[i].maxX =  minWidth; }
    if (rib[i].minX > -minWidth) { rib[i].minX = -minWidth; }
  }


  // ****************************************************************
  // Low-pass filter the width of the tongue with a recursive filter
  // applied forwards and backwards. The width may only get smaller, 
  // but never greater.
  // ****************************************************************

  const double LOWPASS_COEFF = 0.25;

  for (i=1; i < NUM_TONGUE_RIBS; i++)
  {
    z = rib[i-1].maxX + LOWPASS_COEFF*(rib[i].maxX - rib[i-1].maxX);
    if (z <= rib[i].maxX) { rib[i].maxX = z; }

    z = rib[i-1].minX + LOWPASS_COEFF*(rib[i].minX - rib[i-1].minX);
    if (z >= rib[i].minX) { rib[i].minX = z; }
  }

  for (i=NUM_TONGUE_RIBS-2; i >= 0; i--)
  {
    z = rib[i+1].maxX + LOWPASS_COEFF*(rib[i].maxX - rib[i+1].maxX);
    if (z <= rib[i].maxX) { rib[i].maxX = z; }

    z = rib[i+1].minX + LOWPASS_COEFF*(rib[i].minX - rib[i+1].minX);
    if (z >= rib[i].minX) { rib[i].minX = z; }
  }

  // ****************************************************************
  // Calc. the rib points and restrict them.
  // ****************************************************************

  for (i=1; i < NUM_TONGUE_RIBS; i++)
  {
    for (k=0; k < NUM_TONGUE_POINTS; k++)
    {
      z = rib[i].minX + (rib[i].maxX - rib[i].minX)*(double)k/(double)(NUM_TONGUE_POINTS-1);
      if (k <= NUM_TONGUE_POINTS/2)
      {
        t = a2_left[i]*z*z + a3_left[i]*z*z*z;
      }
      else
      {
        t = a2_right[i]*z*z + a3_right[i]*z*z*z;
      }

      // Test if the rib is within the hull.

      index = (int)((z + 0.5*PROFILE_LENGTH) / PROFILE_SAMPLE_LENGTH);
      ratio = (z + 0.5*PROFILE_LENGTH - (double)index*PROFILE_SAMPLE_LENGTH) / PROFILE_SAMPLE_LENGTH;
      if (index < 0) 
      { 
        index = 0; 
        ratio = 0.0;
      }
      if (index > NUM_PROFILE_SAMPLES-2) 
      { 
        index = NUM_PROFILE_SAMPLES-2;
        ratio = 1.0;
      }

      if ((upperProfile[i][index] != INVALID) && (upperProfile[i][index+1] != INVALID))
      {
        y = upperProfile[i][index] + ratio*(upperProfile[i][index+1]-upperProfile[i][index]) - 0.05;
        if (t > y) { t = y; }
      }

      tongue->setVertex(i, k, rib[i].point.x + rib[i].normal.x*t, rib[i].point.y + rib[i].normal.y*t, z);
    }
  }

}


// ****************************************************************************
/// The tongue rib "rigid" is tested for an intersection with the tongue rib
/// "flexible". When they intersect, the normal vector of flexible is adjusted
/// to prevent an intersection.
// ****************************************************************************

void VocalTract::verifyTongueRibNormal(int rigid, int flexible)
{
  double s, t;
  double denominator;
  Point2D P, v, R, w, Q;

  // The middle line is P+t*v with (min <= t <= max).
  
  R = tongueRib[rigid].point + tongueRib[rigid].minY*tongueRib[rigid].normal;
  w = tongueRib[rigid].point + tongueRib[rigid].maxY*tongueRib[rigid].normal - R;

  P = tongueRib[flexible].point;
  v = tongueRib[flexible].normal;

  Q = P - R;
  denominator = v.x*w.y - v.y*w.x;
  if (denominator != 0.0)
  {
    s = (v.x*Q.y - Q.x*v.y) / denominator;
    if ((s >= 0.0) && (s <= 1.0))
    {
      t = (w.x*Q.y - Q.x*w.y) / denominator;
      if ((t <= 0.0) && (t >= tongueRib[flexible].minY))
      {
        tongueRib[flexible].normal = tongueRib[flexible].point - R;
        tongueRib[flexible].normal.normalize();
      }
      else
      if ((t >= 0.0) && (t <= tongueRib[flexible].maxY))
      {
        tongueRib[flexible].normal = R+w - tongueRib[flexible].point;
        tongueRib[flexible].normal.normalize();
      }
    }
  }
}


// ****************************************************************************
/// Calculate the center line of the vocal tract.
// ****************************************************************************

void VocalTract::calcCenterLine()
{
  const int SKIPPED_LIP_POINTS = 1;
  int i, k;

  // Center and radii of the tongue body circle.

  Point3D tongueCenter(param[TCX].limitedX, param[TCY].limitedX, 0.0);
  double  tongueRadiusX = anatomy.tongueCenterRadiusX_cm;
  double  tongueRadiusY = anatomy.tongueCenterRadiusY_cm;

  // ****************************************************************
  // Calc. the length of the three sections of the mu-line.
  // ****************************************************************

  double muSectionLength[3];
  double muLength;
  Point2D R;

  double minY = upperOutline.getControlPoint(0).y;
  R = lowerOutline.getControlPoint(0);
  if (R.y > minY) { minY = R.y; }

  // The anterior termination plane is at one of the rib points of the
  // mid-sagittal lip rib - either the upper or lower lip, depending
  // which one is more posterior (= maxX).
  // For non-retracted lips, this corresponds here to a fixed position
  // between the lip corners and the anterior tangent to the upper and
  // lower lip analogous to Lindblom 2007: "On the acoustics of spread lips".

  double maxX = upperOutline.getControlPoint(upperOutline.getNumPoints() - 1 - SKIPPED_LIP_POINTS).x;
  R = lowerOutline.getControlPoint(lowerOutline.getNumPoints() - 1 - SKIPPED_LIP_POINTS);
  if (R.x < maxX) { maxX = R.x; }

  muSectionLength[0] = tongueCenter.y - minY;
  muSectionLength[1] = 0.5*M_PI*0.5*(tongueRadiusX + tongueRadiusY);
  muSectionLength[2] = maxX - tongueCenter.x;
  muLength = muSectionLength[0] + muSectionLength[1] + muSectionLength[2];

  // ****************************************************************
  // Intersect the hull contours of the vocal tract with the normals
  // of the mu-line at equal distances along the line.
  // ****************************************************************

  // Two points more for the subsequent smoothing
  Point2D nuLine[NUM_CENTERLINE_POINTS+2];    
  Point2D centerPoint;
  double pos;
  double angle;
  double sinus, cosinus;
  Point2D P, v;
  Point2D S[4];
  double t[4];
  bool ok[4];
  double tMax;
  int top, bottom;

  for (i=0; i < NUM_CENTERLINE_POINTS; i++)
  {
    pos = muLength*(double)i / (double)(NUM_CENTERLINE_POINTS-1);

    // Get the point on the mu-line and the corresponding normal vector.
    if (pos <= muSectionLength[0])
    {
      P.set(tongueCenter.x - tongueRadiusX, tongueCenter.y - muSectionLength[0] + pos);
      v.set(-1.0, 0.0);
    }
    else
    if (pos <= muSectionLength[0]+muSectionLength[1])
    {
      // Increase the tongue radius by 4 cm in order to always get the leftmost
      // intersection with the upper cover contour.
      angle = M_PI - 0.5*M_PI*(pos-muSectionLength[0])/muSectionLength[1];
      cosinus = cos(angle);
      sinus   = sin(angle);
      P.set(tongueCenter.x + tongueRadiusX*cosinus, tongueCenter.y + tongueRadiusY*sinus);
      v.set(tongueRadiusX*cosinus, tongueRadiusY*sinus);
      v.normalize();
    }
    else
    {
      P.set(tongueCenter.x + pos - muSectionLength[0] - muSectionLength[1], 
        tongueCenter.y + tongueRadiusY);
      v.set(0.0, 1.0);
    }

    // Calc. the intersection points with the three contour lines.

    ok[0] = upperOutline.getSpecialIntersection(P, v, t[0], S[0]);
    ok[1] = lowerOutline.getClosestIntersection(P, v, t[1], S[1]);
    ok[2] = tongueOutline.getClosestIntersection(P, v, t[2], S[2]);
    ok[3] = epiglottisOutline.getFirstIntersection(P, v, t[3], S[3]);

    centerPoint.set(10.0, 10.0);        // Prevent the error case.

    top = 0;
    bottom = 1;

    if (ok[0])      // The upper contour was intersected.
    {
      tMax = -1000000.0;
      if ((ok[1]) && (t[1] > tMax)) { tMax = t[1]; bottom = 1; }
      if ((ok[2]) && (t[2] > tMax)) { tMax = t[2]; bottom = 2; }
      if ((ok[3]) && (t[3] > tMax)) { tMax = t[3]; bottom = 3; }
    }
    
    centerPoint = 0.5*(S[top] + S[bottom]);
  
    // Fill the global variable for displaying the "rough center line".
    roughCenterLine[i].point = centerPoint;
    roughCenterLine[i].normal = v;
    roughCenterLine[i].max = (t[top] - t[bottom])/2.0;
    roughCenterLine[i].min = -roughCenterLine[i].max;
    roughCenterLine[i].pos = 0.0;

    // The nuLine is used subsequently.
    nuLine[i+1] = centerPoint;
  }

  // ****************************************************************
  // The nu-line is the smoothed mu-line.
  // ****************************************************************

  // Calc. the length and position on the rough center line.

  // Length of the nu-line (without the extra points at the beginning and end).
  double nuLineLength = 0.0;            
  double nu[NUM_CENTERLINE_POINTS+2];   // Two points more.
  Point2D diff;

  for (i=1; i < NUM_CENTERLINE_POINTS; i++)
  {
    nu[i] = nuLineLength;
    diff = nuLine[i+1] - nuLine[i];
    nuLineLength+= diff.magnitude();
  }
  nu[NUM_CENTERLINE_POINTS] = nuLineLength;

  // Variables for smoothing of the nu-line. One value for the 
  // beginning and the end of the smoothing window of each variable.
  
  Point2D rangePoint[2];
  int     rangeIndex[2];    // 0 <= index <= NUM_CENTERLINE_POINTS+1
  double  rangePos[2];
  double  length;
  double  d;

  const double RANGE_LENGTH = 2.0;
  const double HALF_RANGE_LENGTH = 0.5*RANGE_LENGTH;

  // Determine the first and last (helping) point on the nu-line.

  nuLine[0] = nuLine[1] + Point2D(0.0, -HALF_RANGE_LENGTH);
  nuLine[NUM_CENTERLINE_POINTS+1] = nuLine[NUM_CENTERLINE_POINTS] + Point2D(HALF_RANGE_LENGTH, 0.0);
  nu[0] = -HALF_RANGE_LENGTH;
  nu[NUM_CENTERLINE_POINTS+1] = nu[NUM_CENTERLINE_POINTS] + HALF_RANGE_LENGTH;

  // Initialize the range.

  rangeIndex[0] = rangeIndex[1] = 0;

  // ****************************************************************
  // Run through all points of the center line.
  // ****************************************************************

  for (i=0; i < NUM_CENTERLINE_POINTS; i++)
  {
    rangePos[0] = nuLineLength*(double)i / (double)(NUM_CENTERLINE_POINTS-1) - HALF_RANGE_LENGTH;
    rangePos[1] = rangePos[0] + RANGE_LENGTH;

    // Move the three range points forwards until they reached their
    // assigned positions.

    for (k=0; k < 2; k++)
    {
      while ((rangeIndex[k] < NUM_CENTERLINE_POINTS) && 
             (rangePos[k] > nu[rangeIndex[k] + 1])) { rangeIndex[k]++; }
      
      d = (rangePos[k] - nu[rangeIndex[k]]) / (nu[rangeIndex[k]+1] - nu[rangeIndex[k]]);
      rangePoint[k] = nuLine[rangeIndex[k]] + d*(nuLine[rangeIndex[k]+1] - nuLine[rangeIndex[k]]);
    }

    // Sum up all line segments in the current window.

    centerPoint.set(0.0, 0.0);

    if (rangeIndex[1] > rangeIndex[0])
    {
      // First segment
      length = (nuLine[rangeIndex[0] + 1] - rangePoint[0]).magnitude();
      centerPoint+= length*0.5*(rangePoint[0] + nuLine[rangeIndex[0] + 1]);

      // Middle segments
      for (k=rangeIndex[0]+1; k < rangeIndex[1]; k++)
      {
        length = (nuLine[k+1] - nuLine[k]).magnitude();
        centerPoint+= length*0.5*(nuLine[k] + nuLine[k+1]);
      }

      // Last segment
      length = (rangePoint[1] - nuLine[rangeIndex[1]]).magnitude();
      centerPoint+= length*0.5*(rangePoint[1] + nuLine[rangeIndex[1]]);
    }
    else
    {
      // The first and last point of the window are on the same line segment.
      length = (rangePoint[1] - rangePoint[0]).magnitude();
      centerPoint+= length*0.5*(rangePoint[0] + rangePoint[1]);
    }

    centerLine[i].point = centerPoint / RANGE_LENGTH;
  }

  // ****************************************************************
  // Change the positions of the first and last point of the center 
  // line such that they lie "on" the glottis and the lip plane.
  // ****************************************************************

  Point2D Q[2];
  double x, y;

  // At the glottis.

  Q[0] = upperOutline.getControlPoint(0);
  Q[1] = lowerOutline.getControlPoint(0);
  x = centerLine[0].point.x;
  y = Q[0].y + (Q[1].y-Q[0].y)*(x-Q[0].x) / (Q[1].x-Q[0].x);

  centerLine[0].point.y = y;
  for (i=0; i < 4; i++)
  {
    if (centerLine[i].point.y < y) { centerLine[i].point.y = y; }
  }

  // Set the normal vector.
  centerLine[0].normal = Q[0] - Q[1];
  centerLine[0].normal.normalize();

  // At the lips.

  centerLine[NUM_CENTERLINE_POINTS-1].normal.set(0.0, 1.0);
  x = centerLine[NUM_CENTERLINE_POINTS-1].point.x;
  for (i=NUM_CENTERLINE_POINTS-4; i < NUM_CENTERLINE_POINTS; i++)
  {
    if (centerLine[i].point.x > x) { centerLine[i].point.x = x; }
  }

  // ****************************************************************
  // Calc. the position of the center line points along the center line
  // and the whole length of the center line.
  // ****************************************************************

  centerLineLength = 0.0;
  for (i=0; i < NUM_CENTERLINE_POINTS-1; i++)
  {
    centerLine[i].pos = centerLineLength;
    d = (centerLine[i+1].point - centerLine[i].point).magnitude();
    centerLineLength+= d;
  }
  centerLine[NUM_CENTERLINE_POINTS-1].pos = centerLineLength;

  // Calc. the normal vectors.

  for (i=1; i < NUM_CENTERLINE_POINTS-1; i++)
  {
    Q[0] = centerLine[i].point - centerLine[i-1].point;
    Q[1] = centerLine[i+1].point - centerLine[i].point;

    centerLine[i].normal.set(-Q[0].y-Q[1].y, Q[0].x+Q[1].x);
    centerLine[i].normal.normalize();
  }

  // Calc. the min- and max-parameters for each cut vector.

  for (i=0; i < NUM_CENTERLINE_POINTS; i++)
  {
    centerLine[i].reserved = 0.0;   // The cutting line is ok.

    Point2D P = centerLine[i].point;
    Point2D n = centerLine[i].normal;

    ok[0] = upperOutline.getSpecialIntersection(P, n, t[0], S[0]);
    ok[1] = lowerOutline.getClosestIntersection(P, n, t[1], S[1]);
    ok[2] = tongueOutline.getClosestIntersection(P, n, t[2], S[2]);

    // Upper limit.

    if (ok[0])
    {
      centerLine[i].max = t[0];
    }
    else
    {
      centerLine[i].max = 3.0;        // Default value in the error case.
      centerLine[i].reserved = 1.0;   // Is supposed to be newly calculated later.
    }

    // Lower limit.

    if ((ok[1]) && (ok[2]))
    {
      if (t[1] > t[2]) 
      { 
        centerLine[i].min = t[1]; 
      } 
      else 
      { 
        centerLine[i].min = t[2]; 
      }
    }
    else
    if (ok[2])
    {
      centerLine[i].min = t[2];
    }
    else
    if (ok[1])
    {
      centerLine[i].min = t[1];
    }
    else
    {
      centerLine[i].min = -3.0;       // Default value in the error case.
      centerLine[i].reserved = 1.0;   // Is supposed to be newly calculated later.
    }
  }

  // ****************************************************************
  // Run through all normal vectors in a certain order and adjust 
  // their directions such that they don't overlap.
  // ****************************************************************

  int numToCheck;
  int distance;
  int index;

  for (i=0; i < NUM_CENTERLINE_POINTS_EXPONENT; i++)
  {
    numToCheck = 1 << i;
    distance = 1 << (NUM_CENTERLINE_POINTS_EXPONENT-i);
    index = distance/2;

    for (k=0; k < numToCheck; k++)
    {
      verifyCenterLineNormal(index-distance/2, index, index+distance/2);
      index+= distance;
    }
  }

  // ****************************************************************
  // For the normal vectors that have changed, recalculate the min-
  // and max-values.
  // ****************************************************************

  Point2D n;

  for (i=0; i < NUM_CENTERLINE_POINTS; i++)
  {
    if (centerLine[i].reserved != 0.0)
    {
      P = centerLine[i].point;
      n = centerLine[i].normal;
    
      ok[0] = upperOutline.getSpecialIntersection(P, n, t[0], S[0]);
      ok[1] = lowerOutline.getClosestIntersection(P, n, t[1], S[1]);
      ok[2] = tongueOutline.getClosestIntersection(P, n, t[2], S[2]);

      // Upper border.
      if (ok[0])
      {
        centerLine[i].max = t[0];
      }
      else
      {
        centerLine[i].max = 3.0;        // Default-Wert im Fehlerfall
      }

      // Lower border.
      if ((ok[1]) && (ok[2]))
      {
        if (t[1] > t[2]) { centerLine[i].min = t[1]; } else { centerLine[i].min = t[2]; }
      }
      else
      if (ok[2])
      {
        centerLine[i].min = t[2];
      }
      else
      if (ok[1])
      {
        centerLine[i].min = t[1];
      }
      else
      {
        centerLine[i].min = -3.0;       // Default-Wert im Fehlerfall
      }
    }
  }

}

// ****************************************************************************
/// If the normal vector of the center line intersects its left or right 
/// neightbor, its direction is adjusted such that the intersection is just
/// prevented.
// ****************************************************************************

void VocalTract::verifyCenterLineNormal(int l, int m, int r)
{
  double s, t;
  double denominator;
  Point2D P, v, R, w, Q;

  // The middle line is P+t*v with (min <= t <= max)
  
  P = centerLine[m].point;
  v = centerLine[m].normal;

  // ****************************************************************
  // Intersection with the left line ?
  // ****************************************************************

  R = centerLine[l].point + centerLine[l].min*centerLine[l].normal;
  w = centerLine[l].point + centerLine[l].max*centerLine[l].normal - R;

  Q = P - R;
  denominator = v.x*w.y - v.y*w.x;
  if (denominator != 0.0)
  {
    s = (v.x*Q.y - Q.x*v.y) / denominator;
    if ((s >= 0.0) && (s <= 1.0))
    {
      t = (w.x*Q.y - Q.x*w.y) / denominator;
      if ((t <= 0.0) && (t >= centerLine[m].min))
      {
        centerLine[m].normal = centerLine[m].point - R;
        centerLine[m].normal.normalize();
        centerLine[m].reserved = 1.0;     // Because it was changed.
      }
      else
      if ((t >= 0.0) && (t <= centerLine[m].max))
      {
        centerLine[m].normal = R+w - centerLine[m].point;
        centerLine[m].normal.normalize();
        centerLine[m].reserved = 1.0;     // Because it was changed.
      }
    }
  }

  // ****************************************************************
  // Intersection with the right line ?
  // ****************************************************************

  R = centerLine[r].point + centerLine[r].min*centerLine[r].normal;
  w = centerLine[r].point + centerLine[r].max*centerLine[r].normal - R;

  Q = P - R;
  denominator = v.x*w.y - v.y*w.x;
  if (denominator != 0.0)
  {
    s = (v.x*Q.y - Q.x*v.y) / denominator;
    if ((s >= 0.0) && (s <= 1.0))
    {
      t = (w.x*Q.y - Q.x*w.y) / denominator;
      if ((t <= 0.0) && (t >= centerLine[m].min))
      {
        centerLine[m].normal = centerLine[m].point - R;
        centerLine[m].normal.normalize();
        centerLine[m].reserved = 1.0;     // Because it was changed.
      }
      else
      if ((t >= 0.0) && (t <= centerLine[m].max))
      {
        centerLine[m].normal = R+w - centerLine[m].point;
        centerLine[m].normal.normalize();
        centerLine[m].reserved = 1.0;     // Because it was changed.
      }
    }
  }
}


// ****************************************************************************
/// Calculates for each cut vector on the center line the cross-sectional
/// profile.
// ****************************************************************************

void VocalTract::calcCrossSections()
{
  double upperProfile[NUM_PROFILE_SAMPLES];
  double lowerProfile[NUM_PROFILE_SAMPLES];
  double tongueTipRadius = anatomy.tongueTipRadius_cm;
  Tube::Articulator articulator;
  int i;

  // ****************************************************************

  for (i=0; i < NUM_CENTERLINE_POINTS; i++)
  {
    getCrossProfiles(centerLine[i].point, centerLine[i].normal, upperProfile, lowerProfile, true, articulator);
    getCrossSection(upperProfile, lowerProfile, &crossSection[i]);

    crossSection[i].pos = centerLine[i].pos;
    crossSection[i].articulator = articulator;
  }

  // ****************************************************************
  // Calculate data about the nasal port.
  // ****************************************************************

  Surface *s = &surface[UPPER_COVER];
  int portRib = NUM_LARYNX_RIBS + NUM_PHARYNX_RIBS;

  // Area:
  nasalPortArea_cm2 = anatomy.maxNasalPortArea_cm2*param[VO].x;
  // Min. value for VO < 0. Therefore, limit the port area here.
  if (nasalPortArea_cm2 < 0.0) { nasalPortArea_cm2 = 0.0; }

  // Position:
  int bestIndex;
  double bestT;
  nasalPortPos_cm = getCenterLinePos(s->getVertex(portRib, s->numRibPoints/2).toPoint2D(), bestIndex, bestT);

  // ****************************************************************
  // Position of the incisors on the center line.
  // ****************************************************************

  const double EPSILON = 0.000001;
  Point2D P, Q, R, w;
  double teethX = surface[UPPER_TEETH].getVertex(NUM_TEETH_RIBS-1, 2).x;
  double t;

  for (i=0; i < NUM_CENTERLINE_POINTS-1; i++)
  {
    P = centerLine[i].point;
    Q = centerLine[i+1].point;
    if ((P.x < teethX) && (Q.x >= teethX))
    {
      t = Q.x - P.x;
      if (t < EPSILON) { t = EPSILON; }
      incisorPos_cm = centerLine[i].pos + ((teethX - P.x) / t)*(Q-P).magnitude();
    }
  }

  // ****************************************************************
  // For cross-sections that are right from the anterior side of the
  // tongue tip, the associated articulator may not be the tongue,
  // because, if the tongue was cut, it is only the mouth floor.
  // ****************************************************************

  for (i=0; i < NUM_CENTERLINE_POINTS; i++)
  {
    P = centerLine[i].point;
    w = centerLine[i].normal;
    // Q is the center of the tongue tip circle.
    Q.x = param[TTX].limitedX;
    Q.y = param[TTY].limitedX;

    R = Q - P;
    
    // Is the cross-sectional line more than the tongue tip radius
    // right from the tongue tip center point ?

    if ((R.y*w.x - R.x*w.y > tongueTipRadius) && (P.x > param[TCX].limitedX))
    {
      // If so, intersections with the tongue (e.g., at the mouth
      // floor) do not count as tongue any more.
      if (crossSection[i].articulator == Tube::TONGUE)
      {
        crossSection[i].articulator = Tube::OTHER_ARTICULATOR;
      }
    }
  }

  // ****************************************************************
  // Make sure that the minimal areas in different regions are not 
  // undershot.
  // ****************************************************************

  // Directly around the upper incisors, the cross-sectional area
  // is never below 15 mm^2, because of the gaps between the teeth.
  // This is required for /f,v/ and partly for "th" in English.

  const double LEFT_INCISOR_MARGIN_CM = 0.5;   // 5 mm left from the incisor pos.
  const double RIGHT_INCISOR_MARGIN_CM = 0.3;   // 3 mm right from the incisor pos.
  const double MIN_INCISOR_AREA_CM2 = 0.15;

  for (i = 0; i < NUM_CENTERLINE_POINTS; i++)
  {
    if ((crossSection[i].pos >= incisorPos_cm - LEFT_INCISOR_MARGIN_CM) &&
      (crossSection[i].pos <= incisorPos_cm + RIGHT_INCISOR_MARGIN_CM) &&
      (crossSection[i].area < MIN_INCISOR_AREA_CM2))
    {
      crossSection[i].area = MIN_INCISOR_AREA_CM2;
      crossSection[i].circ = 2.0 * sqrt(MIN_INCISOR_AREA_CM2 * M_PI);
    }
  }

  // ****************************************************************
  // Process the regions of the tongue back an the tongue tip.
  // ****************************************************************

  double minArea[2] = { 0.0 };
  double minCirc[2] = { 0.0 };
  
  minArea[0] = tongueSideParamToMinArea_cm2(param[TS2].x);    // For the tongue back region
  // Lateral passages in the tongue back region are not possible.
  if (param[TS2].x < 0.0)
  {
    minArea[0] = 0.0;
  }
  minArea[1] = tongueSideParamToMinArea_cm2(param[TS3].x);    // For the tongue tip region

  minCirc[0] = 2.0 * sqrt(minArea[0] * M_PI);
  minCirc[1] = 2.0 * sqrt(minArea[1] * M_PI);

  // Where along the center line start and end the tongue tip (circle)
  // and the lower lip?

  double pos_cm = 0.0;
  double tongueTipRight_cm = 0.0;
  double lowerLipLeft_cm = 1000000.0;

  for (i = 0; i < NUM_CENTERLINE_POINTS; i++)
  {
    pos_cm = crossSection[i].pos;
    if (crossSection[i].articulator == Tube::TONGUE)
    {
      tongueTipRight_cm = pos_cm;
    }
    
    if ((crossSection[i].articulator == Tube::LOWER_LIP) && (pos_cm < lowerLipLeft_cm))
    {
      lowerLipLeft_cm = pos_cm;
    }
  }

  double tongueTipLeft_cm = tongueTipRight_cm - 2.0;

  // ****************************************************************

  CrossSection* cs = NULL;

  for (i = 0; i < NUM_CENTERLINE_POINTS; i++)
  {
    cs = &crossSection[i];

    if (cs->pos <= tongueTipLeft_cm)
    {
      if (cs->area < minArea[0]) { cs->area = minArea[0]; }
      if (cs->circ < minCirc[0]) { cs->circ = minCirc[0]; }
    }

    // The "real" tongue tip often extends somewhat further towards 
    // the right. Therefore take the left end of the lower lips as 
    // the right boundary of the tongue tip region.

    if ((cs->pos >= tongueTipLeft_cm) && (cs->pos <= lowerLipLeft_cm))
    {
      if (cs->area < minArea[1]) { cs->area = minArea[1]; }
      if (cs->circ < minCirc[1]) { cs->circ = minCirc[1]; }
    }
  }

}


// ****************************************************************************
/// Calculates the upper and lower profile of a cross-section defined by the
/// Point P and normal vector v on the center line.
// ****************************************************************************

void VocalTract::getCrossProfiles(Point2D P, Point2D v, double *upperProfile, 
       double *lowerProfile, bool considerTongue, Tube::Articulator &articulator, bool debug)
{
  const double MIN_SQUARED_NORMAL_LENGTH = 0.0000001;
  const double INVALID = INVALID_PROFILE_SAMPLE;
  const int N = NUM_PROFILE_SAMPLES;
  const int N2 = NUM_PROFILE_SAMPLES/2;
  const int TOP = 1;       // = Bit at pos. 0
  const int BOTTOM = 2;    // = Bit at pos. 1
  const int NUM_PROFILE_SURFACES = 10;

  int i, k;
  Surface *s;
  Point2D Q, P0, P1, n;
  Point2D lowestTeethPoint;
  int left, right;
  int localIndex, globalIndex;
  bool rightOrientation;

  // Global surface indices for all pixels in the upper and lower
  // profiles, i.e., which model surfaces made the points in the
  // profiles?

  int upperProfileSurface[NUM_PROFILE_SAMPLES];
  int lowerProfileSurface[NUM_PROFILE_SAMPLES];

  // Handle all upper-posterior surfaces first, and then the other
  // surfaces.
  int profileSurfaceIndex[NUM_PROFILE_SURFACES] = 
  {
    // Surfaces that contribute to the upper profile
    UPPER_COVER, 
    UPPER_TEETH, 
    UPPER_LIP, 
    UVULA, 
    // Surfaces that contribute to the lower profile
    LOWER_COVER,
    LOWER_TEETH, 
    LOWER_LIP, 
    EPIGLOTTIS, 
    // Surfaces that may contribute to the upper and lower profile
    LEFT_COVER, 
    RADIATION
  };

  const int MAX_CUTS = 2048;
  
  struct Cut
  {
    Point2D P0, P1, n;
    int globalSurfaceIndex;
    int localSurfaceIndex;
  };
  
  Cut cut[MAX_CUTS];
  int numCuts;

  // For the fast intersection method:

  const int MAX_LIST_ENTRIES = 1024;
  int indexList[MAX_LIST_ENTRIES];
  int numListEntries;

  // ****************************************************************
  // Initialize all profile values.
  // ****************************************************************

  for (i=0; i < NUM_PROFILE_SAMPLES; i++) 
  { 
    upperProfile[i] = EXTREME_PROFILE_VALUE;
    lowerProfile[i] = -EXTREME_PROFILE_VALUE;
    upperProfileSurface[i] = -1;
    lowerProfileSurface[i] = -1;
  }

  // ****************************************************************
  // Cut all surfaces except the tongue and, if considerTongue == false,
  // except the uvula and epiglottis. The uvuala and epiglottis sur-
  // faces are formed after the tongue, so they are not valid yet.
  // ****************************************************************

  numCuts = 0;
  for (k=0; k < NUM_PROFILE_SURFACES; k++)
  {
    globalIndex = profileSurfaceIndex[k];

    if ((considerTongue) || ((considerTongue == false) && 
        (globalIndex != UVULA) && (globalIndex != EPIGLOTTIS)))
    {
      s = &surface[globalIndex];

      // The fast intersection method.

      if (makeFasterIntersections) 
      {
        if (intersectionsPrepared[globalIndex] == false)
        {
          // Assign all triangles to the tiles (only once !)
          s->prepareIntersections();
          intersectionsPrepared[globalIndex] = true;
        }

        s->prepareIntersection(P, v);
        s->getTriangleList(indexList, numListEntries, MAX_LIST_ENTRIES);

        for (i=0; i < numListEntries; i++)
        {
          if ((s->getTriangleIntersection(indexList[i], P0, P1, n)) && (numCuts < MAX_CUTS) &&
              (P0.y < MAX_PROFILE_VALUE) && (P1.y < MAX_PROFILE_VALUE) &&
              (P1.y > MIN_PROFILE_VALUE) && (P1.y > MIN_PROFILE_VALUE))
          {
            cut[numCuts].P0 = P0;
            cut[numCuts].P1 = P1;
            cut[numCuts].n = n;
            cut[numCuts].globalSurfaceIndex = globalIndex;
            cut[numCuts].localSurfaceIndex = k;
            numCuts++;
          }
        }
      }
      else

      // The "normal", slower intersection method.

      {
        s->prepareIntersection(P, v);

        for (i=0; i < s->numTriangles; i++)
        {
          if ((s->getTriangleIntersection(i, P0, P1, n)) && (numCuts < MAX_CUTS) &&
              (P0.y < MAX_PROFILE_VALUE) && (P1.y < MAX_PROFILE_VALUE) &&
              (P1.y > MIN_PROFILE_VALUE) && (P1.y > MIN_PROFILE_VALUE))
          {
            cut[numCuts].P0 = P0;
            cut[numCuts].P1 = P1;
            cut[numCuts].n = n;
            cut[numCuts].globalSurfaceIndex = globalIndex;
            cut[numCuts].localSurfaceIndex = k;
            numCuts++;
          }
        }
      }

    }
  }

  // ****************************************************************
  // Dertermine how the teeth and lips were cut!
  // ****************************************************************

  double upperTeethMinY = EXTREME_PROFILE_VALUE;
  double lowerTeethMaxY = -EXTREME_PROFILE_VALUE;
  double upperTeethMaxY = -EXTREME_PROFILE_VALUE;
  double lowerTeethMinY = EXTREME_PROFILE_VALUE;
  double upperTeethMinX = 0.0;
  double lowerTeethMinX = 0.0;
  double upperLipMinY = EXTREME_PROFILE_VALUE;
  double lowerLipMaxY = -EXTREME_PROFILE_VALUE;
  bool bothTeethCut;        // Were both lips cut ?
  bool bothLipsCut;         // Were both teeth cut ?

  for (i=0; i < numCuts; i++)
  {
    P0 = cut[i].P0;
    P1 = cut[i].P1;
    globalIndex = cut[i].globalSurfaceIndex;

    // Teeth data
    if (globalIndex == UPPER_TEETH)
    {
      if (P0.y < upperTeethMinY) { upperTeethMinY = P0.y; }
      if (P1.y < upperTeethMinY) { upperTeethMinY = P1.y; }
      if (P0.y > upperTeethMaxY) { upperTeethMaxY = P0.y; }
      if (P1.y > upperTeethMaxY) { upperTeethMaxY = P1.y; }
      if (P0.x < upperTeethMinX) { upperTeethMinX = P0.x; }
      if (P1.x < upperTeethMinX) { upperTeethMinX = P1.x; }
    }
    if (globalIndex == LOWER_TEETH)
    {
      if (P0.y > lowerTeethMaxY) { lowerTeethMaxY = P0.y; }
      if (P1.y > lowerTeethMaxY) { lowerTeethMaxY = P1.y; }
      if (P0.y < lowerTeethMinY) { lowerTeethMinY = P0.y; }
      if (P1.y < lowerTeethMinY) { lowerTeethMinY = P1.y; }
      if (P0.x < lowerTeethMinX) { lowerTeethMinX = P0.x; }
      if (P1.x < lowerTeethMinX) { lowerTeethMinX = P1.x; }
    }

    // Lip data
    if (globalIndex == UPPER_LIP)
    {
      if (P0.y < upperLipMinY) { upperLipMinY = P0.y; }
      if (P1.y < upperLipMinY) { upperLipMinY = P1.y; }
    }
    if (globalIndex == LOWER_LIP)
    {
      if (P0.y > lowerLipMaxY) { lowerLipMaxY = P0.y; }
      if (P1.y > lowerLipMaxY) { lowerLipMaxY = P1.y; }
    }
  }

  // Where both teeth or teeth cut in the right orientation ?

  if (upperTeethMaxY > lowerTeethMinY) 
  { 
    bothTeethCut = true; 
  } 
  else 
  { 
    bothTeethCut = false; 
  }
  
  if ((upperLipMinY != EXTREME_PROFILE_VALUE) && (lowerLipMaxY != -EXTREME_PROFILE_VALUE)) 
  { 
    bothLipsCut = true; 
  } 
  else 
  { 
    bothLipsCut = false; 
  }

  // ****************************************************************
  // Create the upper and lower profile.
  // ****************************************************************

  for (i=0; i < numCuts; i++)
  {
    P0 = cut[i].P0;
    P1 = cut[i].P1;
    n = cut[i].n;
    localIndex = cut[i].localSurfaceIndex;
    globalIndex = cut[i].globalSurfaceIndex;

    // **************************************************************
    // Surfaces that contribute to the upper profile.
    // **************************************************************

    if (globalIndex == UPPER_COVER)
    {
      if (n.y < 0.0) { insertUpperProfileLine(P0, P1, globalIndex, upperProfile, upperProfileSurface); }
    }
    else
    
    if (globalIndex == UPPER_TEETH)
    {
      if (n.y < 0.0) { insertUpperProfileLine(P0, P1, globalIndex, upperProfile, upperProfileSurface); }
    }
    else

    if (globalIndex == UPPER_LIP)
    {
      if (n.y < 0.0) { insertUpperProfileLine(P0, P1, globalIndex, upperProfile, upperProfileSurface); }
    }
    else

    if (globalIndex == UVULA)
    {
      if (n.y < 0.0) { insertUpperProfileLine(P0, P1, globalIndex, upperProfile, upperProfileSurface); }
    }
    else

    // **************************************************************
    // Surfaces that contribute to the lower profile.
    // **************************************************************

    if (globalIndex == LOWER_COVER)
    {
      // The lower cover is special - it may also supplement
      // the upper profile.
      if (n.y > 0.0) 
      { 
        insertLowerCoverProfileLine(P0, P1, globalIndex, upperProfile, 
          upperProfileSurface, lowerProfile, lowerProfileSurface); 
      }
    }
    else

    if (globalIndex == LOWER_TEETH)
    {
      if (n.y > 0.0) { insertLowerProfileLine(P0, P1, globalIndex, lowerProfile, lowerProfileSurface); }
    }
    else

    if (globalIndex == LOWER_LIP)
    {
      if (n.y > 0.0) { insertLowerProfileLine(P0, P1, globalIndex, lowerProfile, lowerProfileSurface); }
    }
    else

    if (globalIndex == EPIGLOTTIS)
    {
      if (n.y > 0.0) { insertLowerProfileLine(P0, P1, globalIndex, lowerProfile, lowerProfileSurface); }
    }
    else

    // **************************************************************
    // All other surfaces (LEFT_COVER and RADIATION) can contribute 
    // to either the upper or the lower profile, depending on the 
    // surface normal.
    // And only if the surface normal is long enough (> epsilon).
    // **************************************************************

    {
      rightOrientation = true;

      // Consider these surfaces only in the oral and upper pharyngeal
      // part of the vocal tract, because it may cause sporadic errors 
      // in the lower pharyngeal part.
      // Therefore, the normal vector v of the cutting line must not
      // point into the bottom left (posterior) quadrant.

      if ((v.x < 0.0) && (v.y < 0.0))
      {
        rightOrientation = false;
      }

      if ((n.x*n.x + n.y*n.y > MIN_SQUARED_NORMAL_LENGTH) && (rightOrientation))
      {
        if (n.y <= 0.0)
        {
          insertUpperProfileLine(P0, P1, globalIndex, upperProfile, upperProfileSurface);
        }
        else
        if (n.y >= 0.0)
        {
          insertLowerProfileLine(P0, P1, globalIndex, lowerProfile, lowerProfileSurface);
        }
      }
    }
  }

  // ****************************************************************
  // Mark invalid profile values as such.
  // ****************************************************************

  for (i=0; i < N2; i++)
  {
    if (upperProfile[i] == EXTREME_PROFILE_VALUE)  { upperProfile[i] = INVALID; }
    if (lowerProfile[i] == -EXTREME_PROFILE_VALUE) { lowerProfile[i] = INVALID; }
  }

  // ****************************************************************
  // When the upper teeth are cut, the profile more lateral from the
  // teeth may not be higher than the lower edge of the upper teeth.
  // Analog for the lower profile.
  // This is to make sure that in the vicinity of the incisors, the
  // area between the lips left and right from the teeth is not
  // considered in the profile (only with the profile height between 
  // the teeth). This ensures a constantly smaller area between the
  // incisors, which is important for transitions from and to /sh/.
  // ****************************************************************


  double upperLimit = 1000000.0;
  double lowerLimit = -1000000.0;

  // Go from median to lateral.
  for (i = N2 - 1; i >= 0; i--)
  {
    if (upperProfile[i] != INVALID)
    {
      if (upperProfileSurface[i] == UPPER_TEETH)
      {
        if (upperProfile[i] < upperLimit)
        {
          upperLimit = upperProfile[i];
        }
      }

      if (upperProfile[i] > upperLimit)
      {
        upperProfile[i] = upperLimit;
      }
    }

    if (lowerProfile[i] != INVALID)
    {
      if (lowerProfileSurface[i] == LOWER_TEETH)
      {
        if (lowerProfile[i] > lowerLimit)
        {
          lowerLimit = lowerProfile[i];
        }
      }

      if (lowerProfile[i] < lowerLimit)
      {
        lowerProfile[i] = lowerLimit;
      }
    }
  }


  // ****************************************************************
  // Fill gaps in the upper profile.
  // ****************************************************************
  
  double lastValue = INVALID;

  for (i=0; i < N2; i++)
  {
    if (upperProfile[i] == INVALID)
    {
      upperProfile[i] = lastValue;
    }
    else
    {
      lastValue = upperProfile[i];
    }
  }  

  // ****************************************************************
  // All invalid samples of the lower profile are filled up with 
  // minimum values starting from the middle, until the first valid
  // value is found.
  // If all samples are invalid, assume the width of the upper profile.
  // ****************************************************************

  int firstUpperValid = -1;
  int firstLowerValid = -1;
  for (i=0; i < N2; i++)
  {
    if ((firstLowerValid == -1) && (lowerProfile[i] != INVALID)) { firstLowerValid = i; }
    if ((firstUpperValid == -1) && (upperProfile[i] != INVALID)) { firstUpperValid = i; }
  }

  if (firstLowerValid == -1)
  {
    if (firstUpperValid == -1) { firstUpperValid = 0; }
    for (i=firstUpperValid; i < N2; i++) { lowerProfile[i] = MIN_PROFILE_VALUE; }
  }
  else
  {
    i = N2 - 1;
    while ((i > 0) && (lowerProfile[i] == INVALID)) { lowerProfile[i--] = MIN_PROFILE_VALUE; }
  }

  // ****************************************************************
  // Interpolate gaps in the lower profile.
  // ****************************************************************

  left = 0;
  right = 0;

  for (i=0; i < N2; i++)
  {
    if (lowerProfile[i] == INVALID)
    {
      if ((i <= left) || (i >= right))
      {
        // Find a new valid left and right sample.
        left = i;
        right = i;
        while ((lowerProfile[left] == INVALID)  && (left > 0))     { left--; }
        while ((lowerProfile[right] == INVALID) && (right < N2-1)) { right++; }
      }

      if ((i > left) && (i < right))
      {
        if ((lowerProfile[left] != INVALID) && (lowerProfile[right] != INVALID))
        {
          lowerProfile[i] = lowerProfile[left] + (double)(i-left)*(lowerProfile[right]-lowerProfile[left])/(right-left);
        }
        else
        if (lowerProfile[left] != INVALID)
        {
          lowerProfile[i] = lowerProfile[left];
        }
      }
    }
  }

  // ****************************************************************
  // The lower profile samples must always lie below the upper 
  // profile samples.
  // ****************************************************************

  for (i=0; i < N2; i++)
  {
    if ((lowerProfile[i] != INVALID) && (upperProfile[i] != INVALID) && (lowerProfile[i] > upperProfile[i]))
    {
      lowerProfile[i] = upperProfile[i];
    }
  }

  // ****************************************************************
  // Find the leftmost valid sample in both profiles.
  // ****************************************************************

  int upperLeft = 0;
  while ((upperProfile[upperLeft] == INVALID) && (upperLeft < N2-1)) { upperLeft++; }

  int lowerLeft = 0;
  while ((lowerProfile[lowerLeft] == INVALID) && (lowerLeft < N2-1)) { lowerLeft++; }

  // ****************************************************************
  // The leftmost valid sample index.
  // ****************************************************************

  int leftmost = N2-1;
  if ((upperProfile[upperLeft] != INVALID) && (upperLeft < leftmost)) { leftmost = upperLeft; }
  if ((lowerProfile[lowerLeft] != INVALID) && (lowerLeft < leftmost)) { leftmost = lowerLeft; }

  // ****************************************************************
  // Make both profiles have an equal width.
  // ****************************************************************

  if (upperProfile[upperLeft] == INVALID) { upperProfile[upperLeft] = MAX_PROFILE_VALUE; }
  if (lowerProfile[lowerLeft] == INVALID) { lowerProfile[lowerLeft] = MIN_PROFILE_VALUE; }

  if (upperLeft < lowerLeft)
  {
    for (i=upperLeft; i < lowerLeft; i++) { lowerProfile[i] = lowerProfile[lowerLeft]; }
  }
  else
  {
    for (i=lowerLeft; i < upperLeft; i++) { upperProfile[i] = upperProfile[upperLeft]; }
  }

  // ****************************************************************
  // When both lips were cut, and they are closed in the midsagittal
  // plane, then the profile is closed completely!
  // ****************************************************************

  if ((bothLipsCut) && (bothTeethCut == false) && (upperProfile[N2-1] == lowerProfile[N2-1]))
  {
    double d;
    for (i=0; i < N2; i++)
    {
      if ((upperProfile[i] != INVALID) && (lowerProfile[i] != INVALID))
      {
        d = 0.5*(lowerProfile[i] + upperProfile[i]);
        lowerProfile[i] = d;
        upperProfile[i] = d;
      }
    }
  }

  // ****************************************************************
  // Mirror the right half of the profiles to the left.
  // ****************************************************************

  for (i=0; i <= N2; i++)
  {
    upperProfile[N-1-i] = upperProfile[i];
    lowerProfile[N-1-i] = lowerProfile[i];
  }
  int rightmost = N-1-leftmost;


  // ****************************************************************
  // Consider the tongue.
  // ****************************************************************

  if (considerTongue)
  {
    double tongueProfile[NUM_PROFILE_SAMPLES];
    int tongueProfileSurface[NUM_PROFILE_SAMPLES];

    // Init. the profile.

    for (i=0; i < N; i++) 
    { 
      tongueProfile[i] = -EXTREME_PROFILE_VALUE; 
      tongueProfileSurface[i] = -1;
    }

    s = &surface[TONGUE];

    // Faster intersection method.

    if (makeFasterIntersections)
    {
      if (intersectionsPrepared[TONGUE] == false)
      {
        // Assign all triangles to the tiles (only once !)
        s->prepareIntersections();
        intersectionsPrepared[TONGUE] = true;
      }

      s->prepareIntersection(P, v);
      s->getTriangleList(indexList, numListEntries, MAX_LIST_ENTRIES);

      for (i=0; i < numListEntries; i++)
      {
        if (s->getTriangleIntersection(indexList[i], P0, P1, n))
        {
          if (n.y >= 0.0) 
          { 
            insertLowerProfileLine(P0, P1, TONGUE, tongueProfile, tongueProfileSurface); 
          }
        }
      }
    }
    else

    // Slower intersection method.

    {
      s->prepareIntersection(P, v);
      for (i=0; i < s->numTriangles; i++)
      {
        if (s->getTriangleIntersection(i, P0, P1, n))
        {
          if (n.y >= 0.0) 
          { 
            insertLowerProfileLine(P0, P1, TONGUE, tongueProfile, tongueProfileSurface); 
          }
        }
      }
    }

    for (i=0; i < N; i++)
    {
      if (tongueProfile[i] == -EXTREME_PROFILE_VALUE) { tongueProfile[i] = INVALID; }
    }

    // **************************************************************
    // Linear interpolation of gaps in the tongue profile.
    // **************************************************************

    left = 0;
    right = 0;

    for (i=0; i < N; i++)
    {
      if (tongueProfile[i] == INVALID)
      {
        if ((i <= left) || (i >= right))
        {
          // Find a new valid left and right sample.
          left = i;
          right = i;
          while ((tongueProfile[left] == INVALID)  && (left > 0))    { left--; }
          while ((tongueProfile[right] == INVALID) && (right < N-1)) { right++; }
        }

        if (right > left)
        {
          if ((tongueProfile[left] != INVALID) && (tongueProfile[right] != INVALID))
          {
            tongueProfile[i] = tongueProfile[left] + (double)(i-left)*(tongueProfile[right]-tongueProfile[left])/(right-left);
          }
          else
          if (tongueProfile[left] != INVALID)
          {
            tongueProfile[i] = tongueProfile[left];
          }
          else
          if (tongueProfile[right] != INVALID)
          {
            tongueProfile[i] = tongueProfile[right];
          }
        }
      }
    }

    // **************************************************************
    // Merge the tongue with the lower profile.
    // **************************************************************

    for (i=leftmost; i <= rightmost; i++)
    {
      if ((tongueProfile[i] != INVALID) && (lowerProfile[i] != INVALID) && 
          (tongueProfile[i] > lowerProfile[i])) 
      { 
        lowerProfile[i] = tongueProfile[i]; 
        lowerProfileSurface[i] = tongueProfileSurface[i];
      }
    }

    while ((leftmost < N2)  && (lowerProfile[leftmost] > upperProfile[leftmost]-0.1)) { leftmost++; }
    while ((rightmost > N2) && (lowerProfile[rightmost] > upperProfile[rightmost]-0.1)) { rightmost--; }

    // **************************************************************
    // The lower profile samples must always lie below the upper 
    // profile samples.
    // **************************************************************

    for (i=0; i < N; i++)
    {
      if ((i >= leftmost) && (i <= rightmost))
      {
        if (lowerProfile[i] > upperProfile[i]) { lowerProfile[i] = upperProfile[i]; }
      }
      else
      {
        lowerProfile[i] = INVALID;
        upperProfile[i] = INVALID;
      }
    }
  }

  // ****************************************************************
  // Go from the middle to the left until the first open part 
  // (greater than a certain minimum distance) is found in the profile. 
  // This is mostly directly in the middle (apart from /l/).
  // Then declare the profile as valid only as long as it remains 
  // open towards the sides. More than one closure within the
  // profile are not allowed.
  // ****************************************************************

  const double OPEN_THRESHOLD = 0.1;   // = 1 mm
  const double CUTOFF_THRESHOLD = 0.01;

  i = N2-1;
  while ((upperProfile[i] - lowerProfile[i] < OPEN_THRESHOLD) && (i > 0)) { i--; }
  while ((upperProfile[i] - lowerProfile[i] >= CUTOFF_THRESHOLD) && (i > 0)) { i--; }
  while (i >= 0)
  {
    upperProfile[i] = INVALID;
    lowerProfile[i] = INVALID;
    i--;
  }

  i = N2-1;
  while ((upperProfile[i] - lowerProfile[i] < OPEN_THRESHOLD) && (i < N-1)) { i++; }
  while ((upperProfile[i] - lowerProfile[i] >= CUTOFF_THRESHOLD) && (i < N-1)) { i++; }
  while (i <= N-1)
  {
    upperProfile[i] = INVALID;
    lowerProfile[i] = INVALID;
    i++;
  }

  // ****************************************************************
  // Make the connection between both profiles.
  // ****************************************************************

  leftmost = 0;
  while (((upperProfile[leftmost] == INVALID) || (lowerProfile[leftmost] == INVALID)) &&
         (leftmost < N2-1)) { leftmost++; }
  if ((upperProfile[leftmost] != INVALID) && (lowerProfile[leftmost] != INVALID))
  {
    upperProfile[leftmost] = lowerProfile[leftmost] = 0.5*(upperProfile[leftmost] + lowerProfile[leftmost]);
  }

  rightmost = N-1;
  while (((upperProfile[rightmost] == INVALID) || (lowerProfile[rightmost] == INVALID)) &&
         (rightmost > N2)) { rightmost--; }
  if ((upperProfile[rightmost] != INVALID) && (lowerProfile[rightmost] != INVALID))
  {
    upperProfile[rightmost] = lowerProfile[rightmost] = 0.5*(upperProfile[rightmost] + lowerProfile[rightmost]);
  }

  // ****************************************************************
  // Determine the articulator that determines the lower profile
  // in the midsagittal plane.
  // Always test the two most median samples.
  // ****************************************************************

  articulator = Tube::OTHER_ARTICULATOR;
  
  bool hasTongue = false;
  bool hasLowerLip = false;
  bool hasLowerTeeth = false;

  double CHECK_RANGE_CM = 0.5;
  int NUM_CHECK_SAMPLES = (int)(CHECK_RANGE_CM / PROFILE_SAMPLE_LENGTH);

  for (i=1; i < NUM_CHECK_SAMPLES; i++)
  {
    if (lowerProfileSurface[N2-i] == TONGUE)
    {
      hasTongue = true;
    }
    else
    if ((lowerProfileSurface[N2-i] == LOWER_LIP) || (lowerProfileSurface[N2-i] == RADIATION))
    {
      hasLowerLip = true;
    }
    else
    if (lowerProfileSurface[N2-i] == LOWER_TEETH)
    {
      hasLowerTeeth = true;
    }

  }

  if (hasLowerTeeth)
  {
    articulator = Tube::LOWER_INCISORS;
  }
  else
  if (hasLowerLip)
  {
    articulator = Tube::LOWER_LIP;
  }
  else
  if (hasTongue)
  {
    articulator = Tube::TONGUE;
  }

}


// ****************************************************************************
// Samples the cutting line between P0 and P1 horizontally and inserts it into 
// the given upper profile.
// ****************************************************************************

void VocalTract::insertUpperProfileLine(Point2D P0, Point2D P1, int surfaceIndex, 
  double *upperProfile, int *upperProfileSurface)
{
  int i;
  Point2D v;

  // If both points have the same x-value, do nothing.
  if (P0.x == P1.x) { return; }

  // Move both points into the middle.
  P0.x+= 0.5*PROFILE_LENGTH;
  P1.x+= 0.5*PROFILE_LENGTH;

  // P0 must always be left of P1.
  if (P0.x > P1.x)
  {
    Point2D tempPoint = P0;
    P0 = P1;
    P1 = tempPoint;
  }

  // Extend the line by a little bit at both ends.
  const double LENGTH_INC = 0.01;
  v = P1 - P0;
  v.normalize();

  P0-= LENGTH_INC*v;
  P1+= LENGTH_INC*v;

  // Do the horizontal sampling.

  int firstSample = (int)(P0.x / PROFILE_SAMPLE_LENGTH);
  int lastSample  = (int)(P1.x / PROFILE_SAMPLE_LENGTH);

  if (firstSample == lastSample) { return; }

  double dx = PROFILE_SAMPLE_LENGTH;
  double dy = (P1.y - P0.y)*dx / (P1.x - P0.x);
  double y  = P0.y + ((firstSample + 1.0)*dx - P0.x)*dy/dx;

  for (i=firstSample+1; i <= lastSample; i++)
  {
    if ((i >= 0) && (i < NUM_PROFILE_SAMPLES))
    {
      if ((y <= upperProfile[i]) && (y >= MIN_PROFILE_VALUE) && (y <= MAX_PROFILE_VALUE))
      { 
        upperProfile[i] = y; 
        upperProfileSurface[i] = surfaceIndex;
      }
    }
    y+= dy;
  }
}


// ****************************************************************************
// Samples the cutting line between P0 and P1 horizontally and inserts it into 
// the given lower profile.
// ****************************************************************************

void VocalTract::insertLowerProfileLine(Point2D P0, Point2D P1, int surfaceIndex, 
  double *lowerProfile, int *lowerProfileSurface)
{
  int i;
  Point2D v;

  // If both points have the same x-value, do nothing.
  if (P0.x == P1.x) { return; }

  // Move both points into the middle.
  P0.x+= 0.5*PROFILE_LENGTH;
  P1.x+= 0.5*PROFILE_LENGTH;

  // P0 must always be left of P1.
  if (P0.x > P1.x)
  {
    Point2D tempPoint = P0;
    P0 = P1;
    P1 = tempPoint;
  }

  // Extend the line by a little bit at both ends.
  const double LENGTH_INC = 0.01;
  v = P1 - P0;
  v.normalize();

  P0-= LENGTH_INC*v;
  P1+= LENGTH_INC*v;

  // Do the horizontal sampling.

  int firstSample = (int)(P0.x / PROFILE_SAMPLE_LENGTH);
  int lastSample  = (int)(P1.x / PROFILE_SAMPLE_LENGTH);

  if (firstSample == lastSample) { return; }

  double dx = PROFILE_SAMPLE_LENGTH;
  double dy = (P1.y - P0.y)*dx / (P1.x - P0.x);
  double y  = P0.y + ((firstSample + 1.0)*dx - P0.x)*dy/dx;

  for (i=firstSample+1; i <= lastSample; i++)
  {
    if ((i >= 0) && (i < NUM_PROFILE_SAMPLES))
    {
      if ((y >= lowerProfile[i]) && (y >= MIN_PROFILE_VALUE) && (y <= MAX_PROFILE_VALUE)) 
      { 
        lowerProfile[i] = y; 
        lowerProfileSurface[i] = surfaceIndex;
      }
    }
    y+= dy;
  }
}



// ****************************************************************************
// Samples the cutting line between P0 and P1 horizontally and inserts it into 
// the given lower profile, and also in the upper profile, if the current
// upper profile is supplemented (but not replaced).
// ****************************************************************************

void VocalTract::insertLowerCoverProfileLine(Point2D P0, Point2D P1, int surfaceIndex,
  double *upperProfile, int *upperProfileSurface, double *lowerProfile, int *lowerProfileSurface)
{
  int i;
  Point2D v;

  // If both points have the same x-value, do nothing.
  if (P0.x == P1.x) { return; }

  // Move both points into the middle.
  P0.x+= 0.5*PROFILE_LENGTH;
  P1.x+= 0.5*PROFILE_LENGTH;

  // P0 must always be left of P1.
  if (P0.x > P1.x)
  {
    Point2D tempPoint = P0;
    P0 = P1;
    P1 = tempPoint;
  }

  // Extend the line by a little bit at both ends.
  const double LENGTH_INC = 0.01;
  v = P1 - P0;
  v.normalize();

  P0-= LENGTH_INC*v;
  P1+= LENGTH_INC*v;

  // Do the horizontal sampling.

  int firstSample = (int)(P0.x / PROFILE_SAMPLE_LENGTH);
  int lastSample  = (int)(P1.x / PROFILE_SAMPLE_LENGTH);

  if (firstSample == lastSample) { return; }

  double dx = PROFILE_SAMPLE_LENGTH;
  double dy = (P1.y - P0.y)*dx / (P1.x - P0.x);
  double y  = P0.y + ((firstSample + 1.0)*dx - P0.x)*dy/dx;

  for (i=firstSample+1; i <= lastSample; i++)
  {
    if ((i >= 0) && (i < NUM_PROFILE_SAMPLES))
    {
      if ((y >= MIN_PROFILE_VALUE) && (y <= MAX_PROFILE_VALUE)) 
      {
        if (y >= lowerProfile[i])
        { 
          lowerProfile[i] = y; 
          lowerProfileSurface[i] = surfaceIndex;
        }
        // Insert into the upper profile only if it was not defined yet at this place.
        if ((upperProfile[i] == EXTREME_PROFILE_VALUE) && (y <= upperProfile[i]))
        { 
          upperProfile[i] = y; 
          upperProfileSurface[i] = surfaceIndex;
        }
      }
    }
    y+= dy;
  }
}


// ****************************************************************************
/// Given the upper and lower profile of a cross-section, the area and
/// perimeter are calculated.
// ****************************************************************************

void VocalTract::getCrossSection(double *upperProfile, double *lowerProfile, CrossSection *section)
{
  const double INVALID = INVALID_PROFILE_SAMPLE;
  int i;
  double a, b;
  double deltaArea;
  const double d2 = PROFILE_SAMPLE_LENGTH*PROFILE_SAMPLE_LENGTH;

  // ****************************************************************
  // Calculate area and circumference.
  // ****************************************************************

  section->area = 0.0;
  section->circ = 0.0;

  for (i=0; i < NUM_PROFILE_SAMPLES-1; i++)
  {
    if ((upperProfile[i] != INVALID) && (upperProfile[i+1] != INVALID) &&
        (lowerProfile[i] != INVALID) && (lowerProfile[i+1] != INVALID))
    {
      a = upperProfile[i] - lowerProfile[i];
      b = upperProfile[i+1] - lowerProfile[i+1];
      deltaArea = 0.5*(a + b)*PROFILE_SAMPLE_LENGTH;
      section->area += deltaArea;

      a = upperProfile[i+1] - upperProfile[i];
      b = lowerProfile[i+1] - lowerProfile[i];
      section->circ+= sqrt(a*a + d2) + sqrt(b*b + d2);
    }
  }
}


// ****************************************************************************
/// Calculates for each cut vector on the center line the cross-sectional
/// profile.
// ****************************************************************************

bool VocalTract::exportCrossSections(const string &fileName)
{
  double upperProfile[NUM_PROFILE_SAMPLES];
  double lowerProfile[NUM_PROFILE_SAMPLES];
  Tube::Articulator articulator;
  int i, k;

  // ****************************************************************

  ofstream os(fileName);
  if (!os)
  {
    return false;
  }
  
  os << "# x, y (coordinates of the point on the centerline in cm)" << endl;
  os << "# n_x, n_y (coordinates of the normal of the point on the centerline in cm)" << endl;
  os << "# u0, 01, ..., u95 (samples of the upper profile in cm; 1000000 means 'invalid')" << endl;
  os << "# l0, l1, ..., l95 (samples of the lower profile in cm; 1000000 means 'invalid')" << endl;
  os << "# There are 129 slices." << endl;

  for (i = 0; i < NUM_CENTERLINE_POINTS; i++)
  {
    getCrossProfiles(centerLine[i].point, centerLine[i].normal, upperProfile, lowerProfile, true, articulator);

    os << centerLine[i].point.x << " " << centerLine[i].point.y << endl;
    os << centerLine[i].normal.x << " " << centerLine[i].normal.y << endl;

    for (k = 0; k < NUM_PROFILE_SAMPLES; k++)
    {
      os << upperProfile[k] << " ";
    }
    os << endl;

    for (k = 0; k < NUM_PROFILE_SAMPLES; k++)
    {
      os << lowerProfile[k] << " ";
    }
    os << endl;
  }

  // ****************************************************************
  // Close the file.
  // ****************************************************************

  os.close();

  return true;
}


// ****************************************************************************
// Save the vocal tract contour lines as SVG file.
// ****************************************************************************

bool VocalTract::exportTractContourSvg(const string &fileName, bool addCenterLine, bool addCutVectors)
{
  const double SCALE = 37.8;    // 1 cm in Corel Draw are 37.8 "default units" (pixels?)
  int indent = 0;
  int i;
  const string standardAttributes = "stroke=\"black\" stroke-width=\"1.5\" stroke-linecap=\"round\" "
    "stroke-linejoin=\"round\" fill=\"none\"";
  const string stippledAttribute = "stroke-dasharray=\"4.158, 4.158\"";
  Surface *s = NULL;
  Point3D Q;
  bool includeTeeth = true;
  bool includeRibs = false;

  ofstream os(fileName);

  if (!os)
  {
    printf("Error: Could not open the SVG file\n");
    return false;
  }

  // Write the "header and open the svg- and the group element.

  os << "<?xml version=\"1.0\" encoding=\"utf-8\"?>" << endl;
  os << "<svg width=\"100%\" height=\"100%\" viewBox=\"-60 -90 300 500\" version=\"1.1\" xmlns=\"http://www.w3.org/2000/svg\">" << endl;
  indent += 2;

  os << string(indent, ' ') << "<g>" << endl;
  indent += 2;

  // ****************************************************************
  // Draw all contour lines.
  // ****************************************************************

  // ****************************************************************
  // The upper contour.
  // ****************************************************************

  // First part of the upper cover
  os << string(indent, ' ') << "<polyline " << standardAttributes << " points=\"";  // Open the polyline tag and points attribute
  s = &surface[UPPER_COVER];
  addRibsSvg(os, s, 0, NUM_LARYNX_RIBS + NUM_PHARYNX_RIBS, s->numRibPoints - 1);
  os << "\"/>" << endl;      // Close the points attribute and the polyline tag

  // Uvula
  os << string(indent, ' ') << "<polyline " << standardAttributes << " points=\"";  // Open the polyline tag and points attribute
  s = &surface[UVULA];
  addRibsSvg(os, s, 0, s->numRibs - 1, s->numRibPoints - 1);
  addRibsSvg(os, s, s->numRibs - 1, 0, 0);
  os << "\"/>" << endl;      // Close the points attribute and the polyline tag

  // Last part of the upper cover
  os << string(indent, ' ') << "<polyline " << standardAttributes << " points=\"";  // Open the polyline tag and points attribute
  s = &surface[UPPER_COVER];
  addRibsSvg(os, s, NUM_LARYNX_RIBS + NUM_PHARYNX_RIBS + 1, s->numRibs - 1, s->numRibPoints - 1);

  // Upper incisors
  s = &surface[UPPER_TEETH];
  addRibPointsSvg(os, s, s->numRibs - 2, 0, 3);

  // Upper lip
  s = &surface[UPPER_LIP];
  addRibPointsSvg(os, s, s->numRibs - 1, 1, s->numRibPoints - 1);
  os << "\"/>" << endl;      // Close the points attribute and the polyline tag


  // ****************************************************************
  // Lower cover
  // ****************************************************************

  // Laryngeal part of the lower cover
  os << string(indent, ' ') << "<polyline " << standardAttributes << " points=\"";  // Open the polyline tag and points attribute
  s = &surface[LOWER_COVER];
  addRibsSvg(os, s, 0, NUM_LARYNX_RIBS - 1, s->numRibPoints - 1);
  os << "\"/>" << endl;      // Close the points attribute and the polyline tag

  // Epiglottis
  os << string(indent, ' ') << "<polyline " << standardAttributes << " points=\"";  // Open the polyline tag and points attribute
  s = &surface[EPIGLOTTIS];
  addRibsSvg(os, s, 0, s->numRibs - 1, s->numRibPoints - 1);
  addRibsSvg(os, s, s->numRibs - 1, 0, 0);
  os << "\"/>" << endl;      // Close the points attribute and the polyline tag

  // Rest of the lower cover
  os << string(indent, ' ') << "<polyline " << standardAttributes << " points=\"";  // Open the polyline tag and points attribute
  s = &surface[LOWER_COVER];
  addRibsSvg(os, s, NUM_LARYNX_RIBS, s->numRibs - 1, s->numRibPoints - 1);

  // Lower incisors
  s = &surface[LOWER_TEETH];
  addRibPointsSvg(os, s, s->numRibs - 2, 0, 3);

  // Lower lip
  s = &surface[LOWER_LIP];
  addRibPointsSvg(os, s, s->numRibs - 1, 1, s->numRibPoints - 1);
  os << "\"/>" << endl;      // Close the points attribute and the polyline tag

  // ****************************************************************
  // The ribs for the upper and lower cover.
  // ****************************************************************

  if (includeRibs)
  {
    // Ribs of the upper cover.
    s = &surface[UPPER_COVER];
    for (i = 0; i < s->numRibs; i++)
    {
      os << string(indent, ' ') << "<polyline " << standardAttributes << " points=\"";  // Open the polyline tag and points attribute
      addRibPointsSvg(os, s, i, 0, s->numRibPoints - 1);
      os << "\"/>" << endl;      // Close the points attribute and the polyline tag
    }

    // Ribs of the lower cover.
    s = &surface[LOWER_COVER];
    for (i = 0; i < s->numRibs; i++)
    {
      os << string(indent, ' ') << "<polyline " << standardAttributes << " points=\"";  // Open the polyline tag and points attribute
      addRibPointsSvg(os, s, i, 0, s->numRibPoints - 1);
      os << "\"/>" << endl;      // Close the points attribute and the polyline tag
    }

    // Border between the upper and lower cover.
    s = &surface[UPPER_COVER];
    os << string(indent, ' ') << "<polyline " << standardAttributes << " points=\"";  // Open the polyline tag and points attribute
    addRibsSvg(os, s, 0, NUM_LARYNX_RIBS + NUM_PHARYNX_RIBS +
      NUM_VELUM_RIBS - 1, 0);
    os << "\"/>" << endl;      // Close the points attribute and the polyline tag
  }

  // ****************************************************************
  // Tongue
  // ****************************************************************

  s = &surface[TONGUE];
  os << string(indent, ' ') << "<polyline " << standardAttributes << " points=\"";  // Open the polyline tag and points attribute
  addRibsSvg(os, s, 0, s->numRibs - 1, s->numRibPoints / 2);
  os << "\"/>" << endl;      // Close the points attribute and the polyline tag

  // Off the center line (1 cm to the right)
  os << string(indent, ' ') << "<polyline " << standardAttributes << " " << stippledAttribute;
  os << " points=\"";  // Open the polyline tag and points attribute
  addRibsSvg(os, s, 0, s->numRibs - 1, 1);
  os << "\"/>" << endl;      // Close the points attribute and the polyline tag

  // ****************************************************************
  // Teeth
  // ****************************************************************

  if (includeTeeth)
  {
    // Upper teeth
    s = &surface[UPPER_TEETH];
    os << string(indent, ' ') << "<polyline " << standardAttributes << " points=\"";  // Open the polyline tag and points attribute
    addRibsSvg(os, s, 0, s->numRibs - 1, 0);    // Upper inner edge
    os << "\"/>" << endl;      // Close the points attribute and the polyline tag

    os << string(indent, ' ') << "<polyline " << standardAttributes << " points=\"";  // Open the polyline tag and points attribute
    addRibsSvg(os, s, 0, s->numRibs - 5, 2);    // Lower outer edge
    os << "\"/>" << endl;      // Close the points attribute and the polyline tag

    os << string(indent, ' ') << "<polyline " << standardAttributes << " points=\"";  // Open the polyline tag and points attribute
    addRibPointsSvg(os, s, 0, 0, 1);          // Draw the most posterior rib completely
    os << "\"/>" << endl;      // Close the points attribute and the polyline tag

    // Lower teeth
    s = &surface[LOWER_TEETH];
    os << string(indent, ' ') << "<polyline " << standardAttributes << " points=\"";  // Open the polyline tag and points attribute
    addRibsSvg(os, s, 0, s->numRibs - 1, 0);    // Lower inner edge
    os << "\"/>" << endl;      // Close the points attribute and the polyline tag

    os << string(indent, ' ') << "<polyline " << standardAttributes << " points=\"";  // Open the polyline tag and points attribute
    addRibsSvg(os, s, 0, s->numRibs - 8, 2);    // Upper outer edge
    os << "\"/>" << endl;      // Close the points attribute and the polyline tag

    os << string(indent, ' ') << "<polyline " << standardAttributes << " points=\"";  // Open the polyline tag and points attribute
    addRibPointsSvg(os, s, 0, 0, 1);          // Draw the most posterior rib completely
    os << "\"/>" << endl;      // Close the points attribute and the polyline tag
  }


  // ****************************************************************
  // Include the center line.
  // ****************************************************************

  CenterLinePoint *cl = &centerLine[0];

  if (addCenterLine)
  {
    Point2D P0, P1;

    os << string(indent, ' ') << "<g>" << endl;
    indent += 2;

    for (i = 0; i < NUM_CENTERLINE_POINTS - 1; i++)
    {
      P0 = cl[i].point;
      P1 = cl[i + 1].point;

      os << string(indent, ' ')
        << "<line x1=\"" << SCALE*P0.x << "\" "
        << "y1=\"" << -SCALE*P0.y << "\" "
        << "x2=\"" << SCALE*P1.x << "\" "
        << "y2=\"" << -SCALE*P1.y << "\" "
        << "stroke=\"rgb(0,0,0)\" stroke-width=\"2.0\" stroke-linecap=\"round\" stroke-linejoin=\"round\"/>"
        << endl;
    }

    indent -= 2;
    os << string(indent, ' ') << "</g>" << endl;
  }

  // ****************************************************************
  // Include the normal vectors.
  // ****************************************************************

  if (addCutVectors)
  {
    Point2D P0, P1;
    double strokeWidth = 0.5;

    os << string(indent, ' ') << "<g>" << endl;
    indent += 2;

    for (i = 0; i < NUM_CENTERLINE_POINTS; i++)
    {
      P0 = cl[i].point + cl[i].min*cl[i].normal;
      P1 = cl[i].point + cl[i].max*cl[i].normal;

      os << string(indent, ' ')
        << "<line x1=\"" << SCALE*P0.x << "\" "
        << "y1=\"" << -SCALE*P0.y << "\" "
        << "x2=\"" << SCALE*P1.x << "\" "
        << "y2=\"" << -SCALE*P1.y << "\" "
        << "stroke=\"rgb(0,0,0)\" stroke-width=\"1.0\" stroke-linecap=\"round\" stroke-linejoin=\"round\"/>"
        << endl;
    }

    indent -= 2;
    os << string(indent, ' ') << "</g>" << endl;
  }

  // ****************************************************************
  // Finished drawing the contour lines.
  // ****************************************************************

  // Close the group element ****************************************
  indent -= 2;
  os << string(indent, ' ') << "</g>" << endl;

  // Close the svg element ******************************************
  indent -= 2;
  os << "</svg>" << endl;

  // Close the file *************************************************
  os.close();

  return true;
}


// ****************************************************************************
// Adds the point coordinates of some rib points to an output stream.
// ****************************************************************************

void VocalTract::addRibPointsSvg(ostream &os, Surface *s, int rib, int firstRibPoint, int lastRibPoint)
{
  if (s == NULL) { return; }
  if ((rib < 0) || (rib >= s->numRibs) || (firstRibPoint < 0) || (firstRibPoint >= s->numRibPoints) ||
    (lastRibPoint < 0) || (lastRibPoint >= s->numRibPoints)) {
    return;
  }

  const double SCALE = 37.8;    // 1 cm in Corel Draw are 37.8 "default units" (pixels?)
  int i;
  Point3D P;
  char st[256];

  if (lastRibPoint >= firstRibPoint)
  {
    for (i = firstRibPoint; i <= lastRibPoint; i++)
    {
      P = s->getVertex(rib, i);
      sprintf(st, "%2.4f %2.4f ", P.x*SCALE, -P.y*SCALE);
      os << st;
    }
  }
  else
  {
    for (i = firstRibPoint; i >= lastRibPoint; i--)
    {
      P = s->getVertex(rib, i);
      sprintf(st, "%2.4f %2.4f ", P.x*SCALE, -P.y*SCALE);
      os << st;
    }
  }
}


// ****************************************************************************
// Adds the point coordinates of some ribs to an output stream.
// ****************************************************************************

void VocalTract::addRibsSvg(ostream &os, Surface *s, int firstRib, int lastRib, int ribPoint)
{
  if (s == NULL) { return; }
  if ((ribPoint < 0) || (ribPoint >= s->numRibPoints) || (firstRib < 0) || (firstRib >= s->numRibs) ||
    (lastRib < 0) || (lastRib >= s->numRibs)) {
    return;
  }

  const double SCALE = 37.8;    // 1 cm in Corel Draw are 37.8 "default units" (pixels?)
  int i;
  Point3D P;
  char st[256];

  if (lastRib >= firstRib)
  {
    for (i = firstRib; i <= lastRib; i++)
    {
      P = s->getVertex(i, ribPoint);
      sprintf(st, "%2.4f %2.4f ", P.x*SCALE, -P.y*SCALE);
      os << st;
    }
  }
  else
  {
    for (i = firstRib; i >= lastRib; i--)
    {
      P = s->getVertex(i, ribPoint);
      sprintf(st, "%2.4f %2.4f ", P.x*SCALE, -P.y*SCALE);
      os << st;
    }
  }
}


// ****************************************************************************
/// Calculates the piecewise constant area function from the piecewise linear
/// area function (using the minimum norm).
// ****************************************************************************

void VocalTract::crossSectionsToTubeSections()
{
  const double EPSILON = 0.000001;
  int i, k, m;

  // ****************************************************************
  // Assingn the position and length to all discrete tube sections.
  // ****************************************************************

  for (i=0; i < NUM_PHARYNX_SECTIONS; i++)
  {
    tubeSection[i].pos = nasalPortPos_cm*(double)i / (double)NUM_PHARYNX_SECTIONS;
  }

  const int NUM_MOUTH_SECTIONS = NUM_TUBE_SECTIONS - NUM_PHARYNX_SECTIONS;
  double mouthLength = centerLineLength - nasalPortPos_cm;
  double firstSectionLength = mouthLength / 16;
  double lastSectionLength  = mouthLength / 32;
  double im = NUM_MOUTH_SECTIONS;
  double im2 = im*NUM_MOUTH_SECTIONS;
  double im3 = im2*NUM_MOUTH_SECTIONS;
  double im4 = im3*NUM_MOUTH_SECTIONS;

  double a1 = firstSectionLength;
  double a2 = ((mouthLength - a1*im)*3.0*im2 - im3*(lastSectionLength - a1)) / im4;
  double a3 = (im2*(lastSectionLength - a1) - 2.0*im*(mouthLength - a1*im)) / im4;
  
  for (i=NUM_PHARYNX_SECTIONS; i < NUM_TUBE_SECTIONS; i++)
  {
    k = i-NUM_PHARYNX_SECTIONS;
    tubeSection[i].pos = nasalPortPos_cm + a1*k + a2*k*k + a3*k*k*k;
  }

  for (i=0; i < NUM_TUBE_SECTIONS; i++)
  {
    if (i < NUM_TUBE_SECTIONS-1)
    {
      tubeSection[i].length = tubeSection[i+1].pos - tubeSection[i].pos;
    }
    else
    {
      tubeSection[i].length = centerLineLength - tubeSection[i].pos;
    }
  }

  // ****************************************************************
  // Run through all discrete tube sections.
  // ****************************************************************

  const double EXTREME = 1000000.0;
  int leftIndex  = 0;
  int rightIndex = 0;
  double leftX   = 0.0;
  double rightX  = 0.0;
  double y, deltaY;
  double length;

  for (i=0; i < NUM_TUBE_SECTIONS; i++)
  {
    leftX = tubeSection[i].pos;
    rightX = leftX + tubeSection[i].length;

    while ((crossSection[leftIndex].pos > leftX) && (leftIndex > 0)) { leftIndex--; }
    while ((crossSection[leftIndex+1].pos < leftX) && (leftIndex < NUM_CENTERLINE_POINTS-2)) { leftIndex++; }

    while ((crossSection[rightIndex].pos > rightX) && (rightIndex > 0)) { rightIndex--; }
    while ((crossSection[rightIndex+1].pos < rightX) && (rightIndex < NUM_CENTERLINE_POINTS-2)) { rightIndex++; }

    // **************************************************************
    // The articulator of a tube section is the one of the cross-
    // section with the smallest area.
    // **************************************************************

    tubeSection[i].articulator = Tube::OTHER_ARTICULATOR;

    // Both left and right index are left of the tube section.
    if (leftIndex == rightIndex)
    {
      tubeSection[i].articulator = crossSection[leftIndex].articulator;
    }
    else
    {
      // m is the index of the first cross-section IN the tube.
      m = leftIndex + 1;
      for (k=leftIndex + 1; k <= rightIndex; k++)
      {
        if (crossSection[k].area < crossSection[m].area)
        {
          m = k;
        }
      }

      tubeSection[i].articulator = crossSection[m].articulator;
    }


    // **************************************************************
    // Run through all corresponding indices on the center line
    // to determine the area and perimenter of the tube section.
    // **************************************************************

    tubeSection[i].area = EXTREME;
    tubeSection[i].circ = EXTREME;

    for (k=leftIndex; k <= rightIndex; k++)
    {
      // Area and circumference need a more precise determination.

      if ((k == leftIndex) && (k == rightIndex))
      {
        length = crossSection[leftIndex+1].pos - crossSection[leftIndex].pos;
        if (length < EPSILON) { length = EPSILON; }

        // Cross-sectional area.
        deltaY = crossSection[leftIndex+1].area - crossSection[leftIndex].area;
        y = crossSection[leftIndex].area + deltaY*(leftX  - crossSection[leftIndex].pos) / length;
        if (y < tubeSection[i].area) { tubeSection[i].area = y; }

        y = crossSection[leftIndex].area + deltaY*(rightX - crossSection[leftIndex].pos) / length;
        if (y < tubeSection[i].area) { tubeSection[i].area = y; }

        // Perimeter.
        deltaY = crossSection[leftIndex+1].circ - crossSection[leftIndex].circ;
        y = crossSection[leftIndex].circ + deltaY*(leftX  - crossSection[leftIndex].pos) / length;
        if (y < tubeSection[i].circ) { tubeSection[i].circ = y; }

        y = crossSection[leftIndex].circ + deltaY*(rightX - crossSection[leftIndex].pos) / length;
        if (y < tubeSection[i].circ) { tubeSection[i].circ = y; }
      }
      else

      if (k == leftIndex)
      {
        length = crossSection[leftIndex+1].pos - crossSection[leftIndex].pos;
        if (length < EPSILON) { length = EPSILON; }

        // Cross-sectional area.
        deltaY = crossSection[leftIndex+1].area - crossSection[leftIndex].area;
        y = crossSection[leftIndex].area + deltaY*(leftX  - crossSection[leftIndex].pos) / length;
        if (y < tubeSection[i].area) { tubeSection[i].area = y; }
       
        y = crossSection[leftIndex+1].area;
        if (y < tubeSection[i].area) { tubeSection[i].area = y; }

        // Perimeter.
        deltaY = crossSection[leftIndex+1].circ - crossSection[leftIndex].circ;
        y  = crossSection[leftIndex].circ + deltaY*(leftX  - crossSection[leftIndex].pos) / length;
        if (y < tubeSection[i].circ) { tubeSection[i].circ = y; }

        y = crossSection[leftIndex+1].circ;
        if (y < tubeSection[i].circ) { tubeSection[i].circ = y; }
      }
      else

      if (k == rightIndex)
      {
        length = crossSection[rightIndex+1].pos - crossSection[rightIndex].pos;
        if (length < EPSILON) { length = EPSILON; }

        // Cross-sectional area.
        deltaY = crossSection[rightIndex+1].area - crossSection[rightIndex].area;
        y  = crossSection[rightIndex].area;
        if (y < tubeSection[i].area) { tubeSection[i].area = y; }

        y = crossSection[rightIndex].area + deltaY*(rightX  - crossSection[rightIndex].pos) / length;
        if (y < tubeSection[i].area) { tubeSection[i].area = y; }

        // Perimeter.
        deltaY = crossSection[rightIndex+1].circ - crossSection[rightIndex].circ;
        y  = crossSection[rightIndex].circ;
        if (y < tubeSection[i].circ) { tubeSection[i].circ = y; }

        y = crossSection[rightIndex].circ + deltaY*(rightX  - crossSection[rightIndex].pos) / length;
        if (y < tubeSection[i].circ) { tubeSection[i].circ = y; }
      }
      else
      {
        // Cross-sectional area.
        if (crossSection[k].area < tubeSection[i].area) { tubeSection[i].area = crossSection[k].area; }
        if (crossSection[k+1].area < tubeSection[i].area) { tubeSection[i].area = crossSection[k+1].area; }

        // Perimeter.
        if (crossSection[k].circ < tubeSection[i].circ) { tubeSection[i].circ = crossSection[k].circ; }
        if (crossSection[k+1].circ < tubeSection[i].circ) { tubeSection[i].circ = crossSection[k+1].circ; }
      }
    }
  }

}


// ****************************************************************************
// Returns the index of the vocal tract shape with the given name.
// ****************************************************************************

int VocalTract::getShapeIndex(const string &name)
{
  int k, index = -1;
  for (k=0; (k < (int)shapes.size()) && (index == -1); k++)
  {
    if (name == shapes[k].name)
    {
      index = k; 
    }
  }
  return index;
}


// ****************************************************************************
/// Returns true, if the given name is that of a vocalic vocal tract shape.
/// that are all names that don't start with "tt-", "tb-", or "ll-".
// ****************************************************************************

bool VocalTract::isVowelShapeName(const string &name)
{
  if (name.length() < 3)
  {
    return true;
  }

  string firstChars = name.substr(0, 3);
  if ((firstChars == "tt-") || (firstChars == "tb-") || (firstChars == "ll-"))
  {
    return false;
  }
  else
  {
    return true;
  }
}


// ****************************************************************************
// Returns the x-coord. of a point lying on the rear pharynx border line at the
// height y.
// ****************************************************************************

double VocalTract::getPharynxBackX(double y)
{
  const double MIN_ANGLE_DEG = -135.0;
  const double MAX_ANGLE_DEG = -45.0;

  Point2D A = anatomy.pharynxFulcrum;
  double angle_deg = anatomy.pharynxRotationAngle_deg;

  if (angle_deg > 0.0) { angle_deg-= 2.0*M_PI; }
  if (angle_deg < MIN_ANGLE_DEG) { angle_deg = MIN_ANGLE_DEG; }
  if (angle_deg > MAX_ANGLE_DEG) { angle_deg = MAX_ANGLE_DEG; }
  double angle_rad = angle_deg*M_PI/180.0;

  return A.x + (y-A.y)*cos(angle_rad) / sin(angle_rad);
}


// ****************************************************************************
/// Writes the geometry of the tube sections, the teeth position and the
/// tongue tip side elevation into the tube object.
// ****************************************************************************

void VocalTract::getTube(Tube *tube)
{
  const int N = Tube::NUM_PHARYNX_MOUTH_SECTIONS;

  double length_cm[N];
  double area_cm2[N];
  Tube::Articulator articulator[N];

  TubeSection *ts = NULL;
  int i;

  tube->initPiriformFossa(anatomy.piriformFossaLength_cm, anatomy.piriformFossaVolume_cm3);
  tube->initSubglottalCavity(anatomy.subglottalCavityLength_cm);
  tube->initNasalCavity(anatomy.nasalCavityLength_cm);

  for (i=0; i < Tube::NUM_PHARYNX_MOUTH_SECTIONS; i++)
  {
    ts = &tubeSection[i];
    length_cm[i] = ts->length;
    area_cm2[i]  = ts->area;
    articulator[i] = ts->articulator;
  }

  double tongueTipSideElevation = param[TS3].x;

  tube->setPharynxMouthGeometry(length_cm, area_cm2, articulator, 
    incisorPos_cm, tongueTipSideElevation);
  tube->setVelumOpening(nasalPortArea_cm2);
}


// ****************************************************************************
/// Returns the position on the center line that is nearest the point Q.
// ****************************************************************************

double VocalTract::getCenterLinePos(Point2D Q, int &bestIndex, double &bestT)
{
  const double EPSILON = 0.000001;
  int i;
  double minDist = 1000000.0;
  double dist;
  bestIndex = 0;
  bestT = 0.0;
  Point2D P0, P1, n0, n1, v, m, R;
  double denominator;
  double t, t0, t1;
  double p, q;

  // Default setting.

  bestIndex = -1;
  bestT = 0.0;

  // ****************************************************************
  // Run through all line segments.
  // ****************************************************************

  for (i=0; i < NUM_CENTERLINE_POINTS-1; i++)
  {
    P0 = centerLine[i].point;
    n0 = centerLine[i].normal;
    P1 = centerLine[i+1].point;
    n1 = centerLine[i+1].normal;

    // Separate P0 and P1 by a little bit to avoid numeric errors.

    v = (P1 - P0)*EPSILON;
    P0-= v;
    P1+= v;

    // Is Q right of P0+t*n0 and left of P1+t*n1 ?

    if (((Q.x-P0.x)*n0.y - n0.x*(Q.y-P0.y) >= 0.0) && ((Q.x-P1.x)*n1.y - n1.x*(Q.y-P1.y) <= 0.0))
    {
      v = P1 - P0;
      R = Q - P0;
      m = n1 - n0;
      
      denominator = v.y*m.x - v.x*m.y;
      if (denominator == 0.0) { denominator = EPSILON; }
      p = (R.x*m.y - R.y*m.x + v.y*n0.x - v.x*n0.y) / denominator;
      q = (R.x*n0.y - R.y*n0.x) / denominator;

      q = 0.25*p*p - q;
      if (q < 0.0) { q = 0.0; }
      q = sqrt(q);
      t0 = -0.5*p + q;
      t1 = -0.5*p - q;

      if ((t0 > -EPSILON) && (t0 < 1.0+EPSILON)) { t = t0; } else { t = t1; }

      // Check the distance between Q and the point on the center line.
      R = P0 + v*t;
      dist = (Q-R).magnitude();

      if (dist < minDist)
      {
        minDist = dist;
        bestIndex = i;
        bestT = t;
      }
    }
  }

  // ****************************************************************
  // If no position was found yet, check the distance to the
  // center line points.
  // ****************************************************************

  if (bestIndex == -1)
  {
    bestIndex = 0;

    for (i=0; i < NUM_CENTERLINE_POINTS; i++)
    {
      dist = (centerLine[i].point - Q).magnitude();

      if (dist < minDist)
      {
        minDist = dist;
        bestIndex = i;
        bestT = 0.0;
      }
    }
  }

  if (bestIndex == NUM_CENTERLINE_POINTS-1)
  {
    bestIndex = NUM_CENTERLINE_POINTS-2;
    bestT = 1.0;
  }

  return centerLine[bestIndex].pos + bestT*(centerLine[bestIndex+1].pos - centerLine[bestIndex].pos);
}


// ****************************************************************************
/// Returns the point and normal vector at an arbitrary position along the
/// center line.
// ****************************************************************************

void VocalTract::getCutVector(double pos, Point2D &P, Point2D &v)
{
  if (pos < 0.0) { pos = 0.0; }
  if (pos > centerLineLength) { pos = centerLineLength; }

  int i = (int)(NUM_CENTERLINE_POINTS*(pos/centerLineLength));   // First estimate.

  if (i < 0) { i = 0; }
  if (i >= NUM_CENTERLINE_POINTS-1) { i = NUM_CENTERLINE_POINTS-2; }

  while ((i < NUM_CENTERLINE_POINTS-2) && (pos > centerLine[i+1].pos)) { i++; }
  while ((i > 0) && (pos < centerLine[i].pos)) { i--; }

  double diff = centerLine[i+1].pos - centerLine[i].pos;
  if (diff == 0.0) { diff = 0.000001; }
  double t = (pos - centerLine[i].pos) / diff;

  P = (1.0-t)*centerLine[i].point + t*centerLine[i+1].point;
  v = (1.0-t)*centerLine[i].normal + t*centerLine[i+1].normal;

  v.normalize();
}


// ****************************************************************************
/// Restricts the value of the given parameter to its value range.
// ****************************************************************************

void VocalTract::restrictParam(int index)
{
  if (param[index].x < param[index].min) 
  { 
    param[index].x = param[index].min; 
  }
  
  if (param[index].x > param[index].max) 
  { 
    param[index].x = param[index].max; 
  }
}

// ****************************************************************************
// ****************************************************************************

bool VocalTract::hasUnsavedChanges()
{
  return false;
}


// ****************************************************************************
// ****************************************************************************

void VocalTract::clearUnsavedChanges()
{
  // ...
}


// ****************************************************************************
/// Temporarily store (cache) the control parameter values, so that they can 
/// restored later.
// ****************************************************************************

void VocalTract::storeControlParams()
{
  hasStoredControlParams = true;
  int i;
  for (i = 0; i < NUM_PARAMS; i++)
  {
    storedControlParams[i] = param[i].x;
  }
}


// ****************************************************************************
/// Restore the temporarily stored (cached) control parameter values.
// ****************************************************************************

void VocalTract::restoreControlParams()
{
  if (hasStoredControlParams == false)
  {
    return;
  }

  int i;
  for (i = 0; i < NUM_PARAMS; i++)
  {
    param[i].x = storedControlParams[i];
  }
  hasStoredControlParams = false;

  calculateAll();
}


// ****************************************************************************

