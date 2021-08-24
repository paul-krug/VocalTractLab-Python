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

#ifndef __GLOTTIS_H__
#define __GLOTTIS_H__

#include <ostream>
#include <string>
#include <vector>

#include "XmlNode.h"

using namespace std;


// ****************************************************************************
// The is the base class for glottis models.
// ****************************************************************************

class Glottis
{
  // **************************************************************************
  // Public data.
  // **************************************************************************

public:
  static const int FREQUENCY = 0;
  static const int PRESSURE = 1;
  static const int MAX_CONTROL_PARAMS = 32;
  static const double DEFAULT_ASPIRATION_STRENGTH_DB;

  struct Parameter
  {
    string name;
    string abbr;
    string cgsUnit;     ///< CGS unit of this parameter (used for its values)
    double factor;      ///< userUnit = factor * cgsUnit
    string userUnit;    ///< Unit displayed for the user
    double min;
    double max;
    double neutral;
    double x;           ///< Current value of the parameter
  };

  // Structure for a shape like "normal phonation", "breathy phonation", or "open glottis"
  struct Shape
  {
    string name;
    vector<double> controlParam;
  };

  // Data of the last saved state.
  struct SavedState
  {
    vector<double> staticParam;
    vector<Shape> shape;
  };

  // **************************************************************************

  vector<Parameter> staticParam;
  vector<Parameter> controlParam;
  vector<Parameter> derivedParam;
  vector<Shape> shape;

  SavedState savedState;

  // **************************************************************************
  // Public functions.
  // **************************************************************************

public:
  virtual ~Glottis() {}
  virtual string getName() = 0;
  virtual void resetMotion() = 0;
  /// Requires four pressure values: subglottal, lower glottis, upper glottis, supraglottal
  virtual void incTime(const double timeIncrement_s, const double pressure_dPa[]) = 0;
  virtual void calcGeometry() = 0;
  virtual void getTubeData(double *length_cm, double *area_cm2) = 0;
  virtual int getApertureParamIndex() = 0;
  virtual double getAspirationStrength_dB();

  Shape *getShape(const string &name);
  bool hasUnsavedChanges();
  void clearUnsavedChanges();
  bool writeToXml(ostream &os, int initialIndent, bool isSelected);
  bool readFromXml(XmlNode &node);

  void printParamNames(ostream &os);
  void printParamValues(ostream& os, double glottalFlow_cm3_s,
    double glottalPressure_dPa[], double mouthFlow_cm3_s, double nostrilFlow_cm3_s,
    double skinFlow_cm3_s, double radiatedPressure_dPa);
  void restrictParams(vector<Parameter> &p);

  void storeControlParams();
  void restoreControlParams();

  // **************************************************************************
  // Private data.
  // **************************************************************************

private:
  bool hasStoredControlParams;
  double storedControlParams[MAX_CONTROL_PARAMS];
};


#endif
