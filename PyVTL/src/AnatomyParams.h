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

#ifndef __ANATOMY_PARAMS_H__
#define __ANATOMY_PARAMS_H__

#include <string>
#include "VocalTract.h"
#include "Geometry.h"

using namespace std;

// ****************************************************************************
// The units for the parameters are cm and degree (for angles).
// ****************************************************************************

class AnatomyParams
{
  // **************************************************************************
  // Public data.
  // **************************************************************************

public:
  enum AnatomyParameters
  {
    LIP_WIDTH,
    MANDIBLE_HEIGHT,
    LOWER_MOLARS_HEIGHT,
    UPPER_MOLARS_HEIGHT,
    PALATE_HEIGHT,
    PALATE_DEPTH,
    HARD_PALATE_LENGTH,
    SOFT_PALATE_LENGTH,
    PHARYNX_LENGTH,
    LARYNX_LENGTH,
    LARYNX_WIDTH,
    VOCAL_FOLD_LENGTH,
    ORAL_PHARYNGEAL_ANGLE,
    NUM_ANATOMY_PARAMS
  };

  struct Param
  {
    string name;
    string abbr;
    string unit;
    double min;
    double max;
    double x;
  };

  Param param[NUM_ANATOMY_PARAMS];

  // **************************************************************************
  // Public functions.
  // **************************************************************************

public:
  AnatomyParams();
  void restrictParams();
  void calcFromAge(int age_month, bool isMale);
  void getFrom(VocalTract *tract);
  void setFor(VocalTract *tract);
  void adjustTongueRootCalculation(VocalTract *targetVocalTract);
  bool loadReferenceVocalTract(const string &fileName);

  double transformX(AnatomyParams *origAnatomy, double origX);
  double transformY(AnatomyParams *origAnatomy, double origY);
  double transformY(AnatomyParams *origAnatomy, double origX, double origY);

  void adaptArticulation(double *oldParams, double *newParams);

  // **************************************************************************
  // Private data.
  // **************************************************************************

private:
  static VocalTract *referenceVocalTract;
};

#endif
