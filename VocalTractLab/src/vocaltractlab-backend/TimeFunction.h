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

#ifndef __TIME_FUNCTION_H__
#define __TIME_FUNCTION_H__

#include <vector>

using namespace std;

// ****************************************************************************
/// This class represents a time function that is defined by a set of time-
/// value nodes. Between the nodes, the function is linearly interpolated.
// ****************************************************************************

class TimeFunction
{
  // **************************************************************************
  // Public data.
  // **************************************************************************

public:
  struct Node
  {
    double time;
    double value;
  };

  // **************************************************************************
  /// Public functions.
  // **************************************************************************

public:
  TimeFunction();
  bool setNodes(const Node n[], const int numNodes);
  bool setNodes(const vector<Node> n);
  void getNodes(vector<Node> &nodes);
  double getValue(double t);
  
  static void test();

  // **************************************************************************
  // Private data.
  // **************************************************************************

private:
  vector<Node> nodes;
};


#endif
