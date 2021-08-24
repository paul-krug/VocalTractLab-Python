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

#include "TimeFunction.h"
#include <cstdio>


// ****************************************************************************
/// Constructor. Create an empty node list.
// ****************************************************************************

TimeFunction::TimeFunction()
{
  nodes.clear();
}


// ****************************************************************************
// ****************************************************************************

bool TimeFunction::setNodes(const Node n[], const int numNodes)
{
  vector<Node> temp;
  int i;

  for (i=0; i < numNodes; i++)
  {
    temp.push_back(n[i]);
  }
  
  return setNodes(temp);
}

// ****************************************************************************
/// Set the given node list. If it is not a valid list, the current node list
/// is cleared.
// ****************************************************************************

bool TimeFunction::setNodes(const vector<Node> n)
{
  nodes = n;

  // Check that the times of the nodes are in ascending order.
  
  int i;
  bool ok = true;
  int N = (int)nodes.size();
  
  for (i=1; i < N; i++)
  {
    if (nodes[i].time < nodes[i-1].time)
    {
      ok = false;
      break;
    }
  }

  if (ok == false)
  {
    nodes.clear(); 
    printf("ERROR: Invalid node list for time function!\n");
  }

  return ok;
}


// ****************************************************************************
/// Return the time-value pairs of the current node list.
// ****************************************************************************

void TimeFunction::getNodes(vector<Node> &nodes)
{
  nodes = this->nodes;
}


// ****************************************************************************
// ****************************************************************************

double TimeFunction::getValue(double t)
{
  if (nodes.empty())
  {
    return 0.0;
  }

  if (nodes.size() < 2)
  {
    return nodes[0].value;
  }

  if (t < nodes[0].time)
  {
    return nodes[0].value;
  }

  if (t >= nodes[nodes.size()-1].time)
  {
    return nodes[nodes.size()-1].value;
  }

  // ****************************************************************
  // Find the two nodes that form the interval where t lies.
  // ****************************************************************

  int N = (int)nodes.size();
  int leftIndex = 0;
  int rightIndex = N-1;
  int centerIndex;

  while ((rightIndex - leftIndex > 1))
  {
    centerIndex = leftIndex + (rightIndex - leftIndex) / 2;
    if (nodes[centerIndex].time < t)
    {
      leftIndex = centerIndex;
    }
    else
    {
      rightIndex = centerIndex;
    }
  }

  // Make linear interpolation between leftIndex and rightIndex

  const double EPSILON = 0.000000001;
  double deno = nodes[rightIndex].time - nodes[leftIndex].time;
  if (deno < EPSILON)
  {
    deno = EPSILON;
  }

  double a = (t - nodes[leftIndex].time) / deno;
  double value = nodes[leftIndex].value + a*(nodes[rightIndex].value - nodes[leftIndex].value);

  return value;
}


// ****************************************************************************
/// Unit test.
// ****************************************************************************

void TimeFunction::test()
{
  printf("\nTimeFunction unit test\n");
  printf("======================\n");

  TimeFunction func;
  Node n[] = { {-1, -1}, {1, 1}, {2, -2} };
  func.setNodes(n, 3);

  vector<Node> no;
  func.getNodes(no);

  int i;
  for (i=0; i < (int)no.size(); i++)
  {
    printf("i=%d: (%f %f)\n", i, no[i].time, no[i].value);
  }

  double x;
  x = -1.01;
  printf("f(%f) = %f\n", x, func.getValue(x));
  x = -1.0;
  printf("f(%f) = %f\n", x, func.getValue(x));
  x = -0.9;
  printf("f(%f) = %f\n", x, func.getValue(x));
  x = 0.0;
  printf("f(%f) = %f\n", x, func.getValue(x));
  x = 1.5;
  printf("f(%f) = %f\n", x, func.getValue(x));
  x = 2.5;
  printf("f(%f) = %f\n", x, func.getValue(x));

  printf("\n");
}

// ****************************************************************************

