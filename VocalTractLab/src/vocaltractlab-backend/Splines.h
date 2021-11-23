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

#ifndef __SPLINES_H__
#define __SPLINES_H__

#include "Geometry.h"

const int MAX_SPLINE_POINTS = 256;
const int MAX_SPLINE_SEGMENTS = MAX_SPLINE_POINTS - 1;

// ****************************************************************************
// Base class for different types of splines that work with multiple control
// points.
// ****************************************************************************

class Spline3D
{
public:
  Spline3D();
  virtual ~Spline3D() {}
  Spline3D(int newNumPoints, const Point3D *points);
  Spline3D(int newNumPoints, const Point3D *points, const double *weights);

  void reset(int newNumPoints);
  void setPoints(int newNumPoints, const Point3D *points);
  void setPoints(int newNumPoints, const Point3D *points, const double *weights);
  void setPoint(int index, Point3D point, double weight = 1.0);
  void addPoint(Point3D point, double weight = 1.0);

  Point3D getControlPoint(int index);
  Point3D getControlPoint(int index, double &weight);
  virtual Point3D getPoint(double t);
  virtual Point3D getTangent(double t);
  double getUniformParam(double t);
  double getIntersection(const Point3D planePoint, const Point3D planeNormal, double tMin, double tMax);

  int getNumPoints() { return numPoints; }

protected:
  Point3D P[MAX_SPLINE_POINTS];     // control points
  double  w[MAX_SPLINE_POINTS];     // weights
  int  numPoints;
  bool pointsChanged;
};

// ****************************************************************************
// A simple piecewise linear interpolation of control points.
// ****************************************************************************

class LineStrip3D : public Spline3D
{
public:
  LineStrip3D();
  LineStrip3D(int newNumPoints, Point3D *points);

  Point3D getPoint(double t);
  double  getCurveParam(int index);
  double getIntersection(const Point3D planePoint, const Point3D planeNormal);

private:
  void calculateParams();
  double pos[MAX_SPLINE_POINTS];    // 0 <= pos[i] <= 1
};


// ****************************************************************************
// Rational (weighted) Bezier curves.
// ****************************************************************************

class BezierCurve3D : public Spline3D
{
public:
  BezierCurve3D();
  BezierCurve3D(int newNumPoints, const Point3D *points);
  BezierCurve3D(int newNumPoints, const Point3D *points, const double *weights);

  Point3D getPoint(double t);

private:
  void getBernsteinCoeff(int i, int n, double *coeff);
  void calculateCoeff();

  Point3D A[MAX_SPLINE_POINTS];   // Coefficients of the numerator polynomial
  double  B[MAX_SPLINE_POINTS];   // Coefficients of the denominator polynomial
};

// ****************************************************************************
// A simple piecewise linear interpolation of control points in 2D.
// ****************************************************************************

class LineStrip2D
{
public:
  LineStrip2D();
  LineStrip2D(int newNumPoints, Point2D *points);

  void reset(int newNumPoints);
  void setPoints(int newNumPoints, const Point2D *points);
  void setPoint(int index, Point2D point);
  void addPoint(Point2D point);
  void delPoint();

  Point2D getControlPoint(int index);
  Point2D getPoint(double t);
  double  getFunctionValue(double x);
  Point2D getTangent(double t);
  double  getCurveParam(int index);
  bool    getClosestIntersection(const Point2D Q, const Point2D v, double &t, Point2D &intersection);
  bool    getFirstIntersection(const Point2D Q, const Point2D v, double &t, Point2D &intersection);
  bool    getSpecialIntersection(const Point2D Q, const Point2D v, double &t, Point2D &intersection);

  int getNumPoints() { return numPoints; }

private:
  void calculateParams();
  
  Point2D P[MAX_SPLINE_POINTS];     // Control points
  double  pos[MAX_SPLINE_POINTS];    // 0 <= pos[i] <= 1
  int     numPoints;
  bool    pointsChanged;
};

// ****************************************************************************
#endif
