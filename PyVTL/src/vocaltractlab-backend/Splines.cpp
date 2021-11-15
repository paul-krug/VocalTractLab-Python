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

#include "Splines.h"
#include <cstdio>
#include <cmath>

// ----------------------------------------------------------------------------
// Base class for all 3D splines.
// ----------------------------------------------------------------------------

// ****************************************************************************
// ****************************************************************************

Spline3D::Spline3D()
{
  setPoints(0, NULL);
}

// ****************************************************************************
// ****************************************************************************

Spline3D::Spline3D(int newNumPoints, const Point3D *points)
{
  setPoints(newNumPoints, points);
}

// ****************************************************************************
// ****************************************************************************

Spline3D::Spline3D(int newNumPoints, const Point3D *points, const double *weights)
{
  setPoints(newNumPoints, points, weights);
}

// ****************************************************************************
// Set the number of control points and init. them with zeros.
// ****************************************************************************

void Spline3D::reset(int newNumPoints)
{
  numPoints = newNumPoints;
  if (numPoints > MAX_SPLINE_POINTS) { numPoints = MAX_SPLINE_POINTS; }
  if (numPoints < 0) { numPoints = 0; }

  for (int i=0; i < numPoints; i++)
  {
    P[i].set(0.0, 0.0, 0.0);
    w[i] = 1.0;         // Default-Gewicht
  }
  pointsChanged = true;
}

// ****************************************************************************
// Set the number of control points and init. them.
// ****************************************************************************

void Spline3D::setPoints(int newNumPoints, const Point3D *points)
{
  numPoints = newNumPoints;
  if (numPoints > MAX_SPLINE_POINTS) { numPoints = MAX_SPLINE_POINTS; }
  if ((numPoints < 0) || (points == NULL)) { numPoints = 0; }

  for (int i=0; i < numPoints; i++)
  {
    P[i] = points[i];
    w[i] = 1.0;         // Default-Gewicht
  }
  pointsChanged = true;
}


// ****************************************************************************
// Set the number of control points and init. them.
// ****************************************************************************

void Spline3D::setPoints(int newNumPoints, const Point3D *points, const double *weights)
{
  numPoints = newNumPoints;
  if (numPoints > MAX_SPLINE_POINTS) { numPoints = MAX_SPLINE_POINTS; }
  if ((numPoints < 0) || (points == NULL) || (weights == NULL)) { numPoints = 0; }

  for (int i=0; i < numPoints; i++)
  {
    P[i] = points[i];
    w[i] = weights[i];
  }
  pointsChanged = true;
}

// ****************************************************************************
// ****************************************************************************

void Spline3D::setPoint(int index, Point3D point, double weight)
{
  if ((index < 0) || (index >= numPoints)) { return; }

  P[index] = point;
  w[index] = weight;
  pointsChanged = true;
}

// ****************************************************************************
// Add a control point ot the end of the list.
// ****************************************************************************

void Spline3D::addPoint(Point3D point, double weight)
{
  if (numPoints >= MAX_SPLINE_POINTS) { return; }

  P[numPoints] = point;
  w[numPoints] = weight;

  numPoints++;
  pointsChanged = true;
}

// ****************************************************************************
// Returns a single control point.
// ****************************************************************************

Point3D Spline3D::getControlPoint(int index)
{
  if ((index < 0) || (index >= numPoints)) 
    { return Point3D(0.0, 0.0, 0.0); }
  else
    { return P[index]; }
}

// ****************************************************************************
// Returns a single control point.
// ****************************************************************************

Point3D Spline3D::getControlPoint(int index, double &weight)
{
  if ((index < 0) || (index >= numPoints)) 
  { 
    weight = 1.0;
    return Point3D(0.0, 0.0, 0.0); 
  }
  else
  { 
    weight = w[index];
    return P[index]; 
  }
}

// ****************************************************************************
// Returns the point at the relative position t in [0..1].
// ****************************************************************************

Point3D Spline3D::getPoint(double t)
{
  if (numPoints < 2) { return Point3D(0.0, 0.0, 0.0); }

  // Find the right segment.

  int numSegments = numPoints - 1;
  int segment = (int)(t*(double)numSegments);

  if (segment < 0) { segment = 0; }
  if (segment >= numSegments) { segment = numSegments - 1; }

  double deltaT = 1.0 / numSegments;
  double u = (t - segment*deltaT) / deltaT;   // Local parameter for the segment

  if (u < 0.0) { u = 0.0; }
  if (u > 1.0) { u = 1.0; }
  
  double u1 = 1.0 - u;

  return (P[segment]*u1 + P[segment+1]*u);
}

// ****************************************************************************
// The parameter t with 0 <= t <= 1 is the ratio of the curve length between
// the beginning and the seeked position and the total length of the curve.
// Returned is a value u that can be directly passed as curve parameter to 
// getPoint() or getTangent().
// ****************************************************************************

double Spline3D::getUniformParam(double t)
{
  const double EPSILON = 0.000001;
  const int N = 100;
  Point3D P[N];
  double s[N];
  double pos[N];
  int i;
  double l;

  // N Punkte auf der Kurve in nichtgleichfoermigen Abstaenden berechnen.

  for (i=0; i < N; i++)
  {
    s[i] = (double)i / (double)(N-1);
    P[i] = getPoint(s[i]);
    if (i == 0)
    {
      pos[i] = 0.0;
    }
    else
    {
      pos[i] = pos[i-1] + (P[i] - P[i-1]).magnitude();
    }
  }

  // Alle Positionen auf die Gesamtlaenge normieren.

  l = pos[N-1];
  if (l < EPSILON) { l = EPSILON; }

  for (i=0; i < N; i++) { pos[i]/= l; }

  // Zwischen welche beiden Positionen faellt t ?

  if (t < 0.0) { t = 0.0; }
  if (t > 1.0) { t = 1.0; }

  int k = -1;
  double ratio = 0.0;

  for (i=0; i < N-1; i++) 
  { 
    if ((t >= pos[i]) && (t <= pos[i+1])) 
    { 
      k = i; 
      l = pos[i+1] - pos[i];
      if (l < EPSILON) { l = EPSILON; }
      ratio = (t - pos[i]) / l;
    }
  }

  if (k == -1)
  {
    if (t < 0.5) 
    { 
      k = 0; 
      ratio = 0.0;
    }
    else 
    { 
      k = N-2; 
      ratio = 1.0;
    }
  }

  // Den neuen Parameter u ableiten.

  double u = s[k] + ratio*(s[k+1] - s[k]);
  return u;
}

// ****************************************************************************
// Returns the tangent at the curve parameter t.
// ****************************************************************************

Point3D Spline3D::getTangent(double t)
{
  const double DELTA = 0.000001;

  return (getPoint(t+0.5*DELTA) - getPoint(t-0.5*DELTA)) / DELTA;
}

// ****************************************************************************
// Calc. the intersecting point between this curve and the plane that is defined
// by the given plane point and plane normal.
// Returned is the t of the intersection point that is closest to planePoint.
// If there is no intersection, the return value is -1.
// ****************************************************************************

double Spline3D::getIntersection(const Point3D planePoint, const Point3D planeNormal, double tMin, double tMax)
{
  const int NUM_VERTICES = 16;
  double dotProduct[NUM_VERTICES];
  double vertexT[NUM_VERTICES];
  Point3D Q[NUM_VERTICES];
  Point3D C;
  double d;
  int i, k;

  if (numPoints < 2) { return -1.0; }

  // NUM_VERTICES aequidistante Punkte auf der Kurve zwischen tMin und tMax holen
  // und bestimmen, ob sie vor hinter der Ebene liegen.

  for (i=0; i < NUM_VERTICES; i++)
  {
    vertexT[i] = tMin + (tMax - tMin)*(double)i / (double)(NUM_VERTICES-1);
    Q[i] = getPoint(vertexT[i]);
    
    // Skalarprodukt von planeNormal und (Q[i]-planePoint) bilden
    dotProduct[i] = scalarProduct(planeNormal, Q[i] - planePoint);
  }

  // Ueberpruefen, ob zwei benachbarte Stuetzpunkte auf unterschiedl. Seiten der Flaeche sind.

  double minDist = -1.0;
  k = -1;

  for (i=0; i < NUM_VERTICES-1; i++)
  {
    // Liegen die Punkte i und i+1 auf unterschiedlichen Seiten der Ebene ?
    if (((dotProduct[i] >= 0.0) && (dotProduct[i+1] <= 0.0)) ||
        ((dotProduct[i] <= 0.0) && (dotProduct[i+1] >= 0.0)))
    {
      d = 0.0;
      C = (Q[i] + Q[i+1])*0.5 - planePoint;
      d = C.x*C.x + C.y*C.y + C.z*C.z;

      if ((minDist < 0.0) || (d < minDist))
      {
        k = i;
        minDist = d;
      }
    }
  }

  // Wurde ein Schnittpunkt mit der Flaeche gefunden ?

  if (k != -1)
  {
    // Die Unterteilung der Kurve war nun fein genug (< 1%).
    if (tMax - tMin < 0.0001)
    {
      return 0.5*(vertexT[k] + vertexT[k+1]);
    }
    else
    // Den ermittelten Linienabschnitt weiter unterteilen.
    {
      return getIntersection(planePoint, planeNormal, vertexT[k], vertexT[k+1]);
    }
  }

  return -1;
}

// ----------------------------------------------------------------------------
// 2D line strip.
// ----------------------------------------------------------------------------

// ****************************************************************************
// ****************************************************************************

LineStrip2D::LineStrip2D()
{
  reset(0);
}

// ****************************************************************************
// ****************************************************************************

LineStrip2D::LineStrip2D(int newNumPoints, Point2D *points)
{
  setPoints(newNumPoints, points);
}

// ****************************************************************************
// ****************************************************************************

void LineStrip2D::reset(int newNumPoints)
{
  numPoints = newNumPoints;
  if (numPoints > MAX_SPLINE_POINTS) { numPoints = MAX_SPLINE_POINTS; }
  if (numPoints < 0) { numPoints = 0; }

  for (int i=0; i < numPoints; i++) { P[i].set(0.0, 0.0); }
  pointsChanged = true;
}

// ****************************************************************************
// ****************************************************************************

void LineStrip2D::setPoints(int newNumPoints, const Point2D *points)
{
  numPoints = newNumPoints;
  if (numPoints > MAX_SPLINE_POINTS) { numPoints = MAX_SPLINE_POINTS; }
  if ((numPoints < 0) || (points == NULL)) { numPoints = 0; }

  for (int i=0; i < numPoints; i++) { P[i] = points[i]; }
  pointsChanged = true;
}

// ****************************************************************************
// ****************************************************************************

void LineStrip2D::setPoint(int index, Point2D point)
{
  if ((index < 0) || (index >= numPoints)) { return; }
  P[index] = point;
  pointsChanged = true;
}

// ****************************************************************************
// Add a new control point to the end of the list.
// ****************************************************************************

void LineStrip2D::addPoint(Point2D point)
{
  if (numPoints >= MAX_SPLINE_POINTS) { return; }

  P[numPoints++] = point;
  pointsChanged = true;
}

// ****************************************************************************
// Deletes a control point from the end of the list.
// ****************************************************************************

void LineStrip2D::delPoint()
{
  if (numPoints > 0) { numPoints--; }
}

// ****************************************************************************
// ****************************************************************************

Point2D LineStrip2D::getControlPoint(int index)
{
  if ((index < 0) || (index >= numPoints)) 
    { return Point2D(0.0, 0.0); }
  else
    { return P[index]; }
}

// ****************************************************************************
// Returns the point at the curve position t in [0..1].
// ****************************************************************************

Point2D LineStrip2D::getPoint(double t)
{
  if (pointsChanged) { calculateParams(); }
  
  if (numPoints < 1) { return Point2D(0.0, 0.0); }
  if (numPoints == 1) { return P[0]; }
  if (t < 0.0) { t = 0.0; }
  if (t > 1.0) { t = 1.0; }

  const double EPSILON = 0.000001;
  double length;
  double ratio = 0.0;
  int i;
  int k = -1;

  for (i=0; i < numPoints-1; i++)
  {
    if ((t >= pos[i]-EPSILON) && (t <= pos[i+1]+EPSILON))
    {
      k = i;
      length = pos[i+1] - pos[i];
      if (length < EPSILON) { length = EPSILON; }
      ratio = (t - pos[i]) / length;
    }
  }

  if (k == -1) { return Point2D(0.0, 0.0); }
  return P[k] + ratio*(P[k+1]-P[k]);
}

// ****************************************************************************
// Returns the function value y when the line strip is interpreted as
// a function y = f(x).
// ****************************************************************************

double LineStrip2D::getFunctionValue(double x)
{
  if (numPoints < 1) { return 0.0; }
  if (numPoints == 1) { return P[0].y; }

  const double EPSILON = 0.000001;
  double length;
  double result = 0.0;
  int i;

  for (i=0; i < numPoints-1; i++)
  {
    if ((x >= P[i].x-EPSILON) && (x <= P[i+1].x+EPSILON))
    {
      length = P[i+1].x - P[i].x;
      if (length < EPSILON) { length = EPSILON; }
      result = P[i].y + (P[i+1].y - P[i].y)*(x-P[i].x) / length;
    }
  }

  return result;
}

// ****************************************************************************
// Returns the tanget at the curve parameter t in [0..1].
// ****************************************************************************

Point2D LineStrip2D::getTangent(double t)
{
  const double DELTA = 0.000001;
  return (getPoint(t+0.5*DELTA) - getPoint(t-0.5*DELTA)) / DELTA;
}

// ****************************************************************************
// Returns the curve parameter t in [0..1] for the given control point.
// ****************************************************************************

double LineStrip2D::getCurveParam(int index)
{
  if (pointsChanged) { calculateParams(); }
  if (index < 0) { index = 0; }
  if (index > numPoints-1) { index = numPoints-1; }

  return pos[index];  
}

// ****************************************************************************
// Returns the intersection of the line Q+t*v with this line strip.
// Returns true, when the line strip was really intersected.
// When there are multiple potential intersections, the one closest to Q will
// be returned.
// ****************************************************************************

bool LineStrip2D::getClosestIntersection(const Point2D Q, const Point2D v, double &t, Point2D &intersection)
{
  const double EPSILON = 0.000001;
  int i;

  // Einen Normaleneinheitsvektor bilden, der senkrecht (90 Grad nach links
  // gedreht) auf v steht.

  Point2D n(-v.y, v.x);
  n.normalize();

  // Die zwei Punkte P_left und P_right berechnen, die im Abstand EPSILON
  // links bzw. rechts der Gerade Q+t*v auf der Hoehe von Q liegen.

  Point2D P_left  = Q + EPSILON*n;
  Point2D P_right = Q - EPSILON*n;

  // Alle Punkte des Linienzuges durchlaufen und jeweils ueberpruefen,
  // ob sie links, rechts, oder innerhalb der EPSILON-Umgebung der
  // Schnittlinie liegen.

  bool ok = false;              // Es gab bisher keinen Schnittpunkt
  t = 1000000.0;                // Schnittpunkt waere sehr, sehr weit weg
  intersection.set(0.0, 0.0);   // Vorbelegung

  Point2D w, R;
  int section = 0;
  int oldSection = 0;
  double s, d;
  double denominator;

  for (i=0; i < numPoints; i++)
  {
    section = 0;    // In der EPSILON-Umgebung der Schnittlinie

    w = P[i] - P_left;
    d = w.x*v.y - w.y*v.x;
    if (d < 0.0) { section = -1; }

    w = P[i] - P_right;
    d = w.x*v.y - w.y*v.x;
    if (d > 0.0) { section = 1; }

    // Ein Schnitt ist moeglich.

    if ((i > 0) &&
        (((section >= 0) && (oldSection <= 0)) ||
         ((section <= 0) && (oldSection >= 0))))
    {
      R = Q - P[i-1];
      w = P[i] - P[i-1];

      denominator = v.x*w.y - v.y*w.x;
      if (denominator != 0.0)
      {
        s = (v.x*R.y - R.x*v.y) / denominator;
        if ((s >= -EPSILON) && (s <= 1.0+EPSILON))
        {
          d = (w.x*R.y - R.x*w.y) / denominator;

          // Ist der neue Schnittpunkt dichter an P als ein evtl. alter ?
          if (fabs(d) < fabs(t))
          {
            t = d;
            intersection = Q + t*v;
            ok = true;
          }
        }
      }
    }
    
    oldSection = section;
  }

  return ok;
}


// ****************************************************************************
// Returns the intersection of the line Q+t*v with this line strip.
// Returns true, when the line strip was really intersected.
// When there are multiple potential intersections, the first one will be
// returned (in the order of the line strip control points).
// ****************************************************************************

bool LineStrip2D::getFirstIntersection(const Point2D Q, const Point2D v, double &t, Point2D &intersection)
{
  const double EPSILON = 0.000001;
  int i;

  // Einen Normaleneinheitsvektor bilden, der senkrecht (90 Grad nach links
  // gedreht) auf v steht.

  Point2D n(-v.y, v.x);
  n.normalize();

  // Die zwei Punkte P_left und P_right berechnen, die im Abstand EPSILON
  // links bzw. rechts der Gerade Q+t*v auf der Hoehe von Q liegen.

  Point2D P_left  = Q + EPSILON*n;
  Point2D P_right = Q - EPSILON*n;

  // Alle Punkte des Linienzuges durchlaufen und jeweils ueberpruefen,
  // ob sie links, rechts, oder innerhalb der EPSILON-Umgebung der
  // Schnittlinie liegen.

  bool ok = false;              // Es gab bisher keinen Schnittpunkt
  t = 1000000.0;                // Schnittpunkt waere sehr, sehr weit weg
  intersection.set(0.0, 0.0);   // Vorbelegung

  Point2D w, R;
  int section = 0;
  int oldSection = 0;
  double s, d;
  double denominator;

  for (i=0; (i < numPoints) && (ok == false); i++)
  {
    section = 0;    // In der EPSILON-Umgebung der Schnittlinie

    w = P[i] - P_left;
    d = w.x*v.y - w.y*v.x;
    if (d < 0.0) { section = -1; }

    w = P[i] - P_right;
    d = w.x*v.y - w.y*v.x;
    if (d > 0.0) { section = 1; }

    // Ein Schnitt ist moeglich.

    if ((i > 0) &&
        (((section >= 0) && (oldSection <= 0)) ||
         ((section <= 0) && (oldSection >= 0))))
    {
      R = Q - P[i-1];
      w = P[i] - P[i-1];

      denominator = v.x*w.y - v.y*w.x;
      if (denominator != 0.0)
      {
        s = (v.x*R.y - R.x*v.y) / denominator;
        if ((s >= -EPSILON) && (s <= 1.0+EPSILON))
        {
          t = (w.x*R.y - R.x*w.y) / denominator;
          intersection = Q + t*v;
          ok = true;
        }
      }
    }
    
    oldSection = section;
  }

  return ok;
}

// ****************************************************************************
// Returns the intersection of the line Q+t*v with this line strip.
// Returns true, when the line strip was really intersected.
// When there are multiple potential intersections, the one will be
// returned that is expected for an intersection with one of the vocal tract 
// contours.
// ****************************************************************************

bool LineStrip2D::getSpecialIntersection(const Point2D Q, const Point2D v, double &t, Point2D &intersection)
{
  const double EPSILON = 0.000001;
  int i;

  // Einen Normaleneinheitsvektor bilden, der senkrecht (90 Grad nach links
  // gedreht) auf v steht.

  Point2D n(-v.y, v.x);
  n.normalize();

  // Die zwei Punkte P_left und P_right berechnen, die im Abstand EPSILON
  // links bzw. rechts der Gerade Q+t*v auf der Hoehe von Q liegen.

  Point2D P_left  = Q + EPSILON*n;
  Point2D P_right = Q - EPSILON*n;

  // Alle Punkte des Linienzuges durchlaufen und jeweils ueberpruefen,
  // ob sie links, rechts, oder innerhalb der EPSILON-Umgebung der
  // Schnittlinie liegen.

  bool ok = false;              // Es gab bisher keinen Schnittpunkt
  t = 1000000.0;                // Schnittpunkt waere sehr, sehr weit weg
  intersection.set(0.0, 0.0);   // Vorbelegung

  Point2D w, R;
  int section = 0;
  int oldSection = 0;
  double s, d;
  double denominator;

  for (i=0; i < numPoints; i++)
  {
    section = 0;    // In der EPSILON-Umgebung der Schnittlinie

    w = P[i] - P_left;
    d = w.x*v.y - w.y*v.x;
    if (d < 0.0) { section = -1; }

    w = P[i] - P_right;
    d = w.x*v.y - w.y*v.x;
    if (d > 0.0) { section = 1; }

    // Ein Schnitt ist moeglich.

    if ((i > 0) && (section >= 0) && (oldSection <= 0)) 
    {
      R = Q - P[i-1];
      w = P[i] - P[i-1];

      denominator = v.x*w.y - v.y*w.x;
      if (denominator != 0.0)
      {
        s = (v.x*R.y - R.x*v.y) / denominator;
        if ((s >= -EPSILON) && (s <= 1.0+EPSILON))
        {
          d = (w.x*R.y - R.x*w.y) / denominator;

          if ((ok == false) || ((d >= 0.0) && (d < t)) || ((d <= 0.0) && (d > t)))
          {
            t = d;
            intersection = Q + t*v;
            ok = true;
          }
        }
      }
    }
    
    oldSection = section;
  }

  return ok;
}

// ****************************************************************************
// Calculates the positions of the individual control points.
// ****************************************************************************

void LineStrip2D::calculateParams()
{
  int i;
  if (numPoints < 1) { return; }

  // Der Parameter pos[] speichert zu jedem Kontrollpunkt seine Position
  // relativ zur Gesamtlaenge der Kurve.

  pos[0] = 0.0;
  for (i=1; i < numPoints; i++)
  {
    pos[i] = pos[i-1] + (P[i] - P[i-1]).magnitude();
  }

  double length = pos[numPoints-1];
  if (length > 0.0)
  {
    for (i=1; i < numPoints; i++) { pos[i]/= length; }
  }
}

// ----------------------------------------------------------------------------
// 3D line strip.
// ----------------------------------------------------------------------------

// ****************************************************************************
// ****************************************************************************

LineStrip3D::LineStrip3D() : Spline3D() { }

// ****************************************************************************
// ****************************************************************************

LineStrip3D::LineStrip3D(int newNumPoints, Point3D *points) : 
  Spline3D(newNumPoints, points) { }

// ****************************************************************************
// Returns the point at the normalized position t.
// ****************************************************************************
 
Point3D LineStrip3D::getPoint(double t)
{
  if (pointsChanged) { calculateParams(); }
  
  if (numPoints < 1) { return Point3D(0.0, 0.0, 0.0); }
  if (numPoints == 1) { return P[0]; }
  if (t < 0.0) { t = 0.0; }
  if (t > 1.0) { t = 1.0; }

  const double EPSILON = 0.000001;
  double length;
  double ratio = 0.0;
  int i;
  int k = -1;

  for (i=0; i < numPoints-1; i++)
  {
    if ((t >= pos[i]-EPSILON) && (t <= pos[i+1]+EPSILON))
    {
      k = i;
      length = pos[i+1] - pos[i];
      if (length < EPSILON) { length = EPSILON; }
      ratio = (t - pos[i]) / length;
    }
  }

  if (k == -1) { return Point3D(0.0, 0.0, 0.0); }
  return P[k] + ratio*(P[k+1]-P[k]);
}

// ****************************************************************************
// Returns the normalized curve parameter for the point with the given index.
// ****************************************************************************

double LineStrip3D::getCurveParam(int index)
{
  if (pointsChanged) { calculateParams(); }
  if (index < 0) { index = 0; }
  if (index > numPoints-1) { index = numPoints-1; }

  return pos[index];  
}

// ****************************************************************************
// ****************************************************************************

double LineStrip3D::getIntersection(const Point3D planePoint, const Point3D planeNormal)
{
  if (pointsChanged) { calculateParams(); }
  if (numPoints < 2) { return 0.0; }

  const double EPSILON = 0.000001;
  int i;
  Point3D P0, P1, Q, v;
  double length;
  double d[2];
  double result = 0.0;
  double t;
  double denominator;
  double dist;
  double minDist = 1000000.0;


  for (i=0; i < numPoints-1; i++)
  {
    // Das erste Liniensegment wird ein kleines Stueck nach links verlaengert.

    P0 = P[i]; 

    if (i == 0)
    {
      v = P[i+1] - P[i];
      length = v.magnitude();
      if (length > 0.0)
      {
        v/= length;
        P0-= v*EPSILON;
      }
    }

    // Das letzte Liniensegment wird ein kleines Stueck nach rechts verlaengert.

    P1 = P[i+1];

    if (i+1 == numPoints-1)
    {
      v = P[i+1] - P[i];
      length = v.magnitude();
      if (length > 0.0)
      {
        v/= length;
        P1+= v*EPSILON;
      }
    }

    // Liegen P0 und P1 auf unterschiedlichen Seiten der Ebene.
    d[0] = scalarProduct(P0 - planePoint, planeNormal);
    d[1] = scalarProduct(P1 - planePoint, planeNormal);

    if (((d[0] <= 0.0) && (d[1] >= 0.0)) || ((d[1] <= 0.0) && (d[0] >= 0.0)))
    {
      v = P1 - P0;
      denominator = scalarProduct(v, planeNormal);
      if (denominator != 0.0)
      {
        t = -scalarProduct(P0-planePoint, planeNormal) / denominator;
        Q = P0 + t*v;
        dist = (Q - planePoint).magnitude();
        if (dist < minDist)
        {
          result = pos[i] + t*(pos[i+1] - pos[i]);
          minDist = dist;
        }
      }
    }
  }

  return result;
}

// ****************************************************************************
// Calculate the positions of the individual control points.
// ****************************************************************************

void LineStrip3D::calculateParams()
{
  int i;
  if (numPoints < 1) { return; }

  // Der Parameter pos[] speichert zu jedem Kontrollpunkt seine Position
  // relativ zur Gesamtlaenge der Kurve.

  pos[0] = 0.0;
  for (i=1; i < numPoints; i++)
  {
    pos[i] = pos[i-1] + (P[i] - P[i-1]).magnitude();
  }

  double length = pos[numPoints-1];
  if (length > 0.0)
  {
    for (i=1; i < numPoints; i++) { pos[i]/= length; }
  }
}

// ----------------------------------------------------------------------------
// Bezier curve.
// ----------------------------------------------------------------------------

// ****************************************************************************
// ****************************************************************************

BezierCurve3D::BezierCurve3D() : Spline3D() { }

// ****************************************************************************
// ****************************************************************************

BezierCurve3D::BezierCurve3D(int newNumPoints, const Point3D *points) : 
  Spline3D(newNumPoints, points) { }

// ****************************************************************************
// ****************************************************************************

BezierCurve3D::BezierCurve3D(int newNumPoints, const Point3D *points, const double *weights) : 
  Spline3D(newNumPoints, points, weights) { }

// ****************************************************************************
// Returns a point on the curve with 0 <= t <= 1.
// ****************************************************************************

Point3D BezierCurve3D::getPoint(double t)
{
  if (pointsChanged) { calculateCoeff(); }

  Point3D Q(0.0, 0.0, 0.0);
  if (numPoints < 2) { return Q; }

  int i;
  int n = numPoints-1;
  double f = 1.0;       // f = t^i
  double denominator = 0.0;

  for (i=0; i <= n; i++)
  {
    Q+= A[i]*f;
    denominator+= B[i]*f;
    f*= t;
  }

  Q/= denominator;
  return Q;
}

// ****************************************************************************
// Calc. the coefficients in front of the powers of t of the Bernstein polynoms
// B_(k,n)(t) = (n ueber k)*t^k*(1-t)^(n-k) 
// The order of the polynom is always n.
// ****************************************************************************

void BezierCurve3D::getBernsteinCoeff(int k, int n, double *coeff)
{
  int i, j;
  double newCoeff[MAX_SPLINE_POINTS];

  // Alle Koeffizienten sind erstmal Null.
  for (i=0; i <= n; i++) { coeff[i] = 0.0; }

  // Der Term t^k macht genau den Koeffizienten coeff[k] zu Eins.
  coeff[k] = 1.0;

  // (n-k) mal den Term (1-t) dazu multiplizieren.

  for (i=1; i <= n-k; i++)
  {
    // Die neuen Koeffizienten sind erstmal alle Null
    for (j=0; j <= n; j++) { newCoeff[j] = 0.0; }

    // Jeder alte Koeffizient addiert durch die Mult. mit (1-t) zwei.
    // neue Koeff. dazu

    for (j=0; j <= n; j++)
    {
      newCoeff[j]+= coeff[j];
      if (j+1 < MAX_SPLINE_POINTS) { newCoeff[j+1]-= coeff[j]; }
    }

    // Die neuen Koeffizienten werden zu den alten.
    for (j=0; j <= n; j++) { coeff[j] = newCoeff[j]; }
  }

  // Den Vorfaktor (n ueber k) berechnen und zu allen Koeff. multiplizieren.

  double numerator = 1.0;
  double denominator = 1.0;
  
  // Binomialkoeffizient (n ueber k)
  for (j=2; j <= n; j++)   { numerator*= (double)j; }
  for (j=2; j <= k; j++)   { denominator*= (double)j; }
  for (j=2; j <= n-k; j++) { denominator*= (double)j; }
  double f = numerator / denominator;

  for (j=0; j <= n; j++) { coeff[j]*= f; }
}

// ****************************************************************************
// Calc. the coefficients A[] and B[] of the polynoms in the numerator and
// denominator for the efficient calculation of the points on the curve.
// ****************************************************************************

void BezierCurve3D::calculateCoeff()
{
  double coeff[MAX_SPLINE_POINTS];
  int i, j;

  int n = numPoints-1;    // Polynomgrad und obere Summationsgrenze

  // Die Koeffizienten mit Null initialisieren.
  for (j=0; j <= n; j++)
  {
    A[j].x = 0.0;       // Zaehlerkoeffizienten
    A[j].y = 0.0;
    A[j].z = 0.0;
    B[j]   = 0.0;       // Nennerkoeffizienten
  }

  // Die Summen ueber und unter dem Bruchstrich durchlaufen.

  for (i=0; i <= n; i++)
  {
    getBernsteinCoeff(i, n, coeff);

    // Die Bernsteinkoeffizienten gewichtet aufaddieren.
    for (j=0; j <= n; j++)
    {
      A[j]+= w[i]*coeff[j]*P[i];
      B[j]+= w[i]*coeff[j];
    }
  }

  pointsChanged = false;
}

// ****************************************************************************

