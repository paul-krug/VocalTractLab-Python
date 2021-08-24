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

#ifndef __GEOMETRY_H__
#define __GEOMETRY_H__

#include <cmath>

class Point2D;
class Point3D;
class Vector2D;
class Line2D;
class Line3D;
class Circle;
class Ellipse2D;

// ****************************************************************************
// A 2D point.
// ****************************************************************************

class Point2D
{
public:
  Point2D() : x(0), y(0) { }
  Point2D(double X, double Y) : x(X), y(Y) { }
  void set(double X, double Y) { x = X; y = Y; }

  bool isRightFrom(Vector2D& V);
  bool isLeftFrom(Vector2D& V);
  Point3D toPoint3D();
  void leanOn(Vector2D V, Point2D A);
  Point2D normalize();
  Point2D turnRight();
  Point2D turnLeft();
  void turn(double angle);
  double getDistanceFrom(Vector2D V);
  
  double magnitude() { return sqrt(x*x + y*y); }
  double squareMagnitude() { return x*x + y*y; }
  
  Point2D& operator= (const Point2D &Q) { x = Q.x; y = Q.y; return *this; }
  Point2D& operator+=(const Point2D &Q) { x+= Q.x; y+= Q.y; return *this; }
  Point2D& operator-=(const Point2D &Q) { x-= Q.x; y-= Q.y; return *this; }
  Point2D& operator*=(double d)  { x*= d;   y*= d;   return *this; }
  Point2D& operator/=(double d)  { x/= d;   y/= d;   return *this; }

  // ****************************************************************
  
  double x;
  double y;
};

inline bool operator!=(const Point2D &P, const Point2D &Q) { return (P.x != Q.x) || (P.y != Q.y); }
inline bool operator==(const Point2D &P, const Point2D &Q) { return (P.x == Q.x) && (P.y == Q.y); }
inline Point2D operator+(const Point2D &P, const Point2D &Q) { return Point2D(P.x+Q.x, P.y+Q.y); }
inline Point2D operator-(const Point2D &P, const Point2D &Q) { return Point2D(P.x-Q.x, P.y-Q.y); }
inline Point2D operator*(const Point2D &P, double d) { return Point2D(P.x*d, P.y*d); }
inline Point2D operator*(double d, const Point2D &P) { return Point2D(P.x*d, P.y*d); }
inline Point2D operator/(const Point2D &P, double d) { return Point2D(P.x/d, P.y/d); }
inline Point2D operator-(const Point2D &P)           { return Point2D(-P.x, -P.y); }

double scalarProduct(const Point2D &P, const Point2D &Q);

// ****************************************************************************
// A 3D point.
// ****************************************************************************

class Point3D
{
public:
  Point3D() : x(0), y(0), z(0) { }
  Point3D(double X, double Y, double Z) : x(X), y(Y), z(Z) { }
  void set(double X, double Y, double Z) { x = X; y = Y; z = Z; }
  
  bool    isRightFrom(Vector2D& V);
  Point2D toPoint2D();
  Point3D normalize();
  double  magnitude() { return sqrt(x*x + y*y + z*z); }

  Point3D& operator= (const Point3D &Q) { x = Q.x; y = Q.y; z = Q.z; return *this; }
  Point3D& operator+=(const Point3D &Q) { x+= Q.x; y+= Q.y; z+= Q.z; return *this; }
  Point3D& operator-=(const Point3D &Q) { x-= Q.x; y-= Q.y; z-= Q.z; return *this; }
  Point3D& operator*=(double d)  { x*= d;   y*= d;   z*= d;   return *this; }
  Point3D& operator/=(double d)  { x/= d;   y/= d;   z/= d;   return *this; }

  // ****************************************************************

  double x;
  double y;
  double z;
};

inline bool operator!=(const Point3D &P, const Point3D &Q) { return (P.x != Q.x) || (P.y != Q.y) || (P.z != Q.z); }
inline bool operator==(const Point3D &P, const Point3D &Q) { return (P.x == Q.x) && (P.y == Q.y) && (P.z == Q.z); }
inline Point3D operator+(const Point3D &P, const Point3D &Q) { return Point3D(P.x+Q.x, P.y+Q.y, P.z+Q.z); }
inline Point3D operator-(const Point3D &P, const Point3D &Q) { return Point3D(P.x-Q.x, P.y-Q.y, P.z-Q.z); }
inline Point3D operator*(const Point3D &P, double d) { return Point3D(P.x*d, P.y*d, P.z*d); }
inline Point3D operator*(double d, const Point3D &P) { return Point3D(P.x*d, P.y*d, P.z*d); }
inline Point3D operator/(const Point3D &P, double d) { return Point3D(P.x/d, P.y/d, P.z/d); }
inline Point3D operator-(const Point3D &P)           { return Point3D(-P.x, -P.y, -P.z); }

double  scalarProduct(const Point3D &P, const Point3D &Q);
Point3D crossProduct(const Point3D &P, const Point3D &Q);


// ****************************************************************************
// A 2D vector with a basis point and a vector.
// ****************************************************************************

class Vector2D
{
  public:
    Vector2D();
    Vector2D(Point2D _P, Point2D _v);

    void set(Point2D _P, Point2D _v);
    Point2D getIntersection(Vector2D V, double& t);
    Point2D getIntersection(Line2D L, double& t, bool& ok);
    Point3D getIntersection(Line3D L, double& t, bool& ok);
    void normalize();
    Point2D getPoint(double t);
    double getLength(double t);
    bool isNotNull();

    void operator=(Vector2D Q);
    
    // The basis point and the vector.
    
    Point2D P;
    Point2D v;
};


// ****************************************************************************
// A 2D line defined by two points.
// ****************************************************************************

class Line2D
{
  public:
    Line2D();
    Line2D(Point2D P0, Point2D P1);

    void set(Point2D P0, Point2D P1);
    Point2D getIntersection(Vector2D V, double& t, bool& ok);
    Point2D getIntersection(Line2D L, double& t, bool& ok);
    Point2D getPoint(double t);
    double getLength();
    bool encloses(Point2D Q);
    
    // The two points.
    Point2D P[2];
};

// ****************************************************************************
// A 3D line defined by two points.
// ****************************************************************************

class Line3D
{
  public:
    Line3D();
    Line3D(Point3D P0, Point3D P1);

    void set(Point3D P0, Point3D P1);
    Point3D getIntersection(Vector2D V, double& t, bool& ok);
    Point3D getPoint(double t);
    double getLength();
    
    // The two points.

    Point3D P[2];
};


// ****************************************************************************
// A circle.
// ****************************************************************************

class Circle
{
  public:
    Circle();
    Circle(Point2D C, double radius);

    void setValidArc(double angle0, double angle1);
    double getLength();
    double getLength(double a0, double a1);
    Point2D getPoint(double angle);
    Vector2D getNormal(double angle);
    Point2D getIntersection(Vector2D V, double &intersectionAngle, bool &ok);
    Point2D getIntersection(Line2D L, double &intersectionAngle, bool &ok);
    double getTangentContactAngle(Point2D H, bool uhrzeiger);
    double getBend();
    bool isIncluding(Point2D P);
    bool isIncluding(Line2D L);
    bool hasIntersectionWith(Line2D L);
    Line2D getCommonLeftTangentWith(Circle K);

    // **************************************************************

    Point2D M;
    double r;
    double arcAngle[2];
};

double getCircleTangent(Point2D H, Point2D C, double r, bool uhrzeiger);


// ****************************************************************************
// An ellipse. 
// ****************************************************************************

const long NUM_ELLIPSE_PROPS = 33;

class Ellipse2D
{
  public:
    Ellipse2D();
    Ellipse2D(Point2D _M, double _halfWidth, double _halfHeight);

    void set(Point2D _M, double _halfWidth, double _halfHeight);
    Point2D getPoint(double angle);
    double getAngle(double t);
    double getPerimeter();

    double halfWidth;
    double halfHeight;
    Point2D M;

    // **************************************************************

  private:
    double propAngle[NUM_ELLIPSE_PROPS];
    double propPosition[NUM_ELLIPSE_PROPS];
    double perimeter;           // der Umfang
};

double getEllipseTangent(Point2D H, Point2D C, double a, double b, bool uhrzeiger);

// ****************************************************************************

#endif
