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

#ifndef __SURFACE_H__
#define __SURFACE_H__

#include "Geometry.h"
#include <fstream>
#include <string>

using namespace std;

// ****************************************************************************
/// @brief Representation of a surface of the vocal tract model. 
///
/// The surface consists of a 2D-array of vertices with two triangles between 
/// each four neighboring vertices.
// ****************************************************************************

class Surface
{

  // **************************************************************************
  // Public data.
  // **************************************************************************

public:
  /// Maximal number of triangles sharing the same vertex
  static const int NUM_ASSOCIATED_TRIANGLES = 6;  

  /// Max. number of tiles in horizontal direction
  static const int MAX_TILES_X = 15;  
  /// Max. number of tiles in vertical direction
  static const int MAX_TILES_Y = 15;  
  /// Max. number of triangles per tile
  static const int MAX_TILE_TRIANGLES = 150000 / (MAX_TILES_X*MAX_TILES_Y);

  static const double STANDARD_TILE_WIDTH;    ///< The optimal width of a tile
  static const double STANDARD_TILE_HEIGHT;   ///< The optimal height of a tile

  /// The default angle that separates between smooth shading and an edge
  static const double STANDARD_CREASE_ANGLE_DEGREE;

  // ****************************************************************
  /// A vertex of the surface.
  // ****************************************************************
  
  struct Vertex
  {
    Point3D coord;
    int rib;
    int ribPoint;
    int numAssociates;
    int associatedTriangle[NUM_ASSOCIATED_TRIANGLES];
    int associatedCorner[NUM_ASSOCIATED_TRIANGLES];
    int reserved;         // Vertex position in relation to a line: -1=left, +1=right, 0=on the line
	  bool wasTested;       // Was the vertext position in relation to an intersection line tested ?
  };

  // ****************************************************************
  /// A triangle defined by three vertices.
  // ****************************************************************
  
  struct Triangle
  {
    int vertex[3];            ///< Indices of the 3 vertices.
    int edge[3];              ///< Indices of the 3 edges.
    Point3D cornerNormal[3];  ///< Plane normals at the vertex positions.
    Point3D planeNormal;      ///< Normal of the triangle.
    double area;              ///< Area of the triangle.
    double distance;          ///< The minimal z-value of all 3 vertices.
  };

  // ****************************************************************
  /// An edge defined by two vertices.
  // ****************************************************************
 
  struct Edge
  {
    int vertex[2];          ///< Indices of the two vertices.
    bool isIntersected;     ///< Was the edge intersected by the intersecting plane?
	  bool wasTested; 		    ///< Was the edge already tested for an intersection?
    Point2D intersection;   ///< Projection of the intersection point on the intersecting plane.
  };

  // ****************************************************************
  /// @brief A tile represents a rectangular part of the bounding box 
  /// of the surface in the xy-plane.
  ///
  /// A tile contains the indices of all triangles that overlap the
  /// tile region to allow faster intersection point calculations.
  // ****************************************************************

  struct Tile
  {							
		int numTriangles;				///< Number of triangles overlapping the tile.
		int triangle[MAX_TILE_TRIANGLES];	///< Indices of the triangles overplapping the tile.
  };
		
  // ****************************************************************

  int numRibs;        ///< Number of vertices in the first dimension.
  int numRibPoints;   ///< Number of vertices in the second dimension.
  int numTriangles;   ///< Number of triangles.
  int numVertices;    ///< Number of vertices.
  int numEdges;       ///< Number of edges.

  Vertex *vertex;     ///< 1D array of vertices.
  Triangle *triangle; ///< Array of triangles.
  Edge *edge;         ///< Array of edges.
  int *sequence;      ///< Order, in which the triangles must be painted with the painters algorithm.

  // Information about the tiles ************************************

  Tile tile[MAX_TILES_X][MAX_TILES_Y];    ///< 2D array of tiles
  double leftBorder;    ///< Left border of the surface bounding box.
  double rightBorder;   ///< Right border of the surface bounding box.
  double topBorder;     ///< Top border of the surface bounding box.
  double bottomBorder;  ///< Bottom border of the surface bounding box.

  double tileWidth;     ///< With of a tile.
  double tileHeight;    ///< Height of a tile.
  int numTilesX;        ///< Number of tiles in x-direction
  int numTilesY;        ///< Number of tiles in y-direction

  /// The angle that separates between smooth shading and an edge
  double creaseAngle_deg; 

  // **************************************************************************
  // Public functions.
  // **************************************************************************

public:

  Surface();
  Surface(int ribs, int ribPoints);
  ~Surface();
  void init(int ribs, int ribPoints);
  void clear();
  
  /// Set the coordinates of the given vertex.
  inline void setVertex(int rib, int ribPoint, double x, double y, double z)
  {
    int adr = rib*numRibPoints + ribPoint;
    vertex[adr].coord.x = x;
    vertex[adr].coord.y = y;
    vertex[adr].coord.z = z;
  }

  /// Set the coordinates of the given vertex.
  inline void setVertex(int rib, int ribPoint, const Point3D &P)
  {
    vertex[rib*numRibPoints + ribPoint].coord = P;
  }

  /// Get the coordinates of the given vertex.
  inline Point3D getVertex(int rib, int ribPoint)
  {
    return vertex[rib*numRibPoints + ribPoint].coord;
  }

  /// Get the coordinates of the given vertex.
  inline void getVertex(int rib, int ribPoint, double &x, double &y, double &z)
  {
    int adr = rib*numRibPoints + ribPoint;
    x = vertex[adr].coord.x;
    y = vertex[adr].coord.y;
    z = vertex[adr].coord.z;
  }

  /// Get the index of the given vertex.
  inline int getVertexIndex(int rib, int ribPoint) { return rib*numRibPoints + ribPoint; } 

  // Set the normal of the given vertex.
  void setNormal(int rib, int ribPoint, Point3D normal);
  // Get the normal of the given vertex.
  Point3D getNormal(int rib, int ribPoint);

  void swapTriangleOrientation();
  void calculateNormals();
  void flipNormals();
  void calculatePaintSequence(double matrix[16]);
  void reversePaintSequence();

  bool saveAsObjFile(const string &fileName);

  // ****************************************************************
  // New !
  // ****************************************************************

  // Assign the triangles to the tiles.
  void prepareIntersections();
  
  // Prepare the intersection for an individual intersection line.
  void prepareIntersection(Point2D Q, Point2D v);

  // Returns a list with potentially interesected triangles. This
  // function must be called after prepareIntersections().
  bool getTriangleList(int *indexList, int &numEntries, int MAX_ENTRIES);

  // Returns the intersection data for a single triangle.
  bool getTriangleIntersection(int index, Point2D &P0, Point2D &P1, Point2D &n);

  void appendToFile(std::ofstream &file);
  void readFromFile(std::ifstream &file, bool initialize);

  // **************************************************************************
  // Private functions and data.
  // **************************************************************************

private:
  Point2D linePoint;        ///< Origin of the intersecting plane/line (in the xy-plane).
  Point2D leftLinePoint;    ///< The line origin moved to the left (with resprect to the line) by a tiny amount.
  Point2D rightLinePoint;   ///< The line origin moved to the right (with resprect to the line) by a tiny amount.
  Point2D lineVector;       ///< Normalized vector specifying the direction of the intersecting line.

  void quickSort(int firstIndex, int lastIndex);
  bool getEdgeIntersection(int edgeIndex);
};

// ****************************************************************************

#endif
