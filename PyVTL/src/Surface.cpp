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

#include "Surface.h"
#include <cmath>
#include <iostream>
#include <fstream>

using namespace std;


// Borders for the tile wall.

const double Surface::STANDARD_TILE_WIDTH = 0.5;
const double Surface::STANDARD_TILE_HEIGHT = 0.5;

const double Surface::STANDARD_CREASE_ANGLE_DEGREE = 70.0;

// ****************************************************************************
/// @brief The constructor.
///
/// The function init() must be called to initialize the surface.
// ****************************************************************************

Surface::Surface()
{
  vertex   = NULL;
  triangle = NULL;
  edge     = NULL;
  sequence = NULL;

  creaseAngle_deg = STANDARD_CREASE_ANGLE_DEGREE;
  init(0, 0);
}

// ****************************************************************************
/// @brief The constructor.
///
/// Initialization of the surface with the given surface dimensions.
/// @param ribs Number of ribs.
/// @param ribPoints Number of vertices per rib.
// ****************************************************************************

Surface::Surface(int ribs, int ribPoints)
{
  vertex   = NULL;
  triangle = NULL;
  edge     = NULL;
  sequence = NULL;

  creaseAngle_deg = STANDARD_CREASE_ANGLE_DEGREE;
  init(ribs, ribPoints);
}

// ****************************************************************************
/// @brief The destructor.
// ****************************************************************************

Surface::~Surface()
{
  clear();
}

// ****************************************************************************
/// @brief This function never needs to be called explicitely.
// ****************************************************************************

void Surface::clear()
{
  if (vertex != NULL)   { delete [] vertex; }
  if (triangle != NULL) { delete [] triangle; }
  if (edge != NULL)     { delete [] edge; }
  if (sequence != NULL) { delete [] sequence; }

  numRibs      = 0;
  numRibPoints = 0;
  numTriangles = 0;
  numVertices  = 0;
  numEdges     = 0;
}

// ****************************************************************************
/// @brief Initialization of the surface with the given surface dimensions.
///
/// Memory on the heap for the triangles, edges and vertices is allocated here.
/// @param ribs Number of ribs.
/// @param ribPoints Number of vertices per rib.
// ****************************************************************************

void Surface::init(int ribs, int ribPoints)
{
  int i, j, k;

  clear();    // Alten Speicher wieder freigeben

  numRibs      = ribs;
  numRibPoints = ribPoints;

  if ((ribs == 0) || (ribPoints == 0))
  {
    numTriangles = 0;
    numVertices  = 0;
    numEdges     = 0;
    return;
  }
  else
  {
    numTriangles = 2*(numRibs-1)*(numRibPoints-1);
    numVertices  = numRibs*numRibPoints;
    numEdges     = 3*numRibs*numRibPoints - 2*numRibs - 2*numRibPoints + 1;
  }

  // Speicher für die Oberflächenpunkte und Normalen allokieren.
  
  vertex   = new Vertex[numVertices];
  triangle = new Triangle[numTriangles];
  edge     = new Edge[numEdges];
  sequence = new int[numTriangles];

  // Jedem Punkt den Index seiner Rippe und des Rippenpunktes zuweisen.

  int adr;
  for (i=0; i < numRibs; i++)
  {
    for (j=0; j < numRibPoints; j++)
    {
      adr = i*numRibPoints + j;
      vertex[adr].rib = i;
      vertex[adr].ribPoint = j;
    }
  }

  // Den Kanten ihre Endpunkte zuweisen.

  int e = 0;

  // Horizontale Kanten
  for (j=0; j < numRibPoints; j++)
  {
    for (i=0; i < numRibs-1; i++)
    {
      edge[e].vertex[0] = i*numRibPoints + j;
      edge[e].vertex[1] = (i+1)*numRibPoints + j;
      e++;
    }
  }

  // Vertikale Kanten
  for (i=0; i < numRibs; i++)
  {
    for (j=0; j < numRibPoints-1; j++)
    {
      edge[e].vertex[0] = i*numRibPoints + j;
      edge[e].vertex[1] = i*numRibPoints + j + 1;
      e++;
    }
  }

  // Diagonale Kanten
  for (j=0; j < numRibPoints-1; j++)
  {
    for (i=0; i < numRibs-1; i++)
    {
      edge[e].vertex[0] = i*numRibPoints + j;
      edge[e].vertex[1] = (i+1)*numRibPoints + j + 1;
      e++;
    }
  }

  // Den Dreiecken ihre Eckpunkte und Kanten zuweisen.

  int t = 0;
  int A = numRibPoints*(numRibs-1);
  int B = A + (numRibPoints-1)*numRibs;

  for (i=0; i < numRibs-1; i++)
  {
    for (j=0; j < numRibPoints-1; j++)
    {
      triangle[t].vertex[0] = i*numRibPoints + j;
      triangle[t].vertex[1] = (i+1)*numRibPoints + j + 1;
      triangle[t].vertex[2] = i*numRibPoints + j + 1;

      triangle[t].edge[0] = B + j*(numRibs-1) + i;      // Diagonale Kante
      triangle[t].edge[1] = (j+1)*(numRibs-1) + i;      // Horizontale Kante
      triangle[t].edge[2] = A + i*(numRibPoints-1) + j; // Vertikale Kante
      t++;

      triangle[t].vertex[0] = i*numRibPoints + j;
      triangle[t].vertex[1] = (i+1)*numRibPoints + j;
      triangle[t].vertex[2] = (i+1)*numRibPoints + j + 1;

      triangle[t].edge[0] = B + j*(numRibs-1) + i;          // Diagonale Kante;
      triangle[t].edge[1] = j*(numRibs-1) + i;              // Horizontale Kante
      triangle[t].edge[2] = A + (i+1)*(numRibPoints-1) + j; // Vertikale Kante
      t++;
    }
  }

  // Die an einen Vertex grenzenden Dreiecke ermitteln.

  for (i=0; i < numVertices; i++)
  {
    vertex[i].reserved = 0;
    vertex[i].numAssociates = 0;

    for (j=0; j < NUM_ASSOCIATED_TRIANGLES; j++)
    {
      vertex[i].associatedTriangle[j] = -1;
      vertex[i].associatedCorner[j] = -1;
    }
  }

  for (i=0; i < numTriangles; i++)
  {
    for (j=0; j < 3; j++)
    {
      k = triangle[i].vertex[j];
      vertex[k].associatedTriangle[ vertex[k].numAssociates ] = i;
      vertex[k].associatedCorner[ vertex[k].numAssociates ] = j;
      vertex[k].numAssociates++;
    }
  }

  // Abstand und Reihenfolge der Dreiecke vorbelegen.

  for (i=0; i < numTriangles; i++)
  {
  	triangle[i].distance = 0.0;
    sequence[i] = i;
  }

  for (i=0; i < numEdges; i++)
  {
    edge[i].isIntersected = false;
    edge[i].intersection.set(0.0, 0.0);
  }
}

// ****************************************************************************
/// @brief Swaps the orientation of the surface normals by changing the vertex
/// order of the triangle definitions.
/// 
/// The surface normals should be directed to the "outside" of the surface such
/// that shading can be performed correctly.
// ****************************************************************************

void Surface::swapTriangleOrientation()
{
  int i, j, k;

  for (i=0; i < numTriangles; i++)
  {
    k = triangle[i].vertex[0];
    triangle[i].vertex[0] = triangle[i].vertex[2];
    triangle[i].vertex[2] = k;
  }

  // Die an einen Vertex grenzenden Dreiecke neu ermitteln.

  for (i=0; i < numVertices; i++)
  {
    vertex[i].numAssociates = 0;

    for (j=0; j < NUM_ASSOCIATED_TRIANGLES; j++)
    {
      vertex[i].associatedTriangle[j] = -1;
      vertex[i].associatedCorner[j] = -1;
    }
  }

  for (i=0; i < numTriangles; i++)
  {
    for (j=0; j < 3; j++)
    {
      k = triangle[i].vertex[j];
      vertex[k].associatedTriangle[ vertex[k].numAssociates ] = i;
      vertex[k].associatedCorner[ vertex[k].numAssociates ] = j;
      vertex[k].numAssociates++;
    }
  }
}

// ****************************************************************************
/// @brief Calculates the normalized normal vectors of the triangles and at the 
/// vertices.
// ****************************************************************************

void Surface::calculateNormals()
{
  int i, j, k;

  // Die Normalen aller Flächen berechnen und auf die Eckpunkte übertragen.

  int p0, p1, p2;
  Point3D v0;
  Point3D v1;
  double cosAngle;
  double length;

  for (i=0; i < numTriangles; i++)
  {
    p0 = triangle[i].vertex[0];
    p1 = triangle[i].vertex[1];
    p2 = triangle[i].vertex[2];

    // Die Vektoren von p0->p1 und p0->p2
    v0 = vertex[p1].coord - vertex[p0].coord;
    v1 = vertex[p2].coord - vertex[p0].coord;

    // Normalenvektor der Fläche (Länge ist proportional zum Flächeninhalt)

    triangle[i].planeNormal = crossProduct(v0, v1);
    length = triangle[i].planeNormal.magnitude();
    triangle[i].area = 0.5*length;

    // Den Flächennormalenvektor normieren.
    triangle[i].planeNormal.normalize();

    // Die Eckpunktnormalen werden mit der Flächennormalen initialisiert (skaliert mit Flächeninhalt).
    for (k=0; k < 3; k++)
    {
      triangle[i].cornerNormal[k] = triangle[i].area*triangle[i].planeNormal;
    }
  }

  // Die Eckpunktnormalen jedes Punktes in jedem Dreieck berechnen.

  Point3D planeNormal0, planeNormal1;
  Point3D *cornerNormal0, *cornerNormal1;
  double area0, area1;
  double cosCreaseAngle = cos(creaseAngle_deg*3.1415/180.0);

  for (i=0; i < numVertices; i++)
  {
    // Jedes Zweierpärchen von Dreiecken an diesem Punkt überprüfen.
    for (j=0; j < vertex[i].numAssociates-1; j++)
    {
      planeNormal0 = triangle[ vertex[i].associatedTriangle[j] ].planeNormal;
      area0        = triangle[ vertex[i].associatedTriangle[j] ].area;

      for (k=j+1; k < vertex[i].numAssociates; k++)
      {
        planeNormal1 = triangle[ vertex[i].associatedTriangle[k] ].planeNormal;
        area1        = triangle[ vertex[i].associatedTriangle[k] ].area;

        cosAngle = scalarProduct(planeNormal0, planeNormal1);
        
        if (cosAngle > cosCreaseAngle)
        {
          // Die Flächennormalen der Flächen j und k gegenseitig auf die Punkte addieren.

          cornerNormal0 = &triangle[ vertex[i].associatedTriangle[j] ].cornerNormal[ vertex[i].associatedCorner[j] ];
          cornerNormal1 = &triangle[ vertex[i].associatedTriangle[k] ].cornerNormal[ vertex[i].associatedCorner[k] ];

          (*cornerNormal0)+= area1*planeNormal1;
          (*cornerNormal1)+= area0*planeNormal0;
        }
      }
    }
  }

  // Alle Eckpunktnormalen auf die Länge 1 bringen.

  for (i=0; i < numTriangles; i++)
  {
    for (k=0; k < 3; k++)
    {
      triangle[i].cornerNormal[k].normalize();
    }
  }
}

// ****************************************************************************
/// @brief Flip the normal vectors at all vertices inside out.
// ****************************************************************************

void Surface::flipNormals()
{
  int i, j;

  for (i=0; i < numRibs; i++)
  {
    for (j=0; j < numRibPoints; j++)
    {
      setNormal(i, j, -getNormal(i, j));
    }
  }
}

// ****************************************************************************
/// @brief Sets the normal vector for the given vertex to \a normal in all 
/// triangles having this vertex in common.
/// 
/// @param rib The rib index.
/// @param ribPoint Index of the rib point.
/// @param normal The normal vector.
// ****************************************************************************

void Surface::setNormal(int rib, int ribPoint, Point3D normal)
{
  int adr = rib*numRibPoints + ribPoint;

  for (int i=0; i < vertex[adr].numAssociates; i++)
  {
    triangle[ vertex[adr].associatedTriangle[i] ].cornerNormal[ vertex[adr].associatedCorner[i] ] = normal;
  }
}

// ****************************************************************************
/// @brief Returns the average normal for the given vertex.
///
/// @param rib The rib index.
/// @param ribPoint Index of the rib point.
// ****************************************************************************

Point3D Surface::getNormal(int rib, int ribPoint)
{
  Point3D normal(0, 0, 0);
  int adr = rib*numRibPoints + ribPoint;

  for (int i=0; i < vertex[adr].numAssociates; i++)
  {
    normal+= triangle[ vertex[adr].associatedTriangle[i] ].cornerNormal[ vertex[adr].associatedCorner[i] ];
  }

  normal.normalize();
  return normal;
}

// ****************************************************************************
/// @brief Calculates the order, in which the triangles must be painted
/// according to the painters algorithm and writes the results into \a sequence
/// (Triangles with the smallest z-coord. are painted first).
/// 
/// Therefore, the center points of all triangles are transformed with the matrix
/// \a matrix and sorted. Attention: The matrix elements must be passed in 
/// transposed form: Row 1:[m0, m4, m8, m12]; Row 2:[m1, m5, m9, m13]; ect.
/// @param matrix The transposed transformation (model/view) matrix.
// ****************************************************************************

void Surface::calculatePaintSequence(double matrix[16])
{
  int i;
  Point3D M;
  double z, w;
  int v0, v1, v2;

  for (i=0; i < numTriangles; i++)
  {
    // Den Mittelpunkt des Dreiecks berechnen.
    v0 = triangle[i].vertex[0];
    v1 = triangle[i].vertex[1];
    v2 = triangle[i].vertex[2];

    M = (vertex[v0].coord + vertex[v1].coord + vertex[v2].coord) / 3.0;

    // Den z- und den w-Wert des Transformierten Mittelpunktes berechnen.
    z = matrix[2]*M.x + matrix[6]*M.y + matrix[10]*M.z + matrix[14];
    w = matrix[3]*M.x + matrix[7]*M.y + matrix[11]*M.z + matrix[15];

    triangle[i].distance = z/w;
    sequence[i] = i;
  }

  quickSort(0, numTriangles-1);
}

// ****************************************************************************
/// @brief The triangles are sorted by their z-coordinate.
///
/// The triangle with the smallest z-coord. will be the first in the list and
/// must be painted first with the painters algorithm.
// ****************************************************************************

void Surface::quickSort(int firstIndex, int lastIndex)
{
  if ((numTriangles < 1) || (firstIndex > lastIndex)) { return; }

  int left = firstIndex;
  int right = lastIndex;
  double pivot = triangle[ sequence[(firstIndex + lastIndex)/2] ].distance;
  int temp;

  do
  {
    while (triangle[ sequence[left] ].distance < pivot) { left++; }
    while (triangle[ sequence[right] ].distance > pivot) { right--; }
    if (left <= right)
    {
      temp = sequence[left];
      sequence[left] = sequence[right];
      sequence[right] = temp;
      left++;
      right--;
    }
  } while (left <= right);

  if (firstIndex < right) { quickSort(firstIndex, right); }
  if (left < lastIndex)   { quickSort(left, lastIndex); }
}

// ****************************************************************************
/// @brief Reverses the paint sequence that was calculated in 
/// calculatePaintSequence().
///
/// When this function is called right after calculatePaintSequence(), the 
/// triangle with the greatest z-value will be first in the sequence[] array.
/// In contrast to OpenGL, positive z-coordinates are directed into the screen
/// in Direct3D, so that the first triangle to paint is the one with the
/// greatest z-coordinate.
// ****************************************************************************

void Surface::reversePaintSequence()
{
  int i, temp;
  
  for (i=0; i < numTriangles/2; i++)
  {
    temp = sequence[i];
    sequence[i] = sequence[numTriangles-1-i];
    sequence[numTriangles-1-i] = temp;
  }
}


// ****************************************************************************
/// Save the surface in the Wavefront .obj file format.
// ****************************************************************************

bool Surface::saveAsObjFile(const string &fileName)
{
  ofstream file(fileName);
  if (!file) 
  { 
    return false; 
  }

  int i;
  Point3D normal;


  // ****************************************************************
  // Output the vertices.
  // ****************************************************************

  for (i=0; i < numVertices; i++)
  {
    file << "v  "
      << (double)vertex[i].coord.x << "  " 
      << (double)vertex[i].coord.y << "  " 
      << (double)vertex[i].coord.z << endl;
  }

  // ****************************************************************
  // Output the vertex normals.
  // ****************************************************************

  for (i=0; i < numVertices; i++)
  {
    normal = getNormal(vertex[i].rib, vertex[i].ribPoint);
    file << "vn  "
      << (double)normal.x << "  " 
      << (double)normal.y << "  " 
      << (double)normal.z << endl;
  }

  // ****************************************************************
  // Output the faces (triangles).
  // ****************************************************************

  for (i=0; i < numTriangles; i++)
  {
    file << "f  "
      << (int)(triangle[i].vertex[0] + 1) << "//" << (int)(triangle[i].vertex[0] + 1) << "  " 
      << (int)(triangle[i].vertex[1] + 1) << "//" << (int)(triangle[i].vertex[1] + 1) << "  " 
      << (int)(triangle[i].vertex[2] + 1) << "//" << (int)(triangle[i].vertex[2] + 1) << endl;
  }

  file << endl;

  // The file is closed automatically.
  return true;
}

// ****************************************************************************
/// @brief Assign the triangles to the tiles.
///
/// Each triangle is assigned to all tiles that are overlapped by the triangles'
/// bounding box. This function must be called before getTriangleList() !
// ****************************************************************************

void Surface::prepareIntersections()
{
  const double EXTREME = 1000000.0;
  const double EPSILON = 0.1;
  int i, x, y;
  Point3D *Q;

  // Find the bounding box of the surface.

  leftBorder = EXTREME;
  rightBorder = -EXTREME;
  bottomBorder = EXTREME;
  topBorder = -EXTREME;

  for (i=0; i < numVertices; i++)
  {
    Q = &vertex[i].coord;
    if (Q->x < leftBorder)   { leftBorder = Q->x; }
    if (Q->x > rightBorder)  { rightBorder = Q->x; }
    if (Q->y < bottomBorder) { bottomBorder = Q->y; }
    if (Q->y > topBorder)    { topBorder = Q->y; }
  }

  leftBorder-= EPSILON;
  bottomBorder-= EPSILON;
  rightBorder+= EPSILON;
  topBorder+= EPSILON;

  numTilesX = (int)((rightBorder - leftBorder) / STANDARD_TILE_WIDTH) + 1;
  numTilesY = (int)((topBorder - bottomBorder) / STANDARD_TILE_HEIGHT) + 1;

  if (numTilesX < 1) { numTilesX = 1; }
  if (numTilesY < 1) { numTilesY = 1; }
  if (numTilesX > MAX_TILES_X) { numTilesX = MAX_TILES_X; }
  if (numTilesY > MAX_TILES_Y) { numTilesY = MAX_TILES_Y; }

  tileWidth = (rightBorder - leftBorder) / (double)numTilesX;
  tileHeight = (topBorder - bottomBorder) / (double)numTilesY;

  // Make all tiles empty.

  for (x=0; x < numTilesX; x++)
  {
    for (y=0; y < numTilesY; y++) { tile[x][y].numTriangles = 0; }
  }

  // Fill the tiles with the triangles.

  double minX, maxX, minY, maxY;
  int leftTile, rightTile, bottomTile, topTile;
  Point3D P;

  for (i=0; i < numTriangles; i++)
  {
    // Get the bounding box of triangle i.
    minX = minY = EXTREME;
    maxX = maxY = -EXTREME;

    for (x=0; x < 3; x++)
    {
      P = vertex[ triangle[i].vertex[x] ].coord;
      if (P.x < minX) { minX = P.x; }
      if (P.x > maxX) { maxX = P.x; }
      if (P.y < minY) { minY = P.y; }
      if (P.y > maxY) { maxY = P.y; }
    }

    minX-= EPSILON;
    maxX+= EPSILON;
    minY-= EPSILON;
    maxY+= EPSILON;

    // Add the triangle index to all tiles overlapping the bounding box.

    leftTile   = (int)((minX - leftBorder) / tileWidth);
    rightTile  = (int)((maxX - leftBorder) / tileWidth);
    bottomTile = (int)((minY - bottomBorder) / tileHeight);
    topTile    = (int)((maxY - bottomBorder) / tileHeight);

    if (leftTile < 0)   { leftTile = 0; }
    if (rightTile < 0)  { rightTile = 0; }
    if (bottomTile < 0) { bottomTile = 0; }
    if (topTile < 0)    { topTile = 0; }
    if (leftTile >= numTilesX)  { leftTile = numTilesX - 1; }
    if (rightTile >= numTilesX) { rightTile = numTilesX - 1; }
    if (bottomTile >= numTilesY) { bottomTile = numTilesY - 1; }
    if (topTile >= numTilesY)    { topTile = numTilesY - 1; }

    for (x=leftTile; x <= rightTile; x++)
    {
      for (y=bottomTile; y <= topTile; y++)
      {
        if (tile[x][y].numTriangles < MAX_TILE_TRIANGLES)
        {
          tile[x][y].triangle[ tile[x][y].numTriangles ] = i;
          tile[x][y].numTriangles++;
        }
        else
        {
          // Error !
          printf("Too many triangles for one tile !\n");
        }
      }
    }
  }
}
  
// ****************************************************************************
/// @brief Returns a list with all triangles that are possibly interesected 
/// by the intersecting plane defined by the call of 
/// prepareIntersection(Point2D Q, Point2D v).
///
/// Before this function, you must call prepareIntersections() (once
/// for a certain surface geometry) AND prepareIntersection(Point2D Q, Point2D v)
/// with the intersecting plane parameters.
/// The function returns false, when \a MAX_ENTRIES is smaller than the actual
/// necessary number of entries.
/// @param indexList The list to be filled with the triangle indices. Triangles
/// may occur multiple times in the list, when they overlap muliple tiles.
/// @param numEntries Returns the number of list entries.
/// @param MAX_ENTRIES The maximal number of entries to be put in the list.
// ****************************************************************************
  
bool Surface::getTriangleList(int *indexList, int &numEntries, int MAX_ENTRIES)
{
  int i;
  double x, y;
  int tileX, tileY;         // Index of the current tile
  double nextBorderX, nextBorderY;
  double deltaX, deltaY;
  Tile *t = NULL;
  int *source;

  Point2D Q = linePoint;
  Point2D v = lineVector;

  numEntries = 0;

  // ****************************************************************
  // A rather horizontal intersection line.
  // ****************************************************************

  if (fabs(v.y) <= fabs(v.x)*tileHeight / tileWidth)
  {
    // The vector v must be directed to the right !
    if (v.x < 0.0)
    {
      v.x = -v.x;
      v.y = -v.y;
    }

    y = Q.y + (leftBorder - Q.x)*v.y / v.x;
    deltaY = v.y*tileWidth / v.x;
    tileY = (int)((y - bottomBorder) / tileHeight);

    // v.y > 0.0 ****************************************************

    if (v.y > 0.0)
    {
      nextBorderY = bottomBorder + (double)(tileY+1)*tileHeight;

      for (tileX=0; tileX < numTilesX; tileX++)
      {
        if ((tileY >= 0) && (tileY < numTilesY))
        {
          t = &tile[tileX][tileY];
          if (numEntries + t->numTriangles > MAX_ENTRIES) { return false; }
          source = &t->triangle[0];
          for (i=t->numTriangles-1; i >= 0; i--) { *indexList++ = *source++; }
          numEntries+= t->numTriangles;
        }

        if (y + deltaY > nextBorderY)
        {
          tileY++;
          nextBorderY+= tileHeight;

          if ((tileY >= 0) && (tileY < numTilesY))
          {
            t = &tile[tileX][tileY];
            if (numEntries + t->numTriangles > MAX_ENTRIES) { return false; }
            source = &t->triangle[0];
            for (i=t->numTriangles-1; i >= 0; i--) { *indexList++ = *source++; }
            numEntries+= t->numTriangles;
          }
        }

        y+= deltaY;
      }
    }
    else

    // v.y <= 0.0 ***************************************************

    {
      nextBorderY = bottomBorder + (double)tileY*tileHeight;

      for (tileX=0; tileX < numTilesX; tileX++)
      {
        if ((tileY >= 0) && (tileY < numTilesY))
        {
          t = &tile[tileX][tileY];
          if (numEntries + t->numTriangles > MAX_ENTRIES) { return false; }
          source = &t->triangle[0];
          for (i=t->numTriangles-1; i >= 0; i--) { *indexList++ = *source++; }
          numEntries+= t->numTriangles;
        }

        if (y + deltaY < nextBorderY)
        {
          tileY--;
          nextBorderY-= tileHeight;
          t = &tile[tileX][tileY];

          if ((tileY >= 0) && (tileY < numTilesY))
          {
            if (numEntries + t->numTriangles > MAX_ENTRIES) { return false; }
            source = &t->triangle[0];
            for (i=t->numTriangles-1; i >= 0; i--) { *indexList++ = *source++; }
            numEntries+= t->numTriangles;
          }
        }

        y+= deltaY;
      }
    }
  }
  else

  // ****************************************************************
  // A rather vertical line.
  // ****************************************************************

  {
    // The vector v must be directed upwards !
    if (v.y < 0.0)
    {
      v.x = -v.x;
      v.y = -v.y;
    }

    x = Q.x + (bottomBorder - Q.y)*v.x / v.y;
    deltaX = v.x*tileHeight / v.y;
    tileX = (int)((x - leftBorder) / tileWidth);

    // v.x > 0.0 ****************************************************

    if (v.x > 0.0)
    {
      nextBorderX = leftBorder + (double)(tileX+1)*tileWidth;

      for (tileY=0; tileY < numTilesY; tileY++)
      {
        if ((tileX >= 0) && (tileX < numTilesX))
        {
          t = &tile[tileX][tileY];
          if (numEntries + t->numTriangles > MAX_ENTRIES) { return false; }
          source = &t->triangle[0];
          for (i=t->numTriangles-1; i >= 0; i--) { *indexList++ = *source++; }
          numEntries+= t->numTriangles;
        }

        if (x + deltaX > nextBorderX)
        {
          tileX++;
          nextBorderX+= tileWidth;

          if ((tileX >= 0) && (tileX < numTilesX))
          {
            t = &tile[tileX][tileY];
            if (numEntries + t->numTriangles > MAX_ENTRIES) { return false; }
            source = &t->triangle[0];
            for (i=t->numTriangles-1; i >= 0; i--) { *indexList++ = *source++; }
            numEntries+= t->numTriangles;
          }
        }

        x+= deltaX;
      }
    }
    else

    // v.x <= 0.0 ***************************************************

    {
      nextBorderX = leftBorder + (double)tileX*tileWidth;

      for (tileY=0; tileY < numTilesY; tileY++)
      {
        if ((tileX >= 0) && (tileX < numTilesX))
        {
          t = &tile[tileX][tileY];
          if (numEntries + t->numTriangles > MAX_ENTRIES) { return false; }
          source = &t->triangle[0];
          for (i=t->numTriangles-1; i >= 0; i--) { *indexList++ = *source++; }
          numEntries+= t->numTriangles;
        }

        if (x + deltaX < nextBorderX)
        {
          tileX--;
          nextBorderX-= tileWidth;

          if ((tileX >= 0) && (tileX < numTilesX))
          {
            t = &tile[tileX][tileY];
            if (numEntries + t->numTriangles > MAX_ENTRIES) { return false; }
            source = &t->triangle[0];
            for (i=t->numTriangles-1; i >= 0; i--) { *indexList++ = *source++; }
            numEntries+= t->numTriangles;
          }
        }

        x+= deltaX;
      }
    }
  }

  return true;
}

// ****************************************************************************
/// @brief Prepares the intersection of the surface with the given intersecting
/// plane/line. This function must be called before
/// getTriangleIntersection(int index, Point2D &P0, Point2D &P1, Point2D &n)
/// can be called to test the individual triangles for an intersection.
///
/// @param Q A point on the intersecting line.
/// @param v The direction of the intersecting line in the xy-plane.
// ****************************************************************************
  
void Surface::prepareIntersection(Point2D Q, Point2D v)
{
  const double EPSILON = 0.000001;
  int i;

  for (i=0; i < numVertices; i++) { vertex[i].wasTested = false; }
  for (i=0; i < numEdges; i++)    { edge[i].wasTested = false; }

  v.normalize();            // Normalisierung ist wichtig !
  lineVector = v;           // Für getTriangleIntersection merken
  linePoint = Q;

  // Einen Normaleneinheitsvektor bilden, der senkrecht (90° nach links
  // gedreht) auf v steht.

  Point2D n(-v.y, v.x);

  // Die zwei Punkte leftLinePoint und rightLinePoint berechnen, die im 
  // Abstand EPSILON links bzw. rechts der Gerade Q+t*v auf der Höhe von
  // P liegen.

  leftLinePoint  = Q + EPSILON*n;
  rightLinePoint = Q - EPSILON*n;
}

// ****************************************************************************
/// @brief Returns the intersection data for a single triangle.
///
/// The function prepareIntersection(Point2D Q, Point2D v) must be called before
/// this function to define the intersecting plane.
/// @param index The index of the triangle.
/// @param P0 The first point of the intersection line projected onto the 
/// intersecting plane.
/// @param P1 The second point of the intersection line projected onto the 
/// intersecting plane.
/// @param n The triangle normal projected onto the intersecting plane.
// ****************************************************************************
  
bool Surface::getTriangleIntersection(int index, Point2D &P0, Point2D &P1, Point2D &n)
{
  int e0, e1, e2;   // Die 3 Kanten des Dreiecks

  e0 = triangle[index].edge[0];
  e1 = triangle[index].edge[1];
  e2 = triangle[index].edge[2];

  int numIntersections = 0;
  Point2D Q[3];

  if (getEdgeIntersection(e0)) { Q[numIntersections++] = edge[e0].intersection; }
  if (getEdgeIntersection(e1)) { Q[numIntersections++] = edge[e1].intersection; }
  if (getEdgeIntersection(e2)) { Q[numIntersections++] = edge[e2].intersection; }

  if (numIntersections < 2) { return false; }

  // Den Normalenvektor der Fläche berechnen.

  int v0, v1, v2;
  v0 = triangle[index].vertex[0];
  v1 = triangle[index].vertex[1];
  v2 = triangle[index].vertex[2];

  Point3D normal = crossProduct(vertex[v1].coord - vertex[v0].coord, vertex[v2].coord - vertex[v0].coord);

  // ****************************************************************
  // Die Projektion des Normalenvektors auf die Dreicksfläche
  // Dazu muss lineVector die Länge 1 haben !
  // ****************************************************************

  n.x = normal.z;
  n.y = normal.x*lineVector.x + normal.y*lineVector.y;

  // Bei zwei Schnittpunkten P0 und P1 zurückgeben.

  if (numIntersections == 2)
  {
    P0 = Q[0];
    P1 = Q[1];
    return true;
  }

  // numIntersections > 2:

  double l[3];
  l[0] = (Q[0]-Q[1]).squareMagnitude();
  l[1] = (Q[1]-Q[2]).squareMagnitude();
  l[2] = (Q[2]-Q[0]).squareMagnitude();

  if ((l[0] >= l[1]) && (l[0] >= l[2]))
  {
    P0 = Q[0];
    P1 = Q[1];
  }
  else
  {
    if (l[1] >= l[2])
    {
      P0 = Q[1];
      P1 = Q[2];
    }
    else
    {
      P0 = Q[2];
      P1 = Q[0];
    }
  }

  return true;
}

// ****************************************************************************
/// @brief This function returns true, if a triangle edge was interesected by 
/// the intersecting plane. Do not call this function explicitely!
// ****************************************************************************

bool Surface::getEdgeIntersection(int edgeIndex)
{
  if (edge[edgeIndex].wasTested) { return edge[edgeIndex].isIntersected; }

  // The edge must be tested for an intersection.

  Point2D w;
  double d;

  edge[edgeIndex].wasTested = true;

  int v0 = edge[edgeIndex].vertex[0];
  int v1 = edge[edgeIndex].vertex[1];

  // ****************************************************************
  // Is the first vertex left or right from the intersection line ?
  // Im reserved-Feld jedes Punktes vermerken, ob er links von (res = -1), 
  // rechts von (res = +1) oder (innerhalb einer EPSILON-Umgebung) 
  // auf der Schnittebene liegt (res = 0).
  // ****************************************************************

  if (vertex[v0].wasTested == false)
  {
    vertex[v0].wasTested = true;
    vertex[v0].reserved = 0;

    w.x = vertex[v0].coord.x - leftLinePoint.x;
    w.y = vertex[v0].coord.y - leftLinePoint.y;
    d = w.x*lineVector.y - w.y*lineVector.x;
    if (d < 0.0) { vertex[v0].reserved = -1; }

    w.x = vertex[v0].coord.x - rightLinePoint.x;
    w.y = vertex[v0].coord.y - rightLinePoint.y;
    d = w.x*lineVector.y - w.y*lineVector.x;
    if (d > 0.0) { vertex[v0].reserved = 1; }
  }

  // Is the second vertex left or right from the intersection line ?

  if (vertex[v1].wasTested == false)
  {
    vertex[v1].wasTested = true;
    vertex[v1].reserved = 0;

    w.x = vertex[v1].coord.x - leftLinePoint.x;
    w.y = vertex[v1].coord.y - leftLinePoint.y;
    d = w.x*lineVector.y - w.y*lineVector.x;
    if (d < 0.0) { vertex[v1].reserved = -1; }

    w.x = vertex[v1].coord.x - rightLinePoint.x;
    w.y = vertex[v1].coord.y - rightLinePoint.y;
    d = w.x*lineVector.y - w.y*lineVector.x;
    if (d > 0.0) { vertex[v1].reserved = 1; }
  }

  // Test the edge for an intersection.

  edge[edgeIndex].isIntersected = false;

  if (((vertex[v0].reserved >= 0) && (vertex[v1].reserved <= 0)) ||
      ((vertex[v0].reserved <= 0) && (vertex[v1].reserved >= 0)))
  {
    // Den Schnittpunkt der Kante mit der Schnittebene genau bestimmen.
    const double EPSILON = 0.000001;
    Point3D P, u, R;
    double denominator;

    // Geradengleichung der Kante: X = P + t*u.
    P = vertex[v0].coord;
    u = vertex[v1].coord - P;

    R.x = P.x - linePoint.x;
    R.y = P.y - linePoint.y;
    R.z = P.z;

    denominator = -u.x*lineVector.y + u.y*lineVector.x;

    if (denominator != 0.0)
    {
      // Liegt der Parameter d der Kante zwischen -EPSILON und 1+EPSILON ?
      d = (-lineVector.x*R.y + lineVector.y*R.x) / denominator;
      if ((d >= -EPSILON) && (d < 1.0+EPSILON))
      {
        edge[edgeIndex].isIntersected = true;
        edge[edgeIndex].intersection.x = (lineVector.x*(u.y*R.z - u.z*R.y) + lineVector.y*(u.z*R.x - u.x*R.z)) / denominator;
        edge[edgeIndex].intersection.y = (-u.x*R.y + u.y*R.x) / denominator;
      }
    }
  }

  return edge[edgeIndex].isIntersected;
}

// ****************************************************************************
/// @brief Appends the vertex data of this surface in binary format to a file.
/// This function does NOT open or close the file !
///
/// @param hFile Handle to the file.
// ****************************************************************************

void Surface::appendToFile(std::ofstream &file)
{
  int numNumbers = numRibs*numRibPoints*3;
  double *buffer = new double[numNumbers];
  int i, k;
  Point3D Q;

  // Gitterdimension.

  file.write((char*)&numRibPoints, sizeof(int));
  file.write((char*)&numRibs, sizeof(int));

  // Gitterpunkte.

  for (i=0; i < numRibPoints; i++)
  {
    for (k=0; k < numRibs; k++)
    {
      Q = getVertex(k, i);
      buffer[(i*numRibs + k)*3]   = Q.x;
      buffer[(i*numRibs + k)*3+1] = Q.y;
      buffer[(i*numRibs + k)*3+2] = Q.z;
    }
  }

  if (!file.write((char*)buffer, numNumbers*sizeof(double)))
  {
    printf("Error in Surface::appendToFile(): WriteFile() failed!\n");
  }

  delete[] buffer;
}

// ****************************************************************************
/// @brief Reads the vertex data for this surface from a binary file. Opening
/// and closing of the file is not done here !
///
/// @param hFile The file to read from.
/// @param initialize Shall this surface be initialized with the dimensions
/// read from the file or was it already initialized ?
// ****************************************************************************

void Surface::readFromFile(std::ifstream &file, bool initialize)
{
  int i, k;
  Point3D Q;

  // Gitterdimensionen auslesen.

  file.read((char*)&numRibPoints, sizeof(int));
  file.read((char*)&numRibs, sizeof(int));

  // Neue Struktur anlegen.

  if (initialize) { init(numRibs, numRibPoints); }

  // Daten aus der Datei einlesen.

  int numNumbers = numRibs*numRibPoints*3;
  double *buffer = new double[numNumbers];

  if (!file.read((char*)buffer, numNumbers*sizeof(double)))
  {
    printf("Error in Surface::readFromFile(): ReadFile() failed!\n");
  }
  
  // Gitterpunkte.

  for (i=0; i < numRibPoints; i++)
  {
    for (k=0; k < numRibs; k++)
    {
      Q.x = buffer[(i*numRibs + k)*3];
      Q.y = buffer[(i*numRibs + k)*3+1];
      Q.z = buffer[(i*numRibs + k)*3+2];
      setVertex(k, i, Q);
    }
  }

  delete[] buffer;
}

// ****************************************************************************
