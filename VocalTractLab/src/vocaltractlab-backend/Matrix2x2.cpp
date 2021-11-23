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

#include "Matrix2x2.h"

// ****************************************************************************
// ****************************************************************************

Matrix2x2::Matrix2x2(ComplexValue a, ComplexValue b, ComplexValue c, ComplexValue d)
{
  A = a; 
  B = b; 
  C = c; 
  D = d;
}

// ****************************************************************************
// ****************************************************************************

Matrix2x2::Matrix2x2(void)
{
  A = ComplexValue(0, 0);
  B = ComplexValue(0, 0);
  C = ComplexValue(0, 0);
  D = ComplexValue(0, 0);
}

// ****************************************************************************
// ****************************************************************************

void Matrix2x2::unitMatrix(void)
{
  A = ComplexValue(1, 0);
  B = ComplexValue(0, 0);
  C = ComplexValue(0, 0);
  D = ComplexValue(1, 0);
}

// ****************************************************************************
// Invert the matrix. A requirement is that the determinant = 1 !!
// ****************************************************************************

void Matrix2x2::invert()
{
  ComplexValue t;

  t = A;
  A = D;
  D = t;
  B = -B;
  C = -C;
}

// ****************************************************************************
// ****************************************************************************

Matrix2x2 &Matrix2x2::operator+=(const Matrix2x2 &x)
{
  A+= x.A;
  B+= x.B;
  C+= x.C;
  D+= x.D;
  return(*this);
}

// ****************************************************************************
// ****************************************************************************

Matrix2x2 operator+(const Matrix2x2 x, const Matrix2x2 y)
{
  Matrix2x2 z = x;
  return(z+=y);
}

// ****************************************************************************
// ****************************************************************************

Matrix2x2 &Matrix2x2::operator*=(const Matrix2x2 &x)
{
  ComplexValue a, b, c, d;

  a = A*x.A + B*x.C;
  b = A*x.B + B*x.D;
  c = C*x.A + D*x.C;
  d = C*x.B + D*x.D;

  A = a; B = b; C = c; D = d;
  
  return(*this);
}

// ****************************************************************************
// ****************************************************************************

Matrix2x2 operator*(const Matrix2x2 x, const Matrix2x2 y)
{
  Matrix2x2 z = x;
  return(z*=y);
}

// ****************************************************************************
// ****************************************************************************

void Matrix2x2::operator =(const Matrix2x2 &x)
{
  A = x.A;
  B = x.B;
  C = x.C;
  D = x.D;
}

// ****************************************************************************
