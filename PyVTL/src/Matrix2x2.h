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

#ifndef __MATRIX2X2_H__
#define __MATRIX2X2_H__

#include <complex>

typedef std::complex<double> ComplexValue;

// *********************************************************************************
// 2 x 2 matrices with complex numbers.
// *********************************************************************************

class Matrix2x2
{
  public:
    ComplexValue A, B, C, D;

    Matrix2x2(ComplexValue a, ComplexValue b, ComplexValue c, ComplexValue d);
    Matrix2x2(void);
    void unitMatrix(void);
    void invert();
    Matrix2x2 &operator*= (const Matrix2x2 &x);
    Matrix2x2 &operator+= (const Matrix2x2 &x);
    void operator=(const Matrix2x2 &x);
};

Matrix2x2 operator*(const Matrix2x2 x, const Matrix2x2 y);
Matrix2x2 operator+(const Matrix2x2 x, const Matrix2x2 y);


#endif
