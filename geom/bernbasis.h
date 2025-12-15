// --------------------------------------------------------------------------
// bernbasis.h - file containing the definition of the cBernsteinBasis class.
// --------------------------------------------------------------------------
// Copyright (c) 2018 LMCV/UFC
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without 
// modification, are permitted provided that the following conditions
// are met:
//
// 1. Redistributions of source code must retain the above copyright 
//    notice, this list of conditions and the following disclaimer.
// 2. Redistributions in binary form must reproduce the above copyright 
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS 
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT 
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
// PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT 
// HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, 
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED 
// TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR 
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF 
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING 
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS 
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// -------------------------------------------------------------------------
//
// The nBernstein namespace contains routines related to univariate and
// bivariate (used in Bézier triangles) Bernstein basis.
//
// -------------------------------------------------------------------------
// Public methods:
// -------------------------------------------------------------------------
//
// double CompBasVal(int bind, int deg, double par)
//
//   bind - given bernstein basis index                            (in)
//   deg  - given degree                                           (in)
//   par  - given parametric location                              (in)
//
// This method return value of single bernstein basis at given parametric
// location. Input parameter must be in interval [0,1].
// -------------------------------------------------------------------------
//
// double CompBasDerv(int bind, int deg, double par)
//
//   bind - given bernstein basis index                            (in)
//   deg  - given degree                                           (in)
//   par  - given parametric location                              (in)
//
// This method return frist derivative of single bernstein basis at given
// parametric location. Input parameter must be in interval [0,1].
// -------------------------------------------------------------------------
//
// void CompPols(int n, double u, double *B)
//
//  n  -  degree of the polynomials                                (in)
//  u  -  parametric coordinate (0 <= u <= 1)                      (in)
//  B  -  Bernstein polynomials                                    (out)
//
// This function computes n+1 Bernstein polynomials at given parametric location
// (Algorithm A1.3 - pp. 20-21).
// -------------------------------------------------------------------------
//
// void CompDervPols(int n, double u, double *dBdu)
//
//  n    -  degree of the polynomials                              (in)
//  u    -  parametric coordinate (0 <= u <= 1)                    (in)
//  dBdu -  derivative of Bernstein polynomials                    (out)
//
// This function computes first derivatives of n+1 Bernstein polynomials.
// -------------------------------------------------------------------------
//
// void CompSurfPols(int dr, int ds, double r, double s, double *B)
//
//  dr -  degree in r direction                                    (in)
//  ds -  degree in s direction                                    (in)
//  r  -  parametric location in r direction                       (in)
//  s  -  parametric location in s direction                       (in)
//  B  -  Tensor-product surface Bernstein polynomials             (out)
//
// This function computes tensor-product surface Bernstein polynomials.
// -------------------------------------------------------------------------
//
// int CompTrianNumBas(const int &d)
//
//  d -  degree of Bézier triangle                                 (in)
//
// This function return number of basis of Bézier triangle.
// -------------------------------------------------------------------------
//
// int MapTrianBaryToVecID(const int &i, const int &j, const int &d)
//
//  d -  degree of Bézier triangle                                 (in)
//  i -  Triangular barycentric index i                            (in)
//  j -  Triangular barycentric index j                            (in)
//
// This function map the given triangular barycentric indexes to a vector index.
// -------------------------------------------------------------------------
//
// void IncTrianBaryID(const int &d, int &i, int &j, int &k)
//
//  d     -  degree of Bézier triangle                             (in)
//  i,j,k -  Triangular barycentric indexes                        (in/out)
//
// This function incremente barycentric indices according with the scheme
// discussed in Farin (2002), "Curves and Surfaces for CAGD A Pratical Guide".
// -------------------------------------------------------------------------
//
// void CompTrianPols(int p, double r, double s, double *B)
//
//  d -  degree of Bézier triangle                                 (in)
//  r -  Barycentric value in r direction                          (in)
//  s -  Barycentric value in s direction                          (in)
//  B -  Bézier triangle basis                                     (out)
//
// This function computes Bézier triangle basis functions.
// -------------------------------------------------------------------------
//
// void CompTrianDervPols(int p, double r, double s, double *dr, double *ds)
//
//  d  -  degree of Bézier triangle                                (in)
//  r  -  Barycentric value in r direction                         (in)
//  s  -  Barycentric value in s direction                         (in)
//  dr -  Bézier triangle basis derivative in r direction          (out)
//  ds -  Bézier triangle basis derivative in s direction          (out)
//
// This function computes Bézier triangle basis derivatives.
// -------------------------------------------------------------------------

#ifndef _BERNBASIS_H
#define _BERNBASIS_H

// -------------------------------------------------------------------------
// Definition of cKnotVector class:
//

namespace nBernsteinBasis
{
  // Univariate Bernstein basis.
  double CompBasVal(int,int,double);
  double CompBasDerv(int,int,double);
  void   CompPols(int,double,double*&);
  void   CompDervPols(int,double,double*&);

  // Bernstein Bézier tensor product surface.
  void   CompSurfPols(int,int,double,double,double*);

  // Bernstein Bézier Triangle functions.
  int    CompTrianNumBas(const int&);
  int    MapTrianBaryToVecID(const int&,const int&,const int&);
  void   IncTrianBaryID(const int&,int&,int&,int&);
  void   CompTrianPols(int,double,double,double*);
  void   CompTrianDervPols(int,double,double,double*,double*);
}

#endif
