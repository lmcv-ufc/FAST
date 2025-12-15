// -------------------------------------------------------------------------
// bspbasis.h - file containing the definition of the cBSplineBasis class.
// -------------------------------------------------------------------------
// Copyright (c) 2015 LMCV/UFC
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
// The cBSplineBasis class implements the univariate bspline basis object 
// used by BSpline and NURBS curves, surfaces and solids.
//
// -------------------------------------------------------------------------
// Public methods:
// -------------------------------------------------------------------------
//
// cBsplineBasis(void)
//
// Construct a Bspline object with empty knot vector.
// -------------------------------------------------------------------------
//
// cBsplineBasis(const cBsplineBasis &bsp)
//
//   bsp  - given copy object                                      (in)
//
// Construct a BSpline as a copy of given object (knot vector and degree).
// -------------------------------------------------------------------------
//
// int getBasInd(const int &span) const
//
//   span  - given knot span                                       (in)
//
// This method return b-spline basis index for given knot span.
// -------------------------------------------------------------------------
//
// const cKnotVector& getKnotVec(void) const
//
// This method return a const reference of its knot vector object.
// -------------------------------------------------------------------------
//
// int FindSpan(const double &p) const
//
//   p - given parametric location.                                (in)
//
// This method call FindSpan from knot vector object.
// -------------------------------------------------------------------------
//
// void setDegree(const int &deg) 
//
//   deg - input b-spline basis degree                             (in)
//
// This method set basis degree.
// -------------------------------------------------------------------------
//
// void setKnotVec(const cKnotVector &kv)
//
//   kv - given knot vector object.                                (in)
//
// This method set basis knot vector.
// -------------------------------------------------------------------------
//
// double CompBasFunc(int bind, int deg, double t) const
//
//  bind - given bernstein basis index                             (in)
//  deg  - given degree                                            (in)
//  t    - given parametric location                               (in)
//
// This method return value of single b-spline basis at given parametric
// location.
// -------------------------------------------------------------------------
//
// void CompBasFunc(int u, double *N) const
//
//  u  - given parametric location                                 (in)
//  B  - B-Spline basis functions                                  (out)
//
// This function computes B-Spline non-zero basis functions at given parametric
// location.
// -------------------------------------------------------------------------
//
// double CompBasDerv(int bind, int k, double t) const
//
//   bind - given bernstein basis index                            (in)
//   k    - derivative order                                       (in)
//   t    - given parametric location                              (in)
//
// This method return value of single basis derivative at given parametric
// location.
// -------------------------------------------------------------------------
//
// void CompBasDerv(double u, int n, double** &ders) const
//
//  u  - given parametric location                                 (in)
//  n  - derivative order                                          (in)
//  ders  - B-Spline function derivatives.
//
// This function computes all basis functions derivatives and store them in
// matrix ders. The old content of ders is not destroyed.
// -------------------------------------------------------------------------
//
// double EvalAlpha(int bind, int k, int j) const
//
//   bind - given bernstein basis index                            (in)
//   k    - derivative order                                       (in)
//   j    - ??? see in NURBS book!                                 (in)
//
// This method is used in the implementatin of single derivative evaluation.
//
// -------------------------------------------------------------------------
// Friend Methods:
// -------------------------------------------------------------------------
//
// istream& operator>>(istream &in, cBSplineBasis &bsp)
//
//   in  - given istream object                                    (in/out)
//   bsp - given knot vector object                                (in)
//
// This method read BSpline data from given input stream object.
// -------------------------------------------------------------------------
//
// ostream& operator<<(ostream &out, cBSplineBasis &bsp)
//
//   out  - given ostream object                                   (in/out)
//   kvec - given knot vector object                               (in/out)
//
// This method writes BSpline data in given output stream object.
// -------------------------------------------------------------------------

#ifndef _BSPBASIS_H
#define _BSPBASIS_H

#include "knotvec.h"

// -------------------------------------------------------------------------
// Forward declarations:
//
#include <iosfwd>

// -------------------------------------------------------------------------
// Definition of cBSplineBasis class:
//
class cBSplineBasis
{
 protected:
  cKnotVector  KnotVec;     // Knot vector.
  int          Degree;      // Degree of B-Spline basis fuction
  int          NumBasis;    // Number of basis functions

 public:
                      cBSplineBasis(void) { }
                      cBSplineBasis(const cBSplineBasis&);
                     ~cBSplineBasis(void) { }

         int          getBasInd(const int&) const;
         int          getDegree(void) const {return Degree;}
         int          getNumBas(void) const {return NumBasis;}
  const  cKnotVector& getKnotVec(void) const;
	 int          FindSpan(const double &p) const;
         void         setDegree(const int &deg);
	 void         setKnotVec(const cKnotVector&);
         double       CompBasFunc(int,int,double) const;
         void         CompBasFunc(double,double*) const;
         double       CompBasDerv(int,int,double) const;
         void         CompBasDerv(double,int,double**&) const;
         double       EvalAlpha(int,int,int) const;

  friend std::ostream& operator<<(std::ostream&,const cBSplineBasis&);
  friend std::istream& operator>>(std::istream&,cBSplineBasis&);
};

#endif
