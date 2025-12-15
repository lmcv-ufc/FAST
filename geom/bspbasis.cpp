// -------------------------------------------------------------------------
// bspbasis.cpp - implementation of cBSplineBasis class.
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
// Created:      28-Jan-2015    Elias Saraiva Barroso
//
// Modified:     
// -------------------------------------------------------------------------

#include <cstdlib>
#include <cassert>
#include <iostream>
#include <fstream>
#include <cmath>

#include "bspbasis.h"
#include "knotvec.h"

using namespace std;

// -------------------------------------------------------------------------
// Public methods:
//

// ============================= cBSplineBasis =============================

cBSplineBasis :: cBSplineBasis(const cBSplineBasis &bsp) : KnotVec(bsp.KnotVec)
{ 
  Degree   = bsp.Degree;
  NumBasis = bsp.NumBasis;
}

// ============================= getBasInd =================================

int cBSplineBasis :: getBasInd(const int &span) const
{
  return KnotVec.FindSpan(span) - Degree;
}

// ============================= getKnotVec ================================

const cKnotVector&  cBSplineBasis :: getKnotVec( ) const
{
  return KnotVec;
}

// ============================= setDegree =================================

void cBSplineBasis :: setDegree(const int &deg)
{
  Degree   = deg;
  NumBasis = KnotVec.getSize( ) - Degree - 1;
}

// ============================= setKnotVec ================================

void cBSplineBasis :: setKnotVec(const cKnotVector &kv)
{
  KnotVec  = kv;
  NumBasis = KnotVec.getSize( ) - Degree - 1;
}

// ============================= FindSpan ==================================

int cBSplineBasis :: FindSpan(const double &p) const
{
  return KnotVec.FindSpan(NumBasis,Degree,p);
}

// ============================= CompBasFunc ===============================

double cBSplineBasis :: CompBasFunc(int bind, int deg, double t) const
{
  if (deg == 0)  // Basic case
  {
    if(KnotVec[bind] <= t && t <= KnotVec[bind+1])
      return 1;
    else
      return 0;
  }
  else // Special case
  {
    double f1, f2, Denom; 

    // First component.
    f1 = t - KnotVec[bind];
    Denom = KnotVec[bind+deg] - KnotVec[bind];
    f1 = (fabs(Denom) < 10e-6) ? 0.0 : f1/Denom;

    // Second component.
    f2 = KnotVec[bind+deg+1] - t;
    Denom = KnotVec[bind+deg+1] - KnotVec[bind+1];
    f2 = (fabs(Denom) < 10e-6) ? 0.0 : f2/Denom;

    return f1*CompBasFunc(bind,deg-1,t) + f2*CompBasFunc(bind+1,deg-1,t);
  }
}

// ============================= CompBasFunc ===============================

void cBSplineBasis :: CompBasFunc(double u, double* N) const
{
  // Auxiliary variables.
  int    size;
  double temp, saved;
  double *left, *right;

  // Alloc memory for dynamic variables.
  size  = Degree + 1;
  left  = new double [size];
  right = new double [size];

  // Find span which u lies.
  int span = KnotVec.FindSpan(NumBasis,Degree,u);

  // Compute the nonvanishing basis functions.
  N[0] = 1.0;
  for(int j = 1; j <= Degree; ++j)
  {
    left[j]  = u - KnotVec[span+1-j];
    right[j] = KnotVec[span+j] - u;
    saved    = 0.0;
    
    for(int r = 0; r < j; ++r)
    {
      temp  = N[r]/(right[r+1]+left[j-r]);
      N[r]  = saved+right[r+1]*temp;
      saved = left[j-r]*temp;
    }

    N[j] = saved;
  }

  // Release memory.
  delete [] left;
  delete [] right;
}

// ============================= CompBasDerv ==================================

double cBSplineBasis :: CompBasDerv(int bind, int k, double t) const
{ //TODO:Verificar se a funçaõ está correta com o livro nurbs book
  // Derivate order (k) must be less than Degree
  if (k >= Degree) 
  {
    cout << "Error in evaluate process of Bspline basis, (K >= Degre)";
    exit(0);
  }
  
  // Compute the first fraction of expression

  double frac_Num = 1.0, frac_Denom = 1.0;

  for (int i = Degree; i > 0; i--)
    frac_Num *= i;

  for (int i = (Degree - k); i > 0; i--)
    frac_Denom *= i;

  // Compute the sumation

  double sum = 0;

  for (int j = 0; j <= k; j++)
    sum += EvalAlpha(bind,k,j)*CompBasFunc(bind+j,Degree-k,t);
  
  double derv = (frac_Num/frac_Denom) * sum;  

  return derv;
}

// ============================= CompBasDerv ===============================

void cBSplineBasis :: CompBasDerv(double u, int n, double** &ders) const
{
  // Check if the input is valid.
  assert (n > Degree);

  // Auxiliary variables.
  int i,j,k,r;                 // Loop counters.
  int p;
  int size;
  int s1,s2,rk,pk,j1,j2;
  double d;
  double temp, saved;
  double *left, *right;
  double **ndu, **a;

  // Alloc memory for dynamic variables.
  p     = Degree;
  size  = p + 1;
  left  = new double  [size];
  right = new double  [size];
  ders  = new double* [size];
  ndu   = new double* [size];
  a     = new double* [2];

  for(i = 0; i < n+1; ++i)  ders[i] = new double [size];
  for(i = 0; i < size; ++i) ndu[i]  = new double [size];
  for(i = 0; i < size; ++i) a[i]    = new double [size];

  // Find span which u lies.
  int span = KnotVec.FindSpan(NumBasis,p,u);

  ndu[0][0] = 1.0;
  for(j = 1; j <= p; ++j)
  {
    left[j]  = u - KnotVec[span+1-j];
    right[j] = KnotVec[span+j]-u;
    saved = 0.0;
    for(r = 0; r < j; ++r)
    {
      // Lower triangle.
      ndu[j][r] = right[r+1] + left[j-r];
      temp      = ndu[r][j-1] / ndu[j][r];

      // Upper triangle.
      ndu[r][j] = saved + right[r+1] * temp;
      saved = left[j-r] * temp;
    }
    ndu[j][j] = saved;
  }

  // Load the basis functions.
  for(j = 0; j <= p; ++j)
    ders[0][j] = ndu[j][p];

  // This section computes the derivatives (Eq. [2.9]).
  // Loop over function index.
  for(r = 0; r <= p; ++r)
  {
    s1 = 0;   s2 = 1;
    a[0][0] = 1.0;

    // Loop to compute kth derivative.
    for(k = 0; k <= n; ++k)
    {
      d = 0.0;
      rk = r - k;   pk = p - k;

      if (r >= k)
      {
        a[s2][0] = a[s1][0]/ndu[pk+1][rk];
        d = a[s2][0]*ndu[rk][pk];
      }
      if (rk >= -1)
        j1 = 1;
      else
        j1 = -rk;
      if (r - 1 <= pk)
        j2 = k-1;
      else
        j2 = p-r;

      for(j = j1; j <= j2; ++j)
      {
        a[s2][j] = (a[s1][j]-a[s1][j-1])/ndu[pk+1][rk+j];
        d += a[s2][j]*ndu[rk+j][pk];
      }
      if (r <= pk)
      {
        a[s2][k] = -a[s1][k-1]/ndu[pk+1][r];
        d += a[s2][k] * ndu[r][pk];
      }
      ders[k][r] = d;

      // Switch rows.
      j = s1;   s1 = s2;   s2 = j;
    }
  }
  // Multiply through by the correct factors (Eq. [2.9]).
  r = p;
  for (k = 1; k <= n; ++k)
  {
    for(j = 0; j <= p; ++k)
      ders[k][j] *= r;
    r *= (p-k);
  }

  // Release dynamic variables.
  for(i = 0; i < size; ++i)
  {
    delete [] ndu[i];
    delete [] a[i];
  }

  delete [] left;
  delete [] right;
  delete [] ders;
  delete [] ndu;
  delete [] a;
}

// ============================= EvalAlpha =================================

double cBSplineBasis :: EvalAlpha(int bind,int k, int j) const
{
  if (k == 0 && j == 0)
    return 1;
  
  double frac = 0.0, Denom = 1.0;

  if (k != 0 && j == 0)
  {
    frac = EvalAlpha(bind,k-1,0);

    Denom = KnotVec[bind+Degree-k] - KnotVec[bind-1];
  }
  else if (k != 0 && (j >= 1 && j < k))
  {
    frac = EvalAlpha(bind,k-1,j) - EvalAlpha(bind,k-1,j-1);

    Denom = KnotVec[bind+Degree+j-k] - KnotVec[bind+j-1];
  }
  else if (k != 0 && j == k)
  {
    frac = -EvalAlpha(bind,k-1,k-1);
    Denom = KnotVec[bind+Degree] - KnotVec[bind+k-1];
  }

  if (fabs(Denom) < 1e-6)
    return 0.0;
  else
    return frac / Denom;
}


// ============================== operator<< ===================================

ostream& operator<< (ostream &out, const cBSplineBasis &bsp)
{
  out << bsp.getDegree( ) << " "; 
  out << bsp.KnotVec;
  
  return out;
}

// ============================== operator>> ===================================

istream& operator>> (istream &in,cBSplineBasis &bsp)
{
  if (!(in >> bsp.Degree) || !(in >> bsp.KnotVec))
  {
    cout << "Error in the input of B-Spline data." << endl;
    in.setstate(ios::failbit);
    return in;
  }

  bsp.NumBasis = bsp.KnotVec.getSize( ) - bsp.Degree - 1;

  return in;
}

// ======================================================= End of file =====
