// -------------------------------------------------------------------------
// bernbasis.cpp - implementation of nBernsteinBasis namespace routines.
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
// Created:      29-Jan-2015    Elias Saraiva Barroso
//
// Modified:     16-Jul-2016    Elias Saraiva Barroso
//		 Creation of bernstein triangle basis functions
// -------------------------------------------------------------------------


#include <iostream>
#include <stdlib.h>
#include "bernbasis.h"

using namespace std;

// -------------------------------------------------------------------------
// Univariate Bernstein basis:
//

// ============================= CompBasVal ===============================

double nBernsteinBasis :: CompBasVal(int bind,int deg, double par)
{
  if (bind == 1 && deg == 0)
    return 1;
  
  if (bind < 1 || bind > deg+1)
    return 0;

  double b1 = CompBasVal(bind,deg-1,par);
  double b2 = CompBasVal(bind-1,deg-1,par);
  
  return 0.5*(1-par)*b1+0.5*(1+par)*b2;
}

// ============================= CompBasDerv ===============================

double nBernsteinBasis :: CompBasDerv(int bind,int deg, double par)
{
  int basis, degree;

  double b1;
  basis = bind - 1;
  degree = deg - 1;

  if (basis < -1 || basis > deg)
    b1 = 0;
  else
    b1 = CompBasVal(basis, degree, par);

  double b2;
  basis = bind;
  degree = deg - 1;

  if (basis < -1 || basis > deg)
    b2 = 0;
  else
    b2 = CompBasVal(basis, degree, par);

  return deg*(b1-b2);
}

// ========================== CompPols =====================================

void nBernsteinBasis :: CompPols(int n, double u, double* &B)
{
  B[0] = 1.0;
  double u1 = 1.0 - u;
  for (int j = 1; j <= n; j++)
  {
    double saved = 0.0;
    for (int k = 0; k < j; k++)
    {
      double temp = B[k];
      B[k] = saved + u1*temp;
      saved = u*temp;
    }
    B[j] = saved;
  }
}

// ========================= DerivBernsteinPols ============================

void nBernsteinBasis :: CompDervPols(int n, double u, double* &dBdu)
{
  // Compute the Bernstein polynomials of n-1 degree.
  CompPols(n-1, u, dBdu);

  // Compute the derivatives using property P1.7 (pp. 17).
  dBdu[n] = 0.0;      // B_n,n-1  = 0
  double prev = 0.0;  // B_-1,n-1 = 0
  for (int i = 0; i <= n; i++)
  {
    double curr = dBdu[i];
    dBdu[i] = n*(prev - curr);  // dBdu_i,n = J*n*(B_i-1,n-1 - B_i,n-1)
    prev = curr;
  }
}

// ========================= CompSurfPols ==================================

void nBernsteinBasis :: CompSurfPols(int dr, int ds, double r, double s, double *B)
{
  // Compute univariate bernstein polynomials.
  double *basR  = new double [dr+1];
  double *basS  = new double [ds+1];

  CompPols(dr,r,basR);
  CompPols(ds,s,basS);

  // Compute bivariate basis functions.
  for(int j = 0; j < ds+1; ++j)
    for(int i = 0; i < dr+1; ++i)
      B[i+j*(dr+1)] = basR[i] * basS[j];

  // Release memory.
  delete [] basR;
  delete [] basS;
}

// -------------------------------------------------------------------------
// Bézier triangles:
//

// ============================= CompBezTrianNumBas ========================

int nBernsteinBasis :: CompTrianNumBas(const int &d)
{
  return (d+1)*(d+2)/2;
}

// ============================= MapTrianBaryToVecID =======================

int nBernsteinBasis :: MapTrianBaryToVecID(const int &i, const int &j, 
                                                         const int &d)
{
  return d - i - j + (d - i + 1) * (d - i)/2;
}

// ============================= IncTrianBaryID ============================

void nBernsteinBasis :: IncTrianBaryID(const int &d, int &i, int &j, int &k)
{
  if (j == 0)
  {
    --i;
    j = d - i;
    k = 0;
  }
  else
  {
   --j;
   ++k;
  }
}

// ========================= CompTrianPols =================================

void nBernsteinBasis :: CompTrianPols(int p, double r, double s, double *b)
{
  // Auxiliary Data
  int i, j, q, n;
  double t = 1 - r - s;

  b[0] = 1;
  for(int d = 1; d <= p; d++)
  {
    // Evaluate current number of basis functions.
    n = (d+1)*(d+2)/2;

    // Evaluate r edge.
    b[n-d-1] = s * b[n-2*d-1];
    for(q = n-d; q < n-1; ++q)
      b[q] = s * b[q-d] + t * b[q-d-1];
    b[n-1]     = t * b[n-d-2];

    for(q = n-2-d, i = 1; i < d; ++i, --q)
    {
      b[q] = r * b[q] + t * b[q-d+i-1];
        for(j = 1, --q; j < (d-i); ++j, --q)
          b[q] = r * b[q] + s * b[q-d+i] +
                  t * b[q-d+i-1];
      b[q] = r * b[q] + s * b[q-d+i];
    }
    b[0] = r * b[0];
  }
/*
  // Auxiliary Data
  int i,j,k;                       // Barycentric indices.
  double t = 1 - r - s;            // The third barycentric coord
  int size = CompTrianNumBas(p);   // Number of basis functions.
  double *aux = new double [size]; // Auxiliary vector.

  // Evaluate p - 1 basis.
  b[0] = aux[0] = 1;

  for(int deg = 1; deg <= p; deg++)
  {
    // Evaluate current number of basis functions.
    size = CompTrianNumBas(deg);

    // Initialize barycentric indices for the current degree.
    i = deg;   j = 0;   k = 0;

    // Evaluate each basis of the current degree.
    for(int id = 0; id < size; id++)
    {
      // Compute contribution in each parameter coordinate.
      b[id] = 0.0;
      if (i != 0) b[id] += r * aux[MapTrianBaryToVecID(i-1,j  ,deg-1)];
      if (j != 0) b[id] += s * aux[MapTrianBaryToVecID(i  ,j-1,deg-1)];
      if (k != 0) b[id] += t * aux[MapTrianBaryToVecID(i  ,j  ,deg-1)];

      // Increment berycentric indices.
      IncTrianBaryID(deg,i,j,k);
    }

    // Copy to the auxiliary vector
    for(int id = 0; id < size; id++)
      aux[id] = b[id];
  }

  delete [] aux;
  */
}

// ========================= CompTrianDervPols =============================

void nBernsteinBasis :: CompTrianDervPols(int p, double   r, double   s, 
                                                 double *dr, double *ds)
{
  // Auxiliary Data
  int i, j, q, n;
  double *b = ds;

  // Evaluate p-1 basis functions.
  CompTrianPols(p-1,r,s,b);

  // Evaluate current number of basis functions.
  n = (p+1)*(p+2)/2;

  // Evaluate r edge.
  dr[n-p-1] = 0.0;
  ds[n-p-1] = p * b[n-2*p-1];
  for(q = n-p; q < n-1; ++q)
  {
    dr[q] = -p * b[q-p-1];
    ds[q] =  p * (b[q-p] - b[q-p-1]);
  }
  dr[n-1]  = ds[n-1]  = -p * b[n-p-2];

  for(q = n-2-p, i = 1; i < p; ++i, --q)
  {
    dr[q] =  p * (b[q] - b[q-p+i-1]);
    ds[q] = -p * b[q-p+i-1];
    for(j = 1, --q; j < (p-i); ++j, --q)
    {
      dr[q] = p * (b[q] - b[q-p+i-1]);
      ds[q] = p * (b[q-p+i] - b[q-p+i-1]);
    }
    dr[q] = p * b[q];
    ds[q] = p * b[q-p+i];
  }
  dr[0] = p * b[0];
  ds[0] = 0;



/*

  // Auxiliary Data
  int i,j,k;                        // Barycentric indices.
  int sizep = CompTrianNumBas(p);   // Number of basis for degree p. 
  int sizeq = CompTrianNumBas(p-1); // Number of basis for degree q = p-1. 

  // Compute basis for p - 1 degree.
  double *b = new double [sizeq];
  CompTrianPols(p-1,r,s,b);

  // Initialize barycentric indices.
  i = p;   j = 0;    k = 0;

  // Evaluate each basis derivative.
  for(int id = 0; id < sizep; id++)
  {
    // Compute contribution in each parametric direction.
    dr[id] = 0.0;
    ds[id] = 0.0;

    if (i != 0) dr[id] += b[MapTrianBaryToVecID(i-1,j,p-1)];
    if (j != 0) ds[id] += b[MapTrianBaryToVecID(i,j-1,p-1)]; 
    if (k != 0)
    {
      dr[id] -= b[MapTrianBaryToVecID(i, j,p-1)];
      ds[id] -= b[MapTrianBaryToVecID(i, j,p-1)];
    }

    dr[id] *= p;
    ds[id] *= p;

    // Increment berycentric indices.
    IncTrianBaryID(p,i,j,k);
  }

  // Release memory.
  delete [] b;
  */
}

// ======================================================= End of file =====
