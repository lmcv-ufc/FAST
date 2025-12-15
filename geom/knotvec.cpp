// -------------------------------------------------------------------------
// knotvec.cpp - implementation of cKnotVec class.
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
// Created:      11-Oct-2017    Elias Saraiva Barroso
//
// Modified:     
// -------------------------------------------------------------------------

#include <iostream>
#include <cmath>

using namespace std;

#include "knotvec.h"

static double Tol = 1e-5;

// -------------------------------------------------------------------------
// Auxiliary functions:
//

// ======================== operator>> (eKnotIO) ===========================

istream& operator>> (istream &in, eKnotIO &type)
{
  // Read the patch type label.
  string label;
  if (!(in >> label))
  {
    cout << "Error in the input of the pat type label" << endl;
    in.setstate(ios::failbit);
    return in;
  }

  // Set the appropriate knot IO type.
  if (label == "'General'" || label == "'GENERAL'" || label == "'G'")
    type = KNOT_IO_GENERAL;
  else if (label == "'Uniform'" || label == "'UNIFORM'" || label == "'U'")
    type = KNOT_IO_UNIFORM;
  else
  {
    cout << "Unknown knot input/output type: " + label << "\n";
    in.setstate(ios::failbit);
    return in;
  }

  return in;
}

// ======================== operator<< (eKnotIO) ===========================

ostream& operator<< (ostream &out, eKnotIO &type)
{
  if (type == KNOT_IO_GENERAL)
    out << "'General'";
  else if (type == KNOT_IO_UNIFORM)
    out << "'Uniform'";

  return out;
}

// -------------------------------------------------------------------------
// Private methods:
//

// ============================= ReadKnot ==================================

std::istream& cKnotVector :: ReadKnot(const eKnotIO &t, std::istream &in)
{
  // Read knot vector size.
  if (!(in >> NumKnot))
  {
    cout << "Error in reading of knot vector size!";
    in.setstate(ios::failbit);
    return in;
  }

  // Read each knot.
  if (t == KNOT_IO_GENERAL)
  {
    Val = new double [NumKnot];
    for(int i = 0; i < NumKnot; ++i)
    {
      if (!(in >> Val[i]))
      {
        cout << "Error in reading of knot vector value!";
        in.setstate(ios::failbit);
        return in;
      }
    }  
  }

  else if (t == KNOT_IO_UNIFORM)
  {
    Val = new double [2];

    // Read the first knot value and the uniform increment.
    if (!(in >> Val[0]) || !(in >> Val[1]))
    {
      cout << "Error in reading of uniform knot vector data!";
      in.setstate(ios::failbit);
      return in;
    }
  }

  return in;
}

// ============================= ReadMult ==================================

std::istream& cKnotVector :: ReadMult(const eKnotIO &t, std::istream &in)
{
  int *mvec = new int [NumKnot];
  int ks = 0;

  // Read each multiplicity.
  for(int m = 0; m < NumKnot; ++m)
  {
    if (!(in >> mvec[m]))
    {
      cout << "Error in reading of knot multiplicity!";
      in.setstate(ios::failbit);
      return in;
    }

    ks += mvec[m];
  }

  // Create a knot vector with repeated values. 
  double *kv = new double [ks];

  int id = -1;  

  if (t == KNOT_IO_GENERAL)
    for(int knot = 0; knot < NumKnot; ++knot)
      for(int m = 0; m < mvec[knot]; ++m)
        kv[++id] = Val[knot];
  
  else if (t == KNOT_IO_UNIFORM)
    for(int knot = 0; knot < NumKnot; ++knot)
      for(int m = 0; m < mvec[knot]; ++m)
        kv[++id] = Val[0] + Val[1] * knot;


  // Release memory and stores knot vector.
  delete [] Val;
  delete [] mvec;

  NumKnot = ks;
  Val     = kv;

  // Compute the number of knot spans.
  NumSpan = 0;
  for(int i = 1; i < NumKnot; ++i)
    if (fabs(Val[i] - Val[i-1]) > 1e-6)
      ++NumSpan;

  return in;
}

// -------------------------------------------------------------------------
// Public methods:
//

// ============================= cKnotVector ===============================

cKnotVector :: cKnotVector(void)
{
  Val     = 0;
  NumKnot = 0;
  NumSpan = 0;
}

cKnotVector :: cKnotVector(const cKnotVector &kvec)
{
  Val     = 0;
  NumKnot = 0;
  NumSpan = 0;
  
  *this = kvec;
}

// ============================= ~cKnotVector ==============================

cKnotVector :: ~cKnotVector(void)
{
  delete [] Val;
}

// 

// ============================= operator= =================================

void cKnotVector :: operator=(const cKnotVector &kv)
{
  if (Val)
  { 
    delete Val;
    Val = 0;
  }

  if (kv.NumKnot > 0)
  {
    NumKnot = kv.NumKnot;
    NumSpan = kv.NumSpan;
    Val = new double [NumKnot];
    for(int i = 0; i < NumKnot; ++i)
      Val[i] = kv.Val[i];
  }
}

// ============================= getKnotMult ===============================

int cKnotVector :: getKnotMult(const int &n, const int &p, const double &u) const
{
  // Auxiliary data.
  int    mult = 0;

  // Find u span.
  int span = FindSpan(n,p,u);

  while (Tol > fabs(Val[span--] - u))
  {
    ++mult;
    if (span < 0) break;
  }

  return mult;
}

void cKnotVector :: getKnotMult(int* &mvec) const
{
  // Auxiliary data.
  int knot    = Val[0];
  int mult    = 1;
  int numknot = getNumSpan( ) + 1;
  mvec        = new int [numknot];

  int id = -1;
  for(int kid = 1; kid < NumKnot; ++kid)
  {
    if (Val[kid] - knot > Tol)
    {
      mvec[++id] = mult;
      mult = 1;
      knot = Val[kid];
    }
    else
      ++mult;
  }
  mvec[++id] = mult;
}

// ============================= FindSpan ==================================

int cKnotVector :: FindSpan(const int &n, const int &p, const double &u) const
{
  // Handle special case.
  if (fabs(u - Val[n]) < Tol) // (u == Val[n+1])
    return n-1;

  // Do binary search.
  int low = p, upp = n+1;
  int mid = (low+upp)/2;

  while (u < Val[mid] || u >= Val[mid+1])
  {
    if (u < Val[mid])
      upp = mid;
    else
      low = mid;
    mid = (low+upp)/2;
  }

  return mid;
}

int cKnotVector :: FindSpan(const int &sid) const
{
  // Do linear search.
  int currspan = -1;
  for(int id = 1; id < NumKnot; ++id)
    if ((Val[id] - Val[id-1]) > Tol)
      if (++currspan == sid)
        return (id - 1);

  return NumKnot-1;
}

// ============================= FindSpans =================================

void cKnotVector :: FindSpans(int &s, int *&svec) const
{
  s = getNumSpan( );
  svec = new int [s];

  int currid = 0;
  for(int i = 1; i < NumKnot; ++i)
    if (fabs(Val[i] - Val[i-1]) > 1e-6)
      svec[currid++] = i - 1;
}

// ============================= getSpanID =================================

int cKnotVector :: getSpanID(const int &kid) const
{
  // Do linear search.
  int spanID = 0;
  for(int id = 1; id < NumKnot; ++id)
  {
    if ((Val[id] - Val[id-1]) > Tol)
      ++spanID;

    if (id >= kid)
      break;
  }

  return spanID;
}

// ============================= getKnots ==================================

void cKnotVector :: getKnots(int &s, double *&knots)
{
  s = NumKnot;
  knots = new double [s];
  for(int i = 0; i < s; ++i)
    knots[i] = Val[i];
}

void cKnotVector :: getKnots(int &s, float *&knots)
{
  s = NumKnot;
  knots = new float [s];
  for(int i = 0; i < s; ++i)
    knots[i] = static_cast<float>(Val[i]);
}

// ============================= getEndKnots ===============================

void cKnotVector :: getEndKnots(double &init, double &end) const
{
  init = Val[0];
  end  = Val[NumKnot-1];
}

// ============================= getSpanVal ===============================

void cKnotVector :: getSpanVal(const int &sid, double &vi, double &vf) const
{
  int currspan = -1;
  for(int id = 1; id < NumKnot; ++id)
    if ((Val[id] - Val[id-1]) > Tol)
    {
      vi = Val[id-1];
      vf = Val[id];
      if (++currspan == sid) return;
    }
}

// ============================= getNumSpan ================================

int cKnotVector :: getNumSpan( ) const
{
  return NumSpan;
}

// ============================= FindSpan ==================================

void cKnotVector :: setKnots(const int &ks, double *knots)
{
  if (Val)
    delete [] Val;

  Val = new double [ks];
  for(int i = 0; i < ks; ++i)
    Val[i] = knots[i];

  // Compute the number of knot spans.
  NumSpan = 0;
  for(int i = 1; i < NumKnot; ++i)
    if (fabs(Val[i] - Val[i-1]) > 1e-6)
      ++NumSpan;

  NumKnot = ks;
}

// ============================= CompBezExtMat =============================

void cKnotVector :: CompBezExtMat(int deg, double***&C) const
{
  // Initializations.
  int p = deg;
  int m = NumKnot;
  int a = p + 1;
  int b = a + 1;
  int nb = 1;
  int ss = getNumSpan( );

  C = new double** [ss];
  for(int i = 0; i < ss; ++i)
  {
    C[i] = new double* [p+1];
    for(int j = 0; j < p+1; ++j)
    {
      C[i][j] = new double [p+1];
      for(int k = 0; k < p+1; ++k)
        C[i][j][k] = 0.0;
    }
  }

  for(int i = 0; i < p+1; i++)
    C[0][i][i] = 1.0;

  while (b < m)
  {
    // Initialize the next extraction operator.
    if (nb < ss)
      for(int i = 0; i < p+1; i++)
        C[nb][i][i] = 1;

    int i = b;

    // Count multiplicity of the knot at location b.
    while (b < m && Val[b] == Val[b-1])
      ++b;
    int mult = b-i+1;

    if (mult < p)
    {
      // Use Eq(10) to evaluate alphas.
      double numer = Val[b-1] - Val[a-1];
      double *alphas = new double [p-1];

      for (int j = p; j >= mult+1; j--)
        alphas[j-mult-1] = numer / (Val[a+j-1]-Val[a-1]);

      int r = p - mult;

      // Update matrix coefficients for r new knots.

      for(int j=1; j <=r; j++)
      {
        int save = r-j+1;
        int s = mult+j;

        for(int k = p+1; k >= s+1; k--)
        {
          double alpha = alphas[k-s-1];

          // The following line corresponds to Eq(9).
          for(int w = 1; w < p+1; w++)
            C[nb-1][w][k-1] = alpha*C[nb-1][w][k-1] + (1.0-alpha)*C[nb-1][w][k-2];
        }

        // Update overlapping coefficients of the next operator.
        if ((b < m) && (nb < ss))
          for(int w=save,v=p-j+1; w <= j+save; w++,v++)
            C[nb][w-1][save-1] = C[nb-1][v-1][p];
      }
    }
    ++nb; // Finished with the current operator.

    // Update Indices for next operator.
    if(b < m)
      a = b++;
  }
}

// ============================= operator>> ================================

istream& operator>>(istream &in, cKnotVector &kvec)
{
  // Read knot vector input type.
  eKnotIO type;

  if (!(in >> type))
  {
    cout << "Error in reading of knot vector input type!\n";
    in.setstate(ios::failbit);
    return in;
  }

  // Read knot values data.
  if (!(kvec.ReadKnot(type,in)))
  {
    cout << "Error in reading of knot vector values!\n";
    in.setstate(ios::failbit);
    return in;
  }

  // Read knot multiplicity.
  if (!(kvec.ReadMult(type,in)))
  {
    cout << "Error in reading of knot vector multiplicities!\n";
    in.setstate(ios::failbit);
    return in;
  }

  return in;
}

// ============================= operator<< ================================

ostream& operator<< (ostream &out,const cKnotVector &kvec)
{
  // Get all knot spans.
  int     ns; 
  int    *svec;
  kvec.FindSpans(ns,svec);

  // Write knot values
  out << ns + 1 << " ";
  for(int i = 0; i < ns; ++i)
    out << kvec.Val[svec[i]] << " ";
  out << kvec.Val[svec[ns-1]+1] << " ";

  // Compute knots multiplicities.
  int    *mult;
  kvec.getKnotMult(mult);

  // Write knot multiplicities.
  for(int i = 0; i < ns; ++i)
    out << mult[i] << " ";
  out << mult[ns] << " ";

  // Release memory.
  delete [] svec;
  delete [] mult;

  return out;
}

// ======================================================= End of file =====
