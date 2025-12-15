// -------------------------------------------------------------------------
// knotvec.h - file containing the definition of the cKnotVector class.
// -------------------------------------------------------------------------
// Copyright (c) 2017 LMCV/UFC
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
// The cKnotVector class store knot vector data and implements some algoritms
// related to them, as univariatee Bézier extraction operator. 
//
// -------------------------------------------------------------------------
// Private methods:
// -------------------------------------------------------------------------
//
// istream& ReadKnot(const eKnotIO &t, std::istream &in)
//
//   t  - given knot input/output rule.                            (in)
//   in - given input stream object.                               (in)
//
// This method read knot vector values from given input stream object.
// -------------------------------------------------------------------------
//
// istream& ReadKnot(const eKnotIO &t, std::istream &in)
//
//   t  - given knot input/output rule.                            (in)
//   in - given input stream object.                               (in)
//
// This method read knot vector multiplicuty from given input stream object.
//
// -------------------------------------------------------------------------
// Public methods:
// -------------------------------------------------------------------------
//
// cKnotVector(void)
//
// Construct an empty knot vector object.
// -------------------------------------------------------------------------
//
// cKnotVector(const cKnotVector &kv)
//
//   kv - given knot vector                                        (in)
//
// Construct a knot vector as a copy of given object.
//
// -------------------------------------------------------------------------
//
// void getSize(void)
//
// This method return the number of knot values in knot vector.
// -------------------------------------------------------------------------
//
// void getKnots(int &ks, double *&kv) 
// void getKnots(int &ks, float  *&kv) 
//
//   ks - given knot vector size.                                  (out)
//   kv - given knot vector values.                                (out)
//
// This method create a copy of current knot vector array and store its 
// pointer in kv and its size in ks. (The old content of kv is not deleted).
// -------------------------------------------------------------------------
//
// double& operator[](const int &i) const
//
//   id - Knot index                                               (in)
//
// This method return the reference of idth knot.
// -------------------------------------------------------------------------
//
// void operator=(const cKntVec &kv)
//
//   kv - given knot vector                                        (in)
//
// This method copy content of given knot vector.
// -------------------------------------------------------------------------
//
// void getEndKnots(double &init, double &end) 
//
//   init - Initial knot value.                                    (out)
//   end  - Final knot value.                                      (out)
//
// This method return the first and final knot values of this knot vector.
// -------------------------------------------------------------------------
//
// void getSpanVal(const int &sid, double &vi, double &vf) 
//
//   sid - Span index.                                             (in)
//   vi  - Initial knot value.                                     (out)
//   vf  - Final knot value.                                       (out)
//
// This method return the first and final knot values of the sidth span of this
// knot vector.
// -------------------------------------------------------------------------
//
// void getNumSpan(void) 
//
// This method return the number of knot spans.
// -------------------------------------------------------------------------
//
// int getKnotMult(const int &n, const int &p, const double &u) const
//
//   n - Number of basis in b-spline.                              (in)
//   p - Degree of b-spline.                                       (in)
//   u - parametric value used to find knot index.                 (in)
//
// This method return the multplicity of the knot span containing the given
// parametric value.
// -------------------------------------------------------------------------
//
// void getKnotMult(int* &mvec) const
//
//   mvec - Number of basis in b-spline.                           (out)
//
// This method return the multplicity of each knot span of this knot vector in
// mvec array. The old content mvec is not deleted.
// -------------------------------------------------------------------------
//
// void getSpanID(const int &kid) const
//
//   kid - knot vector index.                                      (out)
//
// This method return knot span index of given knot vector component.
// -------------------------------------------------------------------------
//
// void setKnots(const int &ks, double *kv) 
//
//   ks - given knot vector size.                                  (in)
//   kv - given knot vector values.                                (in)
//
// This method delete the current knot vector content and creates a new knots
// array as a copy of kv.
// -------------------------------------------------------------------------
//
// int FindSpan(const int &n, const int &p, const double &u) const
//
//   n - Number of basis in b-spline.                              (in)
//   p - Degree of b-spline.                                       (in)
//   u - parametric value used to find knot index.                 (in)
//
// This method return a knot vector index "i" satisfying the following
// condiction: KnotVal[i] < u <= KnotVal[i+1]. 
// -------------------------------------------------------------------------
//
// int FindSpan(const int &sid) const
//
//   sid - Span index.                                             (in)
//
// This method return a knot vector index of sidth knot span. A knot vector
// index "i" is considered at the begining of a knot span if it satisfying
// following condiction: KnotVal[i] < KnotVal[i+1].
// -------------------------------------------------------------------------
//
// int FindSpans(int &s, int *&svec) const
//
//   s    - Number of spans.                                       (in)
//   svec - Vector of knot spans.                                  (in)
//
// This method stores all knot spans of this knot vector in svec. The number of
// knot spans is stored in s. The old content of svec is not deleted.
// -------------------------------------------------------------------------
//
// void CompBezExtMat(int deg, double*** &C) const
//
//   deg - Degree of parent b-Spline.                              (in)
//   C   - Vector of knot spans.                                   (out)
//
// This method compute the univariate Bézier extraction matrix of each knot span
// and store in tensor C. The old content of C is not deleted.
//
// Ref: BORDEN, M. J.; SCOTT, M. A.; EVANS, J. A.; HUGHES, T. J. R. Isogeometric
// finite element data structures based on Bezier extraction of NURBS.
// International Journal for Numerical Methods in Engineering, v. 87, n. 1-5, p.
// 15–47, 2011.
//
// -------------------------------------------------------------------------
// Friend Methods:
// -------------------------------------------------------------------------
//
// istream& operator>>(istream &in, cKnotVector &kvec)
//
//   in   - given istream object                                   (in/out)
//   kvec - given knot vector object                               (in)
//
// This method read knot vector data from given input stream object.
// -------------------------------------------------------------------------
//
// ostream& operator<<(ostream &out, cKnotVector &kvec)
//
//   out  - given ostream object                                   (in/out)
//   kvec - given knot vector object                               (in/out)
//
// This method writes knot vector data in given output stream object.
// -------------------------------------------------------------------------

#ifndef _KNOTVEC_H
#define _KNOTVEC_H

// -------------------------------------------------------------------------
// Forward declarations:
//
#include <iosfwd>

// -------------------------------------------------------------------------
// Knot input/output rule.
//
enum eKnotIO
{
  KNOT_IO_GENERAL,
  KNOT_IO_UNIFORM
};

std::istream& operator>> (std::istream&,eKnotIO&);
std::ostream& operator<< (std::ostream&,eKnotIO&);

// -------------------------------------------------------------------------
// Definition of cKnotVector class:
//
class cKnotVector
{
 private:
  std::istream& ReadKnot(const eKnotIO&,std::istream&);
  std::istream& ReadMult(const eKnotIO&,std::istream&);

 protected:
  double        *Val;            // Knots.
  int            NumKnot;        // Number of knots.
  int            NumSpan;        // Number of knot spans.

 public:
            cKnotVector(void);
            cKnotVector(const cKnotVector&);
           ~cKnotVector(void);

  double&   operator[](const int &i) const { return Val[i];  }
  void      operator=(const cKnotVector&);
  int       getSize(void) const            { return NumKnot; }
  void      getKnots(int&,float*&);
  void      getKnots(int&,double*&);
  void      getEndKnots(double&,double&) const;
  void      getSpanVal(const int&,double&,double&) const;
  int       getNumSpan(void) const;
  int       getKnotMult(const int&,const int&,const double&) const;
  void      getKnotMult(int*&) const;
  int       getSpanID(const int&) const;
  void      setKnots(const int&,double*);
  int       FindSpan(const int&,const int&,const double&) const;
  int       FindSpan(const int&) const;
  void      FindSpans(int&,int*&) const;
  void      CompBezExtMat(int,double***&) const;

  friend  std::ostream&     operator<< (std::ostream&,const cKnotVector&);
  friend  std::istream&     operator>> (std::istream&,cKnotVector&);
};

#endif
