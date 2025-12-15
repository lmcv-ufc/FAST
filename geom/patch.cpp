// -------------------------------------------------------------------------
// patch.cpp - implementation of patch class.
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
// Created:      10-Mai-2015     Elias Saraiva Barroso
//
// -------------------------------------------------------------------------

#include <fstream>
#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <assert.h>
#include <cmath>

using namespace std;

#include "patch.h"
#include "bernbasis.h"
#include "bspbasis.h"
#include "knotvec.h"
#include "ctrlpnt.h"

// -------------------------------------------------------------------------
// Local functions:
//
static std::ostream& DefaultWriteCP(std::ostream&,const cPatch&);
static std::istream& DefaultReadCP(std::istream&,cPatch&);

// -------------------------------------------------------------------------
// Static variables:
//
PatWriteFunc cPatch :: WriteCtrlPntFunc = &DefaultWriteCP;
PatReadFunc  cPatch :: ReadCtrlPntFunc  = &DefaultReadCP;
sPatCell*    cPatch :: ListPatHead = 0;
sPatCell*    cPatch :: ListPatTail = 0; 
int          cPatch :: NumPatch    = 0;
int          cPatch :: NewPatLabel = 0;

//vector<cPatch*> cPatch :: VecPatch;
//vector<cNode*>  cPatch :: VecEqNode;
//vector<sPatPnt> cPatch :: VecDisplPts;
//vector<sPatPnt> cPatch :: VecStrPts;

static double Tol = 1e8;

// -------------------------------------------------------------------------
// Auxiliary methods:
//

// ======================== DefaultWriteCP =================================

std::ostream& DefaultWriteCP(std::ostream &out, const cPatch &pat)
{
  // Write the number of control points.
  out << pat.getNumCP( ) << " ";

  // Write each homogeneous coordinate point.
  for(int i = 0; i < pat.getNumCP( ); i++)
    out << *pat.getCtrlPnt(i) << " ";

  return out;
}

// ======================== DefaultReadCP ==================================

std::istream& DefaultReadCP(std::istream &in, cPatch &pat)
{
  // Read each homogeneus coordinate.
  int size         = pat.getNumBas( );
  sCtrlPnt** InpCP = new sCtrlPnt* [size];
  sCtrlPnt coord;
  for(int i = 0; i < size; ++i)
  {
    if (!(in >> coord))
    {
      cout << "Error in read of number of control points!" << endl;
      in.setstate(ios::failbit);
      return in;
    }

    InpCP[i] = new sCtrlPnt(coord);
  }
  pat.setCtrlPnt(InpCP,size);

  // Release memory
  delete [] InpCP;

  return in;
}

// ======================== operator>> (ePatType) ==========================

istream& operator>> (istream &in, ePatType &type)
{
  // Read the patch type label.
  string label;
  if (!(in >> label))
  {
    cout << "Error in the input of the pat type label" << endl;
    in.setstate(ios::failbit);
    return in;
  }

  // set the appropriate patch type.
  if (label == "'bspcurv'")
    type = BSPLINE_CURVE;
  else if (label == "'bspsurf'")
    type = BSPLINE_SURFACE;
  else if (label == "'bspsol'")
    type = BSPLINE_SOLID;
  else if (label == "'bezcurv'")
    type = BEZ_CURVE;
  else if (label == "'bezsurf'")
    type = BEZ_SURFACE;
  else if (label == "'beztrian'")
    type = BEZ_TRIAN;
  else
  {
    cout << "Unknown patch type: " + label << "\n";
    in.setstate(ios::failbit);
    return in;
  }

  return in;
}

// ======================== operator<< (ePatType) ==========================

ostream& operator<< (ostream &out, ePatType &type)
{
  // Write patch type.
  if (type == BSPLINE_CURVE)
    out << "'bspcurv'";
  else if (type == BSPLINE_SURFACE)
    out << "'bspsurf'";
  else if (type == BSPLINE_SOLID)
    out << "'bspsol'";
  else if (type == BEZ_CURVE)
    out << "'bezcurv'";
  else if (type == BEZ_SURFACE)
    out << "'bezsurf'";
  else if (type == BEZ_TRIAN)
    out << "'beztrian'";

  return out;
}

// -------------------------------------------------------------------------
// sPatCell structure:
//

// =============================== sPatCell ================================

sPatCell :: sPatCell( )
{ 
  nxt = 0; 
}

sPatCell :: sPatCell(cPatch *p)
{ 
  nxt = 0; 
  pat = p; 
}

// =============================== Turn ====================================

void sPatCell :: Turn(sPatCell* e)
{
  if (nxt) nxt->Turn(this); 
  nxt = e;
}

// -------------------------------------------------------------------------
// Private methods:
//

// =============================== PopBack =================================

void cPatch :: PopBack( )
{
  if (ListPatTail)
  {
    sPatCell *e = ListPatTail;
    if (ListPatHead == ListPatTail) ListPatHead = e->nxt;
    ListPatTail = e->nxt;
    delete e->pat;
    delete e;
    --NumPatch;
  }
}

void cPatch :: setCtrlPntWriteFunc(PatWriteFunc wfunc)
{
  WriteCtrlPntFunc = wfunc;
}

void cPatch :: setCtrlPntReadFunc(PatReadFunc rfunc)
{
  ReadCtrlPntFunc = rfunc;
}

// =============================== ReadPatch ===============================

void cPatch :: ReadPatch(std::istream &in)
{
  // Destroy current data.
  Destroy( );

  // Read the number of patches.
  int InpNumPat;
  if (!(in >> InpNumPat))
  {
    cout << "Error in the input of the number of patches!\n";
    exit(0);
  }
   
  int id;
  ePatType type;
  cPatch **vecpat = new cPatch* [InpNumPat];
  for(int i = 0; i < InpNumPat; ++i)
  {
    if (!(in >> id) || !(in >> type) || (id <= 0) || (id > InpNumPat))
    {
      cout << "Error in the input of patch label!\n";
      exit(0);
    }

    vecpat[id-1] = cPatch :: CreatePatch(type);
    vecpat[id-1]->setLabel(id);
    vecpat[id-1]->Read(in);
  }

  // Insert each created patch;
  for(int i = 0; i < InpNumPat; ++i)
    InsertPatch(vecpat[i]);

  // Release memory.
  delete [] vecpat;
}

// =============================== WritePatch ==============================

void cPatch :: WritePatch(std::ostream &out)
{
  if (!NumPatch) return;

  // Write the number of patches.
  out << NumPatch << "\n";

  // Turn linked-list to head-tail direction.
  if (ListPatTail) ListPatTail->Turn(0);

  // Re-assign patch label.
  int id = 0;
  for (sPatCell *c = ListPatHead; c; c = c->nxt)
    c->pat->setLabel(++id);

  // Write each patch.
  for (sPatCell *c = ListPatHead; c; c = c->nxt)
  {
    out << c->pat->getLabel( ) << " ";
    out << c->pat->type << " " << *(c->pat) << "\n";
  }
  out << "\n";

  // Turn linked-list back to tail-head direction.
  if (ListPatHead) ListPatHead->Turn(0);
}

// -------------------------------------------------------------------------
// Public methods:
//

// =============================== CreatePatch ===============================

cPatch* cPatch :: CreatePatch(ePatType t)
{
  cPatch *pat;

  switch (t)
  {
    case BEZ_TRIAN:
      pat = new cBezTriangle( );
    break;

    case BSPLINE_CURVE:
      pat = new cBspCurve( );
    break;

    case BSPLINE_SURFACE:
      pat = new cBspSurface( );
    break;
    
    case BSPLINE_SOLID:
      pat = new cBspSolid( );
    break;

    default:
      pat = 0;
    break;
  }

  return pat;
}

// =============================== InsertPatch =============================

void cPatch :: InsertPatch(cPatch *pat)
{
  // Create new cell in linked-list.
  sPatCell *e = new sPatCell(pat);
  ++NumPatch;
  pat->setLabel(++NewPatLabel);

  // Update Head-Tail pointers.
  if (!ListPatTail)
    ListPatTail = ListPatHead = e;
  else
  {
    e->nxt = ListPatTail;
    ListPatTail = e;
  }
}

// =============================== DeletePatch =============================

void cPatch :: DeletePatch(int label)
{
  if (!NumPatch)
    return;
    
  // Find element.
  sPatCell *e;
  for (e = ListPatTail; e ; e = e->nxt)
    if (e->pat->getLabel( ) == label)
      break;

  // Do not find cell with requested label.
  if (!e) return;

  // Move targed patch to end of list.
  if (e != ListPatTail)
  {
    e->pat->setLabel(ListPatTail->pat->getLabel( ));
    ListPatTail->pat->setLabel(label);

    swap(e,ListPatTail);
  }

  // Destroy the last element.
  PopBack( );
}

// =============================== Destroy =================================

void cPatch :: Destroy( )
{
  // Destroy global patch list.
  while (ListPatTail)
    PopBack( );
}

// ================================ FindPatch ==============================

cPatch* cPatch :: FindPatch(int label)
{
  for(sPatCell *e = ListPatTail; e; e = e->nxt)
    if (e->pat->getLabel( ) == label)
      return e->pat;
  
  return 0;
}

// =============================== cPatch ==================================

cPatch :: cPatch( )
{
  NumCP        = 0;
  CP           = 0;
  Rational     = false;
  RelCtrlPntMem = true;
}

// =============================== ~cPatch =================================

cPatch :: ~cPatch(void)
{
  // Release the array of control points.
  if (RelCtrlPntMem)
    for(int i = 0; i < NumCP; ++i)
      delete CP[i];

  if (CP) delete [] CP;
}

// =============================== getCtrlPnt ==============================

sCtrlPnt* cPatch :: getCtrlPnt(const int &id) const
{
  assert(id >= 0 && id < NumCP);
  return CP[id];
}

void cPatch :: getCtrlPnt(sCtrlPnt* &cpvec,  bool alloc) const
{
  if (alloc) cpvec = new sCtrlPnt[NumCP];

  for(int i = 0; i < NumCP; ++i)
    cpvec[i].copy(*CP[i]);
}

void cPatch :: getCtrlPnt(sCtrlPnt** &cpvec,  bool alloc) const
{
  if (alloc) cpvec = new sCtrlPnt*[NumCP];

  for(int i = 0; i < NumCP; ++i)
    cpvec[i] = CP[i];
}

// =============================== setCtrlPnt ==============================

void cPatch :: setCtrlPnt(const int &id, sCtrlPnt *cp, bool destroy)
{
  if (destroy)
    delete CP[id];

  CP[id] = cp;
}

void cPatch :: setCtrlPnt(sCtrlPnt** cpvec, const int &s, bool dest, bool inv)
{
  // Destroy current control points.
  if (dest)
    for(int i = 0; i < NumCP; ++i)
      delete CP[i];

  // Resize control point vector.
  if (NumCP != s)
  {
    delete [] CP;
    NumCP  = s;
    CP     = new sCtrlPnt*[NumCP];
  }

  // Assign new control points.
  if (inv)
    for(int i = 0; i < NumCP; ++i)
      CP[i] = cpvec[s -1 -i];
  else
    for(int i = 0; i < NumCP; ++i)
      CP[i] = cpvec[i];
}

void cPatch :: setCtrlPnt(sCtrlPnt *cp, const int &s, bool destroy, bool inv)
{
  // Copy input control points.
  sCtrlPnt **newcp = new sCtrlPnt* [s];
  for(int i = 0; i < s; ++i)
    newcp[i] = (new sCtrlPnt(cp[i]));

  // Assign control points.
  setCtrlPnt(newcp,s,destroy,inv);

  // Release memory.
  delete [] cp;
}

// =============================== getParVarLimits =========================

void cPatch :: getParVarLimits(const int &dim, double &inf, double &sup) const
{
  inf = 0.0;
  sup = 1.0;
}

// =============================== Evaluate ================================

sCtrlPnt cPatch :: Evaluate(double) const
{
  cout << "This patch object do not support this function!" << endl;
  exit(0);
  return sCtrlPnt( );
}

sCtrlPnt cPatch :: Evaluate(double,double) const
{
  cout << "This patch object do not support this function!" << endl;
  exit(0);
  return sCtrlPnt( );
}

sCtrlPnt cPatch :: Evaluate(double,double,double) const
{
  cout << "This patch object do not support this function!" << endl;
  exit(0);
  return sCtrlPnt( );
}

void cPatch :: getSpanLimits(int sid ,double ParVar[][2]) const
{
  for(int pvar = 0; pvar < NumParVar; ++pvar)
    getParVarLimits(pvar,ParVar[pvar][0],ParVar[pvar][1]);
}

int cPatch :: getSpanID(int spanR, int spanS, int spanT) const
{
  return 0;
}

int cPatch :: getSpanID(double r, double s, double t) const
{
  return 0;
}

void cPatch :: getSpanID(double r, double s, double t,int &sz, int i[]) const
{
  sz = 1;
  i[0] = 0;
}

void cPatch :: getSpanID(double r, double s, int &size, int inds[]) const
{
  getSpanID(r,s,0,size,inds);
}

void cPatch :: getSpanID(double r, int &size, int inds[]) const
{
  getSpanID(r,0,0,size,inds);
}

// -------------------------------------------------------------------------
// I/O operators:
//

// =============================== operator<< ==============================

std::istream& operator>>(std::istream &in, cPatch &pat)
{
  pat.Read(in);
  return in;
}

// =============================== operator>> ==============================

std::ostream& operator<<(std::ostream &out, const cPatch &pat)
{
  pat.Write(out);

  return out;
}

// -------------------------------------------------------------------------
// cBezCurvPat class:
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Public methods:
//

// =============================== cBezCurvPat ===============================

cBezCurve :: cBezCurve( ) : cPatch( )
{
  type = BEZ_CURVE;
}

cBezCurve :: cBezCurve(const int &deg, sCtrlPnt **cp) : cPatch( )
{
  Deg = deg;
  type = BEZ_CURVE;

  // Store control points references

  //int NumCP = Deg+1;
  //CP = new sCtrlPnt* [NumCP];
  //cout << "deg " << deg << endl;
//  for(int i = 0; i < NumCP; i++)
//  {
//  cout << "i " << i << endl;
//    CP[i] = cp[i];
//    cout << CP[i]->getLabel( ) << endl;
//    }
}

// =============================== ~cBezCurvPat ==============================

cBezCurve :: ~cBezCurve(void)
{
  // Destroy each basis
}

// =============================== Read ====================================

void cBezCurve :: Read(std::istream &in)
{
    /*
    cout << "Vai ler bspline"<< endl;
  // Write B-Spline basis.
  in >> *(bsp[0]);
    cout << *(bsp[0]) << endl;

  // Write control points.
  ReadCtrlPntFunc(in,*this);
  WriteCtrlPntFunc(cout,*this);
  */
}

// =============================== Write ===================================

void cBezCurve :: Write(std::ostream &out) const
{
  // Write B-Spline basis.
 // out << *(bsp[0]) << " ";

  // Write control points.
 // WriteCtrlPntFunc(out,*this);
}

// =============================== Read ====================================
/*
void cBezCurve :: Read(ifstream &in)
{
  // Reading rational flag.
  if (!(in >> Rational))
  {
    cout << "Error in the input of rational flag in Bézier curve!\n";
    exit(0);
  }

  // Reading degree
  if (!(in >> Deg))
    Utl::Exit("Error in the input of degree of bézier triangle patch.");

  // Alloc memory for Control Points

  NumCP = CompBezTrianNumBas(Deg);
  CP = new cNode*[NumCP];

  // Reading the Control Points
 
  int id;

  for (int i = 0; i < NumCP; i++)
  {
    in >> id;
    CP[i] = cNode :: FindNode(id);
  }
}
*/
// =============================== Evaluate ================================

sCtrlPnt cBezCurve :: Evaluate(double r) const
{
  // Evaluate basis function
  double *BasFunc;
  nBernsteinBasis :: CompPols(Deg,r,BasFunc);

  // Compute denominator
  double denom = 0.0;
  sCtrlPnt pos;

  if (Rational)
    for (int i = 0; i < NumCP; i++)
      denom += BasFunc[i] * CP[i]->getW( );

  assert(denom != 0.0);

  // Sum contribution from each control point

  if (Rational)
    for (int i = 0; i < NumCP; i++)
      pos += (BasFunc[i] * CP[i]->getW( ) / denom) * (*CP[i]);
  else
    for (int i = 0; i < NumCP; i++)
      pos += BasFunc[i] * (*CP[i]);

  return pos;
}

// =============================== DegreeElevation  ========================
/*    TODO: Passar para lib geom alg
void cBezCurve :: DegreeElevation(const int &di, std::vector<sCtrlPnt>& ncp)
{
  int newdeg;
  sCtrlPnt *aux;

  DegreeElevation(di,newdeg,aux);

  int size = newdeg + 1;

  assert(ncp.size( ) >= size);

  for(int i = 0; i < size; ++i)
    ncp[i].copy(aux[i]);

  delete [] aux;
}

void cBezCurve :: DegreeElevation(const int &di, int &nd, sCtrlPnt *&ncp)
{
  assert(di > 0);

  nd = di + Deg;         // Final degree.
  int fsize = nd + 1;    // Target size of control net.
  int csize = Deg + 1;   // Current size of control net.
  vector<sCtrlPnt> CurrCP; // Current control net.
  vector<sCtrlPnt> NewCP;  // New control net of elevated patch.
  double a;              // Auxiliary variable.

  // Resize containers to the final size.
  CurrCP.resize(fsize);
  NewCP.resize(fsize);

  // get the initial control points;
  getCtrlPnt(CurrCP);

  // Rational.
  if (Rational)
    for(int i = 0; i < csize; ++i)
      CurrCP[i] *= CurrCP[i].getW( );

  for(int deg = Deg+1; deg <= nd; ++deg)
  {
    // Compute the number of control points of the elevated curve.
    csize = deg + 1;

    // Store the last control points.
    NewCP[0]       = CurrCP[0];
    NewCP[csize-1] = CurrCP[csize-1];

    // Evaluate control net in current degree.
    for(int id = 1; id < (csize - 1); ++id)
    {
      a = double(id)/double(deg);
      NewCP[id] = a * CurrCP[id-1] + (1.0 - a) * CurrCP[id];
      NewCP[id].setW(a * CurrCP[id-1].getW( ) + (1.0 - a) * CurrCP[id].getW( ));
    }

    // Update current control net.
    for(int id = 1; id < (csize - 1); id++)
      CurrCP[id].copy(NewCP[id]);
  }
   // Rational.
   if (Rational)
     for(int i = 0; i < csize; ++i)
       CurrCP[i] /= CurrCP[i].getW( );

  ncp = new sCtrlPnt[csize];
  for(int id = 0; id < csize; ++id)
    ncp[id].copy(CurrCP[id]);
}
*/

// -------------------------------------------------------------------------
// cBezSurface class:
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Public methods:
//

// =============================== cBezSurface =============================

cBezSurface :: cBezSurface(void) : cPatch( )
{
  type   = BEZ_SURFACE;
  Deg[0] = 0;
  Deg[1] = 0;
}

cBezSurface :: cBezSurface(int dr, int ds, sCtrlPnt **cp) : cPatch( )
{
  Deg[0] = dr;
  Deg[1] = ds;
  type = BEZ_SURFACE;

  // Store control points references

  //int NumCP = Deg+1;
  //CP = new sCtrlPnt* [NumCP];
  //cout << "deg " << deg << endl;
//  for(int i = 0; i < NumCP; i++)
//  {
//  cout << "i " << i << endl;
//    CP[i] = cp[i];
//    cout << CP[i]->getLabel( ) << endl;
//    }
}

// =============================== ~cBezSurface ============================

cBezSurface :: ~cBezSurface(void)
{
  // Destroy each basis
}

// =============================== getNumBas ===============================

int cBezSurface :: getNumBas( ) const
{
  return (Deg[0] + 1) * (Deg[1] + 1);
}

// =============================== Read ====================================

void cBezSurface :: Read(std::istream &in)
{
    /*
    cout << "Vai ler bspline"<< endl;
  // Write B-Spline basis.
  in >> *(bsp[0]);
    cout << *(bsp[0]) << endl;

  // Write control points.
  ReadCtrlPntFunc(in,*this);
  WriteCtrlPntFunc(cout,*this);
  */
}

// =============================== Write ===================================

void cBezSurface :: Write(std::ostream &out) const
{
  // Write B-Spline basis.
 // out << *(bsp[0]) << " ";

  // Write control points.
 // WriteCtrlPntFunc(out,*this);
}

// =============================== Read ====================================
/*
void cBezCurve :: Read(ifstream &in)
{
  // Reading rational flag.
  if (!(in >> Rational))
  {
    cout << "Error in the input of rational flag in Bézier curve!\n";
    exit(0);
  }

  // Reading degree
  if (!(in >> Deg))
    Utl::Exit("Error in the input of degree of bézier triangle patch.");

  // Alloc memory for Control Points

  NumCP = CompBezTrianNumBas(Deg);
  CP = new cNode*[NumCP];

  // Reading the Control Points

  int id;

  for (int i = 0; i < NumCP; i++)
  {
    in >> id;
    CP[i] = cNode :: FindNode(id);
  }
}
*/

// =============================== Evaluate ================================

sCtrlPnt cBezSurface :: Evaluate(double r, double s) const
{
  // Evaluate basis function
//  double *BasFunc;
  //nBernsteinBasis :: CompBernsteinPols(Deg,r,BasFunc);

  // Compute denominator
  //double denom = 0.0;
  sCtrlPnt pos;
  return pos;
  /*

  if (Rational)
    for (int i = 0; i < NumCP; i++)
      denom += BasFunc[i] * CP[i]->getW( );

  assert(denom != 0.0);

  // Sum contribution from each control point

  if (Rational)
    for (int i = 0; i < NumCP; i++)
      pos += (BasFunc[i] * CP[i]->getW( ) / denom) * (*CP[i]);
  else
    for (int i = 0; i < NumCP; i++)
      pos += BasFunc[i] * (*CP[i]);

  return pos;
  */
}

// -------------------------------------------------------------------------
// cBezTrianPat class:
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Public methods:
//

// =============================== cBezTriangle =============================

cBezTriangle :: cBezTriangle(void)
{
  type = BEZ_TRIAN;
  Deg  = 0;
}

// =============================== ~cBezTriangle ============================

cBezTriangle :: ~cBezTriangle(void)
{

}

// =============================== getNumBas ===============================

int cBezTriangle :: getNumBas( ) const
{
  return ((Deg + 1) * (Deg + 2)) / 2;
}

// =============================== Read ====================================
/*
void cBezTrianPat :: Read(ifstream &in)
{
  // Reading rational flag.

  if (!(in >> Rational))
  {
    cout << "Error in the input of rational flag in Bézier triangle!\n";
    exit(0);
  }

  // Reading degree

  if (!(in >> Deg))
    Utl::Exit("Error in the input of degree of bézier triangle patch.");

  // Alloc memory for Control Points

  NumCP = CompBezTrianNumBas(Deg);
  CP = new cNode*[NumCP];

  // Reading the Control Points
 
  int id;

  for (int i = 0; i < NumCP; i++)
  {
    in >> id;
    CP[i] = cNode :: FindNode(id);
  }
}
*/
// =============================== Evaluate ================================

sCtrlPnt cBezTriangle :: Evaluate(double r, double s) const
{
  // Auxiliary data.
  int numcp = getNumBas( );                  // Number of control points.
  double w = 0.0;                            // Point weigth.
  sCtrlPnt pos;                              // Point position.

  // Evaluate basis function
  int     nb      = nBernsteinBasis :: CompTrianNumBas(Deg);
  double *BasFunc = new double [nb];
  nBernsteinBasis :: CompTrianPols(Deg,r,s,BasFunc);

  if (Rational)
  {
    for(int i = 0; i < numcp; ++i)
    {
      pos += BasFunc[i] * CP[i]->getW( ) * (*CP[i]);
      w   += BasFunc[i] * CP[i]->getW( );
    }
    pos.setW(w);
    pos /= w;
  }
  else
    for(int i = 0; i <= numcp; ++i)
      pos += BasFunc[i] * (*CP[i]);

  // Release memory.
  delete [] BasFunc;

  return pos;
}

// ===================== CompBernTrianSideCtrlPntID ========================

void cBezTriangle :: getSideCtrlPntID(int sc, int *inds) const
{
  assert(sc >= 0 && sc <= 2);
  // Initialize barycentric indices according with side code.
  int i = 0, j = 0;

  // Evaluate side control points.
  if (sc == 0)                         // side i -> j.
  {
    i = Deg;
    for(int id = 0; id < Deg + 1; ++id)
      inds[id] = nBernsteinBasis :: MapTrianBaryToVecID(i--,j++,Deg);
  }
  else if (sc == 1)                    // side j -> k.
  {
    j = Deg;
    for(int id = 0; id < Deg + 1; ++id)
      inds[id] = nBernsteinBasis :: MapTrianBaryToVecID(i,j--,Deg);
  }
  else if (sc == 2)                    // side k -> i.
    for(int id = 0; id < Deg + 1; ++id)
      inds[id] = nBernsteinBasis :: MapTrianBaryToVecID(i++,j,Deg);
}

// ===================== CompInnerCtrlPntID ================================

void cBezTriangle :: getInnerCtrlPntID(int &size, int *inds) const
{
  // Return null vector in the case of degree less then three.
  if (Deg < 3)
  {
    size = 0;
    inds = 0;
    return;
  }

  // Auxiliary data.
  int nl;             // Number of control net lines.
  int nil;            // Number of control net lines with inner cps.
  int id;             // Current inner control point identify.
  int i,j;            // Loop counters.

  // Compute number of control net lines.
  nl = Deg + 1;

  // Compute number of internal control points.
  nil = nl - 3;
  size = (nil * (nil+1))/2;
  inds = new int[size];

  j = 0;

  // Store each internal control points.
  for(int l = 2; l < nl - 1; ++l)
  {
    id = (l * (l+1))/2 + 1;
    for(i = 0; i < (l-2); ++i, ++id)
      inds[j++] = id;
  }
}


// =============================== DegreeElevation =========================
/*    TODO: Passar para lib geom alg
void cBezTriangle :: DegreeElevation(const int &di, std::vector<sCtrlPnt>& ncp)
{
  int newdeg;
  sCtrlPnt *aux;

  DegreeElevation(di,newdeg,aux);

  int size = nBernsteinBasis :: CompBezTrianNumBas(newdeg);

  assert(ncp.size( ) >= size);

  for(int i = 0; i < size; ++i)
    ncp[i].copy(aux[i]);

  delete [] aux;
}

void cBezTriangle :: DegreeElevation(const int &di, int &nd, sCtrlPnt *&ncp)
{
  assert(di > 0);

  int i,j,k;             // Barycentric indices.
  int vecid;             // Equivalent vector indice of barycentric indice.
  nd = di + Deg;         // Final degree.
  int fsize;             // Target size of control net.
  int csize;             // Current size of control net.
  vector<sCtrlPnt> CurrCP; // Current control net.
  vector<sCtrlPnt> NewCP;  // New control net of elevated patch.

  // Compute the initial and final sizes of the control net.
  csize = nBernsteinBasis :: CompBezTrianNumBas(Deg);
  fsize = nBernsteinBasis :: CompBezTrianNumBas(nd);

  // Resize containers to the final size.
  CurrCP.resize(fsize);
  NewCP.resize(fsize);

  // get the initial control points;
  getCtrlPnt(CurrCP);

  // Rational.
  if (Rational)
    for(int i = 0; i < csize; ++i)
      CurrCP[i] *= CurrCP[i].getW( );

  for(int deg = Deg+1; deg <= nd; ++deg)
  {
    // Compute the number of control points of the elevated patch.
    csize = nBernsteinBasis :: CompBezTrianNumBas(deg);

    // Initialize barycentric indices for the current degree.
    i = deg;   j = 0;    k = 0;

    // Evaluate control net in current degree.
    for(int id = 0; id < csize; id++)
    {
      // Reset current coordinate.
      NewCP[id] *= 0.0;   NewCP[id].setW(0.0);

      if (i != 0)
      {
        vecid      = nBernsteinBasis :: MapBaryToVecID(i-1,j,deg-1);
        NewCP[id] += i * CurrCP[vecid];
        NewCP[id].setW(NewCP[id].getW( ) + i * CurrCP[vecid].getW( ));
      }
      if (j != 0)
      {
        vecid = nBernsteinBasis::MapBaryToVecID(i,j-1,deg-1);
        NewCP[id] += j * CurrCP[vecid];
        NewCP[id].setW(NewCP[id].getW( ) + j * CurrCP[vecid].getW( ));
      }
      if (k != 0)
      {
        vecid = nBernsteinBasis::MapBaryToVecID(i,j,deg-1);
        NewCP[id] += k * CurrCP[vecid];
        NewCP[id].setW(NewCP[id].getW( ) + k * CurrCP[vecid].getW( ));
      }

      NewCP[id] *= (1.0/(deg));
      NewCP[id].setW(NewCP[id].getW( )*(1.0/(deg)));

      // Increment berycentric indices.
      nBernsteinBasis :: IncBezTrianBaryID(deg,i,j,k);
    }

    // Update current control net.
    for(int id = 0; id < csize; id++)
      CurrCP[id].copy(NewCP[id]);
  }
   // Rational.
   if (Rational)
     for(int i = 0; i < csize; ++i)
       CurrCP[i] /= CurrCP[i].getW( );

  ncp = new sCtrlPnt[csize];
  for(int id = 0; id < csize; ++id)
    ncp[id].copy(CurrCP[id]);
}
*/

// -------------------------------------------------------------------------
// cBspPat class:
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Public methods:
//

// =============================== cBspPat =================================

cBspPat :: cBspPat( )
{
  bsp          = 0;
  ExtMat       = 0;
}

// =============================== ~cBspPat ================================

cBspPat :: ~cBspPat(void)
{
  // Release the array of basis.
  delete [] bsp;

  // Release extraction matrices.
  ReleaseExtMat( );
}

// =============================== getBasis ================================

const cBSplineBasis& cBspPat :: getBasis(const int &bspId) const
{
  assert (bspId >= 0 && bspId < NumParVar);
  return *bsp[bspId];
}

// =============================== CPInd ===================================

int cBspPat :: getCPInd(int i, int j, int k) const
{
  int result = i;

  if (j != 0)
    result += bsp[0]->getNumBas( )*j;
  if (k != 0)
    result += (bsp[0]->getNumBas( ) * bsp[1]->getNumBas( ))*k;

  return result;
}

// =============================== setDegree ===============================

void cBspPat :: setDegree(const int &bspId, const int &deg)
{
  assert (bspId >= 0 && bspId < NumParVar);
  bsp[bspId]->setDegree(deg);
}

// =============================== getSpanLimits ===========================

void cBspPat :: getSpanLimits(int spanid, double SpanLim[][2]) const
{
  int aux;
  int Span[3] = {0,0,0};

  // Compute span index in each dimension.
  switch (NumParVar)
  {
    case (3):
      aux      = getBasis(0).getKnotVec( ).getNumSpan( ); 
      aux     *= getBasis(1).getKnotVec( ).getNumSpan( ); 
      Span[2]  = spanid / aux;
      spanid  -= aux * Span[2];

    case (2):
   
      aux      = getBasis(0).getKnotVec( ).getNumSpan( ); 
      Span[1]  = spanid / aux;
      spanid  -= aux * Span[1];

    case (1):
      Span[0]  = spanid;
  }
 
  for(int pvar = 0; pvar < NumParVar; ++pvar)
  {
    spanid = getBasis(pvar).getKnotVec( ).FindSpan(Span[pvar]);
    SpanLim[pvar][0] = getBasis(pvar).getKnotVec( )[spanid  ];
    SpanLim[pvar][1] = getBasis(pvar).getKnotVec( )[spanid+1];
  }
}

// =============================== getParVarLimits =========================

void cBspPat :: getParVarLimits(const int &bspId, double &inf, 
                                                  double &sup) const
{
  assert (bspId >= 0 && bspId < NumParVar);
  bsp[bspId]->getKnotVec( ).getEndKnots(inf,sup);
}

// =============================== getParVarLimits =========================

int cBspPat :: getNumSpan( ) const
{
  int size = 1;
  for(int i = 0; i < getNumParVar( ); ++i)
    size *= getBasis(i).getKnotVec( ).getNumSpan( );

  return size;
}

// =============================== getBasInd ===============================

void cBspPat :: getBasInd(const int &bspId, const int &span, int *bi) const
{
  assert (bspId >= 0 && bspId < NumParVar);
  bi[0] = bsp[bspId]->getBasInd(span);
  int p = bsp[bspId]->getDegree( );
  for(int i = 1; i < p + 1; ++i)
    bi[i] = bi[i-1] + 1;
}

int cBspPat :: getSpanID(double r, double s, double t) const
{
  int    kid;
  double parvals[3] = {r,s,t};
  int    currIDS[3] = {0,0,0};

  for(int pvar = 0; pvar < NumParVar; ++pvar)
  {
    kid           = getBasis(pvar).FindSpan(parvals[pvar]);
    currIDS[pvar] = getBasis(pvar).getKnotVec( ).getSpanID(kid);
  }

  return getSpanID(currIDS[0],currIDS[1],currIDS[2]);
}

// =============================== getSpanID ===============================

int  cBspPat :: getSpanID(int spanR, int spanS, int spanT) const
{

  int id = 0;
  switch (getNumParVar( ))
  {
    case 3:
      id += spanT * getBasis(0).getKnotVec( ).getNumSpan( )
                  * getBasis(1).getKnotVec( ).getNumSpan( );

    case 2:
      id += spanS * getBasis(0).getKnotVec( ).getNumSpan( );

    case 1:
      id += spanR;
  }

  return id;
}

void cBspPat :: getSpanID(double r, double s, double t,int &size, 
                                                       int inds[]) const
{
  double parvals[3] = {r,s,t};
  int SpanVec[][2] = {{-1,-1},{-1,-1},{-1,-1}}; 
  int SpanSize[3] = {1,1,1};
  int kid;

  size = 1;
  for(int pvar = 0; pvar < NumParVar; ++pvar)
  {
    kid              = getBasis(pvar).FindSpan(parvals[pvar]);
    SpanVec[pvar][0] = getBasis(pvar).getKnotVec( ).getSpanID(kid);

    // Verify if the input value is at boundary of two spans.
    if (fabs(getBasis(pvar).getKnotVec( )[kid] - parvals[pvar]) < Tol)
      if (SpanVec[pvar][0] < getBasis(pvar).getKnotVec( ).getNumSpan( ) - 1) 
        if (SpanVec[pvar][0] > 0) 
        {
          SpanVec[pvar][1]  = SpanVec[pvar][0] + 1;
          SpanSize[pvar]    = 2;
          size             *= 2;
        }
  }

  // Store each knot span index.
  int id = 0;
  int currIDS[3];

  for(int i = 0; i < SpanSize[0]; ++i)
  {
    currIDS[0] = currIDS[1] = currIDS[2] = 0;
    currIDS[0] = SpanVec[0][i];

    if (NumParVar == 1)
      inds[id++] = getSpanID(currIDS[0],currIDS[1],currIDS[2]);
    else
      for(int j = 0; j < SpanSize[1]; ++j)
      {
        currIDS[1] = SpanVec[1][j];

        if (NumParVar == 2)
          inds[id++] = getSpanID(currIDS[0],currIDS[1],currIDS[2]);
	else
          for(int k = 0; k < SpanSize[2]; ++k)
          {
            currIDS[2] = SpanVec[2][k];
            inds[id++] = getSpanID(currIDS[0],currIDS[1],currIDS[2]);
          }
      
      }
  }
}

// =============================== getExtMat ===============================

double** cBspPat :: getExtMat(const int &bind, const int &span) const 
{
  if (!ExtMat)
    return 0;

  assert (bind >= 0  && bind < NumParVar);
  assert (span  >= 0 && span  < ExtMat[bind].nspan);

  return ExtMat[bind].C[span];
}

// =============================== LoadExtMat ==============================

void cBspPat :: LoadExtMat( )
{
  // Release allocated extraction matrix.
  if (ExtMat)
    ReleaseExtMat( );

  // Alloc a vector of matrices for each parametric variable.
  ExtMat = new sExtMatData [NumParVar];

  // Loop over parametric variables and alloc each matrix set.
  for(int pvar = 0; pvar < NumParVar; ++pvar)
  {
    ExtMat[pvar].deg   = bsp[pvar]->getDegree( ); 
    ExtMat[pvar].nspan = bsp[pvar]->getKnotVec( ).getNumSpan( ); 

    // Load extraction matrices.
    bsp[pvar]->getKnotVec( ).CompBezExtMat(ExtMat[pvar].deg,ExtMat[pvar].C);
  }
}

// =============================== ReleaseExtMat ===========================

void cBspPat :: ReleaseExtMat( )
{
  // Return if no data is stored.
  if (!ExtMat)
    return;

  // Loop over parametric variable and release extraction matrices.
  for(int pvar = 0; pvar < NumParVar; ++pvar)
  {
    for(int i = 0; i < ExtMat[pvar].nspan; ++i)
    {
      for(int j = 0; j < ExtMat[pvar].deg + 1; ++j)
        delete [] ExtMat[pvar].C[i][j];
    
      delete [] ExtMat[pvar].C[i];
    }

    delete [] ExtMat[pvar].C;
  } 

  delete [] ExtMat;
  ExtMat = 0;
}

// =============================== Read ====================================

void cBspPat :: Read(istream &in)
{
  // Reading rational flag.
  if(!(in >> Rational))
  {
    cout << "Error in read of patch rational flag!\n";
    exit(0);
  }

  // Reading the bspline bases
  bsp = new cBSplineBasis* [NumParVar];

  for(int pvar = 0; pvar < NumParVar; ++pvar)
  {
    bsp[pvar] = new cBSplineBasis;
    if(!(in >> *(bsp[pvar])))
    {
      cout << "Error in read of B-Spline basis!\n";
      exit(0);
    }
  }

  // Read control points.
  ReadCtrlPntFunc(in,*this);
}


// =============================== getNumBas ===============================

int cBspPat :: getNumBas( ) const
{
  int nb = 1;
  for(int pvar = 0; pvar < NumParVar; ++pvar)
    nb *= bsp[pvar]->getNumBas( );

  return nb;
}

// -------------------------------------------------------------------------
// cBspCurve class:
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Public methods:
//

// =============================== cBspCurve ===============================

cBspCurve :: cBspCurve(void)
{
  Rational = false;
  type = BSPLINE_CURVE;
  NumParVar = 1;
  bsp = new cBSplineBasis* [1];
  bsp[0] = new cBSplineBasis;
}

cBspCurve :: cBspCurve(cBSplineBasis *b1)
{
  type = BSPLINE_CURVE;
  NumParVar = 1;
  bsp = new cBSplineBasis* [1];
  bsp[0] = b1;
}

// =============================== ~cBspCurve ==============================

cBspCurve :: ~cBspCurve(void)
{
  // Destroy basis

  delete bsp[0];  
  delete [] bsp; 
}

// =============================== getBasis ================================

const cBSplineBasis& cBspCurve :: getBasis( ) const
{
  return *bsp[0];
}

// =============================== getNumBas ===============================
int cBspCurve :: getNumBas( ) const
{
  return bsp[0]->getNumBas( );
}

// =============================== KnotInsertion ===========================
/*    TODO: Passar para lib geom alg
bool cBspCurve :: KnotInsertion(double u, int t, sCtrlPnt* &Qw, int &ps, double* &kv, int &ks)
{
  // Auxiliary data.
  int i, j, L;
  int mp, nq;
  int p  = getDegree( );
  int np = getNumBas( );
  double alpha;

  // get correct multiplicity of the insertion.
  int s = getBasis( )->getKnotMult(u);
  int r;

  if (s + t > getDegree( ))
    r = getDegree( ) - s;
  else
    r = t;

  // Abort the operation if the knot have maximum multiplicity.
  if (r <= 0) return false;

  // Create new knot vector data.
  ks = getNumBas( ) + getDegree( ) + 1 + r;
  kv = new double [ks];

  // Create new control point data.
          ps = getNumBas( ) + r;
          Qw = new sCtrlPnt [ps];
  sCtrlPnt *Rw = new sCtrlPnt [getDegree( ) + 1];

  // get knot span.
  int k = getBasis( )->FindSpan(u);

  mp = getNumBas( ) + getDegree( ) + 1;
  nq = getNumBas( ) + r;

  // Load new knot vector.
  for(i = 0  ; i <= k ; ++i)
      kv[i]   = getBasis( )->getKnot(i);
  for(i = 1  ; i <= r ; ++i)
    kv[k+i] = u;
  for(i = k+1; i < mp; ++i)
    kv[i+r] = getBasis( )->getKnot(i);

  // Save unaltered control points.
  for(i = 0  ; i <= k - p ; ++i) Qw[i].copy(*CP[i]);
  for(i = k-s; i <  np ; ++i)    Qw[i+r].copy(*CP[i]);
  for(i = 0  ; i <= p-s ; ++i)   Rw[i].copy(*(CP[k-p+i]));

  // Project control points into 4th dimension in case of rational curve.
  if (Rational)
  {
    for(i = 0; i < ps; ++i) Qw[i]  *= Qw[i].getW( );
    for(i = 0; i < p+1; ++i) Rw[i] *= Rw[i].getW( );
  }

  // Insert the knot r times.
  for(j = 1  ; j <= r ; ++j)
  {
    L = k - p + j;
    for(i = 0; i <= p-j-s; ++i)
    {
      alpha = (u - getKnot(L+i))/(getKnot(i+k+1) - getKnot(L+i));
      Rw[i] = alpha * Rw[i+1] + (1.0-alpha) * Rw[i];
      if (Rational)
        Rw[i].setW(alpha * Rw[i+1].getW( ) + (1.0-alpha)*Rw[i].getW( ));
    }

    Qw[L].copy(Rw[0]);
    Qw[k+r-j-s].copy(Rw[p-j-s]);
  }

  // Load remaining control points.
  for(i=L+1; i<k-s; ++i)
    Qw[i].copy(Rw[i-L]);

  // Project control points back to 3th dimension.
  if (Rational)
    for(i = 0; i < ps; ++i) Qw[i] /= Qw[i].getW( );

  return true;
}

// =============================== DegreeElevation =========================

bool cBspCurve :: DegreeElevation(int di, int &nks, double* &nkv, int &s,sCtrlPnt **&cp)
{
  assert(di > 0);

  int p = bsp[0]->getDegree( );
  int deg;
  int i,j,k,m;  // Loop iterators.
  double coeff;
  double w;
  sCtrlPnt pnt1, pnt2;

  int ks;
  double *kv;

  int nc, cs;
  sCtrlPnt ***BezCP, **NBCP;
  int fdeg = this->getDegree( ) + di;

  // get initial Bezier segments.
  CompBezCtrlPnt(nc,cs,BezCP);

  // Process Degree elevation in each Bezier segment.
  NBCP = new sCtrlPnt* [nc];

  for(i = 0; i < nc; ++i)
  {
    // Alloc memory for a new bezier curve.
    NBCP[i] = new sCtrlPnt[fdeg+1];

    // Load current bezier control polygon.
    for(j = 0; j < cs; ++j)
    {
      NBCP[i][j].copy(*BezCP[i][j]);
      NBCP[i][j] *= NBCP[i][j].getW( );
    }

    // Do degree elevation.
    for(deg = bsp[0]->getDegree( ) + 1; deg <= fdeg; ++deg)
    {
      coeff = 1.0/deg;
      pnt2.copy(NBCP[i][0]);
      for(k = 1; k < deg+1; ++k)
      {
        pnt1.copy(NBCP[i][k]);
        NBCP[i][k] = coeff * k * pnt2 + (1.0-coeff*k) * NBCP[i][k];
        w          = coeff * k * pnt2.getW( ) + (1.0-coeff*k) * NBCP[i][k].getW( );
        NBCP[i][k].setW(w);
        pnt2.copy(pnt1);
      }
    }
  }

  // Modify current b-spline basis data.
  bsp[0]->getKnots(ks,kv);
  nks = ks + (bsp[0]->getNumKnotSpan( ) + 1) * di;
  nkv = new double [nks];
  for(i = 1, j = 0; i < ks; i++)
  {
    nkv[j++] = kv[i-1];
    if ((kv[i] - kv[i-1]) > 1e-8)
      for(k = 0; k < di; ++k)
        nkv[j++] = kv[i-1];
  }
  for(k = 0; k < di+1; ++k)
    nkv[j++] = kv[i-1];

  // Modify current b-spline basis data.
  bsp[0]->setDegree(fdeg);
  bsp[0]->setKnots(nks,nkv);

  // Compute new extraction matrix.
  double ***C;
  CompBezExtMat(C);

  // Transpose each matrix and process LU decomposition.
  for(i = 0; i < nc; ++i)
  {
    nMatLib :: SqrMatTranspose(fdeg+1,C[i]);
    if (!nMatLib :: DecompLU(fdeg+1,C[i]))
    {
      cout << "Singular matrix found at cBspCurve :: DegreeElevation!" << endl;
      return (0);
    }
  }

  // Alloc memory for the new control polygon.
  s  = bsp[0]->getNumBas( );
  cp = new sCtrlPnt*[s];
  for(i = 0; i < s; ++i)
    cp[i] = 0;

  sCtrlPnt *localcp = new sCtrlPnt[fdeg+1];

  // get each knot span.
  int spansize;
  int *svec;
  getBasis( )->getAllKnotSpan(spansize,svec);

  // Loop over each bezier segment, evaluating NURBS b-spline control polygon.
  for(int i = 0; i < nc; ++i)
  {
    if (cp[svec[i] - p] && (i < nc-1))
      continue;

      // Clear local bsp control polygon.
      for(j = 0; j < fdeg+1; ++j)
        localcp[j].zero( );

      // Top-down substitution to find {y} in [L] {y} = {b}.
      for(j = 0; j < fdeg+1; ++j)
      {
        pnt1.zero( ); w = 0;
        for(k = 0; k < j; ++k)
        {
          pnt1 += C[i][j][k] * localcp[k];
          w    += C[i][j][k] * localcp[k].getW( );
        }
        localcp[j] = NBCP[i][j] - pnt1;
        localcp[j].setW(NBCP[i][j].getW( ) - w);
      }

      // Bottom-up substitution to find {x} in [U] {x} = {y}.
      for(j = fdeg; j > -1; --j)
      {
        pnt1.zero( ); w = 0;
        for(k = fdeg; k > j; --k)
        {
          pnt1 += C[i][j][k] * localcp[k];
          w    += C[i][j][k] * localcp[k].getW( );
        }
        localcp[j] = (localcp[j] - pnt1) / C[i][j][j];
        localcp[j].setW((localcp[j].getW( ) - w) / C[i][j][j]);
      }

      // Store nurbs control point.
      for(int k = 0; k < fdeg+1; ++k)
      {
        if (!cp[svec[i]-fdeg+k])
        {
          cp[svec[i]-fdeg+k] = new sCtrlPnt;
          cp[svec[i]-fdeg+k]->copy(localcp[k]);
          *cp[svec[i]-fdeg+k] /= cp[svec[i]-fdeg+k]->getW( );
        }
      }
  }

  // Reset old knot vector and degree.
  bsp[0]->setDegree(fdeg - di);
  bsp[0]->setKnots(ks,kv);

  // Release memory.
  delete [] kv;
  delete [] localcp;

  for(int i = 0; i < nc; ++i)
  {
    for(int j = 0; j < (cs+di); ++j)
      delete [] C[i][j];

    delete [] C[i];
  }
  delete [] C;


  delete BezCP[0][0];
  for(int i = 0; i < nc; ++i)
  {
    for(int j = 1; j < cs; ++j) // j = 0 was deleted in
      delete BezCP[i][j];       // the last iteration.

    delete [] BezCP[i];
    delete [] NBCP[i];
  }
  delete [] BezCP;
  delete [] NBCP;

  return true;
}
*/

// =============================== Read ====================================
/*
void cBspCurve :: Read(std::istream &in)
{
  // Write B-Spline basis.
  in >> Rational;
  in >> *(bsp[0]);

  // Write control points.
  ReadCtrlPntFunc(in,*this);
  WriteCtrlPntFunc(cout,*this);
}
*/
// =============================== Write ===================================

void cBspCurve :: Write(std::ostream &out) const
{
  // Write B-Spline basis.
  out << Rational << " ";
  out << *(bsp[0]) << " ";

  // Write control points.
  WriteCtrlPntFunc(out,*this);
}

// =============================== Evaluate ================================

sCtrlPnt cBspCurve :: Evaluate(double r) const
{
  // Auxiliary data.
  int span = bsp[0]->FindSpan(r);                // Current knot span.
  int deg  = bsp[0]->getDegree( );               // Curve degree.
  double w = 0.0;                                // Point weigth.
  sCtrlPnt pos;                                  // Point position.

  // Evaluate basis functions.
  double *BasFunc = new double [deg+1];
  bsp[0]->CompBasFunc(r,BasFunc);

  if (Rational)
  {
    for(int i = 0; i <= deg; ++i)
    {
      pos += (BasFunc[i] * CP[span-deg+i]->getW( )) * (*CP[span-deg+i]);
      w   += BasFunc[i] * CP[span-deg+i]->getW( );
    }
    pos.setW(w);
    pos /= w;
  }
  else
    for(int i = 0; i <= deg; ++i)
      pos += BasFunc[i] * (*CP[span-deg+i]);

  // Release memory.
  delete [] BasFunc;

  return pos;
}

// =============================== CompDervFD ==============================
/*    TODO: Passar para lib geom alg
sCtrlPnt cBspCurve :: CompDervFD(double r, double dr)
{
  // Compute curve derivative at r using finite diferences method.
  sCtrlPnt p1,p2,res;

  p1 = Evaluate(r-dr/2.0);
  p2 = Evaluate(r+dr/2.0);

  res.setX((p2.getX( ) - p1.getX( ))/dr);
  res.setY((p2.getY( ) - p1.getY( ))/dr);

  return res;
}

// =============================== Length ==================================


double cBspCurve :: Length(double t1, double t2)
{
  using namespace nMatLib;

  double tol = 1e-4;
  double currlength, length = -1;

  int maxpnt = nQuadTable :: getSize(GAUSS);
  int numpnts;

  sCtrlPnt derv;
  double intpnt, w;
  double dx2, dy2;
  double jacob;


  cout << "vai avaliar " << t1 << " " << t2 << endl;

  for(int order = 3; order < maxpnt; ++order)
  {
    // Evaluate current length.
    numpnts = order;
    currlength = 0;
    for(int pntid = 0; pntid < numpnts; ++pntid)
    {
      nQuadTable :: getParVal(GAUSS,order,pntid,intpnt);
      nQuadTable :: getWeight(GAUSS,order,pntid,w);

      // Write gauss point into t1-t2 interval.
      intpnt = (intpnt + 1.0)/2.0 * (t2 - t1) + t1;

      derv  = CompDervFD(intpnt);
      dx2   = derv.getX( ) * derv.getX( );
      dy2   = derv.getY( ) * derv.getY( );
      jacob = (t2-t1)/2.0;
      currlength += sqrt(dx2 + dy2) * w * jacob;
    }

    if (order > 3)
      if (fabs((length - currlength)/length) < tol)
        break;

    length = currlength;
  }

  return length;
}
*/

// =============================== CompBezCtrlPnt ==========================
/*    TODO: Passar para lib geom alg
void cBspCurve :: CompBezCtrlPnt(int &nc, int &cs, sCtrlPnt*** &cpvec)
{
  // setting the number of curves and control points in each.
  int p  = getDegree( );
      nc = getNumSpan( );
      cs = p + 1;

  // get each knot span.
  int s;
  int *svec;
  getBasis( )->getAllKnotSpan(s,svec);
  assert(s == nc);

  // Alloc memory for control points table.
  cpvec = new sCtrlPnt** [nc];
  for(int i = 0; i < nc; ++i)
    cpvec[i] = new sCtrlPnt* [cs];

  // Evaluate Bezier extraction operator (matrix).
  double ***C;
  CompBezExtMat(C);

  // Project each NURBS control points into 4th dimension.
  if (Rational)
    for(int i = 0; i < NumCP; ++i)
      (*CP[i]) *= CP[i]->getW( );       // Pb = (1/Wb) * Ct * Wn * Pn

  // Evaluate each bézier point;
  for(int i = 0; i < nc; ++i)
  {
    for(int j = 0; j < cs; ++j)
    {
      sCtrlPnt coord;
      if (i != 0 && j == 0)
      {
        cpvec[i][j] = cpvec[i-1][cs-1];
        coord.copy(*cpvec[i][j]);
      }
      else
      {
        if (Rational)
        {
          coord *= 0.0;
          coord.setW(0.0);
        }
        for(int k = 0; k < cs; ++k)
        {
          assert(i < s);
          assert(svec[i] - p + k < getNumCP( ));

          coord += C[i][k][j] * (*CP[svec[i] - p + k]);
          if (Rational)
            coord.setW(coord.getW( ) + C[i][k][j] * CP[svec[i]-p+k]->getW( ));
        }

        // Project control point back to 3th dimension.
        if (Rational)
          coord /= coord.getW( );

        cpvec[i][j] = new sCtrlPnt( );
        cpvec[i][j]->copy(coord);
      }
    }
  }

  // Project each B-Spline control points to 3th dimension.
  if (Rational)
    for(int i = 0; i < NumCP; ++i)
      (*CP[i]) /= CP[i]->getW( );

  // Release memory.
  for(int i = 0; i < nc; ++i)
  {
    for(int j = 0; j < cs; ++j)
      delete [] C[i][j];

    delete [] C[i];
  }

  delete [] C;
  delete [] svec;
}
*/

// -------------------------------------------------------------------------
// cBspSurface class:
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Public methods:
//

// =============================== cBspSurface =============================

cBspSurface :: cBspSurface(void)
{
  type = BSPLINE_SURFACE;
  NumParVar = 2;
}

cBspSurface :: cBspSurface(cBSplineBasis *b1, cBSplineBasis *b2)
{
  type = BSPLINE_SURFACE;
  bsp = new cBSplineBasis* [2];
  bsp[0] = b1;
  bsp[1] = b2;
  NumParVar = 2;
}

// =============================== ~cBspSurface ============================

cBspSurface :: ~cBspSurface(void)
{
  // Destroy each basis

  delete bsp[0];
  delete bsp[1];
}

/*
// =============================== Evaluate ================================

sNodeCoord cBspSurface :: Evaluate(double r, double s, bool deform)
{
  // Compute denominator

  double denom = 0.0, aux;

  if (Rational)
    for (int i = 0; i < bsp[1]->getNumBas( ); i++)
      for (int j = 0; j < bsp[0]->getNumBas( ); j++)
      {
        aux = bsp[1]->getBasValue(i,s)*bsp[0]->getBasValue(j,r);
        denom += aux * CP[getCPInd(j,i)]->getW( );
      }
  else
    denom = 1.0;

  sNodeCoord coord = {0.0, 0.0, 0.0};

  if (denom == 0)
    return coord;

  // Compute numerator

  double resX = 0, resY = 0, resZ = 0, val;

  for (int i = 0; i < bsp[1]->getNumBas( ); i++)
    for (int j = 0; j < bsp[0]->getNumBas( ); j++)
    {
      coord = CP[getCPInd(j,i)]->getCoord( );

      if (deform)
      {
        coord.x += CP[getCPInd(j,i)]->getDispl(0);
        coord.y += CP[getCPInd(j,i)]->getDispl(1);
        coord.z += CP[getCPInd(j,i)]->getDispl(2);
      }
        
      val = bsp[1]->getBasValue(i,s)*bsp[0]->getBasValue(j,r);

      if (Rational)
      {
        resX += val*CP[getCPInd(j,i)]->getW( )*coord.x/denom;
        resY += val*CP[getCPInd(j,i)]->getW( )*coord.y/denom;
        resZ += val*CP[getCPInd(j,i)]->getW( )*coord.z/denom;
      }
      else
      {
        resX += val*coord.x;
        resY += val*coord.y;
        resZ += val*coord.z;
      }
    }

  coord.x = resX;
  coord.y = resY;
  coord.z = resZ;

  return coord;
}
*/

// -------------------------------------------------------------------------
// cBspSolid class:
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Public methods:
//

// =============================== cBspSolid ===============================

cBspSolid :: cBspSolid(void)
{
  type = BSPLINE_SOLID;
  NumParVar = 3;
  bsp = new cBSplineBasis* [3];
  bsp[0] = new cBSplineBasis;
  bsp[1] = new cBSplineBasis;
  bsp[2] = new cBSplineBasis;
}

cBspSolid :: cBspSolid(cBSplineBasis *b1, cBSplineBasis *b2, cBSplineBasis *b3)
{
  type = BSPLINE_SOLID;
  NumParVar = 3;
  bsp = new cBSplineBasis* [3];
  bsp[0] = b1;
  bsp[1] = b2;
  bsp[2] = b3;
}

// =============================== ~cBspSolid ==============================

cBspSolid :: ~cBspSolid(void)
{
  // Destroy each basis

  delete bsp[0];
  delete bsp[1];
  delete bsp[2];
}


/*
// =============================== Evaluate ================================

sNodeCoord cBspSolid :: Evaluate(double r, double s, double t, bool deform)
{
   // Compute denominator

  double denom = 0.0, aux;

  if (Rational)
    for (int i = 0; i < bsp[2]->getNumBas( ); i++)
      for (int j = 0; j < bsp[1]->getNumBas( ); j++)
        for (int k = 0; k < bsp[0]->getNumBas( ); k++)
        {
          aux = bsp[0]->getBasValue(k,r);
          aux *= bsp[1]->getBasValue(j,s);
          aux *= bsp[2]->getBasValue(i,t);
    
          denom += aux * CP[getCPInd(k,j,i)]->getW( );
        } 
  else
    denom = 1.0;

  sNodeCoord coord = {0.0,0.0,0.0};

  if (denom == 0.0)
    return coord;

  // Compute numerator

  double resX = 0, resY = 0, resZ = 0, val;

  for (int i = 0; i < bsp[2]->getNumBas( ); i++)
    for (int j = 0; j < bsp[1]->getNumBas( ); j++)
      for (int k = 0; k < bsp[0]->getNumBas( ); k++)
      {
        coord = CP[getCPInd(j,i)]->getCoord( );
      
        if (deform)
        {
          coord.x += CP[getCPInd(k,j,i)]->getDispl(0);
          coord.y += CP[getCPInd(k,j,i)]->getDispl(1);
          coord.z += CP[getCPInd(k,j,i)]->getDispl(2);
        }
          
        val  = bsp[0]->getBasValue(k,r);
        val *= bsp[1]->getBasValue(j,s);
        val *= bsp[2]->getBasValue(i,t);
      
        if (Rational)
        {
          resX += val*CP[getCPInd(k,j,i)]->getW( )*coord.x/denom;
          resY += val*CP[getCPInd(k,j,i)]->getW( )*coord.y/denom;
          resZ += val*CP[getCPInd(k,j,i)]->getW( )*coord.z/denom;
        }
        else
        {
          resX += val*coord.x;
          resY += val*coord.y;
          resZ += val*coord.z;
        }
      }

  coord.x = resX;
  coord.y = resY;
  coord.z = resZ;

  return coord;
}
*/
// ======================================================= End of file =====
