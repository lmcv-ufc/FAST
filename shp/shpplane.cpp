// -------------------------------------------------------------------------
// shpplane.cpp - implementation of the Shape class.
// -------------------------------------------------------------------------
// Created:      06-Sep-2005     Aurea Silva de Holanda
//
// Modified:     17-May-2011     Evandro Parente Junior
//               Implementation of non-isoparametric elements.
//
// Modified:     14-Nov-2023     Renan Melo Barros/Pedro Ygor Mesquita
//               Computation of second derivatives.
// -------------------------------------------------------------------------

#include <stdlib.h>
#include <iostream>
#include <iomanip>

using namespace std;

#include "shpplane.h"
#include "gblvar.h"
#include "ctrl.h"
#include "vec.h"

#ifdef _OMP_
#include "omp.h"
#endif

// -------------------------------------------------------------------------
// Static variables:
//
static double        _J[3][3]; // Jacobian matrix
static double     _InvJ[3][3]; // Inverse of Jacobian matrix
static double       _Ja[3][3]; // Auxiliar matrix to compute second derivs
static double    _InvJa[3][3]; // Inverse of auxiliar matrix
static sNatCoord _dNrst[500];  // Shape func derivatives/parametric coords
static sNatCoord _dMrst[500];  // Mapping func derivatives/parametric coords
static sNatDrv   _ddNrs[500];  // Shape functions second derivatives.

#pragma omp threadprivate(_J,_InvJ,_dNrst,_dMrst,_A,_InvA,_ddNrs)

// -------------------------------------------------------------------------
// Public methods:
//

// ============================== cShapePlane ==============================

cShapePlane :: cShapePlane(void)
{
}

// ============================= ~cShapePlane ==============================

cShapePlane :: ~cShapePlane(void)
{
}

// =============================== DrvShpXYZ ===============================

void cShapePlane :: DrvShpXYZ(sNatCoord pt, sNodeCoord *coord,
                              double *detJ, sNodeCoord *dNxyz)
{
  int i;

  // Compute dNrs.

  DrvShpRST(pt, _dNrst);

  // Compute [J].

  ShpMatZero(2, _J);
  if (Isoparametric( )) // Isoparametric mapping
  {
    for (i = 0; i < NumNode; i++)
    {
      _J[0][0] += _dNrst[i].r*coord[i].x;
      _J[0][1] += _dNrst[i].r*coord[i].y;
      _J[1][0] += _dNrst[i].s*coord[i].x;
      _J[1][1] += _dNrst[i].s*coord[i].y;
    }
  }
  else                  // Non-isoparametric mapping
  {
    DrvMapRST(pt, _dMrst);
    for (i = 0; i < NumNode; i++)
    {
      _J[0][0] += _dMrst[i].r*coord[i].x;
      _J[0][1] += _dMrst[i].r*coord[i].y;
      _J[1][0] += _dMrst[i].s*coord[i].x;
      _J[1][1] += _dMrst[i].s*coord[i].y;
    }
  }

  // Compute |J| and inv(J).

  ShpMatInv(2, _J, _InvJ, detJ);

  // Compute Ni,xyz.

  for (i = 0; i < NumNode; i++)
  {
    dNxyz[i].x = _InvJ[0][0]*_dNrst[i].r + _InvJ[0][1]*_dNrst[i].s;
    dNxyz[i].y = _InvJ[1][0]*_dNrst[i].r + _InvJ[1][1]*_dNrst[i].s;
  }
}

// =============================== DrvShpXYZ ===============================

void cShapePlane :: DrvShpXYZ(sNatCoord pt, sNodeCoord *coord, double *detJ, 
                              sNodeCoord *dNxyz, sNodeDrv *ddNxy)
{
  // Evaluate first derivatives, Jacobian matrix, and detJ.

  DrvShpXYZ(pt, coord, detJ, dNxyz);

  // Evaluate auxiliary matrix [Ja].

  _Ja[0][0] = _J[0][0] * _J[0][0];                        // x,r^2
  _Ja[0][1] = _J[0][1] * _J[0][1];                        // y,r^2
  _Ja[0][2] = 2.0 * (_J[0][0] * _J[0][1]);                // 2(x,r y,r)
  _Ja[1][0] = _J[1][0] * _J[1][0];                        // x,s^2
  _Ja[1][1] = _J[1][1] * _J[1][1];                        // y,s^2
  _Ja[1][2] = 2.0 * (_J[1][0] * _J[1][1]);                // 2(x,s y,s)
  _Ja[2][0] = _J[0][0] * _J[1][0];                        // x,r x,s
  _Ja[2][1] = _J[0][1] * _J[1][1];                        // y,r y,s
  _Ja[2][2] = _J[0][0] * _J[1][1] + _J[0][1] * _J[1][0];  // (x,r y,s) + (x,s y,r)

  // Invert the matrix [Ja].

  double d;
  ShpMatInv(3, _Ja, _InvJa, &d);

  // Compute the second derivatives w.r.t. the parametric coords (ddNrs).

  DrvShpRST(pt, _dNrst, _ddNrs);

  // Compute x and y derivatives w.r.t. rr, ss, and rs.

  sNatDrv x, y;
  x.rr = x.rs = x.ss = y.rr = y.rs = y.ss = 0.0;
#if 1  
  for (int i = 0; i < NumNode; i++)
  {
    x.rr += _ddNrs[i].rr * coord[i].x;
    x.rs += _ddNrs[i].rs * coord[i].x;
    x.ss += _ddNrs[i].ss * coord[i].x;

    y.rr += _ddNrs[i].rr * coord[i].y;
    y.rs += _ddNrs[i].rs * coord[i].y;
    y.ss += _ddNrs[i].ss * coord[i].y;
  }
#endif  
  // cout << "x.rr: " << x.rr <<  "  " << "x.rs: " << x.rs << "  " << "x.ss: " << x.ss  << "\n";
  // cout << "y.rr: " << y.rr <<  "  " << "y.rs: " << y.rs << "  " << "y.ss: " << y.ss  << "\n\n";

  // Compute the second derivatives w.r.t. the cartesian coords (ddNxy).

  double v[3];
  for (int i = 0; i < NumNode; i++)
  {
    v[0] = _ddNrs[i].rr - (x.rr * dNxyz[i].x + y.rr * dNxyz[i].y);
    v[1] = _ddNrs[i].ss - (x.ss * dNxyz[i].x + y.ss * dNxyz[i].y);
    v[2] = _ddNrs[i].rs - (x.rs * dNxyz[i].x + y.rs * dNxyz[i].y);
    ddNxy[i].xx = _InvJa[0][0]*v[0] + _InvJa[0][1]*v[1] + _InvJa[0][2]*v[2];
    ddNxy[i].yy = _InvJa[1][0]*v[0] + _InvJa[1][1]*v[1] + _InvJa[1][2]*v[2];
    ddNxy[i].xy = _InvJa[2][0]*v[0] + _InvJa[2][1]*v[1] + _InvJa[2][2]*v[2];
  }
}

// ================================ GetFace ================================

int cShapePlane :: GetFace(int *corner, eShpType *type, int *nnode,
                            cNode **conn)
{
  *type   = Type;
  *nnode  = NumNode;

  for(int n = 0; n < NumNode; ++n) conn[n] = ElmNode[n];

  return 1;
}

// -------------------------------------------------------------------------
// Class cShapeT3:
// -------------------------------------------------------------------------

// =============================== cShapeT3 ================================

cShapeT3 :: cShapeT3(void)
{
  Type    = SHAPE_T3;
  TopType = TRIANGULAR_TOPOLOGY;
  NumNode = 3;
  ElmNode = new cNode*[NumNode];
}

// ============================== ~cShapeT3 ================================

cShapeT3 :: ~cShapeT3(void)
{
}

// ============================== NodeNatCoord =============================

void cShapeT3 :: NodeNatCoord(sNatCoord *c)
{
  c[0].r = 0.0; c[0].s = 0.0;
  c[1].r = 1.0; c[1].s = 0.0;
  c[2].r = 0.0; c[2].s = 1.0;
}

// ================================= ShpFunc ===============================

void cShapeT3 :: ShpFunc(sNatCoord p, double *N)
{
  double r = p.r;
  double s = p.s;

  N[0] = r;
  N[1] = s;
  N[2] = 1.0 - r - s;
}

// =============================== DrvShpRST ===============================

void cShapeT3 :: DrvShpRST(sNatCoord p, sNatCoord *dN)
{
  dN[0].r = 1.0;
  dN[0].s = 0.0;
  dN[1].r = 0.0;
  dN[1].s = 1.0;
  dN[2].r = -1.0;
  dN[2].s = -1.0;
}

// ================================ GetEdge ================================

int cShapeT3 :: GetEdge(int *corner, eShpType *type, int *nnode,
                        cNode **conn)
{
  static int NumEdge = 3;
  static tEdge Edge[] = { {SHAPE_L2_2D, 2, {0, 1}},
                          {SHAPE_L2_2D, 2, {1, 2}},
                          {SHAPE_L2_2D, 2, {2, 0}}};

  return(GetShapeEdge(NumEdge, Edge, corner, type, nnode, conn));
}


// -------------------------------------------------------------------------
// Class cShapeT6:
// -------------------------------------------------------------------------

// =============================== cShapeT6 ================================

cShapeT6 :: cShapeT6(void)
{
  Type    = SHAPE_T6;
  TopType = TRIANGULAR_TOPOLOGY;
  NumNode = 6;
  ElmNode = new cNode*[NumNode];
}

// ============================== ~cShapeT6 ================================

cShapeT6 :: ~cShapeT6(void)
{
}

// ============================== NodeNatCoord =============================

void cShapeT6 :: NodeNatCoord(sNatCoord *c)
{
  c[0].r = 0.0; c[0].s = 0.0;
  c[1].r = 0.5; c[1].s = 0.0;
  c[2].r = 1.0; c[2].s = 0.0;
  c[3].r = 0.5; c[3].s = 0.5;
  c[4].r = 0.0; c[4].s = 1.0;
  c[5].r = 0.0; c[5].s = 0.5;
}

// ================================= ShpFunc ===============================

void cShapeT6 :: ShpFunc(sNatCoord p, double *N)
{
  double r = p.r;
  double s = p.s;

  N[0] = 1.0 - 3.0*r - 3.0*s + 4.0*r*s + 2.0*r*r + 2.0*s*s;
  N[1] = 4.0*r - 4.0*r*r - 4.0*r*s;
  N[2] = 2.0*r*r - r;
  N[3] = 4.0*r*s;
  N[4] = 2.0*s*s - s;
  N[5] = 4.0*s - 4.0*r*s - 4.0*s*s;
}

// =============================== DrvShpRST ===============================

void cShapeT6 :: DrvShpRST(sNatCoord p, sNatCoord *dN)
{
  double r = p.r;
  double s = p.s;

  dN[0].r =  4.0*r + 4.0*s - 3.0;
  dN[1].r =  4.0 - 8.0*r - 4.0*s;
  dN[2].r =  4.0*r - 1.0;
  dN[3].r =  4.0*s;
  dN[4].r =  0.0;
  dN[5].r = -4.0*s;

  dN[0].s =  4.0*r + 4.0*s - 3.0;
  dN[1].s = -4.0*r;
  dN[2].s =  0.0;
  dN[3].s =  4.0*r;
  dN[4].s =  4.0*s - 1.0;
  dN[5].s =  4.0 - 4.0*r - 8.0*s;
}

// ================================ GetEdge ================================

int cShapeT6 :: GetEdge(int *corner, eShpType *type, int *nnode,
                        cNode **conn)
{
  static int NumEdge = 3;
  static tEdge Edge[] = { {SHAPE_L3_2D, 3, {0, 1, 2}},
                          {SHAPE_L3_2D, 3, {2, 3, 4}},
                          {SHAPE_L3_2D, 3, {4, 5, 0}}};

  return(GetShapeEdge(NumEdge, Edge, corner, type, nnode, conn));
}

// -------------------------------------------------------------------------
// Class cShapeQ4:
// -------------------------------------------------------------------------

// =============================== cShapeQ4 ================================

cShapeQ4 :: cShapeQ4(void)
{
  Type    = SHAPE_Q4;
  TopType = QUADRILATERAL_TOPOLOGY;
  NumNode = 4;
  ElmNode = new cNode*[NumNode];
}

// ============================== ~cShapeQ4 ================================

cShapeQ4 :: ~cShapeQ4(void)
{
}

// ============================== NodeNatCoord =============================

void cShapeQ4 :: NodeNatCoord(sNatCoord *c)
{
  c[0].r = -1.0; c[0].s = -1.0;
  c[1].r =  1.0; c[1].s = -1.0;
  c[2].r =  1.0; c[2].s =  1.0;
  c[3].r = -1.0; c[3].s =  1.0;
}

// ================================= ShpFunc ===============================

void cShapeQ4 :: ShpFunc(sNatCoord p, double *N)
{
  double r = p.r;
  double s = p.s;

  N[0] = 0.25*(1.0 - r)*(1.0 - s);
  N[1] = 0.25*(1.0 + r)*(1.0 - s);
  N[2] = 0.25*(1.0 + r)*(1.0 + s);
  N[3] = 0.25*(1.0 - r)*(1.0 + s);
}

// =============================== DrvShpRST ===============================

void cShapeQ4 :: DrvShpRST(sNatCoord p, sNatCoord *dN)
{
  double r = p.r;
  double s = p.s;

  dN[0].r = -0.25*(1.0 - s);
  dN[1].r =  0.25*(1.0 - s);
  dN[2].r =  0.25*(1.0 + s);
  dN[3].r = -0.25*(1.0 + s);

  dN[0].s = -0.25*(1.0 - r);
  dN[1].s = -0.25*(1.0 + r);
  dN[2].s =  0.25*(1.0 + r);
  dN[3].s =  0.25*(1.0 - r);
}

// ================================ GetEdge ================================

int cShapeQ4 :: GetEdge(int *corner, eShpType *type, int *nnode,
                        cNode **conn)
{
  static int NumEdge = 4;
  static tEdge Edge[] = { {SHAPE_L2_2D, 2, {0, 1}},
                          {SHAPE_L2_2D, 2, {1, 2}},
                          {SHAPE_L2_2D, 2, {2, 3}},
                          {SHAPE_L2_2D, 2, {3, 0}}};

  return(GetShapeEdge(NumEdge, Edge, corner, type, nnode, conn));
}

// -------------------------------------------------------------------------
// Class cShapeQ8:
// -------------------------------------------------------------------------

// =============================== cShapeQ8 ================================

cShapeQ8 :: cShapeQ8(void)
{
  Type    = SHAPE_Q8;
  TopType = QUADRILATERAL_TOPOLOGY;
  NumNode = 8;
  ElmNode = new cNode*[NumNode];
}

// ============================== ~cShapeQ8 ================================

cShapeQ8 :: ~cShapeQ8(void)
{
}

// ============================== NodeNatCoord =============================

void cShapeQ8 :: NodeNatCoord(sNatCoord *c)
{
  c[0].r = -1.0; c[0].s = -1.0;
  c[1].r =  0.0; c[1].s = -1.0;
  c[2].r =  1.0; c[2].s = -1.0;
  c[3].r =  1.0; c[3].s =  0.0;
  c[4].r =  1.0; c[4].s =  1.0;
  c[5].r =  0.0; c[5].s =  1.0;
  c[6].r = -1.0; c[6].s =  1.0;
  c[7].r = -1.0; c[7].s =  0.0;
}

// ================================= ShpFunc ===============================

void cShapeQ8 :: ShpFunc(sNatCoord p, double *N)
{
  double r = p.r;
  double s = p.s;

  N[0] = 0.25*(1.0 - r)*(1.0 - s) - 0.25*(1.0 - r)*(1.0 - s*s) -
         0.25*(1.0 - r*r)*(1.0 - s);
  N[1] = 0.50*(1.0 - r*r)*(1.0 - s);
  N[2] = 0.25*(1.0 + r)*(1.0 - s) - 0.25*(1.0 - r*r)*(1.0 - s) -
         0.25*(1.0 - s*s)*(1.0 + r);
  N[3] = 0.50*(1.0 - s*s)*(1.0 + r);
  N[4] = 0.25*(1.0 + r)*(1.0 + s) - 0.25*(1.0 - s*s)*(1.0 + r) -
         0.25*(1.0 + s)*(1.0 - r*r);
  N[5] = 0.50*(1.0 - r*r)*(1.0 + s);
  N[6] = 0.25*(1.0 - r)*(1.0 + s) - 0.25*(1.0 - r*r)*(1.0 + s) -
         0.25*(1.0 - s*s)*(1.0 - r);
  N[7] = 0.50*(1.0 - s*s)*(1.0 - r);
}

// =============================== DrvShpRST ===============================

void cShapeQ8 :: DrvShpRST(sNatCoord p, sNatCoord *dN)
{
  double r = p.r;
  double s = p.s;

  dN[0].r = (2.0*r - 2.0*r*s - s*s + s)/4.0;
  dN[1].r = r*s - r;
  dN[2].r = (2.0*r - 2.0*r*s + s*s - s)/4.0;
  dN[3].r = (1.0 - s*s)/2.0;
  dN[4].r = (2.0*r + 2.0*r*s + s*s + s)/4.0;
  dN[5].r = -r - r*s;
  dN[6].r = (2.0*r + 2.0*r*s - s*s - s)/4.0;
  dN[7].r = (s*s - 1.0)/2.0;

  dN[0].s = (2.0*s - r*r - 2.0*r*s + r)/4.0;
  dN[1].s = (r*r - 1.0)/2.0;
  dN[2].s = (2.0*s - r*r + 2.0*r*s - r)/ 4.0;
  dN[3].s = -s - r*s;
  dN[4].s = (2.0*s + r*r + 2.0*r*s + r)/4.0;
  dN[5].s = (1.0 - r*r)/2.0;
  dN[6].s = (2.0*s + r*r - 2.0*r*s - r)/4.0;
  dN[7].s = r*s - s;
}

// ================================ GetEdge ================================

int cShapeQ8 :: GetEdge(int *corner, eShpType *type, int *nnode,
                        cNode **conn)
{
  static int NumEdge = 4;
  static tEdge Edge[] = { {SHAPE_L3_2D, 3, {0, 1, 2}},
                          {SHAPE_L3_2D, 3, {2, 3, 4}},
                          {SHAPE_L3_2D, 3, {4, 5, 6}},
                          {SHAPE_L3_2D, 3, {6, 7, 0}}};

  return(GetShapeEdge(NumEdge, Edge, corner, type, nnode, conn));
}

// -------------------------------------------------------------------------
// Class cShapeQ9:
// -------------------------------------------------------------------------

// =============================== cShapeQ9 ================================

cShapeQ9 :: cShapeQ9(void)
{
  Type    = SHAPE_Q9;
  TopType = QUADRILATERAL_TOPOLOGY;
  NumNode = 9;
  ElmNode = new cNode*[NumNode];
}

// ============================== ~cShapeQ9 ================================

cShapeQ9 :: ~cShapeQ9(void)
{
}

// ============================== NodeNatCoord =============================

void cShapeQ9 :: NodeNatCoord(sNatCoord *c)
{
  c[0].r = -1.0; c[0].s = -1.0;
  c[1].r =  0.0; c[1].s = -1.0;
  c[2].r =  1.0; c[2].s = -1.0;
  c[3].r =  1.0; c[3].s =  0.0;
  c[4].r =  1.0; c[4].s =  1.0;
  c[5].r =  0.0; c[5].s =  1.0;
  c[6].r = -1.0; c[6].s =  1.0;
  c[7].r = -1.0; c[7].s =  0.0;
  c[8].r =  0.0; c[8].s =  0.0;
}

// ================================= ShpFunc ===============================

void cShapeQ9 :: ShpFunc(sNatCoord p, double *N)
{
  double r = p.r;
  double s = p.s;

  N[0] = 0.25*r*(r - 1.0)*s*(s - 1.0);
  N[1] = 0.50*(1.0 - r*r)*s*(s - 1.0);
  N[2] = 0.25*r*(r + 1.0)*s*(s - 1.0);
  N[3] = 0.50*(1.0 - s*s)*r*(r + 1.0);
  N[4] = 0.25*r*(r + 1.0)*s*(s + 1.0);
  N[5] = 0.50*(1.0 - r*r)*s*(s + 1.0);
  N[6] = 0.25*r*(r - 1.0)*s*(s + 1.0);
  N[7] = 0.50*(1.0 - s*s)*r*(r - 1.0);
  N[8] = (1.0 - s*s)*(1.0 - r*r);
}

// =============================== DrvShpRST ===============================

void cShapeQ9 :: DrvShpRST(sNatCoord p, sNatCoord *dN)
{
  double r = p.r;
  double s = p.s;

  dN[0].r =  0.25*s*(s - 1.0)*(2.0*r - 1.0);
  dN[1].r = -r*s*(s - 1.0);
  dN[2].r =  0.25*s*(s - 1.0)*(2.0*r + 1.0);
  dN[3].r =  0.50*(1.0 - s*s)*(2.0*r + 1.0);
  dN[4].r =  0.25*s*(s + 1.0)*(2.0*r + 1.0);
  dN[5].r = -r*s*(s + 1.0);
  dN[6].r =  0.25*s*(s + 1.0)*(2.0*r - 1.0);
  dN[7].r =  0.50*(1.0 - s*s)*(2.0*r - 1.0);
  dN[8].r =  -2.0*r*(1.0 - s*s);

  dN[0].s =  0.25*r*(r - 1.0)*(2.0*s - 1.0);
  dN[1].s =  0.50*(1.0 - r*r)*(2.0*s - 1.0);
  dN[2].s =  0.25*r*(r + 1.0)*(2.0*s - 1.0);
  dN[3].s = -s*r*(r + 1.0);
  dN[4].s =  0.25*r*(r + 1.0)*(2.0*s + 1.0);
  dN[5].s =  0.50*(1.0 - r*r)*(2.0*s + 1.0);
  dN[6].s =  0.25*r*(r - 1.0)*(2.0*s + 1.0);
  dN[7].s = -s*r*(r - 1.0);
  dN[8].s =  -2.0*s*(1.0 - r*r);
}

// ================================ GetEdge ================================

int cShapeQ9 :: GetEdge(int *corner, eShpType *type, int *nnode,
                        cNode **conn)
{
  static int NumEdge = 4;
  static tEdge Edge[] = { {SHAPE_L3_2D, 3, {0, 1, 2}},
                          {SHAPE_L3_2D, 3, {2, 3, 4}},
                          {SHAPE_L3_2D, 3, {4, 5, 6}},
                          {SHAPE_L3_2D, 3, {6, 7, 0}}};

  return(GetShapeEdge(NumEdge, Edge, corner, type, nnode, conn));
}

// ======================================================= End of file =====
