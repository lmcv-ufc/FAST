// -------------------------------------------------------------------------
// shpinf.cpp  - this file contains methods for the 2D mapped
//               infinite elements shape.
// -------------------------------------------------------------------------
// Created:      23-Apr-2006     Aurea Silva de Holanda
//
// Modified:     17-May-2011     Evandro Parente Junior
//               Implementation of non-isoparametric elements.
// -------------------------------------------------------------------------

#include <stdio.h>

#include "shpinf.h"

// -------------------------------------------------------------------------
// Public methods:
//

// =========================== cShapePlaneInf ==============================

cShapePlaneInf :: cShapePlaneInf(void)
{
}

// ============================ ~cShapePlaneInf ============================

cShapePlaneInf :: ~cShapePlaneInf(void)
{
}


// -------------------------------------------------------------------------
// Class cShapeL6Inf:
// -------------------------------------------------------------------------

// ============================== cShapeL6Inf ==============================

cShapeL6Inf :: cShapeL6Inf(void)
{
  Type = SHAPE_L6_INF;
  NumNode = 6;
  ElmNode = new cNode*[NumNode];
}

// ============================== ~cShapeL6Inf =============================

cShapeL6Inf :: ~cShapeL6Inf(void)
{
}

// ============================== NodeNatCoord =============================

void cShapeL6Inf :: NodeNatCoord(sNatCoord *c)
{
  c[0].r = -1.0; c[0].s = -1.0;
  c[1].r =  0.0; c[1].s = -1.0;
  c[2].r =  1.0; c[2].s = -1.0;
  c[3].r =  1.0; c[3].s =  0.0;
  c[4].r =  0.0; c[4].s =  0.0;
  c[5].r = -1.0; c[5].s =  0.0;
}

// ================================= MapFunc ===============================

void cShapeL6Inf :: MapFunc(sNatCoord p, double *M)
{
  double r = p.r;
  double s = p.s;

  M[0] = (1.0 - r)*r*s/(1.0 - s);
  M[1] = 2.0*(r*r - 1.0)*s/(1.0 - s);
  M[2] = - (1.0 + r)*r*s/(1.0 - s);
  M[3] = 0.5*r*(r + 1.0)*(1.0 + s)/(1.0 - s);
  M[4] = (1.0 - r*r)*(1.0 + s)/(1.0 - s);
  M[5] = 0.5*r*(r - 1.0)*(1.0 + s)/(1.0 - s);
}

// =============================== DrvMapRST ===============================

void cShapeL6Inf :: DrvMapRST(sNatCoord p, sNatCoord *dM)
{
  double r = p.r;
  double s = p.s;

  dM[0].r = ((1.0 - 2.0*r)*s)/(1.0 - s);
  dM[0].s =((1.0 - r)*r)/((1.0 - s)*(1.0 - s));
  dM[1].r =(4.0*r*s)/(1 - s);
  dM[1].s =(2.0*(r*r - 1.0))/((1.0 - s)*(1.0 - s));
  dM[2].r =((-s)*(1.0 + 2.0*r))/(1.0 - s);
  dM[2].s =((-r)*(1.0 + r))/((1.0 - s)*(1.0 - s));
  dM[3].r =((1.0 + s)*(1 + 2.0*r))/(2.0*(1.0 - s));
  dM[3].s =(r*(r + 1.0))/((1.0 - s)*(1.0 - s));
  dM[4].r =((-2.0)*r*(1.0 + s))/(1.0 - s);
  dM[4].s =(2.0*(1.0 - r*r))/((1.0 - s)*(1.0 - s));
  dM[5].r =((1.0 + s)*(2.0*r - 1))/(2.0*(1.0 - s));
  dM[5].s =(r*(r - 1.0))/((1.0 - s)*(1.0 - s));
}

// ================================= ShpFunc ===============================

void cShapeL6Inf :: ShpFunc(sNatCoord p, double *N)
{
  double r = p.r;
  double s = p.s;

  N[0] = 0.25*r*(r - 1.0)*s*(s - 1.0);
  N[1] = 0.50*(1.0 - r*r)*s*(s - 1.0);
  N[2] = 0.25*r*(r + 1.0)*s*(s - 1.0);
  N[3] = 0.50*(1.0 - s*s)*r*(r + 1.0);
  N[4] = (1.0 - s*s)*(1.0 - r*r);
  N[5] = 0.50*(1.0 - s*s)*r*(r - 1.0);
}

// =============================== DrvShpRST ===============================

void cShapeL6Inf :: DrvShpRST(sNatCoord p, sNatCoord *dN)
{
  double r = p.r;
  double s = p.s;

  dN[0].r = 0.25*s*(s - 1.0)*(2.0*r - 1.0);
  dN[0].s = 0.25*r*(r - 1.0)*(2.0*s - 1.0);
  dN[1].r = (-r)*s*(s - 1.0);
  dN[1].s = 0.50*(1.0 - r*r)*(2.0*s - 1.0);
  dN[2].r = 0.25*s*(s - 1.0)*(2.0*r + 1.0);
  dN[2].s = 0.25*r*(r + 1.0)*(2.0*s - 1.0);
  dN[3].r = 0.50*(1.0 - s*s)*(2.0*r + 1.0);
  dN[3].s = (-s)*r*(r + 1.0);
  dN[4].r = (-2.0)*r*(1.0 - s*s);
  dN[4].s = (-2.0)*s*(1.0 - r*r);
  dN[5].r = 0.50*(1.0 - s*s)*(2.0*r - 1.0);
  dN[5].s = (-s)*r*(r - 1.0);
}

// ================================ GetEdge ================================

int cShapeL6Inf :: GetEdge(int *corner, eShpType *type, int *nnode,
                           cNode **conn)
{
  static int NumEdge = 3;
  static tEdge Edge[] = { {SHAPE_L2_2D, 2, {0, 1}},
                          {SHAPE_L2_2D, 2, {1, 2}},
                          {SHAPE_L2_2D, 2, {2, 1}}};

  return(GetShapeEdge(NumEdge, Edge, corner, type, nnode, conn));
}

/*
// -------------------------------------------------------------------------
// Class cShapeL4Inf:
// -------------------------------------------------------------------------

// ============================== cShapeL4Inf ==============================

cShapeL4Inf :: cShapeL4Inf (void)
{
  Type = SHAPE_L4_INF;
  NumNode = 4;
  ElmNode = new cNode*[NumNode];
}

// ============================== ~cShapeL4Inf =============================

cShapeL4Inf :: ~cShapeL4Inf (void)
{
}

// ============================== NodeNatCoord =============================

void cShapeL4Inf :: NodeNatCoord(sNatCoord *c)
{
  c[0].r = 0.0; c[0].s = 0.0;
  c[1].r = 0.5; c[1].s = 0.0;
  c[2].r = 1.0; c[2].s = 0.0;
  c[3].r = 0.5; c[3].s = 0.5;
}

// ================================= ShpFunc ===============================

void cShapeL4Inf :: ShpFunc(sNatCoord p, double *N)
{
  double r = p.r;
  double s = p.s;

  N[0] = 1.0 - 3.0*r - 3.0*s + 4.0*r*s + 2.0*r*r + 2.0*s*s;
  N[1] = 4.0*r - 4.0*r*r - 4.0*r*s;
  N[2] = 2.0*r*r - r;
  N[3] = 4.0*r*s;
}

// =============================== DrvShpRST ===============================

void cShapeL4Inf :: DrvShpRST(sNatCoord p, sNatCoord *dN)
{
  double r = p.r;
  double s = p.s;

  dN[0].r =  4.0*r + 4.0*s - 3.0;
  dN[1].r =  4.0 - 8.0*r - 4.0*s;
  dN[2].r =  4.0*r - 1.0;
  dN[3].r =  4.0*s;

  dN[0].s =  4.0*r + 4.0*s - 3.0;
  dN[1].s = -4.0*r;
  dN[2].s =  0.0;
  dN[3].s =  4.0*r;
}

// ================================ GetEdge ================================

int cShapeL4Inf :: GetEdge(int *corner, eShpType *type, int *nnode,
                        cNode **conn)
{
  static int NumEdge = 3;
  static tEdge Edge[] = { {SHAPE_L3_2D, 3, {0, 1, 2}},
                          {SHAPE_L3_2D, 3, {2, 3, 4}},
                          {SHAPE_L3_2D, 3, {4, 5, 0}}};

  return(GetShapeEdge(NumEdge, Edge, corner, type, nnode, conn));
}
*/
// ============================================================= End of File
