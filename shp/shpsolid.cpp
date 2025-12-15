// -------------------------------------------------------------------------
// shpsolid.cpp - implementation of the Shape class.
// -------------------------------------------------------------------------
// Created:      06-Sep-2005     Aurea Silva de Holanda
//
// Modified:     17-May-2011     Evandro Parente Junior
//               Implementation of non-isoparametric elements.
//
// Modified:     03-Oct-2013     Evandro Parente Junior
//               Implementation of tVector.
//
// Modified:     22-May-2023     Evandro Parente Junior
//               Implementation of Tet4 and Tet10 shapes.
// -------------------------------------------------------------------------

#include <vector>

using namespace std;

#include "shpsolid.h"
#include "vec.h"

#ifdef _OMP_
#include "omp.h"
#endif

// -------------------------------------------------------------------------
// Static variables:
//
static double       _J[3][3];    // Jabobian matrix
static double    _InvJ[3][3];    // Inverse of Jacobian matrix
static sNatCoord _dNrst[500];    // Shape func derivatives/parametric coords
static sNatCoord _dMrst[500];    // Mapping func derivatives/parametric coords

//static vector<sNatCoord> _dNrst;
//static vector<sNatCoord> _dMrst;

#pragma omp threadprivate(_J,_InvJ,_dNrst,_dMrst)

// -------------------------------------------------------------------------
// Public methods:
//

// ============================== cShapeSolid ==============================

cShapeSolid :: cShapeSolid(void)
{
}

// ============================= ~cShapeSolid ==============================

cShapeSolid :: ~cShapeSolid(void)
{
}

// =============================== DrvShpXYZ ===============================

void cShapeSolid ::  DrvShpXYZ(sNatCoord pt, sNodeCoord *coord,
                               double *detJ, sNodeCoord *dNxyz)
{
  int i;

  // Compute dNrs.

  DrvShpRST(pt, _dNrst);

  // Compute [J].

  ShpMatZero(3, _J);
  if (Isoparametric( )) // Isoparametric mapping
  {
    for (i = 0; i < NumNode; i++)
    {
      _J[0][0] += _dNrst[i].r*coord[i].x;
      _J[0][1] += _dNrst[i].r*coord[i].y;
      _J[0][2] += _dNrst[i].r*coord[i].z;

      _J[1][0] += _dNrst[i].s*coord[i].x;
      _J[1][1] += _dNrst[i].s*coord[i].y;
      _J[1][2] += _dNrst[i].s*coord[i].z;

      _J[2][0] += _dNrst[i].t*coord[i].x;
      _J[2][1] += _dNrst[i].t*coord[i].y;
      _J[2][2] += _dNrst[i].t*coord[i].z;
    }
  }
  else                  // Non-isoparametric mapping
  {
    DrvMapRST(pt, _dMrst);
    for (i = 0; i < NumNode; i++)
    {
      _J[0][0] += _dMrst[i].r*coord[i].x;
      _J[0][1] += _dMrst[i].r*coord[i].y;
      _J[0][2] += _dMrst[i].r*coord[i].z;

      _J[1][0] += _dMrst[i].s*coord[i].x;
      _J[1][1] += _dMrst[i].s*coord[i].y;
      _J[1][2] += _dMrst[i].s*coord[i].z;

      _J[2][0] += _dMrst[i].t*coord[i].x;
      _J[2][1] += _dMrst[i].t*coord[i].y;
      _J[2][2] += _dMrst[i].t*coord[i].z;
    }
  }

  // Compute |J| and inv(J).

  ShpMatInv(3, _J, _InvJ, detJ);

  // Compute Ni,xyz.

  for (i = 0; i < NumNode; i++)
  {
    dNxyz[i].x = _InvJ[0][0] * _dNrst[i].r +
                 _InvJ[0][1] * _dNrst[i].s +
                 _InvJ[0][2] * _dNrst[i].t;

    dNxyz[i].y = _InvJ[1][0] * _dNrst[i].r +
                 _InvJ[1][1] * _dNrst[i].s +
                 _InvJ[1][2] * _dNrst[i].t;

    dNxyz[i].z = _InvJ[2][0] * _dNrst[i].r +
                 _InvJ[2][1] * _dNrst[i].s +
                 _InvJ[2][2] * _dNrst[i].t;
  }
}

// ============================== tVector ==================================

void cShapeSolid :: tVector(sNatCoord pt, sNodeCoord *coord, cVector &w)
{
  // Compute dMrs.

  DrvMapRST(pt, _dMrst);

  // Compute {w} = d{x}/dt.

  w.Zero( );
  for (int i = 0; i < NumNode; i++)
  {
    w[0] += _dMrst[i].t*coord[i].x;
    w[1] += _dMrst[i].t*coord[i].y;
    w[2] += _dMrst[i].t*coord[i].z;
  }

  // Return the unit vector {w}/|w|.

  w.Normalize( );
}

// -------------------------------------------------------------------------
// Class cShapeBrick8:
// -------------------------------------------------------------------------

// ============================== cShapeBrick8 =============================

cShapeBrick8 :: cShapeBrick8(void)
{
  Type    = SHAPE_BRICK8;
  TopType = HEXAHEDRAL_TOPOLOGY;
  NumNode = 8;
  ElmNode = new cNode*[NumNode];
}

// ============================= ~cShapeBrick8 =============================

cShapeBrick8 :: ~cShapeBrick8(void)
{
}

// ============================== NodeNatCoord =============================

void cShapeBrick8 :: NodeNatCoord(sNatCoord *c)
{
   c[0].r = -1.0; c[0].s = -1.0; c[0].t = 1.0;
   c[1].r =  1.0; c[1].s = -1.0; c[1].t = 1.0;
   c[2].r =  1.0; c[2].s =  1.0; c[2].t = 1.0;
   c[3].r = -1.0; c[3].s =  1.0; c[3].t = 1.0;

   c[4].r = -1.0; c[4].s = -1.0; c[4].t = -1.0;
   c[5].r =  1.0; c[5].s = -1.0; c[5].t = -1.0;
   c[6].r =  1.0; c[6].s =  1.0; c[6].t = -1.0;
   c[7].r = -1.0; c[7].s =  1.0; c[7].t = -1.0;
}

// ================================= ShpFunc ===============================

void cShapeBrick8 :: ShpFunc(sNatCoord p, double *N)
{
  double r = p.r;
  double s = p.s;
  double t = p.t;

  N[0] = 0.125*(1.0-r)*(1.0-s)*(1.0+t);
  N[1] = 0.125*(1.0+r)*(1.0-s)*(1.0+t);
  N[2] = 0.125*(1.0+r)*(1.0+s)*(1.0+t);
  N[3] = 0.125*(1.0-r)*(1.0+s)*(1.0+t);
  N[4] = 0.125*(1.0-r)*(1.0-s)*(1.0-t);
  N[5] = 0.125*(1.0+r)*(1.0-s)*(1.0-t);
  N[6] = 0.125*(1.0+r)*(1.0+s)*(1.0-t);
  N[7] = 0.125*(1.0-r)*(1.0+s)*(1.0-t);
}

// =============================== DrvShpRST ===============================

void cShapeBrick8 :: DrvShpRST(sNatCoord p, sNatCoord *dN)
{
  double r = p.r;
  double s = p.s;
  double t = p.t;

   dN[0].r = -0.125*(1.0-s)*(1.0+t);
   dN[1].r =  0.125*(1.0-s)*(1.0+t);
   dN[2].r =  0.125*(1.0+s)*(1.0+t);
   dN[3].r = -0.125*(1.0+s)*(1.0+t);
   dN[4].r = -0.125*(1.0-s)*(1.0-t);
   dN[5].r =  0.125*(1.0-s)*(1.0-t);
   dN[6].r =  0.125*(1.0+s)*(1.0-t);
   dN[7].r = -0.125*(1.0+s)*(1.0-t);

   dN[0].s = -0.125*(1.0-r)*(1.0+t);
   dN[1].s = -0.125*(1.0+r)*(1.0+t);
   dN[2].s =  0.125*(1.0+r)*(1.0+t);
   dN[3].s =  0.125*(1.0-r)*(1.0+t);
   dN[4].s = -0.125*(1.0-r)*(1.0-t);
   dN[5].s = -0.125*(1.0+r)*(1.0-t);
   dN[6].s =  0.125*(1.0+r)*(1.0-t);
   dN[7].s =  0.125*(1.0-r)*(1.0-t);

   dN[0].t =  0.125*(1.0-r)*(1.0-s);
   dN[1].t =  0.125*(1.0+r)*(1.0-s);
   dN[2].t =  0.125*(1.0+r)*(1.0+s);
   dN[3].t =  0.125*(1.0-r)*(1.0+s);
   dN[4].t = -0.125*(1.0-r)*(1.0-s);
   dN[5].t = -0.125*(1.0+r)*(1.0-s);
   dN[6].t = -0.125*(1.0+r)*(1.0+s);
   dN[7].t = -0.125*(1.0-r)*(1.0+s);
}

// ================================ GetEdge ================================

int cShapeBrick8 :: GetEdge(int *corner, eShpType *type, int *nnode,
                            cNode **conn)
{
  static int NumEdge = 12;
  static tEdge Edge[] = { {SHAPE_L2_3D, 2, {0, 1}},
                          {SHAPE_L2_3D, 2, {1, 2}},
                          {SHAPE_L2_3D, 2, {2, 3}},
                          {SHAPE_L2_3D, 2, {3, 0}},
                          {SHAPE_L2_3D, 2, {0, 4}},
                          {SHAPE_L2_3D, 2, {1, 5}},
                          {SHAPE_L2_3D, 2, {2, 6}},
                          {SHAPE_L2_3D, 2, {3, 7}},
                          {SHAPE_L2_3D, 2, {4, 5}},
                          {SHAPE_L2_3D, 2, {5, 6}},
                          {SHAPE_L2_3D, 2, {6, 7}},
                          {SHAPE_L2_3D, 2, {7, 4}}};

  return(GetShapeEdge(NumEdge, Edge, corner, type, nnode, conn));
}

// ================================ GetFace ================================

int cShapeBrick8 :: GetFace(int *corner, eShpType *type, int *nnode,
                            cNode **conn)
{
  static int NumFace = 6;
  static tFace Face[] = {{SHAPE_Q4_SURF, 4, {0, 1, 2, 3}, 4, {0, 1, 2, 3}},
                         {SHAPE_Q4_SURF, 4, {3, 2, 6, 7}, 4, {3, 2, 6, 7}},
                         {SHAPE_Q4_SURF, 4, {1, 5, 6, 2}, 4, {1, 5, 6, 2}},
                         {SHAPE_Q4_SURF, 4, {0, 4, 5, 1}, 4, {0, 4, 5, 1}},
                         {SHAPE_Q4_SURF, 4, {0, 3, 7, 4}, 4, {0, 3, 7, 4}},
                         {SHAPE_Q4_SURF, 4, {4, 7, 6, 5}, 4, {4, 7, 6, 5}}};

  return(GetShapeFace(NumFace, Face, corner, type, nnode, conn));
}


// -------------------------------------------------------------------------
// Class cShapeBrick20:
// -------------------------------------------------------------------------

// ============================= cShapeBrick20 =============================

cShapeBrick20 :: cShapeBrick20(void)
{
  Type    = SHAPE_BRICK20;
  TopType = HEXAHEDRAL_TOPOLOGY;
  NumNode = 20;
  ElmNode = new cNode*[NumNode];
}

// ============================ ~cShapeBrick20 =============================

cShapeBrick20 :: ~cShapeBrick20(void)
{
}

// ============================== NodeNatCoord =============================

void cShapeBrick20 :: NodeNatCoord(sNatCoord *c)
{
   c[0].r = -1.0; c[0].s = -1.0; c[0].t =  1.0;
   c[1].r =  0.0; c[1].s = -1.0; c[1].t =  1.0;
   c[2].r =  1.0; c[2].s = -1.0; c[2].t =  1.0;
   c[3].r =  1.0; c[3].s =  0.0; c[3].t =  1.0;
   c[4].r =  1.0; c[4].s =  1.0; c[4].t =  1.0;
   c[5].r =  0.0; c[5].s =  1.0; c[5].t =  1.0;
   c[6].r = -1.0; c[6].s =  1.0; c[6].t =  1.0;
   c[7].r = -1.0; c[7].s =  0.0; c[7].t =  1.0;

   c[ 8].r = -1.0; c[ 8].s = -1.0; c[ 8].t =  0.0;
   c[ 9].r =  1.0; c[ 9].s = -1.0; c[ 9].t =  0.0;
   c[10].r =  1.0; c[10].s =  1.0; c[10].t =  0.0;
   c[11].r = -1.0; c[11].s =  1.0; c[11].t =  0.0;

   c[12].r = -1.0; c[12].s = -1.0; c[12].t = -1.0;
   c[13].r =  0.0; c[13].s = -1.0; c[13].t = -1.0;
   c[14].r =  1.0; c[14].s = -1.0; c[14].t = -1.0;
   c[15].r =  1.0; c[15].s =  0.0; c[15].t = -1.0;
   c[16].r =  1.0; c[16].s =  1.0; c[16].t = -1.0;
   c[17].r =  0.0; c[17].s =  1.0; c[17].t = -1.0;
   c[18].r = -1.0; c[18].s =  1.0; c[18].t = -1.0;
   c[19].r = -1.0; c[19].s =  0.0; c[19].t = -1.0;
}

// ================================= ShpFunc ===============================

void cShapeBrick20 :: ShpFunc(sNatCoord p, double *N)
{
  double r = p.r;
  double s = p.s;
  double t = p.t;
  double r2 = r*r;
  double s2 = s*s;
  double t2 = t*t;

  N[ 0] = 0.125 * (1.0-r) * (1.0-s) * (1.0+t)*(-r-s+t-2.0);
  N[ 1] = 0.250 * (1.0-r2)* (1.0-s) * (1.0+t);
  N[ 2] = 0.125 * (1.0+r) * (1.0-s) * (1.0+t)*(r-s+t-2.0);
  N[ 3] = 0.250 * (1.0+r) * (1.0-s2)* (1.0+t);
  N[ 4] = 0.125 * (1.0+r) * (1.0+s) * (1.0+t)*(r+s+t-2.0);
  N[ 5] = 0.250 * (1.0-r2)* (1.0+s) * (1.0+t);
  N[ 6] = 0.125 * (1.0-r) * (1.0+s) * (1.0+t)*(-r+s+t-2.0);
  N[ 7] = 0.250 * (1.0-r) * (1.0-s2)* (1.0+t);
  N[ 8] = 0.250 * (1.0-r) * (1.0-s) * (1.0-t2);
  N[ 9] = 0.250 * (1.0+r) * (1.0-s) * (1.0-t2);
  N[10] = 0.250 * (1.0+r) * (1.0+s) * (1.0-t2);
  N[11] = 0.250 * (1.0-r) * (1.0+s) * (1.0-t2);
  N[12] = 0.125 * (1.0-r) * (1.0-s) * (1.0-t)*(-r-s-t-2.0);
  N[13] = 0.250 * (1.0-r2)* (1.0-s) * (1.0-t);
  N[14] = 0.125 * (1.0+r) * (1.0-s) * (1.0-t)*(r-s-t-2.0);
  N[15] = 0.250 * (1.0+r) * (1.0-s2)* (1.0-t);
  N[16] = 0.125 * (1.0+r) * (1.0+s) * (1.0-t)*(r+s-t-2.0);
  N[17] = 0.250 * (1.0-r2)* (1.0+s) * (1.0-t);
  N[18] = 0.125 * (1.0-r) * (1.0+s) * (1.0-t)*(-r+s-t-2.0);
  N[19] = 0.250 * (1.0-r) * (1.0-s2)* (1.0-t);
}

// =============================== DrvShpRST ===============================

void cShapeBrick20 :: DrvShpRST(sNatCoord p, sNatCoord *dN)
{
  double r = p.r;
  double s = p.s;
  double t = p.t;
  double r2 = r*r;
  double s2 = s*s;
  double t2 = t*t;

  dN[ 0].r = -0.125 * (1.0-s) * (1.0+t) * (-2.0*r-s+t-1.0);
  dN[ 1].r = -0.500 * r * (1.0-s) * (1.0+t);
  dN[ 2].r =  0.125 * (1.0-s) * (1.0+t) * (2.0*r-s+t-1.0);
  dN[ 3].r =  0.250 * (1.0-s2) * (1.0+t);
  dN[ 4].r =  0.125 * (1.0+s) * (1.0+t) * (2.0*r+s+t-1.0);
  dN[ 5].r = -0.500 * r * (1.0+s) * (1.0+t);
  dN[ 6].r = -0.125 * (1.0+s) * (1.0+t) * (-2.0*r+s+t-1.0);
  dN[ 7].r = -0.250 * (1.0-s2) * (1.0+t);
  dN[ 8].r = -0.250 * (1.0-s) * (1.0-t2);
  dN[ 9].r =  0.250 * (1.0-s) * (1.0-t2);
  dN[10].r =  0.250 * (1.0+s) * (1.0-t2);
  dN[11].r = -0.250 * (1.0+s) * (1.0-t2);
  dN[12].r = -0.125 * (1.0-s) * (1.0-t) * (-2.0*r-s-t-1.0);
  dN[13].r = -0.500 * r * (1.0-s) * (1.0-t);
  dN[14].r =  0.125 * (1.0-s) * (1.0-t) * (2.0*r-s-t-1.0);
  dN[15].r =  0.250 * (1.0-s2) * (1.0-t);
  dN[16].r =  0.125 * (1.0+s) * (1.0-t) * (2.0*r+s-t-1.0);
  dN[17].r = -0.500 * r * (1.0+s) * (1.0-t);
  dN[18].r = -0.125 * (1.0+s) * (1.0-t) * (-2.0*r+s-t-1.0);
  dN[19].r = -0.250 * (1.0-s2) * (1.0-t);

  dN[ 0].s = -0.125 * (1.0-r) * (1.0+t) * (-r-2.0*s+t-1.0);
  dN[ 1].s = -0.250 * (1.0-r2) * (1.0+t);
  dN[ 2].s = -0.125 * (1.0+r) * (1.0+t) * (r-2.0*s+t-1.0);
  dN[ 3].s = -0.500 * (1.0+r) * s * (1.0+t);
  dN[ 4].s =  0.125 * (1.0+r) * (1.0+t) * (r+2.0*s+t-1.0);
  dN[ 5].s =  0.250 * (1.0-r2) * (1.0+t);
  dN[ 6].s =  0.125 * (1.0-r) * (1.0+t) * (-r+2.0*s+t-1.0);
  dN[ 7].s = -0.500 * (1.0-r) * s * (1.0+t);
  dN[ 8].s = -0.250 * (1.0-r) * (1-t2);
  dN[ 9].s = -0.250 * (1.0+r) * (1.0-t2);
  dN[10].s =  0.250 * (1.0+r) * (1.0-t2);
  dN[11].s =  0.250 * (1.0-r) * (1.0-t2);
  dN[12].s = -0.125 * (1.0-r) * (1.0-t) * (-r-2.0*s-t-1.0);
  dN[13].s = -0.250 * (1.0-r2) * (1.0-t);
  dN[14].s = -0.125 * (1.0+r) * (1.0-t) * (r-2.0*s-t-1.0);
  dN[15].s = -0.500 * (1.0+r) * s * (1.0-t);
  dN[16].s =  0.125 * (1.0+r) * (1.0-t) * (r+2.0*s-t-1.0);
  dN[17].s =  0.250 * (1.0-r2) * (1.0-t);
  dN[18].s =  0.125 * (1.0-r) * (1.0-t) * (-r+2.0*s-t-1.0);
  dN[19].s = -0.500 * (1.0-r) * s * (1.0-t);

  dN[ 0].t =  0.125 * (1.0-r) * (1.0-s) * (-r-s+2.0*t-1.0);
  dN[ 1].t =  0.250 * (1.0-r2) * (1.0-s);
  dN[ 2].t =  0.125 * (1.0+r) * (1.0-s) * (r-s+2.0*t-1.0);
  dN[ 3].t =  0.250 * (1.0+r) * (1.0-s2);
  dN[ 4].t =  0.125 * (1.0+r) * (1.0+s) * (r+s+2.0*t-1.0);
  dN[ 5].t =  0.250 * (1.0-r2) * (1.0+s);
  dN[ 6].t =  0.125 * (1.0-r) * (1.0+s) * (-r+s+2.0*t-1.0);
  dN[ 7].t =  0.250 * (1.0-r) * (1.0-s2);
  dN[ 8].t = -0.500 * (1.0-r) * (1.0-s) * t;
  dN[ 9].t = -0.500 * (1.0+r) * (1.0-s) * t;
  dN[10].t = -0.500 * (1.0+r) * (1.0+s) * t;
  dN[11].t = -0.500 * (1.0-r) * (1.0+s) * t;
  dN[12].t = -0.125 * (1.0-r) * (1.0-s) * (-r-s-2.0*t-1.0);
  dN[13].t = -0.250 * (1.0-r2) * (1.0-s);
  dN[14].t = -0.125 * (1.0+r) * (1.0-s) * (r-s-2.0*t-1.0);
  dN[15].t = -0.250 * (1.0+r) * (1.0-s2);
  dN[16].t = -0.125 * (1.0+r) * (1.0+s) * (r+s-2.0*t-1.0);
  dN[17].t = -0.250 * (1.0-r2) * (1.0+s);
  dN[18].t = -0.125 * (1.0-r) * (1.0+s) * (-r+s-2.0*t-1.0);
  dN[19].t = -0.250 * (1.0-r) * (1.0-s2);
}

// ================================ GetEdge ================================

int cShapeBrick20 :: GetEdge(int *corner, eShpType *type, int *nnode,
                             cNode **conn)
{
  static int NumEdge = 12;
  static tEdge Edge[] = { {SHAPE_L3_3D, 3, { 0,  1,  2}},
                          {SHAPE_L3_3D, 3, { 2,  3,  4}},
                          {SHAPE_L3_3D, 3, { 4,  5,  6}},
                          {SHAPE_L3_3D, 3, { 6,  7,  0}},
                          {SHAPE_L3_3D, 3, { 0,  8, 12}},
                          {SHAPE_L3_3D, 3, { 2,  9, 14}},
                          {SHAPE_L3_3D, 3, { 4, 10, 16}},
                          {SHAPE_L3_3D, 3, { 6, 11, 18}},
                          {SHAPE_L3_3D, 3, {12, 13, 14}},
                          {SHAPE_L3_3D, 3, {14, 15, 16}},
                          {SHAPE_L3_3D, 3, {16, 17, 18}},
                          {SHAPE_L3_3D, 3, {18, 19, 12}}};

  return(GetShapeEdge(NumEdge, Edge, corner, type, nnode, conn));
}

// ================================ GetFace ================================

int cShapeBrick20 :: GetFace(int *corner, eShpType *type, int *nnode,
                             cNode **conn)
{
  static int NumFace = 6;
  static tFace Face[] = { {SHAPE_Q8_SURF, 4, { 0,  2,  4,  6}, 8, { 0,  1,  2,  3,  4,  5,  6,  7}},
                          {SHAPE_Q8_SURF, 4, { 2,  0, 12, 14}, 8, { 2,  1,  0,  8, 12, 13, 14,  9}},
                          {SHAPE_Q8_SURF, 4, { 4,  2, 14, 16}, 8, { 4,  3,  2,  9, 14, 15, 16, 10}},
                          {SHAPE_Q8_SURF, 4, { 6,  4, 16, 18}, 8, { 6,  5,  4, 10, 16, 17, 18, 11}},
                          {SHAPE_Q8_SURF, 4, {12, 18, 16, 14}, 8, {12, 19, 18, 17, 16, 15, 14, 13}},
                          {SHAPE_Q8_SURF, 4, {18,  6,  0, 12}, 8, {18, 11,  6,  7,  0,  8, 12, 19}}};

  return(GetShapeFace(NumFace, Face, corner, type, nnode, conn));
}


// -------------------------------------------------------------------------
// Class cShapeTet4:
// -------------------------------------------------------------------------

// ============================== cShapeTet4 ===============================

cShapeTet4 :: cShapeTet4(void)
{
  Type = SHAPE_TET4;
  TopType = TETRAHEDRAL_TOPOLOGY;
  NumNode = 4;
  ElmNode = new cNode*[NumNode];
}

// ============================= ~cShapeTet4 ===============================

cShapeTet4 :: ~cShapeTet4(void)
{
}

// ============================== NodeNatCoord =============================

void cShapeTet4 :: NodeNatCoord(sNatCoord *c)
{
   c[0].r = 0.0; c[0].s = 0.0; c[0].t = 0.0;
   c[1].r = 1.0; c[1].s = 0.0; c[1].t = 0.0;
   c[2].r = 0.0; c[2].s = 1.0; c[2].t = 0.0;
   c[3].r = 0.0; c[3].s = 0.0; c[3].t = 1.0;
}

// ================================= ShpFunc ===============================

void cShapeTet4 :: ShpFunc(sNatCoord p, double *N)
{
  double r = p.r;
  double s = p.s;
  double t = p.t;

  N[0] = 1.0 - r - s - t;
  N[1] = r;
  N[2] = s;
  N[3] = t;
}

// =============================== DrvShpRST ===============================

void cShapeTet4 :: DrvShpRST(sNatCoord p, sNatCoord *dN)
{
  double r = p.r;
  double s = p.s;
  double t = p.t;

  dN[0].r = -1.0;
  dN[1].r =  1.0;
  dN[2].r =  0.0;
  dN[3].r =  0.0;

  dN[0].s = -1.0;
  dN[1].s =  0.0;
  dN[2].s =  1.0;
  dN[3].s =  0.0;

  dN[0].t = -1.0;
  dN[1].t =  0.0;
  dN[2].t =  0.0;
  dN[3].t =  1.0;
}

// ================================ GetEdge ================================

int cShapeTet4 :: GetEdge(int *corner, eShpType *type, int *nnode,
                            cNode **conn)
{
  static int NumEdge = 6;
  static tEdge Edge[] = {{SHAPE_L2_3D, 2, {0, 1}},
                         {SHAPE_L2_3D, 2, {1, 2}},
                         {SHAPE_L2_3D, 2, {2, 0}},
                         {SHAPE_L2_3D, 2, {0, 3}},
                         {SHAPE_L2_3D, 2, {1, 3}},
                         {SHAPE_L2_3D, 2, {2, 3}}};

  return(GetShapeEdge(NumEdge, Edge, corner, type, nnode, conn));
}

// ================================ GetFace ================================

int cShapeTet4 :: GetFace(int *corner, eShpType *type, int *nnode,
                            cNode **conn)
{
  static int NumFace = 4;
  static tFace Face[] = {{SHAPE_T3_SURF, 3, {0, 1, 2}, 3, {0, 1, 2}},  // Bottom
                         {SHAPE_T3_SURF, 3, {0, 2, 3}, 3, {0, 2, 3}},  // Left
                         {SHAPE_T3_SURF, 3, {0, 3, 1}, 3, {0, 3, 1}},  // Back
                         {SHAPE_T3_SURF, 3, {1, 3, 2}, 3, {1, 3, 2}}}; // Right

  return(GetShapeFace(NumFace, Face, corner, type, nnode, conn));
}


// -------------------------------------------------------------------------
// Class cShapeTet10:
// -------------------------------------------------------------------------

// ============================== cShapeTet10 ==============================

cShapeTet10 :: cShapeTet10(void)
{
  Type = SHAPE_TET10;
  TopType = TETRAHEDRAL_TOPOLOGY;
  NumNode = 10;
  ElmNode = new cNode*[NumNode];
}

// ============================= ~cShapeTet10 ==============================

cShapeTet10 :: ~cShapeTet10(void)
{
}

// ============================== NodeNatCoord =============================

void cShapeTet10 :: NodeNatCoord(sNatCoord *c)
{
   c[0].r = 0.0; c[0].s = 0.0; c[0].t = 0.0;
   c[1].r = 0.5; c[1].s = 0.0; c[1].t = 0.0;
   c[2].r = 1.0; c[2].s = 0.0; c[2].t = 0.0;
   c[3].r = 0.5; c[3].s = 0.5; c[3].t = 0.0;
   c[4].r = 0.0; c[4].s = 1.0; c[4].t = 0.0;
   c[5].r = 0.0; c[5].s = 0.5; c[5].t = 0.0;
   c[6].r = 0.0; c[6].s = 0.0; c[6].t = 0.5;
   c[7].r = 0.5; c[7].s = 0.0; c[7].t = 0.5;
   c[8].r = 0.5; c[8].s = 0.5; c[8].t = 0.5;
   c[9].r = 0.0; c[9].s = 0.0; c[9].t = 1.0;
}

// ================================= ShpFunc ===============================

void cShapeTet10 :: ShpFunc(sNatCoord p, double *N)
{
  double r = p.r;
  double s = p.s;
  double t = p.t;

  N[0] = (1-r-s-t) * (1-2*r-2*s-2*t);
  N[1] = 4*r*(1-r-s-t);
  N[2] = r * (2*r-1);
  N[3] = 4*r*s;
  N[4] = s * (2*s-1);
  N[5] = 4*s*(1-r-s-t);
  N[6] = 4*t*(1-r-s-t);
  N[7] = 4*t*r;
  N[8] = 4*s*t;
  N[9] = t * (2*t-1);
}

// =============================== DrvShpRST ===============================

void cShapeTet10 :: DrvShpRST(sNatCoord p, sNatCoord *dN)
{
  double r = p.r;
  double s = p.s;
  double t = p.t;

  dN[0].r =  2*t-2*(-t-s-r+1)+2*s+2*r-1;
  dN[1].r =  4*(1-s-t) - 8*r;
  dN[2].r =  4*r - 1;
  dN[3].r =  4*s;
  dN[4].r =  0.0;
  dN[5].r = -4*s;
  dN[6].r = -4*t;
  dN[7].r =  4*t;
  dN[8].r =  0.0;
  dN[9].r =  0.0;

  dN[0].s =  2*t-2*(-t-s-r+1)+2*s+2*r-1;
  dN[1].s =  -4*r;
  dN[2].s =  0.0;
  dN[3].s =  4*r;
  dN[4].s =  4*s-1;
  dN[5].s =  4*(1-r-t) - 8*s;
  dN[6].s = -4*t;
  dN[7].s =  0.0;
  dN[8].s =  4*t;
  dN[9].s =  0.0;

  dN[0].t =  2*t-2*(-t-s-r+1)+2*s+2*r-1;
  dN[1].t =  -4*r;
  dN[2].t =  0.0;
  dN[3].t =  0.0;
  dN[4].t =  0.0;
  dN[5].t =  -4*s;
  dN[6].t =  4*(1-r-s) - 8*t;
  dN[7].t =  4*r;
  dN[8].t =  4*s;
  dN[9].t =  4*t-1;
}

// ================================ GetEdge ================================

int cShapeTet10 :: GetEdge(int *corner, eShpType *type, int *nnode,
                            cNode **conn)
{
  static int NumEdge = 6;
  static tEdge Edge[] = {{SHAPE_L3_3D, 3, {0, 1, 2}},
                         {SHAPE_L3_3D, 3, {2, 3, 4}},
                         {SHAPE_L3_3D, 3, {4, 5, 0}},
                         {SHAPE_L3_3D, 3, {0, 6, 9}},
                         {SHAPE_L3_3D, 3, {2, 7, 9}},
                         {SHAPE_L3_3D, 3, {4, 8, 9}}};

  return(GetShapeEdge(NumEdge, Edge, corner, type, nnode, conn));
}

// ================================ GetFace ================================

int cShapeTet10 :: GetFace(int *corner, eShpType *type, int *nnode,
                            cNode **conn)
{
  static int NumFace = 4;
  static tFace Face[] = {{SHAPE_T6_SURF, 3, {0, 2, 4}, 6, {0, 1, 2, 3, 4, 5}},  // Bottom
                         {SHAPE_T6_SURF, 3, {0, 4, 9}, 6, {0, 5, 4, 8, 9, 6}},  // Left
                         {SHAPE_T6_SURF, 3, {0, 9, 2}, 6, {0, 6, 9, 7, 2, 1}},  // Back
                         {SHAPE_T6_SURF, 3, {4, 2, 9}, 6, {4, 3, 2, 7, 9, 8}}}; // Right

  return(GetShapeFace(NumFace, Face, corner, type, nnode, conn));
}

// ======================================================= End of file =====
