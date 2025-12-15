// ------------------------------------------------------------------------
// shpsurf.cpp - This file contains methods of Curve shape class.
// -------------------------------------------------------------------------
// Created:      06-Sep-2005     Aurea Silva de Holanda
//
// Modified:     17-May-2011     Evandro Parente Junior
//               Implementation of non-isoparametric elements.
//
// Modified:     02-Oct-2013     Evandro Parente Junior
//               Implementation of tVector. Use of DrvMapRST in LocalSys.
// -------------------------------------------------------------------------

#include <math.h>
#include <stdio.h>
#include <iostream>

using namespace std;

#include "shpsurf.h"
#include "shape.h"
#include "vec.h"
#include "mat.h"

#ifdef _OMP_
#include "omp.h"
#endif

// -------------------------------------------------------------------------
// Static variables:
//
static double       _J[3][3];   // Jabobian matrix
static double    _InvJ[3][3];   // Inverse of Jacobian matrix
static sNatCoord _dNrst[500];    // Shape func derivatives/parametric coords
static sNatCoord _dMrst[500];    // Mapping func derivatives/parametric coords

#pragma omp threadprivate(_J,_InvJ,_dNrst,_dMrst)

// -------------------------------------------------------------------------
// Public methods:
//

// ============================== cShapeSurf ===============================

cShapeSurf :: cShapeSurf(void)
{
}

// ============================= ~cShapeSurf ===============================

cShapeSurf :: ~cShapeSurf(void)
{
}

// ============================== DrvShpXYZ ================================

void cShapeSurf :: DrvShpXYZ(sNatCoord pt, sNodeCoord *coord,
                             double *detJ, sNodeCoord *dNxyz)
{
  int i;

  // Compute dNrs.

  DrvShpRST(pt, _dNrst);

  // Compute [J].
  // Evaluate the first and second rows (tangent vectors).

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
    }
  }

  // Get third row (normal vector) by cross product of the tangent vectors.

  ShpCrossProd(_J[0], _J[1], _J[2]);

  // Numerically invert the Jacobian matrix.

  ShpMatInv(3, _J, _InvJ, detJ);

  // Evaluate the |J| as the length of the normal vector.

  *detJ = ShpVecLen3D(_J[2]);

  // Compute Ni,xyz.

  for (i = 0; i < NumNode; i++)
  {
    dNxyz[i].x = _InvJ[0][0]*_dNrst[i].r + _InvJ[0][1]*_dNrst[i].s;
    dNxyz[i].y = _InvJ[1][0]*_dNrst[i].r + _InvJ[1][1]*_dNrst[i].s;
    dNxyz[i].z = _InvJ[2][0]*_dNrst[i].r + _InvJ[2][1]*_dNrst[i].s;
  }
}

// ============================== LocalSys =================================

void cShapeSurf :: LocalSys(sNatCoord pt, sNodeCoord *coord, cMatrix &R)
{
  int i;

  // Compute dMrs.

  DrvMapRST(pt, _dMrst);

  // Evaluate the tangent vectors: {u} = d{x}/dr and {v} = d{x}/ds.

  cVector u(3),v(3);
  u.Zero( );
  v.Zero( );
  for (i = 0; i < NumNode; i++)
  {
    u[0] += _dMrst[i].r*coord[i].x;
    u[1] += _dMrst[i].r*coord[i].y;
    u[2] += _dMrst[i].r*coord[i].z;

    v[0] += _dMrst[i].s*coord[i].x;
    v[1] += _dMrst[i].s*coord[i].y;
    v[2] += _dMrst[i].s*coord[i].z;
  }
  u.Normalize( );
  v.Normalize( );

  // Evaluate the normal vector {w} = {u} x {v}.

  cVector w(3);
  CrossProd(u, v, w);

  // Ensure that the 3 vectors are orthogonal to each other.

  CrossProd(w, u, v); // {v} = {w} x {u}

  // Assembly the rotation matrix (vectors in columns).

  for (i = 0; i < 3; i++)
  {
    R[i][0] = u[i];
    R[i][1] = v[i];
    R[i][2] = w[i];
  }
}

// ============================== tVector ==================================

void cShapeSurf :: tVector(sNatCoord pt, sNodeCoord *coord, cVector &w)
{
  // Compute dMrs.

  DrvMapRST(pt, _dMrst);

  // Evaluate the tangent vectors ({u} => r and {v} => s).

  cVector u(3),v(3);
  u.Zero( );
  v.Zero( );
  for (int i = 0; i < NumNode; i++)
  {
    u[0] += _dMrst[i].r*coord[i].x;
    u[1] += _dMrst[i].r*coord[i].y;
    u[2] += _dMrst[i].r*coord[i].z;

    v[0] += _dMrst[i].s*coord[i].x;
    v[1] += _dMrst[i].s*coord[i].y;
    v[2] += _dMrst[i].s*coord[i].z;
  }
  u.Normalize( );
  v.Normalize( );

  // Evaluate the normal vector {w} = {u} x {v}.
  
  CrossProd(u, v, w);
  w.Normalize( );       // Return the unit vector {w}/|w|
}

// -------------------------------------------------------------------------
// Public methods (class cShapeT3Surf):
//

// ============================ cShapeT3Surf ===============================

cShapeT3Surf :: cShapeT3Surf(void)
{
  Type = SHAPE_T3_SURF;
}

// ============================ ~cShapeT3Surf ==============================

cShapeT3Surf :: ~cShapeT3Surf(void)
{
}

// ================================ GetEdge ================================

int cShapeT3Surf :: GetEdge(int *corner, eShpType *type, int *nnode,
                            cNode **conn)
{
  static int NumEdge = 3;
  static tEdge Edge[] = { {SHAPE_L2_3D, 2, {0, 1}},
                          {SHAPE_L2_3D, 2, {1, 2}},
                          {SHAPE_L2_3D, 2, {2, 1}}};

  return(GetShapeEdge(NumEdge, Edge, corner, type, nnode, conn));
}

// -------------------------------------------------------------------------
// Public methods (class cShapeT6Surf):
//

// ============================ cShapeT6Surf ===============================

cShapeT6Surf :: cShapeT6Surf(void)
{
  Type = SHAPE_T6_SURF;
}

// ============================ ~cShapeT6Surf ==============================

cShapeT6Surf :: ~cShapeT6Surf(void)
{

}

// ================================ GetEdge ================================

int cShapeT6Surf :: GetEdge(int *corner, eShpType *type, int *nnode,
                        cNode **conn)
{
  static int NumEdge = 3;
  static tEdge Edge[] = { {SHAPE_L3_3D, 3, {0, 1, 2}},
                          {SHAPE_L3_3D, 3, {2, 3, 4}},
                          {SHAPE_L3_3D, 3, {4, 5, 0}}};

  return(GetShapeEdge(NumEdge, Edge, corner, type, nnode, conn));
}

// -------------------------------------------------------------------------
// Public methods (class cShapeQ4Surf):
//

// ============================ cShapeQ4Surf ===============================

cShapeQ4Surf :: cShapeQ4Surf(void)
{
  Type = SHAPE_Q4_SURF;
}

// ============================ ~cShapeQ4Surf ==============================

cShapeQ4Surf :: ~cShapeQ4Surf(void)
{

}

// ================================ GetEdge ================================

int cShapeQ4Surf :: GetEdge(int *corner, eShpType *type, int *nnode,
                            cNode **conn)
{
  static int NumEdge = 4;
  static tEdge Edge[] = { {SHAPE_L2_3D, 2, {0, 1}},
                          {SHAPE_L2_3D, 2, {1, 2}},
                          {SHAPE_L2_3D, 2, {2, 3}},
                          {SHAPE_L2_3D, 2, {3, 0}}};

  return(GetShapeEdge(NumEdge, Edge, corner, type, nnode, conn));
}

// -------------------------------------------------------------------------
// Public methods (class cShapeQ8Surf):
//

// ============================ cShapeQ8Surf ===============================

cShapeQ8Surf :: cShapeQ8Surf(void)
{
  Type = SHAPE_Q8_SURF;
}

// =========================== ~cShapeQ8Surf ===============================

cShapeQ8Surf :: ~cShapeQ8Surf(void)
{

}

// ================================ GetEdge ================================

int cShapeQ8Surf :: GetEdge(int *corner, eShpType *type, int *nnode,
                            cNode **conn)
{
  static int NumEdge = 4;
  static tEdge Edge[] = { {SHAPE_L3_3D, 3, {0, 1, 2}},
                          {SHAPE_L3_3D, 3, {2, 3, 4}},
                          {SHAPE_L3_3D, 3, {4, 5, 6}},
                          {SHAPE_L3_3D, 3, {6, 7, 0}}};

  return(GetShapeEdge(NumEdge, Edge, corner, type, nnode, conn));
}

// -------------------------------------------------------------------------
// Public methods (class cShapeQ9Surf):
//

// ============================ cShapeQ9Surf ===============================

cShapeQ9Surf :: cShapeQ9Surf(void)
{
  Type = SHAPE_Q9_SURF;
}

// =========================== ~cShapeQ9Surf ===============================

cShapeQ9Surf :: ~cShapeQ9Surf(void)
{

}

// ================================ GetEdge ================================

int cShapeQ9Surf :: GetEdge(int *corner, eShpType *type, int *nnode,
                            cNode **conn)
{
  static int NumEdge = 4;
  static tEdge Edge[] = { {SHAPE_L3_3D, 3, {0, 1, 2}},
                          {SHAPE_L3_3D, 3, {2, 3, 4}},
                          {SHAPE_L3_3D, 3, {4, 5, 6}},
                          {SHAPE_L3_3D, 3, {6, 7, 0}}};

  return(GetShapeEdge(NumEdge, Edge, corner, type, nnode, conn));
}

// ====================================================== End of File ======
