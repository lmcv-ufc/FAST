// -------------------------------------------------------------------------
// shpline.cpp - implementation of the Line Shape class.
// -------------------------------------------------------------------------
// Created:      06-Sep-2005     Aurea Silva de Holanda
//
// Modified:     17-May-2011     Evandro Parente Junior
//               Implementation of non-isoparametric elements.
// -------------------------------------------------------------------------

#include <math.h>

#include "shpline.h"
#include "shape.h"
#include "mat.h"

#ifdef _OMP_
#include "omp.h"
#endif

#include <iostream>

// -------------------------------------------------------------------------
// Local constants:
//
const double LINE_TOL = 1.0e-12;

// -------------------------------------------------------------------------
// Static variables:
//
static double       _J[3][3];   // Jabobian matrix
static double    _InvJ[3][3];   // Inverse of Jacobian matrix
static sNatCoord _dNrst[50];    // Shape func derivatives/parametric coords
static sNatCoord _dMrst[50];    // Mapping func derivatives/parametric coords

#pragma omp threadprivate(_J,_InvJ,_dNrst,_dMrst)

// -------------------------------------------------------------------------
// Class cShapeBar:
// -------------------------------------------------------------------------

// =============================== cShapeBar ===============================

cShapeBar :: cShapeBar(void)
{
  Type    = SHAPE_BAR;
  TopType = LINE_TOPOLOGY;
  NumNode = 2;
  ElmNode = new cNode*[NumNode];
}

// ============================== ~cShapeBar ===============================

cShapeBar :: ~cShapeBar(void)
{
}

// -------------------------------------------------------------------------
// Class cShapeLine2D:
// -------------------------------------------------------------------------

// ============================== cShapeLine2D =============================

cShapeLine2D :: cShapeLine2D(void)
{
  TopType = LINE_TOPOLOGY;
}

// ============================= ~cShapeLine2D =============================

cShapeLine2D :: ~cShapeLine2D(void)
{

}

// =============================== DrvShpXYZ ===============================

void cShapeLine2D :: DrvShpXYZ(sNatCoord pt, sNodeCoord *coord,
                               double *detJ, sNodeCoord *dNxyz)
{
  int i;

  // Compute dNr.

  DrvShpRST(pt, _dNrst);

  // Compute the tangent vector.

  double tx = 0.0;
  double ty = 0.0;
  if (Isoparametric( )) // Isoparametric mapping
  {
    for (i = 0; i < NumNode; i++)
    {
      tx += _dNrst[i].r*coord[i].x;
      ty += _dNrst[i].r*coord[i].y;
    }
  }
  else                  // Non-isoparametric mapping
  {
    DrvMapRST(pt, _dMrst);
    for (i = 0; i < NumNode; i++)
    {
      tx += _dMrst[i].r*coord[i].x;
      ty += _dMrst[i].r*coord[i].y;
    }
  }

  // Compute the [J] matrix.

  _J[0][0] =  tx;
  _J[0][1] =  ty;
  _J[1][0] = -ty;
  _J[1][1] =  tx;

  // Compute |J| and inv(J).

  ShpMatInv(2, _J, _InvJ, detJ);
  *detJ = sqrt(tx*tx + ty*ty);   // |J| = ||t||

  // Compute Ni,xy.

  for (i = 0; i < NumNode; i++)
  {
    dNxyz[i].x = _InvJ[0][0]*_dNrst[i].r;
    dNxyz[i].y = _InvJ[1][0]*_dNrst[i].r;
  }
}

// ================================ LocalSys ===============================

void cShapeLine2D :: LocalSys(sNatCoord pt, sNodeCoord *coord, cMatrix &R)
{
  // Compute dMr.

  DrvMapRST(pt, _dMrst);

  // Compute the unit tangent vector.

  double tx = 0.0;
  double ty = 0.0;
  for (int i = 0; i < NumNode; i++)
  {
    tx += _dMrst[i].r*coord[i].x;
    ty += _dMrst[i].r*coord[i].y;
  }
  double len = sqrt(tx*tx + ty*ty);
  if (len != 0.0)
  {
    tx /= len;
    ty /= len;
  }

  // Assembly the rotation matrix (vectors in columns).

  if (R.NRow( ) == 2)
  {
    R[0][0] =  tx;  R[0][1] = -ty;
    R[1][0] =  ty;  R[1][1] =  tx;
  }
  else
  {
    R[0][0] =  tx;  R[0][1] = -ty;  R[0][2] = 0.0;
    R[1][0] =  ty;  R[1][1] =  tx;  R[1][2] = 0.0;
    R[2][0] = 0.0;  R[2][1] = 0.0;  R[2][2] = 1.0;
  }
}

// ================================ GetEdge ================================

int cShapeLine2D :: GetEdge(int *corner, eShpType *type, int *nnode,
                            cNode **conn)
{
  int last = NumNode - 1;

  // Direct order.

  if (corner[0] == ElmNode[0]->GetLabel( ) &&
      corner[1] == ElmNode[last]->GetLabel( ))
  {
    *type  = Type;
    *nnode = NumNode;
    for (int i = 0; i < NumNode; i++) conn[i] = ElmNode[i];
    return(1);
  }

  // Inverse order.

  if (corner[0] == ElmNode[last]->GetLabel( ) &&
      corner[1] == ElmNode[0]->GetLabel( ))
  {
    *type  = Type;
    *nnode = NumNode;
    for (int i = 0; i < NumNode; i++) conn[i] = ElmNode[last-i];
    return(1);
  }

  // Not found.

  return(0);
}


// -------------------------------------------------------------------------
// Class cShapeLine2:
// -------------------------------------------------------------------------

// ============================== cShapeLine2 ==============================

cShapeLine2 :: cShapeLine2(void)
{
  Type = SHAPE_L2_2D;
  NumNode = 2;
  ElmNode = new cNode*[NumNode];
}

// ============================= ~cShapeLine2 ==============================

cShapeLine2 :: ~cShapeLine2(void)
{
}

// ============================= NodeNatCoord ==============================

void cShapeLine2 :: NodeNatCoord(sNatCoord *c)
{
  c[0].r = -1.0;
  c[1].r =  1.0;
}

// ================================= ShpFunc ===============================

void cShapeLine2 :: ShpFunc(sNatCoord p, double *N)
{
  double r = p.r;

  N[0] = 0.5*(1.0 - r);
  N[1] = 0.5*(1.0 + r);
}

// =============================== DrvShpRST ===============================

void cShapeLine2 :: DrvShpRST(sNatCoord p, sNatCoord *dN)
{
  dN[0].r = -0.5;
  dN[1].r =  0.5;
}


// -------------------------------------------------------------------------
// Class cShapeLine3:
// -------------------------------------------------------------------------

// ============================== cShapeLine3 ==============================

cShapeLine3 :: cShapeLine3(void)
{
  Type = SHAPE_L3_2D;
  NumNode = 3;
  ElmNode = new cNode*[NumNode];
}

// ============================= ~cShapeLine3 ==============================

cShapeLine3 :: ~cShapeLine3(void)
{
}

// ============================== NodeNatCoord =============================

void cShapeLine3 :: NodeNatCoord(sNatCoord *c)
{
  c[0].r = -1.0;
  c[1].r =  0.0;
  c[2].r =  1.0;
}

// ================================= ShpFunc ===============================

void cShapeLine3 :: ShpFunc(sNatCoord p, double *N)
{
  double r = p.r;

  N[0] = 0.5*r*(r - 1.0);
  N[1] = 1.0 - r*r;
  N[2] = 0.5*r*(r + 1.0);
}

// =============================== DrvShpRST ===============================

void cShapeLine3 :: DrvShpRST(sNatCoord p, sNatCoord *dN)
{
  double r = p.r;

  dN[0].r =  0.5*(2.0*r - 1.0);
  dN[1].r = -2.0*r;
  dN[2].r =  0.5*(2.0*r + 1.0);
}

// -------------------------------------------------------------------------
// Class cShapeLine3D:
// -------------------------------------------------------------------------

// ============================== cShapeLine3D =============================

cShapeLine3D :: cShapeLine3D(void)
{
  TopType = LINE_TOPOLOGY;
}

// ============================= ~cShapeLine3D =============================

cShapeLine3D :: ~cShapeLine3D(void)
{
}

// =============================== DrvShpXYZ ===============================

void cShapeLine3D :: DrvShpXYZ(sNatCoord pt, sNodeCoord *coord,
                               double *detJ, sNodeCoord *dNxyz)
{
  int i;

  // Compute dNr.

  DrvShpRST(pt, _dNrst);

  // Evaluate the first line (tangent vector) => {u}.

  ShpMatZero(3, _J);
  if (Isoparametric( )) // Isoparametric mapping
  {
    for (i = 0; i < NumNode; i++)
    {
      _J[0][0] += _dNrst[i].r*coord[i].x;
      _J[0][1] += _dNrst[i].r*coord[i].y;
      _J[0][2] += _dNrst[i].r*coord[i].z;
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
    }
  }

  // Verify if the jacobian valid.

  double len = ShpVecLen3D(_J[0]);
  if (len < LINE_TOL)
  {
    *detJ = 0.0;
    return;
  }

  // Determine a vector orthogonal to the tangent vector => {v}.

  double ez[3] = { 0.0, 0.0, 1.0 };
  ShpCrossProd(ez, _J[0], _J[1]);
  if ((ShpVecLen3D(_J[1])/len) < LINE_TOL)
  {
    double ey[3] = { 0.0, 1.0, 0.0 };
    ShpCrossProd(ey, _J[0], _J[1]);
  }

  // Evaluate the third line => {w} = {u} x {v}.

  ShpCrossProd(_J[0], _J[1], _J[2]);

  // Numerically invert the Jacobian.

  ShpMatInv(3, _J, _InvJ, detJ);

  // Evaluate the determinant as the length of the tangent vector.

  *detJ = len;

  // Compute Ni,xyz.

  for (i = 0; i < NumNode; i++)
  {
    dNxyz[i].x = _InvJ[0][0]*_dNrst[i].r;
    dNxyz[i].y = _InvJ[1][0]*_dNrst[i].r;
    dNxyz[i].z = _InvJ[2][0]*_dNrst[i].r;
  }
}

// ================================ LocalSys ===============================

void cShapeLine3D :: LocalSys(sNatCoord pt, sNodeCoord *coord, cMatrix &R)
{
  int i;
  double u[3] = { 0.0, 0.0, 0.0 };
  double v[3] = { 0.0, 0.0, 0.0 };
  double w[3] = { 0.0, 0.0, 0.0 };

  // Compute dMr.

  DrvMapRST(pt, _dMrst);

  // Evaluate the first line (tangent vector) => {u}.

  for (i = 0; i < NumNode; i++)
  {
    u[0] += _dMrst[i].r*coord[i].x;
    u[1] += _dMrst[i].r*coord[i].y;
    u[2] += _dMrst[i].r*coord[i].z;
  }

  // Normalize the tangent vector.

  double len = ShpVecLen3D(u);
  if (len < LINE_TOL) return;
  ShpNormVec3D(u);

  // Determine a vector orthogonal to the tangent vector => {v}.

  double ez[3] = { 0.0, 0.0, 1.0 };
  ShpCrossProd(ez, u, v);
  if ((ShpVecLen3D(v)/len) < LINE_TOL)
  {
    double ey[3] = { 0.0, 1.0, 0.0 };
    ShpCrossProd(ey, u, v);
  }

  // Evaluate the third line => {w} = {u} x {v}.

  ShpCrossProd(u, v, w);

  // Assembly the rotation matrix (vectors in columns).

  for (i = 0; i < 3; i++)
  {
    R[i][0] = u[i];
    R[i][1] = v[i];
    R[i][2] = w[i];
  }
}

// -------------------------------------------------------------------------
// Class cShape3DLine2:
// -------------------------------------------------------------------------

// ============================= cShape3DLine2 =============================

cShape3DLine2 :: cShape3DLine2(void)
{
  Type = SHAPE_L2_3D;
}

// ============================ ~cShape3DLine2 =============================

cShape3DLine2 :: ~cShape3DLine2(void)
{
}


// -------------------------------------------------------------------------
// Class cShape3DLine3:
// -------------------------------------------------------------------------

// ============================= cShape3DLine3 =============================

cShape3DLine3 :: cShape3DLine3(void)
{
  Type = SHAPE_L3_3D;
}

// ============================ ~cShape3DLine3 =============================

cShape3DLine3 :: ~cShape3DLine3(void)
{
}

// ======================================================= End of file =====
