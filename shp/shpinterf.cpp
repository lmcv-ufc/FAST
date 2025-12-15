// -------------------------------------------------------------------------
// shpinterf.cpp - implementation of the Shape class.
// -------------------------------------------------------------------------
// Created:      16-Dez-2013     Edson Moreira Dantas Junior
//
// Modified:     06-May-2014     Evandro Parente Junior
//               Creation of two-dimensional interface shapes.
// -------------------------------------------------------------------------

#include <math.h>
#include <stdio.h>

#include "shpinterf.h"
#include "vec.h"

// -------------------------------------------------------------------------
// Static variables:
//
static sNatCoord _dMrst[50];    // Mapping func derivatives/parametric coords

// -------------------------------------------------------------------------
// Class cShapeInterf2D:
// -------------------------------------------------------------------------

// ============================= cShapeInterf2D ============================

cShapeInterf2D :: cShapeInterf2D(void)
{
  TopType = QUADRILATERAL_TOPOLOGY;
}

// ============================ ~cShapeInterf2D ============================

cShapeInterf2D :: ~cShapeInterf2D(void)
{
}

// =============================== RMatrix =================================

void cShapeInterf2D :: RMatrix(sNatCoord pt, sNodeCoord *coord,
                               double *detJ, cMatrix &R)
{
  int i;

  // Compute dMrs.

  DrvMapRST(pt, _dMrst);

  // Evaluate the tangent vector {t} = d{x}/dr.

  double tx = 0.0;
  double ty = 0.0;
  for (i = 0; i < NumNode; i++)
  {
    tx += _dMrst[i].r*coord[i].x/2.0;
    ty += _dMrst[i].r*coord[i].y/2.0;
  }

  // Normalize the tangent vector.

  double len = sqrt(tx*tx + ty*ty);
  if (len != 0.0)
  {
    tx /= len;
    ty /= len;
  }
  *detJ = len;

  // Assembly the rotation matrix (vectors in rows).

  R[0][0] =  tx;  R[0][1] = ty;
  R[1][0] = -ty;  R[1][1] = tx;
}


// -------------------------------------------------------------------------
// Class cShapeInterfL2:
// -------------------------------------------------------------------------

// ============================ cShapeInterfL2 =============================

cShapeInterfL2 :: cShapeInterfL2(void)
{
  Type = SHAPE_INTERF_L2;
  NumNode = 4;
  ElmNode = new cNode*[NumNode];
}

// =========================== ~cShapeInterfL2 =============================

cShapeInterfL2 :: ~cShapeInterfL2(void)
{
}

// ============================= NodeNatCoord ==============================

void cShapeInterfL2 :: NodeNatCoord(sNatCoord *c)
{
  c[0].r = c[2].r = -1.0;
  c[1].r = c[3].r =  1.0;
}

// ================================= ShpFunc ===============================

void cShapeInterfL2 :: ShpFunc(sNatCoord p, double *N)
{
  double r = p.r;

  N[0] = N[2] = 0.5*(1.0 - r);
  N[1] = N[3] = 0.5*(1.0 + r);
}

// =============================== DrvShpRST ===============================

void cShapeInterfL2 :: DrvShpRST(sNatCoord p, sNatCoord *dN)
{
  dN[0].r = dN[2].r = -0.5;
  dN[1].r = dN[3].r =  0.5;
}


// -------------------------------------------------------------------------
// Class cShapeInterfL3:
// -------------------------------------------------------------------------

// ============================= cShapeInterfL3 ============================

cShapeInterfL3 :: cShapeInterfL3(void)
{
  Type = SHAPE_INTERF_L3;
  NumNode = 6;
  ElmNode = new cNode*[NumNode];
}

// ============================ ~cShapeInterfL3 ============================

cShapeInterfL3 :: ~cShapeInterfL3(void)
{
  TopType = QUADRILATERAL_TOPOLOGY;
}

// ============================== NodeNatCoord =============================

void cShapeInterfL3 :: NodeNatCoord(sNatCoord *c)
{
  c[0].r = c[3].r = -1.0;
  c[1].r = c[4].r =  0.0;
  c[2].r = c[5].r =  1.0;
}

// ================================= ShpFunc ===============================

void cShapeInterfL3 :: ShpFunc(sNatCoord p, double *N)
{
  double r = p.r;

  N[0] = N[3] = 0.5*r*(r - 1.0);
  N[1] = N[4] = 1.0 - r*r;
  N[2] = N[5] = 0.5*r*(r + 1.0);
}

// =============================== DrvShpRST ===============================

void cShapeInterfL3 :: DrvShpRST(sNatCoord p, sNatCoord *dN)
{
  double r = p.r;

  dN[0].r = dN[3].r =  0.5*(2.0*r - 1.0);
  dN[1].r = dN[4].r = -2.0*r;
  dN[2].r = dN[5].r =  0.5*(2.0*r + 1.0);
}


// -------------------------------------------------------------------------
// Class cShapeInterf3D:
// -------------------------------------------------------------------------

// ============================= cShapeInterf3D ============================

cShapeInterf3D :: cShapeInterf3D(void)
{
  TopType = HEXAHEDRAL_TOPOLOGY;
}

// ============================ ~cShapeInterf3D ============================

cShapeInterf3D :: ~cShapeInterf3D(void)
{
}

// =============================== RMatrix =================================

void cShapeInterf3D :: RMatrix(sNatCoord pt, sNodeCoord *coord,
                               double *detJ, cMatrix &R)
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
    u[0] += _dMrst[i].r*coord[i].x/2.0;
    u[1] += _dMrst[i].r*coord[i].y/2.0;
    u[2] += _dMrst[i].r*coord[i].z/2.0;

    v[0] += _dMrst[i].s*coord[i].x/2.0;
    v[1] += _dMrst[i].s*coord[i].y/2.0;
    v[2] += _dMrst[i].s*coord[i].z/2.0;
  }

  // Evaluate the normal vector {w} = {u} x {v}.

  cVector w(3);
  CrossProd(u, v, w);
  *detJ = w.Length( );

  // Ensure that the 3 vectors are orthogonal to each other.

  CrossProd(w, u, v); // {v} = {w} x {u}

  // Assembly the rotation matrix (vectors in rows).

  u.Normalize( );
  v.Normalize( );
  w.Normalize( );
  for (i = 0; i < 3; i++)
  {
    R[0][i] = u[i];
    R[1][i] = v[i];
    R[2][i] = w[i];
  }
}


// -------------------------------------------------------------------------
// Class cShapeInterfQ4:
// -------------------------------------------------------------------------

// ============================ cShapeInterfQ4 =============================

cShapeInterfQ4 :: cShapeInterfQ4(void)
{
  Type = SHAPE_INTERF_Q4;
  NumNode = 8;
  ElmNode = new cNode*[NumNode];
}

// =========================== ~cShapeInterfQ4 =============================

cShapeInterfQ4  :: ~cShapeInterfQ4(void)
{
}

// ============================== NodeNatCoord =============================

void cShapeInterfQ4 :: NodeNatCoord(sNatCoord *c)
{
  c[0].r = c[4].r = -1.0;  c[0].s = c[4].s = -1.0;
  c[1].r = c[5].r =  1.0;  c[1].s = c[5].s = -1.0;
  c[2].r = c[6].r =  1.0;  c[2].s = c[6].s =  1.0;
  c[3].r = c[7].r = -1.0;  c[3].s = c[7].s =  1.0;
}

// ================================= ShpFunc ===============================

void cShapeInterfQ4 :: ShpFunc(sNatCoord p, double *N)
{
  double r = p.r;
  double s = p.s;

  N[0] = N[4] = 0.25*(1.0 - r)*(1.0 - s);
  N[1] = N[5] = 0.25*(1.0 + r)*(1.0 - s);
  N[2] = N[6] = 0.25*(1.0 + r)*(1.0 + s);
  N[3] = N[7] = 0.25*(1.0 - r)*(1.0 + s);
}

// =============================== DrvShpRST ===============================

void cShapeInterfQ4 :: DrvShpRST(sNatCoord p, sNatCoord *dN)
{
  double r = p.r;
  double s = p.s;

  dN[0].r = dN[4].r = -0.25*(1.0 - s);
  dN[1].r = dN[5].r =  0.25*(1.0 - s);
  dN[2].r = dN[6].r =  0.25*(1.0 + s);
  dN[3].r = dN[7].r = -0.25*(1.0 + s);

  dN[0].s = dN[4].s = -0.25*(1.0 - r);
  dN[1].s = dN[5].s = -0.25*(1.0 + r);
  dN[2].s = dN[6].s =  0.25*(1.0 + r);
  dN[3].s = dN[7].s =  0.25*(1.0 - r);
}


// -------------------------------------------------------------------------
// Class cShapeInterfQ8:
// -------------------------------------------------------------------------

// ============================ cShapeInterfQ8 =============================

cShapeInterfQ8 :: cShapeInterfQ8(void)
{
  Type = SHAPE_INTERF_Q8;
  NumNode = 16;
  ElmNode = new cNode*[NumNode];
}

// =========================== ~cShapeInterfQ8 =============================

cShapeInterfQ8  :: ~cShapeInterfQ8(void)
{
}

// ============================== NodeNatCoord =============================

void cShapeInterfQ8 :: NodeNatCoord(sNatCoord *c)
{
  c[0].r = c[ 8].r = -1.0;  c[0].s = c[ 8].s = -1.0;
  c[1].r = c[ 9].r =  0.0;  c[1].s = c[ 9].s = -1.0;
  c[2].r = c[10].r =  1.0;  c[2].s = c[10].s = -1.0;
  c[3].r = c[11].r =  1.0;  c[3].s = c[11].s =  0.0;
  c[4].r = c[12].r =  1.0;  c[4].s = c[12].s =  1.0;
  c[5].r = c[13].r =  0.0;  c[5].s = c[13].s =  1.0;
  c[6].r = c[14].r = -1.0;  c[6].s = c[14].s =  1.0;
  c[7].r = c[15].r = -1.0;  c[7].s = c[15].s =  0.0;
}

// ================================= ShpFunc ===============================

void cShapeInterfQ8 :: ShpFunc(sNatCoord p, double *N)
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

  N[ 8] = N[0];
  N[ 9] = N[1];
  N[10] = N[2];
  N[11] = N[3];
  N[12] = N[4];
  N[13] = N[5];
  N[14] = N[6];
  N[15] = N[7];
}

// =============================== DrvShpRST ===============================

void cShapeInterfQ8 :: DrvShpRST(sNatCoord p, sNatCoord *dN)
{
  double r = p.r;
  double s = p.s;

  dN[0].r = dN[ 8].r = (2.0*r - 2.0*r*s - s*s + s)/4.0;
  dN[1].r = dN[ 9].r = r*s - r;
  dN[2].r = dN[10].r = (2.0*r - 2.0*r*s + s*s - s)/4.0;
  dN[3].r = dN[11].r = (1.0 - s*s)/2.0;
  dN[4].r = dN[12].r = (2.0*r + 2.0*r*s + s*s + s)/4.0;
  dN[5].r = dN[13].r = -r - r*s;
  dN[6].r = dN[14].r = (2.0*r + 2.0*r*s - s*s - s)/4.0;
  dN[7].r = dN[15].r = (s*s - 1.0)/2.0;

  dN[0].s = dN[ 8].s = (2.0*s - r*r - 2.0*r*s + r)/4.0;
  dN[1].s = dN[ 9].s = (r*r - 1.0)/2.0;
  dN[2].s = dN[10].s = (2.0*s - r*r + 2.0*r*s - r)/ 4.0;
  dN[3].s = dN[11].s = -s - r*s;
  dN[4].s = dN[12].s = (2.0*s + r*r + 2.0*r*s + r)/4.0;
  dN[5].s = dN[13].s = (1.0 - r*r)/2.0;
  dN[6].s = dN[14].s = (2.0*s + r*r - 2.0*r*s - r)/4.0;
  dN[7].s = dN[15].s = r*s - s;
}

// ======================================================= End of file =====
