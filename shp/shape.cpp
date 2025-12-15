// -------------------------------------------------------------------------
// shape.cpp - implementation of the Shape class.
// -------------------------------------------------------------------------
// Created:      30-Apr-2005     Evandro Parente Junior
//
// Modified:     02-Oct-2013     Evandro Parente Junior
//               Creation of tVector method.
//
// Modified:     08-Apr-2015     Elias Saraiva Barroso
//               Input/output operations using c++ streams.
// -------------------------------------------------------------------------

#include <iostream>
#include <math.h>
#include <stdlib.h>

using namespace std;

#include "shape.h"
#include "shpline.h"
#include "shpplane.h"
#include "shpsolid.h"
#include "shpsurf.h"
#include "shpshell.h"
#include "shpinterf.h"
#include "shpinf.h"
#include "shpbsp.h"
#include "shpbez.h"
#include "shpinterfiga.h"
#include "node.h"
#include "vec.h"
#include "mat.h"
#include "gblvar.h"

// -------------------------------------------------------------------------
// Public methods:
//

// ============================== CreateShape ==============================

cShape *cShape :: CreateShape(eShpType type)
{
  cShape *shp = 0;

  switch(type)
  {
    case SHAPE_BAR:
     shp = new cShapeBar( );
    break;

    case SHAPE_L2_2D:
     shp = new cShapeLine2( );
    break;

    case SHAPE_L3_2D:
     shp = new cShapeLine3( );
    break;

    case SHAPE_L2_3D:     //  2-noded line in 3D space
     shp = new cShape3DLine2( );
    break;

    case SHAPE_L3_3D:     //  3-noded line in 3D space
     shp = new cShape3DLine3( );
    break;

    case SHAPE_BEZ_PLCURV:
     shp = new cShapePlBezCurve( );
    break;

    case SHAPE_BEZ_CURV:
     shp = new cShapeBezCurve( );
    break;

    case SHAPE_BSP_PLCURV:
     shp = new cShapePlBspCurve( );
    break;

    case SHAPE_BSP_CURV:
     shp = new cShapeBspCurve( );
    break;

    case SHAPE_T3:
     shp = new cShapeT3( );
    break;

    case SHAPE_T3_SURF:
     shp = new cShapeT3Surf( );
    break;

    case SHAPE_T6:
     shp = new cShapeT6( );
    break;

    case SHAPE_T6_SURF:
     shp = new cShapeT6Surf( );
    break;

    case SHAPE_Q4:
     shp = new cShapeQ4( );
    break;

    case SHAPE_Q4_SURF:
     shp = new cShapeQ4Surf( );
    break;

    case SHAPE_Q8:
     shp = new cShapeQ8( );
    break;

    case SHAPE_Q8_SURF:
     shp = new cShapeQ8Surf( );
    break;

    case SHAPE_Q9_SURF:
     shp = new cShapeQ9Surf( );
    break;

    case SHAPE_Q9:
     shp = new cShapeQ9( );
    break;

    case SHAPE_BSP_PLSURF:
     shp = new cShapePlBspSurf( );
    break;

    case SHAPE_BSP_SURF:
     shp = new cShapeBspSurf( );
    break;

    case SHAPE_BEZ_PLSURF:
     shp = new cShapePlBezSurf( );
    break;

    case SHAPE_BEZ_SURF:
     shp = new cShapeBezSurf( );
    break;

    case SHAPE_BEZTRI_PLSURF:
     shp = new cShapePlBezTrianSurf( );
    break;

    case SHAPE_BEZTRI_SURF:
     shp = new cShapeBezTrianSurf( );
    break;

    case SHAPE_BSP_SOLID:
     shp = new cShapeBspSolid( );
    break;

    case SHAPE_BRICK8:
     shp = new cShapeBrick8( );
    break;

    case SHAPE_BRICK20:
     shp = new cShapeBrick20( );
    break;

    case SHAPE_TET4:
     shp = new cShapeTet4( );
    break;

    case SHAPE_TET10:
     shp = new cShapeTet10( );
    break;

    case SHAPE_INTERF_L2:
     shp = new cShapeInterfL2( );
    break;

    case SHAPE_BSP_INTERF_LINE:
     shp = new cShapeInterflineBsp( );
    break;

    case SHAPE_INTERF_L3:
     shp = new cShapeInterfL3( );
    break;

    case SHAPE_INTERF_Q4:
     shp = new cShapeInterfQ4( );
    break;

    case SHAPE_INTERF_Q8:
     shp = new cShapeInterfQ8( );
    break;

    case SHAPE_L6_INF:
     shp = new cShapeL6Inf( );
    break;

    case SHAPE_T6_SHELL:
    case SHAPE_Q8_SHELL:
    case SHAPE_Q9_SHELL:
      shp = new cShapeShell(type);
    break;

    case SHAPE_BEZTRI_SHELL:
    case SHAPE_BEZ_SHELL:
    case SHAPE_BSP_SHELL:
      shp = new cShapeShellIGA(type);
    break;

/*  case SHAPE_L4_INF:
   shp = new cShapeL4Inf( );
    break;*/
  }

  return(shp);
}

// ================================= cShape ================================

cShape :: cShape(void)
{
}

// ================================ ~cShape ================================

cShape :: ~cShape(void)
{
  delete []ElmNode;
}

// ================================= Read ==================================

void cShape :: Read(void)
{
  int nid;
  for (int i = 0; i < NumNode; i++)
  {
    in >> nid;
    ElmNode[i] = cNode :: FindNode(nid);
    if (!ElmNode[i])
    {
      cout << "Error in the input of element node = " << nid << "!\n";
      exit(0);
    }
  }
}

// =============================== SetNodes ================================

void cShape :: SetNodes(cNode **node, int size)
{
  if (size == 0)
    size = NumNode;

  for (int i = 0; i < size; i++)
  {
    ElmNode[i] = node[i];
    if (!ElmNode[i])
    {
      cout << "Error: Invalid node in SetNodes function!\n";
      exit(0);
    }
  }
}

// ============================== NodalCoord ===============================

void cShape :: NodalCoord(sNodeCoord *coord)
{
  for (int i = 0; i < NumNode; i++) coord[i] = ElmNode[i]->GetCoord( );
}

// ============================== UpdatedNodalCoord ========================

void cShape :: UpdatedNodalCoord(sNodeCoord *coord)
{
  for (int i = 0; i < NumNode; i++)
  {
    coord[i].x = ElmNode[i]->GetCoord( ).x + ElmNode[i]->GetDispl(0); 
    coord[i].y = ElmNode[i]->GetCoord( ).y + ElmNode[i]->GetDispl(1); 
    coord[i].z = ElmNode[i]->GetCoord( ).z + ElmNode[i]->GetDispl(2);
  }
}

// ================================== Evaluate =============================

void cShape :: Evaluate(sNatCoord p, double *N, sNodeCoord *coord, sNodeCoord &c)
{
  // Compute shape functions.
  ShpFunc(p,N);

  // Get nodal coordinates.
  NodalCoord(coord);

  // Evaluate parametric point in physical space.
  c.x = c.y = c.z = 0.0;
  for(int n = 0; n < NumNode; ++n)
  {
    c.x += N[n] * coord[n].x;
    c.y += N[n] * coord[n].y;
    c.z += N[n] * coord[n].z;
  }
}

// ================================ NMatrix ================================

void cShape :: NMatrix(int ndofn, double *shpfunc, cMatrix &N)
{
  N.Zero( );
  for (int i = 0; i < ndofn; i++)
    for (int j = 0; j < NumNode; j++) 
      N[i][ndofn*j+i] = shpfunc[j];
}

// ================================ LocalSys ===============================

void cShape :: LocalSys(sNatCoord p, sNodeCoord *coord, cMatrix &R)
{
  // Identity matrix.

  if (R.NRow( ) == 2)
  {
    R[0][0] = 1.0;  R[0][1] = 0.0;
    R[1][0] = 0.0;  R[1][1] = 1.0;
  }
  else
  {
    R[0][0] = 1.0;  R[0][1] = 0.0;  R[0][2] = 0.0;
    R[1][0] = 0.0;  R[1][1] = 1.0;  R[1][2] = 0.0;
    R[2][0] = 0.0;  R[2][1] = 0.0;  R[2][2] = 1.0;
  }
}

// ================================ RMatrix ================================

void cShape :: RMatrix(sNatCoord p, sNodeCoord *coord, double *detJ, cMatrix &R)
{
  cout << "Invalid method - cShape::LocalSys\n";
  exit(0);
}

// ================================ tVector ================================

void cShape :: tVector(sNatCoord p, sNodeCoord *coord, cVector &v)
{
  // Default out-of-plane vector parallel to z-axis.

  v[0] = 0.0;
  v[1] = 0.0;
  v[2] = 1.0;
}

// -------------------------------------------------------------------------
// Protected methods:
//

// ============================== ShpMatZero ===============================

void cShape :: ShpMatZero(int n, double A[3][3])
{
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++) A[i][j] = 0.0;
}

// ============================== ShpMatInv ================================

int cShape :: ShpMatInv(int n, double A[3][3], double B[3][3], double *det)
{
  if (n == 2)
  {
    (*det) = A[0][0]*A[1][1] - A[0][1]*A[1][0];
    if (fabs(*det) == 0.0) return(0);

    B[0][0] =  A[1][1]/(*det);
    B[0][1] = -A[0][1]/(*det);
    B[1][0] = -A[1][0]/(*det);
    B[1][1] =  A[0][0]/(*det);
  }
  else
  {
    (*det) = A[0][0]*(A[1][1]*A[2][2] - A[1][2]*A[2][1]) -
             A[0][1]*(A[1][0]*A[2][2] - A[1][2]*A[2][0]) +
             A[0][2]*(A[1][0]*A[2][1] - A[1][1]*A[2][0]);
    if (fabs(*det) == 0.0) return(0);

    B[0][0] =  (A[1][1]*A[2][2] - A[1][2]*A[2][1])/(*det);
    B[0][1] = -(A[0][1]*A[2][2] - A[0][2]*A[2][1])/(*det);
    B[0][2] =  (A[0][1]*A[1][2] - A[0][2]*A[1][1])/(*det);
    B[1][0] = -(A[1][0]*A[2][2] - A[1][2]*A[2][0])/(*det);
    B[1][1] =  (A[0][0]*A[2][2] - A[0][2]*A[2][0])/(*det);
    B[1][2] = -(A[0][0]*A[1][2] - A[0][2]*A[1][0])/(*det);
    B[2][0] =  (A[1][0]*A[2][1] - A[1][1]*A[2][0])/(*det);
    B[2][1] = -(A[0][0]*A[2][1] - A[0][1]*A[2][0])/(*det);
    B[2][2] =  (A[0][0]*A[1][1] - A[0][1]*A[1][0])/(*det);
  }

  return(1);
}

// ============================== ShpVecLen3D ==============================

double cShape :: ShpVecLen3D(double *v)
{
  return(sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]));
}

// ============================= ShpNormVec3D ==============================

void cShape :: ShpNormVec3D(double *v)
{
  double norm = ShpVecLen3D(v);
  v[0] /= norm;
  v[1] /= norm;
  v[2] /= norm;
}

// ============================== ShpCrossProd =============================

void cShape :: ShpCrossProd(double *u, double *v, double *w)
{
  w[0] = u[1]*v[2] - u[2]*v[1];
  w[1] = u[2]*v[0] - u[0]*v[2];
  w[2] = u[0]*v[1] - u[1]*v[0];
}

// ============================== GetShapeEdge =============================

int cShape :: GetShapeEdge(int NumEdge, tEdge *Edge, int *corner,
                           eShpType *type, int *nnode, cNode **conn)
{
  for (int i = 0; i < NumEdge; i++)
  {
    int ni = 0;                   // First node
    int nf = Edge[i].nnode - 1;   // Last node

    // Direct order.

    if (corner[0] == ElmNode[Edge[i].node[ni]]->GetLabel( ) &&
        corner[1] == ElmNode[Edge[i].node[nf]]->GetLabel( ))
    {
      *type  = Edge[i].type;
      *nnode = Edge[i].nnode;
      for (int j = 0; j < Edge[i].nnode; j++)
        conn[j] = ElmNode[Edge[i].node[j]];
      return(1);
    }

    // Inverse order.

    if (corner[0] == ElmNode[Edge[i].node[nf]]->GetLabel( ) &&
        corner[1] == ElmNode[Edge[i].node[ni]]->GetLabel( ))
    {
      *type  = Edge[i].type;
      *nnode = Edge[i].nnode;
      for (int j = 0; j < Edge[i].nnode; j++)
        conn[j] = ElmNode[Edge[i].node[nf-j]];
      return(1);
    }
  }

  // Not found.

  return(0);
}

// ============================== GetShapeFace =============================

int cShape :: GetShapeFace(int NumFace, tFace *Face, int *corner,
                           eShpType *type, int *nnode, cNode **conn)
{
  for (int i = 0; i < NumFace; i++)
  {
    int found0 = 0;
    int found1 = 0;
    int found2 = 0;
    for (int j = 0; j < Face[i].ncorner; j++)
    {
      if (corner[0] == ElmNode[Face[i].corner[j]]->GetLabel( )) found0 = 1;
      if (corner[1] == ElmNode[Face[i].corner[j]]->GetLabel( )) found1 = 1;
      if (corner[2] == ElmNode[Face[i].corner[j]]->GetLabel( )) found2 = 1;
    }
    int found = found0 + found1 + found2;
    if (found == 3)
    {
      *type  = Face[i].type;
      *nnode = Face[i].nnode;
      for (int k = 0; k < Face[i].nnode; k++)
        conn[k] = ElmNode[Face[i].node[k]];
      return(1);
    }
  }

  // Not found.

  return(0);
}

// ======================================================= End of file =====
