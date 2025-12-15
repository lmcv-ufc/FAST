// ------------------------------------------------------------------------
// shpshell.cpp - This file contains methods of Shell shape class.
// -------------------------------------------------------------------------
// Created:      18-Sep-2020     Elias Saraiva Barroso
//
// Modified:     17-May-2024     Elias Saraiva Barroso
//               Nodal vectors (axes and iaxes) for nonlinear analysis.
//
// Modified:     24-May-2024     Evandro Parente Junior
//               Method UpdNodalAxes for nonlinear analysis (TL).
// -------------------------------------------------------------------------

#include <math.h>
#include <stdio.h>
#include <iostream>

using namespace std;

#include "shpshell.h"
#include "vec.h"
#include "mat.h"
#include "gblvar.h"

#ifdef _OMP_
#include "omp.h"
#endif


// -------------------------------------------------------------------------
// Static variables:
//
double cShapeShell :: RefVec[3] = {0.0, 1.0, 0.0}; // {ey}

// -------------------------------------------------------------------------
// Local variables:
//
static double       _J[3][3];     // Jabobian matrix
static double    _InvJ[3][3];     // Inverse of Jacobian matrix
static double     _Mrst[500];     // Mapping functions
static sNatCoord _dMrst[500];     // Mapping func derivatives/parametric coords
static sNatCoord _dNrst[500];     // Shape func derivatives/parametric coords

#pragma omp threadprivate(_J,_InvJ,_Mrst,_dNrst,_dMrst)

// -------------------------------------------------------------------------
// Public methods:
//

// ============================== cShapeShell ===============================

cShapeShell :: cShapeShell(eShpType type)
{
  Type = type;

  // Create base shape.
  switch (type)
  {
    case SHAPE_T3_SHELL:
      baseShp = CreateShape(SHAPE_T3_SURF);
    break;

    case SHAPE_T6_SHELL:
      baseShp = CreateShape(SHAPE_T6_SURF);
    break;

    case SHAPE_Q8_SHELL:
      baseShp = CreateShape(SHAPE_Q8_SURF);
    break;

    case SHAPE_Q9_SHELL:
      baseShp = CreateShape(SHAPE_Q9_SURF);
    break;

    default:
    break;
  }

  if (baseShp->GetTopologyType( ) == TRIANGULAR_TOPOLOGY)
    TopType = WEDGE_TOPOLOGY; 
  else
    TopType = HEXAHEDRAL_TOPOLOGY;
   
  NumNode  = baseShp->GetNumNode( );
}

// ============================= ~cShapeShell ===============================

cShapeShell :: ~cShapeShell(void)
{
  delete baseShp;
  delete []Thk;
  delete []Axes;
  delete []InitAxes;
  delete []pa;
  delete []pb;
  delete []paa;
  delete []pab;
  delete []pbb;
}

// ================================= Read ==================================

void cShapeShell :: Read(void)
{
  baseShp->Read( );
  NumNode  = baseShp->GetNumNode( );
  InitAxes = 0;

  // Store node references from base shape.
  ElmNode = new cNode*[NumNode];
  for(int i = 0; i < NumNode; ++i)
    ElmNode[i] = baseShp->GetNode(i);

  // Store nodal thickness.
  Thk = new double [NumNode];
  for(int i = 0; i < NumNode; ++i)
    Thk[i] = GetNode(i)->GetThk( );
 
  // Alloc memory for nodal axes.
  Axes = new sNodeAxes[NumNode];
  for(int i = 0; i < NumNode; ++i)
  {
    Axes[i].V1.Resize(3);
    Axes[i].V2.Resize(3);
    Axes[i].V3.Resize(3);
  }
  
  if (cNode :: GetNumNormal( ) > 0) // Copy normal vector from node data.
    for(int i = 0; i < NumNode; ++i)
      Axes[i].V3 = GetNode(i)->GetNormal( );
  else
    CompNodalNormal( ); // Compute using tVector().
}

// ============================== LoadInitAxes =============================

void cShapeShell :: LoadInitAxes( )
{
  // Alloc memory for current nodal axes.
  InitAxes = new sNodeAxes[NumNode];
  for(int i = 0; i < NumNode; ++i)
  {
    InitAxes[i].V1 = Axes[i].V1;
    InitAxes[i].V2 = Axes[i].V2;
    InitAxes[i].V3 = Axes[i].V3;
  }

  pb  = new cVector[NumNode];
  pa  = new cVector[NumNode];
  paa = new cVector[NumNode];
  pab = new cVector[NumNode];
  pbb = new cVector[NumNode];

  for(int i = 0; i < NumNode; ++i)
  {
    pa[i].Resize(3);
    pb[i].Resize(3);
    paa[i].Resize(3);
    pab[i].Resize(3);
    pbb[i].Resize(3);
    pa[i] = (-Thk[i]/2.0)*InitAxes[i].V2; 
    pb[i] = ( Thk[i]/2.0)*InitAxes[i].V1;
    
    paa[i].Zero( );
    pab[i].Zero( );
    pbb[i].Zero( );
  }
}

// =========================== ReadReferenceVector =========================

void cShapeShell :: ReadReferenceVector(void)
{
  // Read the orientation vector.

  double a1,a2,a3;
  if (!(in >> a1) || !(in >> a2) || !(in >> a3))
  {
    cout << "Error in the input of section orientation vector\n";
    exit(0);
  }
  RefVec[0] = a1;
  RefVec[1] = a2;
  RefVec[2] = a3;
}

// ============================== CompNodalNormal ==========================

void cShapeShell :: CompNodalNormal( ) 
{
  // Get Node coord.
  sNodeCoord *coord  = new sNodeCoord[NumNode];
  sNatCoord  *ncoord = new sNatCoord[NumNode];
  baseShp->NodalCoord(coord);
  baseShp->NodeNatCoord(ncoord);

  for (int i = 0; i < NumNode; i++) 
  {
    Axes[i].V3.Resize(3);
    baseShp->tVector(ncoord[i],coord,Axes[i].V3);
    //Axes[i].V3 *= -1.0;
  }

  // Release memory.

  delete []coord;
  delete []ncoord;
}

// ============================= CompNodalTangents =========================

void cShapeShell :: CompNodalTangents(void)
{
  // Compute each node directors V1 and V2.

  double ex[] = {1.0, 0.0, 0.0};
  cVector e1(3, ex);
  cVector refvec(3, RefVec);
  refvec.Normalize( );
  for (int i = 0; i < NumNode; i++) 
  {
    // Normalize V3 and compute V1 = rv x V3.

    Axes[i].V3.Normalize( );
    CrossProd(refvec, Axes[i].V3, Axes[i].V1);

    // Check if V1 is a valid vector and compute V2.

    if (Axes[i].V1.Length( ) > 1.0e-6)
    {
      CrossProd(Axes[i].V3, Axes[i].V1, Axes[i].V2);  // V2 = V3 x V1
    } 
    else // Compute V2 and V1
    {
      CrossProd(Axes[i].V3, e1, Axes[i].V2);          // V2 = V3 x e1 
      CrossProd(Axes[i].V2, Axes[i].V3, Axes[i].V1);  // V1 = V2 x V3 
    } 

    // Normalize V1 and V2.
      
    Axes[i].V1.Normalize( );
    Axes[i].V2.Normalize( );
  }
}

// ============================== NodalCoord ===============================

void cShapeShell :: NodalCoord(sNodeCoord *coord)
{
  for (int i = 0; i < NumNode; i++) 
  {
    coord[i]       = GetNode(i)->GetCoord( );
    coord[i].thk   = &(Thk[i]); 
    coord[i].axes  = &(Axes[i]); 
    coord[i].iaxes = &(InitAxes[i]); 
    coord[i].pa    = &(pa[i]);
    coord[i].pb    = &(pb[i]);
    coord[i].paa   = &(paa[i]);
    coord[i].pab   = &(pab[i]);
    coord[i].pbb   = &(pbb[i]);
  }
}

// ============================== NodeNatCoord =============================

void cShapeShell :: NodeNatCoord(sNatCoord *c)
{
  baseShp->NodeNatCoord(c);
}

// ================================= ShpFunc ===============================

void cShapeShell :: ShpFunc(sNatCoord p, double *N)
{
  baseShp->ShpFunc(p, N);  // 2D shape functions
  N[NumNode] = p.t;        // Zeta 
}

// =============================== DrvShpRST ===============================

void cShapeShell :: DrvShpRST(sNatCoord p, sNatCoord *dN)
{
  baseShp->DrvShpRST(p, dN);  // Derivatives of 2D shape functions
  dN[NumNode].r = 0.0;        // dzeta/dr
  dN[NumNode].s = 0.0;        // dzeta/ds
  dN[NumNode].t = 1.0;        // dzeta/dt
}

// ============================== DrvShpXYZ ================================

void cShapeShell :: DrvShpXYZ(sNatCoord p, sNodeCoord *coord,
                             double *detJ, sNodeCoord *dNxyz)
{
  int i;
  double     thk2;
  sNodeAxes *axes;

  // Compute dNrs.

  DrvShpRST(p, _dNrst);

  // Compute [J].

  ShpMatZero(3, _J);
  DrvMapRST(p, _dMrst);
  MapFunc(p, _Mrst);
  double zeta = p.t;
  for (i = 0; i < NumNode; i++)
  {
    // Get current node thickness and normal vector.
    
    thk2 = *(coord[i].thk)/2.0;
    axes = coord[i].iaxes;        // Evandro (24-May-2024)

    // First column. (x,rst)

   _J[0][0] += _dMrst[i].r*(coord[i].x + zeta*thk2*axes->V3[0]);
   _J[1][0] += _dMrst[i].s*(coord[i].x + zeta*thk2*axes->V3[0]);
   _J[2][0] += _Mrst[i]*thk2*axes->V3[0];

    // Second column. (y,rst)

   _J[0][1] += _dMrst[i].r*(coord[i].y + zeta*thk2*axes->V3[1]);
   _J[1][1] += _dMrst[i].s*(coord[i].y + zeta*thk2*axes->V3[1]);
   _J[2][1] += _Mrst[i]*thk2*axes->V3[1];

    // Third column. (z,rst)

   _J[0][2] += _dMrst[i].r*(coord[i].z + zeta*thk2*axes->V3[2]);
   _J[1][2] += _dMrst[i].s*(coord[i].z + zeta*thk2*axes->V3[2]);
   _J[2][2] += _Mrst[i]*thk2*axes->V3[2];
  }

  // Numerically invert the Jacobian matrix.

  ShpMatInv(3, _J, _InvJ, detJ);

  // Compute Ni,xyz.

  for (i = 0; i < NumNode; i++)
  {
    dNxyz[i].x = _InvJ[0][0]*_dNrst[i].r + _InvJ[0][1]*_dNrst[i].s;
    dNxyz[i].y = _InvJ[1][0]*_dNrst[i].r + _InvJ[1][1]*_dNrst[i].s;
    dNxyz[i].z = _InvJ[2][0]*_dNrst[i].r + _InvJ[2][1]*_dNrst[i].s;
  }

  // Store zeta,xyz.

  dNxyz[NumNode].x = _InvJ[0][2];
  dNxyz[NumNode].y = _InvJ[1][2];
  dNxyz[NumNode].z = _InvJ[2][2];
}

// ================================ NMatrix ================================

void cShapeShell :: NMatrix(int, double *shpfunc, cMatrix &N)
{
  double zeta = shpfunc[NumNode];
  for (int i = 0, col = 0; i < NumNode; i++, col+=5) 
  {
    N[0][col+0] = shpfunc[i];
    N[1][col+0] = 0.0;
    N[2][col+0] = 0.0;
    N[3][col+0] = 0.0;
    N[4][col+0] = 0.0;

    N[0][col+1] = 0.0;
    N[1][col+1] = shpfunc[i];
    N[2][col+1] = 0.0;
    N[3][col+1] = 0.0;
    N[4][col+1] = 0.0;

    N[0][col+2] = 0.0;
    N[1][col+2] = 0.0;
    N[2][col+2] = shpfunc[i];
    N[3][col+2] = 0.0;
    N[4][col+2] = 0.0;

    N[0][col+3] = -0.5*zeta*shpfunc[i]*Thk[i]*Axes[i].V2[0];
    N[1][col+3] = -0.5*zeta*shpfunc[i]*Thk[i]*Axes[i].V2[1];
    N[2][col+3] = -0.5*zeta*shpfunc[i]*Thk[i]*Axes[i].V2[2];
    N[3][col+3] = 0.0;
    N[4][col+3] = 0.0;

    N[0][col+4] = 0.5*zeta*shpfunc[i]*Thk[i]*Axes[i].V1[0];
    N[1][col+4] = 0.5*zeta*shpfunc[i]*Thk[i]*Axes[i].V1[1];
    N[2][col+4] = 0.5*zeta*shpfunc[i]*Thk[i]*Axes[i].V1[2];
    N[3][col+4] = 0.0;
    N[4][col+4] = 0.0;
  }
}

// ============================== LocalSys =================================

void cShapeShell :: LocalSys(sNatCoord p, sNodeCoord *coord, cMatrix &R)
{
  int i;

  // Compute Shape functions of middle surface.

  baseShp->ShpFunc(p, _Mrst);

  // Evaluate vectos at curret position.

  cVector e1(3),e2(3), e3(3);
  e1.Zero( );
  e2.Zero( );
  e3.Zero( );
  for (i = 0; i < NumNode; i++)
  {
    e1 += _Mrst[i]*coord[i].iaxes->V1;   // Evandro (24-May-2024)
    e2 += _Mrst[i]*coord[i].iaxes->V2;   // Evandro (24-May-2024)
    e3 += _Mrst[i]*coord[i].iaxes->V3;   // Evandro (24-May-2024)
  }
  e1.Normalize( );
  e2.Normalize( );
  e3.Normalize( );

  // Assembly the rotation matrix (vectors in columns).

  for (i = 0; i < 3; i++)
  {
    R[i][0] = e1[i];
    R[i][1] = e2[i];
    R[i][2] = e3[i];
  }
}

// ================================ GetEdge ================================

int cShapeShell :: GetEdge(int *corner, eShpType *type, int *nnode,
                           cNode **conn)
{
  return(baseShp->GetEdge(corner,type,nnode,conn));
}

// ================================ GetFace ================================

int cShapeShell :: GetFace(int *corner, eShpType *type, int *nnode,
                           cNode **conn)
{
  return(baseShp->GetFace(corner,type,nnode,conn));
}

// ============================== tVector ==================================

void cShapeShell :: tVector(sNatCoord p, sNodeCoord *coord, cVector &w)
{
  return baseShp->tVector(p,coord,w);
}

// ============================== SetNodalAxes =============================

void cShapeShell :: SetNodalAxes(sNodeAxes *axes)
{
  // Set normal vector.
  for (int i = 0; i < NumNode; i++)
    Axes[i].V3 = axes[i].V3;
  
  // Evaluate nodal tangent directors.
  CompNodalTangents( );
}

// ============================== UpdNodalAxes =============================

void cShapeShell :: UpdNodalAxes(sNodeAxes *axes)
{
  // Update nodal vectors.

  for (int i = 0; i < NumNode; i++)
  {
//    cout << "Before update:\n";	  
//    Axes[i].V1.Print();
//    Axes[i].V2.Print();
//    Axes[i].V3.Print();

    Axes[i].V1 = axes[i].V1;
    Axes[i].V2 = axes[i].V2;
    Axes[i].V3 = axes[i].V3;

//    cout << "After update:\n";	  
//    Axes[i].V1.Print();
//    Axes[i].V2.Print();
//    Axes[i].V3.Print();
  }
}

// -------------------------------------------------------------------------
// Classe shape shell IGA:
//

// ============================ cShapeShellIGA ==============================

cShapeShellIGA :: cShapeShellIGA(eShpType type)
{
  Type = type;

  // Create base shape.
  switch (type)
  {
    case SHAPE_BEZTRI_SHELL:
      baseShp = CreateShape(SHAPE_BEZTRI_SURF);
      TopType = WEDGE_TOPOLOGY;
    break;

    case SHAPE_BEZ_SHELL:
      baseShp = CreateShape(SHAPE_BEZ_SURF);
      TopType = HEXAHEDRAL_TOPOLOGY;
    break;

    case SHAPE_BSP_SHELL:
      baseShp = CreateShape(SHAPE_BSP_SURF);
      TopType = HEXAHEDRAL_TOPOLOGY;
    break;

    default:
    break;
  }

  NumNode = baseShp->GetNumNode( );
}

// =========================== ~cShapeShellIGA ==============================

cShapeShellIGA :: ~cShapeShellIGA(void)
{
}

// ============================== SetNodalAxes =============================

void cShapeShellIGA :: SetNodalAxes(sNodeAxes *axes)
{
  // TODO: Calcular normais em NURBS dentro do FAST!	
  if (Type == SHAPE_BSP_SHELL) 
  {
    this->cShapeShell::SetNodalAxes(axes);
    return;
  }	  

  // Auxiliary variables. 
  double    *N      = new double [GetNumMap( )];    // Mapping functions.
  sNatCoord *ncoord = new sNatCoord[NumNode];       // Natural coordinates.

  // Evaluate mapping function matrix at a set of uniformly distributed points
  // in parametric space.
  NodeNatCoord(ncoord);
  cMatrix B(NumNode,NumNode); 
  for (int i = 0; i < NumNode; i++) 
  {
    MapFunc(ncoord[i],N);
    for (int j = 0; j < NumNode; j++) B[i][j] = N[j]; 
  }

  // LU decomposition.
  B.DecompLU( );

  // Evaluate control point normal vectors.
  cVector X(NumNode), Xcp(NumNode);
  for(int comp = 0; comp < 3; ++comp)
  {
    for (int i = 0; i < NumNode; i++) X[i] = axes[i].V3[comp];  
    B.SolveLU(X,Xcp);
    for (int i = 0; i < NumNode; i++) Axes[i].V3[comp] = Xcp[i];
  }

  // Normalize normal vectors.
  for (int i = 0; i < NumNode; i++) Axes[i].V3.Normalize( );

  // Evaluate nodal tangent directors.
  CompNodalTangents( );

  // Release memory.

  delete []N;
  delete []ncoord;
}

// ====================================================== End of File ======
