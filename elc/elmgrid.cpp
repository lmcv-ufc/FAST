// -------------------------------------------------------------------------
// elmgrid.cpp - implementation of grid class.
// -------------------------------------------------------------------------
// Created:      15-Jun-2012     Elias Saraiva Barroso
// -------------------------------------------------------------------------

#include <math.h>
#include <iostream>
#include <cstdlib>

#include "elmgrid.h"
#include "node.h"
#include "material.h"
#include "section.h"
#include "secbar.h"
#include "shpline.h"
#include "intpoint.h"
#include "ctrl.h"
#include "load.h"
#include "vec.h"
#include "mat.h"
#include "gblvar.h"
#include "gbldef.h"

// -------------------------------------------------------------------------
// cGrid class:
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Public methods:
//

// ============================= ValidSection ==============================

bool cGrid :: ValidSection(eSecType type)
{
  switch (type)
  {
    case SEC_BAR_GEN:
    case SEC_BAR_B3DCOUPLED:
    case SEC_BAR_HOM_RECT:
    case SEC_BAR_HOM_CIRCLE:
    case SEC_BAR_HOM_TUBE:
      return(true);
    break;

    default:
    break;
  }

  return(false);
}

// =============================== cGrid ================================

cGrid :: cGrid(int id, cSection *sec) : cElement(id)
{
  Type = GRID;
  Section = (cSecBar *)sec;
  Shape = new cShapeBar( );
}

// ============================== ~cGrid ================================

cGrid :: ~cGrid(void)
{
}

// ================================= Read ==================================

void cGrid :: Read(void)
{
  // Read beam orientations.
  
  int ort;
  if (!(in >> ort))
  {
    cout << "Error in the input of element " << Label << "!\n";
    exit(0);
  }
  
  if ((ort < 1) || (ort > NumBeamOrt))
  {
    cout << "Error in the input of element " << Label << "!\n";
    exit(0);
  }
  BeamOrtIdx = ort - 1;
      
  // Read element incidence.

  Shape->Read( );
}

// ============================== GetActDir ================================

void cGrid :: GetActDir(int *dir)
{
  dir[0] = 0;   // u
  dir[1] = 0;   // v
  dir[2] = 1;   // w
  dir[3] = 1;   // rx
  dir[4] = 1;   // ry
  dir[5] = 0;   // rz
}

// ============================= GetStrLabels ==============================

void cGrid :: GetStrLabels(int *label)
{
  label[0] = FORCE_Z;
  label[1] = MOMENT_X;
  label[2] = MOMENT_Y;
}

// =============================== IntForce ================================

int cGrid :: IntForce(cVector &g)
{
  // Get the nodal coordinates.
  
  sNodeCoord coord[2];
  Shape->NodalCoord(coord);
    
  // Compute the element length.
  
  double L = CalcLength(coord);
  
  // Assembly the transformation matrix.
  
  cMatrix T(6, 6);
  GetTrnMat(T);
  
  // Get the nodal displacements.
  
  cVector u(6);
  NodalDispl(u);
  
  // Compute the local displacements => {ul} = [T]{u}.
  
  cVector ul(6);
  ul = T*u;
  
  // Assembly the local stiffness matrix.
  
  cMatrix Ke(6, 6);
  LocStiffMat(L, Ke);
  
  // Compute the local forces => {gl} = [Ke]{ul}.
  
  cVector gl(6);
  gl = Ke*ul;
  
  // Transform to the global system => {g} = [T]t{gl}
  
  g = t(T)*gl;
  
  return (1);
}

// =============================== StiffMat ================================

void cGrid :: StiffMat(cMatrix &K)
{
  // Get te nodal coordinates.
  
  sNodeCoord coord[2];
  Shape->NodalCoord(coord);
  
  // Compute the element length.
  
  double L = CalcLength(coord);
  
  // Assembly the transformation matrix.
  
  cMatrix T(6, 6);
  GetTrnMat(T);
  
  // Assembly the local stiffness matrix.
  
  cMatrix Ke(6, 6);
  LocStiffMat(L, Ke);
  
  // Compute the stiffness matrix => [K] = [T]t[Ke][T].
    
  K.Zero( );
  MatTripBtCB(T, Ke, 1.0, K);
}

// ============================== NodalStress ==============================

void cGrid :: NodalStress(cMatrix &S)
{
  // Get the nodal coordinates.
  
  sNodeCoord coord[2];
  Shape->NodalCoord(coord);
  
  // Assembly the transformation matrix.
  
  cMatrix T(6, 6);
  GetTrnMat(T);
  
  // Compute the internal force vector.
  
  cVector g(6);
  IntForce(g);
  
  // Transform to the local system => {gl} = [T]{g}.
  
  cVector gl(6);
  gl = T*g;
  
  // Add the fixed-end forces (- equivforces).
  
  cVector fl(6);
  EqvForces(fl);
  gl -= fl;
  
  // Return nodal forces.

  if (OutputConv == DIRECT_STIFFNESS)
  {
    S[0][0] =  gl[0];  // Initial node
    S[0][1] =  gl[1];
    S[0][2] =  gl[2];
    S[1][0] =  gl[3];  // Final node
    S[1][1] =  gl[4];
    S[1][2] =  gl[5];
  }
  else
  {
    S[0][0] =  gl[0];  // Initial node
    S[0][1] = -gl[1];
    S[0][2] =  gl[2];
    S[1][0] = -gl[3];  // Final node
    S[1][1] =  gl[4];
    S[1][2] = -gl[5];
  }
}

/*
// ============================== CalcRotMat ==============================

void cGrid :: CalcRotMat(cMatrix &R)
{
  double e0x[3];  // vector of the direct cosine of the bar
  double e0z[3];  // vector of the direct cosine of the z-axis
  double e0y[3];  // vector of the direct cosine of the y-axis
  
  // Get the nodal coordinates.
  
  sNodeCoord coord[2];
  Shape->NodalCoord(coord);
  
  // Computer the element length.
  
  double L = CalcLength(coord);
  
  // Computer the direct cosines.
  
  double dx = coord[1].x - coord[0].x;
  double dy = coord[1].y - coord[0].y;
  double dz = coord[1].z - coord[0].z;
  
  e0x[0] = dx/L;  // lx - direct cosine between bar and global x-axis
  e0x[1] = dy/L;  // mx - direct cosine between bar and global y-axis
  e0x[2] = dz/L;  // nx - direct cosine between bar and global z-axis
  
  VecCrossProd(e0x, VecBeamOrt[BeamOrtIdx].e0y, e0z);
  VecCrossProd(e0z, e0x, e0y);
  
  // Assembly the rotation matrix.
  
  R.Zero( );
  R[0][0] = e0x[0];  R[0][1] = e0x[1];  R[0][2] = e0x[2];
  R[1][0] = e0y[0];  R[1][1] = e0y[1];  R[1][2] = e0y[2];
  R[2][0] = e0z[0];  R[2][1] = e0z[1];  R[2][2] = e0z[2];
}
*/

// -------------------------------------------------------------------------
// Protected methods:
//

// ================================ CalcLength =============================
//
// This method computes the length of a given bar.
//
//   coord  - vector of nodal coordinates                              (in)
//   L      - element length                                          (out)
//
double cGrid :: CalcLength(sNodeCoord *coord)
{
  double dx = coord[1].x - coord[0].x;
  double dy = coord[1].y - coord[0].y;
  double dz = coord[1].z - coord[0].z;
  double L = sqrt(dx*dx + dy*dy + dz*dz);
  
  return (L);
}

// =============================== GetTrnMat ===============================
//
// This method computes the transformation matrix between the local and the
// global coordinate system.
//
//   T  - transformation matrix (6x6)                               (out)
//
void cGrid :: GetTrnMat(cMatrix &T)
{
  double s;  // vector of the direct cosine of the bar
  double c;  // vector of the direct cosine of the z-axis
  
  // Get nodal coordinates.
  
  sNodeCoord coord[2];
  Shape->NodalCoord(coord);
  
  // Computer the direct cosines.

  double dx = coord[1].x - coord[0].x;
  double dy = coord[1].y - coord[0].y;
  double dz = coord[1].z - coord[0].z;
  double l  = sqrt(dx*dx + dy*dy + dz*dz);  // element length
  
  c = dx/l;  // lx - direct cosine between bar and global x-axis
  s = dy/l;  // mx - direct cosine between bar and global y-axis
    
  T.Zero( );
  T[0][0] = T[3][3] = 1;  
  T[0][1] = T[3][4] = 0;  
  T[0][2] = T[3][5] = 0;  
    
  T[1][0] = T[4][3] = 0;  
  T[1][1] = T[4][4] = c;  
  T[1][2] = T[4][5] = s;  
    
  T[2][0] = T[5][3] = 0;  
  T[2][1] = T[5][4] = -s;  
  T[2][2] = T[5][5] = c;  

/*
  cMatrix R(3, 3);
  CalcRotMat(R);
  
  T.Zero( );
  T[0][0] = T[3][3] = T[6][6] = T[ 9][ 9] = R[0][0];  // lx
  T[0][1] = T[3][4] = T[6][7] = T[ 9][10] = R[0][1];  // mx
  T[0][2] = T[3][5] = T[6][8] = T[ 9][11] = R[0][2];  // nx
  
  T[1][0] = T[4][3] = T[7][6] = T[10][ 9] = R[1][0];  // ly
  T[1][1] = T[4][4] = T[7][7] = T[10][10] = R[1][1];  // my
  T[1][2] = T[4][5] = T[7][8] = T[10][11] = R[1][2];  // ny
  
  T[2][0] = T[5][3] = T[8][6] = T[11][ 9] = R[2][0];  // lz
  T[2][1] = T[5][4] = T[8][7] = T[11][10] = R[2][1];  // mz
  T[2][2] = T[5][5] = T[8][8] = T[11][11] = R[2][2];  // nz
*/
}

// ============================== LocSitffMat ==============================
//
// This method computes the element sitffness matrix the local system.
//
//   L  - element length                                               (in)
//   Kl - elastic matrix (6x6)                                        (out)
//
void cGrid :: LocStiffMat(double L, cMatrix &Kl)
{
  Kl.Zero( );
  
  // Classical stiffness matrix (diagonal section matrix).

  if (Section->GetType( ) != SEC_BAR_B3DCOUPLED)
  {
    cMatMec *mat = Section->GetMaterial( )->GetMecProp( );
    double matpar[2];
    mat->GetParam(matpar);
    double E  = matpar[0];
    double nu = matpar[1];
    double G = 0.5*E/(1.0 + nu);

    double EIy = E*Section->GetIy( );
    double GJ  = G*Section->GetJt( );

    double l2 = L*L;
    double l3 = l2*L;
        
    Kl.Zero( );
      
    Kl[0][0] =  12*EIy/l3;
    Kl[0][2] = Kl[ 2][ 0] = -6*EIy/l2;
    Kl[0][3] = Kl[ 3][ 0] = -12*EIy/l3;
    Kl[0][5] = Kl[5][ 0] =  -6*EIy/l2;
    Kl[1][1] = Kl[ 4][ 4] =  GJ/L;
    Kl[1][4] = Kl[ 4][ 1] = -GJ/L;
    Kl[2][2] =  4*EIy/L;
    Kl[2][3] = Kl[ 3][ 2] = 6*EIy/l2;
    Kl[2][5] = Kl[5][ 2] = 2*EIy/L;
    Kl[3][3] = 12*EIy/l3;
    Kl[3][5] = Kl[5][3] = 6*EIy/l2;       
    Kl[5][5] = 4*EIy/L;
  }
  else // Coupled section constitutive matrix
  {
    cout << "Error: Grid with coupled section not implemented yet!";
    exit(0); 
  }
}

// =============================== EqvForces ===============================
//
// This method returns the equivalent nodal forces due to the loads
// applied at the element.
//
//   fl - equivalent nodal forces in the local system                 (out)
//
void cGrid :: EqvForces(cVector &fl)
{
  // Get the equivalent nodal forces in the global system.
  
  cVector fg(6);
  fg.Zero( );
  double t = cControl::GetTotTime( );
  cLoad :: EqvForces(this, t, fg);
  
  // Assembly the transformation matrix.
  
  cMatrix T(6, 6);
  GetTrnMat(T);
  
  // Rotate fix end forces to local system => {fl} = [T]{fg}.
  
  fl = T*fg;
}

// Código antigo!! FAST 3.0!!

/*
// -------------------------------------------------------------------------
// cGrid class:
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Public methods:
//

// =============================== cGrid ================================

cGrid :: cGrid(int id, eSecType sec) : cElmBar(id)
{
  Type = GRID;
  SecType = sec;
}

// =============================== ~cGrid ===============================

cGrid :: ~cGrid(void)
{
}

// =============================== GetActDir ===============================

void cGrid :: GetActDir(int *dir)
{
  dir[0] = 0;   // u
  dir[1] = 0;   // v
  dir[2] = 1;   // w
  dir[3] = 1;   // rx
  dir[4] = 1;   // ry
  dir[5] = 0;   // rz
}

// ============================== GetStrLabels ==============================

void cGrid :: GetStrLabels(int *label)
{
  label[0] = FORCE_Z;
  label[1] = MOMENT_X;
  label[2] = MOMENT_Y;
}

// ================================ IntForce ================================

int cGrid :: IntForce(cVector &g)
{
  // Get the node coordinates
  sNodeCoord coord[2];
  Shape->NodalCoord(coord);
    
  // Compute the element length.
    
  double l0;
  CalcLength(coord, &l0);
    
  // Assembly the transformation matrix.
    
  cMatrix T(6, 6);
  GetTrnMat(T);
    
  // Get the nodal displacements.
    
  cVector u(6);
  NodalDispl(u);
    
  // Compute the local displacements.
    
  cVector ul(6);
  ul = T*u;
    
  // Assembly the element matrix.
    
  cMatrix Ke(6,6);
  GetElastMat(l0, Ke);
    
  // Compute the local forces.
    
  cVector gl(6);
  gl = Ke*ul;
    
  // Trasform to the global system => {g} = [T]t {gl}.
    
  g = t(T)*gl;
    
  return (1);
}

// =============================== StiffMat ================================

void cGrid :: StiffMat(cMatrix &K)
{
  // Get the nodal coordinates.

  sNodeCoord coord[2];
  Shape->NodalCoord(coord);
    
  // Compute the element length.
    
  double l0;
  CalcLength(coord, &l0);
        
  // Assembly the transformation matrix.
    
  cMatrix T(6, 6);
  GetTrnMat(T);
  
  // Assembly the element matrix.
    
  cMatrix Ke(6,6);
  GetElastMat(l0, Ke);
  
  // Compute the stiffness matrix => [K] = [T]t[Ke][T].
    
  K.Zero( );
  MatTripBtCB(T, Ke, 1.0, K);
}

// ============================== NodalStress ==============================

void cGrid :: NodalStress(cMatrix &S)
{
  // Get the nodal coordinates.
  sNodeCoord coord[2];
  Shape->NodalCoord(coord);
    
  // Compute the element length.
    
  double l0;
  CalcLength(coord, &l0);
    
  // Assembly the transformation matrix.
    
  cMatrix T(6,6);
  GetTrnMat(T);
    
  // Get the nodal displacement.
    
  cVector u(6);
  NodalDispl(u);
    
  // Compute the local displacements.
    
  cVector ul(6);
  ul = T*u;
    
  // Assembly the element matrix.
    
  cMatrix Ke(6,6);
  GetElastMat(l0, Ke);
    
  // Compute the local forces.
    
  cVector gl(6);
  gl = Ke*ul;
  
  // Add the fixed-end forces (- equivforces).
    
  cVector f(6);
  EqvForces(f);
  gl -= f;
      
  // Store in the matrix format (isostatic convention).
    
  S[0][0] =  gl[0];  // initial nodal
  S[0][1] = -gl[1];
  S[0][2] =  gl[2];
  
  S[1][0] = -gl[3];  // final nodal
  S[1][1] =  gl[4];
  S[1][2] = -gl[5];
}

// // --------------------------------------------------------------------------
// Protected methods:
//

// ================================ CalcLength ==============================
//
// This method computes the length of a given bar.
//
//   coord  - vector of nodal coordinates                               (in)
//   l      - element length                                           (out)
//
void cGrid :: CalcLength(sNodeCoord *coord, double *l)
{
  double dx = coord[1].x - coord[0].x;
  double dy = coord[1].y - coord[0].y;
  double dz = coord[1].z - coord[0].z;
    
  *l = sqrt(dx*dx + dy*dy + dz*dz);
}

// =============================== GetTrnMat ===============================
//
// This method computes the transformation matrix between the natural
// and the global coordinate system.
//
//   T  - transformation matrix (6x6)                                (out)
//
void cGrid :: GetTrnMat(cMatrix &T)
{
  double s;  // vector of the direct cosine of the bar
  double c;  // vector of the direct cosine of the z-axis
  
  // Get nodal coordinates.
  
  sNodeCoord coord[2];
  Shape->NodalCoord(coord);
  
  // Computer the direct cosines.

  double dx = coord[1].x - coord[0].x;
  double dy = coord[1].y - coord[0].y;
  double dz = coord[1].z - coord[0].z;
  double l  = sqrt(dx*dx + dy*dy + dz*dz);  // element length
  
  c = dx/l;  // lx - direct cosine between bar and global x-axis
  s = dy/l;  // mx - direct cosine between bar and global y-axis
    
  T.Zero( );
  T[0][0] = T[3][3] = 1;  
  T[0][1] = T[3][4] = 0;  
  T[0][2] = T[3][5] = 0;  
    
  T[1][0] = T[4][3] = 0;  
  T[1][1] = T[4][4] = c;  
  T[1][2] = T[4][5] = s;  
    
  T[2][0] = T[5][3] = 0;  
  T[2][1] = T[5][4] = -s;  
  T[2][2] = T[5][5] = c;  
}

// ============================== GetElastMat ==============================
//
// This method computes the elastic constitutive matrix for the natural
// coordinate system.
//
//   l - element length                                                (in)
//   C - elastic matrix (6x6)                                       (out)
//
void cGrid :: GetElastMat(double l, cMatrix &C)
{
  double matpar[2];
  Material->GetParam(matpar);
  double E  = matpar[0];
  double nu = matpar[1];
  double G = 0.5*E/(1.0 + nu);
        
// double EA  = E*VecSec[Section->GetCrossSecIdx( )].Ax; nao precisa?
// double EIz = E*VecSec[Section->GetCrossSecIdx( )].Iz; nao precisa?
  double EIy = E*VecSec[Section->GetCrossSecIdx( )].Iy;
  double GJ  = G*VecSec[Section->GetCrossSecIdx( )].Ix;
  std::cout << "Ix = " << VecSec[Section->GetCrossSecIdx( )].Ix << "\n";
  std::cout << "Iy = " << VecSec[Section->GetCrossSecIdx( )].Iy << "\n";

  double l2 = l*l;
  double l3 = l2*l;
      
  C.Zero( );
    
  C[0][0] =  12*EIy/l3;
  C[0][2] = C[ 2][ 0] = -6*EIy/l2;
  C[0][3] = C[ 3][ 0] = -12*EIy/l3;
  C[0][5] = C[5][ 0] =  -6*EIy/l2;
  C[1][1] = C[ 4][ 4] =  GJ/l;
  C[1][4] = C[ 4][ 1] = -GJ/l;
  C[2][2] =  4*EIy/l;
  C[2][3] = C[ 3][ 2] = 6*EIy/l2;
  C[2][5] = C[5][ 2] = 2*EIy/l;
  C[3][3] = 12*EIy/l3;
  C[3][5] = C[5][3] = 6*EIy/l2;       
  C[5][5] = 4*EIy/l;
  
  
}

// ============================== EqvForces =============================
//
// This method returns the equivalent nodal forces due to the loads
// applied at the element.
//
//   f - equivalent nodal forces in the local system                (out)
//
void cGrid :: EqvForces(cVector &f)
{
  // Get the equivalent nodal forces in the global system.
  cVector fg(6);
  fg.Zero( );
  double t = cControl :: GetTotTime( );
  cLoad :: EqvForces(this, t, fg);
    
  // Rotate fix end forces to local system => {fl} = [T] {fg}
  
  cMatrix T(6,6);
  GetTrnMat(T);
    
  f = T*fg;
}
*/

// ======================================================= End of file =====
