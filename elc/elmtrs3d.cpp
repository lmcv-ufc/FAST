// -------------------------------------------------------------------------
// elmtrs3d.cpp - implementation of 3D truss elements.
// -------------------------------------------------------------------------
// Created:      06-Sep-2015     Evandro Parente Junior
//               Based on plane truss elements (elmpltrs.cpp).
//
// Modified:
// -------------------------------------------------------------------------

#include <iostream>
#include <math.h>
#include <stdlib.h>

using namespace std;

#include "elmtrs3d.h"
#include "node.h"
#include "material.h"
#include "cmodel.h"
#include "secbar.h"
#include "shpline.h"
#include "vec.h"
#include "gblvar.h"
#include "gbldef.h"

// -------------------------------------------------------------------------
// Class cTruss3D:
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Public methods:
//

// =============================== cTruss3D ================================

cTruss3D :: cTruss3D(int id) : cElement(id)
{
  Type = TRUSS3D;
  Shape = new cShapeBar( );
}

// ============================== ~cTruss3D ================================

cTruss3D :: ~cTruss3D(void)
{
  delete ConstMod;
}

// ================================= Read ==================================

void cTruss3D :: Read(void)
{
  int sec;
  if (!(in >> sec))
  {
    cout << "Error in the input of truss element " << Label << "!\n";
    cout << "Invalid section!\n";
    exit(0);
  }

  Section = (cSecBar *)cSection::GetSection(sec);
  if (!Section)
  {
    cout << "Error in the input of truss element " << Label << "!\n";
    cout << "Invalid section!\n";
    exit(0);
  }

  // Create the constitutive model.

  cMaterial *mat = Section->GetMaterial( );
  ConstMod = cConstModel::Create(this, mat);

  // Read element incidence.

  Shape->Read( );
}

// ============================== GetActDir ================================

void cTruss3D :: GetActDir(int *dir)
{
  dir[0] = 1;   // u
  dir[1] = 1;   // v
  dir[2] = 1;   // w
  dir[3] = 0;   // rx
  dir[4] = 0;   // ry
  dir[5] = 0;   // rz
}

// ============================= GetStrLabels ==============================

void cTruss3D :: GetStrLabels(int *label)
{
  label[0] = FORCE_X;
}

// =============================== IntForce ================================

int cTruss3D :: IntForce(cVector &g)
{
  // Get the nodal coordinates.

  sNodeCoord coord[2];
  Shape->NodalCoord(coord);

  // Get the nodal displacements.

  cVector u(6);
  NodalDispl(u);

  // Evaluate the strain-displacement matrix.

  double L;
  cVector B0(6);
  B0Matrix(coord, &L, B0);

  // Evaluate the normal force.

  cVector eps(1);
  cVector sig(1);
  eps[0] = B0*u;
  if (!ConstMod->Stress(eps, sig)) return(0);
  double N = Section->GetA( )*sig[0];

  // Evaluate {g} = [B]t*N*L;

  g = (N*L)*B0;

  return(1);
}

// =============================== StiffMat ================================

void cTruss3D :: StiffMat(cMatrix &K)
{
  // Get the nodal coordinates.

  sNodeCoord coord[2];
  Shape->NodalCoord(coord);

  // Evaluate the strain-displacement matrix.

  double L;
  cVector B0(6);
  B0Matrix(coord, &L, B0);

  // Evaluate the tangent elastic modulus.

  cMatrix Ct(1, 1);
  ConstMod->TangMat(Ct);

  // Evaluate [Kt] = A*L*[B]t*E*[B]

  double EAL = Ct[0][0]*Section->GetA( )*L;
  for (int i = 0; i < 6; i++)
    for (int j = 0; j < 6; j++) K[i][j] = EAL*B0[i]*B0[j];
}

// =============================== MassMat =================================

void cTruss3D :: MassMat(cMatrix &M)
{
  // Compute the bar length.

  sNodeCoord coord[2];
  Shape->NodalCoord(coord);
  double dx = coord[1].x - coord[0].x;
  double dy = coord[1].y - coord[0].y;
  double dz = coord[1].z - coord[0].z;
  double L  = sqrt(dx*dx + dy*dy + dz*dz);

  // Assembly the mass matrix.

  cMaterial *mat = Section->GetMaterial( );
  double rho = mat->GetDensity( );
  double A = Section->GetA( );
  double coeff = rho*A*L/6.0;
  M.Zero( );
  M[0][0] = M[1][1] = M[2][2] = M[3][3] = M[4][4] = M[5][5] = 2.0*coeff;
  M[3][0] = M[4][1] = M[5][2] = coeff;
  M[0][3] = M[1][4] = M[2][5] = coeff;
}

// =============================== GeomStiff ===============================

void cTruss3D :: GeomStiff(cMatrix &G)
{
  // Get the nodal coordinates.

  sNodeCoord coord[2];
  Shape->NodalCoord(coord);

  // Get the nodal displacements.

  cVector u(6);
  NodalDispl(u);

  // Evaluate the strain-displacement matrix.

  double L;
  cVector B0(6);
  B0Matrix(coord, &L, B0);

  // Evaluate the normal force.

  cVector eps(1);
  cVector sig(1);
  eps[0] = B0*u;
  ConstMod->Stress(eps, sig);
  double N = Section->GetA( )*sig[0];

  // Assembly the geometric stiffness matrix.

  G.Zero( );
  G[0][0] = G[1][1] = G[2][2] = G[3][3] = G[4][4] = G[5][5] = N/L;
  G[3][0] = G[4][1] = G[5][2] = -N/L;
  G[0][3] = G[1][4] = G[2][5] = -N/L;
}

// ============================== NodalStress ==============================

void cTruss3D :: NodalStress(cMatrix &S)
{
  // Get the nodal coordinates.

  sNodeCoord coord[2];
  Shape->NodalCoord(coord);

  // Get the nodal displacements.

  cVector u(6);
  NodalDispl(u);

  // Evaluate the strain-displacement matrix.

  double L;
  cVector B0(6);
  B0Matrix(coord, &L, B0);

  // Evaluate the normal force.

  cVector eps(1);
  cVector sig(1);
  eps[0] = B0*u;
  ConstMod->Stress(eps, sig);
  double N = Section->GetA( )*sig[0];

  // Store in the given matrix.

  if (OutputConv == DIRECT_STIFFNESS)
  {   
    S[0][0] = -N;
    S[1][0] =  N;
  }
  else
  {
    S[0][0] = S[1][0] = N;
  }
}

// ============================== UpdateState ==============================

void cTruss3D :: UpdateState(void)
{
  ConstMod->UpdateState( );
}

// -------------------------------------------------------------------------
// Protected methods:
//

// ================================ B0Matrix ===============================
//
// This method computes the element linear strain-displacement matrix.
//
//   coord - vector of nodal coordinates                               (in)
//   L     - element length                                           (out)
//   B0    - linear strain displacement matrix                        (out)
//
void cTruss3D :: B0Matrix(sNodeCoord *coord, double *L, cVector &B0)
{
  // Evaluate element length.

  double dx = coord[1].x - coord[0].x;
  double dy = coord[1].y - coord[0].y;
  double dz = coord[1].z - coord[0].z;
  double L2 = dx*dx + dy*dy + dz*dz;
  *L = sqrt(L2);

  // Evaluate [B0].

  B0[0] = -dx/L2;
  B0[1] = -dy/L2;
  B0[2] = -dz/L2;
  B0[3] = -B0[0];
  B0[4] = -B0[1];
  B0[5] = -B0[2];
}


// -------------------------------------------------------------------------
// Class cTruss3DTL:
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Public methods:
//

// ============================== cTruss3DTL ===============================

cTruss3DTL :: cTruss3DTL(int id) : cTruss3D(id)
{
  Type = TRUSS3DTL;
}

// ============================= ~cTruss3DTL ===============================

cTruss3DTL :: ~cTruss3DTL(void)
{
}

// =============================== IntForce ================================

int cTruss3DTL :: IntForce(cVector &g)
{
  // Get the nodal coordinates.

  sNodeCoord coord[2];
  Shape->NodalCoord(coord);

  // Get the nodal displacements.

  cVector u(6);
  NodalDispl(u);

  // Evaluate the strain-displacement matrices.

  double L;
  cVector B0(6),Bl(6),B(6),Bi(6);
  B0Matrix(coord, &L, B0);
  BlMatrix(coord, u, &L, Bl);
  B  = B0 + 0.5*Bl;
  Bi = B0 + Bl;

  // Evaluate the normal force.

  cVector eps(1);
  cVector sig(1);
  eps[0] = B*u;
  if (!ConstMod->Stress(eps, sig)) return(0);
  double N = Section->GetA( )*sig[0];

  // Evaluate {g} = [Bi]t*N*L;

  g = (N*L)*Bi;

  return(1);
}

// =============================== StiffMat ================================

void cTruss3DTL :: StiffMat(cMatrix &Kt)
{
  // Get the nodal coordinates.

  sNodeCoord coord[2];
  Shape->NodalCoord(coord);

  // Get the nodal displacements.

  cVector u(6);
  NodalDispl(u);

  // Evaluate the strain-displacement matrices.

  double L;
  cVector B0(6),Bl(6),B(6),Bi(6);
  B0Matrix(coord, &L, B0);
  BlMatrix(coord, u, &L, Bl);
  B  = B0 + 0.5*Bl;
  Bi = B0 + Bl;

  // Evaluate the tangent elastic modulus.

  cMatrix Ct(1, 1);
  ConstMod->TangMat(Ct);

  // Evaluate [Ke] = (E*A*L)*[Bi]t*[Bi]

  int i,j;
  double EA  = Ct[0][0]*Section->GetA( );
  double EAL = EA*L;
  for (i = 0; i < 6; i++)
    for (j = 0; j < 6; j++) Kt[i][j] = EAL*Bi[i]*Bi[j];

  // Evaluate the normal force.

  cVector eps(1);
  cVector sig(1);
  eps[0] = B*u;
  ConstMod->Stress(eps, sig);
  double N = Section->GetA( )*sig[0];

  // Add the geometric stiffness matrix: [Kg] = (N/L)*[G]t*[G].

  double NL = N/L;
  for (i = 0; i < 6; i++) Kt[i][i] += NL;
  Kt[0][3] -= NL;
  Kt[3][0] -= NL;
  Kt[1][4] -= NL;
  Kt[4][1] -= NL;
  Kt[2][5] -= NL;
  Kt[5][2] -= NL;
}

// ============================== NodalStress ==============================

void cTruss3DTL :: NodalStress(cMatrix &S)
{
  // Get the nodal coordinates.

  sNodeCoord coord[2];
  Shape->NodalCoord(coord);

  // Get the nodal displacements.

  cVector u(6);
  NodalDispl(u);

  // Evaluate the strain-displacement matrices.

  double L;
  cVector B0(6),Bl(6),B(6);
  B0Matrix(coord, &L, B0);
  BlMatrix(coord, u, &L, Bl);
  B = B0 + 0.5*Bl;

  // Evaluate the normal force.

  cVector eps(1);
  cVector sig(1);
  eps[0] = B*u;
  ConstMod->Stress(eps, sig);
  double N = Section->GetA( )*sig[0];

  // Store in the given matrix.

  if (OutputConv == DIRECT_STIFFNESS)
  {   
    S[0][0] = -N;
    S[1][0] =  N;
  }
  else
  {
    S[0][0] = S[1][0] = N;
  }
}

// -------------------------------------------------------------------------
// Protected methods:
//

// ================================ BlMatrix ===============================
//
// This method computes the element nonlinear strain-displacement matrix.
//
//   coord - vector of nodal coordinates                               (in)
//   u     - vector of nodal displacements                             (in)
//   L     - element length                                           (out)
//   Bl    - nonlinear strain displacement matrix                     (out)
//
void cTruss3DTL :: BlMatrix(sNodeCoord *coord, cVector &u,
                            double *L, cVector &Bl)
{
  // Evaluate element length.

  double dx = coord[1].x - coord[0].x;
  double dy = coord[1].y - coord[0].y;
  double dz = coord[1].z - coord[0].z;
  double L2 = dx*dx + dy*dy + dz*dz;
  *L = sqrt(L2);

  // Evaluate the displacement differences.

  double du = u[3] - u[0];
  double dv = u[4] - u[1];
  double dw = u[5] - u[2];

  // Evaluate [Bl].

  Bl[0] = -du/L2;
  Bl[1] = -dv/L2;
  Bl[2] = -dw/L2;
  Bl[3] = -Bl[0];
  Bl[4] = -Bl[1];
  Bl[5] = -Bl[2];
}


// -------------------------------------------------------------------------
// Class cTruss3DCR:
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Public methods:
//

// ============================== cTruss3DCR ===============================

cTruss3DCR :: cTruss3DCR(int id) : cTruss3D(id)
{
  Type = TRUSS3DCR;
}

// ============================= ~cTruss3DCR ===============================

cTruss3DCR :: ~cTruss3DCR(void)
{
}

// =============================== IntForce ================================

int cTruss3DCR :: IntForce(cVector &g)
{
  // Get the nodal coordinates.

  sNodeCoord coord[2];
  Shape->NodalCoord(coord);

  // Get the nodal displacements.

  cVector u(6);
  NodalDispl(u);

  // Evaluate l0, ln and {r}.

  double l0,ln;
  cVector r(6);
  RVector(coord, u, &l0, &ln, r);

  // Evaluate the normal force.

  cVector eps(1);
  cVector sig(1);
  eps[0] = ln/l0 - 1.0;
  if (!ConstMod->Stress(eps, sig)) return(0);
  double N = Section->GetA( )*sig[0];

  // Evaluate {g} = N*{r}.

  g = N*r;

  return(1);
}

// =============================== StiffMat ================================

void cTruss3DCR :: StiffMat(cMatrix &Kt)
{
  // Get the nodal coordinates.

  sNodeCoord coord[2];
  Shape->NodalCoord(coord);

  // Get the nodal displacements.

  cVector u(6);
  NodalDispl(u);

  // Evaluate l0, ln and {r}.

  double l0,ln;
  cVector r(6);
  RVector(coord, u, &l0, &ln, r);

  // Evaluate the tangent elastic modulus.

  cMatrix Ct(1, 1);
  ConstMod->TangMat(Ct);

  // Evaluate [Ke] = (EA/l0)*{r}*{r}t.

  int i,j;
  double EA = Ct[0][0]*Section->GetA( );
  double EAl0 = EA/l0;
  for (i = 0; i < 6; i++)
    for (j = 0; j < 6; j++) Kt[i][j] = EAl0*r[i]*r[j];

  // Evaluate the normal force.

  cVector eps(1);
  cVector sig(1);
  eps[0] = ln/l0 - 1.0;
  ConstMod->Stress(eps, sig);
  double N = Section->GetA( )*sig[0];

  // Add the first geometric stiffness matrix: [Kg1] = (N/ln)*[G]t*[G].

  double Nln = N/ln;
  for (i = 0; i < 6; i++) Kt[i][i] += Nln;
  Kt[0][3] -= Nln;
  Kt[3][0] -= Nln;
  Kt[1][4] -= Nln;
  Kt[4][1] -= Nln;
  Kt[2][5] -= Nln;
  Kt[5][2] -= Nln;

  // Add the second geometric stiffness matrix [Kg2] = -(N/ln)*{r}*{r}t

  for (i = 0; i < 6; i++)
    for (j = 0; j < 6; j++) Kt[i][j] -= Nln*r[i]*r[j];
}

// ============================== NodalStress ==============================

void cTruss3DCR :: NodalStress(cMatrix &S)
{
  // Get the nodal coordinates.

  sNodeCoord coord[2];
  Shape->NodalCoord(coord);

  // Get the nodal displacements.

  cVector u(6);
  NodalDispl(u);

  // Evaluate {r}.

  double l0,ln;
  cVector r(6);
  RVector(coord, u, &l0, &ln, r);

  // Evaluate the normal force.

  cVector eps(1);
  cVector sig(1);
  eps[0] = ln/l0 - 1.0;
  ConstMod->Stress(eps, sig);
  double N = Section->GetA( )*sig[0];

  // Store in the given matrix.

  if (OutputConv == DIRECT_STIFFNESS)
  {   
    S[0][0] = -N;
    S[1][0] =  N;
  }
  else
  {
    S[0][0] = S[1][0] = N;
  }
}

// -------------------------------------------------------------------------
// Protected methods:
//

// ================================ RVector ================================

void cTruss3DCR :: RVector(sNodeCoord *coord, cVector &u, double *l0,
                           double *ln, cVector &r)
{
  // Evaluate the initial length.

  double dx = coord[1].x - coord[0].x;
  double dy = coord[1].y - coord[0].y;
  double dz = coord[1].z - coord[0].z;
  *l0 = sqrt(dx*dx + dy*dy + dz*dz);

  // Evaluate the displacement differences.

  double du = u[3] - u[0];
  double dv = u[4] - u[1];
  double dw = u[5] - u[2];

  // Evaluate the current length.

  double dxl = dx + du;
  double dyl = dy + dv;
  double dzl = dz + dw;
  *ln = sqrt(dxl*dxl + dyl*dyl + dzl*dzl);

  // Evaluate {r}.

  r[0] = -dxl/(*ln);
  r[1] = -dyl/(*ln);
  r[2] = -dzl/(*ln);
  r[3] = -r[0];
  r[4] = -r[1];
  r[5] = -r[2];
}

// ======================================================= End of file =====
