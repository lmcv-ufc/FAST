// -------------------------------------------------------------------------
// elmpltrs.cpp - implementation of plane truss elements.
// -------------------------------------------------------------------------
// Created:      16-Nov-2000     Evandro Parente Junior
//
// Modified:     26-Nov-2000     Evandro Parente Junior
//               Implementation of Total Lagrangian element.
//
// Modified:     03-Oct-2001     Evandro Parente Junior
//               Use of vector/matrix classes.
//
// Modified:     28-Jan-2002     Evandro Parente Junior
//               Implementation of the corotational formulation.
//
// Modified:     02-Aug-2014     Evandro Parente Junior
//               Implementation of material nonlinearity.
//
// Modified:     08-Apr-2015     Elias Saraiva Barroso
//               Input/output operations using c++ streams.
// -------------------------------------------------------------------------

#include <iostream>
#include <math.h>
#include <stdlib.h>

using namespace std;

#include "elmpltrs.h"
#include "node.h"
#include "material.h"
#include "cmodel.h"
#include "secbar.h"
#include "shpline.h"
#include "vec.h"
#include "gblvar.h"
#include "gbldef.h"

// -------------------------------------------------------------------------
// Class cPlTruss:
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Public methods:
//

// =============================== cPlTruss ================================

cPlTruss :: cPlTruss(int id) : cElement(id)
{
  Type = PLTRUSS;
  Shape = new cShapeBar( );
}

// ============================== ~cPlTruss ================================

cPlTruss :: ~cPlTruss(void)
{
  delete ConstMod;
}

// ================================= Read ==================================

void cPlTruss :: Read(void)
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

void cPlTruss :: GetActDir(int *dir)
{
  dir[0] = 1;   // u
  dir[1] = 1;   // v
  dir[2] = 0;   // w
  dir[3] = 0;   // rx
  dir[4] = 0;   // ry
  dir[5] = 0;   // rz
}

// ============================= GetStrLabels ==============================

void cPlTruss :: GetStrLabels(int *label)
{
  label[0] = FORCE_X;
}

// =============================== IntForce ================================

int cPlTruss :: IntForce(cVector &g)
{
  // Get the nodal coordinates.

  sNodeCoord coord[2];
  Shape->NodalCoord(coord);

  // Get the nodal displacements.

  cVector u(4);
  NodalDispl(u);

  // Evaluate the strain-displacement matrix.

  double L;
  cVector B0(4);
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

void cPlTruss :: StiffMat(cMatrix &K)
{
  // Get the nodal coordinates.

  sNodeCoord coord[2];
  Shape->NodalCoord(coord);

  // Evaluate the strain-displacement matrix.

  double L;
  cVector B0(4);
  B0Matrix(coord, &L, B0);

  // Evaluate the tangent elastic modulus.

  cMatrix Ct(1, 1);
  ConstMod->TangMat(Ct);

  // Evaluate [Kt] = A*L*[B]t*E*[B]

  double EAL = Ct[0][0]*Section->GetA( )*L;
  for (int i = 0; i < 4; i++)
    for (int j = 0; j < 4; j++) K[i][j] = EAL*B0[i]*B0[j];
}

// =============================== MassMat =================================

void cPlTruss :: MassMat(cMatrix &M)
{
  // Compute the bar length.

  sNodeCoord coord[2];
  Shape->NodalCoord(coord);
  double dx = coord[1].x - coord[0].x;
  double dy = coord[1].y - coord[0].y;
  double L  = sqrt(dx*dx + dy*dy);

  // Assembly the lower matrix.

  cMaterial *mat = Section->GetMaterial( );
  double rho = mat->GetDensity( );
  double A = Section->GetA( );
  double coeff = rho*A*L/6.0;
  M.Zero( );
  M[0][0] = M[1][1] = M[2][2] = M[3][3] = 2.0*coeff;
  M[2][0] = M[0][2] = M[3][1] = M[1][3] = coeff;
}

// =============================== GeomStiff ===============================

void cPlTruss :: GeomStiff(cMatrix &G)
{
  cout << "Error: method cPlTruss :: GeomStiff not implemented!\n";
  exit(0);
}

// ============================== NodalStress ==============================

void cPlTruss :: NodalStress(cMatrix &S)
{
  // Get the nodal coordinates.

  sNodeCoord coord[2];
  Shape->NodalCoord(coord);

  // Get the nodal displacements.

  cVector u(4);
  NodalDispl(u);

  // Evaluate the strain-displacement matrix.

  double L;
  cVector B0(4);
  B0Matrix(coord, &L, B0);

  // Evaluate the normal force.

  cVector eps(1);
  cVector sig(1);
  eps[0] = B0*u;
  ConstMod->Stress(eps, sig);
  double N = Section->GetA( )*sig[0];

  // Store in the given matrix.

  S[0][0] = S[1][0] = N;
}

// ============================== UpdateState ==============================

void cPlTruss :: UpdateState(void)
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
void cPlTruss :: B0Matrix(sNodeCoord *coord, double *L, cVector &B0)
{
  // Evaluate element length.

  double dx = coord[1].x - coord[0].x;
  double dy = coord[1].y - coord[0].y;
  double L1 = sqrt(dx*dx + dy*dy);
  double L2 = L1*L1;
  *L = L1;

  // Evaluate [B0].

  B0[0] = -dx/L2;
  B0[1] = -dy/L2;
  B0[2] = -B0[0];
  B0[3] = -B0[1];
}


// -------------------------------------------------------------------------
// Class cPlTrussTL:
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Public methods:
//

// ============================== cPlTrussTL ===============================

cPlTrussTL :: cPlTrussTL(int id) : cPlTruss(id)
{
  Type = PLTRUSSTL;
}

// ============================= ~cPlTrussTL ===============================

cPlTrussTL :: ~cPlTrussTL(void)
{
}

// =============================== IntForce ================================

int cPlTrussTL :: IntForce(cVector &g)
{
  // Get the nodal coordinates.

  sNodeCoord coord[2];
  Shape->NodalCoord(coord);

  // Get the nodal displacements.

  cVector u(4);
  NodalDispl(u);

  // Evaluate the strain-displacement matrices.

  double L;
  cVector B0(4),Bl(4),B(4),Bi(4);
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

void cPlTrussTL :: StiffMat(cMatrix &Kt)
{
  // Get the nodal coordinates.

  sNodeCoord coord[2];
  Shape->NodalCoord(coord);

  // Get the nodal displacements.

  cVector u(4);
  NodalDispl(u);

  // Evaluate the strain-displacement matrices.

  double L;
  cVector B0(4),Bl(4),B(4),Bi(4);
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
  for (i = 0; i < 4; i++)
    for (j = 0; j < 4; j++) Kt[i][j] = EAL*Bi[i]*Bi[j];

  // Evaluate the normal force.

  cVector eps(1);
  cVector sig(1);
  eps[0] = B*u;
  ConstMod->Stress(eps, sig);
  double N = Section->GetA( )*sig[0];

  // Add the geometric stiffness matrix: [Kg] = (N/L)*[G]t*[G].

  double NL = N/L;
  for (i = 0; i < 4; i++) Kt[i][i] += NL;
  Kt[0][2] -= NL;
  Kt[2][0] -= NL;
  Kt[1][3] -= NL;
  Kt[3][1] -= NL;
}

// ============================== NodalStress ==============================

void cPlTrussTL :: NodalStress(cMatrix &S)
{
  // Get the nodal coordinates.

  sNodeCoord coord[2];
  Shape->NodalCoord(coord);

  // Get the nodal displacements.

  cVector u(4);
  NodalDispl(u);

  // Evaluate the strain-displacement matrices.

  double L;
  cVector B0(4),Bl(4),B(4);
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

  S[0][0] = S[1][0] = N;
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
void cPlTrussTL :: BlMatrix(sNodeCoord *coord, cVector &u,
                            double *L, cVector &Bl)
{
  // Evaluate element length.

  double dx = coord[1].x - coord[0].x;
  double dy = coord[1].y - coord[0].y;
  double L1 = sqrt(dx*dx + dy*dy);
  double L2 = L1*L1;
  *L = L1;

  // Evaluate the displacement differences.

  double du = u[2] - u[0];
  double dv = u[3] - u[1];

  // Evaluate [Bl].

  Bl[0] = -du/L2;
  Bl[1] = -dv/L2;
  Bl[2] = -Bl[0];
  Bl[3] = -Bl[1];
}


// -------------------------------------------------------------------------
// Class cPlTrussCR:
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Public methods:
//

// ============================== cPlTrussCR ===============================

cPlTrussCR :: cPlTrussCR(int id) : cPlTruss(id)
{
  Type = PLTRUSSCR;
}

// ============================= ~cPlTrussCR ===============================

cPlTrussCR :: ~cPlTrussCR(void)
{
}

// =============================== IntForce ================================

int cPlTrussCR :: IntForce(cVector &g)
{
  // Get the nodal coordinates.

  sNodeCoord coord[2];
  Shape->NodalCoord(coord);

  // Get the nodal displacements.

  cVector u(4);
  NodalDispl(u);

  // Evaluate l0, ln and {r}.

  double l0,ln;
  cVector r(4);
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

void cPlTrussCR :: StiffMat(cMatrix &Kt)
{
  // Get the nodal coordinates.

  sNodeCoord coord[2];
  Shape->NodalCoord(coord);

  // Get the nodal displacements.

  cVector u(4);
  NodalDispl(u);

  // Evaluate l0, ln and {r}.

  double l0,ln;
  cVector r(4);
  RVector(coord, u, &l0, &ln, r);

  // Evaluate the tangent elastic modulus.

  cMatrix Ct(1, 1);
  ConstMod->TangMat(Ct);

  // Evaluate [Ke] = (EA/l0)*{r}*{r}t.

  int i,j;
  double EA = Ct[0][0]*Section->GetA( );
  double EAl0 = EA/l0;
  for (i = 0; i < 4; i++)
    for (j = 0; j < 4; j++) Kt[i][j] = EAl0*r[i]*r[j];

  // Evaluate the normal force.

  cVector eps(1);
  cVector sig(1);
  eps[0] = ln/l0 - 1.0;
  ConstMod->Stress(eps, sig);
  double N = Section->GetA( )*sig[0];

  // Add the first geometric stiffness matrix: [Kg1] = (N/ln)*[G]t*[G].

  double Nln = N/ln;
  for (i = 0; i < 4; i++) Kt[i][i] += Nln;
  Kt[0][2] -= Nln;
  Kt[2][0] -= Nln;
  Kt[1][3] -= Nln;
  Kt[3][1] -= Nln;

  // Add the second geometric stiffness matrix [Kg2] = -(N/ln)*{r}*{r}t

  for (i = 0; i < 4; i++)
    for (j = 0; j < 4; j++) Kt[i][j] -= Nln*r[i]*r[j];
}

// ============================== NodalStress ==============================

void cPlTrussCR :: NodalStress(cMatrix &S)
{
  // Get the nodal coordinates.

  sNodeCoord coord[2];
  Shape->NodalCoord(coord);

  // Get the nodal displacements.

  cVector u(4);
  NodalDispl(u);

  // Evaluate {r}.

  double l0,ln;
  cVector r(4);
  RVector(coord, u, &l0, &ln, r);

  // Evaluate the normal force.

  cVector eps(1);
  cVector sig(1);
  eps[0] = ln/l0 - 1.0;
  ConstMod->Stress(eps, sig);
  double N = Section->GetA( )*sig[0];

  // Store in the given matrix.

  S[0][0] = S[1][0] = N;
}

// -------------------------------------------------------------------------
// Protected methods:
//

// ================================ RVector ================================

void cPlTrussCR :: RVector(sNodeCoord *coord, cVector &u, double *l0,
                           double *ln, cVector &r)
{
  // Evaluate the initial length.

  double dx = coord[1].x - coord[0].x;
  double dy = coord[1].y - coord[0].y;
  *l0 = sqrt(dx*dx + dy*dy);

  // Evaluate the displacement differences.

  double du = u[2] - u[0];
  double dv = u[3] - u[1];

  // Evaluate the current length.

  double dxl = dx + du;
  double dyl = dy + dv;
  *ln = sqrt(dxl*dxl + dyl*dyl);

  // Evaluate {r}.

  r[0] = -dxl/(*ln);
  r[1] = -dyl/(*ln);
  r[2] = -r[0];
  r[3] = -r[1];
}

// ======================================================= End of file =====
