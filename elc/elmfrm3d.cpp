// -------------------------------------------------------------------------
// elmplfrm3d.cpp - implementation of three dimensional frame classes.
// -------------------------------------------------------------------------
// Created:      12-Nov-2012     Luiz Antonio Taumaturgo Mororo
//
// Modified:     
//             
// -------------------------------------------------------------------------

#include <math.h>
#include <stdlib.h>
#include <iostream>

using namespace std;

#include "elmfrm3d.h"
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
// cFrame3D class:
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Public methods:
//

// ============================= ValidSection ==============================

bool cFrame3D :: ValidSection(eSecType type)
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

// =============================== cFrame3D ================================

cFrame3D :: cFrame3D(int id, cSection *sec) : cElement(id)
{
  Type = FRAME3D;
  Section = (cSecBar *)sec;
  Shape = new cShapeBar( );
}

// ============================== ~cFrame3D ================================

cFrame3D :: ~cFrame3D(void)
{
}

// ================================= Read ==================================

void cFrame3D :: Read(void)
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

void cFrame3D :: GetActDir(int *dir)
{
  dir[0] = 1;   // u
  dir[1] = 1;   // v
  dir[2] = 1;   // w
  dir[3] = 1;   // rx
  dir[4] = 1;   // ry
  dir[5] = 1;   // rz
}

// ============================= GetStrLabels ==============================

void cFrame3D :: GetStrLabels(int *label)
{
  label[0] = FORCE_X;
  label[1] = FORCE_Y;
  label[2] = FORCE_Z;
  label[3] = MOMENT_X;
  label[4] = MOMENT_Y;
  label[5] = MOMENT_Z;
}

// =============================== IntForce ================================

int cFrame3D :: IntForce(cVector &g)
{
  // Get the nodal coordinates.
  
  sNodeCoord coord[2];
  Shape->NodalCoord(coord);
    
  // Compute the element length.
  
  double L = CalcLength(coord);
  
  // Assembly the transformation matrix.
  
  cMatrix T(12, 12);
  GetTrnMat(T);
  
  // Get the nodal displacements.
  
  cVector u(12);
  NodalDispl(u);
  
  // Compute the local displacements => {ul} = [T]{u}.
  
  cVector ul(12);
  ul = T*u;
  
  // Assembly the local stiffness matrix.
  
  cMatrix Ke(12, 12);
  LocStiffMat(L, Ke);
  
  // Compute the local forces => {gl} = [Ke]{ul}.
  
  cVector gl(12);
  gl = Ke*ul;
  
  // Transform to the global system => {g} = [T]t{gl}
  
  g = t(T)*gl;
  
  return (1);
}

// =============================== StiffMat ================================

void cFrame3D :: StiffMat(cMatrix &K)
{
  // Get te nodal coordinates.
  
  sNodeCoord coord[2];
  Shape->NodalCoord(coord);
  
  // Compute the element length.
  
  double L = CalcLength(coord);
  
  // Assembly the transformation matrix.
  
  cMatrix T(12, 12);
  GetTrnMat(T);
  
  // Assembly the local stiffness matrix.
  
  cMatrix Ke(12, 12);
  LocStiffMat(L, Ke);
  
  // Compute the stiffness matrix => [K] = [T]t[Ke][T].
    
  K.Zero( );
  MatTripBtCB(T, Ke, 1.0, K);
}

// ============================== NodalStress ==============================

void cFrame3D :: NodalStress(cMatrix &S)
{
  // Get the nodal coordinates.
  
  sNodeCoord coord[2];
  Shape->NodalCoord(coord);
  
  // Assembly the transformation matrix.
  
  cMatrix T(12, 12);
  GetTrnMat(T);
  
  // Compute the internal force vector.
  
  cVector g(12);
  IntForce(g);
  
  // Transform to the local system => {gl} = [T]{g}.
  
  cVector gl(12);
  gl = T*g;
  
  // Add the fixed-end forces (- equivforces).
  
  cVector fl(12);
  EqvForces(fl);
  gl -= fl;
  
  // Return nodal forces.

  if (OutputConv == DIRECT_STIFFNESS)
  {
    S[0][0] = gl[ 0];  // Initial node
    S[0][1] = gl[ 1];
    S[0][2] = gl[ 2];
    S[0][3] = gl[ 3];
    S[0][4] = gl[ 4];
    S[0][5] = gl[ 5];
    S[1][0] = gl[ 6];  // Final node
    S[1][1] = gl[ 7];
    S[1][2] = gl[ 8];
    S[1][3] = gl[ 9];
    S[1][4] = gl[10];
    S[1][5] = gl[11];
  }
  else
  {
    S[0][0] = -gl[ 0];  // Initial node
    S[0][1] =  gl[ 1];
    S[0][2] =  gl[ 2];
    S[0][3] = -gl[ 3];
    S[0][4] = -gl[ 4];
    S[0][5] =  gl[ 5];
    S[1][0] =  gl[ 6];  // Final node
    S[1][1] = -gl[ 7];
    S[1][2] = -gl[ 8];
    S[1][3] =  gl[ 9];
    S[1][4] =  gl[10];
    S[1][5] = -gl[11];
  }
}

// ============================== CalcRotMat ==============================

void cFrame3D :: CalcRotMat(cMatrix &R)
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
double cFrame3D :: CalcLength(sNodeCoord *coord)
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
//   T  - transformation matrix (12x12)                               (out)
//
void cFrame3D :: GetTrnMat(cMatrix &T)
{
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
}

// ============================== LocSitffMat ==============================
//
// This method computes the element sitffness matrix the local system.
//
//   L  - element length                                               (in)
//   Kl - elastic matrix (12x12)                                      (out)
//
void cFrame3D :: LocStiffMat(double L, cMatrix &Kl)
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

    double EA  = E*Section->GetA( );
    double EIy = E*Section->GetIy( );
    double EIz = E*Section->GetIz( );
    double GJ  = G*Section->GetJt( );
  
    double L2 = L*L;
    double L3 = L2*L;
  
    Kl[ 0][ 0] = Kl[ 6][ 6] =  EA/L;
    Kl[ 0][ 6] = Kl[ 6][ 0] = -EA/L;
  
    Kl[ 1][ 1] =  12*EIz/L3;
    Kl[ 1][ 5] = Kl[ 5][ 1] =   6*EIz/L2;
    Kl[ 1][ 7] = Kl[ 7][ 1] = -12*EIz/L3;
    Kl[ 1][11] = Kl[11][ 1] =   6*EIz/L2;
  
    Kl[ 2][ 2] =  12*EIy/L3;
    Kl[ 2][ 4] = Kl[ 4][ 2] =  -6*EIy/L2;
    Kl[ 2][ 8] = Kl[ 8][ 2] = -12*EIy/L3;
    Kl[ 2][10] = Kl[10][ 2] =  -6*EIy/L2;
  
    Kl[ 3][ 3] = Kl[ 9][ 9] =  GJ/L;
    Kl[ 3][ 9] = Kl[ 9][ 3] = -GJ/L;
  
    Kl[ 4][ 4] =  4*EIy/L;
    Kl[ 4][ 8] = Kl[ 8][ 4] = 6*EIy/L2;
    Kl[ 4][10] = Kl[10][ 4] = 2*EIy/L;
  
    Kl[ 5][ 5] =  4*EIz/L;
    Kl[ 5][ 7] = Kl[ 7][ 5] = -6*EIz/L2;
    Kl[ 5][11] = Kl[11][ 5] =  2*EIz/L;
  
    Kl[ 7][ 7] =  12*EIz/L3;
    Kl[ 7][11] = Kl[11][7] = -6*EIz/L2;
  
    Kl[ 8][ 8] =  12*EIy/L3;
    Kl[ 8][10] = Kl[10][8] = 6*EIy/L2;
  
    Kl[10][10] =  4*EIy/L;
  
    Kl[11][11] =  4*EIz/L;
  }
  else // Coupled section constitutive matrix
  {
    cout << "Error: Frame3D with coupled section not implemented yet!";
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
void cFrame3D :: EqvForces(cVector &fl)
{
  // Get the equivalent nodal forces in the global system.
  
  cVector fg(12);
  fg.Zero( );
  double t = cControl::GetTotTime( );
  cLoad :: EqvForces(this, t, fg);
  
  // Assembly the transformation matrix.
  
  cMatrix T(12, 12);
  GetTrnMat(T);
  
  // Rotate fix end forces to local system => {fl} = [T]{fg}.
  
  fl = T*fg;
}


// -------------------------------------------------------------------------
// cFrame3DTL class:
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Public methods:
//

// ============================== cFrame3DTL ===============================

cFrame3DTL :: cFrame3DTL(int id, cSection *sec) : cFrame3D(id,sec)
{
  Type = FRAME3DTL;
}

// ============================= ~cFrame3DTL ===============================

cFrame3DTL :: ~cFrame3DTL(void)
{
}

// =============================== IntForce ================================

int cFrame3DTL :: IntForce(cVector &g)
{
  // Get the nodal coordinates.
  
  sNodeCoord coord[2];
  Shape->NodalCoord(coord);
  
  // Compute the element length.
  
  double L = CalcLength(coord);
  
  // Assembly the transformation matrix.
  
  cMatrix T(12, 12);
  GetTrnMat(T);
    
  // Get the nodal displacements.
  
  cVector u(12);
  NodalDispl(u);
  
  // Compute the local displacements => {ul} = [T]{u}.
  
  cVector ul(12);
  ul = T*u;
  
  // Compute the stiffness and geometric properties.

  cMatMec *mat = Section->GetMaterial( )->GetMecProp( );
  double matpar[2];
  mat->GetParam(matpar);
  double A  = Section->GetA( );
  double Ip = Section->GetJt( );

  // Compute the internal local force.
  
  cVector gl(12);
  LocIntForce(L, Ip, A, ul, gl);
  
  // Transform to the global system => {g} = [T]t{gl}.
  
  g = t(T)*gl;
  
  return (1);
}

// =============================== StiffMat ================================

void cFrame3DTL :: StiffMat(cMatrix &K)
{
  // Compute the stiffness and geometric properties.
  
  cMatMec *mat = Section->GetMaterial( )->GetMecProp( );
  double matpar[2];
  mat->GetParam(matpar);
  
  double A  = Section->GetA( );
  double Ip = Section->GetJt( );

  // Get the nodal coordinates.
  
  sNodeCoord coord[2];
  Shape->NodalCoord(coord);
  
  // Compute the element length.
  
  double L = CalcLength(coord);
  
  // Assembly the transformation matrix.
  
  cMatrix T(12, 12);
  GetTrnMat(T);
  
  // Get the nodal displacements.
  
  cVector u(12);
  NodalDispl(u);
  
  // Compute the local displacements => {ul} = [T]{u}.
  
  cVector ul(12);
  ul = T*u;
  
  // Compute the "geometric" stiffness matrix => [kg] = NL*[A]t.
  
  cMatrix At(12, 12);
  AMat(L, Ip, A, At);
  
  cMatrix kg(12, 12);
  double NL = NormalForce(L, Ip, A, ul);  // Esta correto ???
  kg = NL*At;
  
  // Compute the [kl] matrix => [Kl] = int([B]t[C][B]).
  
  cMatrix kl(12, 12);
  LocStiffMat(L, Ip, A, ul, kl);
  
  // Compute the [ke] matrix => [ke] = [kl] + [kg].
  
  cMatrix ke(12, 12);
  ke = kl + kg;
  
  // Compute the stiffness matrix => [K] = [T]t[ke][T].
  
  K.Zero( );
  MatTripBtCB(T, ke, 1.0, K);
}

// -------------------------------------------------------------------------
// Protected methods:
//

// ================================ LocIntForce ============================
//
// This method computes the internal force vector.
//
//   L  - element length                                               (in)
//   Ip - polar moment of inertia of the cross section                 (in)
//   A  - area of the cross section                                    (in)
//   u  - local displacement vector (12)                               (in)
//   gl - local internal force vector (12)                            (out)
//
void cFrame3DTL :: LocIntForce(double L, double Ip, double A, cVector &u, 
                               cVector &gl)
{
  // Classical formulation (diagonal section matrix).

  if (Section->GetType( ) != SEC_BAR_B3DCOUPLED)
  {
    cout << "Error: Frame3DTL with uncoupled section not implemented yet!";
    exit(0); 
  }

  // Section stiffness - fully coupled 4x4 matrix.
  
  cMatrix Cb(4, 4);
  Section->SecStiffness3D(Cb);
  double C11 = Cb[0][0]; double C12 = Cb[0][1]; double C13 = Cb[0][2]; double C14 = Cb[0][3];
  double C21 = Cb[1][0]; double C22 = Cb[1][1]; double C23 = Cb[1][2]; double C24 = Cb[1][3];
  double C31 = Cb[2][0]; double C32 = Cb[2][1]; double C33 = Cb[2][2]; double C34 = Cb[2][3];
  double C41 = Cb[3][0]; double C42 = Cb[3][1]; double C43 = Cb[3][2]; double C44 = Cb[3][3];
  
  // Auxiliary variables (Maple code).

  double t1 = L * L;
  double t2 = t1 * L;
  double t3 = 0.1e1 / t2;
  double t4 = t3 * C13;
  double t7 = t1 * t1;
  double t8 = 0.1e1 / t7;
  double t9 = t8 * C13;
  double t12 = t3 * C12;
  double t15 = t8 * C12;
  double t26 = 0.6e1 * t4 * u[4] - 0.12e2 * t9 * u[2] - 0.6e1 * t12 * u[5] - 0.12e2 * t15 * u[1] + 0.12e2 * t15 * u[7] + 0.6e1 * t4 * u[10] + 0.12e2 * t9 * u[8] - 0.6e1 * t12 * u[11];
  double t29 = 0.1e1 / L;
  double t30 = t29 * C11;
  double t31 = t30 * u[0];
  double t32 = 0.1e1 / t1;
  double t33 = u[1] * t32;
  double t35 = u[5] * t29;
  double t37 = u[7] * t32;
  double t39 = u[11] * t29;
  double t41 = 0.3e1 / 0.5e1 * t33 + t35 / 0.20e2 - 0.3e1 / 0.5e1 * t37 + t39 / 0.20e2;
  double t43 = 0.6e1 * t12;
  double t44 = -t30 * t41 + t43;
  double t49 = Ip / A;
  double t53 = -u[3] * t32 * t49 + u[9] * t32 * t49;
  double t55 = t32 * C14;
  double t56 = -t30 * t53 / 0.2e1 - t55;
  double t59 = u[2] * t32;
  double t61 = u[4] * t29;
  double t63 = u[8] * t32;
  double t65 = u[10] * t29;
  double t67 = 0.3e1 / 0.5e1 * t59 - t61 / 0.20e2 - 0.3e1 / 0.5e1 * t63 - t65 / 0.20e2;
  double t69 = 0.6e1 * t4;
  double t70 = -t30 * t67 + t69;
  double t73 = u[1] * t29;
  double t74 = t73 / 0.20e2;
  double t76 = u[7] * t29;
  double t77 = t76 / 0.20e2;
  double t79 = t74 + u[5] / 0.15e2 - t77 - u[11] / 0.60e2;
  double t81 = t32 * C12;
  double t83 = -t30 * t79 + 0.4e1 * t81;
  double t87 = t30 * t53 / 0.2e1 + t55;
  double t90 = u[2] * t29;
  double t91 = t90 / 0.20e2;
  double t93 = u[8] * t29;
  double t94 = t93 / 0.20e2;
  double t96 = -t91 + u[4] / 0.15e2 + t94 - u[10] / 0.60e2;
  double t98 = t32 * C13;
  double t100 = -t30 * t96 - 0.4e1 * t98;
  double t104 = t30 * t41 - t43;
  double t109 = -t91 - u[4] / 0.60e2 + t94 + u[10] / 0.15e2;
  double t112 = -t30 * t109 - 0.2e1 * t98;
  double t116 = t30 * t67 - t69;
  double t119 = t30 * u[6];
  double t122 = t74 - u[5] / 0.60e2 - t77 + u[11] / 0.15e2;
  double t125 = -t30 * t122 + 0.2e1 * t81;
  double t128 = t26 * t1 / 0.2e1 + t31 + t44 * u[1] * L + t56 * u[9] * L + t70 * u[2] * L + t83 * u[5] * L + t87 * u[3] * L + t100 * u[4] * L + t104 * u[7] * L + t112 * u[10] * L + t116 * u[8] * L - t119 + t125 * u[11] * L;
  double t130 = 0.1e1 / t7 / t1;
  double t131 = t130 * C23;
  double t135 = 0.1e1 / t7 / L;
  double t136 = t135 * C22;
  double t139 = t130 * C22;
  double t144 = t135 * C23;
  double t153 = 0.144e3 * t131 * u[2] + 0.72e2 * t136 * u[5] - 0.144e3 * t139 * u[7] - 0.144e3 * t131 * u[8] - 0.72e2 * t144 * u[10] - 0.72e2 * t144 * u[4] + 0.144e3 * t139 * u[1] + 0.72e2 * t136 * u[11];
  double t156 = t3 * C21;
  double t158 = -0.12e2 * t156 * t41;
  double t163 = 0.6e1 / 0.5e1 * t33 + t35 / 0.10e2 - 0.6e1 / 0.5e1 * t37 + t39 / 0.10e2;
  double t166 = 0.6e1 * t32 * C22;
  double t167 = t163 * C12 - t166;
  double t169 = 0.12e2 * t167 * t3;
  double t170 = 0.72e2 * t136;
  double t173 = t8 * C21;
  double t175 = 0.12e2 * t173 * u[6];
  double t177 = 0.12e2 * t173 * u[0];
  double t179 = 0.12e2 * t156 * t109;
  double t182 = 0.6e1 * t32 * C23;
  double t183 = t163 * C13 - t182;
  double t185 = 0.6e1 * t183 * t32;
  double t186 = t8 * C23;
  double t187 = 0.24e2 * t186;
  double t191 = t8 * C24;
  double t192 = t156 * t53 / 0.2e1 + t191;
  double t195 = -t156 * t53 / 0.2e1 - t191;
  double t198 = 0.12e2 * t156 * t67;
  double t200 = 0.12e2 * t183 * t3;
  double t201 = 0.72e2 * t144;
  double t205 = 0.12e2 * t156 * t96;
  double t206 = 0.48e2 * t186;
  double t210 = 0.12e2 * t156 * t122;
  double t212 = 0.6e1 * t167 * t32;
  double t213 = t8 * C22;
  double t214 = 0.24e2 * t213;
  double t218 = 0.12e2 * t156 * t79;
  double t219 = 0.48e2 * t213;
  double t223 = 0.12e2 * t156 * t41;
  double t227 = -0.12e2 * t156 * t67;
  double t230 = (t158 - t169 + t170) * u[7] + t175 - t177 + (t179 - t185 + t187) * u[10] + 0.12e2 * t192 * u[9] + 0.12e2 * t195 * u[3] + (t198 + t200 - t201) * u[2] + (t205 - t185 + t206) * u[4] + (t210 + t212 - t214) * u[11] + (t218 + t212 - t219) * u[5] + (t223 + t169 - t170) * u[1] + (t227 - t200 + t201) * u[8];
  double t234 = t32 * C21;
  double t235 = 0.6e1 * t234;
  double t236 = t163 * C11 - t235;
  double t250 = 0.6e1 * t32 * C24;
  double t252 = (t163 * C14 - t250) * t29;
  double t257 = t167 * t29;
  double t263 = t183 * t29;
  double t290 = t153 * t2 / 0.3e1 + t230 * t1 / 0.2e1 - t236 * u[0] + (t236 * t67 - t185) * u[2] * L + (t236 * t41 - t212) * u[1] * L + t236 * u[6] + (-t236 * t53 / 0.2e1 - t252) * u[3] * L + (t236 * t79 - 0.4e1 * t257) * u[5] * L + (t236 * t96 + 0.4e1 * t263) * u[4] * L + (t236 * t53 / 0.2e1 + t252) * u[9] * L + (t236 * t122 - 0.2e1 * t257) * u[11] * L + (-t236 * t67 + t185) * u[8] * L + (-t236 * t41 + t212) * u[7] * L + (t236 * t109 + 0.2e1 * t263) * u[10] * L;
  double t291 = t130 * C33;
  double t294 = t135 * C32;
  double t297 = t130 * C32;
  double t302 = t135 * C33;
  double t311 = 0.144e3 * t291 * u[2] + 0.72e2 * t294 * u[5] - 0.144e3 * t297 * u[7] - 0.144e3 * t291 * u[8] - 0.72e2 * t302 * u[10] - 0.72e2 * t302 * u[4] + 0.144e3 * t297 * u[1] + 0.72e2 * t294 * u[11];
  double t314 = t3 * C31;
  double t316 = -0.12e2 * t314 * t41;
  double t321 = 0.6e1 / 0.5e1 * t59 - t61 / 0.10e2 - 0.6e1 / 0.5e1 * t63 - t65 / 0.10e2;
  double t324 = 0.6e1 * t32 * C32;
  double t325 = t321 * C12 - t324;
  double t327 = 0.12e2 * t325 * t3;
  double t328 = 0.72e2 * t294;
  double t331 = t8 * C31;
  double t333 = 0.12e2 * t331 * u[6];
  double t335 = 0.12e2 * t331 * u[0];
  double t337 = 0.12e2 * t314 * t109;
  double t340 = 0.6e1 * t32 * C33;
  double t341 = t321 * C13 - t340;
  double t343 = 0.6e1 * t341 * t32;
  double t344 = t8 * C33;
  double t345 = 0.24e2 * t344;
  double t349 = t8 * C34;
  double t350 = t314 * t53 / 0.2e1 + t349;
  double t353 = -t314 * t53 / 0.2e1 - t349;
  double t356 = 0.12e2 * t314 * t67;
  double t358 = 0.12e2 * t341 * t3;
  double t359 = 0.72e2 * t302;
  double t363 = 0.12e2 * t314 * t96;
  double t364 = 0.48e2 * t344;
  double t368 = 0.12e2 * t314 * t122;
  double t370 = 0.6e1 * t325 * t32;
  double t371 = t8 * C32;
  double t372 = 0.24e2 * t371;
  double t376 = 0.12e2 * t314 * t79;
  double t377 = 0.48e2 * t371;
  double t381 = 0.12e2 * t314 * t41;
  double t385 = -0.12e2 * t314 * t67;
  double t388 = (t316 - t327 + t328) * u[7] + t333 - t335 + (t337 - t343 + t345) * u[10] + 0.12e2 * t350 * u[9] + 0.12e2 * t353 * u[3] + (t356 + t358 - t359) * u[2] + (t363 - t343 + t364) * u[4] + (t368 + t370 - t372) * u[11] + (t376 + t370 - t377) * u[5] + (t381 + t327 - t328) * u[1] + (t385 - t358 + t359) * u[8];
  double t392 = t32 * C31;
  double t393 = 0.6e1 * t392;
  double t394 = t321 * C11 - t393;
  double t408 = 0.6e1 * t32 * C34;
  double t410 = (t321 * C14 - t408) * t29;
  double t415 = t325 * t29;
  double t421 = t341 * t29;
  double t448 = t311 * t2 / 0.3e1 + t388 * t1 / 0.2e1 - t394 * u[0] + (t394 * t67 - t343) * u[2] * L + (t394 * t41 - t370) * u[1] * L + t394 * u[6] + (-t394 * t53 / 0.2e1 - t410) * u[3] * L + (t394 * t79 - 0.4e1 * t415) * u[5] * L + (t394 * t96 + 0.4e1 * t421) * u[4] * L + (t394 * t53 / 0.2e1 + t410) * u[9] * L + (t394 * t122 - 0.2e1 * t415) * u[11] * L + (-t394 * t67 + t343) * u[8] * L + (-t394 * t41 + t370) * u[7] * L + (t394 * t109 + 0.2e1 * t421) * u[10] * L;
  double t450 = t29 * C43;
  double t451 = -t53 * C13 - t450;
  double t452 = t451 * t32;
  double t455 = t451 * t3;
  double t459 = t29 * C42;
  double t460 = -t53 * C12 - t459;
  double t461 = t460 * t32;
  double t464 = t460 * t3;
  double t479 = t29 * C41;
  double t480 = -t53 * C11 - t479;
  double t483 = 0.6e1 * t461;
  double t489 = t29 * C44;
  double t491 = (-t53 * C14 - t489) * t29;
  double t496 = 0.6e1 * t452;
  double t501 = t460 * t29;
  double t511 = t451 * t29;
  double t535 = (-0.6e1 * t452 * u[4] + 0.12e2 * t455 * u[2] + 0.6e1 * t461 * u[5] + 0.12e2 * t464 * u[1] - 0.12e2 * t464 * u[7] - 0.6e1 * t452 * u[10] - 0.12e2 * t455 * u[8] + 0.6e1 * t461 * u[11]) * t1 / 0.2e1 - t480 * u[0] + (t480 * t41 - t483) * u[1] * L + (t480 * t53 / 0.2e1 + t491) * u[9] * L + (t480 * t67 - t496) * u[2] * L + (t480 * t79 - 0.4e1 * t501) * u[5] * L + (-t480 * t53 / 0.2e1 - t491) * u[3] * L + (t480 * t96 + 0.4e1 * t511) * u[4] * L + (-t480 * t41 + t483) * u[7] * L + (t480 * t109 + 0.2e1 * t511) * u[10] * L + (-t480 * t67 + t496) * u[8] * L + t480 * u[6] + (t480 * t122 - 0.2e1 * t501) * u[11] * L;
  double t554 = (-0.72e2 * t302 * u[2] - 0.36e2 * t371 * u[5] + 0.72e2 * t294 * u[7] + 0.72e2 * t302 * u[8] + 0.36e2 * t344 * u[10] + 0.36e2 * t344 * u[4] - 0.72e2 * t294 * u[1] - 0.36e2 * t371 * u[11]) * t2 / 0.3e1;
  double t556 = -0.6e1 * t392 * t41;
  double t557 = t90 / 0.10e2;
  double t559 = t93 / 0.10e2;
  double t561 = -t557 + 0.2e1 / 0.15e2 * u[4] + t559 - u[10] / 0.30e2;
  double t563 = t29 * C32;
  double t565 = t561 * C12 + 0.4e1 * t563;
  double t567 = 0.12e2 * t565 * t3;
  double t568 = 0.36e2 * t371;
  double t572 = 0.6e1 * t314 * u[6];
  double t574 = 0.6e1 * t314 * u[0];
  double t576 = 0.6e1 * t392 * t109;
  double t578 = t29 * C33;
  double t580 = t561 * C13 + 0.4e1 * t578;
  double t582 = 0.6e1 * t580 * t32;
  double t583 = t3 * C33;
  double t584 = 0.12e2 * t583;
  double t588 = t3 * C34;
  double t590 = 0.6e1 * (-t392 * t53 / 0.2e1 - t588) * u[9];
  double t593 = 0.6e1 * (t392 * t53 / 0.2e1 + t588) * u[3];
  double t595 = 0.6e1 * t392 * t67;
  double t597 = 0.12e2 * t580 * t3;
  double t598 = 0.36e2 * t344;
  double t602 = 0.6e1 * t392 * t96;
  double t603 = 0.24e2 * t583;
  double t607 = 0.6e1 * t392 * t122;
  double t609 = 0.6e1 * t565 * t32;
  double t610 = t3 * C32;
  double t611 = 0.12e2 * t610;
  double t615 = 0.6e1 * t392 * t79;
  double t616 = 0.24e2 * t610;
  double t620 = 0.6e1 * t392 * t41;
  double t624 = -0.6e1 * t392 * t67;
  double t627 = (-t556 - t567 - t568) * u[7] - t572 + t574 + (-t576 - t582 - t584) * u[10] + t590 + t593 + (-t595 + t597 + t598) * u[2] + (-t602 - t582 - t603) * u[4] + (-t607 + t609 + t611) * u[11] + (-t615 + t609 + t616) * u[5] + (-t620 + t567 + t568) * u[1] + (-t624 - t597 - t598) * u[8];
  double t631 = t29 * C31;
  double t633 = t561 * C11 + 0.4e1 * t631;
  double t646 = t29 * C34;
  double t649 = (t561 * C14 + 0.4e1 * t646) * t29;
  double t654 = t565 * t29;
  double t660 = t580 * t29;
  double t687 = t554 + t627 * t1 / 0.2e1 - t633 * u[0] + (t633 * t67 - t582) * u[2] * L + (t633 * t41 - t609) * u[1] * L + t633 * u[6] + (-t633 * t53 / 0.2e1 - t649) * u[3] * L + (t633 * t79 - 0.4e1 * t654) * u[5] * L + (t633 * t96 + 0.4e1 * t660) * u[4] * L + (t633 * t53 / 0.2e1 + t649) * u[9] * L + (t633 * t122 - 0.2e1 * t654) * u[11] * L + (-t633 * t67 + t582) * u[8] * L + (-t633 * t41 + t609) * u[7] * L + (t633 * t109 + 0.2e1 * t660) * u[10] * L;
  double t706 = (0.72e2 * t144 * u[2] + 0.36e2 * t213 * u[5] - 0.72e2 * t136 * u[7] - 0.72e2 * t144 * u[8] - 0.36e2 * t186 * u[10] - 0.36e2 * t186 * u[4] + 0.72e2 * t136 * u[1] + 0.36e2 * t213 * u[11]) * t2 / 0.3e1;
  double t708 = -0.6e1 * t234 * t41;
  double t709 = t73 / 0.10e2;
  double t711 = t76 / 0.10e2;
  double t713 = t709 + 0.2e1 / 0.15e2 * u[5] - t711 - u[11] / 0.30e2;
  double t715 = t29 * C22;
  double t717 = t713 * C12 - 0.4e1 * t715;
  double t719 = 0.12e2 * t717 * t3;
  double t720 = 0.36e2 * t213;
  double t724 = 0.6e1 * t156 * u[6];
  double t726 = 0.6e1 * t156 * u[0];
  double t728 = 0.6e1 * t234 * t109;
  double t730 = t29 * C23;
  double t732 = t713 * C13 - 0.4e1 * t730;
  double t734 = 0.6e1 * t732 * t32;
  double t735 = t3 * C23;
  double t736 = 0.12e2 * t735;
  double t740 = t3 * C24;
  double t742 = 0.6e1 * (t234 * t53 / 0.2e1 + t740) * u[9];
  double t745 = 0.6e1 * (-t234 * t53 / 0.2e1 - t740) * u[3];
  double t747 = 0.6e1 * t234 * t67;
  double t749 = 0.12e2 * t732 * t3;
  double t750 = 0.36e2 * t186;
  double t754 = 0.6e1 * t234 * t96;
  double t755 = 0.24e2 * t735;
  double t759 = 0.6e1 * t234 * t122;
  double t761 = 0.6e1 * t717 * t32;
  double t762 = t3 * C22;
  double t763 = 0.12e2 * t762;
  double t767 = 0.6e1 * t234 * t79;
  double t768 = 0.24e2 * t762;
  double t772 = 0.6e1 * t234 * t41;
  double t776 = -0.6e1 * t234 * t67;
  double t779 = (t708 - t719 + t720) * u[7] + t724 - t726 + (t728 - t734 + t736) * u[10] + t742 + t745 + (t747 + t749 - t750) * u[2] + (t754 - t734 + t755) * u[4] + (t759 + t761 - t763) * u[11] + (t767 + t761 - t768) * u[5] + (t772 + t719 - t720) * u[1] + (t776 - t749 + t750) * u[8];
  double t783 = t29 * C21;
  double t785 = t713 * C11 - 0.4e1 * t783;
  double t798 = t29 * C24;
  double t801 = (t713 * C14 - 0.4e1 * t798) * t29;
  double t806 = t717 * t29;
  double t812 = t732 * t29;
  double t839 = t706 + t779 * t1 / 0.2e1 - t785 * u[0] + (t785 * t67 - t734) * u[2] * L + (t785 * t41 - t761) * u[1] * L + t785 * u[6] + (-t785 * t53 / 0.2e1 - t801) * u[3] * L + (t785 * t79 - 0.4e1 * t806) * u[5] * L + (t785 * t96 + 0.4e1 * t812) * u[4] * L + (t785 * t53 / 0.2e1 + t801) * u[9] * L + (t785 * t122 - 0.2e1 * t806) * u[11] * L + (-t785 * t67 + t734) * u[8] * L + (-t785 * t41 + t761) * u[7] * L + (t785 * t109 + 0.2e1 * t812) * u[10] * L;
  double t862 = -t26 * t1 / 0.2e1 - t31 - t44 * u[1] * L - t56 * u[9] * L - t70 * u[2] * L - t83 * u[5] * L - t87 * u[3] * L - t100 * u[4] * L - t104 * u[7] * L - t112 * u[10] * L - t116 * u[8] * L + t119 - t125 * u[11] * L;
  double t866 = -t163 * C12 + t166;
  double t868 = 0.12e2 * t866 * t3;
  double t872 = -t163 * C13 + t182;
  double t874 = 0.6e1 * t872 * t32;
  double t880 = 0.12e2 * t872 * t3;
  double t886 = 0.6e1 * t866 * t32;
  double t895 = (-t158 - t868 - t170) * u[7] - t175 + t177 + (-t179 - t874 - t187) * u[10] - 0.12e2 * t192 * u[9] - 0.12e2 * t195 * u[3] + (-t198 + t880 + t201) * u[2] + (-t205 - t874 - t206) * u[4] + (-t210 + t886 + t214) * u[11] + (-t218 + t886 + t219) * u[5] + (-t223 + t868 + t170) * u[1] + (-t227 - t880 - t201) * u[8];
  double t899 = -t163 * C11 + t235;
  double t913 = (-t163 * C14 + t250) * t29;
  double t918 = t866 * t29;
  double t924 = t872 * t29;
  double t951 = -t153 * t2 / 0.3e1 + t895 * t1 / 0.2e1 - t899 * u[0] + (t899 * t67 - t874) * u[2] * L + (t899 * t41 - t886) * u[1] * L + t899 * u[6] + (-t899 * t53 / 0.2e1 - t913) * u[3] * L + (t899 * t79 - 0.4e1 * t918) * u[5] * L + (t899 * t96 + 0.4e1 * t924) * u[4] * L + (t899 * t53 / 0.2e1 + t913) * u[9] * L + (t899 * t122 - 0.2e1 * t918) * u[11] * L + (-t899 * t67 + t874) * u[8] * L + (-t899 * t41 + t886) * u[7] * L + (t899 * t109 + 0.2e1 * t924) * u[10] * L;
  double t955 = -t321 * C12 + t324;
  double t957 = 0.12e2 * t955 * t3;
  double t961 = -t321 * C13 + t340;
  double t963 = 0.6e1 * t961 * t32;
  double t969 = 0.12e2 * t961 * t3;
  double t975 = 0.6e1 * t955 * t32;
  double t984 = (-t316 - t957 - t328) * u[7] - t333 + t335 + (-t337 - t963 - t345) * u[10] - 0.12e2 * t350 * u[9] - 0.12e2 * t353 * u[3] + (-t356 + t969 + t359) * u[2] + (-t363 - t963 - t364) * u[4] + (-t368 + t975 + t372) * u[11] + (-t376 + t975 + t377) * u[5] + (-t381 + t957 + t328) * u[1] + (-t385 - t969 - t359) * u[8];
  double t988 = -t321 * C11 + t393;
  double t1002 = (-t321 * C14 + t408) * t29;
  double t1007 = t955 * t29;
  double t1013 = t961 * t29;
  double t1040 = -t311 * t2 / 0.3e1 + t984 * t1 / 0.2e1 - t988 * u[0] + (t988 * t67 - t963) * u[2] * L + (t988 * t41 - t975) * u[1] * L + t988 * u[6] + (-t988 * t53 / 0.2e1 - t1002) * u[3] * L + (t988 * t79 - 0.4e1 * t1007) * u[5] * L + (t988 * t96 + 0.4e1 * t1013) * u[4] * L + (t988 * t53 / 0.2e1 + t1002) * u[9] * L + (t988 * t122 - 0.2e1 * t1007) * u[11] * L + (-t988 * t67 + t963) * u[8] * L + (-t988 * t41 + t975) * u[7] * L + (t988 * t109 + 0.2e1 * t1013) * u[10] * L;
  double t1042 = t53 * C13 + t450;
  double t1043 = t1042 * t32;
  double t1046 = t1042 * t3;
  double t1050 = t53 * C12 + t459;
  double t1051 = t1050 * t32;
  double t1054 = t1050 * t3;
  double t1069 = t53 * C11 + t479;
  double t1072 = 0.6e1 * t1051;
  double t1079 = (t53 * C14 + t489) * t29;
  double t1084 = 0.6e1 * t1043;
  double t1089 = t1050 * t29;
  double t1099 = t1042 * t29;
  double t1123 = (-0.6e1 * t1043 * u[4] + 0.12e2 * t1046 * u[2] + 0.6e1 * t1051 * u[5] + 0.12e2 * t1054 * u[1] - 0.12e2 * t1054 * u[7] - 0.6e1 * t1043 * u[10] - 0.12e2 * t1046 * u[8] + 0.6e1 * t1051 * u[11]) * t1 / 0.2e1 - t1069 * u[0] + (t1069 * t41 - t1072) * u[1] * L + (t1069 * t53 / 0.2e1 + t1079) * u[9] * L + (t1069 * t67 - t1084) * u[2] * L + (t1069 * t79 - 0.4e1 * t1089) * u[5] * L + (-t1069 * t53 / 0.2e1 - t1079) * u[3] * L + (t1069 * t96 + 0.4e1 * t1099) * u[4] * L + (-t1069 * t41 + t1072) * u[7] * L + (t1069 * t109 + 0.2e1 * t1099) * u[10] * L + (-t1069 * t67 + t1084) * u[8] * L + t1069 * u[6] + (t1069 * t122 - 0.2e1 * t1089) * u[11] * L;
  double t1126 = -t557 - u[4] / 0.30e2 + t559 + 0.2e1 / 0.15e2 * u[10];
  double t1129 = t1126 * C12 + 0.2e1 * t563;
  double t1131 = 0.12e2 * t1129 * t3;
  double t1136 = t1126 * C13 + 0.2e1 * t578;
  double t1138 = 0.6e1 * t1136 * t32;
  double t1142 = 0.12e2 * t1136 * t3;
  double t1148 = 0.6e1 * t1129 * t32;
  double t1157 = (-t556 - t1131 - t568) * u[7] - t572 + t574 + (-t576 - t1138 - t584) * u[10] + t590 + t593 + (-t595 + t1142 + t598) * u[2] + (-t602 - t1138 - t603) * u[4] + (-t607 + t1148 + t611) * u[11] + (-t615 + t1148 + t616) * u[5] + (-t620 + t1131 + t568) * u[1] + (-t624 - t1142 - t598) * u[8];
  double t1162 = t1126 * C11 + 0.2e1 * t631;
  double t1177 = (t1126 * C14 + 0.2e1 * t646) * t29;
  double t1182 = t1129 * t29;
  double t1188 = t1136 * t29;
  double t1215 = t554 + t1157 * t1 / 0.2e1 - t1162 * u[0] + (t1162 * t67 - t1138) * u[2] * L + (t1162 * t41 - t1148) * u[1] * L + t1162 * u[6] + (-t1162 * t53 / 0.2e1 - t1177) * u[3] * L + (t1162 * t79 - 0.4e1 * t1182) * u[5] * L + (t1162 * t96 + 0.4e1 * t1188) * u[4] * L + (t1162 * t53 / 0.2e1 + t1177) * u[9] * L + (t1162 * t122 - 0.2e1 * t1182) * u[11] * L + (-t1162 * t67 + t1138) * u[8] * L + (-t1162 * t41 + t1148) * u[7] * L + (t1162 * t109 + 0.2e1 * t1188) * u[10] * L;
  double t1218 = t709 - u[5] / 0.30e2 - t711 + 0.2e1 / 0.15e2 * u[11];
  double t1221 = t1218 * C12 - 0.2e1 * t715;
  double t1223 = 0.12e2 * t1221 * t3;
  double t1228 = t1218 * C13 - 0.2e1 * t730;
  double t1230 = 0.6e1 * t1228 * t32;
  double t1234 = 0.12e2 * t1228 * t3;
  double t1240 = 0.6e1 * t1221 * t32;
  double t1249 = (t708 - t1223 + t720) * u[7] + t724 - t726 + (t728 - t1230 + t736) * u[10] + t742 + t745 + (t747 + t1234 - t750) * u[2] + (t754 - t1230 + t755) * u[4] + (t759 + t1240 - t763) * u[11] + (t767 + t1240 - t768) * u[5] + (t772 + t1223 - t720) * u[1] + (t776 - t1234 + t750) * u[8];
  double t1254 = t1218 * C11 - 0.2e1 * t783;
  double t1269 = (t1218 * C14 - 0.2e1 * t798) * t29;
  double t1274 = t1221 * t29;
  double t1280 = t1228 * t29;
  double t1307 = t706 + t1249 * t1 / 0.2e1 - t1254 * u[0] + (t1254 * t67 - t1230) * u[2] * L + (t1254 * t41 - t1240) * u[1] * L + t1254 * u[6] + (-t1254 * t53 / 0.2e1 - t1269) * u[3] * L + (t1254 * t79 - 0.4e1 * t1274) * u[5] * L + (t1254 * t96 + 0.4e1 * t1280) * u[4] * L + (t1254 * t53 / 0.2e1 + t1269) * u[9] * L + (t1254 * t122 - 0.2e1 * t1274) * u[11] * L + (-t1254 * t67 + t1230) * u[8] * L + (-t1254 * t41 + t1240) * u[7] * L + (t1254 * t109 + 0.2e1 * t1280) * u[10] * L;

  // Final vector.

  gl[0] = t128;
  gl[1] = t290;
  gl[2] = t448;
  gl[3] = t535;
  gl[4] = t687;
  gl[5] = t839;
  gl[6] = t862;
  gl[7] = t951;
  gl[8] = t1040;
  gl[9] = t1123;
  gl[10] = t1215;
  gl[11] = t1307;
}

// ================================= LocStiffMat ===========================
//
// This method computes the material tangent stiffness matrix considering
// the existing couplings.
//
//   L  - element length                                               (in)
//   Ip - polar moment of inertia of the cross section                 (in)
//   A  - area of the cross section                                    (in)
//   u  - local displacement vector (12)                               (in)
//   kl - tangent stiffness matrix (12x12)                            (out)
//
void cFrame3DTL :: LocStiffMat(double L, double Ip, double A, cVector &u, 
                               cMatrix &kfull)
{
  // Classical formulation (diagonal section matrix).

  if (Section->GetType( ) != SEC_BAR_B3DCOUPLED)
  {
    cout << "Error: Frame3DTL with uncoupled section not implemented yet!";
    exit(0); 
  }

  // Section stiffness - fully coupled 4x4 matrix.
  
  cMatrix Cb(4, 4);
  Section->SecStiffness3D(Cb);
  double C11 = Cb[0][0]; double C12 = Cb[0][1]; double C13 = Cb[0][2]; double C14 = Cb[0][3];
  double C21 = Cb[1][0]; double C22 = Cb[1][1]; double C23 = Cb[1][2]; double C24 = Cb[1][3];
  double C31 = Cb[2][0]; double C32 = Cb[2][1]; double C33 = Cb[2][2]; double C34 = Cb[2][3];
  double C41 = Cb[3][0]; double C42 = Cb[3][1]; double C43 = Cb[3][2]; double C44 = Cb[3][3];
  
  // Auxiliary variables (Maple code).  

  double t1 = 0.1e1 / L;
  double t2 = t1 * C11;
  double t3 = L * L;
  double t4 = 0.1e1 / t3;
  double t13 = 0.6e1 / 0.5e1 * u[1] * t4 + u[5] * t1 / 0.10e2 - 0.6e1 / 0.5e1 * u[7] * t4 + u[11] * t1 / 0.10e2;
  double t14 = t13 * C11;
  double t23 = 0.6e1 / 0.5e1 * u[2] * t4 - u[4] * t1 / 0.10e2 - 0.6e1 / 0.5e1 * u[8] * t4 - u[10] * t1 / 0.10e2;
  double t24 = t23 * C11;
  double t27 = Ip / A;
  double t31 = u[3] * t4 * t27 - u[9] * t4 * t27;
  double t32 = t31 * C11;
  double t33 = t1 * C14;
  double t34 = -t32 + t33;
  double t35 = t1 * C13;
  double t37 = u[2] * t1 / 0.10e2;
  double t40 = u[8] * t1 / 0.10e2;
  double t42 = -t37 + 0.2e1 / 0.15e2 * u[4] + t40 - u[10] / 0.30e2;
  double t43 = t42 * C11;
  double t44 = -t35 - t43;
  double t45 = t1 * C12;
  double t47 = u[1] * t1 / 0.10e2;
  double t50 = u[7] * t1 / 0.10e2;
  double t52 = t47 + 0.2e1 / 0.15e2 * u[5] - t50 - u[11] / 0.30e2;
  double t53 = t52 * C11;
  double t54 = t45 - t53;
  double t55 = -t13 * C11;
  double t56 = -t23 * C11;
  double t57 = -t31 * C11;
  double t58 = -t57 - t33;
  double t61 = -t37 - u[4] / 0.30e2 + t40 + 0.2e1 / 0.15e2 * u[10];
  double t62 = t61 * C11;
  double t63 = t35 - t62;
  double t66 = t47 - u[5] / 0.30e2 - t50 + 0.2e1 / 0.15e2 * u[11];
  double t67 = t66 * C11;
  double t68 = -t45 - t67;
  double t70 = 0.1e1 / t3 / L;
  double t71 = t70 * C22;
  double t72 = 0.48e2 * t71;
  double t73 = t70 * C21;
  double t75 = 0.12e2 * t73 * t13;
  double t76 = t13 * C12;
  double t77 = t4 * C22;
  double t78 = 0.6e1 * t77;
  double t79 = t76 - t78;
  double t81 = 0.12e2 * t79 * t70;
  double t82 = t3 * t3;
  double t84 = 0.1e1 / t82 / L;
  double t86 = 0.72e2 * t84 * C22;
  double t90 = t4 * C21;
  double t91 = 0.6e1 * t90;
  double t92 = t14 - t91;
  double t96 = 0.6e1 * t79 * t1;
  double t98 = t70 * C23;
  double t99 = 0.48e2 * t98;
  double t101 = 0.12e2 * t73 * t23;
  double t102 = t13 * C13;
  double t103 = t4 * C23;
  double t104 = 0.6e1 * t103;
  double t105 = t102 - t104;
  double t107 = 0.12e2 * t105 * t70;
  double t109 = 0.72e2 * t84 * C23;
  double t116 = 0.6e1 * t105 * t1;
  double t119 = 0.1e1 / t82;
  double t120 = t119 * C24;
  double t121 = t73 * t31 - t120;
  double t126 = t13 * C14;
  double t128 = 0.6e1 * t4 * C24;
  double t130 = 0.48e2 * t103;
  double t132 = 0.12e2 * t73 * t42;
  double t134 = 0.6e1 * t105 * t4;
  double t135 = t119 * C23;
  double t136 = 0.48e2 * t135;
  double t144 = 0.48e2 * t77;
  double t146 = 0.12e2 * t73 * t52;
  double t148 = 0.6e1 * t79 * t4;
  double t149 = t119 * C22;
  double t150 = 0.48e2 * t149;
  double t159 = -0.12e2 * t73 * t13;
  double t167 = -0.12e2 * t73 * t23;
  double t175 = -t73 * t31 + t120;
  double t181 = 0.36e2 * t103;
  double t183 = 0.12e2 * t73 * t61;
  double t184 = 0.24e2 * t135;
  double t192 = 0.36e2 * t77;
  double t194 = 0.12e2 * t73 * t66;
  double t195 = 0.24e2 * t149;
  double t203 = t70 * C32;
  double t204 = 0.48e2 * t203;
  double t205 = t70 * C31;
  double t207 = 0.12e2 * t205 * t13;
  double t208 = t23 * C12;
  double t209 = t4 * C32;
  double t210 = 0.6e1 * t209;
  double t211 = t208 - t210;
  double t213 = 0.12e2 * t211 * t70;
  double t215 = 0.72e2 * t84 * C32;
  double t219 = t4 * C31;
  double t220 = 0.6e1 * t219;
  double t221 = t24 - t220;
  double t225 = 0.6e1 * t211 * t1;
  double t227 = t70 * C33;
  double t228 = 0.48e2 * t227;
  double t230 = 0.12e2 * t205 * t23;
  double t231 = t23 * C13;
  double t232 = t4 * C33;
  double t233 = 0.6e1 * t232;
  double t234 = t231 - t233;
  double t236 = 0.12e2 * t234 * t70;
  double t238 = 0.72e2 * t84 * C33;
  double t245 = 0.6e1 * t234 * t1;
  double t248 = t119 * C34;
  double t249 = t205 * t31 - t248;
  double t254 = t23 * C14;
  double t256 = 0.6e1 * t4 * C34;
  double t258 = 0.48e2 * t232;
  double t260 = 0.12e2 * t205 * t42;
  double t262 = 0.6e1 * t234 * t4;
  double t263 = t119 * C33;
  double t264 = 0.48e2 * t263;
  double t272 = 0.48e2 * t209;
  double t274 = 0.12e2 * t205 * t52;
  double t276 = 0.6e1 * t211 * t4;
  double t277 = t119 * C32;
  double t278 = 0.48e2 * t277;
  double t287 = -0.12e2 * t205 * t13;
  double t295 = -0.12e2 * t205 * t23;
  double t303 = -t205 * t31 + t248;
  double t309 = 0.36e2 * t232;
  double t311 = 0.12e2 * t205 * t61;
  double t312 = 0.24e2 * t263;
  double t320 = 0.36e2 * t209;
  double t322 = 0.12e2 * t205 * t66;
  double t323 = 0.24e2 * t277;
  double t331 = t1 * C41;
  double t332 = -t32 + t331;
  double t339 = t31 * C14;
  double t340 = t1 * C44;
  double t342 = t31 * C13;
  double t343 = t1 * C43;
  double t347 = t31 * C12;
  double t348 = t1 * C42;
  double t365 = t1 * C31;
  double t366 = -t365 - t43;
  double t367 = 0.24e2 * t209;
  double t369 = 0.6e1 * t219 * t13;
  double t370 = t42 * C12;
  double t371 = t1 * C32;
  double t373 = t370 + 0.4e1 * t371;
  double t375 = 0.12e2 * t373 * t70;
  double t376 = 0.36e2 * t277;
  double t381 = t43 + 0.4e1 * t365;
  double t385 = 0.6e1 * t373 * t1;
  double t387 = 0.24e2 * t232;
  double t389 = 0.6e1 * t219 * t23;
  double t390 = t42 * C13;
  double t391 = t1 * C33;
  double t393 = t390 + 0.4e1 * t391;
  double t395 = 0.12e2 * t393 * t70;
  double t396 = 0.36e2 * t263;
  double t403 = 0.6e1 * t393 * t1;
  double t406 = t70 * C34;
  double t409 = 0.3e1 * (-t219 * t31 + t406) * t3;
  double t412 = t42 * C14;
  double t413 = t1 * C34;
  double t414 = 0.4e1 * t413;
  double t418 = 0.6e1 * t219 * t42;
  double t420 = 0.6e1 * t393 * t4;
  double t421 = 0.24e2 * t227;
  double t431 = 0.6e1 * t219 * t52;
  double t433 = 0.6e1 * t373 * t4;
  double t434 = 0.24e2 * t203;
  double t443 = -0.6e1 * t219 * t13;
  double t451 = -0.6e1 * t219 * t23;
  double t461 = 0.3e1 * (t219 * t31 - t406) * t3;
  double t465 = 0.20e2 * t391;
  double t467 = 0.6e1 * t219 * t61;
  double t468 = 0.12e2 * t227;
  double t476 = 0.20e2 * t371;
  double t478 = 0.6e1 * t219 * t66;
  double t479 = 0.12e2 * t203;
  double t487 = t1 * C21;
  double t488 = t487 - t53;
  double t489 = 0.24e2 * t77;
  double t491 = 0.6e1 * t90 * t13;
  double t492 = t52 * C12;
  double t493 = t1 * C22;
  double t495 = t492 - 0.4e1 * t493;
  double t497 = 0.12e2 * t495 * t70;
  double t498 = 0.36e2 * t149;
  double t503 = t53 - 0.4e1 * t487;
  double t507 = 0.6e1 * t495 * t1;
  double t509 = 0.24e2 * t103;
  double t511 = 0.6e1 * t90 * t23;
  double t512 = t52 * C13;
  double t513 = t1 * C23;
  double t515 = t512 - 0.4e1 * t513;
  double t517 = 0.12e2 * t515 * t70;
  double t518 = 0.36e2 * t135;
  double t525 = 0.6e1 * t515 * t1;
  double t528 = t70 * C24;
  double t531 = 0.3e1 * (t90 * t31 - t528) * t3;
  double t534 = t52 * C14;
  double t535 = t1 * C24;
  double t536 = 0.4e1 * t535;
  double t540 = 0.6e1 * t90 * t42;
  double t542 = 0.6e1 * t515 * t4;
  double t543 = 0.24e2 * t98;
  double t553 = 0.6e1 * t90 * t52;
  double t555 = 0.6e1 * t495 * t4;
  double t556 = 0.24e2 * t71;
  double t565 = -0.6e1 * t90 * t13;
  double t573 = -0.6e1 * t90 * t23;
  double t583 = 0.3e1 * (-t90 * t31 + t528) * t3;
  double t587 = 0.20e2 * t513;
  double t589 = 0.6e1 * t90 * t61;
  double t590 = 0.12e2 * t98;
  double t598 = 0.20e2 * t493;
  double t600 = 0.6e1 * t90 * t66;
  double t601 = 0.12e2 * t71;
  double t609 = -t13 * C12;
  double t610 = t609 + t78;
  double t612 = 0.12e2 * t610 * t70;
  double t616 = t55 + t91;
  double t620 = 0.6e1 * t610 * t1;
  double t622 = -t13 * C13;
  double t623 = t622 + t104;
  double t625 = 0.12e2 * t623 * t70;
  double t632 = 0.6e1 * t623 * t1;
  double t638 = -t13 * C14;
  double t641 = 0.6e1 * t623 * t4;
  double t650 = 0.6e1 * t610 * t4;
  double t689 = -t23 * C12;
  double t690 = t689 + t210;
  double t692 = 0.12e2 * t690 * t70;
  double t696 = t56 + t220;
  double t700 = 0.6e1 * t690 * t1;
  double t702 = -t23 * C13;
  double t703 = t702 + t233;
  double t705 = 0.12e2 * t703 * t70;
  double t712 = 0.6e1 * t703 * t1;
  double t718 = -t23 * C14;
  double t721 = 0.6e1 * t703 * t4;
  double t730 = 0.6e1 * t690 * t4;
  double t769 = -t57 - t331;
  double t776 = -t31 * C14;
  double t778 = -t31 * C13;
  double t782 = -t31 * C12;
  double t799 = t365 - t62;
  double t800 = t61 * C12;
  double t802 = t800 + 0.2e1 * t371;
  double t804 = 0.12e2 * t802 * t70;
  double t809 = t62 + 0.2e1 * t365;
  double t813 = 0.6e1 * t802 * t1;
  double t815 = t61 * C13;
  double t817 = t815 + 0.2e1 * t391;
  double t819 = 0.12e2 * t817 * t70;
  double t826 = 0.6e1 * t817 * t1;
  double t830 = t61 * C14;
  double t831 = 0.2e1 * t413;
  double t834 = 0.6e1 * t817 * t4;
  double t843 = 0.6e1 * t802 * t4;
  double t882 = -t487 - t67;
  double t883 = t66 * C12;
  double t885 = t883 - 0.2e1 * t493;
  double t887 = 0.12e2 * t885 * t70;
  double t892 = t67 - 0.2e1 * t487;
  double t896 = 0.6e1 * t885 * t1;
  double t898 = t66 * C13;
  double t900 = t898 - 0.2e1 * t513;
  double t902 = 0.12e2 * t900 * t70;
  double t909 = 0.6e1 * t900 * t1;
  double t913 = t66 * C14;
  double t914 = 0.2e1 * t535;
  double t917 = 0.6e1 * t900 * t4;
  double t926 = 0.6e1 * t885 * t4;
  
  // Assembly the [kfull] matrix (Maple code).
  
  kfull[0][0] = t2;
  kfull[0][1] = -t14;
  kfull[0][2] = -t24;
  kfull[0][3] = t34;
  kfull[0][4] = t44;
  kfull[0][5] = t54;
  kfull[0][6] = -t2;
  kfull[0][7] = -t55;
  kfull[0][8] = -t56;
  kfull[0][9] = t58;
  kfull[0][10] = t63;
  kfull[0][11] = t68;
  kfull[1][0] = -t14;
  kfull[1][1] = t72 + (t75 + t81 - t86) * t3 / 0.2e1 + t92 * t13 * L - t96;
  kfull[1][2] = t99 + (t101 + t107 - t109) * t3 / 0.2e1 + t92 * t23 * L - t116;
  kfull[1][3] = 0.6e1 * t121 * t3 + t92 * t31 * L - t126 + t128;
  kfull[1][4] = -t130 + (t132 - t134 + t136) * t3 / 0.2e1 + t92 * t42 * L + 0.4e1 * t102;
  kfull[1][5] = t144 + (t146 + t148 - t150) * t3 / 0.2e1 + t92 * t52 * L - 0.4e1 * t76;
  kfull[1][6] = t14;
  kfull[1][7] = -t72 + (t159 - t81 + t86) * t3 / 0.2e1 - t92 * t13 * L + t96;
  kfull[1][8] = -t99 + (t167 - t107 + t109) * t3 / 0.2e1 - t92 * t23 * L + t116;
  kfull[1][9] = 0.6e1 * t175 * t3 - t92 * t31 * L + t126 - t128;
  kfull[1][10] = -t181 + (t183 - t134 + t184) * t3 / 0.2e1 + t92 * t61 * L + 0.2e1 * t102;
  kfull[1][11] = t192 + (t194 + t148 - t195) * t3 / 0.2e1 + t92 * t66 * L - 0.2e1 * t76;
  kfull[2][0] = -t24;
  kfull[2][1] = t204 + (t207 + t213 - t215) * t3 / 0.2e1 + t221 * t13 * L - t225;
  kfull[2][2] = t228 + (t230 + t236 - t238) * t3 / 0.2e1 + t221 * t23 * L - t245;
  kfull[2][3] = 0.6e1 * t249 * t3 + t221 * t31 * L - t254 + t256;
  kfull[2][4] = -t258 + (t260 - t262 + t264) * t3 / 0.2e1 + t221 * t42 * L + 0.4e1 * t231;
  kfull[2][5] = t272 + (t274 + t276 - t278) * t3 / 0.2e1 + t221 * t52 * L - 0.4e1 * t208;
  kfull[2][6] = t24;
  kfull[2][7] = -t204 + (t287 - t213 + t215) * t3 / 0.2e1 - t221 * t13 * L + t225;
  kfull[2][8] = -t228 + (t295 - t236 + t238) * t3 / 0.2e1 - t221 * t23 * L + t245;
  kfull[2][9] = 0.6e1 * t303 * t3 - t221 * t31 * L + t254 - t256;
  kfull[2][10] = -t309 + (t311 - t262 + t312) * t3 / 0.2e1 + t221 * t61 * L + 0.2e1 * t231;
  kfull[2][11] = t320 + (t322 + t276 - t323) * t3 / 0.2e1 + t221 * t66 * L - 0.2e1 * t208;
  kfull[3][0] = t332;
  kfull[3][1] = -t332 * t13 * L;
  kfull[3][2] = -t332 * t23 * L;
  kfull[3][3] = -t332 * t31 * L - t339 + t340;
  kfull[3][4] = t342 - t343 - t332 * t42 * L;
  kfull[3][5] = -t347 + t348 - t332 * t52 * L;
  kfull[3][6] = -t332;
  kfull[3][7] = t332 * t13 * L;
  kfull[3][8] = t332 * t23 * L;
  kfull[3][9] = t332 * t31 * L + t339 - t340;
  kfull[3][10] = -t342 + t343 - t332 * t61 * L;
  kfull[3][11] = t347 - t348 - t332 * t66 * L;
  kfull[4][0] = t366;
  kfull[4][1] = -t367 + (-t369 + t375 + t376) * t3 / 0.2e1 + t381 * t13 * L - t385;
  kfull[4][2] = -t387 + (-t389 + t395 + t396) * t3 / 0.2e1 + t381 * t23 * L - t403;
  kfull[4][3] = t409 + t381 * t31 * L - t412 - t414;
  kfull[4][4] = 0.28e2 * t391 + (-t418 - t420 - t421) * t3 / 0.2e1 + t381 * t42 * L + 0.4e1 * t390;
  kfull[4][5] = -0.28e2 * t371 + (-t431 + t433 + t434) * t3 / 0.2e1 + t381 * t52 * L - 0.4e1 * t370;
  kfull[4][6] = -t366;
  kfull[4][7] = t367 + (-t443 - t375 - t376) * t3 / 0.2e1 - t381 * t13 * L + t385;
  kfull[4][8] = t387 + (-t451 - t395 - t396) * t3 / 0.2e1 - t381 * t23 * L + t403;
  kfull[4][9] = t461 - t381 * t31 * L + t412 + t414;
  kfull[4][10] = t465 + (-t467 - t420 - t468) * t3 / 0.2e1 + t381 * t61 * L + 0.2e1 * t390;
  kfull[4][11] = -t476 + (-t478 + t433 + t479) * t3 / 0.2e1 + t381 * t66 * L - 0.2e1 * t370;
  kfull[5][0] = t488;
  kfull[5][1] = t489 + (t491 + t497 - t498) * t3 / 0.2e1 + t503 * t13 * L - t507;
  kfull[5][2] = t509 + (t511 + t517 - t518) * t3 / 0.2e1 + t503 * t23 * L - t525;
  kfull[5][3] = t531 + t503 * t31 * L - t534 + t536;
  kfull[5][4] = -0.28e2 * t513 + (t540 - t542 + t543) * t3 / 0.2e1 + t503 * t42 * L + 0.4e1 * t512;
  kfull[5][5] = 0.28e2 * t493 + (t553 + t555 - t556) * t3 / 0.2e1 + t503 * t52 * L - 0.4e1 * t492;
  kfull[5][6] = -t488;
  kfull[5][7] = -t489 + (t565 - t497 + t498) * t3 / 0.2e1 - t503 * t13 * L + t507;
  kfull[5][8] = -t509 + (t573 - t517 + t518) * t3 / 0.2e1 - t503 * t23 * L + t525;
  kfull[5][9] = t583 - t503 * t31 * L + t534 - t536;
  kfull[5][10] = -t587 + (t589 - t542 + t590) * t3 / 0.2e1 + t503 * t61 * L + 0.2e1 * t512;
  kfull[5][11] = t598 + (t600 + t555 - t601) * t3 / 0.2e1 + t503 * t66 * L - 0.2e1 * t492;
  kfull[6][0] = -t2;
  kfull[6][1] = t14;
  kfull[6][2] = t24;
  kfull[6][3] = -t34;
  kfull[6][4] = -t44;
  kfull[6][5] = -t54;
  kfull[6][6] = t2;
  kfull[6][7] = t55;
  kfull[6][8] = t56;
  kfull[6][9] = -t58;
  kfull[6][10] = -t63;
  kfull[6][11] = -t68;
  kfull[7][0] = -t55;
  kfull[7][1] = -t72 + (-t75 + t612 + t86) * t3 / 0.2e1 + t616 * t13 * L - t620;
  kfull[7][2] = -t99 + (-t101 + t625 + t109) * t3 / 0.2e1 + t616 * t23 * L - t632;
  kfull[7][3] = -0.6e1 * t121 * t3 + t616 * t31 * L - t638 - t128;
  kfull[7][4] = t130 + (-t132 - t641 - t136) * t3 / 0.2e1 + t616 * t42 * L + 0.4e1 * t622;
  kfull[7][5] = -t144 + (-t146 + t650 + t150) * t3 / 0.2e1 + t616 * t52 * L - 0.4e1 * t609;
  kfull[7][6] = t55;
  kfull[7][7] = t72 + (-t159 - t612 - t86) * t3 / 0.2e1 - t616 * t13 * L + t620;
  kfull[7][8] = t99 + (-t167 - t625 - t109) * t3 / 0.2e1 - t616 * t23 * L + t632;
  kfull[7][9] = -0.6e1 * t175 * t3 - t616 * t31 * L + t638 + t128;
  kfull[7][10] = t181 + (-t183 - t641 - t184) * t3 / 0.2e1 + t616 * t61 * L + 0.2e1 * t622;
  kfull[7][11] = -t192 + (-t194 + t650 + t195) * t3 / 0.2e1 + t616 * t66 * L - 0.2e1 * t609;
  kfull[8][0] = -t56;
  kfull[8][1] = -t204 + (-t207 + t692 + t215) * t3 / 0.2e1 + t696 * t13 * L - t700;
  kfull[8][2] = -t228 + (-t230 + t705 + t238) * t3 / 0.2e1 + t696 * t23 * L - t712;
  kfull[8][3] = -0.6e1 * t249 * t3 + t696 * t31 * L - t718 - t256;
  kfull[8][4] = t258 + (-t260 - t721 - t264) * t3 / 0.2e1 + t696 * t42 * L + 0.4e1 * t702;
  kfull[8][5] = -t272 + (-t274 + t730 + t278) * t3 / 0.2e1 + t696 * t52 * L - 0.4e1 * t689;
  kfull[8][6] = t56;
  kfull[8][7] = t204 + (-t287 - t692 - t215) * t3 / 0.2e1 - t696 * t13 * L + t700;
  kfull[8][8] = t228 + (-t295 - t705 - t238) * t3 / 0.2e1 - t696 * t23 * L + t712;
  kfull[8][9] = -0.6e1 * t303 * t3 - t696 * t31 * L + t718 + t256;
  kfull[8][10] = t309 + (-t311 - t721 - t312) * t3 / 0.2e1 + t696 * t61 * L + 0.2e1 * t702;
  kfull[8][11] = -t320 + (-t322 + t730 + t323) * t3 / 0.2e1 + t696 * t66 * L - 0.2e1 * t689;
  kfull[9][0] = t769;
  kfull[9][1] = -t769 * t13 * L;
  kfull[9][2] = -t769 * t23 * L;
  kfull[9][3] = -t769 * t31 * L - t776 - t340;
  kfull[9][4] = t778 + t343 - t769 * t42 * L;
  kfull[9][5] = -t782 - t348 - t769 * t52 * L;
  kfull[9][6] = -t769;
  kfull[9][7] = t769 * t13 * L;
  kfull[9][8] = t769 * t23 * L;
  kfull[9][9] = t769 * t31 * L + t776 + t340;
  kfull[9][10] = -t778 - t343 - t769 * t61 * L;
  kfull[9][11] = t782 + t348 - t769 * t66 * L;
  kfull[10][0] = t799;
  kfull[10][1] = -t367 + (-t369 + t804 + t376) * t3 / 0.2e1 + t809 * t13 * L - t813;
  kfull[10][2] = -t387 + (-t389 + t819 + t396) * t3 / 0.2e1 + t809 * t23 * L - t826;
  kfull[10][3] = t409 + t809 * t31 * L - t830 - t831;
  kfull[10][4] = t465 + (-t418 - t834 - t421) * t3 / 0.2e1 + t809 * t42 * L + 0.4e1 * t815;
  kfull[10][5] = -t476 + (-t431 + t843 + t434) * t3 / 0.2e1 + t809 * t52 * L - 0.4e1 * t800;
  kfull[10][6] = -t799;
  kfull[10][7] = t367 + (-t443 - t804 - t376) * t3 / 0.2e1 - t809 * t13 * L + t813;
  kfull[10][8] = t387 + (-t451 - t819 - t396) * t3 / 0.2e1 - t809 * t23 * L + t826;
  kfull[10][9] = t461 - t809 * t31 * L + t830 + t831;
  kfull[10][10] = 0.16e2 * t391 + (-t467 - t834 - t468) * t3 / 0.2e1 + t809 * t61 * L + 0.2e1 * t815;
  kfull[10][11] = -0.16e2 * t371 + (-t478 + t843 + t479) * t3 / 0.2e1 + t809 * t66 * L - 0.2e1 * t800;
  kfull[11][0] = t882;
  kfull[11][1] = t489 + (t491 + t887 - t498) * t3 / 0.2e1 + t892 * t13 * L - t896;
  kfull[11][2] = t509 + (t511 + t902 - t518) * t3 / 0.2e1 + t892 * t23 * L - t909;
  kfull[11][3] = t531 + t892 * t31 * L - t913 + t914;
  kfull[11][4] = -t587 + (t540 - t917 + t543) * t3 / 0.2e1 + t892 * t42 * L + 0.4e1 * t898;
  kfull[11][5] = t598 + (t553 + t926 - t556) * t3 / 0.2e1 + t892 * t52 * L - 0.4e1 * t883;
  kfull[11][6] = -t882;
  kfull[11][7] = -t489 + (t565 - t887 + t498) * t3 / 0.2e1 - t892 * t13 * L + t896;
  kfull[11][8] = -t509 + (t573 - t902 + t518) * t3 / 0.2e1 - t892 * t23 * L + t909;
  kfull[11][9] = t583 - t892 * t31 * L + t913 - t914;
  kfull[11][10] = -0.16e2 * t513 + (t589 - t917 + t590) * t3 / 0.2e1 + t892 * t61 * L + 0.2e1 * t898;
  kfull[11][11] = 0.16e2 * t493 + (t600 + t926 - t601) * t3 / 0.2e1 + t892 * t66 * L - 0.2e1 * t883;
}


// =================================== AMat ================================
//
// This method computes the auxiliary matrix [At] due to transverse
// displacements v and w, and rotation of the cross section effects.
//
//   L  - element length                                               (in)
//   Ip - polar moment of inertia of the cross section                 (in)
//   A  - area of the cross section                                    (in)
//   At - auxiliary matrix (12x12)                                    (out)
//
void cFrame3DTL :: AMat(double L, double Ip, double A, cMatrix &At)
{
  // Auxiliary variables (Maple code).
  
  double t1 = L * L;
  double t2 = 0.1e1 / t1;
  double t3 = 0.6e1 / 0.5e1 * t2;
  double t5 = 0.1e1 / L / 0.10e2;
  double t8 = t2 * Ip / A;
  
  // Assembly the [A] matrix (Maple code).
  
  At[0][0] = 0.0;
  At[0][1] = 0.0;
  At[0][2] = 0.0;
  At[0][3] = 0.0;
  At[0][4] = 0.0;
  At[0][5] = 0.0;
  At[0][6] = 0.0;
  At[0][7] = 0.0;
  At[0][8] = 0.0;
  At[0][9] = 0.0;
  At[0][10] = 0.0;
  At[0][11] = 0.0;
  At[1][0] = 0.0;
  At[1][1] = t3;
  At[1][2] = 0.0;
  At[1][3] = 0.0;
  At[1][4] = 0.0;
  At[1][5] = t5;
  At[1][6] = 0.0;
  At[1][7] = -t3;
  At[1][8] = 0.0;
  At[1][9] = 0.0;
  At[1][10] = 0.0;
  At[1][11] = t5;
  At[2][0] = 0.0;
  At[2][1] = 0.0;
  At[2][2] = t3;
  At[2][3] = 0.0;
  At[2][4] = -t5;
  At[2][5] = 0.0;
  At[2][6] = 0.0;
  At[2][7] = 0.0;
  At[2][8] = -t3;
  At[2][9] = 0.0;
  At[2][10] = -t5;
  At[2][11] = 0.0;
  At[3][0] = 0.0;
  At[3][1] = 0.0;
  At[3][2] = 0.0;
  At[3][3] = t8;
  At[3][4] = 0.0;
  At[3][5] = 0.0;
  At[3][6] = 0.0;
  At[3][7] = 0.0;
  At[3][8] = 0.0;
  At[3][9] = -t8;
  At[3][10] = 0.0;
  At[3][11] = 0.0;
  At[4][0] = 0.0;
  At[4][1] = 0.0;
  At[4][2] = -t5;
  At[4][3] = 0.0;
  At[4][4] = 0.2e1 / 0.15e2;
  At[4][5] = 0.0;
  At[4][6] = 0.0;
  At[4][7] = 0.0;
  At[4][8] = t5;
  At[4][9] = 0.0;
  At[4][10] = -0.1e1 / 0.30e2;
  At[4][11] = 0.0;
  At[5][0] = 0.0;
  At[5][1] = t5;
  At[5][2] = 0.0;
  At[5][3] = 0.0;
  At[5][4] = 0.0;
  At[5][5] = 0.2e1 / 0.15e2;
  At[5][6] = 0.0;
  At[5][7] = -t5;
  At[5][8] = 0.0;
  At[5][9] = 0.0;
  At[5][10] = 0.0;
  At[5][11] = -0.1e1 / 0.30e2;
  At[6][0] = 0.0;
  At[6][1] = 0.0;
  At[6][2] = 0.0;
  At[6][3] = 0.0;
  At[6][4] = 0.0;
  At[6][5] = 0.0;
  At[6][6] = 0.0;
  At[6][7] = 0.0;
  At[6][8] = 0.0;
  At[6][9] = 0.0;
  At[6][10] = 0.0;
  At[6][11] = 0.0;
  At[7][0] = 0.0;
  At[7][1] = -t3;
  At[7][2] = 0.0;
  At[7][3] = 0.0;
  At[7][4] = 0.0;
  At[7][5] = -t5;
  At[7][6] = 0.0;
  At[7][7] = t3;
  At[7][8] = 0.0;
  At[7][9] = 0.0;
  At[7][10] = 0.0;
  At[7][11] = -t5;
  At[8][0] = 0.0;
  At[8][1] = 0.0;
  At[8][2] = -t3;
  At[8][3] = 0.0;
  At[8][4] = t5;
  At[8][5] = 0.0;
  At[8][6] = 0.0;
  At[8][7] = 0.0;
  At[8][8] = t3;
  At[8][9] = 0.0;
  At[8][10] = t5;
  At[8][11] = 0.0;
  At[9][0] = 0.0;
  At[9][1] = 0.0;
  At[9][2] = 0.0;
  At[9][3] = -t8;
  At[9][4] = 0.0;
  At[9][5] = 0.0;
  At[9][6] = 0.0;
  At[9][7] = 0.0;
  At[9][8] = 0.0;
  At[9][9] = t8;
  At[9][10] = 0.0;
  At[9][11] = 0.0;
  At[10][0] = 0.0;
  At[10][1] = 0.0;
  At[10][2] = -t5;
  At[10][3] = 0.0;
  At[10][4] = -0.1e1 / 0.30e2;
  At[10][5] = 0.0;
  At[10][6] = 0.0;
  At[10][7] = 0.0;
  At[10][8] = t5;
  At[10][9] = 0.0;
  At[10][10] = 0.2e1 / 0.15e2;
  At[10][11] = 0.0;
  At[11][0] = 0.0;
  At[11][1] = t5;
  At[11][2] = 0.0;
  At[11][3] = 0.0;
  At[11][4] = 0.0;
  At[11][5] = -0.1e1 / 0.30e2;
  At[11][6] = 0.0;
  At[11][7] = -t5;
  At[11][8] = 0.0;
  At[11][9] = 0.0;
  At[11][10] = 0.0;
  At[11][11] = 0.2e1 / 0.15e2;
}

// ================================ NormalForce ============================
//
// This method compute the normal force considering the existing couplings.
//
//   L  - element length                                               (in)
//   Ip - polar moment of inertia of the cross section                 (in)
//   A  - area of the cross section                                    (in)
//   u  - local displacement vector (12)                               (in)
//   N  - normal force                                                (out)
//
double cFrame3DTL :: NormalForce(double L, double Ip, double A, cVector &u)
{
  // Classical formulation (diagonal section matrix).

  if (Section->GetType( ) != SEC_BAR_B3DCOUPLED)
  {
    cout << "Error: Frame3DTL with uncoupled section not implemented yet!";
    exit(0); 
  }

  // First line of the fully coupled section stiffness matrix.
  
  cMatrix Cb(4, 4);
  Section->SecStiffness3D(Cb);
  double C11 = Cb[0][0]; double C12 = Cb[0][1]; double C13 = Cb[0][2]; double C14 = Cb[0][3];
  
  // Auxiliary variables (Maple code).

  double t1 = L * L;
  double t3 = 0.1e1 / t1 / L;
  double t6 = 0.1e1 / t1;
  double t28 = 0.1e1 / L;
  double t30 = u[1] * t6;
  double t32 = u[5] * t28;
  double t34 = u[7] * t6;
  double t36 = u[11] * t28;
  double t38 = 0.3e1 / 0.5e1 * t30 + t32 / 0.20e2 - 0.3e1 / 0.5e1 * t34 + t36 / 0.20e2;
  double t40 = u[2] * t6;
  double t42 = u[4] * t28;
  double t44 = u[8] * t6;
  double t46 = u[10] * t28;
  double t48 = 0.3e1 / 0.5e1 * t40 - t42 / 0.20e2 - 0.3e1 / 0.5e1 * t44 - t46 / 0.20e2;
  double t52 = Ip / A;
  double t56 = u[3] * t6 * t52 - u[9] * t6 * t52;
  double t59 = u[2] * t28 / 0.20e2;
  double t62 = u[8] * t28 / 0.20e2;
  double t67 = u[1] * t28 / 0.20e2;
  double t70 = u[7] * t28 / 0.20e2;
  double t86 = -t28 * u[0] + t38 * u[1] + t48 * u[2] + t56 * u[3] / 0.2e1 + (-t59 + u[4] / 0.15e2 + t62 - u[10] / 0.60e2) * u[4] + (t67 + u[5] / 0.15e2 - t70 - u[11] / 0.60e2) * u[5] + t28 * u[6] - t38 * u[7] - t48 * u[8] - t56 * u[9] / 0.2e1 + (-t59 - u[4] / 0.60e2 + t62 + u[10] / 0.15e2) * u[10] + (t67 - u[5] / 0.60e2 - t70 + u[11] / 0.15e2) * u[11];
  
  // Compute the Nfull.
  
  double Nfull = (C12 * (0.12e2 * t3 * u[1] + 0.6e1 * t6 * u[5] - 0.12e2 * t3 * u[7] + 0.6e1 * t6 * u[11]) + C13 * (0.12e2 * t3 * u[2] - 0.6e1 * t6 * u[4] - 0.12e2 * t3 * u[8] - 0.6e1 * t6 * u[10])) * t1 / 0.2e1 + C11 * t86 * L + C12 * (-0.6e1 * t30 - 0.4e1 * t32 + 0.6e1 * t34 - 0.2e1 * t36) * L + C13 * (-0.6e1 * t40 + 0.4e1 * t42 + 0.6e1 * t44 + 0.2e1 * t46) * L + C14 * (-t28 * u[3] + t28 * u[9]) * L;
  
  return (Nfull);
}

// ======================================================= End of file =====
