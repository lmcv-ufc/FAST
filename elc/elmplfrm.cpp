// -------------------------------------------------------------------------
// elmplfrm.cpp - implementation of plane frame classes.
// -------------------------------------------------------------------------
// Created:      21-Nov-2000     Evandro Parente Junior
//
// Modified:     03-Oct-2001     Evandro Parente Junior
//               Use of vector/matrix classes.
//
// Modified:     17-Oct-2001     Evandro Parente Junior
//               Implementation of the corotational element.
//
// Modified:     04-Dec-2011     Iuri Barcelos Rocha
//               Added the geometric stiffness matrix of the linear element.
//
// Modified:     12-Oct-2012     Luiz Antonio Taumaturgo Mororo
//               Implementation of the Total Lagrangian elements.
//
// Modified:     03-Aug-2014     Evandro Parente Junior
//               Use of SecAnalysis (path dependent materials).
// -------------------------------------------------------------------------

#include <math.h>

#include "elmplfrm.h"
#include "node.h"
#include "material.h"
#include "section.h"
#include "secbar.h"
#include "secanalysis.h"
#include "shpline.h"
#include "intpoint.h"
#include "ctrl.h"
#include "load.h"
#include "vec.h"
#include "mat.h"
#include "gblvar.h"
#include "gbldef.h"

// -------------------------------------------------------------------------
// Auxiliary functions:
//

// =============================== GetAtan =================================

inline double GetAtan(double y, double x)
{
 double at = atan2(y,x);

 if (at < 0) at += 2*PI;

 return(at);
}

// -------------------------------------------------------------------------
// cPlFrame class:
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Public methods:
//

// ============================= ValidSection ==============================

bool cPlFrame :: ValidSection(eSecType type)
{
  switch (type)
  {
    case SEC_BAR_GEN:
    case SEC_BAR_HOM_RECT:
    case SEC_BAR_HOM_CIRCLE:
    case SEC_BAR_HOM_TUBE:
    case SEC_BAR_RC_RECT:
      return(true);
    break;

    default:
    break;
  }

  return(false);
}

// =============================== cPlFrame ================================

cPlFrame :: cPlFrame(int id, cSection *sec) : cElement(id)
{
  Type    = PLFRAME;
  Section = sec;
  Shape   = new cShapeBar( );
}

// ============================== ~cPlFrame ================================

cPlFrame :: ~cPlFrame(void)
{
}

// ============================== GetActDir ================================

void cPlFrame :: GetActDir(int *dir)
{
  dir[0] = 1;   // u
  dir[1] = 1;   // v
  dir[2] = 0;   // w
  dir[3] = 0;   // rx
  dir[4] = 0;   // ry
  dir[5] = 1;   // rz
}

// ============================= GetStrLabels ==============================

void cPlFrame :: GetStrLabels(int *label)
{
  label[0] = FORCE_X;
  label[1] = FORCE_Y;
  label[2] = MOMENT_Z;
}

// ================================= Read ==================================

void cPlFrame :: Read(void)
{
  // Read element incidence.
					    
  Shape->Read( );
}

// =============================== IntForce ================================

int cPlFrame :: IntForce(cVector &g)
{
  // Get the nodal coordinates.
  
  sNodeCoord coord[2];
  Shape->NodalCoord(coord);
  
  // Compute the bar configuration.
  
  double L, c, s;
  CalcConfig(coord, &L, &c, &s);
  
  // Assembly the transformation matrix.
  
  cMatrix T(6, 6);
  GetTrnMat(c, s, T);
  
  // Get nodal displacements.
  
  cVector u(6);
  NodalDispl(u);
  
  // Compute the local displacements.
  
  cVector ul(6);
  ul = T*u;
  
  // Compute the local stiffness matrix.
  
  cMatrix Kl(6, 6);
  LocStiffMat(L, Kl);
  
  // Compute the local forces.
  
  cVector gl(6);
  gl = Kl*ul;
  
  // Transform to the global system => {g} = [T]t{gl}.
  
  g = t(T)*gl;
  
  return (1);
}

// =============================== StiffMat ================================

void cPlFrame :: StiffMat(cMatrix &K)
{
  // Get the nodal coordinates.
  
  sNodeCoord coord[2];
  Shape->NodalCoord(coord);
  
  // Compute the bar configuration.
  
  double L, c, s;
  CalcConfig(coord, &L, &c, &s);
  
  // Assembly the transformation matrix.
  
  cMatrix T(6, 6);
  GetTrnMat(c, s, T);
  
  // Compute the local stiffness matrix.
  
  cMatrix Kl(6, 6);
  LocStiffMat(L, Kl);
  
  // Compute the global stiffness matrix => [K] = [T]t[Kl][T].
  
  K.Zero( );
  MatTripBtCB(T, Kl, 1.0, K);
}

// =============================== MassMat =================================

void cPlFrame :: MassMat(cMatrix &M)
{
  // Get the nodal coordinates.
  
  sNodeCoord coord[2];
  Shape->NodalCoord(coord);
  
  // Compute the bar configuration.
  
  double L, c, s;
  CalcConfig(coord, &L, &c, &s);

  // Assembly the rotation matrix.

  cMatrix T(6, 6);
  GetTrnMat(c, s, T);
  
  // Assembly the local mass matrix.

  cMaterial *mat = Section->GetMaterial( );
  double rho = mat->GetDensity( );
  double A = ((cSecBar *)Section)->GetA( );
  double coeff = rho*A*L/420.0;

  cMatrix Ml(6, 6);
  Ml[0][0] = 140.0*coeff;     // 1st line
  Ml[1][0] = 0.0;             // 2nd line
  Ml[1][1] = 156.0*coeff;
  Ml[2][0] = 0.0;             // 3rd line
  Ml[2][1] = 22.0*L*coeff;
  Ml[2][2] = 4.0*L*L*coeff;
  Ml[3][0] = 70.0*coeff;      // 4th line
  Ml[3][1] = 0.0;
  Ml[3][2] = 0.0;
  Ml[3][3] = 140.0*coeff;
  Ml[4][0] = 0.0;             // 5th line
  Ml[4][1] = 54.0*coeff;
  Ml[4][2] = 13.0*L*coeff;
  Ml[4][3] = 0.0;
  Ml[4][4] = 156.0*coeff;
  Ml[5][0] = 0.0;             // 6th line
  Ml[5][1] = -13.0*L*coeff;
  Ml[5][2] = -3.0*L*L*coeff;
  Ml[5][3] = 0.0;
  Ml[5][4] = -22.0*L*coeff;
  Ml[5][5] = 4.0*L*L*coeff;

  // Compute the global mass matrix => [M] = [T]t[Ml][T].

  M.Zero( );
  MatTripBtCB(T, Ml, 1.0, M);
}

// ============================== GeomStiff ================================

void cPlFrame :: GeomStiff(cMatrix &G)
{
  // Get the nodal coordinates.

  sNodeCoord coord[2];
  Shape->NodalCoord(coord);
  
  // Get the nodal displaements.

  cVector u(6);
  NodalDispl(u);

  // Compute the bar configuration.

  double L, c, s;
  CalcConfig(coord, &L, &c, &s);

  // Assembly the transformation matrix.

  cMatrix T(6, 6);
  GetTrnMat(c, s, T);

  // Compute the local displacements.

  cVector ul(6);
  ul = T*u;

  // Compute the normal force.

  cMatMec *mat = Section->GetMaterial( )->GetMecProp( );
  double matpar[2];
  mat->GetParam(matpar);
  double EA = matpar[0]*((cSecBar *)Section)->GetA( );
  double e = (ul[3] - ul[0])/L;   // Compute the axial strain
  double N = EA*e;                // Compute the normal force

  // Assembly the local geometric stiffness.

  cMatrix Gl(6, 6);
  Gl.Zero( );

  Gl[1][1] = 36.0;
  Gl[1][2] = Gl[2][1] =  3.0*L;
  Gl[1][4] = Gl[4][1] = -36.0;
  Gl[1][5] = Gl[5][1] =  3.0*L;

  Gl[2][2] = 4.0*L*L;
  Gl[2][4] = Gl[4][2] = -3.0*L;
  Gl[2][5] = Gl[5][2] = -L*L;

  Gl[4][4] = 36.0;
  Gl[4][5] = Gl[5][4] = -3.0*L;

  Gl[5][5] = 4.0*L*L;

  // Compute the global geometric stiffness [G] = [T]t[Gl][T].

  G.Zero( );
  MatTripBtCB(T, Gl, N/(30.0*L), G);
}

// ============================== NodalStress ==============================

void cPlFrame :: NodalStress(cMatrix &S)
{
  // Get the nodal coordiantes.
  
  sNodeCoord coord[2];
  Shape->NodalCoord(coord);
  
  // Compute the bar configuration.
  
  double L, c, s;
  CalcConfig(coord, &L, &c, &s);
  
  // Assembly the transformation matrix.
  
  cMatrix T(6, 6);
  GetTrnMat(c, s, T);
  
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
  
  // Store in the matrix format.
  
  if (OutputConv == DIRECT_STIFFNESS)
  {
    S[0][0] = gl[0];
    S[0][1] = gl[1];
    S[0][2] = gl[2];
    S[1][0] = gl[3];
    S[1][1] = gl[4];
    S[1][2] = gl[5];
  }
  else
  {
    S[0][0] = -gl[0];
    S[0][1] =  gl[1];
    S[0][2] = -gl[2];
    S[1][0] =  gl[3];
    S[1][1] = -gl[4];
    S[1][2] =  gl[5];
  }
}

// ============================= IntPntStress ==============================

void cPlFrame :: IntPntStress(cMatrix &Strain, cMatrix &Stress)
{
  Strain.Zero( );
  Stress.Zero( );
}

// -------------------------------------------------------------------------
// Protected methods:
//

// =============================== CalcConfig ==============================
//
// This method computes the geometric configuration (angle and length) of a
// given bar.
//
//   coord - vector of nodal coordinates                               (in)
//   L     - element length                                           (out)
//   c     - cosine of bar angle                                      (out)
//   s     - sine of bar angle                                        (out)
//
void cPlFrame :: CalcConfig(sNodeCoord *coord, double *L, double *c, double *s)
{
  double dx = coord[1].x - coord[0].x;
  double dy = coord[1].y - coord[0].y;

  *L = sqrt(dx*dx + dy*dy);
  *c = dx/(*L);
  *s = dy/(*L);
}

// =============================== GetTrnMat ===============================
//
// This method computes the transformation matrix between the local and
// global coordinate systems.
//
//   c - cosine of the element angle with the x axis                   (in)
//   s - sine of the element angle with the x axis                     (in)
//   T - transformation matrix (6x6)                                  (out)
//
void cPlFrame :: GetTrnMat(double c, double s, cMatrix &T)
{
  T.Zero( );
  
  T[0][0] = T[1][1] = c;
  T[0][1] = s;
  T[1][0] = -s;
  T[2][2] = 1.0;
  
  T[3][3] = T[4][4] = c;
  T[3][4] = s;
  T[4][3] = -s;
  T[5][5] = 1.0;
}

// ============================== LocStiffMat ==============================
//
// This method computes the stiffness matrix in the local system.
//
//   L - element length                                                (in)
//   Kl - local stiffness matrix (6x6)                                (out)
//
void cPlFrame :: LocStiffMat(double L, cMatrix &Kl)
{
  cMatMec *mat = Section->GetMaterial( )->GetMecProp( );
  double matpar[2];
  mat->GetParam(matpar);
  double EA = matpar[0]*((cSecBar *)Section)->GetA( );
  double EI = matpar[0]*((cSecBar *)Section)->GetIz( );

  double L2 = L*L;    // auxiliary variable
  double L3 = L2*L;   // auxiliary variable
  
  Kl.Zero( );
  
  Kl[0][0] = Kl[3][3] =  EA/L;
  Kl[0][3] = Kl[3][0] = -EA/L;
  
  Kl[1][1] = 12*EI/L3;
  Kl[1][2] = Kl[2][1] =   6*EI/L2;
  Kl[1][4] = Kl[4][1] = -12*EI/L3;
  Kl[1][5] = Kl[5][1] =   6*EI/L2;
  
  Kl[2][2] = 4*EI/L;
  Kl[2][4] = Kl[4][2] = -6*EI/L2;
  Kl[2][5] = Kl[5][2] =  2*EI/L;
  
  Kl[4][4] = 12*EI/L3;
  Kl[4][5] = Kl[5][4] = -6*EI/L2;
  
  Kl[5][5] = 4*EI/L;
}

// =============================== EqvForces ===============================
//
// This method returns the equivalent nodal forces due to the loads
// applied at the element.
//
//   fl - equivalent nodal forces in the local system                 (out)
//
void cPlFrame :: EqvForces(cVector &fl)
{
  // Get the equivalent nodal forces in the global system.

  cVector fg(6);
  fg.Zero( );
  double t = cControl::GetTotTime( );
  cLoad :: EqvForces(this, t, fg);

  // Get the nodal coordinates.

  sNodeCoord coord[2];
  Shape->NodalCoord(coord);

  // Compute the bar configuration.
  
  double L, c, s;
  CalcConfig(coord, &L, &c, &s);
  
  // Assembly the transformation matrix.
  
  cMatrix T(6, 6);
  GetTrnMat(c, s, T);
  
  // Rotate fix end forces to local system => {fl} = [T]{fg}.
  
  fl = T*fg;
}


// -------------------------------------------------------------------------
// cPlFrameNI class:
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Public methods:
//

// ============================== cPlFrameNI ===============================

cPlFrameNI :: cPlFrameNI(int id, cSection *sec) : cPlFrame(id, sec)
{
  Type = PLFRAMENI;
  IntPnt = cIntPoint :: CreateLinePoints(2, GAUSS, &NumIntPnt);
//  IntPnt = cIntPoint :: CreateLobattoLinePoints(3, &NumIntPnt);

  // Create objects for section analysis at each integration point.

  SecAn = new cSecAnalysis*[NumIntPnt];
  for (int i = 0; i < NumIntPnt; i++)
    SecAn[i] = cSecAnalysis::Create(this, &IntPnt[i]);
}

// ============================= ~cPlFrameNI ===============================

cPlFrameNI :: ~cPlFrameNI(void)
{
  delete []IntPnt;
  for (int i = 0; i < NumIntPnt; i++) delete SecAn[i];
  delete []SecAn; 
}

// =============================== IntForce ================================

int cPlFrameNI :: IntForce(cVector &g)
{
  // Get the nodal coordinates.
  
  sNodeCoord coord[2];
  Shape->NodalCoord(coord);
  
  // Compute the bar configuration.
  
  double L,c,s;
  CalcConfig(coord, &L, &c, &s);
  
  // Assembly the transformation matrix.
  
  cMatrix T(6, 6);
  GetTrnMat(c, s, T);
  
  // Get the nodal displacements.
  
  cVector u(6);
  NodalDispl(u);
  
  // Compute the local displacements.
  
  cVector ul(6);
  ul = T*u;
  
  // Create auxiliary vectors and matrices.
  
  cVector eps(2);       // Strain vector
  cVector sig(2);       // Stress vector
  cVector gl(6);        // Internal force vector in local system
  cMatrix B(2, 6);      // Strain-displacement matrix
  
  // Compute the local internal forces by Gauss integration.
  
  gl.Zero( );
  for (int i = 0; i < NumIntPnt; i++)
  {
    // Get integration point data.
      
    sNatCoord p = IntPnt[i].GetCoord( );
    double wgt  = IntPnt[i].GetWeight( );
    
    // Evaluate the cartesian coordinate of current integration point.
    
    double x = L*(1.0 + p.r)/2.0; 
        
    // Compute the strain-displacement matrix.
      
    BMatrix(x, L, B);
    
    // Evaluate the strains and stresses.
    
    eps = B*ul;
    SecAn[i]->Stress(eps, sig);
    
    // Compute {gl} += coeff*[B]t{sig}.
  
    double coeff = wgt*L/2.0;
    MultTAcc(coeff, B, sig, gl);
  }    
  
  // Transform to the global system => {g} = [T]t{gl}.
  
  g = t(T)*gl;
  
  return(1);
}

// =============================== StiffMat ================================

void cPlFrameNI :: StiffMat(cMatrix &Kt)
{
  // Get the nodal coordinates.
  
  sNodeCoord coord[2];
  Shape->NodalCoord(coord);
  
  // Compute the bar configuration.
  
  double L,c,s;
  CalcConfig(coord, &L, &c, &s);
  
  // Assembly the transformation matrix.
  
  cMatrix T(6, 6);
  GetTrnMat(c, s, T);
  
  // Get the nodal displacements.
  
  cVector u(6);
  NodalDispl(u);
  
  // Compute the local displacements.

  cVector ul(6);
  ul = T*u;
  
  // Create auxiliary vectors and matrices.
  
  cVector eps(2);       // Strain vector
  cMatrix B(2, 6);      // Strain-displacement matrix
  cMatrix C(2, 2);      // Constitutive matrix
  cMatrix Kl(6, 6);     // Stiffness matrix in the local system
  
  // Compute the local stiffness matrix by Gauss integration.
  
  Kl.Zero( );
  for (int i = 0; i < NumIntPnt; i++)
  {
    // Get integration point data.
    
    sNatCoord p = IntPnt[i].GetCoord( );
    double wgt  = IntPnt[i].GetWeight( );
    
    // Evaluate the cartesian coordinate of current integration point.
    
    double x = L*(1.0 + p.r)/2.0; 
    
    // Compute the strain-displacement matrix.
    
    BMatrix(x, L, B);
  
    // Compute the strains and the tangent constitutive matrix.
    
    eps = B*ul;
    SecAn[i]->CMatrix(C);
  
    // [Kl] += [B]t[C][B]*coeff.
    
    double coeff = wgt*L/2.0;
    MatTripBtCB(B, C, coeff, Kl);
  }    
  
  // Transform to the global system => [Kt] = [T]t[Kl][T].
  
  Kt.Zero( );
  MatTripBtCB(T, Kl, 1.0, Kt);
}

// ============================== UpdateState ==============================

void cPlFrameNI :: UpdateState(void)
{
  for (int i = 0; i < NumIntPnt; i++) SecAn[i]->UpdateState( );
}

// ============================= IntPntStress ==============================

void cPlFrameNI :: IntPntStress(cMatrix &Strain, cMatrix &Stress)
{
  // Get the nodal coordinates.
  
  sNodeCoord coord[2];
  Shape->NodalCoord(coord);
  
  // Compute the bar configuration.
  
  double L,c,s;
  CalcConfig(coord, &L, &c, &s);
  
  // Assembly the transformation matrix.
  
  cMatrix T(6, 6);
  GetTrnMat(c, s, T);
  
  // Get the nodal displacements.
  
  cVector u(6);
  NodalDispl(u);
  
  // Compute the local displacements.
  
  cVector ul(6);
  ul = T*u;
  
  // Create auxiliary vectors and matrices.

  cVector eps(2);       // Strain vector
  cVector sig(2);       // Stress vector
  cMatrix B(2, 6);      // Strain-displacement matrix
  
  // Loop over integration points.
  
  Strain.Zero( );
  Stress.Zero( );
  double x,x0 = 0.0;
  for (int i = 0; i < NumIntPnt; i++)
  {
    // Get integration point data.
      
    sNatCoord p = IntPnt[i].GetCoord( );
    
    // Evaluate the cartesian coordinate of current integration point.
    
    x = L*(1.0 + p.r)/2.0; 
    if (i == 0) x0 = x;
        
    // Compute the strain-displacement matrix.
      
    BMatrix(x, L, B);
    
    // Evaluate the strains and stresses.
    
    eps = B*ul;
    SecAn[i]->Stress(eps, sig);
    
    // Store the computed values in the given matrices.

    Strain[i][0] = eps[0];
    Strain[i][1] = 0.0;     // Shear
    Strain[i][2] = eps[1];
    Stress[i][0] = sig[0];
    Stress[i][1] = 0.0;     // Shear
    Stress[i][2] = sig[1];
  }    
  
  // Evaluate an approximation of the shear force: Q = dM/dx.

  double dx = x - x0;
  double Q  = (Stress[NumIntPnt-1][2] - Stress[0][2])/dx;
  for (int i = 0; i < NumIntPnt; i++)
    Stress[i][1] = Q;
}

// -------------------------------------------------------------------------
// Protected methods:
//

// ================================ BMatrix ================================
//
// This method computes the strain-displacement matrix
//
//   x - cartesian coordinate                                          (in)
//   L - element length                                                (in)
//   B - strain-displacement matrix (2x6)                             (out)
//
void cPlFrameNI :: BMatrix(double x, double L, cMatrix &B)
{
  double L2 = L*L;
  double L3 = L2*L;
  
  B[0][0] = -1.0/L;
  B[0][1] =  0.0;
  B[0][2] =  0.0;
  B[0][3] =  1.0/L;
  B[0][4] =  0.0;
  B[0][5] =  0.0;
  
  B[1][0] =  0.0;
  B[1][1] = -6.0/L2 + 12.0*x/L3;
  B[1][2] = -4.0/L + 6.0*x/L2;
  B[1][3] =  0.0;
  B[1][4] =  6.0/L2 - 12.0*x/L3;
  B[1][5] = -2.0/L + 6.0*x/L2;
}

// -------------------------------------------------------------------------
// cPlFrameCR class:
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Public methods:
//

// ============================== cPlFrameCR ===============================

cPlFrameCR :: cPlFrameCR(int id, cSection *sec) : cPlFrame(id, sec)
{
  RigRot = 0.0;
}

// ============================= ~cPlFrameCR ===============================

cPlFrameCR :: ~cPlFrameCR(void)
{
}

// =============================== IntForce ================================

int cPlFrameCR :: IntForce(cVector &g)
{
  // Get the nodal coordinates.
  
  sNodeCoord coord[2];
  Shape->NodalCoord(coord);
  
  // Get the nodal displacements.
  
  cVector u(6);
  NodalDispl(u);
    
  // Compute the initial and current configurations.
  
  UpdateGeom(coord, u);
  double Ln, c, s;
  CalcConfig(coord, &Ln, &c, &s);
    
  // Get the local internal force.
  
  cVector gl(3);
  LocIntForce(gl);
  
  // Assembly the transformation matrix => [T].
  
  cMatrix T(3, 6);
  GetTrnMat(Ln, c, s, T);
  
  // Transform to the global system => {g} = [T]t{gl}.
  
  g = t(T)*gl;
  
  return (1);
}

// =============================== StiffMat ================================

void cPlFrameCR :: StiffMat(cMatrix &Kt)
{
  // Get the nodal coordinates.
  
  sNodeCoord coord[2];
  Shape->NodalCoord(coord);
  
  // Get the nodal displacements.
  
  cVector u(6);
  NodalDispl(u);
    
  // Compute the initial and current configuration.
  
  UpdateGeom(coord, u);
  double Ln, c, s;
  CalcConfig(coord, &Ln, &c, &s);
    
  // Get local stiffness matrix.
  
  cMatrix Kl(3, 3);
  LocStiffMat(Kl);
  
  // Assembly the transformation matrix => [T].
  
  cMatrix T(3, 6);
  GetTrnMat(Ln, c, s, T);

  // Transform to the global system => [Kt] = [T]t[Kl][T].
  
  Kt.Zero( );
  MatTripBtCB(T, Kl, 1.0, Kt);
  
  // Get the local internal force.
  
  cVector gl(3);
  LocIntForce(gl);
  
  // Assembly the auxiliary vectors {r} and {z}.
  
  cVector r(6), z(6);
  GetTrnVec(c, s, r, z);
  
  // Add the global geometric stiffness matrix.
  
  double a1 = gl[0]/Ln;
  double a2 = (gl[1] + gl[2])/(Ln*Ln);
  for (int i = 0; i < 6; i++)
    for (int j = 0; j < 6; j++)
      Kt[i][j] += a1*z[i]*z[j] + a2*(z[i]*r[j] + z[j]*r[i]);
}

// ============================== NodalStress ==============================

void cPlFrameCR :: NodalStress(cMatrix &S)
{
  // Compute the internal forces in the local (3 dofs) system.
  
  cVector gl(3);
  LocIntForce(gl);

  // Compute the initial length.
  
  sNodeCoord coord[2];
  Shape->NodalCoord(coord);
  double L0 = CalcLength(coord);
  
  // Compute the internal force in local system with 6 dofs.
  
  cVector g(6);
  g[0] = -gl[0];
  g[1] = (gl[1] + gl[2])/L0;
  g[2] =  gl[1];
  g[3] =  gl[0];
  g[4] = -g[1];
  g[5] =  gl[2];

  // Add the fixed-end forces (- equivforces).
  
  cVector fl(6);
  EqvForces(fl);
  g -= fl;
  
  // Store in the matrix format.
  
  if (OutputConv == DIRECT_STIFFNESS)
  {
    S[0][0] = g[0];
    S[0][1] = g[1];
    S[0][2] = g[2];
    S[1][0] = g[3];
    S[1][1] = g[4];
    S[1][2] = g[5];
  }
  else
  {
    S[0][0] = -g[0];
    S[0][1] =  g[1];
    S[0][2] = -g[2];
    S[1][0] =  g[3];
    S[1][1] = -g[4];
    S[1][2] =  g[5];
  }
}

// ============================== UpdateState ==============================

void cPlFrameCR :: UpdateState(void)
{
  // Get the nodal coordinates.
  
  sNodeCoord coord[2];
  Shape->NodalCoord(coord);
  
  // Get the nodal displacements.
  
  cVector u(6);
  NodalDispl(u);
  
  // Update the rigid rotation.
  
  RigRot = GetCurrRot(coord, u);

  // Update section data.

  UpdateSection( ); 
}

// -------------------------------------------------------------------------
// Protected methods:
//

// =============================== CalcLength ==============================
//
// This method computes the element length.
//
//   coord - vector of nodal coordinates                               (in)
//
double cPlFrameCR :: CalcLength(sNodeCoord *coord)
{
  double dx = coord[1].x - coord[0].x;
  double dy = coord[1].y - coord[0].y;
  double L  = sqrt(dx*dx + dy*dy);
  
  return(L);
}

// =============================== GetCurrRot ==============================
//
// This method computes the current rigid rotation.
//
//   coord - vector of nodal coordinates                               (in)
//   u     - vector fo nodal displacements                             (in)
//
double cPlFrameCR :: GetCurrRot(sNodeCoord *coord, cVector &u)
{
  // Compute initial vector => {u} = {x0}_21
  
  double ux = coord[1].x - coord[0].x;
  double uy = coord[1].y - coord[0].y;
  
  // Compute current vector => {v} = [R] {u}
  
  double s  = sin(RigRot);
  double c  = cos(RigRot);
  double vx = ux*c - uy*s;
  double vy = ux*s + uy*c;
  
  // Compute the displacement difference vector => {d}_21
  
  double dx = u[3] - u[0];   // u2 - u1
  double dy = u[4] - u[1];   // v2 - v1
  
  // Compute new vector => {w} = {xn}_21 = {x0}_21 + {d}_21
  
  double wx = ux + dx;
  double wy = uy + dy;
  
  // Compute incremental rotation (c = v . w and s = v x w)
  
  c = vx*wx + vy*wy;
  s = vx*wy - vy*wx;
  double da = atan2(s, c);
  
  // Return the total rotation.
  
  return(RigRot + da);
}

// =============================== UpdateGeom ==============================
//
// This method updates the coordinates.
//
//   coord - vector of nodal coordinates                               (in)
//   u     - vector fo nodal displacements                             (in)
//
void cPlFrameCR :: UpdateGeom(sNodeCoord *coord, cVector &u)
{
  coord[0].x += u[0];   // x1 += u1
  coord[0].y += u[1];   // y1 += v1
  coord[1].x += u[3];   // x2 += u2
  coord[1].y += u[4];   // y2 += v2
}

// =============================== GetTrnMat ===============================
//
// This method computes the transformation matrix between the local
// and the global coordinate system.
//
//   L - element length                                                (in)
//   c - cosine of the element angle with the x axis                   (in)
//   s - sine of the element angle with the x axis                     (in)
//   T - transformation matrix (3x6)                                  (out)
//
void cPlFrameCR :: GetTrnMat(double L, double c, double s, cMatrix &T)
{
  T[0][0] = -c;
  T[0][1] = -s;
  T[0][2] =  0.0;
  T[0][3] =  c;
  T[0][4] =  s;
  T[0][5] =  0.0;
  
  T[1][0] = -s/L;
  T[1][1] =  c/L;
  T[1][2] =  1.0;
  T[1][3] =  s/L;
  T[1][4] = -c/L;
  T[1][5] =  0.0;
  
  T[2][0] = -s/L;
  T[2][1] =  c/L;
  T[2][2] =  0.0;
  T[2][3] =  s/L;
  T[2][4] = -c/L;
  T[2][5] =  1.0;
}

// =============================== GetTrnVec ===============================
//
// This method computes the auxiliary vectors {r} and {z}.
//
//   c - cosine of the element angle with the x axis                   (in)
//   s - sine of the element angle with the x axis                     (in)
//   r - auxiliary vector (6)                                         (out)
//   z - auxiliary vector (6)                                         (out)
//
void cPlFrameCR :: GetTrnVec(double c, double s, cVector &r, cVector &z)
{
  r[0] = -c;
  r[1] = -s;
  r[2] =  0.0;
  r[3] =  c;
  r[4] =  s;
  r[5] =  0.0;
  
  z[0] =  s;
  z[1] = -c;
  z[2] =  0.0;
  z[3] = -s;
  z[4] =  c;
  z[5] =  0.0;
}

// ================================ EqvForces ==============================
//
// This method returns the equivalent nodal forces due to the loads
// applied at the element.
//
//   fl - equivalent nodal forces in the local system                 (out)
//
void cPlFrameCR :: EqvForces(cVector &fl)
{
  // Get the equivalent nodal forces in the global system.
  
  cVector fg(6);
  fg.Zero( );
  double t = cControl::GetTotTime( );
  cLoad :: EqvForces(this, t, fg);
  
  // Get the nodal coordinates.
  
  sNodeCoord coord[2];
  Shape->NodalCoord(coord);
  
  // Compute the bar configuration.
  
  double L, c, s;
  CalcConfig(coord, &L, &c, &s);

  // Assembly the transformation matrix.
  
  cMatrix T(6,6);
  T.Zero( );
  T[0][0] = T[1][1] = c;
  T[0][1] = s;
  T[1][0] = -s;
  T[2][2] = 1.0;
  
  T[3][3] = T[4][4] = c;
  T[3][4] = s;
  T[4][3] = -s;
  T[5][5] = 1.0;

  // Rotate the fixed-end forces to the local system => {fl} = [T]{fg}.
  
  fl = T*fg;
}


// -------------------------------------------------------------------------
// cPlFrameCR1 class:
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Public methods:
//

// ============================= cPlFrameCR1 ===============================

cPlFrameCR1 :: cPlFrameCR1(int id, cSection *sec) : cPlFrameCR(id, sec)
{
  Type = PLFRAMECR1;
}

// ============================ ~cPlFrameCR1 ===============================

cPlFrameCR1 :: ~cPlFrameCR1(void)
{
}

// ============================= LocIntForce ===============================

void cPlFrameCR1 :: LocIntForce(cVector &gl)
{
  // Get the nodal coordinates.
  
  sNodeCoord coord[2];
  Shape->NodalCoord(coord);
  
  // Get the nodal displacements.
  
  cVector u(6);
  NodalDispl(u);
  
  // Get the current rigid rotation.
  
  double rot = GetCurrRot(coord, u);
  
  // Compute the initial and current configuration.
  
  double L0 = CalcLength(coord);
  UpdateGeom(coord, u);
  double Ln = CalcLength(coord);
  
  // Compute the local displacements.
  
  cVector ul(3);
  ul[0] = Ln - L0;
  ul[1] = u[2] - rot;
  ul[2] = u[5] - rot;
  
  // Assembly the local stiffness matrix.
  
  cMatMec *mat = Section->GetMaterial( )->GetMecProp( );
  double matpar[2];
  mat->GetParam(matpar);
  double EA = matpar[0]*((cSecBar *)Section)->GetA( );
  double EI = matpar[0]*((cSecBar *)Section)->GetIz( );
  
  cMatrix Kl(3, 3);
  Kl.Zero( );
  Kl[0][0] = EA/L0;
  Kl[1][1] = Kl[2][2] = 4.0*EI/L0;
  Kl[1][2] = Kl[2][1] = 2.0*EI/L0;

  // Compute the local internal force vector => {gl} = [Kl]{ul}.
  
  gl = Kl*ul;
}

// ============================= LocStiffMat ===============================

void cPlFrameCR1 :: LocStiffMat(cMatrix &Kl)
{
  // Get nodal coordinates.
  
  sNodeCoord coord[2];
  Shape->NodalCoord(coord);
  
  // Compute the element length.
  
  double L0 = CalcLength(coord);
  
  // Assembly the local stiffness matrix.
  
  cMatMec *mat = Section->GetMaterial( )->GetMecProp( );
  double matpar[2];
  mat->GetParam(matpar);
  double EA = matpar[0]*((cSecBar *)Section)->GetA( );
  double EI = matpar[0]*((cSecBar *)Section)->GetIz( );
  
  Kl.Zero( );
  Kl[0][0] = EA/L0;
  Kl[1][1] = Kl[2][2] = 4.0*EI/L0;
  Kl[1][2] = Kl[2][1] = 2.0*EI/L0;
}


// -------------------------------------------------------------------------
// cPlFrameCR1NI class:
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Public methods:
//

// ============================ cPlFrameCR1NI ==============================

cPlFrameCR1NI :: cPlFrameCR1NI(int id, cSection *sec) : 
                 cPlFrameCR1(id, sec)
{
  Type = PLFRAMECR1NI;
  IntPnt = cIntPoint :: CreateLinePoints(2, GAUSS, &NumIntPnt);
//  IntPnt = cIntPoint :: CreateLobattoLinePoints(3, &NumIntPnt);

  // Create objects for section analysis at each integration point.

  SecAn = new cSecAnalysis*[NumIntPnt];
  for (int i = 0; i < NumIntPnt; i++)
    SecAn[i] = cSecAnalysis::Create(this, &IntPnt[i]);
}

// =========================== ~cPlFrameCR1NI ==============================

cPlFrameCR1NI :: ~cPlFrameCR1NI(void)
{
  delete []IntPnt;
  for (int i = 0; i < NumIntPnt; i++) delete SecAn[i];
  delete []SecAn; 
}

// ============================= LocIntForce ===============================

void cPlFrameCR1NI :: LocIntForce(cVector &gl)
{
  // Get the nodal coordinates.
  
  sNodeCoord coord[2];
  Shape->NodalCoord(coord);
  
  // Get the nodal displacements.
  
  cVector u(6);
  NodalDispl(u);
  
  // Get the current rigid rotation.
  
  double rot = GetCurrRot(coord, u);
  
  // Compute the initial and current configuration.
  
  double L0 = CalcLength(coord);
  UpdateGeom(coord, u);
  double Ln = CalcLength(coord);
  
  // Compute the local displacements.
  
  cVector ul(3);
  ul[0] = Ln - L0;
  ul[1] = u[2] - rot;
  ul[2] = u[5] - rot;
  
  // Create auxiliary vectors and matrices.
  
  cVector eps(2);       // Strain vector
  cVector sig(2);       // Stress vector
  cMatrix B(2, 3);      // Strain-displacement matrix
  
  // Compute the local internal forces by Gauss integration.
  
  gl.Zero( );
  for (int i = 0; i < NumIntPnt; i++)
  {
    // Get integration point data.
      
    sNatCoord p = IntPnt[i].GetCoord( );
    double wgt  = IntPnt[i].GetWeight( );
    
    // Evaluate the cartesian coordinate of current integration point.
    
    double x = L0*(1.0 + p.r)/2.0; 
        
    // Compute the strain-displacement matrix.
      
    BMatrix(x, L0, B);
    
    // Evaluate the strains and stresses.
    
    eps = B*ul;
    SecAn[i]->Stress(eps, sig);
    
    // Compute {gl} += coeff*[B]t{sig}.
  
    double coeff = wgt*L0/2.0;
    MultTAcc(coeff, B, sig, gl);
  }    
}

// ============================= LocStiffMat ===============================

void cPlFrameCR1NI :: LocStiffMat(cMatrix &Kl)
{
  // Get the nodal coordinates.
  
  sNodeCoord coord[2];
  Shape->NodalCoord(coord);
  
  // Get the nodal displacements.
  
  cVector u(6);
  NodalDispl(u);
  
  // Get the current rigid rotation.
  
  double rot = GetCurrRot(coord, u);
  
  // Compute the initial and current configuration.
  
  double L0 = CalcLength(coord);
  UpdateGeom(coord, u);
  double Ln = CalcLength(coord);
  
  // Compute the local displacements.
  
  cVector ul(3);
  ul[0] = Ln - L0;
  ul[1] = u[2] - rot;
  ul[2] = u[5] - rot;
  
  // Create auxiliary vectors and matrices.
  
  cVector eps(2);       // Strain vector
  cMatrix B(2, 3);      // Strain-displacement matrix
  cMatrix C(2, 2);      // Constitutive matrix
  
  // Compute the local stiffness matrix by Gauss integration.
  
  Kl.Zero( );
  for (int i = 0; i < NumIntPnt; i++)
  {
    // Get integration point data.
    
    sNatCoord p = IntPnt[i].GetCoord( );
    double wgt  = IntPnt[i].GetWeight( );
    
    // Evaluate the cartesian coordinate of current integration point.
    
    double x = L0*(1.0 + p.r)/2.0; 
    
    // Compute the strain-displacement matrix.
    
    BMatrix(x, L0, B);
  
    // Compute the strains and the tangent constitutive matrix.
    
 //   eps = B*ul;
    SecAn[i]->CMatrix(C);
  
    // [Kl] += [B]t[C][B]*coeff.
    
    double coeff = wgt*L0/2.0;
    MatTripBtCB(B, C, coeff, Kl);
  }    
}

// ============================= UpdateSection =============================

void cPlFrameCR1NI :: UpdateSection(void)
{
  for (int i = 0; i < NumIntPnt; i++) SecAn[i]->UpdateState( );
}

// -------------------------------------------------------------------------
// Protected methods:
//

// ================================ BMatrix ================================
//
// This method computes the strain-displacement matrix
//
//   x - cartesian coordinate                                          (in)
//   L - element length                                                (in)
//   B - strain-displacement matrix (2x3)                             (out)
//
void cPlFrameCR1NI :: BMatrix(double x, double L, cMatrix &B)
{
  double L2 = L*L;
  
  B[0][0] =  1.0/L;
  B[0][1] =  0.0;
  B[0][2] =  0.0;
  
  B[1][0] =  0.0;
  B[1][1] = -4.0/L + 6.0*x/L2;
  B[1][2] = -2.0/L + 6.0*x/L2; 
}


// -------------------------------------------------------------------------
// cPlFrameCR2 class:
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Public methods:
//

// ============================= cPlFrameCR2 ===============================

cPlFrameCR2 :: cPlFrameCR2(int id, cSection *sec) : cPlFrameCR(id, sec)
{
  Type = PLFRAMECR2;
}

// ============================ ~cPlFrameCR2 ===============================

cPlFrameCR2 :: ~cPlFrameCR2(void)
{
}

// ============================= LocIntForce ===============================

void cPlFrameCR2 :: LocIntForce(cVector &gl)
{
  // Get the nodal coordinates.
  
  sNodeCoord coord[2];
  Shape->NodalCoord(coord);
  
  // Get the nodal displacements.
  
  cVector u(6);
  NodalDispl(u);
  
  // Get the current rotation.
  
  double rot = GetCurrRot(coord, u);
  
  // Compute the initial and current configuration.
  
  double L0 = CalcLength(coord);
  UpdateGeom(coord, u);
  double Ln = CalcLength(coord);
  
  // Compute the local rotations.
  
  double r1 = u[2] - rot;
  double r2 = u[5] - rot;
  
  // Compute the local displacements.
  
  cVector ul(3);
  ul[0] = Ln - L0;
  ul[1] = r1;
  ul[2] = r2;
  
  // Get the incremental strain-displacement matrix.
  
  cMatrix Bim(1, 3);
  BimMatrix(L0, r1, r2, Bim);
  
  // Compute the membrane strain.
  
  double em = (Ln/L0 - 1) + (2.0*r1*r1 - r1*r2 + 2.0*r2*r2)/30;
  
  // Compute the internal force due to membrane effects => {gm} = [Bim]t*N*L.
  
  cMatMec *mat = Section->GetMaterial( )->GetMecProp( );
  double matpar[2];
  mat->GetParam(matpar);
  double E  = matpar[0];
  double EA = E*((cSecBar *)Section)->GetA( );
  double N  = EA*em;
  double NL = N*L0;
  
  cVector gm(3);
  for (int i = 0; i < 3; i++)
    gm[i] = Bim[0][i]*NL;
  
  // Compute the internal force due to bending effects => {gb} = [Kbe]{ul}.

  double EI = matpar[0]*((cSecBar *)Section)->GetIz( );
  
  cMatrix Kb(3, 3);
  Kb.Zero( );
  Kb[1][1] = Kb[2][2] = 4.0*EI/L0;
  Kb[1][2] = Kb[2][1] = 2.0*EI/L0;
  
  cVector gb(3);
  gb = Kb*ul;
  
  // Compute the local internal force => {gl} = {gm} + {gb}.
  
  gl = gm + gb;
}

// ============================= LocStiffMat ===============================

void cPlFrameCR2 :: LocStiffMat(cMatrix &Kl)
{
  // Get the nodal coordinates.
  
  sNodeCoord coord[2];
  Shape->NodalCoord(coord);
  
  // Get the nodal displacements.
  
  cVector u(6);
  NodalDispl(u);
  
  // Get the current rotation.
  
  double rot = GetCurrRot(coord, u);
  
  // Compute the initial and current configuration.
  
  double L0 = CalcLength(coord);
  UpdateGeom(coord, u);
  double Ln = CalcLength(coord);
  
  // Compute the local rotations.
  
  double r1 = u[2] - rot;
  double r2 = u[5] - rot;
  
  // Get the incremental strain-displacement matrix.
  
  cMatrix Bim(1, 3);
  BimMatrix(L0, r1, r2, Bim);
  
  // Compute the membrane strain.
  
  double em = (Ln/L0 - 1) + (2.0*r1*r1 - r1*r2 + 2.0*r2*r2)/30;
  
  // Compute the local displacements.
  
  cVector ul(3);
  ul[0] = Ln - L0;
  ul[1] = r1;
  ul[2] = r2;

  // Compute the local geometric stiffness matrix => [Kg] = NL*[A].
  
  cMatMec *mat = Section->GetMaterial( )->GetMecProp( );
  double matpar[2];
  mat->GetParam(matpar);
  double E  = matpar[0];
  double EA = E*((cSecBar *)Section)->GetA( );
  double N  = EA*em;
  double NL = N*L0;
  
  cMatrix A(3,3);
  A.Zero( );
  A[1][1] = A[2][2] =  4.0/30.0;
  A[1][2] = A[2][1] = -1.0/30.0;
  
  cMatrix Kg(3,3);
  Kg = NL*A;
  
  // Compute the bending stiffness matrix.
  
  double EI = matpar[0]*((cSecBar *)Section)->GetIz( );
  
  cMatrix Kb(3, 3);
  Kb.Zero( );
  Kb[1][1] = Kb[2][2] = 4.0*EI/L0;
  Kb[1][2] = Kb[2][1] = 2.0*EI/L0;
  
  // Compute the membrane stiffness matrix => [Km] = [Bim]tEAL[Bim].
    
  cMatrix Km(3, 3);
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      Km[i][j] = EA*L0*Bim[0][j]*Bim[0][i];
    
  // Compute the local stiffness matrix => [Kl] = [Kb] + [Km] + [Kg].
  
  Kl  = Kb + Km;
  Kl += Kg;
}

// -------------------------------------------------------------------------
// Protected methods:
//

// =============================== BimMatrix ===============================
//
// This method computes the incremental strain-displacement matrix due to
// local quadratic Green-Lagrange ("shallow arch") term (1/2*v,x^2).
//
//   Bim - incremental strain-displacement matrix                     (out)
//
void cPlFrameCR2 :: BimMatrix(double L, double r1, double r2, cMatrix &Bim)
{
  Bim.Zero( );
  Bim[0][0] = 1.0/L;
  Bim[0][1] = (4.0*r1 - r2)/30.0;
  Bim[0][2] = (4.0*r2 - r1)/30.0;
}


// -------------------------------------------------------------------------
// cPlFrameCR2NI class:
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Public methods:
//

// ============================ cPlFrameCR2NI ==============================

cPlFrameCR2NI :: cPlFrameCR2NI(int id, cSection *sec) : 
                 cPlFrameCR2(id, sec)
{
  Type = PLFRAMECR2NI;
  IntPnt = cIntPoint :: CreateLinePoints(2, GAUSS, &NumIntPnt);
//  IntPnt = cIntPoint :: CreateLobattoLinePoints(3, &NumIntPnt);

  // Create objects for section analysis at each integration point.

  SecAn = new cSecAnalysis*[NumIntPnt];
  for (int i = 0; i < NumIntPnt; i++)
    SecAn[i] = cSecAnalysis::Create(this, &IntPnt[i]);
}

// =========================== ~cPlFrameCR2NI ==============================

cPlFrameCR2NI :: ~cPlFrameCR2NI(void)
{
  delete []IntPnt;
  for (int i = 0; i < NumIntPnt; i++) delete SecAn[i];
  delete []SecAn; 
}

// ============================= LocIntForce ===============================

void cPlFrameCR2NI :: LocIntForce(cVector &gl)
{
  // Get the nodal coordinates.
  
  sNodeCoord coord[2];
  Shape->NodalCoord(coord);
  
  // Get the nodal displacements.
  
  cVector u(6);
  NodalDispl(u);
  
  // Get the current rigid rotation.
  
  double rot = GetCurrRot(coord, u);
  
  // Compute the initial and current configuration.
  
  double L0 = CalcLength(coord);
  UpdateGeom(coord, u);
  double Ln = CalcLength(coord);
  
  // Compute the local rotations.
  
  double r1 = u[2] - rot;
  double r2 = u[5] - rot;
  
  // Create auxiliary vectors and matrices.
  
  cVector eps(2);       // Strain vector
  cVector sig(2);       // Stress vector
  cMatrix Bbar(2, 3);   // Incremental strain-displacement matrix
  
  // Compute the membrane strain.
  
  eps[0] = (Ln/L0 - 1) + (2.0*r1*r1 - r1*r2 + 2.0*r2*r2)/30;
  
  // Compute the local internal forces by Gauss integration.
  
  gl.Zero( );
  double L02 = L0*L0;
  for (int i = 0; i < NumIntPnt; i++)
  {
    // Get integration point data.
      
    sNatCoord p = IntPnt[i].GetCoord( );
    double wgt  = IntPnt[i].GetWeight( );
    
    // Evaluate the cartesian coordinate of current integration point.
    
    double x = L0*(1.0 + p.r)/2.0; 
        
    // Compute the bending strain (curvature).
      
    eps[1] = (-4.0/L0 + 6.0*x/L02)*r1 + (-2.0/L0 + 6.0*x/L02)*r2;
    
    // Compute the incremental strain-displacement matrix.
    
    BbarMatrix(x, L0, r1, r2, Bbar);

    // Evaluate the section stresses.
    
    SecAn[i]->Stress(eps, sig);
    
    // Compute {gl} += coeff*[B]t{sig}.
  
    double coeff = wgt*L0/2.0;
    MultTAcc(coeff, Bbar, sig, gl);
  }    
}

// ============================= LocStiffMat ===============================

void cPlFrameCR2NI :: LocStiffMat(cMatrix &Kl)
{
  // Get the nodal coordinates.
  
  sNodeCoord coord[2];
  Shape->NodalCoord(coord);
  
  // Get the nodal displacements.
  
  cVector u(6);
  NodalDispl(u);
  
  // Get the current rigid rotation.
  
  double rot = GetCurrRot(coord, u);
  
  // Compute the initial and current configuration.
  
  double L0 = CalcLength(coord);
  UpdateGeom(coord, u);
  double Ln = CalcLength(coord);
  
  // Compute the local rotations.
  
  double r1 = u[2] - rot;
  double r2 = u[5] - rot;
  
  // Create auxiliary vectors and matrices.
  
  cVector eps(2);       // Strain vector
  cVector sig(2);       // Stress vector
  cMatrix C(2, 2);      // Constitutive matrix
  cMatrix Bbar(2, 3);   // Incremental strain-displacement matrix
  
  // Compute the membrane strain.
  
  eps[0] = (Ln/L0 - 1) + (2.0*r1*r1 - r1*r2 + 2.0*r2*r2)/30;
  
  // Compute the local stiffness matrix by Gauss integration.
  
  Kl.Zero( );
  double L02 = L0*L0;
  double NL  = 0.0;
  for (int i = 0; i < NumIntPnt; i++)
  {
    // Get integration point data.
    
    sNatCoord p = IntPnt[i].GetCoord( );
    double wgt  = IntPnt[i].GetWeight( );
    
    // Evaluate the cartesian coordinate of current integration point.
    
    double x = L0*(1.0 + p.r)/2.0; 
    
    // Compute the bending strain (curvature).
      
    eps[1] = (-4.0/L0 + 6.0*x/L02)*r1 + (-2.0/L0 + 6.0*x/L02)*r2;
    
    // Compute the section stresses and tangent constitutive matrix.
    
    SecAn[i]->Stress(eps, sig);
    SecAn[i]->CMatrix(C);
  
    // Compute the incremental strain-displacement matrix.
    
    BbarMatrix(x, L0, r1, r2, Bbar);

    // [Kl] += [Bbar]t[C][Bbar]*coeff.
    
    double coeff = wgt*L0/2.0;
    MatTripBtCB(Bbar, C, coeff, Kl);

    // Integration of normal force.

    NL += coeff*sig[0];
  }    

  // Add the local geometric stiffness matrix.

  Kl(1, 1) += NL*4.0/30.0;
  Kl(1, 2) -= NL/30.0;
  Kl(2, 1) -= NL/30.0;
  Kl(2, 2) += NL*4.0/30.0;
}

// ============================= UpdateSection =============================

void cPlFrameCR2NI :: UpdateSection(void)
{
  for (int i = 0; i < NumIntPnt; i++) SecAn[i]->UpdateState( );
}

// -------------------------------------------------------------------------
// Protected methods:
//

// =============================== BbarMatrix ==============================
//
// This method computes the incremental strain-displacement matrix.
//
//   Bbar - incremental strain-displacement matrix (2x3)              (out)
//
void cPlFrameCR2NI :: BbarMatrix(double x, double L, double r1, double r2, 
                                 cMatrix &Bbar)
{
  double L2 = L*L;

  Bbar[0][0] = 1.0/L;
  Bbar[0][1] = (4.0*r1 - r2)/30.0;
  Bbar[0][2] = (4.0*r2 - r1)/30.0;
  
  Bbar[1][0] =  0.0;
  Bbar[1][1] = -4.0/L + 6.0*x/L2;
  Bbar[1][2] = -2.0/L + 6.0*x/L2; 
}


// -------------------------------------------------------------------------
// cPlFrameTL1 class:
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Public methods:
//

// ============================= cPlFrameTL1 ===============================

cPlFrameTL1 :: cPlFrameTL1(int id, cSection *sec) : cPlFrame(id, sec)
{
  Type = PLFRAMETL1;
}

// ============================ ~cPlFrameTL1 ===============================

cPlFrameTL1 :: ~cPlFrameTL1(void)
{
}

// =============================== IntForce ================================

int cPlFrameTL1 :: IntForce(cVector &g)
{
  // Get the nodal coordinates.
  
  sNodeCoord coord[2];
  Shape->NodalCoord(coord);
  
  // Compute the bar configuration.
  
  double L, c, s;
  CalcConfig(coord, &L, &c, &s);
  
  // Assembly the transformation matrix.
  
  cMatrix T(6, 6);
  GetTrnMat(c, s, T);
  
  // Get the nodal displacements.
  
  cVector u(6);
  NodalDispl(u);
  
  // Compute the local displacements.
  
  cVector ul(6);
  ul = T*u;
  
  // Compute the membrane strain-displacement matrices.
  
  cMatrix B0m(1, 6);
  cMatrix Blm1(1, 6);
  cMatrix Blm2(1, 6);
  B0mMatrix(L, B0m);
  Blm1Matrix(L, ul, Blm1);
  Blm2Matrix(L, ul, Blm2);

  cMatrix Blm(1, 6);
  cMatrix Bim(1, 6);
  Blm = Blm1 + Blm2;
  Bim = B0m + Blm;

  // Compute the membrane strain => em = ([B0m] + 1/2*[Blm]){u}.
  
  double em = 0.0;
  for (int i = 0; i < 6; i++)
    em += (B0m[0][i] + 0.5*Blm[0][i])*ul[i];

  // Compute the membrane internal force vector => {gm} = [Bim]t*N*L.
  
  cMatMec *mat = Section->GetMaterial( )->GetMecProp( );
  double matpar[2];
  mat->GetParam(matpar);
  double E  = matpar[0];
  double EA = E*((cSecBar *)Section)->GetA( );
  double N  = EA*em;
  double NL = N*L;
  
  cVector gm(6);
  gm.Zero( );
  for (int i = 0; i < 6; i++)
    gm[i] = Bim[0][i]*NL;

  // Compute the bending internal force vector => {gb} = [Kbe]{ul}.
  
  cMatrix Kbe(6, 6);
  KbeMatrix(L, Kbe);
  cVector gb(6);
  gb = Kbe*ul;
  
  // Compute the local internal force vector => {gl} = {gm} + {gb}.
  
  cVector gl(6);
  gl = gm + gb;
  
  // Transform to the global system => {g} = [T]t{gl}.
  
  g = t(T)*gl;
  
  return (1);  
}

// =============================== StiffMat ================================

void cPlFrameTL1 :: StiffMat(cMatrix &K)
{
  // Get the nodal coordinates.
  
  sNodeCoord coord[2];
  Shape->NodalCoord(coord);
  
  // Compute the bar configuration.
  
  double L, c, s;
  CalcConfig(coord, &L, &c, &s);
  
  // Assembly the transformation matrix
  
  cMatrix T(6,6);
  GetTrnMat(c, s, T);
  
  // Get the nodal displacements.
  
  cVector u(6);
  NodalDispl(u);
  
  // Compute the local displacements.
  
  cVector ul(6);
  ul = T*u;
  
  // Compute the membrane strain-displacement matrices.
  
  cMatrix B0m(1, 6);
  cMatrix Blm1(1, 6);
  cMatrix Blm2(1, 6);
  B0mMatrix(L, B0m);
  Blm1Matrix(L, ul, Blm1);
  Blm2Matrix(L, ul, Blm2);

  cMatrix Blm(1, 6);
  cMatrix Bim(1, 6);
  Blm = Blm1 + Blm2;
  Bim = B0m + Blm;
  
  // Compute the membrane strain => em = ([B0m] + 1/2*[Blm]){u}.
  
  double em = 0.0;
  for (int i = 0; i < 6; i++)
    em += (B0m[0][i] + 0.5*Blm[0][i])*ul[i];

  // Compute [A] matrix.
  
  cMatrix A(6, 6);
  AMatrix(L, A);

  // Compute the membrane geometric stiffness matrix.
  
  cMatMec *mat = Section->GetMaterial( )->GetMecProp( );
  double matpar[2];
  mat->GetParam(matpar);
  double E  = matpar[0];
  double EA = E*((cSecBar *)Section)->GetA( );
  double N  = EA*em;
  double NL = N*L;

  cMatrix Km1(6, 6);         // [Km1] = ([B0m]t[B0m] + [A])*N*L
  Km1.Zero( );
  Km1 += NL*A;               // [Km1] += [A]*N*L
  int i, j;                  // auxiliary variable
  for (i = 0; i < 6; i++)    // [Km1] += [B0m]t[B0m]*N*L
    for (j = 0; j < 6; j++)
      Km1[i][j] += NL*B0m[0][i]*B0m[0][j];
  
  // Compute the membrane elastic stiffness matrix.

  cMatrix Km2(6, 6);         // [Km2] = [Bim]t*EA*[Bim]*L
  Km2.Zero( );
  for (i = 0; i < 6; i++)
    for (j = 0; j < 6; j++)
      Km2[i][j] += EA*L*Bim[0][i]*Bim[0][j];
  
  // Compute the membrane stiffness matrix => [Kme] = [Km1] + [Km2].

  cMatrix Kme(6, 6);
  Kme = Km1 + Km2;

  // Compute the bending stiffness matrix => [Kbe]=[Bl]t*EI*[Bl].
  
  cMatrix Kbe(6, 6);
  KbeMatrix(L, Kbe);

  // Compute the element stiffness matrix => [Kte] = [Kme] + [Kbe].
  
  cMatrix Kte(6, 6);
  Kte = Kme + Kbe;
  
  // Compute the stiffness matrix => [K] = [T]t[Kte][T].
  
  K.Zero( );
  MatTripBtCB(T, Kte, 1.0, K);
}


// -------------------------------------------------------------------------
// Protected methods:
//

// ================================ AMatrix ================================
//
// This method computes the auxiliary matrix.
//
//   L  -  element length                                              (in)
//   A  -  auxiliary matrix                                           (out)
//
void cPlFrameTL1 :: AMatrix(double L, cMatrix &A)
{
  double L2 = L*L;  // auxialiary variable
  
  // Assembly the [A] matrix.
  
  A.Zero( );
  A[1][1] = 6.0/(5.0*L2);
  A[1][2] = A[1][5] = A[2][1] = A[5][1] = 1.0/(10.0*L);
  A[1][4] = A[4][1] = - 6/(5.0*L2);
  
  A[2][2] = 2.0/15.0;
  A[2][4] = A[4][2] = -1.0/(10.0*L);
  A[2][5] = A[5][2] = -2.0/15.0;
  
  A[4][4] = 6.0/(5.0*L2);
  A[4][5] = A[5][4] = -1.0/(10.0*L);
  
  A[5][5] = 2.0/15.0;
}

// =============================== B0mMatrix ===============================
//
// This method computes the B0m matrix (linear) related to membrane strains.
//
//   L    -  element length                                            (in)
//   B0m  -  B0m matrix                                               (out)
//
void cPlFrameTL1 :: B0mMatrix(double L, cMatrix &B0m)
{
  B0m[0][0] = -1.0/L;
  B0m[0][1] = 0.0;
  B0m[0][2] = 0.0;
  B0m[0][3] = 1.0/L;
  B0m[0][4] = 0.0;
  B0m[0][5] = 0.0;
}

// =============================== Blm1Matrix ==============================
//
// This method computes the Blm1 matrix (nonlinear) related to membrane 
// strains.
//
//   L     -  element length                                           (in)
//   u     -  nodal displacement vector                                (in)
//   Blm1  -  Blm1 matrix                                             (out)
//
void cPlFrameTL1 :: Blm1Matrix(double L, cVector &u, cMatrix &Blm1)
{
  double L2 = L*L;

  Blm1[0][0] = u[0]/L2 - u[3]/L2;
  Blm1[0][1] = 0.0;
  Blm1[0][2] = 0.0;
  Blm1[0][3] = -u[0]/L2 + u[3]/L2;
  Blm1[0][4] = 0.0;
  Blm1[0][5] = 0.0;
}

// =============================== Blm2Matrix ==============================
//
// This method computes the Blm2 matrix (nonlinear) related to membrane
// strains.
//
//   L     -  element length                                           (in)
//   u     -  nodal displacement vector                                (in)
//   Blm2  -  Blm2 matrix                                             (out)
//
void cPlFrameTL1 :: Blm2Matrix(double L, cVector &u, cMatrix &Blm2)
{
  double L2 = L*L;
  
  Blm2[0][0] = 0.0;
  Blm2[0][1] = 6*u[1]/(5*L2) + u[2]/(10*L) - 6*u[4]/(5*L2) + u[5]/(10*L);
  Blm2[0][2] = u[1]/(10*L) + 2*u[2]/15 - u[4]/(10*L) - u[5]/30;
  Blm2[0][3] = 0.0;
  Blm2[0][4] = -6*u[1]/(5*L2) - u[2]/(10*L) + 6*u[4]/(5*L2) - u[5]/(10*L);
  Blm2[0][5] = u[1]/(10*L) - u[2]/30 - u[4]/(10*L) + 2*u[5]/15;
}

// =============================== KbeMatrix ===============================
//
// This method computes the bending stiffness matrix.
//
//   L    -  element length                                            (in)
//   Kbe  -  Kbe matrix                                               (out)
//
void cPlFrameTL1 :: KbeMatrix(double L, cMatrix &Kbe)
{
  // Get the material parameter.
  
  cMatMec *mat = Section->GetMaterial( )->GetMecProp( );
  double matpar[2];
  mat->GetParam(matpar);
  double E  = matpar[0];
  double EI = E*((cSecBar *)Section)->GetIz( );
  double L2 = L*L;
  double L3 = L2*L;
  
  // Assembly the bending matrix.
  
  Kbe.Zero( );
  Kbe[1][1] =  12.0*EI/L3;
  Kbe[1][2] = Kbe[2][1] =   6.0*EI/L2;
  Kbe[1][4] = Kbe[4][1] = -12.0*EI/L3;
  Kbe[1][5] = Kbe[5][1] =   6.0*EI/L2;
  
  Kbe[2][2] =  4.0*EI/L;
  Kbe[2][4] = Kbe[4][2] = -6.0*EI/L2;
  Kbe[2][5] = Kbe[5][2] =  2.0*EI/L;
  
  Kbe[4][4] = 12.0*EI/L3;
  Kbe[4][5] = Kbe[5][4] = -6.0*EI/L2;
  
  Kbe[5][5] = 4.0*EI/L;
}

// -------------------------------------------------------------------------
// cPlFrameTL2 class:
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Public methods:
//

// ============================= cPlFrameTL2 ===============================

cPlFrameTL2 :: cPlFrameTL2(int id, cSection *sec) : cPlFrameTL1(id, sec)
{
  Type = PLFRAMETL2;
}

// ============================= ~cPlFrameTL2 ==============================

cPlFrameTL2 :: ~cPlFrameTL2(void)
{
}

// =============================== IntForce ================================

int cPlFrameTL2 :: IntForce(cVector &g)
{
  // Get the nodal coordinates.
  
  sNodeCoord coord[2];
  Shape->NodalCoord(coord);
  
  // Compute the bar configuration.
  
  double L, c, s;
  CalcConfig(coord, &L, &c, &s);
  
  // Assembly the transformation matrix.
  
  cMatrix T(6, 6);
  GetTrnMat(c, s, T);
  
  // Get the nodal coordinates.
  
  cVector u(6);
  NodalDispl(u);
  
  // Compute the local displacements.
  
  cVector ul(6);
  ul = T*u;
  
  // Compute the membrane strain-displacement matrices.
  
  cMatrix Bim(1, 6);
  cMatrix B0m(1, 6);
  cMatrix Blm(1, 6);
  B0mMatrix(L, B0m);
  Blm2Matrix(L, ul, Blm);
  Bim = B0m + Blm;

  // Compute the membrane strain => em = ([B0m] + 1/2*[Blm]){u}.
  
  double em = 0.0;
  for (int i = 0; i < 6; i++)
    em += (B0m[0][i] + 0.5*Blm[0][i])*ul[i];

  // Compute the membrane internal force vector => {gm} = [Bim]t*N*L.
  
  cMatMec *mat = Section->GetMaterial( )->GetMecProp( );
  double matpar[2];
  mat->GetParam(matpar);
  double E  = matpar[0];
  double EA = E*((cSecBar *)Section)->GetA( );
  double N  = EA*em;
  double NL = N*L;

  cVector gm(6);
  gm.Zero( );
  for (int i = 0; i < 6; i++)
    gm[i] = Bim[0][i]*NL;

  // Compute the internal force due to bending effects => {gb} = [Kbe]{ul}
  
  cMatrix Kbe(6, 6);
  KbeMatrix(L, Kbe);
  cVector gb(6);
  gb = Kbe*ul;

  // Compute the local internal force => {gl} = {gm} + {gb}.
  
  cVector gl(6);
  gl = gm + gb;
  
  // Transform to the global system => {g} = [T]t{gl}.
  
  g = t(T)*gl;
  
  return (1);
}

// =============================== StiffMat ================================

void cPlFrameTL2 :: StiffMat(cMatrix &K)
{
  // Get the nodal coordinates.
  
  sNodeCoord coord[2];
  Shape->NodalCoord(coord);
  
  // Compute the bar configuration.
  
  double L, c, s;
  CalcConfig(coord, &L, &c, &s);
  
  // Assembly the transformation matrix
  
  cMatrix T(6, 6);
  GetTrnMat(c, s, T);
  
  // Get the nodal displacements.
  
  cVector u(6);
  NodalDispl(u);
  
  // Compute the local displacements.
  
  cVector ul(6);
  ul = T*u;
  
  // Compute the membrane strain-displacement matrices.
  
  cMatrix Bim(1, 6);
  cMatrix B0m(1, 6);
  cMatrix Blm(1, 6);
  B0mMatrix(L, B0m);
  Blm2Matrix(L, ul, Blm);
  Bim  = B0m + Blm;
  
  // Compute the membrane strain => em = ([B0m] + 1/2*[Blm]){u}.
  
  double em = 0.0;  // em = (B0m[0][0] + 0.5*Blm[0][0])*ul[0] + ...
  for (int i = 0; i < 6; i++)
    em += (B0m[0][i] + 0.5*Blm[0][i])*ul[i];

  // Get [A] matrix.
  
  cMatrix A(6, 6);
  AMatrix(L, A);
  
  // Compute the membrane geometric stiffness matrix.
  
  cMatMec *mat = Section->GetMaterial( )->GetMecProp( );
  double matpar[2];
  mat->GetParam(matpar);
  double E  = matpar[0];
  double EA = E*((cSecBar *)Section)->GetA( );
  double N  = EA*em;
  double NL = N*L;
  
  cMatrix Km1(6, 6);     // [Km1] = [A]*N*L
  Km1 = NL*A;
  
  // Compute the membrane elastic stiffness matrix.

  cMatrix Km2(6, 6);     // [Km2] = [Bim]t*EA*[Bim]*L
  Km2.Zero( );
  int i,j;
  for (i = 0; i < 6; i++)
    for (j = 0; j < 6; j++)
      Km2[i][j] += EA*L*Bim[0][i]*Bim[0][j];
  
  cMatrix Kme(6, 6);
  Kme = Km1 + Km2;
  
  // Compute the bending elastic stiffness matrix.
  
  cMatrix Kbe(6, 6);
  KbeMatrix(L, Kbe);
  
  // Compute the element stiffness matrix => [Kte] = [Kme] + [Kbe].

  cMatrix Kte(6, 6);
  Kte = Kme + Kbe;
  
  // Compute the stiffness matrix => [K] = [T]t[Kte][T].
  
  K.Zero( );
  MatTripBtCB(T, Kte, 1.0, K);
}

// ======================================================= End of file =====
