// -------------------------------------------------------------------------
// anmodel.cpp - implementation of the Analysis Model class.
// -------------------------------------------------------------------------
// Created:      29-Apr-2005     Evandro Parente Junior
//
// Modified:     15-Jun-2011     Iuri Barcelos Rocha
//               Implementation of ThickPlate and ShallowShell models.
//
// Modified:     21-Oct-2013     Evandro Parente Junior
//               Implementation of CMatrixOrtho.
//
// Modified:     28-Oct-2013     Iuri Barcelos Rocha
//               Implementation of ExpandStress function.
//
// Modified:     08-Apr-2015     Elias Saraiva Barroso
//               Input/output operations using c++ streams.
//
// Modified:     05-Feb-2019     Samir Parente Auad
//               Implementation of DonnellShell model.
//
// Modified:     16-Sep-2022     Pedro Ygor Rodrigues Mesquita
//               Implementation of HSDTPlate model.
//
// Modified:     01-Feb-2023     Renan Melo Barros 
//               Implementation of CalcDev, CalcTen, CalcFOID and CalcSOID.
//
// Modified:     10-Mar-2023     Evandro Parente Junior
//               Implementation of DimBMatrix, DimQMatrix, and QMatrix.
//
// Modified:     16-Feb-2024     Evandro Parente Junior
//               Implementation of DimNMatrix and NMatrix for HSDTPlate.
// -------------------------------------------------------------------------

#include <math.h>
#include <iostream>
#include <iomanip>
#include <stdlib.h>

using namespace std;

#include "anmodel.h"
#include "vec.h"
#include "mat.h"
#include "gbldef.h"

// -------------------------------------------------------------------------
// Public methods:
//

// ============================== CreateModel ==============================

cAnModel *cAnModel :: CreateModel(eAnmType type)
{
  cAnModel *anm = 0;

  switch(type)
  {
    case PLANE_STRESS:
     anm = new cPlaneStress( );
    break;

    case PLANE_STRAIN:
     anm = new cPlaneStrain( );
    break;

    case AXISYMMETRIC:
     anm = new cAxisymmetric( );
    break;

    case SOLID:
     anm = new cSolid( );
    break;

    case THICK_PLATE:
      anm = new cThickPlate( );
    break;

    case HSDT_PLATE:
    anm = new cHSDTPlate( );
    break;

    case SHALLOW_SHELL:
      anm = new cShallowShell( );
    break;

    case DONNELL_SHELL:
      anm = new cDonnellShell( );
    break;

    case SHELL:
      anm = new cShell( );
    break;

    case INTERFACE_2D:
      anm = new cInterface2D( );
    break;

    case INTERFACE_3D:
      anm = new cInterface3D( );
    break;

    case PLANE_HEAT_TRANSFER:
      anm = new cPlaneHeatTransfer( );
    break;
  }

  return(anm);
}

// ================================ cAnModel ===============================

cAnModel :: cAnModel(void)
{
}

// =============================== ~cAnModel ===============================

cAnModel :: ~cAnModel(void)
{
}

// ================================ GetMecMod ==============================

cMecModel* cAnModel :: GetMecMod(void)
{
  cMecModel *mod = dynamic_cast<cMecModel*>(this);

  if (!mod)
  {
    cout << "Invalid cast to mechanical analysis model!" << endl;
    exit(0);
  }

  return mod;
}

// ================================ GetHeatMod =============================

cHeatModel* cAnModel :: GetHeatMod(void)
{
  cHeatModel *mod = dynamic_cast<cHeatModel*>(this);

  if (!mod)
  {
    cout << "Invalid cast to thermal analysis model!" << endl;
    exit(0);
  }

  return mod;
}


// -------------------------------------------------------------------------
// Class cMecModel:
// -------------------------------------------------------------------------

// ============================= cMecModel =================================

cMecModel :: cMecModel(void)
{
}

// ============================= ~cMecModel =================================

cMecModel :: ~cMecModel(void)
{
}

// ============================= GetDimBnlMatrix ===========================

int cMecModel :: GetDimBnlMatrix(void)
{
  cout << "\ncMecModel :: GetDimBnlMatrix not implemented.  ";
  cout << "Geometric nonlinear analysis not available for this element!\n\n";
  return(0);
}

// ================================= NMatrix ===============================

void cMecModel :: NMatrix(int nn, double *shpfunc, cMatrix &N) 
{
  cout << "\ncMecModel :: NMatrix not implemented.  ";
  cout << "Dynamic analysis not available for this element!\n\n";
  N.Zero( );
}

// ================================= NMatrix ===============================

void cMecModel :: NMatrix(int nn, double *shpfunc, sNodeCoord *drvshp, cMatrix &N) 
{
  cout << "\ncMecModel :: NMatrix not implemented.  ";
  cout << "Dynamic analysis not available for this element!\n\n";
  N.Zero( );
}

// ================================= AlphaVec ==============================

void cMecModel :: AlphaVec(double *param, cVector &Alpha)
{
  cout << "\ncMecModel :: AlphaVec not implemented.  ";
  cout << "Thermal loading not available for this element!\n\n";
  exit(0);
}

// =============================== AlphaVecOrtho ===========================

void cMecModel :: AlphaVecOrtho(double *param, cVector &Alpha)
{
  cout << "\ncMecModel :: AlphaVecOrtho not implemented.  ";
  cout << "Thermal loading not available for this element!\n\n";
  exit(0);
}


// -------------------------------------------------------------------------
// Class cHeatMod:
// -------------------------------------------------------------------------

// ============================= cHeatModel ================================

cHeatModel :: cHeatModel(void)
{
}

// ============================= ~cHeatModel ===============================

cHeatModel :: ~cHeatModel(void)
{
}

// -------------------------------------------------------------------------
// Class cPlaneStress:
// -------------------------------------------------------------------------

// ============================= cPlaneStress ==============================

cPlaneStress :: cPlaneStress(void)
{
  Type = PLANE_STRESS;
}

// ============================ ~cPlaneStress ==============================

cPlaneStress :: ~cPlaneStress(void)
{
}

// ============================== GetActDir ================================

void cPlaneStress :: GetActDir(int *dir)
{
  dir[0] = 1;   // u
  dir[1] = 1;   // v
}

// ============================== GetStrLabels =============================

void cPlaneStress :: GetStrLabels(int *label)
{
  label[0] = STRESS_XX;
  label[1] = STRESS_YY;
  label[2] = STRESS_XY;
}

// ============================== PrincStress ==============================

void cPlaneStress :: PrincStress(cVector &str, cVector &prs)
{
  double sxx = str[0];
  double syy = str[1];
  double sxy = str[2];
  double cen = (sxx + syy)/2.0;        // Mohr's circle center
  double ds  = (sxx - syy)/2.0;
  double rad = sqrt(ds*ds + sxy*sxy);  // Mohr's circle radius
  prs[0] = cen + rad;
  prs[1] = cen - rad;
  prs[2] = 0.0;
  prs.Sort(DESCENDING);
}

// ================================ CalcI1 =================================

double cPlaneStress :: CalcI1(cVector &str)
{
  double I1 = str[0] + str[1];
  return(I1);
}

// ================================ CalcJ2 =================================

double cPlaneStress :: CalcJ2(cVector &str)
{
  double sxx = str[0];
  double syy = str[1];
  double sxy = str[2];
  double dxy = sxx - syy;
  double J2 = (dxy*dxy + sxx*sxx + syy*syy)/6.0 + sxy*sxy;
  return(J2);
}

// ================================ GradI1 =================================

void cPlaneStress :: GradI1(cVector &str, cVector &g)
{
  g.Zero( );
  g[0] = 1.0;
  g[1] = 1.0;
}

// ================================ GradJ2 =================================

void cPlaneStress :: GradJ2(cVector &str, cVector &g)
{
  g.Zero( );
  g[0] = (2.0*str[0] - str[1])/3.0;
  g[1] = (2.0*str[1] - str[0])/3.0;
  g[2] =  2.0*str[2];
}

// ================================= QMatrix ===============================

void cPlaneStress :: QMatrix(double *param, cMatrix &Q)
{
  double E  = param[0];
  double nu = param[1];

  Q.Zero( );
  Q[0][0] = Q[1][1] = E/(1.0 - nu*nu);
  Q[0][1] = Q[1][0] = nu*E/(1.0 - nu*nu);
  Q[2][2] = 0.5*E/(1.0 + nu);
}

// ============================== QMatrixOrtho =============================

void cPlaneStress :: QMatrixOrtho(double *param, cMatrix &Q)
{
  double E1   = param[0];
  double E2   = param[1];
//  double E3   = param[2];
  double Nu12 = param[3];
//  double Nu13 = param[4];
//  double Nu23 = param[5];
  double G12  = param[6];

  double Nu21 = Nu12*E2/E1;
  double d = 1.0/(1.0 - Nu12*Nu21);

  Q.Zero( );
  Q[0][0] = d*E1;
  Q[1][1] = d*E2;
  Q[0][1] = Q[1][0] = d*Nu12*E2;
  Q[2][2] = G12;
}

// ================================= DMatrix ===============================

void cPlaneStress :: DMatrix(double *param, cMatrix &D)
{
  double E  = param[0];
  double nu = param[1];

  D.Zero( );
  D[0][0] = D[1][1] = 1.0/E;
  D[0][1] = D[1][0] = -nu/E;
  D[2][2] = 2.0*(1.0 + nu)/E;
}

// ================================ TMatrix ================================

void cPlaneStress :: TMatrix(double angle, cMatrix &T)
{
  T.Zero( );

  // Trigonometric constants.

  double c  = cos(angle);
  double s  = sin(angle);
  double c2 = c*c;
  double s2 = s*s;

  // Transformation matrix.

  T[0][0] = T[1][1] = c2;
  T[0][1] = T[1][0] = s2;
  T[0][2] =  c*s;
  T[1][2] = -c*s;
  T[2][0] = -2.0*c*s;
  T[2][1] =  2.0*c*s;
  T[2][2] = c2 - s2;
}

// ================================= BMatrix ===============================

void cPlaneStress :: BMatrix(int nn, double *shpfunc, double *mapfunc,
                             sNodeCoord *coord, sNodeCoord *drvshp, cMatrix &B)
{
  for (int i = 0; i < nn; i++)
  {
    B[0][2*i  ] = drvshp[i].x;
    B[0][2*i+1] = 0.0;
    B[1][2*i  ] = 0.0;
    B[1][2*i+1] = drvshp[i].y;
    B[2][2*i  ] = drvshp[i].y;
    B[2][2*i+1] = drvshp[i].x;
  }
}

// ================================ BlMatrix ===============================

void cPlaneStress :: BlMatrix(int nn, double *shpfunc, double *mapfunc,
                              cVector &u, sNodeCoord *coord,
                              sNodeCoord *drvshp, cMatrix &Bl)
{
  // Compute the displacement derivatives.

  double l11 = 0.0;
  double l12 = 0.0;
  double l21 = 0.0;
  double l22 = 0.0;
  for (int i = 0; i < nn; i++)
  {
    l11 += u[2*i  ]*drvshp[i].x;
    l12 += u[2*i  ]*drvshp[i].y;
    l21 += u[2*i+1]*drvshp[i].x;
    l22 += u[2*i+1]*drvshp[i].y;
  }

  // Assembly [Bl].

  for (int i = 0; i < nn; i++)
  {
    Bl[0][2*i]   = l11*drvshp[i].x;
    Bl[0][2*i+1] = l21*drvshp[i].x;
    Bl[1][2*i]   = l12*drvshp[i].y;
    Bl[1][2*i+1] = l22*drvshp[i].y;
    Bl[2][2*i]   = l11*drvshp[i].y + l12*drvshp[i].x;
    Bl[2][2*i+1] = l21*drvshp[i].y + l22*drvshp[i].x;
  }
}

// ================================ BnlMatrix ===============================

void cPlaneStress :: BnlMatrix(int nn, double *shpfunc, double *mapfunc,
                               sNodeCoord *coord, sNodeCoord *drvshp,
                               cMatrix &Bnl)
{
  for (int i = 0; i < nn; i++)
  {
    Bnl[0][2*i]   = drvshp[i].x;
    Bnl[0][2*i+1] = 0.0;
    Bnl[1][2*i]   = drvshp[i].y;
    Bnl[1][2*i+1] = 0.0;
    Bnl[2][2*i]   = 0.0;
    Bnl[2][2*i+1] = drvshp[i].x;
    Bnl[3][2*i]   = 0.0;
    Bnl[3][2*i+1] = drvshp[i].y;
  }
}

// ================================ SMatrix =================================

void cPlaneStress :: SMatrix(cVector &str, cMatrix &S)
{
  S.Zero( );
  S[0][0] = str[0];
  S[1][0] = str[2];
  S[0][1] = str[2];
  S[1][1] = str[1];

  S[2][2] = str[0];
  S[3][2] = str[2];
  S[2][3] = str[2];
  S[3][3] = str[1];
}

// ================================= NMatrix ===============================

void cPlaneStress :: NMatrix(int nn, double *shpfunc, cMatrix &N)
{
  for (int i = 0; i < nn; i++)
  {
    N[0][2*i  ] = shpfunc[i];
    N[0][2*i+1] = 0.0;
    N[1][2*i  ] = 0.0;
    N[1][2*i+1] = shpfunc[i];
  }
}

// ============================== ExpandStress =============================

void cPlaneStress :: ExpandStress(cVector &sig, cVector &expsig)
{
  expsig[0] = sig[0];  // Sxx
  expsig[1] = sig[1];  // Syy
  expsig[2] = 0.0;     // Szz
  expsig[3] = sig[2];  // Sxy
  expsig[4] = 0.0;     // Sxz
  expsig[5] = 0.0;     // Syz
}

// ================================= AlphaVec ==============================

void cPlaneStress :: AlphaVec(double *param, cVector &Alpha)
{
  double alpha = param[2];   // Coefficient of Thermal Expansion

  Alpha.Zero( );
  Alpha[0] = alpha;
  Alpha[1] = alpha;
}

// =============================== AlphaVecOrtho ===========================

void cPlaneStress :: AlphaVecOrtho(double *param, cVector &Alpha)
{
  double alpha1 = param[9];  // Coefficient of Thermal Expansion - dir 1
  double alpha2 = param[10]; // Coefficient of Thermal Expansion - dir 2 

  Alpha.Zero( );
  Alpha[0] = alpha1;
  Alpha[1] = alpha2;
}


// -------------------------------------------------------------------------
// Class cPlaneStrain:
// -------------------------------------------------------------------------

// ============================= cPlaneStrain ==============================

cPlaneStrain :: cPlaneStrain(void)
{
  Type = PLANE_STRAIN;
}

// ============================ ~cPlaneStrain ==============================

cPlaneStrain :: ~cPlaneStrain (void)
{
}

// ============================== GetActDir ================================

void cPlaneStrain :: GetActDir(int *dir)
{
  dir[0] = 1;   // u
  dir[1] = 1;   // v
}

// ============================== GetStrLabels =============================

void cPlaneStrain :: GetStrLabels(int *label)
{
  label[0] = STRESS_XX;
  label[1] = STRESS_YY;
  label[2] = STRESS_ZZ;
  label[3] = STRESS_XY;
}

// ============================== PrincStress ==============================

void cPlaneStrain :: PrincStress(cVector &str, cVector &prs)
{
  double sxx = str[0];
  double syy = str[1];
  double szz = str[2];
  double sxy = str[3];
  double cen = (sxx + syy)/2.0;        // Mohr's circle center
  double ds  = (sxx - syy)/2.0;
  double rad = sqrt(ds*ds + sxy*sxy);  // Mohr's circle radius
  prs[0] = cen + rad;
  prs[1] = cen - rad;
  prs[2] = szz;
  prs.Sort(DESCENDING);
}

// ================================ CalcI1 =================================

double cPlaneStrain :: CalcI1(cVector &str)
{
  double I1 = str[0] + str[1] + str[2];
  return(I1);
}

// ================================ CalcJ2 =================================

double cPlaneStrain :: CalcJ2(cVector &str)
{
  double dxy = str[0] - str[1];
  double dxz = str[0] - str[2];
  double dyz = str[1] - str[2];
  double J2 = (dxy*dxy + dxz*dxz + dyz*dyz)/6.0 + str[3]*str[3];
  return(J2);
}

// ================================ GradI1 =================================

void cPlaneStrain :: GradI1(cVector &str, cVector &g)
{
  g.Zero( );
  g[0] = 1.0;
  g[1] = 1.0;
  g[2] = 1.0;
}

// ================================ GradJ2 =================================

void cPlaneStrain :: GradJ2(cVector &str, cVector &g)
{
  g.Zero( );
  g[0] = (2.0*str[0] - str[1] - str[2])/3.0;
  g[1] = (2.0*str[1] - str[0] - str[2])/3.0;
  g[2] = (2.0*str[2] - str[0] - str[1])/3.0;
  g[3] =  2.0*str[3];
}

// ================================= QMatrix ===============================

void cPlaneStrain :: QMatrix(double *param, cMatrix &Q)
{
  double E  = param[0];
  double nu = param[1];
  double d  = (1.0 + nu)*(1.0 - 2.0*nu);

  Q.Zero( );
  Q[0][0] = Q[1][1] = Q[2][2] = E*(1.0 - nu)/d;
  Q[0][1] = Q[1][0] =
  Q[0][2] = Q[2][0] =
  Q[1][2] = Q[2][1] = E*nu/d;
  Q[3][3] = 0.5*E/(1.0 + nu);
}

// ============================== QMatrixOrtho =============================

void cPlaneStrain :: QMatrixOrtho(double *param, cMatrix &Q)
{
  double E1   = param[0];
  double E2   = param[1];
  double E3   = param[2];
  double Nu12 = param[3];
  double Nu13 = param[4];
  double Nu23 = param[5];
  double G12  = param[6];

  double Nu21 = Nu12*E2/E1;
  double Nu31 = Nu13*E3/E1;
  double Nu32 = Nu23*E3/E2;

  double d = 1.0 - Nu12*Nu21 - Nu13*Nu31 - Nu23*Nu32 -
             Nu12*Nu23*Nu31 - Nu13*Nu21*Nu32;

  Q.Zero( );
  Q[0][0] = E1*(1.0 - Nu23*Nu32)/d;
  Q[1][1] = E2*(1.0 - Nu13*Nu31)/d;
  Q[2][2] = E3*(1.0 - Nu12*Nu21)/d;

  Q[0][1] = Q[1][0] = E2*(Nu12 + Nu13*Nu32)/d;
  Q[0][2] = Q[2][0] = E3*(Nu13 + Nu12*Nu23)/d;
  Q[1][2] = Q[2][1] = E3*(Nu23 + Nu21*Nu13)/d;

  Q[3][3] = G12;
}


// ================================= DMatrix ===============================

void cPlaneStrain :: DMatrix(double *param, cMatrix &D)
{
  double E  = param[0];
  double nu = param[1];

  D.Zero( );
  D[0][0] = D[1][1] = D[2][2] = 1.0/E;
  D[0][1] = D[1][0] =
  D[0][2] = D[2][0] =
  D[1][2] = D[2][1] = -nu/E;
  D[3][3] = 2.0*(1.0 + nu)/E;
}

// ================================= BMatrix ===============================

void cPlaneStrain :: BMatrix(int nn, double *shpfunc, double *mapfunc,
                             sNodeCoord *coord, sNodeCoord *drvshp, cMatrix &B)
{
  for (int i = 0; i < nn; i++)
  {
    B[0][2*i  ] = drvshp[i].x;
    B[0][2*i+1] = 0.0;
    B[1][2*i  ] = 0.0;
    B[1][2*i+1] = drvshp[i].y;
    B[2][2*i  ] = 0.0;
    B[2][2*i+1] = 0.0;
    B[3][2*i  ] = drvshp[i].y;
    B[3][2*i+1] = drvshp[i].x;

  }
}

// ================================ BlMatrix ===============================

void cPlaneStrain :: BlMatrix(int nn, double *shpfunc, double *mapfunc,
                              cVector &u, sNodeCoord *coord,
                              sNodeCoord *drvshp, cMatrix &Bl)
{
  // Compute the displacement derivatives.

  double l11 = 0.0;
  double l12 = 0.0;
  double l21 = 0.0;
  double l22 = 0.0;
  for (int i = 0; i < nn; i++)
  {
    l11 += u[2*i  ]*drvshp[i].x;
    l12 += u[2*i  ]*drvshp[i].y;
    l21 += u[2*i+1]*drvshp[i].x;
    l22 += u[2*i+1]*drvshp[i].y;
  }

  // Assembly [Bl].

  for (int i = 0; i < nn; i++)
  {
    Bl[0][2*i  ] = l11*drvshp[i].x;
    Bl[0][2*i+1] = l21*drvshp[i].x;
    Bl[1][2*i  ] = l12*drvshp[i].y;
    Bl[1][2*i+1] = l22*drvshp[i].y;
    Bl[2][2*i  ] = l11*drvshp[i].y + l12*drvshp[i].x;
    Bl[2][2*i+1] = l21*drvshp[i].y + l22*drvshp[i].x;
    Bl[3][2*i  ] = 0.0;
    Bl[3][2*i+1] = 0.0;
  }
}

// ================================ BnlMatrix ===============================

void cPlaneStrain :: BnlMatrix(int nn, double *shpfunc, double *mapfunc,
                               sNodeCoord *coord, sNodeCoord *drvshp,
                               cMatrix &Bnl)
{
  for (int i = 0; i < nn; i++)
  {
    Bnl[0][2*i  ] = Bnl[2][2*i+1] = drvshp[i].x;
    Bnl[1][2*i  ] = Bnl[3][2*i+1] = drvshp[i].y;
    Bnl[0][2*i+1] = Bnl[1][2*i+1] = 0.0;
    Bnl[2][2*i  ] = Bnl[3][2*i  ] = 0.0;
    Bnl[4][2*i  ] = Bnl[4][2*i+1] = 0.0;
  }
}

// ================================= SMatrix ===============================

void cPlaneStrain :: SMatrix(cVector &str, cMatrix &S)
{
  S.Zero( );
  S[0][0] = str[0];
  S[1][0] = str[3];
  S[0][1] = str[3];
  S[1][1] = str[1];
  S[2][2] = str[0];
  S[3][2] = str[3];
  S[2][3] = str[3];
  S[3][3] = str[1];
  S[4][4] = str[2];
}

// ================================= NMatrix ===============================

void cPlaneStrain :: NMatrix(int nn, double *shpfunc, cMatrix &N)
{
  for (int i = 0; i < nn; i++)
  {
    N[0][2*i  ] = shpfunc[i];
    N[0][2*i+1] = 0.0;
    N[1][2*i  ] = 0.0;
    N[1][2*i+1] = shpfunc[i];
  }
}

// ============================== ExpandStress =============================

void cPlaneStrain :: ExpandStress(cVector &sig, cVector &expsig)
{
  expsig[0] = sig[0];  // Sxx
  expsig[1] = sig[1];  // Syy
  expsig[2] = sig[2];  // Szz
  expsig[3] = sig[3];  // Sxy
  expsig[4] = 0.0;     // Sxz
  expsig[5] = 0.0;     // Syz
}

// ============================== CalcDev =============================

void cPlaneStrain :: CalcDev(double p, cVector &sig, cVector &sdev)
{
  sdev[0] = sig[0] - p;  // Sxx
  sdev[1] = sig[1] - p;  // Syy
  sdev[2] = sig[2] - p;  // Szz
  sdev[3] = sig[3];  // Sxy
}

// ============================== CalcTen =============================

void cPlaneStrain :: CalcTen(double p, cVector &sdev, cVector &ten)
{
  ten[0] = sdev[0] + p;  // Sxx
  ten[1] = sdev[1] + p;  // Syy
  ten[2] = sdev[2] + p;  // Szz
  ten[3] = sdev[3];  // Sxy
}

// ================================ CalcFOID =================================

void cPlaneStrain :: CalcFOID(cMatrix &FO)
{
  FO.Zero();
  FO[0][0] = FO[1][1] = FO[2][2] = 1.0;
  FO[3][3] = 0.5;
}

// ================================ CalcSOID =================================

void cPlaneStrain :: CalcSOID(cVector &SO)
{
  SO.Zero();
  SO[0] = SO[1] = SO[2] = 1.0;
}

// -------------------------------------------------------------------------
// Class cAxisymmetric:
// -------------------------------------------------------------------------

// ============================= cAxisymmetric =============================

cAxisymmetric :: cAxisymmetric(void)
{
  Type = AXISYMMETRIC;
}

// ============================ ~cAxisymmetric =============================

cAxisymmetric :: ~cAxisymmetric(void)
{
}

// ============================== GetActDir ================================

void cAxisymmetric :: GetActDir(int *dir)
{
  dir[0] = 1;   // u
  dir[1] = 1;   // v
}

// ============================== GetStrLabels =============================

void cAxisymmetric :: GetStrLabels(int *label)
{
  label[0] = STRESS_XX;
  label[1] = STRESS_YY;
  label[2] = STRESS_ZZ;
  label[3] = STRESS_XY;
}

// ============================== PrincStress ==============================

void cAxisymmetric :: PrincStress(cVector &str, cVector &prs)
{
  double sxx = str[0];
  double syy = str[1];
  double szz = str[2];
  double sxy = str[3];
  double cen = (sxx + syy)/2.0;        // Mohr's circle center
  double ds  = (sxx - syy)/2.0;
  double rad = sqrt(ds*ds + sxy*sxy);  // Mohr's circle radius
  prs[0] = cen + rad;
  prs[1] = cen - rad;
  prs[2] = szz;
  prs.Sort(DESCENDING);
}

// ================================ CalcI1 =================================

double cAxisymmetric :: CalcI1(cVector &str)
{
  double I1 = str[0] + str[1] + str[2];
  return(I1);
}

// ================================ CalcJ2 =================================

double cAxisymmetric :: CalcJ2(cVector &str)
{
  double dxy = str[0] - str[1];
  double dxz = str[0] - str[2];
  double dyz = str[1] - str[2];
  double J2 = (dxy*dxy + dxz*dxz + dyz*dyz)/6.0 + str[3]*str[3];
  return(J2);
}

// ================================ GradI1 =================================

void cAxisymmetric :: GradI1(cVector &str, cVector &g)
{
  g.Zero( );
  g[0] = 1.0;
  g[1] = 1.0;
  g[2] = 1.0;
}

// ================================ GradJ2 =================================

void cAxisymmetric :: GradJ2(cVector &str, cVector &g)
{
  g.Zero( );
  g[0] = (2.0*str[0] - str[1] - str[2])/3.0;
  g[1] = (2.0*str[1] - str[0] - str[2])/3.0;
  g[2] = (2.0*str[2] - str[0] - str[1])/3.0;
  g[3] =  2.0*str[3];
}

// ================================= QMatrix ===============================

void cAxisymmetric :: QMatrix(double *param, cMatrix &Q)
{
  double E  = param[0];
  double nu = param[1];
  double d  = (1.0 + nu)*(1.0 - 2.0*nu);

  Q.Zero( );
  Q[0][0] = Q[1][1] = Q[2][2] = E*(1.0 - nu)/d;
  Q[0][1] = Q[1][0] =
  Q[0][2] = Q[2][0] =
  Q[1][2] = Q[2][1] = E*nu/d;
  Q[3][3] = 0.5*E/(1.0 + nu);
}

// ============================== QMatrixOrtho =============================

void cAxisymmetric :: QMatrixOrtho(double *param, cMatrix &Q)
{
  double E1   = param[0];
  double E2   = param[1];
  double E3   = param[2];
  double Nu12 = param[3];
  double Nu13 = param[4];
  double Nu23 = param[5];
  double G12  = param[6];

  double Nu21 = Nu12*E2/E1;
  double Nu31 = Nu13*E3/E1;
  double Nu32 = Nu23*E3/E2;

  double d = 1.0 - Nu12*Nu21 - Nu13*Nu31 - Nu23*Nu32 -
             Nu12*Nu23*Nu31 - Nu13*Nu21*Nu32;

  Q.Zero( );
  Q[0][0] = E1*(1.0 - Nu23*Nu32)/d;
  Q[1][1] = E2*(1.0 - Nu13*Nu31)/d;
  Q[2][2] = E3*(1.0 - Nu12*Nu21)/d;

  Q[0][1] = Q[1][0] = E2*(Nu12 + Nu13*Nu32)/d;
  Q[0][2] = Q[2][0] = E3*(Nu13 + Nu12*Nu23)/d;
  Q[1][2] = Q[2][1] = E3*(Nu23 + Nu21*Nu13)/d;

  Q[3][3] = G12;
}


// ================================= DMatrix ===============================

void cAxisymmetric :: DMatrix(double *param, cMatrix &D)
{
  double E  = param[0];
  double nu = param[1];

  D.Zero( );
  D[0][0] = D[1][1] = D[2][2] = 1.0/E;
  D[0][1] = D[1][0] = 
  D[0][2] = D[2][0] =
  D[1][2] = D[2][1] = -nu/E;
  D[3][3] = 2.0*(1.0 + nu)/E;
}

// ============================= VolCoeff ==================================

double cAxisymmetric :: VolCoeff(double t, int nn, double *mapfunc,
                                sNodeCoord *coord)
{
  // Compute the radial coordinate.

  double r = 0.0;
  for (int i = 0; i < nn; i++) r += mapfunc[i]*coord[i].x;

  return(2.0*PI*r);
}

// ================================= BMatrix ===============================

void cAxisymmetric :: BMatrix(int nn, double *shpfunc, double *mapfunc,
                              sNodeCoord *coord, sNodeCoord *drvshp,
                              cMatrix &B)
{
  // Compute the radial coordinate.

  int i;
  double r = 0.0;
  for (i = 0; i < nn; i++) r += mapfunc[i]*coord[i].x;

  // Assembly [B].

  for (i = 0; i < nn; i++)
  {
    B[0][2*i  ] = drvshp[i].x;
    B[0][2*i+1] = 0.0;
    B[1][2*i  ] = 0.0;
    B[1][2*i+1] = drvshp[i].y;
    B[2][2*i  ] = shpfunc[i]/r;
    B[2][2*i+1] = 0.0;
    B[3][2*i  ] = drvshp[i].y;
    B[3][2*i+1] = drvshp[i].x;
  }
}

// ================================ BlMatrix ===============================

void cAxisymmetric :: BlMatrix(int nn, double *shpfunc, double *mapfunc,
                               cVector &u, sNodeCoord *coord,
                               sNodeCoord *drvshp, cMatrix &Bl)
{
  // Compute the radial coordinate and the displacement derivatives.

  int i;
  double r = 0.0;
  double l11 = 0.0;
  double l12 = 0.0;
  double l21 = 0.0;
  double l22 = 0.0;
  double l33 = 0.0;
  for (i = 0; i < nn; i++)
  {
    r += mapfunc[i]*coord[i].x;
    l11 += u[2*i  ]*drvshp[i].x;
    l12 += u[2*i  ]*drvshp[i].y;
    l21 += u[2*i+1]*drvshp[i].x;
    l22 += u[2*i+1]*drvshp[i].y;
    l33 += u[2*i  ]*shpfunc[i];
  }
  l33 /= r;

  // Assembly [Bl].

  for (i = 0; i < nn; i++)
  {
    Bl[0][2*i  ] = l11*drvshp[i].x;
    Bl[0][2*i+1] = l21*drvshp[i].x;

    Bl[1][2*i  ] = l12*drvshp[i].y;
    Bl[1][2*i+1] = l22*drvshp[i].y;

    Bl[2][2*i  ] = l11*drvshp[i].y + l12*drvshp[i].x;
    Bl[2][2*i+1] = l21*drvshp[i].y + l22*drvshp[i].x;

    Bl[3][2*i  ] = l33*shpfunc[i]/r;
    Bl[3][2*i+1] = 0.0;
  }
}

// ================================ BnlMatrix ===============================

void cAxisymmetric :: BnlMatrix(int nn, double *shpfunc, double *mapfunc,
                                sNodeCoord *coord, sNodeCoord *drvshp,
                                cMatrix &Bnl)
{
  // Compute the radial coordinate.

  int i;
  double r = 0.0;
  for (i = 0; i < nn; i++) r += mapfunc[i]*coord[i].x;

  // Assembly [Bnl].

  for (i = 0; i < nn; i++)
  {
    Bnl[0][2*i  ] = drvshp[i].x;
    Bnl[0][2*i+1] = 0.0;

    Bnl[1][2*i  ] = drvshp[i].y;
    Bnl[1][2*i+1] = 0.0;

    Bnl[2][2*i  ] = 0.0;
    Bnl[2][2*i+1] = drvshp[i].x;

    Bnl[3][2*i  ] = 0.0;
    Bnl[3][2*i+1] = drvshp[i].y;

    Bnl[4][2*i  ] = shpfunc[i]/r;
    Bnl[4][2*i+1] = 0.0;
  }
}

// ================================= SMatrix ===============================

void cAxisymmetric :: SMatrix(cVector &str, cMatrix &S)
{
  S.Zero( );
  S[0][0] = str[0];
  S[1][0] = str[3];
  S[0][1] = str[3];
  S[1][1] = str[1];
  S[2][2] = str[0];
  S[3][2] = str[3];
  S[2][3] = str[3];
  S[3][3] = str[1];
  S[4][4] = str[2];
}

// ================================= NMatrix ===============================

void cAxisymmetric :: NMatrix(int nn, double *shpfunc, cMatrix &N)
{
  for (int i = 0; i < nn; i++)
  {
    N[0][2*i  ] = shpfunc[i];
    N[0][2*i+1] = 0.0;
    N[1][2*i  ] = 0.0;
    N[1][2*i+1] = shpfunc[i];
  }
}

// ============================== ExpandStress =============================

void cAxisymmetric :: ExpandStress(cVector &sig, cVector &expsig)
{
  expsig[0] = sig[0];  // Sxx
  expsig[1] = sig[1];  // Syy
  expsig[2] = sig[2];  // Szz
  expsig[3] = sig[3];  // Sxy
  expsig[4] = 0.0;     // Sxz
  expsig[5] = 0.0;     // Syz
}

// ============================== CalcDev =============================

void cAxisymmetric :: CalcDev(double p, cVector &sig, cVector &sdev)
{
  sdev[0] = sig[0] - p;  // Sxx
  sdev[1] = sig[1] - p;  // Syy
  sdev[2] = sig[2] - p;  // Szz
  sdev[3] = sig[3];  // Sxy
}

// ============================== CalcTen =============================

void cAxisymmetric :: CalcTen(double p, cVector &sdev, cVector &ten)
{
  ten[0] = sdev[0] + p;  // Sxx
  ten[1] = sdev[1] + p;  // Syy
  ten[2] = sdev[2] + p;  // Szz
  ten[3] = sdev[3];  // Sxy
}

// ================================ CalcFOID =================================

void cAxisymmetric :: CalcFOID(cMatrix &FO)
{
  FO.Zero();
  FO[0][0] = FO[1][1] = FO[2][2] = 1.0;
  FO[3][3] = 0.5;
}

// ================================ CalcSOID =================================

void cAxisymmetric :: CalcSOID(cVector &SO)
{
  SO.Zero();
  SO[0] = SO[1] = SO[2] = 1.0;
}

// -------------------------------------------------------------------------
// Class cSolid:
// -------------------------------------------------------------------------

// ================================ cSolid =================================

cSolid :: cSolid(void)
{
  Type = SOLID;
}

// =============================== ~cSolid =================================

cSolid :: ~cSolid(void)
{
}

// ============================== GetActDir ================================

void cSolid :: GetActDir(int *dir)
{
  dir[0] = 1;   // u
  dir[1] = 1;   // v
  dir[2] = 1;   // w
}

// ============================== GetStrLabels =============================

void cSolid :: GetStrLabels(int *label)
{
  label[0] = STRESS_XX;
  label[1] = STRESS_YY;
  label[2] = STRESS_ZZ;
  label[3] = STRESS_XY;
  label[4] = STRESS_XZ;
  label[5] = STRESS_YZ;
}

// ============================== PrincStress ==============================

void cSolid :: PrincStress(cVector &str, cVector &prs)
{
  // Compute the principal stresses solving the cubic equation
  // s^3 - I1*s^2 + I2*s - I3 = 0. Since stress tensor is symmetric, the
  // three eigenvalues (princ stresses) are real.
  //
  double I1 = str[0] + str[1] + str[2];
  double I2 = str[0]*str[1] + str[0]*str[2] + str[1]*str[2] -
              str[3]*str[3] - str[4]*str[4] - str[5]*str[5];
  double I3 = str[0]*str[1]*str[2] - str[0]*str[5]*str[5] -
              str[1]*str[4]*str[4] - str[2]*str[3]*str[3] +
              2.0*str[3]*str[4]*str[5];

  double xn = I1/3.0;
  double yn = -2.0*I1*I1*I1/27.0 + I1*I2/3.0 - I3;
  double dl = sqrt(I1*I1 - 3.0*I2)/3.0;
  double h  = 2.0*dl*dl*dl;
  double tet = acos(-yn/h)/3.0;
  prs[0] = xn + 2.0*dl*cos(tet);
  prs[1] = xn + 2.0*dl*cos(tet + 2.0*PI/3.0);
  prs[2] = xn + 2.0*dl*cos(tet + 4.0*PI/3.0);
  prs.Sort(DESCENDING);
}

// ================================ CalcI1 =================================

double cSolid :: CalcI1(cVector &str)
{
  double I1 = str[0] + str[1] + str[2];
  return(I1);
}

// ================================ CalcJ2 =================================

double cSolid :: CalcJ2(cVector &str)
{
  double dxy = str[0] - str[1];
  double dxz = str[0] - str[2];
  double dyz = str[1] - str[2];
  double J2 = (dxy*dxy + dxz*dxz + dyz*dyz)/6.0 +
              str[3]*str[3] + str[4]*str[4] + str[5]*str[5];
  return(J2);
}

// ================================ GradI1 =================================

void cSolid :: GradI1(cVector &str, cVector &g)
{
  g.Zero( );
  g[0] = 1.0;
  g[1] = 1.0;
  g[2] = 1.0;
}

// ================================ GradJ2 =================================

void cSolid :: GradJ2(cVector &str, cVector &g)
{
  g[0] = (2.0*str[0] - str[1] - str[2])/3.0;
  g[1] = (2.0*str[1] - str[0] - str[2])/3.0;
  g[2] = (2.0*str[2] - str[0] - str[1])/3.0;
  g[3] =  2.0*str[3];
  g[4] =  2.0*str[4];
  g[5] =  2.0*str[5];
}

// ================================= QMatrix ===============================

void cSolid :: QMatrix(double *param, cMatrix &Q)
{
  double E  = param[0];
  double nu = param[1];
  double d  = (1.0 + nu)*(1.0 - 2.0*nu);

  Q.Zero( );
  Q[0][0] = Q[1][1] = Q[2][2] = E*(1.0 - nu)/d;
  Q[0][1] = Q[1][0] =
  Q[0][2] = Q[2][0] =
  Q[1][2] = Q[2][1] = E*nu/d;
  Q[3][3] = Q[4][4] = Q[5][5] = 0.5*E/(1.0 + nu);
}

// ============================== QMatrixOrtho =============================

void cSolid :: QMatrixOrtho(double *param, cMatrix &Q)
{
  double E1   = param[0];
  double E2   = param[1];
  double E3   = param[2];
  double Nu12 = param[3];
  double Nu13 = param[4];
  double Nu23 = param[5];
  double G12  = param[6];
  double G13  = param[7];
  double G23  = param[8];

  double Nu21 = Nu12*E2/E1;
  double Nu31 = Nu13*E3/E1;
  double Nu32 = Nu23*E3/E2;

  double d = 1.0 - Nu12*Nu21 - Nu13*Nu31 - Nu23*Nu32 -
             Nu12*Nu23*Nu31 - Nu13*Nu21*Nu32;

  Q.Zero( );

  Q[0][0] = E1*(1.0 - Nu23*Nu32)/d;
  Q[1][1] = E2*(1.0 - Nu13*Nu31)/d;
  Q[2][2] = E3*(1.0 - Nu12*Nu21)/d;

  Q[0][1] = Q[1][0] = E2*(Nu12 + Nu13*Nu32)/d;
  Q[0][2] = Q[2][0] = E3*(Nu13 + Nu12*Nu23)/d;
  Q[1][2] = Q[2][1] = E3*(Nu23 + Nu21*Nu13)/d;

  Q[3][3] = G12;
  Q[4][4] = G13;
  Q[5][5] = G23;
}

// ================================= DMatrix ===============================

void cSolid :: DMatrix(double *param, cMatrix &D)
{
  double E  = param[0];
  double nu = param[1];

  D.Zero( );
  D[0][0] = D[1][1] = D[2][2] = 1.0/E;
  D[0][1] = D[1][0] =
  D[0][2] = D[2][0] =
  D[1][2] = D[2][1] = -nu/E;
  D[3][3] = D[4][4] = D[5][5] = 2.0*(1.0 + nu)/E;
}

// ================================= TMatrix ===============================

void cSolid :: TMatrix(cVector &e1, cVector &e2, cVector &e3, cMatrix &T)
{
  T.Zero( );

  // Trigonometric constants.

  double l1 = e1[0];
  double m1 = e1[1];
  double n1 = e1[2];

  double l2 = e2[0];
  double m2 = e2[1];
  double n2 = e2[2];

  double l3 = e3[0];
  double m3 = e3[1];
  double n3 = e3[2];

  // Transformation Matrix.

  T[0][0] = l1*l1;
  T[0][1] = m1*m1;
  T[0][2] = n1*n1;
  T[0][3] = l1*m1;
  T[0][4] = n1*l1;
  T[0][5] = m1*n1;

  T[1][0] = l2*l2;
  T[1][1] = m2*m2;
  T[1][2] = n2*n2;
  T[1][3] = l2*m2;
  T[1][4] = n2*l2;
  T[1][5] = m2*n2;

  T[2][0] = l3*l3;
  T[2][1] = m3*m3;
  T[2][2] = n3*n3;
  T[2][3] = l3*m3;
  T[2][4] = n3*l3;
  T[2][5] = m3*n3;

  T[3][0] = 2*l1*l2;
  T[3][1] = 2*m1*m2;
  T[3][2] = 2*n1*n2;
  T[3][3] = l1*m2+l2*m1;
  T[3][4] = n1*l2+n2*l1;
  T[3][5] = m1*n2+m2*n1;

  T[4][0] = 2*l3*l1;
  T[4][1] = 2*m3*m1;
  T[4][2] = 2*n3*n1;
  T[4][3] = l3*m1+l1*m3;
  T[4][4] = n3*l1+n1*l3;
  T[4][5] = m3*n1+m1*n3;

  T[5][0] = 2*l2*l3;
  T[5][1] = 2*m2*m3;
  T[5][2] = 2*n2*n3;
  T[5][3] = l2*m3+l3*m2;
  T[5][4] = n2*l3+n3*l2;
  T[5][5] = m2*n3+m3*n2;
}

// ================================= BMatrix ===============================

void cSolid :: BMatrix(int nn, double *shpfunc, double *mapfunc,
                       sNodeCoord *coord, sNodeCoord *drvshp, cMatrix &B)
{
  for (int i = 0; i < nn; i++)
  {
    B[0][3*i  ] = drvshp[i].x;
    B[0][3*i+1] = 0.0;
    B[0][3*i+2] = 0.0;

    B[1][3*i  ] = 0.0;
    B[1][3*i+1] = drvshp[i].y;
    B[1][3*i+2] = 0.0;

    B[2][3*i  ] = 0.0;
    B[2][3*i+1] = 0.0;
    B[2][3*i+2] = drvshp[i].z;

    B[3][3*i  ] = drvshp[i].y;
    B[3][3*i+1] = drvshp[i].x;
    B[3][3*i+2] = 0.0;

    B[4][3*i  ] = drvshp[i].z;
    B[4][3*i+1] = 0.0;
    B[4][3*i+2] = drvshp[i].x;

    B[5][3*i  ] = 0.0;
    B[5][3*i+1] = drvshp[i].z;
    B[5][3*i+2] = drvshp[i].y;
  }
}

// ================================ BlMatrix ===============================

void cSolid :: BlMatrix(int nn, double *shpfunc, double *mapfunc, cVector &u,
                        sNodeCoord *coord, sNodeCoord *drvshp, cMatrix &Bl)
{
  // Compute the displacement derivatives.

  int i;
  double l11 = 0.0;
  double l12 = 0.0;
  double l13 = 0.0;
  double l21 = 0.0;
  double l22 = 0.0;
  double l23 = 0.0;
  double l31 = 0.0;
  double l32 = 0.0;
  double l33 = 0.0;
  for (i = 0; i < nn; i++)
  {
    l11 += u[3*i  ]*drvshp[i].x;
    l12 += u[3*i  ]*drvshp[i].y;
    l13 += u[3*i  ]*drvshp[i].z;

    l21 += u[3*i+1]*drvshp[i].x;
    l22 += u[3*i+1]*drvshp[i].y;
    l23 += u[3*i+1]*drvshp[i].z;

    l31 += u[3*i+2]*drvshp[i].x;
    l32 += u[3*i+2]*drvshp[i].y;
    l33 += u[3*i+2]*drvshp[i].z;
  }

  // Assembly [Bl].

  Bl.Zero( );
  for (i = 0; i < nn; i++)
  {
    Bl[0][3*i  ] = l11*drvshp[i].x;
    Bl[0][3*i+1] = l21*drvshp[i].x;
    Bl[0][3*i+2] = l31*drvshp[i].x;

    Bl[1][3*i  ] = l12*drvshp[i].y;
    Bl[1][3*i+1] = l22*drvshp[i].y;
    Bl[1][3*i+2] = l32*drvshp[i].y;

    Bl[2][3*i  ] = l13*drvshp[i].z;
    Bl[2][3*i+1] = l23*drvshp[i].z;
    Bl[2][3*i+2] = l33*drvshp[i].z;

    Bl[3][3*i  ] = l11*drvshp[i].y + l12*drvshp[i].x;
    Bl[3][3*i+1] = l21*drvshp[i].y + l22*drvshp[i].x;
    Bl[3][3*i+2] = l31*drvshp[i].y + l32*drvshp[i].x;

    Bl[4][3*i  ] = l11*drvshp[i].z + l13*drvshp[i].x;
    Bl[4][3*i+1] = l21*drvshp[i].z + l23*drvshp[i].x;
    Bl[4][3*i+2] = l31*drvshp[i].z + l33*drvshp[i].x;

    Bl[5][3*i  ] = l12*drvshp[i].z + l13*drvshp[i].y;
    Bl[5][3*i+1] = l22*drvshp[i].z + l23*drvshp[i].y;
    Bl[5][3*i+2] = l32*drvshp[i].z + l33*drvshp[i].y;
  }
}

// ================================ BnlMatrix ===============================

void cSolid :: BnlMatrix(int nn, double *shpfunc, double *mapfunc,
                         sNodeCoord *coord, sNodeCoord *drvshp, cMatrix &Bnl)
{
  Bnl.Zero( );
  for (int i = 0; i < nn; i++)
  {
    Bnl[0][3*i  ] = drvshp[i].x;
    Bnl[1][3*i  ] = drvshp[i].y;
    Bnl[2][3*i  ] = drvshp[i].z;

    Bnl[3][3*i+1] = drvshp[i].x;
    Bnl[4][3*i+1] = drvshp[i].y;
    Bnl[5][3*i+1] = drvshp[i].z;

    Bnl[6][3*i+2] = drvshp[i].x;
    Bnl[7][3*i+2] = drvshp[i].y;
    Bnl[8][3*i+2] = drvshp[i].z;
  }
}

// ================================= SMatrix ===============================

void cSolid :: SMatrix(cVector &str, cMatrix &S)
{
  S.Zero( );
  S[0][0] = str[0];
  S[1][0] = str[3];
  S[2][0] = str[4];
  S[0][1] = str[3];
  S[1][1] = str[1];
  S[2][1] = str[5];
  S[0][2] = str[4];
  S[1][2] = str[5];
  S[2][2] = str[2];

  S[3][3] = str[0];
  S[4][3] = str[3];
  S[5][3] = str[4];
  S[3][4] = str[3];
  S[4][4] = str[1];
  S[5][4] = str[5];
  S[3][5] = str[4];
  S[4][5] = str[5];
  S[5][5] = str[2];

  S[6][6] = str[0];
  S[7][6] = str[3];
  S[8][6] = str[4];
  S[6][7] = str[3];
  S[7][7] = str[1];
  S[8][7] = str[5];
  S[6][8] = str[4];
  S[7][8] = str[5];
  S[8][8] = str[2];
}

// ================================= NMatrix ===============================

void cSolid :: NMatrix(int nn, double *shpfunc, cMatrix &N)
{
  for (int i = 0; i < nn; i++)
  {
    N[0][3*i  ] = shpfunc[i];
    N[0][3*i+1] = 0.0;
    N[0][3*i+2] = 0.0;

    N[1][3*i  ] = 0.0;
    N[1][3*i+1] = shpfunc[i];
    N[1][3*i+2] = 0.0;

    N[2][3*i  ] = 0.0;
    N[2][3*i+1] = 0.0;
    N[2][3*i+2] = shpfunc[i];
  }
}

// ============================== ExpandStress =============================

void cSolid :: ExpandStress(cVector &sig, cVector &expsig)
{
  expsig = sig;
}

// ============================== CalcDev ============================

void cSolid :: CalcDev(double p, cVector &sig, cVector &sdev)
{
  sdev[0] = sig[0] - p;  // Sxx
  sdev[1] = sig[1] - p;  // Syy
  sdev[2] = sig[2] - p;  // Szz
  sdev[3] = sig[3];  // Sxy
  sdev[4] = sig[4];  // Sxz
  sdev[5] = sig[5];  // Syz
}

// ============================== CalcTen =============================

void cSolid :: CalcTen(double p, cVector &sdev, cVector &ten)
{
  ten[0] = sdev[0] + p;  // Sxx
  ten[1] = sdev[1] + p;  // Syy
  ten[2] = sdev[2] + p;  // Szz
  ten[3] = sdev[3];  // Sxy
  ten[4] = sdev[4];  // Sxz
  ten[5] = sdev[5];  // Syz
}

// ================================ CalcFOID =================================

void cSolid :: CalcFOID(cMatrix &FO)
{
  FO.Zero();
  FO[0][0] = FO[1][1] = FO[2][2] = 1.0;
  FO[3][3] = FO[4][4] = FO[5][5] = 0.5;
}

// ================================ CalcSOID =================================

void cSolid :: CalcSOID(cVector &SO)
{
  SO.Zero();
  SO[0] = SO[1] = SO[2] = 1.0;
}
// -------------------------------------------------------------------------
// Class cThickPlate:
// -------------------------------------------------------------------------

// ============================== cThickPlate ==============================

cThickPlate :: cThickPlate(void)
{
  Type = THICK_PLATE;
}

// ============================= ~cThickPlate ==============================

cThickPlate :: ~cThickPlate(void)
{
}

// ============================== GetActDir ================================

void cThickPlate :: GetActDir(int *dir)
{
  dir[2] = 1;   // w
  dir[3] = 1;   // rx
  dir[4] = 1;   // ry
}

// ============================== GetStrLabels =============================

void cThickPlate :: GetStrLabels(int *label)
{
  label[0] = MOMENT_Y;
  label[1] = MOMENT_X;
  label[2] = MOMENT_XY;
  label[3] = FORCE_XZ;
  label[4] = FORCE_YZ;
}

// ============================== PrincStress ==============================

void cThickPlate :: PrincStress(cVector &str, cVector &prs)
{
  cout << "cThickPlate :: PrincStress not implemented!\n";
}

// ================================ QMatrix ================================

void cThickPlate :: QMatrix(double *param, cMatrix &Q)
{
  double E  = param[0];
  double nu = param[1];
  double G  = 0.5*E/(1.0 + nu);
  double d  = 1.0/(1.0 - nu*nu);

  Q.Zero( );

  // Bending terms.

  Q[0][0] = Q[1][1] = d*E;
  Q[0][1] = Q[1][0] = d*nu*E;
  Q[2][2] = G;

  // Transverse shear terms.

  Q[3][3] = Q[4][4] = G;
}

// ============================== QMatrixOrtho =============================

void cThickPlate :: QMatrixOrtho(double *param, cMatrix &Q)
{
  double E1   = param[0];
  double E2   = param[1];
//  double E3   = param[2];
  double Nu12 = param[3];
//  double Nu13 = param[4];
//  double Nu23 = param[5];
  double G12  = param[6];
  double G13  = param[7];
  double G23  = param[8];

  double Nu21 = Nu12*E2/E1;
  double d = 1.0/(1.0 - Nu12*Nu21);

  Q.Zero( );

  // Bending terms.

  Q[0][0] = d*E1;
  Q[1][1] = d*E2;
  Q[0][1] = Q[1][0] = d*Nu12*E2;
  Q[2][2] = G12;

  // Transverse shear terms.

  Q[3][3] = G13;
  Q[4][4] = G23;
}

// ================================= DMatrix ===============================

void cThickPlate :: DMatrix(double *param, cMatrix &D)
{
  cout << "cThickPlate :: DMatrix not implemented!\n";
}

// ================================ TMatrix ================================

void cThickPlate :: TMatrix(double angle, cMatrix &T)
{
  T.Zero( );

  // Trigonometric constants.

  double c  = cos(angle);
  double s  = sin(angle);
  double c2 = c*c;
  double s2 = s*s;

  // Bending transformation matrix.

  T[0][0] = T[1][1] = c2;
  T[0][1] = T[1][0] = s2;
  T[0][2] =  c*s;
  T[1][2] = -c*s;
  T[2][0] = -2.0*c*s;
  T[2][1] =  2.0*c*s;
  T[2][2] = c2 - s2;

  // Shear transformation matrix.

  T[3][3] = T[4][4] = c;
  T[3][4] = -s;
  T[4][3] =  s;
}

// ================================= BMatrix ===============================

void cThickPlate :: BMatrix(int nn, double *shpfunc, double *mapfunc,
                            sNodeCoord *coord, sNodeCoord *drvshp, cMatrix &B)
{
  for (int i = 0; i < nn; i++)
  {
    B[0][3*i  ] =  0.0;
    B[0][3*i+1] =  0.0;
    B[0][3*i+2] = -drvshp[i].x;

    B[1][3*i  ] =  0.0;
    B[1][3*i+1] =  drvshp[i].y;
    B[1][3*i+2] =  0.0;

    B[2][3*i  ] =  0.0;
    B[2][3*i+1] = -drvshp[i].x;
    B[2][3*i+2] =  drvshp[i].y;

    B[3][3*i  ] =  drvshp[i].x;
    B[3][3*i+1] =  0.0;
    B[3][3*i+2] =  shpfunc[i];

    B[4][3*i  ] =  drvshp[i].y;
    B[4][3*i+1] = -shpfunc[i];
    B[4][3*i+2] =  0.0;
  }
}

// ================================ BlMatrix ===============================

void cThickPlate :: BlMatrix(int nn, double *shpfunc, double *mapfunc,
                             cVector &u, sNodeCoord *coord,
                             sNodeCoord *drvshp, cMatrix &Bl)
{
  cout << "cThickPlate :: BlMatrix not implemented!\n";
  exit(0);
}

// =============================== BnlMatrix ===============================

void cThickPlate :: BnlMatrix(int nn, double *shpfunc, double *mapfunc,
                              sNodeCoord *coord, sNodeCoord *drvshp,
                              cMatrix &Bnl)
{
  cout << "cThickPlate :: BnlMatrix not implemented!\n";
  exit(0);
}

// ================================= SMatrix ===============================

void cThickPlate :: SMatrix(cVector &str, cMatrix &S)
{
  cout << "cThickPlate :: SMatrix not implemented!\n";
  exit(0);
}

// ================================= NMatrix ===============================

void cThickPlate :: NMatrix(int nn, double *shpfunc, cMatrix &N)
{
  for (int i = 0; i < nn; i++)
  {
    N[0][3*i  ] = shpfunc[i];
    N[0][3*i+1] = 0.0;
    N[0][3*i+2] = 0.0;

    N[1][3*i  ] = 0.0;
    N[1][3*i+1] = shpfunc[i];
    N[1][3*i+2] = 0.0;

    N[2][3*i  ] = 0.0;
    N[2][3*i+1] = 0.0;
    N[2][3*i+2] = shpfunc[i];
  }
}

// ============================== ExpandStress =============================

void cThickPlate :: ExpandStress(cVector &sig, cVector &expsig)
{
  expsig[0] = sig[0];  // Sxx
  expsig[1] = sig[1];  // Syy
  expsig[2] = 0.0;     // Szz
  expsig[3] = sig[2];  // Sxy
  expsig[4] = sig[3];  // Sxz
  expsig[5] = sig[4];  // Syz
}


// -------------------------------------------------------------------------
// Class cHSDTPlate:
// -------------------------------------------------------------------------

// =========================== cHSDTPlate ===============================

cHSDTPlate :: cHSDTPlate(void)
{
  Type = HSDT_PLATE;
}

// =========================== ~cHSDTPlate ==============================

cHSDTPlate :: ~cHSDTPlate(void)
{
}

// ============================== GetActDir ================================

void cHSDTPlate :: GetActDir(int *dir)
{
  dir[0] = 1;   // u
  dir[1] = 1;   // v
  dir[2] = 1;   // w
  dir[3] = 1;   // rx
  dir[4] = 1;   // ry
  dir[5] = 0;   // rz
}

// ============================== GetStrLabels =============================

void cHSDTPlate :: GetStrLabels(int *label)
{
  label[0] = FORCE_X;
  label[1] = FORCE_Y;
  label[2] = FORCE_XY;
  label[3] = MOMENT_X;
  label[4] = MOMENT_Y;
  label[5] = MOMENT_XY;
  label[6] = P_XX;
  label[7] = P_YY;
  label[8] = P_XY;
  label[9] = R_XZ;
  label[10] = R_YZ;
}

// ============================== PrincStress ==============================

void cHSDTPlate :: PrincStress(cVector &str, cVector &prs)
{
  cout << "cHSDTPlate :: PrincStress not implemented!\n";
}

// ================================ QMatrix ================================

void cHSDTPlate :: QMatrix(double *param, cMatrix &Q)
{
  double E  = param[0];
  double nu = param[1];
  double G  = 0.5*E/(1.0 + nu);
  double d  = 1.0/(1.0 - nu*nu);

  Q.Zero( );

  // Plane-stress terms.

  Q[0][0] = d*E;
  Q[1][1] = d*E;
  Q[0][1] = Q[1][0] = d*nu*E;
  Q[2][2] = G;

  // Transverse shear terms.

  Q[3][3] = G;
  Q[4][4] = G;
}

// ============================== QMatrixOrtho =============================

void cHSDTPlate :: QMatrixOrtho(double *param, cMatrix &Q)
{
  cout << "cHSDTPlate :: MatrixOrtho not implemented!\n";
}

// ================================= DMatrix ===============================

void cHSDTPlate :: DMatrix(double *param, cMatrix &D)
{
  cout << "cHSDTPlate :: DMatrix not implemented!\n";
}

// ================================ TMatrix ================================

void cHSDTPlate :: TMatrix(double angle, cMatrix &T)
{
  cout << "cHSDTPlate :: TMatrix not implemented!\n";
}

// ================================= BMatrix ===============================

void cHSDTPlate :: BMatrix(double thk, int nn, double *shpfunc, double *mapfunc,
                              sNodeCoord *coord, sNodeCoord *drvshp,
                              sNodeDrv *drvshp2, cMatrix &B)
{
  // Assembly [B].

  for (int i = 0; i < nn; i++)
  {
    B[0][5*i  ] = drvshp[i].x;
    B[0][5*i+1] = 0.0;
    B[0][5*i+2] = 0.0;
    B[0][5*i+3] = 0.0;
    B[0][5*i+4] = 0.0;

    B[1][5*i  ] = 0.0;
    B[1][5*i+1] = drvshp[i].y;
    B[1][5*i+2] = 0.0;
    B[1][5*i+3] = 0.0;
    B[1][5*i+4] = 0.0;

    B[2][5*i  ] = drvshp[i].y;
    B[2][5*i+1] = drvshp[i].x;
    B[2][5*i+2] = 0.0;
    B[2][5*i+3] = 0.0;
    B[2][5*i+4] = 0.0;

    B[3][5*i  ] = 0.0;
    B[3][5*i+1] = 0.0;
    B[3][5*i+2] = -drvshp2[i].xx;
    B[3][5*i+3] = 0.0;
    B[3][5*i+4] = 0.0;

    B[4][5*i  ] = 0.0;
    B[4][5*i+1] = 0.0;
    B[4][5*i+2] = -drvshp2[i].yy;
    B[4][5*i+3] = 0.0;
    B[4][5*i+4] = 0.0;

    B[5][5*i  ] = 0.0;
    B[5][5*i+1] = 0.0;
    B[5][5*i+2] = -2.0*drvshp2[i].xy; // drv2shp.xy
    B[5][5*i+3] = 0.0;
    B[5][5*i+4] = 0.0;

    B[6][5*i  ] = 0.0;
    B[6][5*i+1] = 0.0;
    B[6][5*i+2] = drvshp2[i].xx; // drv2shp.xx
    B[6][5*i+3] = -drvshp[i].x;
    B[6][5*i+4] = 0.0;

    B[7][5*i  ] = 0.0;
    B[7][5*i+1] = 0.0;
    B[7][5*i+2] = drvshp2[i].yy; // drv2shp.yy
    B[7][5*i+3] = 0.0;
    B[7][5*i+4] = -drvshp[i].y;

    B[8][5*i  ] = 0.0;
    B[8][5*i+1] = 0.0;
    B[8][5*i+2] = 2.0*drvshp2[i].xy; // drv2shp.xy
    B[8][5*i+3] = -drvshp[i].y;
    B[8][5*i+4] = -drvshp[i].x;

    B[9][5*i  ] = 0.0;
    B[9][5*i+1] = 0.0;
    B[9][5*i+2] = drvshp[i].x;
    B[9][5*i+3] = -shpfunc[i];
    B[9][5*i+4] = 0.0;

    B[10][5*i  ] = 0.0;
    B[10][5*i+1] = 0.0;
    B[10][5*i+2] = drvshp[i].y;
    B[10][5*i+3] = 0.0;
    B[10][5*i+4] = -shpfunc[i];
  }
}

// ================================ BlMatrix ===============================

void cHSDTPlate :: BlMatrix(int nn, double *shpfunc, double *mapfunc,
                            cVector &u, sNodeCoord *coord,
                            sNodeCoord *drvshp, cMatrix &Bl)
{
  // Compute the displacement derivatives w,x and w,y.

  int i;
  double wx = 0;
  double wy = 0;
  for (i = 0; i < nn; i++)
  {
    double w = u[5*i+2];
    wx += w*drvshp[i].x;
    wy += w*drvshp[i].y;
  }

  // Assembly [Bl].

  for (i = 0; i < nn; i++)
  {
    Bl[0][5*i  ] = 0.0;
    Bl[0][5*i+1] = 0.0;
    Bl[0][5*i+2] = wx*drvshp[i].x;
    Bl[0][5*i+3] = 0.0;
    Bl[0][5*i+4] = 0.0;

    Bl[1][5*i  ] = 0.0;
    Bl[1][5*i+1] = 0.0;
    Bl[1][5*i+2] = wy*drvshp[i].y;
    Bl[1][5*i+3] = 0.0;
    Bl[1][5*i+4] = 0.0;

    Bl[2][5*i  ] = 0.0;
    Bl[2][5*i+1] = 0.0;
    Bl[2][5*i+2] = wy*drvshp[i].x + wx*drvshp[i].y;
    Bl[2][5*i+3] = 0.0;
    Bl[2][5*i+4] = 0.0;

    Bl[3][5*i  ] = 0.0;
    Bl[3][5*i+1] = 0.0;
    Bl[3][5*i+2] = 0.0;
    Bl[3][5*i+3] = 0.0;
    Bl[3][5*i+4] = 0.0;

    Bl[4][5*i  ] = 0.0;
    Bl[4][5*i+1] = 0.0;
    Bl[4][5*i+2] = 0.0;
    Bl[4][5*i+3] = 0.0;
    Bl[4][5*i+4] = 0.0;

    Bl[5][5*i  ] = 0.0;
    Bl[5][5*i+1] = 0.0;
    Bl[5][5*i+2] = 0.0;
    Bl[5][5*i+3] = 0.0;
    Bl[5][5*i+4] = 0.0;

    Bl[6][5*i  ] = 0.0;
    Bl[6][5*i+1] = 0.0;
    Bl[6][5*i+2] = 0.0;
    Bl[6][5*i+3] = 0.0;
    Bl[6][5*i+4] = 0.0;

    Bl[7][5*i  ] = 0.0;
    Bl[7][5*i+1] = 0.0;
    Bl[7][5*i+2] = 0.0;
    Bl[7][5*i+3] = 0.0;
    Bl[7][5*i+4] = 0.0;

    Bl[8][5*i  ] = 0.0;
    Bl[8][5*i+1] = 0.0;
    Bl[8][5*i+2] = 0.0;
    Bl[8][5*i+3] = 0.0;
    Bl[8][5*i+4] = 0.0;

    Bl[9][5*i  ] = 0.0;
    Bl[9][5*i+1] = 0.0;
    Bl[9][5*i+2] = 0.0;
    Bl[9][5*i+3] = 0.0;
    Bl[9][5*i+4] = 0.0;

    Bl[10][5*i  ] = 0.0;
    Bl[10][5*i+1] = 0.0;
    Bl[10][5*i+2] = 0.0;
    Bl[10][5*i+3] = 0.0;
    Bl[10][5*i+4] = 0.0;
  }
}

// =============================== BnlMatrix ===============================

void cHSDTPlate :: BnlMatrix(int nn, double *shpfunc, double *mapfunc,
                             sNodeCoord *coord, sNodeCoord *drvshp,
                             cMatrix &Bnl)
{
  for (int i = 0; i < nn; i++)
  {
    Bnl[0][5*i  ] = 0.0;
    Bnl[0][5*i+1] = 0.0;
    Bnl[0][5*i+2] = drvshp[i].x;
    Bnl[0][5*i+3] = 0.0;
    Bnl[0][5*i+4] = 0.0;

    Bnl[1][5*i  ] = 0.0;
    Bnl[1][5*i+1] = 0.0;
    Bnl[1][5*i+2] = drvshp[i].y;
    Bnl[1][5*i+3] = 0.0;
    Bnl[1][5*i+4] = 0.0;
  }
}

// ================================= SMatrix ===============================

void cHSDTPlate :: SMatrix(cVector &str, cMatrix &S)
{
  S.Zero( );
  S[0][0] = str[0];
  S[1][0] = str[2];
  S[0][1] = str[2];
  S[1][1] = str[1];
}

// ================================= NNMatrix ===============================

void cHSDTPlate :: NMatrix(int nn, double *shpfunc, sNodeCoord *drvshp, cMatrix &N)
{
  for (int i = 0; i < nn; i++)
  {
    N[0][5*i  ] = shpfunc[i];
    N[0][5*i+1] = 0.0;
    N[0][5*i+2] = 0.0;
    N[0][5*i+3] = 0.0;
    N[0][5*i+4] = 0.0;

    N[1][5*i  ] = 0.0;
    N[1][5*i+1] = shpfunc[i];
    N[1][5*i+2] = 0.0;
    N[1][5*i+3] = 0.0;
    N[1][5*i+4] = 0.0;

    N[2][5*i  ] = 0.0;
    N[2][5*i+1] = 0.0;
    N[2][5*i+2] = shpfunc[i];
    N[2][5*i+3] = 0.0;
    N[2][5*i+4] = 0.0;

    N[3][5*i  ] = 0.0;
    N[3][5*i+1] = 0.0;
    N[3][5*i+2] = 0.0;
    N[3][5*i+3] = shpfunc[i];
    N[3][5*i+4] = 0.0;

    N[4][5*i  ] = 0.0;
    N[4][5*i+1] = 0.0;
    N[4][5*i+2] = 0.0;
    N[4][5*i+3] = 0.0;
    N[4][5*i+4] = shpfunc[i];

    N[5][5*i  ] = 0.0;
    N[5][5*i+1] = 0.0;
    N[5][5*i+2] = drvshp[i].x;
    N[5][5*i+3] = 0.0;
    N[5][5*i+4] = 0.0;

    N[6][5*i  ] = 0.0;
    N[6][5*i+1] = 0.0;
    N[6][5*i+2] = drvshp[i].y;
    N[6][5*i+3] = 0.0;
    N[6][5*i+4] = 0.0;
  }
}

// ============================== ExpandStress =============================

void cHSDTPlate :: ExpandStress(cVector &sig, cVector &expsig)
{
  expsig[0] = sig[0];  // Sxx
  expsig[1] = sig[1];  // Syy
  expsig[2] = 0.0;     // Szz
  expsig[3] = sig[2];  // Sxy
  expsig[4] = sig[3];  // Sxz
  expsig[5] = sig[4];  // Syz
}


// -------------------------------------------------------------------------
// Class cShallowShell:
// -------------------------------------------------------------------------

// =========================== cShallowShell ===============================

cShallowShell :: cShallowShell(void)
{
  Type = SHALLOW_SHELL;
}

// =========================== ~cShallowShell ==============================

cShallowShell :: ~cShallowShell(void)
{
}

// ============================== GetActDir ================================

void cShallowShell :: GetActDir(int *dir)
{
  dir[0] = 1;   // u
  dir[1] = 1;   // v
  dir[2] = 1;   // w
  dir[3] = 1;   // rx
  dir[4] = 1;   // ry
}

// ============================== GetStrLabels =============================

void cShallowShell :: GetStrLabels(int *label)
{
  label[0] = FORCE_X;
  label[1] = FORCE_Y;
  label[2] = FORCE_XY;
  label[3] = MOMENT_Y;
  label[4] = MOMENT_X;
  label[5] = MOMENT_XY;
  label[6] = FORCE_XZ;
  label[7] = FORCE_YZ;
}

// ============================== PrincStress ==============================

void cShallowShell :: PrincStress(cVector &str, cVector &prs)
{
  cout << "cShallowShell :: PrincStress not implemented!\n";
}

// ================================= QMatrix ===============================

void cShallowShell :: QMatrix(double *param, cMatrix &Q)
{
  double E  = param[0];
  double nu = param[1];
  double G  = 0.5*E/(1.0 + nu);
  double d  = 1.0/(1.0 - nu*nu);

  Q.Zero( );

  // Membrane terms.

  Q[0][0] = d*E;
  Q[1][1] = d*E;
  Q[0][1] = Q[1][0] = d*nu*E;
  Q[2][2] = G;

  // Transverse shear terms.

  Q[3][3] = Q[4][4] = G;
}

// ============================== QMatrixOrtho =============================

void cShallowShell :: QMatrixOrtho(double *param, cMatrix &Q)
{
  double E1   = param[0];
  double E2   = param[1];
//  double E3   = param[2];
  double Nu12 = param[3];
//  double Nu13 = param[4];
//  double Nu23 = param[5];
  double G12  = param[6];
  double G13  = param[7];
  double G23  = param[8];
  double Nu21 = Nu12*E2/E1;

  Q.Zero( );

  // Evaluate the local plane stress constitutive matrix (lQb).

  Q[0][0] = E1/(1.0 - Nu12*Nu21);
  Q[1][1] = E2/(1.0 - Nu12*Nu21);
  Q[0][1] = Q[1][0] = (Nu12*E2)/(1.0 - Nu12*Nu21);
  Q[2][2] = G12;

  // Evaluate the local shear constitutive matrix (lQs).

  Q[3][3] = G13;
  Q[4][4] = G23;
}

// ================================= DMatrix ===============================

void cShallowShell :: DMatrix(double *param, cMatrix &D)
{
  cout << "cShallowShell :: DMatrix not implemented!\n";
}

// ================================ TMatrix ================================

void cShallowShell :: TMatrix(double angle, cMatrix &T)
{
  T.Zero( );

  // Trigonometric constants.

  double c  = cos(angle);
  double s  = sin(angle);
  double c2 = c*c;
  double s2 = s*s;

  // Bending transformation matrix.

  T[0][0] = T[1][1] = c2;
  T[0][1] = T[1][0] = s2;
  T[0][2] =  c*s;
  T[1][2] = -c*s;
  T[2][0] = -2.0*c*s;
  T[2][1] =  2.0*c*s;
  T[2][2] = c2 - s2;

  // Shear transformation matrix.

  T[3][3] = T[4][4] = c;
  T[3][4] = -s;
  T[4][3] =  s;
}

// ================================= BMatrix ===============================

void cShallowShell :: BMatrix(int nn, double *shpfunc, double *mapfunc,
                              sNodeCoord *coord, sNodeCoord *drvshp,
                              cMatrix &B)
{
  // Compute the derivatives z,x and z,y.

  int i;
  double zx = 0.0;
  double zy = 0.0;
  for (i = 0; i < nn; i++)
  {
    zx += drvshp[i].x*coord[i].z;
    zy += drvshp[i].y*coord[i].z;
  }

  // Assembly [B].

  for (i = 0; i < nn; i++)
  {
    B[0][5*i  ] = drvshp[i].x;
    B[0][5*i+1] = 0.0;
    B[0][5*i+2] = zx*drvshp[i].x;
    B[0][5*i+3] = 0.0;
    B[0][5*i+4] = 0.0;

    B[1][5*i  ] = 0.0;
    B[1][5*i+1] = drvshp[i].y;
    B[1][5*i+2] = zy*drvshp[i].y;
    B[1][5*i+3] = 0.0;
    B[1][5*i+4] = 0.0;

    B[2][5*i  ] = drvshp[i].y;
    B[2][5*i+1] = drvshp[i].x;
    B[2][5*i+2] = zx*drvshp[i].y + zy*drvshp[i].x;
    B[2][5*i+3] = 0.0;
    B[2][5*i+4] = 0.0;

    B[3][5*i  ] = 0.0;
    B[3][5*i+1] = 0.0;
    B[3][5*i+2] = 0.0;
    B[3][5*i+3] = 0.0;
    B[3][5*i+4] = drvshp[i].x;

    B[4][5*i  ] = 0.0;
    B[4][5*i+1] = 0.0;
    B[4][5*i+2] = 0.0;
    B[4][5*i+3] = -drvshp[i].y;
    B[4][5*i+4] = 0.0;

    B[5][5*i  ] = 0.0;
    B[5][5*i+1] = 0.0;
    B[5][5*i+2] = 0.0;
    B[5][5*i+3] = -drvshp[i].x;
    B[5][5*i+4] =  drvshp[i].y;

    B[6][5*i  ] = 0.0;
    B[6][5*i+1] = 0.0;
    B[6][5*i+2] = drvshp[i].x;
    B[6][5*i+3] = 0.0;
    B[6][5*i+4] = shpfunc[i];

    B[7][5*i  ] = 0.0;
    B[7][5*i+1] = 0.0;
    B[7][5*i+2] =  drvshp[i].y;
    B[7][5*i+3] = -shpfunc[i];
    B[7][5*i+4] = 0.0;
  }
}

// ================================ BlMatrix ===============================

void cShallowShell :: BlMatrix(int nn, double *shpfunc, double *mapfunc,
                               cVector &u, sNodeCoord *coord,
                               sNodeCoord *drvshp, cMatrix &Bl)
{
  // Compute the displacement derivatives w,x and w,y.

  int i;
  double wx = 0;
  double wy = 0;
  for (i = 0; i < nn; i++)
  {
    double w = u[5*i+2];
    wx += w*drvshp[i].x;
    wy += w*drvshp[i].y;
  }

  // Assembly [Bl].

  for (i = 0; i < nn; i++)
  {
    Bl[0][5*i  ] = 0.0;
    Bl[0][5*i+1] = 0.0;
    Bl[0][5*i+2] = wx*drvshp[i].x;
    Bl[0][5*i+3] = 0.0;
    Bl[0][5*i+4] = 0.0;

    Bl[1][5*i  ] = 0.0;
    Bl[1][5*i+1] = 0.0;
    Bl[1][5*i+2] = wy*drvshp[i].y;
    Bl[1][5*i+3] = 0.0;
    Bl[1][5*i+4] = 0.0;

    Bl[2][5*i  ] = 0.0;
    Bl[2][5*i+1] = 0.0;
    Bl[2][5*i+2] = wy*drvshp[i].x + wx*drvshp[i].y;
    Bl[2][5*i+3] = 0.0;
    Bl[2][5*i+4] = 0.0;

    Bl[3][5*i  ] = 0.0;
    Bl[3][5*i+1] = 0.0;
    Bl[3][5*i+2] = 0.0;
    Bl[3][5*i+3] = 0.0;
    Bl[3][5*i+4] = 0.0;

    Bl[4][5*i  ] = 0.0;
    Bl[4][5*i+1] = 0.0;
    Bl[4][5*i+2] = 0.0;
    Bl[4][5*i+3] = 0.0;
    Bl[4][5*i+4] = 0.0;

    Bl[5][5*i  ] = 0.0;
    Bl[5][5*i+1] = 0.0;
    Bl[5][5*i+2] = 0.0;
    Bl[5][5*i+3] = 0.0;
    Bl[5][5*i+4] = 0.0;

    Bl[6][5*i  ] = 0.0;
    Bl[6][5*i+1] = 0.0;
    Bl[6][5*i+2] = 0.0;
    Bl[6][5*i+3] = 0.0;
    Bl[6][5*i+4] = 0.0;

    Bl[7][5*i  ] = 0.0;
    Bl[7][5*i+1] = 0.0;
    Bl[7][5*i+2] = 0.0;
    Bl[7][5*i+3] = 0.0;
    Bl[7][5*i+4] = 0.0;
  }
}

// =============================== BnlMatrix ===============================

void cShallowShell :: BnlMatrix(int nn, double *shpfunc, double *mapfunc,
                                sNodeCoord *coord, sNodeCoord *drvshp,
                                cMatrix &Bnl)
{
  for (int i = 0; i < nn; i++)
  {
    Bnl[0][5*i  ] = 0.0;
    Bnl[0][5*i+1] = 0.0;
    Bnl[0][5*i+2] = drvshp[i].x;
    Bnl[0][5*i+3] = 0.0;
    Bnl[0][5*i+4] = 0.0;

    Bnl[1][5*i  ] = 0.0;
    Bnl[1][5*i+1] = 0.0;
    Bnl[1][5*i+2] = drvshp[i].y;
    Bnl[1][5*i+3] = 0.0;
    Bnl[1][5*i+4] = 0.0;
  }
}

// ================================= SMatrix ===============================

void cShallowShell :: SMatrix(cVector &str, cMatrix &S)
{
  S.Zero( );
  S[0][0] = str[0];
  S[1][0] = str[2];
  S[0][1] = str[2];
  S[1][1] = str[1];
}

// ================================= NMatrix ===============================

void cShallowShell :: NMatrix(int nn, double *shpfunc, cMatrix &N)
{
  for (int i = 0; i < nn; i++)
  {
    N[0][5*i  ] = shpfunc[i];
    N[0][5*i+1] = 0.0;
    N[0][5*i+2] = 0.0;
    N[0][5*i+3] = 0.0;
    N[0][5*i+4] = 0.0;

    N[1][5*i  ] = 0.0;
    N[1][5*i+1] = shpfunc[i];
    N[1][5*i+2] = 0.0;
    N[1][5*i+3] = 0.0;
    N[1][5*i+4] = 0.0;

    N[2][5*i  ] = 0.0;
    N[2][5*i+1] = 0.0;
    N[2][5*i+2] = shpfunc[i];
    N[2][5*i+3] = 0.0;
    N[2][5*i+4] = 0.0;

    N[3][5*i  ] = 0.0;
    N[3][5*i+1] = 0.0;
    N[3][5*i+2] = 0.0;
    N[3][5*i+3] = shpfunc[i];
    N[3][5*i+4] = 0.0;

    N[4][5*i  ] = 0.0;
    N[4][5*i+1] = 0.0;
    N[4][5*i+2] = 0.0;
    N[4][5*i+3] = 0.0;
    N[4][5*i+4] = shpfunc[i];
  }
}

// ============================== ExpandStress =============================

void cShallowShell :: ExpandStress(cVector &sig, cVector &expsig)
{
  expsig[0] = sig[0];  // Sxx
  expsig[1] = sig[1];  // Syy
  expsig[2] = 0.0;     // Szz
  expsig[3] = sig[2];  // Sxy
  expsig[4] = sig[3];  // Sxz
  expsig[5] = sig[4];  // Syz
}

// ============================== AlphaVec =============================

void cShallowShell :: AlphaVec(double *param, cVector &Alpha)
{
//  cout << "cShallowShell :: AlphaVec\n";
  double alpha = param[2];  // Coefficient of Thermal Expansion
//  cout << "alpha = " << alpha << endl;

  Alpha.Zero( );
  Alpha[0] = alpha;
  Alpha[1] = alpha;
}

// ============================== AlphaVecOrtho =============================

void cShallowShell :: AlphaVecOrtho(double *param, cVector &Alpha)
{
  double alpha1 = param[9];  // Coefficient of Thermal Expansion - dir 1
  double alpha2 = param[10]; // Coefficient of Thermal Expansion - dir 2 

  Alpha.Zero( );
  Alpha[0] = alpha1;
  Alpha[1] = alpha2;
}


// -------------------------------------------------------------------------
// Class cDonnellShell:
// -------------------------------------------------------------------------

// =========================== cDonnellShell ===============================

cDonnellShell :: cDonnellShell(void)
{
  Type = DONNELL_SHELL;
}

// ========================== ~cDonnellShell ===============================

cDonnellShell :: ~cDonnellShell(void)
{
}

// ============================== GetActDir ================================

void cDonnellShell :: GetActDir(int *dir)
{
  dir[0] = 1;   // u
  dir[1] = 1;   // v
  dir[2] = 1;   // w
  dir[3] = 1;   // rx
  dir[4] = 1;   // ry
  dir[5] = 0;   // rz
}

// ============================== GetStrLabels =============================

void cDonnellShell :: GetStrLabels(int *label)
{
  label[0] = FORCE_X;
  label[1] = FORCE_Y;
  label[2] = FORCE_XY;
  label[3] = MOMENT_Y;
  label[4] = MOMENT_X;
  label[5] = MOMENT_XY;
  label[6] = FORCE_XZ;
  label[7] = FORCE_YZ;
}

// ============================== PrincStress ==============================

void cDonnellShell :: PrincStress(cVector &str, cVector &prs)
{
  cout << "cDonnellShell :: PrincStress not implemented!\n";
}

// ================================ QMatrix ================================

void cDonnellShell :: QMatrix(double *param, cMatrix &Q)
{
  double E  = param[0];
  double nu = param[1];
  double G  = 0.5*E/(1.0 + nu);
  double d  = 1.0/(1.0 - nu*nu);

  Q.Zero( );

  // Membrane terms.

  Q[0][0] = d*E;
  Q[1][1] = d*E;
  Q[0][1] = Q[1][0] = d*nu*E;
  Q[2][2] = G;

  // Transverse shear terms.

  Q[3][3] = Q[4][4] = G;
}

// ============================== QMatrixOrtho =============================

void cDonnellShell :: QMatrixOrtho(double *param, cMatrix &Q)
{
  double E1   = param[0];
  double E2   = param[1];
//  double E3   = param[2];
  double Nu12 = param[3];
//  double Nu13 = param[4];
//  double Nu23 = param[5];
  double G12  = param[6];
  double G13  = param[7];
  double G23  = param[8];
  double Nu21 = Nu12*E2/E1;

  Q.Zero( );

  // Evaluate the local bending constitutive matrix (lQb).

  Q[0][0] = E1/(1.0 - Nu12*Nu21);
  Q[1][1] = E2/(1.0 - Nu12*Nu21);
  Q[0][1] = Q[1][0] = (Nu12*E2)/(1.0 - Nu12*Nu21);
  Q[2][2] = G12;

  // Evaluate the local shear constitutive matrix (lQs).

  Q[3][3] = G13;
  Q[4][4] = G23;
}

// ================================= DMatrix ===============================

void cDonnellShell :: DMatrix(double *param, cMatrix &D)
{
  cout << "cDonnellShell :: DMatrix not implemented!\n";
}

// ================================ TMatrix ================================

void cDonnellShell :: TMatrix(double angle, cMatrix &T)
{
  T.Zero( );

  // Trigonometric constants.

  double c  = cos(angle);
  double s  = sin(angle);
  double c2 = c*c;
  double s2 = s*s;

  // Bending transformation matrix.

  T[0][0] = T[1][1] = c2;
  T[0][1] = T[1][0] = s2;
  T[0][2] =  c*s;
  T[1][2] = -c*s;
  T[2][0] = -2.0*c*s;
  T[2][1] =  2.0*c*s;
  T[2][2] = c2 - s2;

  // Shear transformation matrix.

  T[3][3] = T[4][4] = c;
  T[3][4] = -s;
  T[4][3] =  s;
}

// ================================= BMatrix ===============================

void cDonnellShell :: BMatrix(int nn, double *shpfunc, double *mapfunc,
                              sNodeCoord *coord, sNodeCoord *drvshp,
                              cMatrix &B)
{
  // Compute the derivatives z,x and z,y.

  int i;
  double r = 0.0;
  for (i = 0; i < nn; i++) r += mapfunc[i]*coord[i].z;

  // Assembly [B].

  for (i = 0; i < nn; i++)
  {
    B[0][5*i  ] = drvshp[i].x;
    B[0][5*i+1] = 0.0;
    B[0][5*i+2] = 0.0;
    B[0][5*i+3] = 0.0;
    B[0][5*i+4] = 0.0;

    B[1][5*i  ] = 0.0;
    B[1][5*i+1] = drvshp[i].y;
    B[1][5*i+2] = shpfunc[i]/r;
    B[1][5*i+3] = 0.0;
    B[1][5*i+4] = 0.0;

    B[2][5*i  ] = drvshp[i].y;
    B[2][5*i+1] = drvshp[i].x;
    B[2][5*i+2] = 0.0;
    B[2][5*i+3] = 0.0;
    B[2][5*i+4] = 0.0;

    B[3][5*i  ] = 0.0;
    B[3][5*i+1] = 0.0;
    B[3][5*i+2] = 0.0;
    B[3][5*i+3] = 0.0;
    B[3][5*i+4] = -drvshp[i].x;

    B[4][5*i  ] = 0.0;
    B[4][5*i+1] = 0.0;
    B[4][5*i+2] = 0.0;
    B[4][5*i+3] = drvshp[i].y;
    B[4][5*i+4] = 0.0;

    B[5][5*i  ] = 0.0;
    B[5][5*i+1] = 0.0;
    B[5][5*i+2] = 0.0;
    B[5][5*i+3] = -drvshp[i].x;
    B[5][5*i+4] =  drvshp[i].y;

    B[6][5*i  ] = 0.0;
    B[6][5*i+1] = 0.0;
    B[6][5*i+2] = drvshp[i].x;
    B[6][5*i+3] = 0.0;
    B[6][5*i+4] = shpfunc[i];

    B[7][5*i  ] = 0.0;
    B[7][5*i+1] = -shpfunc[i]/r;
    B[7][5*i+2] =  drvshp[i].y;
    B[7][5*i+3] = -shpfunc[i];
    B[7][5*i+4] = 0.0;
  }
}

// ================================ BlMatrix ===============================

void cDonnellShell :: BlMatrix(int nn, double *shpfunc, double *mapfunc,
                               cVector &u, sNodeCoord *coord,
                               sNodeCoord *drvshp, cMatrix &Bl)
{
  // Compute the displacement derivatives w,x and w,y.

  int i;
  double wx = 0;
  double wy = 0;
  for (i = 0; i < nn; i++)
  {
    double w = u[5*i+2];
    wx += w*drvshp[i].x;
    wy += w*drvshp[i].y;
  }

  // Assembly [Bl].

  for (i = 0; i < nn; i++)
  {
    Bl[0][5*i  ] = 0.0;
    Bl[0][5*i+1] = 0.0;
    Bl[0][5*i+2] = wx*drvshp[i].x;
    Bl[0][5*i+3] = 0.0;
    Bl[0][5*i+4] = 0.0;

    Bl[1][5*i  ] = 0.0;
    Bl[1][5*i+1] = 0.0;
    Bl[1][5*i+2] = wy*drvshp[i].y;
    Bl[1][5*i+3] = 0.0;
    Bl[1][5*i+4] = 0.0;

    Bl[2][5*i  ] = 0.0;
    Bl[2][5*i+1] = 0.0;
    Bl[2][5*i+2] = wy*drvshp[i].x + wx*drvshp[i].y;
    Bl[2][5*i+3] = 0.0;
    Bl[2][5*i+4] = 0.0;

    Bl[3][5*i  ] = 0.0;
    Bl[3][5*i+1] = 0.0;
    Bl[3][5*i+2] = 0.0;
    Bl[3][5*i+3] = 0.0;
    Bl[3][5*i+4] = 0.0;

    Bl[4][5*i  ] = 0.0;
    Bl[4][5*i+1] = 0.0;
    Bl[4][5*i+2] = 0.0;
    Bl[4][5*i+3] = 0.0;
    Bl[4][5*i+4] = 0.0;

    Bl[5][5*i  ] = 0.0;
    Bl[5][5*i+1] = 0.0;
    Bl[5][5*i+2] = 0.0;
    Bl[5][5*i+3] = 0.0;
    Bl[5][5*i+4] = 0.0;

    Bl[6][5*i  ] = 0.0;
    Bl[6][5*i+1] = 0.0;
    Bl[6][5*i+2] = 0.0;
    Bl[6][5*i+3] = 0.0;
    Bl[6][5*i+4] = 0.0;

    Bl[7][5*i  ] = 0.0;
    Bl[7][5*i+1] = 0.0;
    Bl[7][5*i+2] = 0.0;
    Bl[7][5*i+3] = 0.0;
    Bl[7][5*i+4] = 0.0;
  }
}

// =============================== BnlMatrix ===============================

void cDonnellShell :: BnlMatrix(int nn, double *shpfunc, double *mapfunc,
                                sNodeCoord *coord, sNodeCoord *drvshp,
                                cMatrix &Bnl)
{
  for (int i = 0; i < nn; i++)
  {
    Bnl[0][5*i  ] = 0.0;
    Bnl[0][5*i+1] = 0.0;
    Bnl[0][5*i+2] = drvshp[i].x;
    Bnl[0][5*i+3] = 0.0;
    Bnl[0][5*i+4] = 0.0;

    Bnl[1][5*i  ] = 0.0;
    Bnl[1][5*i+1] = 0.0;
    Bnl[1][5*i+2] = drvshp[i].y;
    Bnl[1][5*i+3] = 0.0;
    Bnl[1][5*i+4] = 0.0;
  }
}

// ================================= SMatrix ===============================

void cDonnellShell :: SMatrix(cVector &str, cMatrix &S)
{
  S.Zero( );
  S[0][0] = str[0];
  S[1][0] = str[2];
  S[0][1] = str[2];
  S[1][1] = str[1];
}

// ================================= NMatrix ===============================

void cDonnellShell :: NMatrix(int nn, double *shpfunc, cMatrix &N)
{
  for (int i = 0; i < nn; i++)
  {
    N[0][5*i  ] = shpfunc[i];
    N[0][5*i+1] = 0.0;
    N[0][5*i+2] = 0.0;
    N[0][5*i+3] = 0.0;
    N[0][5*i+4] = 0.0;

    N[1][5*i  ] = 0.0;
    N[1][5*i+1] = shpfunc[i];
    N[1][5*i+2] = 0.0;
    N[1][5*i+3] = 0.0;
    N[1][5*i+4] = 0.0;

    N[2][5*i  ] = 0.0;
    N[2][5*i+1] = 0.0;
    N[2][5*i+2] = shpfunc[i];
    N[2][5*i+3] = 0.0;
    N[2][5*i+4] = 0.0;

    N[3][5*i  ] = 0.0;
    N[3][5*i+1] = 0.0;
    N[3][5*i+2] = 0.0;
    N[3][5*i+3] = shpfunc[i];
    N[3][5*i+4] = 0.0;

    N[4][5*i  ] = 0.0;
    N[4][5*i+1] = 0.0;
    N[4][5*i+2] = 0.0;
    N[4][5*i+3] = 0.0;
    N[4][5*i+4] = shpfunc[i];
  }
}

// ============================== ExpandStress =============================

void cDonnellShell :: ExpandStress(cVector &sig, cVector &expsig)
{
  expsig[0] = sig[0];  // Sxx
  expsig[1] = sig[1];  // Syy
  expsig[2] = 0.0;     // Szz
  expsig[3] = sig[2];  // Sxy
  expsig[4] = sig[3];  // Sxz
  expsig[5] = sig[4];  // Syz
}


// -------------------------------------------------------------------------
// Class cShell:
// -------------------------------------------------------------------------

// =============================== cShell ==================================

cShell :: cShell(void)
{
  Type = SHELL;
}

// =============================== ~cShell =================================

cShell :: ~cShell(void)
{
}

// ============================== GetActDir ================================

void cShell :: GetActDir(int *dir)
{
  dir[0] = 1;   // u
  dir[1] = 1;   // v
  dir[2] = 1;   // w
  dir[3] = 1;   // rx
  dir[4] = 1;   // ry
}

// ============================== GetStrLabels =============================

void cShell :: GetStrLabels(int *label)
{
  label[0] = STRESS_XX;
  label[1] = STRESS_YY;
  label[2] = STRESS_ZZ;
  label[3] = STRESS_XY;
  label[4] = STRESS_XZ;
  label[5] = STRESS_YZ;
}

// ============================== GetNodStrLabel ===========================

void cShell :: GetNodStrLabels(int *label)
{
  label[0] = FORCE_X;
  label[1] = FORCE_Y;
  label[2] = FORCE_XY;
  label[3] = MOMENT_Y;
  label[4] = MOMENT_X;
  label[5] = MOMENT_XY;
  label[6] = FORCE_XZ;
  label[7] = FORCE_YZ;
}

// ============================== PrincStress ==============================

void cShell :: PrincStress(cVector &str, cVector &prs)
{
  cout << "cShell :: PrincStress not implemented!\n";
}

// ================================ QMatrix ================================

void cShell :: QMatrix(double *param, cMatrix &Q)
{
  double E  = param[0];
  double nu = param[1];
  double G  = 0.5*E/(1.0 + nu);
  double d  = 1.0/(1.0 - nu*nu);

  Q.Zero( );

  Q[0][0] = d*E;
  Q[1][1] = d*E;
  Q[0][1] = Q[1][0] = d*nu*E;
  Q[3][3] = G;
  Q[4][4] = Q[5][5] = 5.0/6.0*G; // Evandro: ks should be in SecAnalysis
}

// ============================== QMatrixOrtho =============================

void cShell :: QMatrixOrtho(double *param, cMatrix &Q)
{
  double E1   = param[0];
  double E2   = param[1];
//  double E3   = param[2];
  double Nu12 = param[3];
//  double Nu13 = param[4];
//  double Nu23 = param[5];
  double G12  = param[6];
  double G13  = param[7];
  double G23  = param[8];
  double Nu21 = Nu12*E2/E1;

  Q.Zero( );

  // Evaluate the local bending constitutive matrix (lQb).

  Q[0][0] = E1/(1.0 - Nu12*Nu21);
  Q[1][1] = E2/(1.0 - Nu12*Nu21);
  Q[0][1] = Q[1][0] = (Nu12*E2)/(1.0 - Nu12*Nu21);
  Q[3][3] = G12;

  // Evaluate the local shear constitutive matrix (lQs).

  Q[4][4] = 5.0/6.0*G13; // Evandro: ks should be in SecAnalysis
  Q[5][5] = 5.0/6.0*G23; // Evandro: ks should be in SecAnalysis
}

// ================================= DMatrix ===============================

void cShell :: DMatrix(double *param, cMatrix &D)
{
  cout << "cShell :: DMatrix not implemented!\n";
}

// ================================ TMatrix ================================

void cShell :: TMatrix(cVector &e1, cVector &e2, cVector &e3, cMatrix &T)
{
  T.Zero( );

  // Trigonometric constants.

  double l1 = e1[0];
  double m1 = e1[1];
  double n1 = e1[2];

  double l2 = e2[0];
  double m2 = e2[1];
  double n2 = e2[2];

  double l3 = e3[0];
  double m3 = e3[1];
  double n3 = e3[2];

  // Transformation Matrix.

  T[0][0] = l1*l1;
  T[0][1] = m1*m1;
  T[0][2] = n1*n1;
  T[0][3] = l1*m1;
  T[0][4] = n1*l1;
  T[0][5] = m1*n1;

  T[1][0] = l2*l2;
  T[1][1] = m2*m2;
  T[1][2] = n2*n2;
  T[1][3] = l2*m2;
  T[1][4] = n2*l2;
  T[1][5] = m2*n2;

  T[2][0] = l3*l3;
  T[2][1] = m3*m3;
  T[2][2] = n3*n3;
  T[2][3] = l3*m3;
  T[2][4] = n3*l3;
  T[2][5] = m3*n3;

  T[3][0] = 2*l1*l2;
  T[3][1] = 2*m1*m2;
  T[3][2] = 2*n1*n2;
  T[3][3] = l1*m2+l2*m1;
  T[3][4] = n1*l2+n2*l1;
  T[3][5] = m1*n2+m2*n1;

  T[4][0] = 2*l3*l1;
  T[4][1] = 2*m3*m1;
  T[4][2] = 2*n3*n1;
  T[4][3] = l3*m1+l1*m3;
  T[4][4] = n3*l1+n1*l3;
  T[4][5] = m3*n1+m1*n3;

  T[5][0] = 2*l2*l3;
  T[5][1] = 2*m2*m3;
  T[5][2] = 2*n2*n3;
  T[5][3] = l2*m3+l3*m2;
  T[5][4] = n2*l3+n3*l2;
  T[5][5] = m2*n3+m3*n2;

/*
  T[0][0] = l1 * l1;
  T[0][1] = m1 * m1;
  T[0][2] = n1 * n1;
  T[0][3] = l1 * m1;
  T[0][4] = n1 * l1;
  T[0][5] = m1 * n1;
  
  T[1][0] = l2 * l2;
  T[1][1] = m2 * m2;
  T[1][2] = n2 * n2;
  T[1][3] = l2 * m2;                                                                
  T[1][4] = n2 * l2;                                                                
  T[1][5] = m2 * n2;                                                                
  
  T[2][0] = 2.0 * l1 * l2;                                                          
  T[2][1] = 2.0 * m1 * m2;                                                          
  T[2][2] = 2.0 * n1 * n2; 
  T[2][3] = l1 * m2 + l2 * m1;                                                      
  T[2][4] = n1 * l2 + n2 * l1;                                                      
  T[2][5] = m1 * n2 + m2 * n1;                                                      
  
  T[3][0] = 2.0 * l3 * l1;                       
  T[3][1] = 2.0 * m3 * m1;
  T[3][2] = 2.0 * n3 * n1;
  T[3][3] = l3 * m1 + l1 * m3;
  T[3][4] = n3 * l1 + n1 * l3;
  T[3][5] = m3 * n1 + m1 * n3;
 
  T[4][0] = 2.0 * l2 * l3;
  T[4][1] = 2.0 * m2 * m3;
  T[4][2] = 2.0 * n2 * n3;
  T[4][3] = l2 * m3 + l3 * m2;
  T[4][4] = n2 * l3 + n3 * l2;
  T[4][5] = m2 * n3 + m3 * n2;
  */
}

void cShell :: TMatrix(double angle, cMatrix &T)
{
  T.Zero( );

  // Trigonometric constants.

  double c  = cos(angle);
  double s  = sin(angle);
  double c2 = c*c;
  double s2 = s*s;

  // Bending transformation matrix.

  T[0][0] = T[1][1] = c2;
  T[0][1] = T[1][0] = s2;
  T[0][2] =  c*s;
  T[1][2] = -c*s;
  T[2][0] = -2.0*c*s;
  T[2][1] =  2.0*c*s;
  T[2][2] = c2 - s2;

  // Shear transformation matrix.

  T[3][3] = T[4][4] = c;
  T[3][4] = -s;
  T[4][3] =  s;
}

// ================================= BMatrix ===============================

void cShell :: BMatrix(int nn, double *shpfunc, double *mapfunc,
                       sNodeCoord *coord, sNodeCoord *drvshp, cMatrix &B)
{
  // Get the through-thickness coordinate and its derivatives.

  double zeta = shpfunc[nn];
  double dzetadx = drvshp[nn].x;
  double dzetady = drvshp[nn].y;
  double dzetadz = drvshp[nn].z;

  // Assembly [B].

  for (int i = 0; i < nn; i++)
  {
    double thk2 = *(coord[i].thk)/2.0;
    double dHdx = drvshp[i].x*zeta + shpfunc[i]*dzetadx;
    double dHdy = drvshp[i].y*zeta + shpfunc[i]*dzetady;
    double dHdz = drvshp[i].z*zeta + shpfunc[i]*dzetadz;

    /*
    double pax = -coord[i].axes->V2[0]*thk2;
    double pay = -coord[i].axes->V2[1]*thk2;
    double paz = -coord[i].axes->V2[2]*thk2;
    double pbx =  coord[i].axes->V1[0]*thk2;
    double pby =  coord[i].axes->V1[1]*thk2;
    double pbz =  coord[i].axes->V1[2]*thk2;
    */

    double pax = (*coord[i].pa)[0];
    double pay = (*coord[i].pa)[1];
    double paz = (*coord[i].pa)[2];
    double pbx = (*coord[i].pb)[0];
    double pby = (*coord[i].pb)[1];
    double pbz = (*coord[i].pb)[2];

    B[0][5*i  ] = drvshp[i].x;
    B[0][5*i+1] = 0.0;
    B[0][5*i+2] = 0.0;
    B[0][5*i+3] = dHdx*pax;
    B[0][5*i+4] = dHdx*pbx;

    B[1][5*i  ] = 0.0;
    B[1][5*i+1] = drvshp[i].y;
    B[1][5*i+2] = 0.0;
    B[1][5*i+3] = dHdy*pay;
    B[1][5*i+4] = dHdy*pby;

    B[2][5*i  ] = 0.0;
    B[2][5*i+1] = 0.0;
    B[2][5*i+2] = drvshp[i].z;
    B[2][5*i+3] = dHdz*paz;
    B[2][5*i+4] = dHdz*pbz;

    B[3][5*i  ] = drvshp[i].y;
    B[3][5*i+1] = drvshp[i].x;
    B[3][5*i+2] = 0.0;
    B[3][5*i+3] = dHdx*pay + dHdy*pax;
    B[3][5*i+4] = dHdx*pby + dHdy*pbx;

    B[4][5*i  ] = drvshp[i].z;
    B[4][5*i+1] = 0.0;
    B[4][5*i+2] = drvshp[i].x;
    B[4][5*i+3] = dHdx*paz + dHdz*pax;
    B[4][5*i+4] = dHdx*pbz + dHdz*pbx;

    B[5][5*i  ] = 0.0;
    B[5][5*i+1] = drvshp[i].z;
    B[5][5*i+2] = drvshp[i].y;
    B[5][5*i+3] = dHdy*paz + dHdz*pay;
    B[5][5*i+4] = dHdy*pbz + dHdz*pby;
  }
}

// ================================ BlMatrix ===============================

void cShell :: BlMatrix(int nn, double *shpfunc, double *mapfunc,
                        cVector &disp, sNodeCoord *coord, sNodeCoord *drvshp,
                        cMatrix &Bl)
{
  // Get the through-thickness coordinate and its derivatives.

  double zeta = shpfunc[nn];
  double dzetadx = drvshp[nn].x;
  double dzetady = drvshp[nn].y;
  double dzetadz = drvshp[nn].z;

  // Compute the displacement derivatives.

  double ux = 0.0;
  double uy = 0.0;
  double uz = 0.0;
  double vx = 0.0;
  double vy = 0.0;
  double vz = 0.0;
  double wx = 0.0;
  double wy = 0.0;
  double wz = 0.0;
  for (int i = 0; i < nn; i++)
  {
    double u = disp[5*i];
    double v = disp[5*i+1];
    double w = disp[5*i+2];

    double thk2 = *(coord[i].thk)/2.0;
    double dHdx = drvshp[i].x*zeta + shpfunc[i]*dzetadx;
    double dHdy = drvshp[i].y*zeta + shpfunc[i]*dzetady;
    double dHdz = drvshp[i].z*zeta + shpfunc[i]*dzetadz;
//    cout << scientific << setprecision(6) << showpos;
//    cout << "{dH} = ";
//    cout << dHdx << " " << dHdy << " " << dHdz << "\n";

    double px = (coord[i].axes->V3[0] - coord[i].iaxes->V3[0])*thk2;
    double py = (coord[i].axes->V3[1] - coord[i].iaxes->V3[1])*thk2;
    double pz = (coord[i].axes->V3[2] - coord[i].iaxes->V3[2])*thk2;
//    cout << "{p} = ";
//    cout << px << " " << py << " " << pz << "\n";
//    cout << "{v1} = ";
//    coord[i].axes->V1.Print();
//    cout << "{v2} = ";
//    coord[i].axes->V2.Print();
//    cout << "{v3} = ";
//    coord[i].axes->V3.Print();

    ux += drvshp[i].x*u + dHdx*px;
    uy += drvshp[i].y*u + dHdy*px;
    uz += drvshp[i].z*u + dHdz*px;

    vx += drvshp[i].x*v + dHdx*py;
    vy += drvshp[i].y*v + dHdy*py;
    vz += drvshp[i].z*v + dHdz*py;

    wx += drvshp[i].x*w + dHdx*pz;
    wy += drvshp[i].y*w + dHdy*pz;
    wz += drvshp[i].z*w + dHdz*pz;
  }
//  cout << "{beta} = " << scientific << setprecision(6) << showpos;
//  cout << ux << " " << uy << " " << uz << " ";
//  cout << vx << " " << vy << " " << vz << " ";
//  cout << wx << " " << wy << " " << wz << "\n";

  // Assembly [Bl].

  for (int i = 0; i < nn; i++)
  {
    double thk2 = *(coord[i].thk)/2.0;
    double dHdx = drvshp[i].x*zeta + shpfunc[i]*dzetadx;
    double dHdy = drvshp[i].y*zeta + shpfunc[i]*dzetady;
    double dHdz = drvshp[i].z*zeta + shpfunc[i]*dzetadz;

    /*
    double pax = -coord[i].axes->V2[0]*thk2;
    double pay = -coord[i].axes->V2[1]*thk2;
    double paz = -coord[i].axes->V2[2]*thk2;
    double pbx =  coord[i].axes->V1[0]*thk2;
    double pby =  coord[i].axes->V1[1]*thk2;
    double pbz =  coord[i].axes->V1[2]*thk2;
    */

    double pax = (*coord[i].pa)[0];
    double pay = (*coord[i].pa)[1];
    double paz = (*coord[i].pa)[2];
    double pbx = (*coord[i].pb)[0];
    double pby = (*coord[i].pb)[1];
    double pbz = (*coord[i].pb)[2];

    double phi11 = pax*ux + pay*vx + paz*wx;
    double phi12 = pax*uy + pay*vy + paz*wy;
    double phi13 = pax*uz + pay*vz + paz*wz;
    double phi21 = pbx*ux + pby*vx + pbz*wx;
    double phi22 = pbx*uy + pby*vy + pbz*wy;
    double phi23 = pbx*uz + pby*vz + pbz*wz;

    Bl[0][5*i  ] = drvshp[i].x*ux;
    Bl[0][5*i+1] = drvshp[i].x*vx;
    Bl[0][5*i+2] = drvshp[i].x*wx;
    Bl[0][5*i+3] = dHdx*phi11;
    Bl[0][5*i+4] = dHdx*phi21;

    Bl[1][5*i  ] = drvshp[i].y*uy;
    Bl[1][5*i+1] = drvshp[i].y*vy;
    Bl[1][5*i+2] = drvshp[i].y*wy;
    Bl[1][5*i+3] = dHdy*phi12;
    Bl[1][5*i+4] = dHdy*phi22;

    Bl[2][5*i  ] = drvshp[i].z*uz;
    Bl[2][5*i+1] = drvshp[i].z*vz;
    Bl[2][5*i+2] = drvshp[i].z*wz;
    Bl[2][5*i+3] = dHdz*phi13;
    Bl[2][5*i+4] = dHdz*phi23;

    Bl[3][5*i  ] = drvshp[i].x*uy + drvshp[i].y*ux;
    Bl[3][5*i+1] = drvshp[i].x*vy + drvshp[i].y*vx;
    Bl[3][5*i+2] = drvshp[i].x*wy + drvshp[i].y*wx;
    Bl[3][5*i+3] = dHdx*phi12 + dHdy*phi11;
    Bl[3][5*i+4] = dHdx*phi22 + dHdy*phi21;

    Bl[4][5*i  ] = drvshp[i].x*uz + drvshp[i].z*ux;
    Bl[4][5*i+1] = drvshp[i].x*vz + drvshp[i].z*vx;
    Bl[4][5*i+2] = drvshp[i].x*wz + drvshp[i].z*wx;
    Bl[4][5*i+3] = dHdx*phi13 + dHdz*phi11;
    Bl[4][5*i+4] = dHdx*phi23 + dHdz*phi21;

    Bl[5][5*i  ] = drvshp[i].y*uz + drvshp[i].z*uy;
    Bl[5][5*i+1] = drvshp[i].y*vz + drvshp[i].z*vy;
    Bl[5][5*i+2] = drvshp[i].y*wz + drvshp[i].z*wy;
    Bl[5][5*i+3] = dHdy*phi13 + dHdz*phi12;
    Bl[5][5*i+4] = dHdy*phi23 + dHdz*phi22;
  }
}

// =============================== BnlMatrix ===============================

void cShell :: BnlMatrix(int nn, double *shpfunc, double *mapfunc,
                         sNodeCoord *coord, sNodeCoord *drvshp, cMatrix &Bnl)
{
  // Get the through-thickness coordinate and its derivatives.

  double zeta = shpfunc[nn];
  double dzetadx = drvshp[nn].x;
  double dzetady = drvshp[nn].y;
  double dzetadz = drvshp[nn].z;

  // Assembly [Bnl].

  for (int i = 0; i < nn; i++)
  {
    double thk2 = *(coord[i].thk)/2.0;
    double dHdx = drvshp[i].x*zeta + shpfunc[i]*dzetadx;
    double dHdy = drvshp[i].y*zeta + shpfunc[i]*dzetady;
    double dHdz = drvshp[i].z*zeta + shpfunc[i]*dzetadz;

    /*
    double pax = -coord[i].axes->V2[0]*thk2;
    double pay = -coord[i].axes->V2[1]*thk2;
    double paz = -coord[i].axes->V2[2]*thk2;
    double pbx =  coord[i].axes->V1[0]*thk2;
    double pby =  coord[i].axes->V1[1]*thk2;
    double pbz =  coord[i].axes->V1[2]*thk2;
    */

    double pax = (*coord[i].pa)[0];
    double pay = (*coord[i].pa)[1];
    double paz = (*coord[i].pa)[2];
    double pbx = (*coord[i].pb)[0];
    double pby = (*coord[i].pb)[1];
    double pbz = (*coord[i].pb)[2];

    Bnl[0][5*i  ] = drvshp[i].x;
    Bnl[0][5*i+1] = 0.0;
    Bnl[0][5*i+2] = 0.0;
    Bnl[0][5*i+3] = dHdx*pax;
    Bnl[0][5*i+4] = dHdx*pbx;

    Bnl[1][5*i  ] = drvshp[i].y;
    Bnl[1][5*i+1] = 0.0;
    Bnl[1][5*i+2] = 0.0;
    Bnl[1][5*i+3] = dHdy*pax;
    Bnl[1][5*i+4] = dHdy*pbx;

    Bnl[2][5*i  ] = drvshp[i].z;
    Bnl[2][5*i+1] = 0.0;
    Bnl[2][5*i+2] = 0.0;
    Bnl[2][5*i+3] = dHdz*pax;
    Bnl[2][5*i+4] = dHdz*pbx;

    Bnl[3][5*i  ] = 0.0;
    Bnl[3][5*i+1] = drvshp[i].x;
    Bnl[3][5*i+2] = 0.0;
    Bnl[3][5*i+3] = dHdx*pay;
    Bnl[3][5*i+4] = dHdx*pby;

    Bnl[4][5*i  ] = 0.0;
    Bnl[4][5*i+1] = drvshp[i].y;
    Bnl[4][5*i+2] = 0.0;
    Bnl[4][5*i+3] = dHdy*pay;
    Bnl[4][5*i+4] = dHdy*pby;

    Bnl[5][5*i  ] = 0.0;
    Bnl[5][5*i+1] = drvshp[i].z;
    Bnl[5][5*i+2] = 0.0;
    Bnl[5][5*i+3] = dHdz*pay;
    Bnl[5][5*i+4] = dHdz*pby;

    Bnl[6][5*i  ] = 0.0;
    Bnl[6][5*i+1] = 0.0;
    Bnl[6][5*i+2] = drvshp[i].x;
    Bnl[6][5*i+3] = dHdx*paz;
    Bnl[6][5*i+4] = dHdx*pbz;

    Bnl[7][5*i  ] = 0.0;
    Bnl[7][5*i+1] = 0.0;
    Bnl[7][5*i+2] = drvshp[i].y;
    Bnl[7][5*i+3] = dHdy*paz;
    Bnl[7][5*i+4] = dHdy*pbz;

    Bnl[8][5*i  ] = 0.0;
    Bnl[8][5*i+1] = 0.0;
    Bnl[8][5*i+2] = drvshp[i].z;
    Bnl[8][5*i+3] = dHdz*paz;
    Bnl[8][5*i+4] = dHdz*pbz;
  }
}

// ================================= SMatrix ===============================

void cShell :: SMatrix(cVector &str, cMatrix &S)
{
  S.Zero( );
  S[0][0] = S[3][3] = S[6][6] = str[0];  // Sxx
  S[1][0] = S[4][3] = S[7][6] = str[3];  // Sxy
  S[2][0] = S[5][3] = S[8][6] = str[4];  // Sxz
  S[0][1] = S[3][4] = S[6][7] = str[3];  // Sxy
  S[1][1] = S[4][4] = S[7][7] = str[1];  // Syy
  S[2][1] = S[5][4] = S[8][7] = str[5];  // Syz
  S[0][2] = S[3][5] = S[6][8] = str[4];  // Sxz
  S[1][2] = S[4][5] = S[7][8] = str[5];  // Syz
  S[2][2] = S[5][5] = S[8][8] = str[2];  // Szz
}

// ================================= NMatrix ===============================

#if 0
void cShell :: NMatrix(int nn, double *shpfunc, sNodeCoord *coord, cMatrix &N)
{
  double zeta = shpfunc[nn];
  for (int i = 0; i < nn; i++)
  {
    double T = shpfunc[i]*zeta;
    double thk2 = *(coord[i].thk)/2.0;
    double pax = -coord[i].axes->V2[0]*thk2;
    double pay = -coord[i].axes->V2[1]*thk2;
    double paz = -coord[i].axes->V2[2]*thk2;
    double pbx =  coord[i].axes->V1[0]*thk2;
    double pby =  coord[i].axes->V1[1]*thk2;
    double pbz =  coord[i].axes->V1[2]*thk2;

    N[0][5*i  ] = shpfunc[i];
    N[0][5*i+1] = 0.0;
    N[0][5*i+2] = 0.0;
    N[0][5*i+3] = pax*T;
    N[0][5*i+4] = pbx*T;

    N[1][5*i  ] = 0.0;
    N[1][5*i+1] = shpfunc[i];
    N[1][5*i+2] = 0.0;
    N[1][5*i+3] = pay*T;
    N[1][5*i+4] = pby*T;

    N[2][5*i  ] = 0.0;
    N[2][5*i+1] = 0.0;
    N[2][5*i+2] = shpfunc[i];
    N[2][5*i+3] = paz*T;
    N[2][5*i+4] = pbz*T;
}
#endif

// ================================= NMatrix ===============================

void cShell :: NMatrix(int nn, double *shpfunc, cMatrix &N)
{
  for (int i = 0; i < nn; i++)
  {
    N[0][3*i  ] = shpfunc[i];
    N[0][3*i+1] = 0.0;
    N[0][3*i+2] = 0.0;
    N[0][3*i+3] = 0.0;
    N[0][3*i+4] = 0.0;

    N[1][3*i  ] = 0.0;
    N[1][3*i+1] = shpfunc[i];
    N[1][3*i+2] = 0.0;
    N[1][3*i+3] = 0.0;
    N[1][3*i+4] = 0.0;

    N[2][3*i  ] = 0.0;
    N[2][3*i+1] = 0.0;
    N[2][3*i+2] = shpfunc[i];
    N[2][3*i+3] = 0.0;
    N[2][3*i+4] = 0.0;

    N[3][3*i  ] = 0.0;
    N[3][3*i+1] = 0.0;
    N[3][3*i+2] = 0.0;
    N[3][3*i+3] = shpfunc[i];
    N[3][3*i+4] = 0.0;

    N[4][3*i  ] = 0.0;
    N[4][3*i+1] = 0.0;
    N[4][3*i+2] = 0.0;
    N[4][3*i+3] = 0.0;
    N[4][3*i+4] = shpfunc[i];
  }
}

// ============================== ExpandStress =============================

void cShell :: ExpandStress(cVector &sig, cVector &expsig)
{
  expsig[0] = sig[0];  // Sxx
  expsig[1] = sig[1];  // Syy
  expsig[2] = sig[2];  // Szz
  expsig[3] = sig[3];  // Sxy
  expsig[4] = sig[4];  // Sxz
  expsig[5] = sig[5];  // Syz
}

// -------------------------------------------------------------------------
// Class cInterface_2D:
// -------------------------------------------------------------------------

// ============================ cInterface2D ===============================

cInterface2D :: cInterface2D(void)
{
  Type = INTERFACE_2D;
}

// ============================ ~cInterface2D ==============================

cInterface2D :: ~cInterface2D(void)
{
}

// ============================== GetActDir ================================

void cInterface2D :: GetActDir(int *dir)
{
  dir[0] = 1;   // u
  dir[1] = 1;   // v
}

// ============================== GetStrLabels =============================

void cInterface2D :: GetStrLabels(int *label)
{
  label[0] = STRESS_XX;
  label[1] = STRESS_YY;
}

// ================================= BMatrix ===============================

void cInterface2D :: BMatrix(int nn, double *shpfunc, double *,
                             sNodeCoord *, sNodeCoord *, cMatrix &B)
{
  B.Zero( );

  int n = nn/2;  // Number of nodes in each surface (top/bottom)
  int m = 2*n;   // Number of degrees of each surface (top/bottom)
  for (int j = 0; j < n; j++)
  {
    int col = 2*j;  // u
    B[0][col    ] =  shpfunc[j];  // Top
    B[0][col + m] = -shpfunc[j];  // Bottom

    col++;          // v
    B[1][col    ] =  shpfunc[j];  // Top
    B[1][col + m] = -shpfunc[j];  // Bottom
  }
}


// -------------------------------------------------------------------------
// Class cInterface_3D:
// -------------------------------------------------------------------------

// ============================ cInterface3D ===============================

cInterface3D :: cInterface3D(void)
{
  Type = INTERFACE_3D;
}

// ============================ ~cInterface3D ==============================

cInterface3D :: ~cInterface3D(void)
{
}

// ============================== GetActDir ================================

void cInterface3D :: GetActDir(int *dir)
{
  dir[0] = 1;   // u
  dir[1] = 1;   // v
  dir[2] = 1;   // w
}

// ============================== GetStrLabels =============================

void cInterface3D :: GetStrLabels(int *label)
{
  label[0] = STRESS_XX;
  label[1] = STRESS_YY;
  label[2] = STRESS_ZZ;
}

// ================================= BMatrix ===============================

void cInterface3D :: BMatrix(int nn, double *shpfunc, double *,
                             sNodeCoord *, sNodeCoord *, cMatrix &B)
{
  B.Zero( );

  int n = nn/2;  // Number of nodes in each surface (top/bottom)
  int m = 3*n;   // Number of degrees of each surface (top/bottom)
  for (int j = 0; j < n; j++)
  {
    int col = 3*j;  // u
    B[0][col    ] =  shpfunc[j];  // Top
    B[0][col + m] = -shpfunc[j];  // Bottom

    col++;          // v
    B[1][col    ] =  shpfunc[j];  // Top
    B[1][col + m] = -shpfunc[j];  // Bottom

    col++;          // w
    B[2][col    ] =  shpfunc[j];  // Top
    B[2][col + m] = -shpfunc[j];  // Bottom
  }
}

// -------------------------------------------------------------------------
// Class cPlaneHeatTransfer:
// -------------------------------------------------------------------------

// =========================== cPlaneHeatTransfer ==========================

cPlaneHeatTransfer :: cPlaneHeatTransfer(void)
{
  Type = PLANE_HEAT_TRANSFER;
}

// =========================== ~cPlaneHeatTransfer =========================

cPlaneHeatTransfer :: ~cPlaneHeatTransfer(void)
{
}

// ============================== GetActDir ================================

void cPlaneHeatTransfer :: GetActDir(int *dir)
{
  dir[7] = 1;   // t
}

// ============================== GetStrLabels =============================

void cPlaneHeatTransfer :: GetStrLabels(int *label)
{
  label[0] = HEAT_FLUX_X;
  label[1] = HEAT_FLUX_Y;
}

// ================================= CondMatrix ============================

void cPlaneHeatTransfer :: CondMatrix(double *param, cMatrix &C)
{
  double k  = param[0];

  C.Zero( );
  C[0][0] = k;
  C[1][1] = k;
 // C[2][2] = k;
}

// ================================= CondMatrix ============================

void cPlaneHeatTransfer :: CondMatrixOrtho(double *param, cMatrix &C)
{
  double k1   = param[0];
  double k2   = param[1];

  C.Zero( );
  C[0][0] = k1;
  C[1][1] = k2;
}

// ================================= BMatrix ===============================

void cPlaneHeatTransfer :: BMatrix(int nn, double *shpfunc, double *mapfunc,
                               sNodeCoord *coord, sNodeCoord *drvshp, cMatrix &B)
{
  for (int i = 0; i < nn; i++)
  {
    B[0][i] = drvshp[i].x;
    B[1][i] = drvshp[i].y;
  }
}

// ======================================================= End of file =====
