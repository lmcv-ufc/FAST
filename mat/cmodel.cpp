// -------------------------------------------------------------------------
// cmodel.cpp - implementation of the Constitutive Model class.
// -------------------------------------------------------------------------
// Created:      15-Jul-2005     Evandro Parente Junior
//
// Modified:     18-Jul-2013     Iuri Barcelos Rocha
//               Creation of Composite Progressive Failure models.
//
// Modified:     22-Oct-2013     Evandro Parente Junior
//               Creation of ModOrthoLinear.
//
// Modified:     28-Oct-2013     Iuri Barcelos Rocha
//               Bug fixes in the progressive failure models.
//
// Modified:     09-Mai-2014     Edson Dantas/Evandro Parente Junior
//               Creation of CZM model.
//
// Modified:     02-Aug-2014     Evandro Parente Junior
//               Redefinition of Create and implementation of cConstModel1D.
//
// Modified:     02-Mar-2020     Bergson Matias/Evandro Parente Junior
//               Creation of Mazars, MuModel and Lee-Fenves materials.
//
// Modified:     01-Feb-2023     Renan Melo Barros/Evandro Parente Junior
//               Creation of elastoplastic model with isotropic hardening.
// -------------------------------------------------------------------------

#include <math.h>
#include <iostream>
#include <iomanip>

using namespace std;

#include <stdio.h>
#include <stdlib.h>

#include "cmodel.h"
#include "cmodel1d.h"
#include "element.h"
#include "elmparam.h"
#include "anmodel.h"
#include "material.h"
#include "ctrl.h"
#include "vec.h"
#include "mat.h"
#include "gblvar.h"

const int MAX_FAILURE_ITER = 10;

// -------------------------------------------------------------------------
// Public methods:
//

// ================================= Create ================================

cConstModel *cConstModel :: Create(cElement *elm, cMaterial *mat, double *p)
{
  cConstModel *cmod = 0;

  switch(mat->GetMecProp( )->GetType( ))
  {
    case MAT_ELASTIC_ISOTROPIC:
     if (elm->GetType( ) == PARAMETRIC || elm->GetType( ) == PARAMETRICTL)
       cmod = new cModLinear(elm, mat);
     else
       cmod = new cMod1DLinear(mat);
    break;

    case MAT_ELASTIC_ORTHOTROPIC:
     cmod = new cModOrthoLinear(elm, mat);
    break;

    case MAT_FGM:
     cmod = new cModFGM(elm, mat, p);
    break;

    case MAT_FGM_TTO:
     cmod = new cModFGMTTO(elm, mat, p);
    break;

    case MAT_VISCOELASTIC:
     if (elm->GetType( ) == PARAMETRIC || elm->GetType( ) == PARAMETRICTL)
       cmod = new cModViscoElastic(elm, mat);
     else
       cmod = new cMod1DViscoElastic(mat);
    break;

    case MAT_PROGFAIL_HASHIN:
     cmod = new cModProgFailHashin(elm, mat);
    break;

    case MAT_PROGFAIL_TSAI:
     cmod = new cModProgFailTsai(elm, mat);
    break;

    case MAT_PROGFAIL_ENGELSTAD:
     cmod = new cModProgFailEngelstad(elm, mat);
    break;

    case MAT_CZM_MODEI_BILINEAR:
     cmod = new cModCZMModeIBilinear(elm, mat);
    break;

    case MAT_CZM_MODEI_EXPONENTIAL:
     cmod = new cModCZMModeIExponential(elm, mat);
    break;

    case MAT_PLASTIC_ISOLINHARD:
     cmod = new cModPlastIsoLinHard(elm, mat);
    break;

    case MAT_ELASTIC_NONLINEAR:
     cmod = new cMod1DElastNonlin(mat);
    break;

    case MAT_CONCRETE_EUROCEB:
     cmod = new cMod1DConcEuroCEB(mat);
    break;

    case MAT_CONCRETE_NBRCEB:
     cmod = new cMod1DConcNBRCEB(mat);
    break;    
    
    case MAT_PLASTIC_LINHARD:
     cmod = new cMod1DPlastLinHard(mat);
    break;

    case MAT_PLASTIC_ISONONLINHARD:
     cmod = new cMod1DPlastIsoNonlinHard(mat);
    break;

    case MAT_PLASTIC_ISOEXPHARD:
     cmod = new cMod1DPlastIsoExpHard(mat);
    break;

    case MAT_DAMAGE_MAZARS:
     cmod = new cMod1DDamageMazars(mat);
    break;

    case MAT_DAMAGE_MU_MODEL:
     cmod = new cMod1DDamageMuModel(mat);
    break;

    case MAT_LEE_FENVES:
     cmod = new cMod1DLeeFenves(mat);
    break;

    default:
    break;
  }

  return(cmod);
}

// ============================= cConstModel ===============================

cConstModel :: cConstModel(cMaterial *mat)
{
  Mat = mat->GetMecProp( );
}

// ============================= cConstModel ===============================

cConstModel :: cConstModel(cElement *elm, cMaterial *mat)
{
  Mat = mat->GetMecProp( );
}

// ============================= cConstModel ===============================

cConstModel :: cConstModel(cElement *elm, cMaterial *mat, double *p)
{
  Mat = mat->GetMecProp( );
}

// ============================= ~cConstModel ==============================

cConstModel :: ~cConstModel(void)
{
}


// -------------------------------------------------------------------------
// Class cModLinear:
// -------------------------------------------------------------------------

// =============================== cModLinear ==============================

cModLinear :: cModLinear(cElement *elm, cMaterial *mat) :
              cConstModel(elm, mat)
{
  Anm = elm->GetAnModel( );
}

// ============================== ~cModLinear ==============================

cModLinear :: ~cModLinear(void)
{
}

// ================================= Stress ================================

int cModLinear :: Stress(cVector &eps, cVector &sig)
{
  // Get the elastic matrix [C].

  int cdim = sig.Dim( );
  cMatrix C(cdim, cdim);
  TangMat(C);

  // Evaluate {sig} = [C]{eps}.

  sig = C*eps;

  return(1);
}

// ================================= Stress ================================

int cModLinear :: Stress(double Tref, cVector &Temp, cVector &eps, cVector &sig)
{
//  cout << "\ncModLinear :: Stress\n";
//  cout << "{eps}\n";
//  eps.Print( );

  // Get the elastic matrix [C].

  int cdim = sig.Dim( );
  cMatrix C(cdim, cdim);
  TangMat(C);
//  cout << "[C]\n";
//  C.Print( );

  // Get the vector of thermal expansion coefficients. 

  cVector Alpha(cdim);
  AlphaVec(Alpha);
//  cout << "{alpha}\n";
//  Alpha.Print( );

  // Evaluate the temperature variation at midplane.

  double Tm = (Temp[0] + Temp[1])/2.0; 
  double dT = Tm - Tref;
//  cout << "Tm = " << Tm << "  Tref = " << Tref << "  dT = " << dT << "\n";;

  // Evaluate the thermal strains.

  cVector eth(cdim);
  eth = dT*Alpha;
//  cout << "{eth}\n";
//  eth.Print( );

  // Evaluate the effective mechanical strains {epsmec} = {eps} - {eth}.

  cVector epsmec(cdim);
  epsmec = eps - eth;
//  cout << "{epsmec}\n";
//  epsmec.Print( );

  // Evaluate {sig} = [C]({eps} - {eth}).

  sig = C*epsmec;
//  cout << "{sig}\n";
//  sig.Print( );

  return(1);
}


// ================================= TangMat ===============================

void cModLinear :: TangMat(cMatrix &C)
{
  // Get the material parameters.

  int np = Mat->NumParam( );
  double *param = new double[np];
  Mat->GetParam(param);

  // Compute the [C] matrix using the Analysis Model class.

  Anm->GetMecMod( )->QMatrix(param, C);

  // Release memory.

  delete []param;
}

// ================================= CondMatrix ============================

void cModLinear :: CondMatrix(cMatrix &C)
{
  // Get the material parameters.

  int np = Mat->NumParam( );
  double *param = new double[np];
  Mat->GetParam(param);

  // Compute the [C] matrix using the Analysis Model class.

  Anm->GetHeatMod( )->CondMatrix(param, C);

  // Release memory.

  delete []param;
}

// ================================ AlphaVec ===============================

void cModLinear :: AlphaVec(cVector &Alpha)
{
  // Get the material parameters.

  int np = Mat->NumParam( );
  double *param = new double[np];
  Mat->GetParam(param);

  // Compute the Alpha vector using the Analysis Model class.

  Anm->GetMecMod( )->AlphaVec(param, Alpha);

  // Release memory.

  delete []param;
}

// -------------------------------------------------------------------------
// Class cModOrthoLinear:
// -------------------------------------------------------------------------

// ============================= cModOrthoLinear ===========================

cModOrthoLinear :: cModOrthoLinear(cElement *elm, cMaterial *mat) :
                   cModLinear(elm, mat)
{
}

// ============================ ~cModOrthoLinear ===========================

cModOrthoLinear :: ~cModOrthoLinear(void)
{
}

// ================================= TangMat ===============================

void cModOrthoLinear :: TangMat(cMatrix &C)
{
  // Get the material parameters.

  int np = Mat->NumParam( );
  double *param = new double[np];
  Mat->GetParam(param);

  // Compute the [C] matrix using the Analysis Model class.

  Anm->GetMecMod( )->QMatrixOrtho(param, C);

  // Release memory.

  delete []param;
}

// ================================= TangMat ===============================

void cModOrthoLinear :: AlphaVec(cVector &Alpha)
{
  // Get the material parameters.

  int np = Mat->NumParam( );
  double *param = new double[np];
  Mat->GetParam(param);

  // Compute the [C] matrix using the Analysis Model class.

  Anm->GetMecMod( )->AlphaVecOrtho(param, Alpha);

  // Release memory.

  delete []param;
}


// -------------------------------------------------------------------------
// Class cModFGM:
// -------------------------------------------------------------------------

// ================================= cModFGM ===============================

cModFGM :: cModFGM(cElement *elm, cMaterial *mat, double *p) :
           cConstModel(elm, mat, p)
{
//  cout << "cModFGM :: cModFGM\n";
  Anm = elm->GetAnModel( );
  V2 = p[0];
//  cout << "V2 = " << V2 << "\n";
//  cout << "End of cModFGM :: cModFGM\n";
}

// ================================ ~cModFGM ===============================

cModFGM :: ~cModFGM(void)
{
//  cout << "cModFGM :: ~cModFGM\n";
}

// ================================= Stress ================================

int cModFGM :: Stress(cVector &eps, cVector &sig)
{
//  cout << "cModFGM :: Stress\n";

  // Get the elastic matrix [C].

  int cdim = sig.Dim( );
  cMatrix C(cdim, cdim);
  TangMat(C);

  // Evaluate {sig} = [C]{eps}.

  sig = C*eps;
//  cout << "End of cModFGM :: Stress\n";
  return(1);
}

// ================================= TangMat ===============================

void cModFGM :: TangMat(cMatrix &C)
{
  //cout << "cModFGM :: TangMat\n";
  // Get the material parameters.

  int np = Mat->NumParam( );
  double *param = new double[np];
  Mat->GetParam(param);
  double E1   = param[0];
  double Nu1  = param[1];
  double E2   = param[2];
  double Nu2  = param[3];
  int method  = int(param[4]);

  // Homogenization.

  double E, Nu;
  if (method == 1) // Rule of Mixtures.
  {
//    cout << "Rule of Mixtures\n";
    E  = (E2 - E1)*V2 + E1;
    Nu = (Nu2 - Nu1)*V2 + Nu1;
  }
  else            // Mori-Tanaka
  {
//    cout << "Mori-Tanaka\n";
    // Bulk and shear moduli of matrix and inclusions.

    double K1 = E1/(3.0 - 6.0*Nu1);
    double K2 = E2/(3.0 - 6.0*Nu2);
    double G1 = E1/(2.0 + 2.0*Nu1);
    double G2 = E2/(2.0 + 2.0*Nu2);

    // Auxiliary factors.

    double V1 = 1.0 - V2;
    double f1 = G1*(9.0*K1 + 8.0*G1)/(6.0*K1 + 12.0*G1);
    double a  = (K2 - K1)/(K1 + 4.0*G1/3.0);
    double b  = (G2 - G1)/(G1 + f1);

    // Mori-Tanaka homogenized bulk and shear moduli.

    double K = K1 + (K2 - K1)*V2/(1.0 + V1*a);
    double G = G1 + (G2 - G1)*V2/(1.0 + V1*b);

    // Transform back into Elastic Modulus and Poisson's ratio.

    E  = 9.0*K*G/(3.0*K + G);
    Nu = (3.0*K - 2.0*G)/(6.0*K + 2.0*G);
  }
//  cout << "V2 = " << fixed << setprecision(5) << V2 << "  ";
//  cout << "E  = " << scientific << setprecision(5) << E  << "  ";
//  cout << "Nu = " << fixed << setprecision(5) << Nu << "\n";

  // Compute the [C] matrix using the Analysis Model class.

  double fgmparam[2];
  fgmparam[0] = E;  // Combined elasticity modulus
  fgmparam[1] = Nu; // Combined Poisson's ratio
  Anm->GetMecMod( )->QMatrix(fgmparam, C);
//  cout << "[C]\n";
//  C.Print( );

  // Release memory.

  delete []param;
//  cout << "End of cModFGM :: TangMat \n";
}

// ================================= AlphaVec ==============================

void cModFGM :: AlphaVec(cVector &Alpha)
{
//  cout << "cModFGM :: AlphaVec \n";
  // Get the material parameters.

  int np = Mat->NumParam( );
  double *param = new double[np];
  Mat->GetParam(param);
  double E1   = param[0];
  double Nu1  = param[1];
  double E2   = param[2];
  double Nu2  = param[3];
  int method  = int(param[4]);
  double alp1 = param[5];
  double alp2 = param[6];

  // Homogenization.

  double E, Nu, Alp;
  if (method == 1) // Rule of Mixtures.
  {
//    cout << "Rule of Mixtures\n";
    E  = (E2 - E1)*V2 + E1;
    Nu = (Nu2 - Nu1)*V2 + Nu1;
    Alp = (alp2 - alp1)*V2 + alp1;
  }
  else            // Mori-Tanaka
  {
//    cout << "Mori-Tanaka\n";
    // Bulk and shear moduli of matrix and inclusions.

    double K1 = E1/(3.0 - 6.0*Nu1);
    double K2 = E2/(3.0 - 6.0*Nu2);
    double G1 = E1/(2.0 + 2.0*Nu1);
    double G2 = E2/(2.0 + 2.0*Nu2);

    // Auxiliary factors.

    double V1 = 1.0 - V2;
    double f1 = G1*(9.0*K1 + 8.0*G1)/(6.0*K1 + 12.0*G1);
    double a  = (K2 - K1)/(K1 + 4.0*G1/3.0);
    double b  = (G2 - G1)/(G1 + f1);

    // Mori-Tanaka homogenized bulk and shear moduli.

    double K = K1 + (K2 - K1)*V2/(1.0 + V1*a);
    double G = G1 + (G2 - G1)*V2/(1.0 + V1*b);

    // Transform back into Elastic Modulus and Poisson's ratio.

    E  = 9.0*K*G/(3.0*K + G);
    Nu = (3.0*K - 2.0*G)/(6.0*K + 2.0*G);
    Alp = alp1 + (alp2 - alp1)*(1.0/K - 1.0/K1)/(1.0/K2 - 1.0/K1);
  }

  // Compute the {alpha} vector using the Analysis Model class.

  double fgmparam[3];
  fgmparam[0] = E;  // Combined elasticity modulus
  fgmparam[1] = Nu; // Combined Poisson's ratio
  fgmparam[2] = Alp;  // Combined Thermal Expasion
  Anm->GetMecMod( )->AlphaVec(fgmparam, Alpha);
//  cout << "{alpha}\n";
//  Alpha.Print( );

  // Release memory.

  delete []param;
//  cout << "End of cModFGM :: AlphaVec \n";
}


// -------------------------------------------------------------------------
// Class cModViscoElastic:
// -------------------------------------------------------------------------

// ============================= cModViscoElastic ==========================

cModViscoElastic :: cModViscoElastic(cElement *elm, cMaterial *mat) :
                    cConstModel(elm, mat)
{
  Anm  = elm->GetAnModel( );
  Time = 0.0;

  int nsc = Anm->GetMecMod( )->GetDimQMatrix( );
  Sig.Resize(nsc);
  Sig.Zero( );
  SigIter.Resize(nsc);
  SigIter.Zero( );
  Eps.Resize(nsc);
  Eps.Zero( );
  EpsIter.Resize(nsc);
  EpsIter.Zero( );
  EpsRate.Resize(nsc);
  EpsRate.Zero( );

  int np = Mat->NumParam( );
  int nprony = (np - 2)/2;
  S = new cVector[nprony];
  for (int i = 0; i < nprony; i++)
  {
    S[i].Resize(nsc);
    S[i].Zero( );
  }
}

// ============================ ~cModViscoElastic ==========================

cModViscoElastic :: ~cModViscoElastic(void)
{
  delete []S;
}

// ================================= Stress ================================

int cModViscoElastic :: Stress(cVector &eps, cVector &sig)
{
  // Start with the stored stresses.

  sig = Sig;
//  double t  = cControl::GetTotTime( );
  double dt = cControl::GetTimeStep( );

  // Get material parameters.

  int np = Mat->NumParam( );
  double *param = new double[np];
  Mat->GetParam(param);
  int nprony = (np - 2)/2;

  // Compute {dsigbar} = [Cbar]*{deps}.

  int i;
  int nsc = Anm->GetMecMod( )->GetDimQMatrix( );
  cMatrix C(nsc, nsc);
  cVector dsigbar(nsc);
  cVector deps(nsc);
  deps = eps - Eps;
  double E,rho;
  double Et = param[0];
  for (i = 0; i < nprony; i++)
  {
    E   = param[2*i+2];
    rho = param[2*i+3];
    Et += E*(rho/dt)*(1.0 - exp(-dt/rho));
  }
  param[0] = Et;
  Anm->GetMecMod( )->QMatrix(param, C);
  dsigbar = C*deps;
  sig += dsigbar;

  // Compute {dsighat} => hereditary terms.

  cVector s(nsc);
  cVector aux(nsc);
  cVector dsighat(nsc);
  dsighat.Zero( );
  param[0] = 1.0;
  Anm->GetMecMod( )->QMatrix(param, C);
  for (i = 0; i < nprony; i++)
  {
    E   = param[2*i+2];
    rho = param[2*i+3];
    aux = C*EpsRate;
    s  = exp(-dt/rho)*S[i];
    s += E*rho*(1.0 - exp(-dt/rho))*aux;
    dsighat += -(1.0 - exp(-dt/rho))*s;
  }
  sig += dsighat;

  // Store the iterative vectors.

  EpsIter = eps;
  SigIter = sig;

  // Release memory.

  delete []param;
  return(1);
}

// ================================= TangMat ===============================

void cModViscoElastic :: TangMat(cMatrix &C)
{
  // Get material parameters.

  int np = Mat->NumParam( );
  double *param = new double[np];
  Mat->GetParam(param);
  int nprony = (np - 2)/2;

  // Evalute the tangent relaxation modulus.

  double dt = cControl::GetTimeStep( );
  double Et = param[0];   // Infinite modulus
  double E,rho;
  for (int i = 0; i < nprony; i++)
  {
    E   = param[2*i+2];
    rho = param[2*i+3];
    Et += E*(rho/dt)*(1.0 - exp(-dt/rho));
  }
  param[0] = Et;          // Set E = Et;

  // Evaluate [C] and release the allocated memory.

  Anm->GetMecMod( )->QMatrix(param, C);
  delete []param;
}

// =============================== UpdateState =============================

void cModViscoElastic :: UpdateState(void)
{
  // Get material parameters.

  int np = Mat->NumParam( );
  double *param = new double[np];
  Mat->GetParam(param);
  int nprony = (np - 2)/2;

  // Update the hereditary vectors {S}.

  int nsc = Anm->GetMecMod( )->GetDimQMatrix( );
  cMatrix C(nsc, nsc);
  param[0] = 1.0;
  Anm->GetMecMod( )->QMatrix(param, C);
  cVector aux(nsc);
  double E,rho;
  double dt = cControl::GetTimeStep( );
  for (int i = 0; i < nprony; i++)
  {
    E   = param[2*i+2];
    rho = param[2*i+3];
    aux = C*EpsRate;
    S[i] *= exp(-dt/rho);
    S[i] += E*rho*(1.0 - exp(-dt/rho))*aux;
  }

  // Update the internal vectors.

  Time = cControl::GetTotTime( );
  aux  = EpsIter - Eps;
  Eps  = EpsIter;
  Sig  = SigIter;
  EpsRate = (1.0/dt)*aux;

  // Release memory.

  delete []param;
}


// -------------------------------------------------------------------------
// Class cModProgFailure:
// -------------------------------------------------------------------------

// ============================ cModProgFailure ============================

cModProgFailure :: cModProgFailure(cElement *elm, cMaterial *mat) :
                   cConstModel(elm, mat)
{
  Anm = elm->GetAnModel( );
}

// =========================== ~cModProgFailure ============================

cModProgFailure :: ~cModProgFailure(void)
{
}

// ================================ Stress =================================

int cModProgFailure :: Stress(cVector &eps, cVector &sig)
{
  // Create [C] matrix.

  int cdim = sig.Dim( );
  cMatrix C(cdim, cdim);

  // Initialize the iterative material parameters (ParamIter = Param).

  InitIterParam( );

  // Evaluate the iterative stresses.

  for (int i = 0; i < MAX_FAILURE_ITER; i++)
  {
    // Evaluate [C] using the current (iterative) parameters.

    Anm->GetMecMod( )->QMatrixOrtho(ParamIter, C);

    // Evaluate the corresponding stresses.

    sig = C*eps;

    // Degrade the iterative parameters when a new failure is detected.

    if (FailCrit(sig))
      Degrade( );
    else
      break;
  }

  return(1);
}

// ================================= TangMat ===============================

void cModProgFailure :: TangMat(cMatrix &C)
{
  // Compute the C Matrix using current material properties.

  Anm->GetMecMod( )->QMatrixOrtho(ParamIter, C);
}


// -------------------------------------------------------------------------
// Class cModProgFailHashin:
// -------------------------------------------------------------------------

// =============================== cModProgFail ============================

cModProgFailHashin :: cModProgFailHashin(cElement *elm, cMaterial *mat) :
                      cModProgFailure(elm, mat)
{
  // Memory allocation.

  Strength = new double[9];
  Param = new double[9];
  ParamIter = new double[9];
  FailedIter = new bool[9];

  // Get material parameters.

  int np = Mat->NumParam( );
  double *param = new double[np];
  Mat->GetParam(param);

  // Store material properties in corresponding vectors.

  for (int i = 0; i < 9; i++)
  {
    FailedIter[i] = false;
    Param[i] = ParamIter[i] = param[i];
    Strength[i] = param[i+9];
  }

  // Get the degradation factor.

  Alpha = param[18];

  delete []param;
}

// ============================== ~cModProgFail ============================

cModProgFailHashin :: ~cModProgFailHashin(void)
{
  delete []Strength;
  delete []FailedIter;
  delete []Param;
  delete []ParamIter;
}

// ================================ FailCrit ===============================

bool cModProgFailHashin :: FailCrit(cVector &sig)
{
  // Flag to determine if failure has occurred.

  bool hasFailed = false;

  for (int i = 0; i < 9; i++)
    FailedIter[i] = false;

  // Get stress values.

  cVector expsig(6);
  Anm->GetMecMod( )->ExpandStress(sig, expsig);
  double sig11 = expsig[0];
  double sig22 = expsig[1];
  double sig33 = expsig[2];
  double sig12 = expsig[3];
  double sig13 = expsig[4];
  double sig23 = expsig[5];

  // Get strength values.

  double Xt = Strength[0];
  double Xc = Strength[1];
  double Yt = Strength[2];
  double Yc = Strength[3];
  double S12 = Strength[6];
  double S23 = Strength[8];

  // Useful variables.

  double sqXt = Xt*Xt;
  double sqYt = Yt*Yt;
  double sqYc = Yc*Yc;
  double sqS12 = S12*S12;
  double sqS23 = S23*S23;

  double sqsig11 = sig11*sig11;
  double sqsig12 = sig12*sig12;
  double sqsig13 = sig13*sig13;
  double sqsig23 = sig23*sig23;

  double sqsigsum = (sig22 + sig33)*(sig22 + sig33);

  // Check for matrix tensile failure.

  double aux;
  if ((sig22 + sig33) >= 0)
  {
    aux = sqsigsum/sqYt + (sqsig12 + sqsig13)/sqS12 + (sqsig23 - sig22*sig33)/sqS23;

    if (aux >= 1.0)
    {
      hasFailed = true;

      // Degradation in E2, G12, G13, G23, nu12, nu13 and nu23.

      if (!FailedIter[1]) FailedIter[1] = true;
      if (!FailedIter[3]) FailedIter[3] = true;
      if (!FailedIter[4]) FailedIter[4] = true;
      if (!FailedIter[5]) FailedIter[5] = true;
      if (!FailedIter[6]) FailedIter[6] = true;
      if (!FailedIter[7]) FailedIter[7] = true;
      if (!FailedIter[8]) FailedIter[8] = true;
    }
  }

  // Check for matrix compressive failure.

  if ((sig22 + sig33) < 0)
  {
    aux = (sqYc/(4.0*sqS23) - 1)*(sig22 + sig33)/Yc + sqsigsum/(4.0*sqS23) + (sqsig23 - sig22*sig33)/sqS23 + (sqsig12 + sqsig13)/sqS12;

    if (aux >= 1.0)
    {
      hasFailed = true;

      // Degradation in E2, G12, G13, G23, nu12, nu13 and nu23.

      if (!FailedIter[1]) FailedIter[1] = true;
      if (!FailedIter[3]) FailedIter[3] = true;
      if (!FailedIter[4]) FailedIter[4] = true;
      if (!FailedIter[5]) FailedIter[5] = true;
      if (!FailedIter[6]) FailedIter[6] = true;
      if (!FailedIter[7]) FailedIter[7] = true;
      if (!FailedIter[8]) FailedIter[8] = true;
    }
  }

  // Check for fiber tensile failure.

  if (sig11 >= 0)
  {
    aux = sqsig11/sqXt + (sqsig12 + sqsig13)/sqS12;

    if (aux >= 1.0)
    {
      hasFailed = true;

      // Degradation in all properties.

      if (!FailedIter[0]) FailedIter[0] = true;
      if (!FailedIter[1]) FailedIter[1] = true;
      if (!FailedIter[2]) FailedIter[2] = true;
      if (!FailedIter[3]) FailedIter[3] = true;
      if (!FailedIter[4]) FailedIter[4] = true;
      if (!FailedIter[5]) FailedIter[5] = true;
      if (!FailedIter[6]) FailedIter[6] = true;
      if (!FailedIter[7]) FailedIter[7] = true;
      if (!FailedIter[8]) FailedIter[8] = true;
    }
  }

  // Check for fiber compressive failure.

  if (sig11 < 0)
  {
    aux = fabs(sig11)/Xc;

    if (aux >= 1.0)
    {
      hasFailed = true;

      // Degradation in all properties.

      if (!FailedIter[0]) FailedIter[0] = true;
      if (!FailedIter[1]) FailedIter[1] = true;
      if (!FailedIter[2]) FailedIter[2] = true;
      if (!FailedIter[3]) FailedIter[3] = true;
      if (!FailedIter[4]) FailedIter[4] = true;
      if (!FailedIter[5]) FailedIter[5] = true;
      if (!FailedIter[6]) FailedIter[6] = true;
      if (!FailedIter[7]) FailedIter[7] = true;
      if (!FailedIter[8]) FailedIter[8] = true;
    }
  }

  // Check for fiber-matrix shear-out failure.

  aux = sqsig11/sqXt + sqsig12/sqS12;

  if (aux >= 1.0)
  {
    hasFailed = true;

    // Degradation in E2, G12 and nu12.

    if (!FailedIter[1]) FailedIter[1] = true;
    if (!FailedIter[3]) FailedIter[3] = true;
    if (!FailedIter[6]) FailedIter[6] = true;
  }

  // Return 1 if failure has occurred. Otherwise, return 0.

  return hasFailed;
}

// ================================ Degrade ================================

void cModProgFailHashin :: Degrade(void)
{
  // Degrade the iterative parameters.

  for (int i = 0; i < 9; i++)
    if (FailedIter[i]) ParamIter[i] *= Alpha;
}

// ============================== InitIterParam ============================

void cModProgFailHashin :: InitIterParam(void)
{
  // Roll back material parameters.

  for (int i = 0; i < 9; i++)
    ParamIter[i] = Param[i];
}

// =============================== UpdateState =============================

void cModProgFailHashin :: UpdateState(void)
{
  // Update material parameters.

  for (int i = 0; i < 9; i++)
    Param[i] = ParamIter[i];
}


// -------------------------------------------------------------------------
// Class cModProgFailTsai:
// -------------------------------------------------------------------------

// ============================ cModProgFailTsai ===========================

cModProgFailTsai :: cModProgFailTsai(cElement *elm, cMaterial *mat) :
                    cModProgFailure(elm, mat)
{
  // Degradation parameters.

  Em =  0.08;
  Ef =  0.01;
  N  =  0.1;
  Ic = IcIter = -0.5;

  // Memory allocation.

  Strength = new double[9];
  Param = new double[9];
  ParamIter = new double[9];
  Failed = new bool[2];
  FailedIter = new bool[2];
  FailMode = new bool[2];

  // Get material parameters.

  int np = Mat->NumParam( );
  double *param = new double[np];
  Mat->GetParam(param);

  // Store material parameters in corresponding vectors.

  for (int i = 0; i < 9; i++)
  {
    Param[i] = ParamIter[i] = param[i];
    Strength[i] = param[i+9];
  }

  // Initialize failure flags.

  for (int i = 0; i < 2; i++)
    Failed[i] = FailedIter[i] = FailMode[i] = false;

  delete []param;
}

// =========================== ~cModProgFailTsai ===========================

cModProgFailTsai :: ~cModProgFailTsai(void)
{
  delete []Strength;
  delete []Param;
  delete []ParamIter;
  delete []Failed;
  delete []FailedIter;
  delete []FailMode;
}

// ================================ FailCrit ===============================

bool cModProgFailTsai :: FailCrit(cVector &sig)
{
  // Get the current (iterative) failure state.

  for (int i = 0; i < 2; i++)
    FailMode[i] = FailedIter[i];

  // Prevent the point from failing more than 2 times.

  if (Failed[1])
    return false;

  // Get stress values.

  cVector expsig(6);
  Anm->GetMecMod( )->ExpandStress(sig, expsig);
  double sig11 = expsig[0];
  double sig22 = expsig[1];
  double sig33 = expsig[2];
  double sig12 = expsig[3];
  double sig13 = expsig[4];
  double sig23 = expsig[5];

  // Auxiliary stress variables.

  double aux1 = sig11*sig11;
  double aux2 = sig22*sig22;
  double aux3 = sig33*sig33;
  double aux4 = sig23*sig23;
  double aux5 = sig13*sig13;
  double aux6 = sig12*sig12;

  double aux12 = sig11*sig22;
  double aux13 = sig11*sig33;
  double aux23 = sig22*sig33;

  // Get strength parameters.

  double Xt = Strength[0];
  double Xc = Strength[1];
  double Yt = Strength[2];
  double Yc = Strength[3];
  double Zt = Strength[4];
  double Zc = Strength[5];
  double S12 = Strength[6];
  double S13 = Strength[7];
  double S23 = Strength[8];

  // Degrade the longitudinal compressive strength, if necessary.

  if (FailedIter[0])
    Xc *= pow(Em, N);

  // Evaluate the Tsai-Wu Safety Factor (SF).

  double f1  = 1.0/Xt - 1.0/Xc;
  double f2  = 1.0/Yt - 1.0/Yc;
  double f3  = 1.0/Zt - 1.0/Zc;
  double f11 = 1.0/(Xt*Xc);
  double f22 = 1.0/(Yt*Yc);
  double f33 = 1.0/(Zt*Zc);
  double f44 = 1.0/(S23*S23);
  double f55 = 1.0/(S13*S13);
  double f66 = 1.0/(S12*S12);
  double f12 = IcIter*sqrt(f11*f22);
  double f13 = IcIter*sqrt(f11*f33);
  double f23 = IcIter*sqrt(f22*f33);

  double a = f11*aux1 + f22*aux2 + f33*aux3 + f44*aux4 + f55*aux5 + f66*aux6 + 2.0*f12*aux12 + 2.0*f13*aux13 + 2.0*f23*aux23;
  double b = (f1*sig11 + f2*sig22 + f3*sig33);

  double aux = 0.5*(-b + sqrt(b*b + 4.0*a))/a;

  // Check for failure using the Failure Index (FI).

  if (1.0/aux >= 1.0)
  {
    if (!FailMode[0])
      FailMode[0] = true; // First failure (matrix).
    else
      FailMode[1] = true; // Second failure (fiber).
    return true;
  }
  else
    return false; // No failure detected.
}

// ================================ Degrade ================================

void cModProgFailTsai :: Degrade(void)
{
  // Degrade the iterative parameters.

  if (FailMode[0] && !FailMode[1]) // First failure (matrix).
  {
    FailedIter[0] = true;
    ParamIter[1] *= Em;
    ParamIter[3] *= Em;
    ParamIter[4] *= Em;
    ParamIter[5] *= Em;
    ParamIter[6] *= Em;
    ParamIter[7] *= Em;
    ParamIter[8] *= Em;
    IcIter *= Em;
  }

  else if (FailMode[1]) // Second failure (fiber).
  {
    FailedIter[1] = true;
    for (int i = 0; i < 9; i++)
      ParamIter[i] *= Ef;
    IcIter *= Ef;
  }
}

// ============================== InitIterParam ============================

void cModProgFailTsai :: InitIterParam(void)
{
  // Roll back material parameters.

  for (int i = 0; i < 9; i++)
    ParamIter[i] = Param[i];

  // Roll back interaction coefficient.

  IcIter = Ic;

  // Roll back failure flags.

  for (int i = 0; i < 2; i++)
    FailedIter[i] = Failed[i];
}

// ============================== UpdateState =============================

void cModProgFailTsai :: UpdateState(void)
{
  // Update material parameters.

  for (int i = 0; i < 9; i++)
    Param[i] = ParamIter[i];

  // Update interaction coefficient.

  Ic = IcIter;

  // Update failure flags.

  for (int i = 0; i < 2; i++)
    Failed[i] = FailedIter[i];
}


// -------------------------------------------------------------------------
// Class cModProgFailEngelstad:
// -------------------------------------------------------------------------

// ========================= cModProgFailEngelstad =========================

cModProgFailEngelstad :: cModProgFailEngelstad(cElement *elm, cMaterial *mat) :
                         cModProgFailure(elm, mat)
{
  // Stress interaction coefficient.

  Ic = -0.0;

  // Memory allocation.

  Strength = new double[9];
  Param = new double[9];
  ParamIter = new double[9];
  Failed = new bool[9];
  FailedIter = new bool[9];

  // Get material parameters.

  int np = Mat->NumParam( );
  double *param = new double[np];
  Mat->GetParam(param);

  // Store material parameters in corresponding vectors.

  for (int i = 0; i < 9; i++)
  {
    Failed[i] = FailedIter[i] = false;
    Param[i] = ParamIter[i] = param[i];
    Strength[i] = param[i+9];
  }

  // Get the degradation factor.

  Alpha = param[18];

  delete []param;
}

// ============================== ~cModProgFail ============================

cModProgFailEngelstad :: ~cModProgFailEngelstad(void)
{
  delete []Strength;
  delete []Param;
  delete []Failed;
  delete []FailedIter;
}

// ================================ FailCrit ===============================

bool cModProgFailEngelstad :: FailCrit(cVector &sig)
{
  // Get stress values.

  cVector expsig(6);
  Anm->GetMecMod( )->ExpandStress(sig, expsig);
  double sig11 = expsig[0];
  double sig22 = expsig[1];
  double sig33 = expsig[2];
  double sig12 = expsig[3];
  double sig13 = expsig[4];
  double sig23 = expsig[5];

  // Auxiliary stress variables.

  double aux1 = sig11*sig11;
  double aux2 = sig22*sig22;
  double aux3 = sig33*sig33;
  double aux4 = sig23*sig23;
  double aux5 = sig13*sig13;
  double aux6 = sig12*sig12;

  double aux12 = sig11*sig22;
  double aux13 = sig11*sig33;
  double aux23 = sig22*sig33;

  // Get strength parameters.

  double Xt = Strength[0];
  double Xc = Strength[1];
  double Yt = Strength[2];
  double Yc = Strength[3];
  double Zt = Strength[4];
  double Zc = Strength[5];
  double S12 = Strength[6];
  double S13 = Strength[7];
  double S23 = Strength[8];

  // Evaluate the Tsai-Wu Safety Factor (SF).

  double f1  = 1.0/Xt - 1.0/Xc;
  double f2  = 1.0/Yt - 1.0/Yc;
  double f3  = 1.0/Zt - 1.0/Zc;
  double f11 = 1.0/(Xt*Xc);
  double f22 = 1.0/(Yt*Yc);
  double f33 = 1.0/(Zt*Zc);
  double f44 = 1.0/(S23*S23);
  double f55 = 1.0/(S13*S13);
  double f66 = 1.0/(S12*S12);
  double f12 = Ic*sqrt(f11*f22);
  double f13 = Ic*sqrt(f11*f33);
  double f23 = Ic*sqrt(f22*f33);

  double a = f11*aux1 + f22*aux2 + f33*aux3 + f44*aux4 + f55*aux5 + f66*aux6 + 2.0*f12*aux12 + 2.0*f13*aux13 + 2.0*f23*aux23;
  double b = (f1*sig11 + f2*sig22 + f3*sig33);

  double aux = 0.5*(-b + sqrt(b*b + 4.0*a))/a;

  // Check for failure using the Failure Index (FI).

  if (1.0/aux >= 1.0)
  {
    // Evaluate the failure mode coefficients.

    double H1 = fabs(f1*sig11 + f11*aux1);
    double H2 = fabs(f2*sig22 + f22*aux2);
    double H6 = f66*aux6;
    double H4 = f44*aux4;
    double H5 = f55*aux5;

    // Determine the failure mode.

    if (H1 > H2 && H1 > H4 && H1 > H5 && H1 > H6)
    {
      // Fiber failure (11).

      if (!FailedIter[0]) FailedIter[0] = true;
      if (!FailedIter[1]) FailedIter[1] = true;
      if (!FailedIter[2]) FailedIter[2] = true;
      if (!FailedIter[3]) FailedIter[3] = true;
      if (!FailedIter[4]) FailedIter[4] = true;
      if (!FailedIter[5]) FailedIter[5] = true;
      if (!FailedIter[6]) FailedIter[6] = true;
      if (!FailedIter[7]) FailedIter[7] = true;
      if (!FailedIter[8]) FailedIter[8] = true;
    }
    else if (H2 > H1 && H2 > H4 && H2 > H5 && H2 > H6)
    {
      // Matrix failure (22).

      if (!FailedIter[1]) FailedIter[1] = true;
      if (!FailedIter[3]) FailedIter[3] = true;
      if (!FailedIter[6]) FailedIter[6] = true;
    }
    else if (H6 > H2 && H6 > H1 && H6 > H4 && H6 > H5)
    {
      // Fiber-matrix shear-out failure (12).

      if (!FailedIter[1]) FailedIter[1] = true;
      if (!FailedIter[3]) FailedIter[3] = true;
      if (!FailedIter[6]) FailedIter[6] = true;
    }
    else if (H4 > H1 && H4 > H2 && H4 > H5 && H4 > H6)
    {
      // Transverse shear failure (13).

      if (!FailedIter[7]) FailedIter[7] = true;
      if (!FailedIter[4]) FailedIter[4] = true;
    }
    else if (H5 > H1 && H5 > H2 && H5 > H4 && H5 > H6)
    {
      // Transverse shear failure (23).

      if (!FailedIter[8]) FailedIter[8] = true;
      if (!FailedIter[5]) FailedIter[5] = true;
    }
    return true;
  }

  // If no failure was detected, return 0.

  return false;
}

// ================================ Degrade ================================

void cModProgFailEngelstad :: Degrade(void)
{
  // Degrade the iterative parameters.

  for (int i = 0; i < 9; i++)
    if (FailedIter[i] && !Failed[i]) ParamIter[i] *= Alpha;
}

// ============================== InitIterParam ============================

void cModProgFailEngelstad :: InitIterParam(void)
{
  // Roll back material parameters and failure flags.

  for (int i = 0; i < 9; i++)
  {
    ParamIter[i] = Param[i];
    FailedIter[i] = Failed[i];
  }
}

// =============================== UpdateState =============================

void cModProgFailEngelstad :: UpdateState(void)
{
  // Update material parameters and failure flags.

  for (int i = 0; i < 9; i++)
  {
    Param[i] = ParamIter[i];
    Failed[i] = FailedIter[i];
  }
}


// -------------------------------------------------------------------------
// Class cModCZMModeIBilinear:
// -------------------------------------------------------------------------

// =========================== cModCZMModeIBilinear ========================

cModCZMModeIBilinear :: cModCZMModeIBilinear(cElement *elm, cMaterial *mat) :
                        cConstModel(elm, mat)
{
  Anm = elm->GetAnModel( );

  // Get material parameters.

  int np = Mat->NumParam( );
  double *param = new double[np];
  Mat->GetParam(param);
  double Ft  = param[0];
//  double GIc = param[1];
  double K0  = param[2];
  delete []param;

  // Initialize the internal variables.

  DamVar  = DamVarIter  = 0.0;
  MaxDisp = MaxDispIter = Ft/K0;
  Status  = StatusIter  = ELASTIC;
}

// ========================== ~cModCZMModeIBilinear ========================

cModCZMModeIBilinear :: ~cModCZMModeIBilinear(void)
{
}

// ================================= TangMat ===============================

void cModCZMModeIBilinear :: TangMat(cMatrix &C)
{
  // Initialization.

  C.Zero( );

  // Get material parameters.

  int np = Mat->NumParam( );
  double *param = new double[np];
  Mat->GetParam(param);
  double Ft  = param[0];
  double GIc = param[1];
  double K0  = param[2];
  double d0  = Ft/K0;       // Initial elastic limit
  double df  = 2.0*GIc/Ft;  // Failure displacement
  delete []param;

  // Evaluate the tangent stiffness.

  double Kntan;
  if (StatusIter == COMPRESSION)
  {
    Kntan = K0;
  }
  else if (StatusIter == ELASTIC) // Elastic loading/unloading
  {
    Kntan = (1.0 - DamVarIter)*K0;
  }
  else if (StatusIter == DAMAGE)  // Damage evolution
  {
    Kntan = -Ft/(df - d0);
  }
  else // Total failure
  {
    Kntan = 0.0;
  }

  // Handle 2D/3D models.

  if (Anm->GetMecMod( )->GetDimQMatrix( ) == 2)
    C[1][1] = Kntan;
  else
    C[2][2] = Kntan;
}


// ================================= Stress ================================

int cModCZMModeIBilinear :: Stress(cVector &eps, cVector &sig)
{
  // Initialization.

  sig.Zero( );

  // Get material parameters.

  int np = Mat->NumParam( );
  double *param = new double[np];
  Mat->GetParam(param);
  double Ft  = param[0];
  double GIc = param[1];
  double K0  = param[2];
  double d0  = Ft/K0;       // Initial elastic limit
  double df  = 2.0*GIc/Ft;  // Failure displacement
  delete []param;

  // Initialize the iterative material parameters (ParamIter = Param).

  DamVarIter  = DamVar;
  MaxDispIter = MaxDisp;
  StatusIter  = Status;

  // Get the opening displacement.

  double dn;
  if (Anm->GetMecMod( )->GetDimQMatrix( ) == 2)
    dn = eps[1];
  else
    dn = eps[2];

  // Stress computation.

  double Knsec;
  if (dn < 0.0) // Penalize negative displacements (penetration).
  {
    StatusIter = COMPRESSION;
    Knsec = K0;
  }
  else          // Tension
  {
    if (dn <= MaxDispIter)  // Elastic loading/unloading
    {
      StatusIter = ELASTIC;
      Knsec = (1.0 - DamVarIter)*K0;
    }
    else if (dn <= df)      // Damage evolution
    {
      StatusIter  = DAMAGE;
      MaxDispIter = dn;
      DamVarIter  = (df/MaxDispIter)*(MaxDispIter - d0)/(df - d0);
      Knsec = (1.0 - DamVarIter)*K0;
    }
    else                    // Total failure
    {
      StatusIter  = FAILURE;
      MaxDispIter = dn;
      DamVarIter  = 1.0;
      Knsec = 0.0;
    }
  }

  // Handle 2D/3D models.

  if (Anm->GetMecMod( )->GetDimQMatrix( ) == 2)
    sig[1] = Knsec*dn;
  else
    sig[2] = Knsec*dn;

  return(1);
}

// =============================== UpdateState =============================

void cModCZMModeIBilinear :: UpdateState(void)
{
  // Update material parameters.

  Status  = StatusIter;
  DamVar  = DamVarIter;
  MaxDisp = MaxDispIter;
}


// -------------------------------------------------------------------------
// Class cModCZMModeIExponential:
// -------------------------------------------------------------------------

// =========================== cModCZMModeIExponential ========================

cModCZMModeIExponential :: cModCZMModeIExponential(cElement *elm, cMaterial *mat) :
                           cConstModel(elm, mat)
{
  Anm = elm->GetAnModel( );

  // Initialize the internal variables.

  DamVar  = DamVarIter  = 0.0;
  MaxDisp = MaxDispIter = 0.0;
  Status  = StatusIter  = DAMAGE;
}

// ========================== ~cModCZMModeIBilinear ========================

cModCZMModeIExponential :: ~cModCZMModeIExponential(void)
{
}

// ================================= TangMat ===============================

void cModCZMModeIExponential :: TangMat(cMatrix &C)
{
  // Initialization.

  C.Zero( );

  // Get material parameters.

  int np = Mat->NumParam( );
  double *param = new double[np];
  Mat->GetParam(param);
  double Ft  = param[0];
  double GIc = param[1];
  double df  = GIc/(exp(1.0)*Ft);  // Separation corresponding to Ft
  double K0  = exp(1.0)*Ft/df;     // Initial stiffness
  delete []param;

  // Evaluate the tangent stiffness.

  double Kntan;
  if (StatusIter == COMPRESSION)
  {
    Kntan = K0;
  }
  else if (StatusIter == DAMAGE)
  {
    Kntan = (1.0 - MaxDispIter/df)*(Ft/df)*exp(1.0 - MaxDispIter/df);
  }
  else  // Elastic unloading/reloading
  {
    Kntan = (1.0 - DamVarIter)*K0;
  }

  if (Anm->GetMecMod( )->GetDimQMatrix( ) == 2)
    C[1][1] = Kntan;
  else
    C[2][2] = Kntan;
}


// ================================= Stress ================================

int cModCZMModeIExponential :: Stress(cVector &eps, cVector &sig)
{
  // Initialization.

  sig.Zero( );

  // Get material parameters.

  int np = Mat->NumParam( );
  double *param = new double[np];
  Mat->GetParam(param);
  double Ft  = param[0];
  double GIc = param[1];
  double df  = GIc/(exp(1.0)*Ft);  // Separation corresponding to Ft
  double K0  = exp(1.0)*Ft/df;     // Initial stiffness
  delete []param;

  // Initialize the iterative material parameters (ParamIter = Param).

  DamVarIter  = DamVar;
  MaxDispIter = MaxDisp;

  // Get the opening displacement.

  double dn;
  if (Anm->GetMecMod( )->GetDimQMatrix( ) == 2)
    dn = eps[1];
  else
    dn = eps[2];

  // Stress computation.

  double Knsec;
  if (dn < 0.0) // Penalize negative displacements (penetration).
  {
    StatusIter = COMPRESSION;
    Knsec = K0;
  }
  else          // Tension
  {
    if (dn >= MaxDispIter)  // Damange evolution (loading)
    {
      StatusIter  = DAMAGE;
      MaxDispIter = dn;
      Knsec = (Ft/df)*exp(1.0 - dn/df);
      DamVarIter = 1.0 - Knsec/K0;
    }
    else                    // Elastic unloading/reloading
    {
      StatusIter = ELASTIC;
//      Knsec = (Ft/df)*exp(1.0 - MaxDispIter/df);
      Knsec = (1.0 - DamVarIter)*K0;
    }
  }

  // Handle 2D/3D models.

  if (Anm->GetMecMod( )->GetDimQMatrix( ) == 2)
    sig[1] = Knsec*dn;
  else
    sig[2] = Knsec*dn;

  return(1);
}

// =============================== UpdateState =============================

void cModCZMModeIExponential :: UpdateState(void)
{
  // Update material parameters.

  Status  = StatusIter;
  MaxDisp = MaxDispIter;
  DamVar  = DamVarIter;
}


// -------------------------------------------------------------------------
// Class cModPlastIsoLinHard:
// -------------------------------------------------------------------------

// =========================== cModPlastIsoLinHard =========================

cModPlastIsoLinHard :: cModPlastIsoLinHard(cElement *elm, cMaterial *mat) :
                       cConstModel(elm, mat)
{
  // Get and check the analysis model.

  Anm = elm->GetAnModel( );
  if (Anm->GetType( ) != PLANE_STRESS && Anm->GetType( ) != PLANE_STRAIN &&
      Anm->GetType( ) != AXISYMMETRIC && Anm->GetType( ) != SOLID)
  {
    cout << "Plasticity not implemented yet for this element!\n";
    exit(0);
  }

  // Initialize the internal variables.

  Status = ELASTIC;
  Dgam = 0.0;
  EqvPlEps = EqvPlEpsIter = 0.0;
  int nsc = Anm->GetMecMod( )->GetDimQMatrix( );
  PlastEps.Resize(nsc);
  PlastEps.Zero( );
  PlastEpsIter.Resize(nsc);
  PlastEpsIter.Zero( );
  StressIter.Resize(nsc);
  StressIter.Zero( );
}

// ========================== ~cModPlastIsoLinHard =========================

cModPlastIsoLinHard :: ~cModPlastIsoLinHard(void)
{
}

// ================================= Stress ================================

int cModPlastIsoLinHard :: Stress(cVector &eps, cVector &sig)
{
  // Get material parameters.

  int np = Mat->NumParam( );
  double *param = new double[np];
  Mat->GetParam(param);

  // Compute stresses.

  int code = StressPlast(param, eps, sig);

  // Release memory.

  delete []param;

  return(code);
}

// ================================= TangMat ===============================

void cModPlastIsoLinHard :: TangMat(cMatrix &C)
{
  // Get material parameters.

  int np = Mat->NumParam( );
  double *param = new double[np];
  Mat->GetParam(param);

  // Compute the tangent constitutive matrix.

  TangMatPlast(param, C);

  // Release memory.

  delete []param;
}

// =============================== UpdateState =============================

void cModPlastIsoLinHard :: UpdateState(void)
{
  EqvPlEps = EqvPlEpsIter;
  PlastEps = PlastEpsIter;
}

// -------------------------------------------------------------------------
// Protected methods:
//

// ============================== StressPlast ==============================

int cModPlastIsoLinHard :: StressPlast(double *param, cVector &eps, cVector &sig)
{
  // Handle the plane stress case.

  if (Anm->GetType( ) == PLANE_STRESS)
    return(StressPlastPS(param, eps, sig));

  // Get material parameters.

  double E     = param[0];
  double Nu    = param[1];
  double SigY0 = param[2];
  double H     = param[3];

  // Get the elastic matrix [C].

  int nsc = Anm->GetMecMod( )->GetDimQMatrix( );
  cMatrix C(nsc, nsc);
  Anm->GetMecMod( )->QMatrix(param, C);

  // Elastic predictor.

  Status = ELASTIC;
  Dgam = 0.0;
  PlastEpsIter = PlastEps;
  EqvPlEpsIter = EqvPlEps;

  // Evaluate the trial stresses.

  cVector eleps(nsc);          // Elastic strains
  eleps = eps - PlastEpsIter;  // Consider an elastic strain increment
  sig = C*eleps;

  // Evaluate the trial von Mises stress (qtrial).

  double J2 = Anm->GetMecMod( )->CalcJ2(sig);
  double qtrial = sqrt(3.0*J2);

  // Evaluate the current yield stress considering isotropic hardening.

  double SigY = SigY0 + H*EqvPlEpsIter; 

  // Check if the trial stresses are admissible (elastic step).

  double ftrial = qtrial - SigY;
  if (ftrial <= 0.0)
    return(1);

  // Plastic step - return mapping.

  Status = PLASTIC;
  double G = E/(2.0*(1.0 + Nu));  // Shear modulus
  Dgam = ftrial/(3.0*G + H);      // Plastic multiplier

  // Evaluate the trial hydrostatic (p) and deviatoric (sdev) trial stresses.

  double I1 = Anm->GetMecMod( )->CalcI1(sig);
  double p  = I1/3.0;
  cVector sdev(nsc);
  Anm->GetMecMod( )->CalcDev(p, sig, sdev);

  // Project the deviatoric stress on the plastic surface.

  sdev *= (1.0 - Dgam*3.0*G/qtrial);

  // Evaluate the elastoplastic stresses.

  Anm->GetMecMod( )->CalcTen(p, sdev, sig);
  StressIter = sig;               // Store the current stresses

  // Evaluate the plastic strains (N = sqrt(3/2)*sdev/||sdev||).

  J2 = Anm->GetMecMod( )->CalcJ2(sig);
  double devnrm = sqrt(2.0*J2);   // ||sdev||
  PlastEpsIter += Dgam*sqrt(3.0/2.0)/devnrm*sdev;
  EqvPlEpsIter += Dgam;

  return(1);
}

// ============================= StressPlastPS =============================

int cModPlastIsoLinHard :: StressPlastPS(double *param, cVector &eps, cVector &sig)
{
  // Get material parameters.

  double E     = param[0];
  double Nu    = param[1];
  double SigY0 = param[2];
  double H     = param[3];

  // Get the elastic matrix [C].

  cMatrix C(3,3);
  Anm->GetMecMod( )->QMatrix(param, C);

  // Elastic predictor.

  Status = ELASTIC;
  Dgam = 0.0;
  PlastEpsIter = PlastEps;
  EqvPlEpsIter = EqvPlEps;

  // Evaluate the trial stresses.

//  cout << "{eps} = " << scientific << endl;
//  eps.Print( );
  cVector eleps(3);            // Elastic strains
  eleps = eps - PlastEpsIter;  // Consider an elastic strain increment
  sig = C*eleps;
//  cout << "{sig} = " << scientific << endl;
//  sig.Print( );

  // Evaluate the trial yield function.

  double a1 = (sig[0] + sig[1])*(sig[0] + sig[1]);
  double a2 = (sig[0] - sig[1])*(sig[0] - sig[1]);
  double a3 = sig[2]*sig[2];
  double xi = a1/6.0 + a2/2.0 + 2.0*a3;
  double SigY = SigY0 + H*EqvPlEpsIter; 
//  cout << "xi = " << scientific << setprecision(4) << xi << "  ";
//  cout << "SigY = " << scientific << SigY << "  ";

  // Check if the trial stresses are admissible (elastic step).

  double f  = xi/2.0 - SigY*SigY/3.0;
  double ft = f/(SigY*SigY);
  double Tol = 1.0e-8;
//  cout << "ft = " << scientific << ft << endl;
  if (ft <= Tol)
    return(1);

  // Plastic step - return mapping.

  Status = PLASTIC;
  int MaxIter = 50;
  int conv = 0;
  double G = E/(2.0*(1.0 + Nu));  // Shear modulus
  double E3nu = E/(3.0*(1.0 - Nu));
  double b1 = 1.0;
  double b2 = 1.0;
  for (int i = 0; i < MaxIter; i++)
  {
    // Evaluate the yield function derivative.

    double dxi = -2.0*a1*E3nu/(6.0*b1*b1*b1) - 2.0*(a2/2.0 + 2.0*a3)*2.0*G/(b2*b2*b2);
    double rxi = sqrt(2.0*xi/3.0);
    double deq = rxi + Dgam*dxi/(3.0*rxi);
    double dY2 = 2.0*SigY*H*deq;
    double df = dxi/2.0 - dY2/3.0;

    // Evaluate the new plastic multiplier.

    Dgam -= f/df;
//    cout << "iter = " << i+1 << "  Dgam = " << scientific << Dgam << "  ";
//    cout << "f = " << f << "  ";

    // Update xi, the equivalent plastic strain and the yield stress.

    b1 = 1.0 + Dgam*E3nu;
    b2 = 1.0 + Dgam*2.0*G;
    xi = a1/(6.0*b1*b1) + (a2/2.0 + 2.0*a3)/(b2*b2);
//    cout << "xi = " << xi << "  ";
    EqvPlEpsIter = EqvPlEps + Dgam*sqrt(2.0*xi/3.0);
    SigY = SigY0 + H*EqvPlEpsIter; 
//    cout << "SigY = " << SigY << "  ";

    // Check for convergence.

    f  = xi/2.0 - SigY*SigY/3.0;
    ft = f/(SigY*SigY);
//    cout << "ft = " << scientific << ft << endl;
    if (fabs(ft) <= Tol)
    {
      conv = 1;
      break;
    }
  }
  if (!conv) return(0);

  // Evaluate the current stresses.

  cMatrix A(3,3);
  A.Zero( );
  double A11 = 1.0/b1;
  double A22 = 1.0/b2;
  A[0][0] = A[1][1] = (A11 + A22)/2.0; 
  A[0][1] = A[1][0] = (A11 - A22)/2.0; 
  A[2][2] = A22;
  StressIter = A*sig;
//  cout << "{StressIter} = " << scientific << endl;
//  StressIter.Print( );

  // Evaluate the plastic strains.

  cMatrix D(3,3);
  Anm->GetMecMod( )->DMatrix(param, D); // Elastic compliance matrix
  eleps = D*StressIter;                 // Elastic strain
  PlastEpsIter = eps - eleps;
//  cout << "{eleps} = " << scientific << endl;
//  eleps.Print( );
//  cout << "{PlastEpsIter} = " << scientific << endl;
//  PlastEpsIter.Print( );
//  cout << endl;

  // Return the computed stresses.

  sig = StressIter;
  return(1);
}

// ============================== TangMatPlast =============================

void cModPlastIsoLinHard :: TangMatPlast(double *param, cMatrix &C)
{
  // Evaluate the tangent constitutive matrix.

  if (Status == ELASTIC)  
  {
    // Compute the elastic [C] matrix using the Analysis Model class.

    Anm->GetMecMod( )->QMatrix(param, C);
  }
  else if (Anm->GetType( ) == PLANE_STRESS)
  {
    // Evaluate the shear (G) modulus.

    double E     = param[0];
    double Nu    = param[1];
    double SigY0 = param[2];
    double H     = param[3];
    double G = E/(2.0*(1.0 + Nu));

    // Assembly the [P] matrix.

    cMatrix P(3,3);
    P.Zero( );
    P[0][0] = P[1][1] =  2.0/3.0;
    P[0][1] = P[1][0] = -1.0/3.0;
    P[2][2] = 2.0;

    // Evaluate the first term of the tangent matrix.

    C.Zero( );
    double C11 = 3.0*E/(3.0*(1.0 - Nu) + E*Dgam);
    double C22 = 2.0*G/(1.0 + 2.0*G*Dgam);
    C[0][0] = C[1][1] = (C11 + C22)/2.0; 
    C[0][1] = C[1][0] = (C11 - C22)/2.0; 
    C[2][2] = C22/2.0;
//    cout << "----------------\n";
//    cout << "[C] = " << scientific << endl;
//    C.Print( );

    // Evaluate xi and {n}.

    cVector psig(3),vecn(3);
    psig = P*StressIter;
    vecn = C*psig;
    double xi = StressIter*psig;
//    cout << "Dgam = " << scientific << Dgam << "\n";
//    cout << "xi = " << scientific << setprecision(4) << xi << "  ";
//    cout << "{n} = " << scientific << endl;
//    vecn.Print( );

    // Evaluate alpha.

    double aux1 = vecn*psig;
    double aux2 = 2.0*xi*H/(3.0 - 2.0*H*Dgam);
    double alpha = 1.0/(aux1 + aux2);

    // Add the second term of the tangent matrix.

    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
        C[i][j] -= alpha*vecn[i]*vecn[j]; 
//    cout << "[C] = " << scientific << endl;
//    C.Print( );
  }
  else
  {
    // Evaluate the shear (G) and bulk (K) moduli.

    double E     = param[0];
    double Nu    = param[1];
    double SigY0 = param[2];
    double H     = param[3];
    double G = E/(2.0*(1.0 + Nu));
    double K = E/(3.0*(1.0 - 2.0*Nu));

    // Get the current stresses.

    int nsc = Anm->GetMecMod( )->GetDimQMatrix( );
    cVector sig(nsc);
    sig = StressIter;
    
    // Evaluate the hydrostatic (p) and deviatoric (sdev) stresses.

    double I1 = Anm->GetMecMod( )->CalcI1(sig);
    double p  = I1/3.0;
    cVector sdev(nsc);
    Anm->GetMecMod( )->CalcDev(p, sig, sdev);

    // Evaluate the trial von Mises stress (qtrial).

    double J2 = Anm->GetMecMod( )->CalcJ2(sig);
    double qtrial = sqrt(3.0*J2) + 3.0*G*Dgam;

    // Auxiliary matrices.

    cVector SO(nsc);
    cMatrix FO(nsc, nsc);
    Anm->GetMecMod( )->CalcSOID(SO);
    Anm->GetMecMod( )->CalcFOID(FO);
    cMatrix DevPrj(nsc, nsc);
    for (int i = 0; i < nsc; i++)
      for (int j = 0; j < nsc; j++)
         DevPrj[i][j] = FO[i][j] - SO[i]*SO[j]/3.0;

    // Consistent elastoplastic tangent matrix.

    double devnrm2 = 2.0*J2;  // ||sdev||^2
    double a = 2.0*G*(1.0 - 3.0*G*Dgam/qtrial);
    double b = 6.0*G*G*(Dgam/qtrial - 1.0/(3.0*G + H))/devnrm2;
    for (int i = 0; i < nsc; i++)
      for (int j = 0; j < nsc; j++)
        C[i][j] = a*DevPrj[i][j] + b*sdev[i]*sdev[j] + K*SO[i]*SO[j];
  }
}


// -------------------------------------------------------------------------
// Class cModFGMTTO:
// -------------------------------------------------------------------------

// =============================== cModFGMTTO ==============================

cModFGMTTO :: cModFGMTTO(cElement *elm, cMaterial *mat, double *p) :
              cModPlastIsoLinHard(elm, mat)
{
  V2 = p[0];
}

// ============================== ~cModFGMTTO ==============================

cModFGMTTO :: ~cModFGMTTO(void)
{
}

// ================================= Stress ================================

int cModFGMTTO :: Stress(cVector &eps, cVector &sig)
{
  // Get material parameters.

  int np = Mat->NumParam( );
  double *param = new double[np];
  Mat->GetParam(param);

  // Compute the equivalent properties (homogenization).

  double *fgmparam = new double[4];
  Homogenize(param, fgmparam);

  // Compute stresses.

  int code = StressPlast(fgmparam, eps, sig);

  // Release memory.

  delete []fgmparam;
  delete []param;

  return(code);
}

// ================================= CMatrix ===============================

void cModFGMTTO :: TangMat(cMatrix &C)
{
  // Get material parameters.

  int np = Mat->NumParam( );
  double *param = new double[np];
  Mat->GetParam(param);

  // FGM Homogenization.
  
  double *fgmparam = new double[4];
  Homogenize(param, fgmparam);

  // Compute the tangent constitutive matrix.

  TangMatPlast(fgmparam, C);

  // Release memory.

  delete []fgmparam;
  delete []param;
}

// -------------------------------------------------------------------------
// Protected methods:
//

// ============================== Homogenize ===============================

void cModFGMTTO :: Homogenize(double *param, double *fgmparam)
{
  // Get the material parameters.

  double E1   = param[0];
  double Nu1  = param[1];
  double E2   = param[2];
  double Nu2  = param[3];
  int method  = int(param[4]);
  double SigYm = param[5];
  double Hm = param[6];
  double qtto = param[7];

  // Homogenization.

  double E, Nu;
  if (method == 1) // Rule of Mixtures.
  {
    E  = (E2 - E1)*V2 + E1;
    Nu = (Nu2 - Nu1)*V2 + Nu1;
  }
  else            // Mori-Tanaka
  {
    // Bulk and shear moduli of matrix and inclusions.

    double K1 = E1/(3.0 - 6.0*Nu1);
    double K2 = E2/(3.0 - 6.0*Nu2);
    double G1 = E1/(2.0 + 2.0*Nu1);
    double G2 = E2/(2.0 + 2.0*Nu2);

    // Auxiliary factors.

    double V1 = 1.0 - V2;
    double f1 = G1*(9.0*K1 + 8.0*G1)/(6.0*K1 + 12.0*G1);
    double a  = (K2 - K1)/(K1 + 4.0*G1/3.0);
    double b  = (G2 - G1)/(G1 + f1);

    // Mori-Tanaka homogenized bulk and shear moduli.

    double K = K1 + (K2 - K1)*V2/(1.0 + V1*a);
    double G = G1 + (G2 - G1)*V2/(1.0 + V1*b);

    // Transform back into Elastic Modulus and Poisson's ratio.

    E  = 9.0*K*G/(3.0*K + G);
    Nu = (3.0*K - 2.0*G)/(6.0*K + 2.0*G);
  }

  // Compute the TTO effective properties.
  
  double V1 = 1.0 - V2;
  double SigY, H, Ep, Epm;
  if (Hm == 0)
  {
    Epm = E1;
  }
  else
  {
    Epm = E1*Hm/(E1 + Hm);
  }
  SigY = SigYm*(V1 + ((qtto + E1)/(qtto + E2))*(E2/E1)*V2);
  Ep = (V1*Epm*((qtto + E2)/(qtto + Epm)) + V2*E2)/(V1*((qtto + E2)/(qtto + Epm)) + V2);
  if (Hm == 0)
  {
    H = 0;
  }
  else
  {
    H = E*Ep/(E - Ep);
  }
  
  // Assembly the FGM param vector.

  fgmparam[0] = E;    // Combined elasticity modulus
  fgmparam[1] = Nu;   // Combined Poisson's ratio
  fgmparam[2] = SigY; // Combined yield stress
  fgmparam[3] = H;    // Combined hardening modulus
}

// ======================================================= End of file =====
