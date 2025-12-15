// -------------------------------------------------------------------------
// cmodel1d.cpp - implementation of uniaxial constitutive models.
// -------------------------------------------------------------------------
// Created:      06-Nov-2012     Evandro Parente Junior
//
// Modified:     04-Aug-2014     Carlos David Rodrigues Melo
//               Implementation of Elastoplastic models.
//
// Modified:     08-Apr-2015     Elias Saraiva Barroso
//               Input/output operations using c++ streams.
//
// Modified:     02-Mar-2020     Bergson Matias/Evandro Parente
//               Creation of Mazars, MuModel and Lee-Fenves models.
// -------------------------------------------------------------------------

#include <math.h>
#include <iostream>
#include <iomanip>
#include <stdlib.h>

using namespace std;

#include "cmodel1d.h"
#include "material.h"
#include "ctrl.h"
#include "vec.h"
#include "mat.h"
#include "gbldef.h"

// Gambiarra berg
#include <fstream>

// -------------------------------------------------------------------------
// Public methods:
//

// -------------------------------------------------------------------------
// Class cMod1DLinear:
// -------------------------------------------------------------------------

// ============================== cMod1DLinear =============================

cMod1DLinear :: cMod1DLinear(cMaterial *mat) : cConstModel(mat)
{
}

// ============================= ~cMod1DLinear =============================

cMod1DLinear :: ~cMod1DLinear(void)
{
}

// ================================= Stress ================================

int cMod1DLinear :: Stress(cVector &eps, cVector &sig)
{
  // Get material parameters.

  int np = Mat->NumParam( );
  double *param = new double[np];
  Mat->GetParam(param);
  double E = param[0];

  // Compute the total stress and return success.

  sig[0] = E*eps[0];

  // Release memory and return success.

  delete []param;
  return(1);
}

// ================================= TangMat ===============================

void cMod1DLinear :: TangMat(cMatrix &C)
{
  // Get material parameters.

  int np = Mat->NumParam( );
  double *param = new double[np];
  Mat->GetParam(param);
  double E = param[0];
  delete []param;

  // Assembly the constitutive matrix.

  C[0][0] = E;
}


// -------------------------------------------------------------------------
// Class cMod1DElastNonlin:
// -------------------------------------------------------------------------

// ============================ cMod1DElastNonlin ==========================

cMod1DElastNonlin :: cMod1DElastNonlin(cMaterial *mat) : cConstModel(mat)
{
  Eps = 0.0;
}

// =========================== ~cMod1DElastNonlin ==========================

cMod1DElastNonlin :: ~cMod1DElastNonlin(void)
{
}

// ============================== NumStrainLim =============================

int cMod1DElastNonlin :: NumStrainLim(void)
{
  // Get material pointer.

  cMatElastNonlin *mat = (cMatElastNonlin *)Mat;

  // Number of divisions = number of points - 1.

  return(mat->NumPnt - 2);
}

// ============================== GetStrainLim =============================

void cMod1DElastNonlin :: GetStrainLim(cVector &lim)
{
  // Get material pointer.

  cMatElastNonlin *mat = (cMatElastNonlin *)Mat;

  // Fill the vector with the strain limits.

  for (int i = 0; i < mat->NumPnt-2; i++)
    lim[i] = mat->Eps[i+1];
}

// ================================= Stress ================================

int cMod1DElastNonlin :: Stress(cVector &eps, cVector &sig)
{
  // Get material pointer.

  cMatElastNonlin *mat = (cMatElastNonlin *)Mat;

  // Find interval.

  int i;
  for (i = 1; i < mat->NumPnt; i++)
    if (eps[0] <= mat->Eps[i]) break;  // Interpolation
  if (i == mat->NumPnt) i--;           // Extrapolation

  // Compute the stress by linear interpolation (or extrapolation).

  double de = mat->Eps[i] - mat->Eps[i-1];
  double ds = mat->Sig[i] - mat->Sig[i-1];
  double Et = 0.0;
  if (de != 0.0) Et = ds/de;
  sig[0] = mat->Sig[i-1] + Et*(eps[0] - mat->Eps[i-1]);

  // Store the current strain.

  Eps = eps[0];
  return(1);
}

// ================================= TangMat ===============================

void cMod1DElastNonlin :: TangMat(cMatrix &C)
{
  // Get material pointer.

  cMatElastNonlin *mat = (cMatElastNonlin *)Mat;

  // Find interval.

  int i;
  for (i = 1; i < mat->NumPnt; i++)
    if (Eps <= mat->Eps[i]) break;  // Interpolation
  if (i == mat->NumPnt) i--;        // Extrapolation

  // Compute the tangent modulus.

  double de = mat->Eps[i] - mat->Eps[i-1];
  double ds = mat->Sig[i] - mat->Sig[i-1];
  double Et = 0.0;
  if (de != 0.0) Et = ds/de;
  C[0][0] = Et;
}


// -------------------------------------------------------------------------
// Class cMod1DViscoElastic:
// -------------------------------------------------------------------------

// ============================ cMod1DViscoElastic =========================

cMod1DViscoElastic :: cMod1DViscoElastic(cMaterial *mat) : cConstModel(mat)
{
  Time = 0.0;
  Sig  = 0.0;
  Eps  = 0.0;
  SigIter = 0.0;
  EpsIter = 0.0;
  EpsRate = 0.0;

  int np = Mat->NumParam( );
  int nprony = (np - 2)/2;
  S = new double[nprony];
  for (int i = 0; i < nprony; i++) S[i] = 0.0;
}

// =========================== ~cMod1DViscoElastic =========================

cMod1DViscoElastic :: ~cMod1DViscoElastic(void)
{
  delete []S;
}

// ================================= Stress ================================

int cMod1DViscoElastic :: Stress(cVector &eps, cVector &sig)
{
  // Start with the stored stresses.

  sig[0] = Sig;
//  double t  = cControl::GetTotTime( );
  double dt = cControl::GetTimeStep( );

  // Get material parameters.

  int np = Mat->NumParam( );
  double *param = new double[np];
  Mat->GetParam(param);
  int nprony = (np - 2)/2;

  // Compute dsigbar = Ebar*deps (Obs: Ebar = Et).

  int i;
  double deps = eps[0] - Eps;
  double E,rho;
  double Et = param[0];
  for (i = 0; i < nprony; i++)
  {
    E   = param[2*i+2];
    rho = param[2*i+3];
    Et += E*(rho/dt)*(1.0 - exp(-dt/rho));
  }
  double dsigbar = Et*deps;
  sig[0] += dsigbar;

  // Compute dsighat => hereditary terms.

  double s;
  double dsighat = 0.0;
  for (i = 0; i < nprony; i++)
  {
    E   = param[2*i+2];
    rho = param[2*i+3];
    s   = exp(-dt/rho)*S[i];
    s  += E*rho*(1.0 - exp(-dt/rho))*EpsRate;
    dsighat += -(1.0 - exp(-dt/rho))*s;
  }
  sig[0] += dsighat;

  // Store the iterative vectors.

  EpsIter = eps[0];
  SigIter = sig[0];

  // Release memory and return success.

  delete []param;
  return(1);
}

// ================================= TangMat ===============================

void cMod1DViscoElastic :: TangMat(cMatrix &C)
{
  // Get material parameters.

  int np = Mat->NumParam( );
  double *param = new double[np];
  Mat->GetParam(param);
  int nprony = (np - 2)/2;

  // Evaluate the tangent relaxation modulus.

  double dt = cControl::GetTimeStep( );
  double Et = param[0];   // Infinite modulus
  double E,rho;
  for (int i = 0; i < nprony; i++)
  {
    E   = param[2*i+2];
    rho = param[2*i+3];
    Et += E*(rho/dt)*(1.0 - exp(-dt/rho));
  }
  delete []param;

  // Evaluate [C].

  C[0][0] = Et;
}

// =============================== UpdateState =============================

void cMod1DViscoElastic :: UpdateState(void)
{
  // Get material parameters.

  int np = Mat->NumParam( );
  double *param = new double[np];
  Mat->GetParam(param);
  int nprony = (np - 2)/2;

  // Update the hereditary terms => S.

  double E,rho;
  double dt = cControl::GetTimeStep( );
  for (int i = 0; i < nprony; i++)
  {
    E   = param[2*i+2];
    rho = param[2*i+3];
    S[i] *= exp(-dt/rho);
    S[i] += E*rho*(1.0 - exp(-dt/rho))*EpsRate;
  }

  // Update the internal variables.

  Time = cControl::GetTotTime( );
  EpsRate = (EpsIter - Eps)/dt;
  Eps = EpsIter;
  Sig = SigIter;

  // Release memory.

  delete []param;
}


// -------------------------------------------------------------------------
// Class cMod1DConcEuroCEB:
// -------------------------------------------------------------------------

// ============================ cMod1DConcEuroCEB ==========================

cMod1DConcEuroCEB :: cMod1DConcEuroCEB(cMaterial *mat) : cConstModel(mat)
{
  Eps = 0.0;
}

// =========================== ~cMod1DConcEuroCEB ==========================

cMod1DConcEuroCEB :: ~cMod1DConcEuroCEB(void)
{
}

// ============================== NumStrainLim =============================

int cMod1DConcEuroCEB :: NumStrainLim(void)
{
  // Get material pointer.

  cMatConcEuroCEB *mat = (cMatConcEuroCEB *)Mat;

  // Number of divisions.

  if (mat->fct <= 0.0)
    return(2);
  else
    return(4);
}

// ============================== GetStrainLim =============================

void cMod1DConcEuroCEB :: GetStrainLim(cVector &lim)
{
  // Get material pointer.

  cMatConcEuroCEB *mat = (cMatConcEuroCEB *)Mat;

  // Fill the vector with the strain limits.

  lim[0] = -3.50e-03;
  lim[1] = 0.0;
  if (mat->fct > 0.0)
  {
    lim[2] = mat->fct/mat->Eci; // Cracking strain
    lim[3] = mat->ey;
  }
}

// ============================== GetCrvDegree =============================

int cMod1DConcEuroCEB :: GetCrvDegree(int range)
{
  // Return the degree of the polynomial describing the stress-strain curve.

  if (range == 1 || range == 3)
    return(4);          // Equivalent degree (non polynomial curve)
  else if (range == 2)
    return(1);
  else
    return(0);
}

// ================================= Stress ================================

int cMod1DConcEuroCEB :: Stress(cVector &eps, cVector &sig)
{
  // Get material pointer.

  cMatConcEuroCEB *mat = (cMatConcEuroCEB *)Mat;

  // Define strain intervals for EuroCEB curve.

  double e1 = -3.50e-03;
  double e2 = 0.0;
  double e3 = mat->fct/mat->Eci;
  double e4 = mat->ey;

  // Compute the stress using the EuroCEB curve.

  double e = eps[0];
  if (e >= e1 && e <= e2)                        // Compression
  {
//    double k = -1.1*mat->Ecm*mat->ec1/mat->fcm;
    double k = -1.05*mat->Ecm*mat->ec1/mat->fcm;
    double r = e/mat->ec1;
    sig[0] = -mat->fcm*((k*r - r*r)/(1.0 + k*r - 2.0*r));
  }
  else if (e >= e2 && e < e3 && mat->fct > 0)    // Tension pre-cracking
  {
    sig[0] = mat->Eci*e;
  }
  else if (e >= e3 && e <= e4 && mat->fct > 0)   // Tension post-cracking
  {
    double a = 0.5*mat->Rho*mat->Es*e;
    double b = mat->fct*mat->fct*(1.0 + mat->Rho*mat->Es/mat->Eci);
    sig[0] = -a + sqrt(a*a + b);
  }
  else
    sig[0] = 0.0;

  // Store the current strain.

  Eps = eps[0];
  return(1);
}

// ================================= TangMat ===============================

void cMod1DConcEuroCEB :: TangMat(cMatrix &C)
{
  // Get material pointer.

  cMatConcEuroCEB *mat = (cMatConcEuroCEB *)Mat;

  // Define strain intervals for EuroCEB curve.

  double e1 = -3.50e-03;
  double e2 = 0.0;
  double e3 = mat->fct/mat->Eci;
  double e4 = mat->ey;

  // Compute the tangent modulus according to EuroCEB model.

  double Et = 0.0;
  if (Eps >= e1 && Eps <= e2)                        // Compression
  {
    double k = -1.05*mat->Ecm*mat->ec1/mat->fcm;
//    double k = -1.1*mat->Ecm*mat->ec1/mat->fcm;
    double r = Eps/mat->ec1;
    double a = 1.0 + (k - 2.0)*r;
    double b = (k - 2.0)*(k*r - r*r)/(a*a);
    double c = (k - 2.0*r)/a;
    Et = (mat->fcm/mat->ec1)*(b - c);
  }
  else if (Eps >= e2 && Eps < e3 && mat->fct > 0)    // Tension pre-cracking
  {
    Et = mat->Eci;
  }
  else if (Eps >= e3 && Eps <= e4 && mat->fct > 0)   // Tension post-cracking
  {
    double da = 0.5*mat->Rho*mat->Es;
    double a  = da*Eps;
    double b  = mat->fct*mat->fct*(1.0 + mat->Rho*mat->Es/mat->Eci);
    Et = -da + a*da/sqrt(a*a + b);
  }

  C[0][0] = Et;
}


// -------------------------------------------------------------------------
// Class cMod1DConcNBRCEB:
// -------------------------------------------------------------------------

// ============================ cMod1DConcNBRCEB ===========================

cMod1DConcNBRCEB :: cMod1DConcNBRCEB(cMaterial *mat) : cConstModel(mat)
{
  Eps = 0.0;
}

// =========================== ~cMod1DConcNBRCEB ===========================

cMod1DConcNBRCEB :: ~cMod1DConcNBRCEB(void)
{
}

// ============================== NumStrainLim =============================

int cMod1DConcNBRCEB :: NumStrainLim(void)
{
  // Get material pointer.

  cMatConcNBRCEB *mat = (cMatConcNBRCEB *)Mat;

  // Number of divisions.

  if (mat->fct <= 0.0)
    return(3);
  else
    return(5);
}

// ============================== GetStrainLim =============================

void cMod1DConcNBRCEB :: GetStrainLim(cVector &lim)
{
  // Get material pointer.

  cMatConcNBRCEB *mat = (cMatConcNBRCEB *)Mat;

  // Fill the vector with the strain limits.

  lim[0] = -3.5e-03;
  lim[1] = -2.0e-03;
  lim[2] = 0.0;
  if (mat->fct > 0.0)
  {
    lim[3] = mat->fct/mat->Eci;
    lim[4] = mat->ey;
  }
}

// ============================== GetCrvDegree =============================

int cMod1DConcNBRCEB :: GetCrvDegree(int range)
{
  // Return the degree of the polynomial describing the stress-strain curve.

  if (range == 2)
    return(2);
  else if (range == 3)
    return(1);
  else if (range == 4)
    return(4);          // Equivalent degree (non polynomial curve)
  else
    return(0);
}

// ================================= Stress ================================

int cMod1DConcNBRCEB :: Stress(cVector &eps, cVector &sig)
{
  // Get material pointer.

  cMatConcNBRCEB *mat = (cMatConcNBRCEB *)Mat;

  // Define strain intervals for NBR-CEB curve.

  double e1 = -3.5e-03;
  double e2 = -2.0e-03;
  double e3 =  0.0;
  double e4 =  mat->fct/mat->Eci;
  double e5 =  mat->ey;

  // Compute stress using NBR-CEB curve.

  double e = eps[0];
  if (e >= e1 && e < e2)                         // Compression (rect.)
  {
    sig[0] = -mat->fc;
  }
  else if (e >= e2 && e < e3)                    // Compression (parab.)
  {
    sig[0] = -mat->fc*(1.0 - pow(1.0 - e/e2, 2));
  }
  else if (e >= e3 && e < e4 && mat->fct > 0)    // Tension pre-cracking
  {
    sig[0] = mat->Eci*e;
  }
  else if (e >= e4 && e <= e5 && mat->fct > 0)   // Tension post-cracking
  {
    double a = 0.5*mat->Rho*mat->Es*e;
    double b = mat->fct*mat->fct*(1.0 + mat->Rho*mat->Es/mat->Eci);
    sig[0] = -a + sqrt(a*a + b);
  }
  else
  {
    sig[0] = 0.0;
  }

  // Store the current strain.

  Eps = eps[0];
  return(1);
}

// ================================= TangMat ===============================

void cMod1DConcNBRCEB :: TangMat(cMatrix &C)
{
  // Get material pointer.

  cMatConcNBRCEB *mat = (cMatConcNBRCEB *)Mat;

  // Define strain intervals for NBR CEB curve.

//  double e1 = -3.5e-03;
  double e2 = -2.0e-03;
  double e3 =  0.0;
  double e4 =  mat->fct/mat->Eci;
  double e5 =  mat->ey;

  // Compute the tangent modulus according to NBR CEB model.

  double Et = 0.0;
  if (Eps >= e2 && Eps < e3)                         // Compression (parab.)
  {
    Et = -2.0*(mat->fc/e2)*(1.0 - Eps/e2);
  }
  else if (Eps >= e3 && Eps < e4 && mat->fct > 0)    // Tension pre-cracking
  {
    Et = mat->Eci;
  }
  else if (Eps >= e4 && Eps <= e5 && mat->fct > 0)   // Tension post-cracking
  {
    double da = 0.5*mat->Rho*mat->Es;
    double a  = da*Eps;
    double b  = mat->fct*mat->fct*(1.0 + mat->Rho*mat->Es/mat->Eci);
    Et = -da + a*da/sqrt(a*a + b);
  }

  C[0][0] = Et;
}


// -------------------------------------------------------------------------
// Class cMod1DPlastLinHard:
// -------------------------------------------------------------------------

// =========================== cMod1DPlastLinHard ==========================

cMod1DPlastLinHard :: cMod1DPlastLinHard(cMaterial *mat) :
                         cConstModel(mat)
{
  // Initialize the internal variables.

  Status = ELASTIC;
  PlasticEpsIter = PlasticEps = 0.0;
  AlphaIter = Alpha = 0.0;
  BackStress = BackStressIter = 0.0;
}

// ========================== ~cMod1DPlastLinHard ==========================

cMod1DPlastLinHard :: ~cMod1DPlastLinHard(void)
{
}

// ================================= Stress ================================

int cMod1DPlastLinHard :: Stress(cVector &eps, cVector &sig)
{
  // Get material parameters.

  int np = Mat->NumParam( );
  double *param = new double[np];
  Mat->GetParam(param);
  double E    = param[0];
//  double Nu   = param[1];
  double SigY = param[2];
  double K    = param[3];
  double H    = param[4];
  delete []param;

  // Elastic preditor.

  PlasticEpsIter = PlasticEps;
  AlphaIter = Alpha;
  BackStressIter = BackStress;
  double sigpred = E*(eps[0] - PlasticEpsIter);
  double xsipred = sigpred - BackStressIter;
  double fpred = fabs(xsipred) - (SigY + K*AlphaIter);

  // Stress computation.

  if (fpred <= 0.0) // Elastic increment
  {
    Status = ELASTIC;
    sig[0] = sigpred;
  }
  else              // Plastic increment
  {
    Status = PLASTIC;
    double Dgam = fpred/(E + K + H);
    PlasticEpsIter += Dgam*Sign(xsipred);
    AlphaIter += Dgam;
    BackStressIter += Dgam*H*Sign(xsipred);
    sig[0] = E*(eps[0] - PlasticEpsIter);
  }

  return(1);
}

// ================================= TangMat ===============================

void cMod1DPlastLinHard:: TangMat(cMatrix &C)
{
  // Get material parameters.

  int np = Mat->NumParam( );
  double *param = new double[np];
  Mat->GetParam(param);
  double E    = param[0];
//  double Nu   = param[1];
//  double SigY = param[2];
  double K    = param[3];
  double H    = param[4];
  delete []param;

  // Compute the tangent modulus.

  C.Zero( );
  if (Status == ELASTIC)
  {
    C[0][0] = E;
  }
  else
  {
    C[0][0] = E*(K + H)/(E + K + H);
  }
}

// =============================== UpdateState =============================

void cMod1DPlastLinHard :: UpdateState(void)
{
  // Update the model parameters.

  Alpha = AlphaIter;
  PlasticEps = PlasticEpsIter;
  BackStress = BackStressIter;
}


// -------------------------------------------------------------------------
// Class cMod1DPlastIsoNonlinHard:
// -------------------------------------------------------------------------

// ========================= cMod1DPlastIsoNonlinHard ======================

cMod1DPlastIsoNonlinHard :: cMod1DPlastIsoNonlinHard(cMaterial *mat) :
                            cConstModel(mat)
{
  // Initialize the internal variables.

  Status = ELASTIC;
  PlasticEpsIter = PlasticEps = 0.0;
  AlphaIter = Alpha = 0.0;
}

// ======================== ~cMod1DPlastIsoNonlinHard ======================

cMod1DPlastIsoNonlinHard :: ~cMod1DPlastIsoNonlinHard(void)
{
}

// ================================= Stress ================================

int cMod1DPlastIsoNonlinHard :: Stress(cVector &eps, cVector &sig)
{
  // Get material parameters.

  int np = Mat->NumParam( );
  double *param = new double[np];
  Mat->GetParam(param);
  double E    = param[0];
  double SigY = param[1];
  double K    = param[2];
  double ExpH = param[3];
  delete []param;

  // Elastic preditor.

  PlasticEpsIter = PlasticEps;
  AlphaIter = Alpha;
  double DeltaGama = 0.0;
  double sigpred = E*(eps[0] - PlasticEpsIter);
  double fpred = fabs(sigpred) - SigY*(1 + K*pow(AlphaIter,ExpH));

  cout << "eps[0] = " << eps[0] << " PlasticEpsIter = " << PlasticEpsIter << "\n";

  if (fpred <= 0.0) // Elastic increment
  {
    Status = ELASTIC;
    sig[0] = sigpred;
  }
  else              // Plastic increment
  {
    Status = PLASTIC;

    // Local N-R algorithm.

    int NITER = 30;            // Maximum number of steps to plastic N-R
    double TOL = 1e-6;         // Tolerance
    for (int i = 0; i < NITER; i++)
    {
      double r = fabs(sigpred) - DeltaGama*E - SigY*(1 + K*pow(AlphaIter + DeltaGama,ExpH)); //Residual.
      if(fabs(r) < TOL) break;       //Evaluate if DdeltaGamaIter satisfact the conditions
      double d;                // correction of the DeltaGama
      double dr = -E - SigY*K*ExpH*pow(AlphaIter + DeltaGama,ExpH-1); //Derived of the residual.
      r *= -1;
      d = r/dr;                //Solver the sistem {dr}{d} = -{r}
      if (d < TOL) break;      //Evaluate the size of the step.
      DeltaGama += d;          // Increment the DdeltaGamaIter
    }

    // Recalculate the variables

    AlphaIter += DeltaGama;               // Alpha
    int sign = 1;
    if (sigpred < 0) sign = -1;
    PlasticEpsIter += DeltaGama*sign;     // PlasticEps
    sig[0] = E*(eps[0] - PlasticEpsIter); //Stress
  }
return (1);
}

// ================================= TangMat ===============================

void cMod1DPlastIsoNonlinHard:: TangMat(cMatrix &C)
{
  // Get material parameters.

  int np = Mat->NumParam( );
  double *param = new double[np];
  Mat->GetParam(param);
  double E    = param[0];
//  double SigY = param[1];
  double K    = param[2];
  double ExpH = param[3];
  delete []param;

  // Compute the tangent modulus.

  C.Zero( );
  if (Status == ELASTIC)
  {
    C[0][0] = E;
  }
  else
  {
    C[0][0] = E*K*ExpH*pow(AlphaIter,ExpH-1) / (E + K*ExpH*pow(AlphaIter,ExpH-1));
  }
}

// =============================== UpdateState =============================

void cMod1DPlastIsoNonlinHard :: UpdateState(void)
{
  // Update the model parameters.

  Alpha = AlphaIter;
  PlasticEps = PlasticEpsIter;
}


// -------------------------------------------------------------------------
// Class cMod1DPlastIsoExpHard:
// -------------------------------------------------------------------------

// ========================= cMod1DPlastIsoExpHard =========================

cMod1DPlastIsoExpHard :: cMod1DPlastIsoExpHard(cMaterial *mat) :
                              cConstModel(mat)
{
  // Initialize the internal variables.

  Status = ELASTIC;
  PlasticEpsIter = PlasticEps = 0.0;
  AlphaIter = Alpha = 0.0;
}

// ======================== ~cMod1DPlastExpLinHard =========================

cMod1DPlastIsoExpHard :: ~cMod1DPlastIsoExpHard(void)
{
}

// ================================= Stress ================================

int cMod1DPlastIsoExpHard :: Stress(cVector &eps, cVector &sig)
{
  // Get material parameters.

  int np = Mat->NumParam( );
  double *param = new double[np];
  Mat->GetParam(param);
  double E    = param[0];
  double SigY = param[1];
  double Mu   = param[2];
//  cout << scientific << setprecision(3);
//  cout << "E = " << E << " SigY = " << SigY << " Mu = " << Mu << "\n";
//  cout << resetiosflags(ios::scientific) << setprecision(6);
  delete []param;

  // Elastic preditor.

  PlasticEpsIter = PlasticEps;
  AlphaIter = Alpha;
  sig[0] = E*(eps[0] - PlasticEpsIter);
  double f = fabs(sig[0]) - SigY*exp(Mu*AlphaIter);
  if (f <= 0.0)
  {
    Status = ELASTIC;
    return(1);
  }

  // Plastic increment.

  Status = PLASTIC;
  int MAXITER = 30;          // Maximum number of steps to plastic N-R
  double TOL = 1.0e-8;       // Tolerance
  double Dgam = 0.0;         // DeltaGamma (step)
  for (int i = 1; i <= MAXITER; i++)
  {
    // Evalute the iterative plastic multiplier dgam.

    double dG   = SigY*Mu*exp(Mu*AlphaIter);
    double df   = -(E + dG);
    double dgam = -f/df;

//    cout << scientific << setprecision(6);
//    cout << "i = " << i << " eps = " << eps[0] << " sig = " << sig[0];
//    cout << " f = " << setprecision(3) << f << " ";
//    cout << resetiosflags(ios::scientific) << setprecision(6);

    // Update the plastic variables.

    Dgam += dgam;
    AlphaIter += dgam;
    PlasticEpsIter += dgam*Sign(sig[0]);

    // Compute the updated stress and yield function.

    sig[0] = E*(eps[0] - PlasticEpsIter);
    f = fabs(sig[0]) - SigY*exp(Mu*AlphaIter);

//    cout << scientific << setprecision(6);
//    cout << "dgam = " << dgam << " sig = " << sig[0];
//    cout << " f = " << showpos << setprecision(3) << f;
//    cout << noshowpos << "err = " << f/SigY << "\n";
//    cout << resetiosflags(ios::scientific) << setprecision(6);

    // Test for convergence.

    if (fabs(f)/SigY <= TOL && fabs(dgam/Dgam) <= TOL) return(1);
  }

  // Convergence not achieved.

  return(0);
}

// ================================= TangMat ===============================

void cMod1DPlastIsoExpHard:: TangMat(cMatrix &C)
{
  // Get material parameters.

  int np = Mat->NumParam( );
  double *param = new double[np];
  Mat->GetParam(param);
  double E    = param[0];
  double SigY = param[1];
  double Mu   = param[2];
  delete []param;

  // Compute the tangent modulus.

  C.Zero( );
  if (Status == ELASTIC)
  {
    C[0][0] = E;
  }
  else
  {
    double dG = SigY*Mu*exp(Mu*AlphaIter);
    C[0][0] = E*dG/(E + dG);
  }
}

// =============================== UpdateState =============================

void cMod1DPlastIsoExpHard :: UpdateState(void)
{
  // Update the model parameters.

  Alpha = AlphaIter;
  PlasticEps = PlasticEpsIter;
}


// -------------------------------------------------------------------------
// Class cMod1DDamageMazars:
// -------------------------------------------------------------------------

// ========================= cMod1DDamageMazars =========================

cMod1DDamageMazars :: cMod1DDamageMazars(cMaterial *mat) :
                      cConstModel(mat)
{
  // Initialize the internal variables.

  Status = ELASTIC_TENSION;
  DtIter = Dt = 0.0;
  DcIter = Dc = 0.0;

  int np = Mat->NumParam( );
  double *param = new double[np];
  Mat->GetParam(param);
  double epsd0= param[5];
  delete []param;

  EpsEqtHist = EpsEqcHist = epsd0;
  EpsEqtHistIter = EpsEqcHistIter = 0.0;
}

// ======================== ~cMod1DDamageMazars =========================

cMod1DDamageMazars :: ~cMod1DDamageMazars(void)
{
}

// ================================= Stress ================================

int cMod1DDamageMazars :: Stress(cVector &eps, cVector &sig)
{
  // Get material parameters.

  int np = Mat->NumParam( );
  double *param = new double[np];
  Mat->GetParam(param);
  double E    = param[0];
  double Ac   = param[1];
  double Bc   = param[2];
  double At   = param[3];
  double Bt   = param[4];
  double epsd0= param[5];
  double Nu   = param[6];
  delete []param;

  // Set the trial values of the internal variables.

  DtIter = Dt;
  DcIter = Dc;
  EpsEqtHistIter = EpsEqtHist;
  EpsEqcHistIter = EpsEqcHist;
  //EpsEqtHistIter = MAX(EpsEqtHist, epsd0); // ESTRANHO
  //EpsEqcHistIter = MAX(EpsEqcHist, epsd0); // ESTRANHO

  // Determine the loading status and update the damage parameters.

  if (eps[0] >= 0.0) // Tension
  {
    double epseqt = eps[0];
    if (eps[0] <= EpsEqtHistIter)  // Elastic loading/unloading
    {
      Status = ELASTIC_TENSION;
    }
    else                          // Damage evolution
    {
      Status = DAMAGE_TENSION;
      DtIter = 1.0 - epsd0*(1.0 - At)/epseqt - At/(exp(Bt*(epseqt - epsd0)));
      EpsEqtHistIter = epseqt;
      if (DtIter < 0.0)
      {
        DtIter = 0.0;
        EpsEqtHistIter = EpsEqtHist;
      }
      if (DtIter > 1.0)
      {
        DtIter = 1.0;
        EpsEqtHistIter = EpsEqtHist;
      }
      //if (DtIter < 0.0 || DtIter > 1.0)  // ESTRANHO
      //{
      //  DtIter = Dt;
      //  EpsEqtHistIter = EpsEqtHist;
      //}
    }
  }
  else               // Compression
  {
    double epseqc = sqrt(2.0)*Nu*fabs(eps[0]);
    //double epseqc = fabs(eps[0]);
    if (epseqc <= EpsEqcHistIter)  // Elastic loading/unloading
    {
      Status = ELASTIC_COMPRESSION;
    }
    else                          // Damage evolution
    {
      Status = DAMAGE_COMPRESSION;
      DcIter = 1.0 - epsd0*(1.0 - Ac)/epseqc - Ac/(exp(Bc*(epseqc - epsd0)));
      EpsEqcHistIter = epseqc;
      if (DcIter < 0.0)
      {
        DcIter = 0.0;
        EpsEqcHistIter = EpsEqcHist;
      }
      if (DcIter > 1.0)
      {
        DcIter = 1.0;
        EpsEqcHistIter = EpsEqcHist;
      }
      //if (DcIter < 0.0 || DcIter > 1.0)  // ESTRANHO
      //{
      //  DcIter = Dc;
      //  EpsEqcHistIter = EpsEqcHist;
      //}
    }
  }

  // Evaluate the current stress using the secant modulus.

  if (eps[0] >= 0.0)
  {
    sig[0] = (1.0 - DtIter)*E*eps[0];
  }
  else
  {
    sig[0] = (1.0 - DcIter)*E*eps[0];
  }

  return(1);
}

// ================================= TangMat ===============================

void cMod1DDamageMazars:: TangMat(cMatrix &C)
{
  // Get material parameters.

  int np = Mat->NumParam( );
  double *param = new double[np];
  Mat->GetParam(param);
  double E    = param[0];
  double Ac   = param[1];
  double Bc   = param[2];
  double At   = param[3];
  double Bt   = param[4];
  double epsd0= param[5];
  double Nu   = param[6];
  delete []param;

  // Compute the tangent modulus.

  C.Zero( );
  if (Status == ELASTIC_TENSION)
  {
    C[0][0] = (1.0 - DtIter)*E;
  }
  else if (Status == ELASTIC_COMPRESSION)
  {
    C[0][0] = (1.0 - DcIter)*E;
  }
  else if (Status == DAMAGE_TENSION)
  {
    double aux1 = epsd0*(1.0 - At)/(EpsEqtHist*EpsEqtHist);
    double aux2 = At*Bt/(exp(Bt*(EpsEqtHist - epsd0)));
    C[0][0] = (1.0 - DtIter)*E - (aux1 + aux2)*E*EpsEqtHist;
  }
  else // DAMAGE_COMPRESSION
  {
    double aux1 = epsd0*(1.0 - Ac)/(EpsEqcHist*EpsEqcHist);
    double aux2 = Ac*Bc/(exp(Bc*(EpsEqcHist - epsd0)));
    C[0][0] = (1.0 + DcIter)*E - (aux1 + aux2)*E*EpsEqcHist;  //VERIFICAFR SINAL
  }
  Ct= C[0][0];
}

// =============================== UpdateState =============================

void cMod1DDamageMazars :: UpdateState(void)
{
  // Update the model parameters.

  Dt = DtIter;
  Dc = DcIter;
  EpsEqtHist = EpsEqtHistIter;
  EpsEqcHist = EpsEqcHistIter;
  cout << "Ct = " << Ct <<  "\n";
}


// -------------------------------------------------------------------------
// Class cMod1DDamageMuModel:
// -------------------------------------------------------------------------

// ========================= cMod1DDamageMuModel =========================

cMod1DDamageMuModel :: cMod1DDamageMuModel(cMaterial *mat) :
                       cConstModel(mat)
{
  // Initialize the internal variables.

  Status = ELASTIC_TENSION;
  DtIter = Dt = 0.0;
  DcIter = Dc = 0.0;

  int np = Mat->NumParam( );
  double *param = new double[np];
  Mat->GetParam(param);
  double epsdt0= param[5];
  double epsdc0= param[6];
  delete []param;

  EpsEqtHistIter = EpsEqtHist = epsdt0;
  EpsEqcHistIter = EpsEqcHist = epsdc0;
}

// ======================== ~cMod1DDamageMuModel =========================

cMod1DDamageMuModel :: ~cMod1DDamageMuModel(void)
{
}

// ================================= Stress ================================

int cMod1DDamageMuModel :: Stress(cVector &eps, cVector &sig)
{
  // Get material parameters.

  int np = Mat->NumParam( );
  double *param = new double[np];
  Mat->GetParam(param);
  double E    = param[0];
  double Ac   = param[1];
  double Bc   = param[2];
  double At   = param[3];
  double Bt   = param[4];
  double epsdt0= param[5];
  double epsdc0= param[6];
  double Nu   = param[7];
  delete []param;

  // Set the trial values of the internal variables.

  DtIter = Dt;
  DcIter = Dc;
  EpsEqtHistIter = EpsEqtHist;
  EpsEqcHistIter = EpsEqcHist;
  //EpsEqtHistIter = MAX(EpsEqtHist, epsd0); // ESTRANHO
  //EpsEqcHistIter = MAX(EpsEqcHist, epsd0); // ESTRANHO

  // Determine the loading status and update the damage parameters.

  if (eps[0] >= 0.0) // Tension
  {
    double epseqt = eps[0];
    if (eps[0] <= EpsEqtHistIter)  // Elastic loading/unloading
    {
      Status = ELASTIC_TENSION;
    }
    else                          // Damage evolution
    {
      Status = DAMAGE_TENSION;
      DtIter = 1.0 - epsdt0*(1.0 - At)/epseqt - At/(exp(Bt*(epseqt - epsdt0)));
      EpsEqtHistIter = epseqt;
      if (DtIter < 0.0)
      {
        DtIter = 0.0;
        EpsEqtHistIter = EpsEqtHist;
      }
      if (DtIter > 1.0)
      {
        DtIter = 1.0;
        EpsEqtHistIter = EpsEqtHist;
      }

      //if (DtIter < 0.0 || DtIter > 1.0)  // ESTRANHO
      //{
      //  DtIter = Dt;
      //  EpsEqtHistIter = EpsEqtHist;
      //}
    }
  }
  else               // Compression
  {
    double epseqc = fabs(eps[0]);
    if (epseqc <= EpsEqcHistIter)  // Elastic loading/unloading
    {
      Status = ELASTIC_COMPRESSION;
    }
    else                          // Damage evolution
    {
      Status = DAMAGE_COMPRESSION;
      DcIter = 1.0 - epsdc0*(1.0 - Ac)/epseqc - Ac/(exp(Bc*(epseqc - epsdc0)));
      EpsEqcHistIter = epseqc;
      if (DcIter < 0.0)
      {
        DcIter = 0.0;
        EpsEqcHistIter = EpsEqcHist;
      }
      if (DcIter > 1.0)
      {
        DcIter = 1.0;
        EpsEqcHistIter = EpsEqcHist;
      }

      //if (DcIter < 0.0 || DcIter > 1.0)  // ESTRANHO
      //{
      //  DcIter = Dc;
      //  EpsEqcHistIter = EpsEqcHist;
      //}
    }
  }

  // Evaluate the current stress using the secant modulus.

  if (eps[0] >= 0.0)
  {
    sig[0] = (1.0 - DtIter)*E*eps[0];
  }
  else
  {
    sig[0] = (1.0 - DcIter)*E*eps[0];
  }

  return(1);
}

// ================================= TangMat ===============================

void cMod1DDamageMuModel:: TangMat(cMatrix &C)
{
  // Get material parameters.

  int np = Mat->NumParam( );
  double *param = new double[np];
  Mat->GetParam(param);
  double E    = param[0];
  double Ac   = param[1];
  double Bc   = param[2];
  double At   = param[3];
  double Bt   = param[4];
  double epsdt0= param[5];
  double epsdc0= param[6];
  double Nu   = param[7];
  delete []param;

  // Compute the tangent modulus.

  C.Zero( );
  if (Status == ELASTIC_TENSION)
  {
    C[0][0] = (1.0 - DtIter)*E;
  }
  else if (Status == ELASTIC_COMPRESSION)
  {
    C[0][0] = (1.0 - DcIter)*E;
  }
  else if (Status == DAMAGE_TENSION)
  {
    double aux1 = epsdt0*(1.0 - At)/(EpsEqtHist*EpsEqtHist);
    double aux2 = At*Bt/(exp(Bt*(EpsEqtHist - epsdt0)));
    C[0][0] = (1.0 - DtIter)*E - (aux1 + aux2)*E*EpsEqtHist;
  }
  else // DAMAGE_COMPRESSION
  {
    double aux1 = epsdc0*(1.0 - Ac)/(EpsEqcHist*EpsEqcHist);
    double aux2 = Ac*Bc/(exp(Bc*(EpsEqcHist - epsdc0)));
    C[0][0] = (1.0 - DcIter)*E - (aux1 + aux2)*E*EpsEqcHist;  //VERIFICAFR SINAL
  }
  Ct= C[0][0];
}

// =============================== UpdateState =============================

void cMod1DDamageMuModel :: UpdateState(void)
{
  // Update the model parameters.

  Dt = DtIter;
  Dc = DcIter;
  EpsEqtHist = EpsEqtHistIter;
  EpsEqcHist = EpsEqcHistIter;
  //cout << "Ct = " << Ct <<  "\n";

    static int d = cControl :: GetTotFactor( );
  //static ofstream outDT;
  //static ofstream outDC;

  static fstream outDT("dt.txt",std::fstream::out);
  static fstream outDC("dc.txt",std::fstream::out);
  static fstream outCT("ct.txt",std::fstream::out);
  int tot = cControl :: GetTotFactor( );

  if (tot != d)
  {
     d = tot;

     outDT << "\n";
     outDC << "\n";
     outCT << "\n";
  }


  outDT << Dt << " ";
  outDC << Dc << " ";
  outCT << Ct << " ";
}


// -------------------------------------------------------------------------
// Class cMod1DLeeFenves:
// -------------------------------------------------------------------------

// ============================= cMod1DLeeFenves ===========================

cMod1DLeeFenves :: cMod1DLeeFenves(cMaterial *mat) :
                   cConstModel(mat)
{
  // Initialize the internal variables.

  Status = ELASTIC;
  Dt = DtIter = 0.0;
  Dc = DcIter = 0.0;
  KapatIter = Kapat = 0.0;
  KapacIter = Kapac = 0.0;
  Lambda = LambdaIter = 0.0;
  SigEf = SigEfIter  = 0.0;
  PlasticEpsIter = PlasticEps = 0.0;
  Ct = 0.0;
}

// ============================ ~cMod1DLeeFenves ===========================

cMod1DLeeFenves :: ~cMod1DLeeFenves(void)
{
}

// ================================= Stress ================================

int cMod1DLeeFenves :: Stress(cVector &eps, cVector &sig)
{
  // Get material parameters.

  int np = Mat->NumParam( );
  double *param = new double[np];
  Mat->GetParam(param);
  double E    = param[0];
  double fb0  = param[1];
  double fc0  = param[2];
  double ft0  = param[3];
  double a_c  = param[4];
  double b_c  = param[5];
  double d_c  = param[6];
  double a_t  = param[7];
  double b_t  = param[8];
  double d_t  = param[9];
  double Gt   = param[10];
  double Gc   = param[11];
  double lt   = param[12];
  double lc   = param[13];
  double s0   = param[14];
  delete []param;

  // Elastic predictor.

  DtIter = Dt;
  DcIter = Dc;
  KapatIter = Kapat;
  KapacIter = Kapac;
  LambdaIter = Lambda;
  PlasticEpsIter = PlasticEps;
  double sigefpred = E*(eps[0] - PlasticEpsIter);

  double alpha = (fb0/fc0-1.0)/(2.0*fb0/fc0-1.0);
  double phic  = 1.0 + a_c*(2.0+a_c)*KapacIter;
  double cc    = -fc0*pow(1/a_c*(1.0+a_c-sqrt(phic)),1.0-d_c/b_c)*sqrt(phic);
  double phit  = 1.0 + a_t*(2+a_t)*KapatIter;
  double ct    = ft0*pow(1.0/a_t*(1.0+a_t-sqrt(phit)),1.0-d_t/b_t)*sqrt(phit);
  double beta  = cc/ct*(1.0-alpha)-(1.0+alpha);
  double r = 0.0;
  double fpred = 0.0;

  if (sigefpred >= 0.0)  // Tension
  {
    fpred = 1.0/(1.0-alpha)*(alpha*sigefpred + sigefpred + beta*sigefpred) - cc;
  }
  else                   // Compression
  {
    fpred = 1.0/(1.0-alpha)*(alpha*sigefpred + fabs(sigefpred)) - cc;
  }

  if (fpred <= 0.0) // Elastic increment
  {
    Status = ELASTIC;
    if (sigefpred >= 0.0)
    {
      r = 1.0;
    }
    else
    {
      r = 0.0;
    }
    double s = s0 + (1.0-s0)*r;  // Cracking closing
    double sigef = sigefpred;
    sig[0] = sigef*(1.0-DcIter)*(1.0-s*DtIter);
  }

  else              // Plastic increment
  {
    Status = PLASTIC;

    // Local N-.0R algorithm.

    double kkc = Kapac;
    double kkt = Kapat;
    double kc  = KapacIter;
    double kt  = KapatIter;
    double phic  = 1.0 + a_c*(2.0+a_c)*kc;
    double fcr = fc0/a_c*((1.0+a_c)*sqrt(phic)-phic);
    double fc  = fc0*pow(1.0/a_c*(1.0+a_c-sqrt(phic)),1.0-d_c/b_c)*sqrt(phic);
    double cc  = -fc;
    double phit = 1.0 + a_t*(2.0+a_t)*kt;
    double ftr = ft0/a_t*((1.0+a_t)*sqrt(phit)-phit);
    double ft  = ft0*pow(1.0/a_t*(1+a_t-sqrt(phit)),1.0-d_t/b_t)*sqrt(phit);
    double ct  = ft;
    double lambda = 0.0;
    double sigef = sigefpred;
    double psn1 = PlasticEps;
    double psn  = PlasticEpsIter;

    //cout << "phic = " << phic << " fcr = " << fcr << " fc = " << fc << " cc = " << cc <<"\n";
    //cout << "phit = " << phit << " ftr = " << ftr << " ft = " << ft << " ct = " << ct <<"\n";

    double gc = Gc/lc;
    double gt = Gt/lt;
    double ap = 0.2;

    int NITER = 30;            // Maximum number of steps to plastic N-R
    double TOL1 = 1e-6;        // Tolerance 1
    double TOL2 = 1e-6;        // Tolerance 2

    if (sigefpred >= 0.0)      // Tension
    {
      double r1 = -psn + psn1 + lambda*(sqrt(2.0/3.0) + ap);
      double r2 = -kt + kkt + ftr/gt*lambda*(sqrt(2.0/3.0) + ap);

      double f1 = 1.0/(1.0-alpha)*(alpha*sigef + sigef + beta*sigef) - cc;

      for (int i = 0; i < NITER; i++)
      {
        cVector z(2);  //auxiliary vector
        z.Zero( );

        cVector r(2);
        r.Zero( );
        r(0) = r1;
        r(1) = r2;

        cMatrix iA(2, 2);
        iA.Zero();
        iA[0][0] = 1.0/E;
        iA[1][1] = -1.0 + ft0*(2.0+a_t)/gt*(1.0/2.0*(1.0+a_t)/phit - 1.0)*lambda*(sqrt(2.0/3.0) + ap);

        cMatrix A(2, 2);
        A.Zero();

        double det = iA[0][0]*iA[1][1] - iA[0][1]*iA[1][0];
        if (fabs(det) == 0.0) return(0);

        A[0][0] =  iA[1][1]/(det);
        A[0][1] = -iA[0][1]/(det);
        A[1][0] = -iA[1][0]/(det);
        A[1][1] =  iA[0][0]/(det);

        #if 0
          cout << " iA = " << scientific << setprecision(6);
          iA.Print( );
        #endif

        cVector b(2);
        b.Zero( );
        b(0) = sqrt(2.0/3.0)+ap;
        b(1) = ftr/gt*(sqrt(2.0/3.0)+ap);

        //cout << "ftr = " << ftr<<"\n";
        //cout << "gt = " << gt<<"\n";

        #if 0
          cout << " b = " << scientific << setprecision(6);
          b.Print( );
        #endif

        cVector c(2);
        c.Zero( );
        c(0) = cc/ct;
        c(1) = a_t*(2.0+a_t)*cc*sqrt(phit)*sigef/(2.0*ct*phit)*((1.0-d_t/b_t)/(1.0+a_t-sqrt(phit))-1.0/sqrt(phit));

        #if 0
          cout << " c = " << scientific << setprecision(6);
          c.Print( );
        #endif

        double dlambda = f1;  //Calcular incremento de lambda dlambda = (f1 - c*A*r)/(c*A*b)
        z = A*r;

        dlambda -= c*z;
        z = A*b;

        dlambda *= 1.0/(c*z);

        //cout << scientific << setprecision(6);
        //cout << "dlambda = " << dlambda <<"\n";

        cVector d(2);    //Calcular d = -A*R -A*B*dlambda;
        d.Zero( );
        d = A*r;
        d *= -1.0;
        z = A*b;
        z *= dlambda;
        z *= -1.0;
        d += z;

        #if 0
          cout << " d = " << scientific << setprecision(6);
          d.Print( );
        #endif

        double dpsn = -d(0)/E;
        psn += dpsn;
        sigef += d(0);
        kt += d(1);

        #if 0
          cout << "kt = " << scientific << setprecision(6) << kt <<"\n";
        #endif

        lambda += dlambda;

        phit = 1.0 + a_t*(2.0+a_t)*kt;
        ftr = ft0/a_t*((1.0+a_t)*sqrt(phit)-phit);
        ft  = ft0*pow(1.0/a_t*(1.0+a_t-sqrt(phit)),1.0-d_t/b_t)*sqrt(phit);
        ct  = ft;
        beta = cc/ct*(1.0-alpha)-(1.0+alpha);

        r1 = -psn + psn1 + lambda*(sqrt(2.0/3.0) + ap);
        r2 = -kt + kkt + ftr/gt*lambda*(sqrt(2.0/3.0) + ap);
        f1 = 1.0/(1.0-alpha)*(alpha*sigef + sigef + beta*sigef) - cc;

        #if 0
        cout << scientific << setprecision(6);
        cout << "i = " << i+1 << " r1 = " << r.Length( ) << " f1 = " << f1 << "\n";
        cout << resetiosflags(ios::scientific) << setprecision(6);
        #endif

        if(r.Length( ) < TOL1 && fabs(f1) < TOL2)
        {
           //cout << "Convergence achieved - Lee and Fenves Model!!!\n";
           //cout << "\n";
           break;
        }
        //#if 0
        //else if (i+1 == NITER)
        //{
        //  if(r.Length( ) >= TOL1 || fabs(f1) >= TOL2)
        //  {
        //    cout << "Convergence NOT achieved - Lee and Fenves Model!!!\n";
            //exit(0);
        //  }
        //}
        //#endif
      }
    }
    else  // Compression
    {
      double r1 = -psn + psn1 + lambda*(-sqrt(2.0/3.0) + ap);
      double r2 = -kc + kkc + fcr/gc*lambda*(-sqrt(2.0/3.0) + ap);

      double f1 = 1.0/(1.0-alpha)*(alpha*sigef + fabs(sigef)) - cc;

      for (int i = 0; i < NITER; i++)
      {
        cVector z(2);  //auxiliary vector
        z.Zero( );

        cVector r(2);
        r.Zero( );
        r(0) = r1;
        r(1) = r2;

        cMatrix iA(2, 2);
        iA.Zero();
        iA[0][0] = 1.0/E;
        iA[1][1] = -1.0 + fc0*(2.0+a_c)/gc*(1.0/2.0*(1.0+a_c)/phic - 1.0)*lambda*(-sqrt(2.0/3.0) + ap);

        cMatrix A(2, 2);
        A.Zero();

        double det = iA[0][0]*iA[1][1] - iA[0][1]*iA[1][0];
        if (fabs(det) == 0.0) return(0);

        A[0][0] =  iA[1][1]/(det);
        A[0][1] = -iA[0][1]/(det);
        A[1][0] = -iA[1][0]/(det);
        A[1][1] =  iA[0][0]/(det);

        #if 0
          cout << " iA = " << scientific << setprecision(6);
          iA.Print( );
        #endif

        cVector b(2);
        b.Zero( );
        b(0) = -sqrt(2.0/3.0)+ap;
        b(1) = fcr/gc*(-sqrt(2.0/3.0)+ap);

        #if 0
          cout << " b = " << scientific << setprecision(6);
          b.Print( );
        #endif

        cVector c(2);
        c.Zero( );
        c(0) = -1.0;
        c(1) = a_c*(2.0+a_c)*cc/2.0*((1.0-d_c/b_c)/(sqrt(phic)*(1.0+a_c-sqrt(phic))) -1.0/phic);

        #if 0
          cout << " c = " << scientific << setprecision(6);
          c.Print( );
        #endif

        double dlambda = f1;  //Calcular incremento de lambda dlambda = (f1 - c*A*r)/(c*A*b)
        z = A*r;

        dlambda -= c*z;
        z = A*b;

        dlambda *= 1.0/(c*z);

        //cout << scientific << setprecision(6);
        //cout << "dlambda = " << dlambda <<"\n";

        cVector d(2);    //Calcular d = -A*R -A*B*dlambda;
        d.Zero( );
        d = A*r;
        d *= -1.0;
        z = A*b;
        z *= dlambda;
        z *= -1.0;
        d += z;

        #if 0
          cout << " d = " << scientific << setprecision(6);
          d.Print( );
        #endif

        double dpsn = -d(0)/E;
        psn += dpsn;
        sigef += d(0);
        kc += d(1);

        lambda += dlambda;

        phic = 1.0 + a_c*(2.0+a_c)*kc;
        fcr = fc0/a_c*((1.0+a_c)*sqrt(phic)-phic);
        fc  = fc0*pow(1.0/a_c*(1.0+a_c-sqrt(phic)),1.0-d_c/b_c)*sqrt(phic);
        cc  = -fc;

        r1 = -psn + psn1 + lambda*(-sqrt(2.0/3.0) + ap);
        r2 = -kc + kkc + fcr/gc*lambda*(-sqrt(2.0/3.0) + ap);
        f1 = 1.0/(1.0-alpha)*(alpha*sigef + fabs(sigef)) - cc;

        #if 0
        cout << scientific << setprecision(6);
        cout << "i = " << i+1 << " r1 = " << r.Length( ) << " f1 = " << f1 << "\n";
        cout << resetiosflags(ios::scientific) << setprecision(6);
        #endif

        if(r.Length( ) < TOL1 && fabs(f1) < TOL2)
        {
           //cout << "Convergence achieved - Lee and Fenves Model!!!\n";
           //cout << "\n";
           break;
        }
        //#if 0
        else if (i+1 == NITER)
        {
          if(r.Length( ) >= TOL1 || fabs(f1) >= TOL2)
          {
            cout << "Convergence NOT achieved - Lee and Fenves Model!!!\n";
            //exit(0);
          }
        }
        //#endif
      }
    }

    KapatIter = kt;
    KapacIter = kc;
    LambdaIter = lambda;
    PlasticEpsIter = psn;
    DtIter = 1.0 - pow(1.0/a_t*(1.0+a_t-sqrt(1.0+a_t*(2.0+a_t)*KapatIter)),d_t/b_t);
    DcIter = 1.0 - pow(1.0/a_c*(1.0+a_c-sqrt(1.0+a_c*(2.0+a_c)*KapacIter)),d_c/b_c);
    SigEfIter = sigef;

    if (SigEfIter >= 0.0)
    {
      r = 1.0;
    }
    else
    {
      r = 0.0;
    }
    double s = s0 + (1.0-s0)*r;  // Cracking closing
    sig[0] = SigEfIter*(1.0-DcIter)*(1.0-s*DtIter);
  }
return (1);
}

// ================================= TangMat ===============================

void cMod1DLeeFenves :: TangMat(cMatrix &C)
{
  // Get material parameters.

  int np = Mat->NumParam( );
  double *param = new double[np];
  Mat->GetParam(param);
  double E    = param[0];
  double fb0  = param[1];
  double fc0  = param[2];
  double ft0  = param[3];
  double a_c  = param[4];
  double b_c  = param[5];
  double d_c  = param[6];
  double a_t  = param[7];
  double b_t  = param[8];
  double d_t  = param[9];
  double Gt   = param[10];
  double Gc   = param[11];
  double lt   = param[12];
  double lc   = param[13];
  double s0   = param[14];
  delete []param;

  double gc = Gc/lc;
  double gt = Gt/lt;
  double ap = 0.2;
  double r = 0.0;

  if (SigEfIter >= 0.0)
  {
    r = 1.0;
  }
  else
  {
    r = 0.0;
  }
  double s = s0 + (1.0-s0)*r;  // Cracking closing

  // Compute the tangent modulus.
  C.Zero( );
  if (Status == ELASTIC)
  {
    C[0][0] = E*(1-DcIter)*(1-s*DtIter);
    Ct= C[0][0];

    //cout << "Ct = " << Ct <<  "\n";
  }
  else
  {
    if (SigEfIter >= 0)     // Tension
    {
      double phit  = 1.0 + a_t*(2.0+a_t)*KapatIter;
      double ftr   = ft0/a_t*((1.0+a_t)*sqrt(phit)-phit);
      double ft    = ft0*pow(1.0/a_t*(1.0+a_t-sqrt(phit)),1.0-d_t/b_t)*sqrt(phit);
      double ct    = ft;
      double phic  = 1.0 + a_c*(2.0+a_c)*KapacIter;
      double fcr   = fc0/a_c*((1.0+a_c)*sqrt(phic)-phic);
      double fc    = fc0*pow(1.0/a_c*(1.0+a_c-sqrt(phic)),1.0-d_c/b_c)*sqrt(phic);
      double cc    = -fc;
      double TDsig = 0.0;
      double TDlbd = -((1.0-DcIter)*(1-s*DtIter)*(2.0+a_t)*d_t*a_t*ftr*(sqrt(2.0/3.0)+ap))/(2.0*b_t*gt*sqrt(phit)*(1.0+a_t-sqrt(phit))*(-1.0+(LambdaIter*ft0*(2.0+a_t)*(sqrt(2.0/3.0)+ap))/gt*((1.0+a_t)/(2.0*sqrt(phit))-1.0)));
      double Tfsig = cc/ct;
      double Tflbd = -1.0/(a_t*gt*(-1.0+(LambdaIter*ft0*(2.0+a_t)*(sqrt(2.0/3.0)+ap))/gt*((1.0+a_t)/(2.0*phit)-1.0)))*(a_t*(2.0+a_t)*cc*sqrt(phit)/(2.0*ct*phit)*((1.0-d_t/b_t)/(1+a_t-sqrt(phit))-1.0/sqrt(phit))*SigEfIter*ft0*((1.0+a_t)*sqrt(phit)-phit)*(sqrt(2.0/3.0)+ap));
      //double Tflbd = -1/(at*gt*(-1+(lambda*ft0*(2+at)*(sqrt(2/3)+ap))/gt*((1+at)/(2*phit)-1)))*(at*(2+at)*cc*sqrt(phit)/(2*ct*phit)*((1-dt/bt)/(1+at-sqrt(phit))-1/sqrt(phit))*sigef*ft0*((1+at)*sqrt(phit)-phit)*(sqrt(2/3)+ap));
      double DPHI  = sqrt(2.0/3.0) + ap;
      double S = E;

      #if 0
        cout << scientific << setprecision(6);
        cout << "Kapat = " << Kapat << " Kapac = " << Kapac << "\n";
        cout << "phit = " << phit << " ftr = " << ftr << " ft = " << ft << " ct = " << ct << "\n";
        cout << "phic = " << phic << " fcr = " << fcr << " fc = " << fc << " cc = " << cc << "\n";
        cout << "TDsig = " << TDsig << " TDlbd = " << TDlbd << " Tfsig = " << Tfsig << " Tflbd = " << Tflbd << " DPHI = " << DPHI << "\n";
        cout << resetiosflags(ios::scientific) << setprecision(6);
      #endif

      C[0][0] = ((1.0-DcIter)*(1.0-s*DtIter) - SigEfIter*TDsig);
      C[0][0] *= (S - (S*DPHI*Tfsig*S)/(Tfsig*S*DPHI-Tflbd));
      C[0][0] += -TDlbd*((SigEfIter*Tfsig*S)/(Tfsig*S*DPHI-Tflbd));
      Ct= C[0][0];

    }
    else  // Compression
    {
      //#if 0
      double phit  = 1.0 + a_t*(2.0+a_t)*KapatIter;
      double ftr   = ft0/a_t*((1.0+a_t)*sqrt(phit)-phit);
      double ft    = ft0*pow(1.0/a_t*(1.0+a_t-sqrt(phit)),1.0-d_t/b_t)*sqrt(phit);
      double ct    = ft;
      double phic  = 1.0 + a_c*(2.0+a_c)*KapacIter;
      double fcr   = fc0/a_c*((1.0+a_c)*sqrt(phic)-phic);
      double fc    = fc0*pow(1.0/a_c*(1.0+a_c-sqrt(phic)),1.0-d_c/b_c)*sqrt(phic);
      double cc    = -fc;
      double TDsig = 0.0;
      double TDlbd = -((1.0-DcIter)*(1.0-s*DtIter)*(2.0+a_c)*d_c*a_c*fcr*(-sqrt(2.0/3.0)+ap))/(2.0*b_c*gc*sqrt(phic)*(1.0+a_c-sqrt(phic))*(-1.0+(LambdaIter*fc0*(2.0+a_c)*(-sqrt(2.0/3.0)+ap))/gc*((1.0+a_c)/(2.0*sqrt(phic))-1.0)));
      double Tfsig = -1.0;
      double Tflbd = -1.0/(a_c*gc*(-1.0+(LambdaIter*fc0*(2.0+a_c)*(-sqrt(2.0/3.0)+ap))/gc*((1.0+a_c)/(2.0*phic)-1.0)))*(a_c*(2.0+a_c)*cc/(2.0*sqrt(phic))*((1.0-d_c/b_c)/(1.0+a_c-sqrt(phic))-1.0/sqrt(phic))*fc0*((1.0+a_c)*sqrt(phic)-phic)*(-sqrt(2.0/3.0)+ap));
      double DPHI  = -sqrt(2.0/3.0) + ap;
      double S = E;
      //#endif

      #if 0
      double phit  = 1.0 + a_t*(2.0+a_t)*Kapat;
      double ftr   = ft0/a_t*((1.0+a_t)*sqrt(phit)-phit);
      double ft    = ft0*pow(1.0/a_t*(1.0+a_t-sqrt(phit)),1.0-d_t/b_t)*sqrt(phit);
      double ct    = ft;
      double phic  = 1.0 + a_c*(2.0+a_c)*Kapac;
      double fcr   = fc0/a_c*((1.0+a_c)*sqrt(phic)-phic);
      double fc    = fc0*pow(1.0/a_c*(1.0+a_c-sqrt(phic)),1.0-d_c/b_c)*sqrt(phic);
      double cc    = -fc;
      double TDsig = 0.0;
      double TDlbd = -((1.0-Dc)*(1.0-s*Dt)*(2.0+a_c)*d_c*a_c*fcr*(-sqrt(2.0/3.0)+ap))/(2.0*b_c*gc*sqrt(phic)*(1.0+a_c-sqrt(phic))*(-1.0+(Lambda*fc0*(2.0+a_c)*(-sqrt(2.0/3.0)+ap))/gc*((1.0+a_c)/(2.0*sqrt(phic))-1.0)));
      double Tfsig = -1.0;
      double Tflbd = -1.0/(a_c*gc*(-1.0+(Lambda*fc0*(2.0+a_c)*(-sqrt(2.0/3.0)+ap))/gc*((1.0+a_c)/(2.0*phic)-1.0)))*(a_c*(2.0+a_c)*cc/(2.0*sqrt(phic))*((1.0-d_c/b_c)/(1.0+a_c-sqrt(phic))-1.0/sqrt(phic))*fc0*((1.0+a_c)*sqrt(phic)-phic)*(-sqrt(2.0/3.0)+ap));
      double DPHI  = -sqrt(2.0/3.0) + ap;
      double S = E;
      #endif

      #if 0
        cout << scientific << setprecision(6);
        cout << "Kapat = " << Kapat << " Kapac = " << Kapac << "\n";
        cout << "phit = " << phit << " ftr = " << ftr << " ft = " << ft << " ct = " << ct << "\n";
        cout << "phic = " << phic << " fcr = " << fcr << " fc = " << fc << " cc = " << cc << "\n";
        cout << "TDsig = " << TDsig << " TDlbd = " << TDlbd << " Tfsig = " << Tfsig << " Tflbd = " << Tflbd << " DPHI = " << DPHI << "\n";
        cout << resetiosflags(ios::scientific) << setprecision(6);
      #endif

      C[0][0] = ((1.0-DcIter)*(1.0-s*DtIter) - SigEfIter*TDsig);
      C[0][0] *= (S - (S*DPHI*Tfsig*S)/(Tfsig*S*DPHI-Tflbd));
      C[0][0] += -TDlbd*((SigEfIter*Tfsig*S)/(Tfsig*S*DPHI-Tflbd));
      Ct= C[0][0];
    }
  }
}

// =============================== UpdateState =============================

void cMod1DLeeFenves :: UpdateState(void)
{
  // Update the model parameters.

  Dt = DtIter;
  Dc = DcIter;
  Kapat = KapatIter;
  Kapac = KapacIter;
  Lambda = LambdaIter;
  PlasticEps = PlasticEpsIter;
  SigEf = SigEfIter;

  //cout << "Ct = " << Ct <<  "\n";
  //cout << "Dt = " << Dt <<  "\n";
  //cout << "Pls = " << PlasticEps <<  "\n";

   static int d = cControl :: GetTotFactor( );
  //static ofstream outDT;
  //static ofstream outDC;
  //static ofstream outPS;

  static fstream outDT("dt.txt",std::fstream::out);
  static fstream outDC("dc.txt",std::fstream::out);
  static fstream outCT("ct.txt",std::fstream::out);
  static fstream outPS("pls.txt",std::fstream::out);
  int tot = cControl :: GetTotFactor( );

  if (tot != d)
  {
     d = tot;

     outDT << "\n";
     outDC << "\n";
     outCT << "\n";
     outPS << "\n";
  }


  outDT << Dt << " ";
  outDC << Dc << " ";
  outCT << Ct << " ";
  outPS << PlasticEps  << " ";

  //outDT.close();
  //outDC.close();

  //cout << "Ct = " << Ct <<  "\n";
  //cout << "Ct = "  << scientific << setprecision(4)<< Ct <<  "\n";
}

// ======================================================= End of file =====
