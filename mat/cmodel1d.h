// ------------------------------------------------------------------------
// cmodel1d.h - file containing the definition of uniaxial constitutive
//              models.
// ------------------------------------------------------------------------

#ifndef _CMODEL1D_H
#define _CMODEL1D_H

#include "cmodel.h"

// ------------------------------------------------------------------------
// Definition of the Linear Constitutive Model class:
//
class cMod1DLinear : public cConstModel
{
 public:
                cMod1DLinear(cMaterial *);
  virtual      ~cMod1DLinear(void);
          int   Stress(cVector &, cVector &);
          void  TangMat(cMatrix &);
          void  UpdateState(void) { }
};


// ------------------------------------------------------------------------
// Definition of the Nonlinear Elastic Constitutive Model class:
//
class cMod1DElastNonlin : public cConstModel
{
 protected:
  double  Eps;  // Current strain

 public:
                cMod1DElastNonlin(cMaterial *);
  virtual      ~cMod1DElastNonlin(void);
          int   NumStrainLim(void);
          void  GetStrainLim(cVector &);
          int   GetCrvDegree(int) { return 1; }
          int   Stress(cVector &, cVector &);
          void  TangMat(cMatrix &);
          void  UpdateState(void) { }

};


// ------------------------------------------------------------------------
// Definition of the Nonlinear Elastic EuroCEB Constitutive Model class:
//
class cMod1DConcEuroCEB : public cConstModel
{
 protected:
  double  Eps;  // Current strain

 public:
                cMod1DConcEuroCEB(cMaterial *);
  virtual      ~cMod1DConcEuroCEB(void);
          int   NumStrainLim(void);
          void  GetStrainLim(cVector &);
          int   GetCrvDegree(int);
          int   Stress(cVector &, cVector &);
          void  TangMat(cMatrix &);
          void  UpdateState(void) { }
};


// ------------------------------------------------------------------------
// Definition of the Nonlinear Elastic NBRCEB Constitutive Model class:
//
class cMod1DConcNBRCEB : public cConstModel
{
 protected:
  double  Eps;  // Current strain

 public:
                cMod1DConcNBRCEB(cMaterial *);
  virtual      ~cMod1DConcNBRCEB(void);
          int   NumStrainLim(void);
          void  GetStrainLim(cVector &);
          int   GetCrvDegree(int);
          int   Stress(cVector &, cVector &);
          void  TangMat(cMatrix &);
          void  UpdateState(void) { }
};


// ------------------------------------------------------------------------
// Definition of the Viscoelastic Constitutive Model class:
//
class cMod1DViscoElastic : public cConstModel
{
 protected:
  double  Time;
  double  Sig,SigIter;  // Converged and interactive stresses
  double  Eps,EpsIter;  // Converged and interactive strains
  double  EpsRate;      // Strain rate at the previous step
  double *S;            // Hereditary vector

 public:
                cMod1DViscoElastic(cMaterial *);
  virtual      ~cMod1DViscoElastic(void);
          int   Stress(cVector &, cVector &);
          void  TangMat(cMatrix &);
          void  UpdateState(void);
};


// ------------------------------------------------------------------------
// Definition of the Elastoplastic Constitutive Model with Linear Isotropic
// and Kinematic Hardening class:
//
class cMod1DPlastLinHard : public cConstModel
{
 typedef enum
 {
  ELASTIC,
  PLASTIC
 } eStatus;

 protected:
  eStatus Status;
  double  PlasticEps,PlasticEpsIter;
  double  Alpha,AlphaIter;
  double  BackStress,BackStressIter;

 public:
                cMod1DPlastLinHard(cMaterial *);
  virtual      ~cMod1DPlastLinHard(void);
          int   Stress(cVector &, cVector &);
          void  TangMat(cMatrix &);
          void  UpdateState(void);
};


// ------------------------------------------------------------------------
// Definition of the Elastoplastic Constitutive Model with Nonlinear
// Isotropic Hardening class:
//
class cMod1DPlastIsoNonlinHard : public cConstModel
{
 typedef enum
 {
  ELASTIC,
  PLASTIC
 } eStatus;

 protected:
  eStatus Status;
  double  PlasticEps,PlasticEpsIter;
  double  Alpha,AlphaIter;

 public:
                cMod1DPlastIsoNonlinHard(cMaterial *);
  virtual      ~cMod1DPlastIsoNonlinHard(void);
          int   Stress(cVector &, cVector &);
          void  TangMat(cMatrix &);
          void  UpdateState(void);
};


// ------------------------------------------------------------------------
// Definition of the Elastoplastic Constitutive Model with Exponential
// Isotropic Hardening class:
//
class cMod1DPlastIsoExpHard : public cConstModel
{
 typedef enum
 {
  ELASTIC,
  PLASTIC
 } eStatus;

 protected:
  eStatus Status;
  double  PlasticEps,PlasticEpsIter;
  double  Alpha,AlphaIter;

 public:
                cMod1DPlastIsoExpHard(cMaterial *);
  virtual      ~cMod1DPlastIsoExpHard(void);
          int   Stress(cVector &, cVector &);
          void  TangMat(cMatrix &);
          void  UpdateState(void);
};


// ------------------------------------------------------------------------
// Definition of the Damage Mazars Model class:
//
class cMod1DDamageMazars: public cConstModel
{
 typedef enum
 {
  ELASTIC_TENSION,
  ELASTIC_COMPRESSION,
  DAMAGE_TENSION,
  DAMAGE_COMPRESSION,
 } eStatus;

 protected:
  eStatus Status;
  double  Dt, DtIter;
  double  Dc, DcIter;
  double  EpsEqtHist, EpsEqtHistIter;
  double  EpsEqcHist, EpsEqcHistIter;
  double  Ct;   // Current tangent modulus

 public:
                cMod1DDamageMazars(cMaterial *);
  virtual      ~cMod1DDamageMazars(void);
          int   Stress(cVector &, cVector &);
          void  TangMat(cMatrix &);
          void  UpdateState(void);
};


// ------------------------------------------------------------------------
// Definition of the Damage Mazars Mu-Model class:
//
class cMod1DDamageMuModel: public cConstModel
{
 typedef enum
 {
  ELASTIC_TENSION,
  ELASTIC_COMPRESSION,
  DAMAGE_TENSION,
  DAMAGE_COMPRESSION,
 } eStatus;

 protected:
  eStatus Status;
  double  Dt, DtIter;
  double  Dc, DcIter;
  double  EpsEqtHist, EpsEqtHistIter;
  double  EpsEqcHist, EpsEqcHistIter;
  double  Ct;   // Current tangent modulus

 public:
                cMod1DDamageMuModel(cMaterial *);
  virtual      ~cMod1DDamageMuModel(void);
          int   Stress(cVector &, cVector &);
          void  TangMat(cMatrix &);
          void  UpdateState(void);
};


// ------------------------------------------------------------------------
// Definition of the Plastic Damage Lee-Fenves Model class:
//
class cMod1DLeeFenves: public cConstModel
{
 typedef enum
 {
  ELASTIC,
  PLASTIC
 } eStatus;

 protected:
  eStatus Status;
  double  Dt, DtIter;
  double  Dc, DcIter;
  double  Kapat, KapatIter;
  double  Kapac, KapacIter;
  double  PlasticEps,PlasticEpsIter;
  double  Lambda, LambdaIter;
  double  SigEf, SigEfIter;
  double  Ct;   // Current tangent modulus

 public:
                cMod1DLeeFenves(cMaterial *);
  virtual      ~cMod1DLeeFenves(void);
          int   Stress(cVector &, cVector &);
          void  TangMat(cMatrix &);
          void  UpdateState(void);
};

#endif
