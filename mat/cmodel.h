// ------------------------------------------------------------------------
// cmodel.h - file containing the definition of the Constitutive Model
//            class.
// ------------------------------------------------------------------------
// The Constitutive Model class is responsible for the computation of
// current stresses from the given total strains. In order to allow the
// use of use of implicit analysis algorithms (e.g. Newton-Raphson) it
// also computes the tangent constitutive matrix. A set of internal
// variables can be stored to describe the stress history, which is
// different for each integration point. It is important to note the
// material parameters are not stored in this class, but are obtained from
// the Material class when they are required.
//
// ConstModel
// |-- ModLinear
// |   |-- ModOrthoLinear
// |-- ModFGM
// |-- ModViscoElastic
// |-- ModProgFailure
// |   |-- ModProgFailHashin
// |   |-- ModProgFailTsai
// |   |-- ModProgFailEngelstad
// |-- ModPlastIsoLinHard
// |   |-- ModFGMTTO
// |-- ModCZMModeIBilinear
// |-- ModCZMModeIExponential
// |-- Mod1DLinear
// |-- Mod1DElastNonlin
// |-- Mod1DPlastLinHard
// |-- Mod1DPlastIsoNonlinHard
// |-- Mod1DPlastIsoExpHard
// |-- Mod1DConcEuroCEB
// |-- Mod1DConcNBRCEB
// |-- Mod1DViscoElastic
//
// -------------------------------------------------------------------------
// ModProgFailHashin
//
// Failure Criterion: Hashin
//
// Maximum Number of Failures per Point: No limit
//
// References: Sleight, D. W. Progressive failure analysis methodology for
//             laminated composite structures. Technical Report, National
//             Aeronautics and Space Administration (NASA), 1999.
//
//             Pietropaoli, E. Progressive failure analysis of composite
//             structures using a constitutive material model (USERMAT)
//             developed and implemented in ANSYS. Applied Composite
//             Materials, v. 19, 2012.
//
// Modifications: Changed the choice of parameter degradation for each
//                failure mode.
// -------------------------------------------------------------------------
// ModProgFailTsai
//
// Failure Criterion: Tsai-Wu
//
// Maximum Number of Failures per Point: 2
//
// Reference: Kuraishi, A.; Tsai, S. W.; Liu, K. K. A progressive quadratic
//            failure criterion, Part B. Composites Science and Technology,
//            v. 62, 2002.
//
// Modifications: Adapted for a three-dimensional stress state.
//
// -------------------------------------------------------------------------
// ModProgFailEngelstad
//
// Failure Criterion: Tsai-Wu
//
// Maximum Number of Failures per Point: 1
//
// Reference: Engelstad, S. P.; Reddy, J. N.; Knight JR, N. F. Postbuckling
//            response and failure prediction of graphite-epoxy plates
//            loaded in compression. AIAA Journal, v. 30, 1992.
//
// Modifications: Here, only the Tsai-Wu criterion is used. Also, an
//                interactive failure strategy was chosen, with the
//                degradation of multiple properties for each failure mode.
// -------------------------------------------------------------------------
// Static public methods:
// -------------------------------------------------------------------------
//
// static cConstModel *Create(cElement *elm, cMaterial *mat)
//
// This method creates a constitutive model based on the material type
// associated with the parent element and returns a pointer to the created
// model. It also store a pointer to the associated material.
//
// -------------------------------------------------------------------------
// Pure virtual methods:
// -------------------------------------------------------------------------
//
// virtual int Stress(cVector &eps, cVector &sig)
//
//  eps    -  strain vector                                           (in)
//  sig    -  stress vector                                          (out)
//
// This method computes the current stresses from the given total strains.
// It returns 1/0 if the stress were sucessfully computed or not.
// -------------------------------------------------------------------------
// virtual void TangMat(cMatrix &C)
//
//  C      -  tangent constitutive matrix                            (out)
//
// This method computes the tangent constitutive matrix based on the
// current stress state and stress history (depending on the model).
// -------------------------------------------------------------------------
//
// virtual void UpdateState(void)
//
// This method updates the internal variables of the model. It should be
// called after the convergence of the global analysis (equilibrium)
// algorithm. These variables are used to store the stress history of the
// the corresponding integration point.
// -------------------------------------------------------------------------

#ifndef _CMODEL_H
#define _CMODEL_H

#include "vec.h"

// ------------------------------------------------------------------------
// Forward declarations:
//
class cElement;
class cAnModel;
class cMaterial;
class cMatMec;
class cMatrix;

// ------------------------------------------------------------------------
// Definition of the Constitutive Model class:
//
class cConstModel
{
 protected:
  cMatMec *Mat;

 public:
  static cConstModel *Create(cElement *, cMaterial *, double *p = 0);

                cConstModel(cMaterial *);
                cConstModel(cElement *, cMaterial *);
                cConstModel(cElement *, cMaterial *, double *);
  virtual      ~cConstModel(void);
  virtual int   Stress(cVector &, cVector &) = 0;
  virtual void  TangMat(cMatrix &) = 0;
  virtual void  CondMatrix(cMatrix&) { }
  virtual void  UpdateState(void) = 0;
  virtual int   NumStrainLim(void) { return 0; }
  virtual void  GetStrainLim(cVector &) { }
  virtual int   GetCrvDegree(int) { return 1; }
  virtual void  AlphaVec(cVector &) { }
  virtual int   Stress(double, cVector &, cVector &, cVector &) { return 1; }
};


// ------------------------------------------------------------------------
// Definition of the (Isotropic) Linear Constitutive Model class:
//
class cModLinear : public cConstModel
{
 protected:
  cAnModel  *Anm;

 public:
                cModLinear(cElement *, cMaterial *);
  virtual      ~cModLinear(void);
  virtual int   Stress(cVector &, cVector &);
  virtual void  CondMatrix(cMatrix&);
  virtual void  TangMat(cMatrix &);
  virtual void  UpdateState(void) { }
  virtual void  AlphaVec(cVector &);
  virtual int   Stress(double, cVector &, cVector &, cVector &);
};


// ------------------------------------------------------------------------
// Definition of the Orthotropic Linear Constitutive Model class:
//
class cModOrthoLinear : public cModLinear
{
 public:
                cModOrthoLinear(cElement *, cMaterial *);
  virtual      ~cModOrthoLinear(void);
          void  TangMat(cMatrix &);
          void  AlphaVec(cVector &);
};


// ------------------------------------------------------------------------
// Definition of the FGM Model class:
//
class cModFGM : public cConstModel
{
 protected:
  cAnModel  *Anm;
  double     V2;   // Volume fraction of material 2

 public:
                cModFGM(cElement *, cMaterial *, double *);
  virtual      ~cModFGM(void);
  virtual int   Stress(cVector &, cVector &);
  virtual void  TangMat(cMatrix &);
  virtual void  UpdateState(void) { }
  virtual void  AlphaVec(cVector &);
};


// ------------------------------------------------------------------------
// Definition of the Viscoelastic Constitutive Model class:
//
class cModViscoElastic : public cConstModel
{
 protected:
  cAnModel *Anm;
  double    Time;
  cVector   Sig,SigIter;  // Converged and interactive stresses
  cVector   Eps,EpsIter;  // Converged and interactive strains
  cVector   EpsRate;      // Strain rate at the previous step
  cVector  *S;            // Array of hereditary vectors

 public:
                cModViscoElastic(cElement *, cMaterial *);
  virtual      ~cModViscoElastic(void);
          int   Stress(cVector &, cVector &);
          void  TangMat(cMatrix &);
          void  UpdateState(void);
};


// ------------------------------------------------------------------------
// Definition of Composite Progressive Failure Constitutive Model class:
//
class cModProgFailure : public cConstModel
{
 protected:
  cAnModel *Anm;
  double   *Strength;    // Material strength parameters
  double   *Param;       // Material parameters
  double   *ParamIter;   // Iterative material parameters

 public:
                cModProgFailure(cElement *, cMaterial *);
  virtual      ~cModProgFailure(void);
          int   Stress(cVector &, cVector &);
          void  TangMat(cMatrix &);
  virtual void  UpdateState(void) = 0;
  virtual void  InitIterParam(void) = 0;
  virtual bool  FailCrit(cVector &) = 0;
  virtual void  Degrade(void) = 0;
};


// ------------------------------------------------------------------------
// Definition of Hashin Progressive Failure Constitutive Model class:
//
class cModProgFailHashin : public cModProgFailure
{
 protected:
  double    Alpha;        // Degradation factor
  bool     *FailedIter;   // Iterative failure flags

 public:
                cModProgFailHashin(cElement *, cMaterial *);
  virtual      ~cModProgFailHashin(void);
          bool  FailCrit(cVector &);
          void  Degrade(void);
          void  UpdateState(void);
          void  InitIterParam(void);
};


// ------------------------------------------------------------------------
// Definition of the Tsai Progressive Failure Constitutive Model class:
//
class cModProgFailTsai : public cModProgFailure
{
 protected:
  double  Em;           // Matrix degradation factor
  double  Ef;           // Fiber degradation factor
  double  N;            // Power law factor for Xc reduction
  double  Ic;           // Stress interaction coefficient
  double  IcIter;       // Iterative stress interaction coefficient
  bool   *Failed;       // Failure flags
  bool   *FailedIter;   // Iterative failure flags
  bool   *FailMode;     // Failure mode identifier

 public:
                cModProgFailTsai(cElement *, cMaterial *);
  virtual      ~cModProgFailTsai(void);
          bool  FailCrit(cVector &);
          void  Degrade(void);
          void  UpdateState(void);
          void  InitIterParam(void);
};


// ------------------------------------------------------------------------
// Definition of the Engelstad Progressive Failure Constitutive Model class:
//
class cModProgFailEngelstad : public cModProgFailure
{
 protected:
  double  Alpha;        // Degradation parameter
  double  Ic;           // Stress interaction coefficient
  bool   *Failed;       // Failure flags 
  bool   *FailedIter;   // Iterative failure flags
  
 public:
                cModProgFailEngelstad(cElement *, cMaterial *);
  virtual      ~cModProgFailEngelstad(void);
          bool  FailCrit(cVector &);
          void  Degrade(void);
          void  UpdateState(void);
          void  InitIterParam(void);
};

// ------------------------------------------------------------------------
// Definition of the CZM ModeI Bilinear Constitutive Model class:
//
class cModCZMModeIBilinear : public cConstModel
{
 typedef enum
 {
   COMPRESSION,
   ELASTIC,
   DAMAGE,
   FAILURE
 } eStatus;

 protected:
  cAnModel *Anm;
  eStatus   Status,StatusIter;
  double    DamVar,DamVarIter;    // Degradation parameter
  double    MaxDisp,MaxDispIter;  // Maximum displacement (elastic limit)

 public:
                cModCZMModeIBilinear(cElement *, cMaterial *);
  virtual      ~cModCZMModeIBilinear(void);
          int   Stress(cVector &, cVector &);
          void  TangMat(cMatrix &);
          void  UpdateState(void);
};

// ------------------------------------------------------------------------
// Definition of the CZM ModeI Exponential Constitutive Model class:
//
class cModCZMModeIExponential : public cConstModel
{
 typedef enum
 {
   COMPRESSION,
   ELASTIC,
   DAMAGE
 } eStatus;

 protected:
  cAnModel *Anm;
  eStatus   Status,StatusIter;
  double    DamVar,DamVarIter;    // Degradation parameter
  double    MaxDisp,MaxDispIter;  // Maximum displacement (elastic limit)

 public:
                cModCZMModeIExponential(cElement *, cMaterial *);
  virtual      ~cModCZMModeIExponential(void);
          int   Stress(cVector &, cVector &);
          void  TangMat(cMatrix &);
          void  UpdateState(void);
};


// ------------------------------------------------------------------------
// Definition of the von Mises Elastoplastic Constitutive Model class:
//
class cModPlastIsoLinHard : public cConstModel
{
 protected:
  typedef enum
  {
   ELASTIC,
   PLASTIC
  } eStatus;

  cAnModel  *Anm;
  eStatus    Status;                  // Loading status
  double     Dgam;                    // Plastic multiplier
  double     EqvPlEps,EqvPlEpsIter;   // Equivalent plastic strains
  cVector    PlastEps,PlastEpsIter;   // Plastic strains
  cVector    StressIter;              // Current stresses

 protected:
  virtual int   StressPlast(double *, cVector &, cVector &);
  virtual int   StressPlastPS(double *, cVector &, cVector &);
  virtual void  TangMatPlast(double *, cMatrix &);

 public:
                cModPlastIsoLinHard(cElement *, cMaterial *);
  virtual      ~cModPlastIsoLinHard(void);
  virtual int   Stress(cVector &, cVector &);
  virtual void  TangMat(cMatrix &);
  virtual void  UpdateState(void);
};


// ------------------------------------------------------------------------
// Definition of the FGM-TTO Model class:
//
class cModFGMTTO : public cModPlastIsoLinHard
{
 protected:
  double V2;    // Volume fraction of material 2

 protected:
  virtual void  Homogenize(double *, double *);

 public:
                cModFGMTTO(cElement *, cMaterial *, double *);
  virtual      ~cModFGMTTO(void);
  virtual int   Stress(cVector &, cVector &);
  virtual void  TangMat(cMatrix &);
};

#endif
