// -------------------------------------------------------------------------
// anmodel.h - file containing the definition of the Analysis Model class.
// -------------------------------------------------------------------------
//
// The Analysis Model class handles the element features related to the
// differential equation/variation statement governing the problem to be
// solved, as the number and type of degrees of freedom of each node,
// the number and label of stress components, and the structure of the
// strain-desplacement (B) and constitutive (Q) matrices, to cite a few
// examples.
//
// Analysis Model
// |-- Mechanical Model
// |   |-- Plane Stress
// |   |-- Plane Strain
// |   |-- Axisymmetric
// |   |-- Solid
// |   |-- Thick Plate (Reissner-Mindlin)
// |   |-- HSDT Plate (von Karman membrane strains + HSDT)
// |   |-- Shallow Shell (Marguerre membrane strains + Reissner-Mindlin shear)
// |-- Donnell Shell (Donnell membrane/curvature + Reissner-Mindlin shear)
// |   |-- Shell (Degenerated)
// |   |-- Interface2D
// |   |-- Interface3D
// |-- Thermal Model
// |   |-- Plane Heat Transfer
//
// -------------------------------------------------------------------------
// Static public methods:
// -------------------------------------------------------------------------
//
// static cAnModel *CreateModel(eAnmType type)
//
//   type  -  analysis model type                                     (in)
//
// This method creates an Analysis Model object and returns a pointer to it.
//
// -------------------------------------------------------------------------
// Shared methods (implemented only in the base class):
// -------------------------------------------------------------------------
//
// eAnmType GetType(void)
//
// This method returns the Analysis Model type.
//
// -------------------------------------------------------------------------
// Pure virtual methods:
// -------------------------------------------------------------------------
//
// virtual int GetNumDofNode(void)
//
// This method returns the number of degrees of freedom per node.
// -------------------------------------------------------------------------
//
// virtual int GetActDir(int *dir)
//
//  dir    -  array of the active directions (dofs)                  (out)
//
// This method returns the active directions (degrees of freedom) of a node.
// As an example, for plane stress problems the active directions are
// dir[0] = dir[1] = 1 and dir[2] = dir[3] = dir[4] = dir[5] = 0,
// since only the displacements u,v are active.
// -------------------------------------------------------------------------
//
// virtual int GetDimQMatrix(void)
//
// This method returns the dimension of the constitutive matrix [Q]. This
// size corresponds to the number of stress and strain components of the
// underlying constitutive model (e.g. 3 for planess stress, 6 for solid, 
// and 5 for shell models).
// -------------------------------------------------------------------------
//
// virtual int GetDimBMatrix(void)
//
// This method returns the number of rows of the strain-displacement matrix
// [B]. This number is identical to the number of strain/stress components
// continuum models (e.g. 3 for plane strain and 6 for solid models) and 
// equal to the number of generalized strain/stresses for structural models 
// (e.g. 5 for shallow shell and 6 for solid).
// -------------------------------------------------------------------------
//
// virtual int GetDimBnlMatrix(void)
//
// This method returns the number of rows of the nonlinear strain-displacement
// matrix [Bnl] used for stability and large displacement analyses.
// -------------------------------------------------------------------------
//
// virtual void GetStrLabels(int *label)
//
//  label  -  array of response labels                               (out)
//
// This method returns the indices of the element stress component labels
// in the global vector VEC_SCL_LAB.
// -------------------------------------------------------------------------
//
// virtual void PrincStress(cVector &str, cVector &prs)
//
//  str    -  stress vector in cartesian coordinates                  (in)
//  prs    -  principal stresses                                     (out)
//
// This method computes the principal stresses from a given stress state in
// cartesian coordinates.
// -------------------------------------------------------------------------
//
// virtual double VolCoeff(double t,int nn,double *mapfnc,sNodeCoord *coord)
//
//  t      -  thickness                                               (in)
//  nn     -  number of element nodes                                 (in)
//  mapfnc -  mapping functions                                       (in)
//  coord  -  nodal coordinates                                       (in)
//
// This method returns the coefficient that multiply the Jacobian
// determinant (|J|) in order to obtain the differential volume (dV) used
// in the integration of the element matrices. Some examples:
// plane stress  dV = t*|J| drds,  it returns the thickness t.
// axisymmetric  dV = 2*PI*R*|J| drds, it returns 2*PI*R.
// solids        dV = |J| drdsdt, it returns 1.
// -------------------------------------------------------------------------
//
// virtual void QMatrix(double *param, cMatrix &Q)
//
//  param  -  material parameters (E, nu)                             (in)
//  Q      -  elastic constitutive matrix                            (out)
//
// This method computes the elastic constitutive matrix using the given
// material parameters.
// -------------------------------------------------------------------------
//
// virtual void QMatrixOrtho(double *param, cMatrix &Q)
//
//  param  -  material parameters                                     (in)
//  Q      -  elastic constitutive matrix                            (out)
//
// This method computes the elastic constitutive matrix for orthotropic
// materials using the given parameters (E1, E2, ...).
// -------------------------------------------------------------------------
//
// virtual void DMatrix(double *param, cMatrix &D)
//
//  param  -  material parameters (E, nu)                             (in)
//  D      -  compliance matrix                                      (out)
//
// This method computes the compliance matrix (D = invQ) using the given
// material parameters.
// -------------------------------------------------------------------------
//
// virtual void TMatrix(double angle, cMatrix &T)
//
//  angle  -  angle between the local and global x-axes               (in)
//  T      -  strain transformation matrix                           (out)
//
// This method computes the strain transformation matrix from the global
// to the local (material) system {epsloc} = [T]*{eps} for plane models.
// -------------------------------------------------------------------------
//
// virtual void TMatrix(cVector &e1, cVector &e2, cVector &e3, cMatrix &T)
//
//  e1,e2,e3  -  material axes in the global system                   (in)
//  T         -  strain transformation matrix                        (out)
//
// This method computes the strain transformation matrix from the global
// to the local (material) system {epsloc} = [T]*{eps}.
// -------------------------------------------------------------------------
//
// virtual void  BMatrix(int nn, double *shpfunc, double *mapfunc,
//                       sNodeCoord *coord, sNodeCoord *drvshp, cMatrix &B)
//
//  nn      - number of nodes                                         (in)
//  shpfunc - shape functions                                         (in)
//  mapfunc - mapping functions                                       (in)
//  coord   - nodal coordinates                                       (in)
//  drvshp  - derivatives of the shape functions w.r.t. xyz           (in)
//  B       - strain-displacement matrix                             (out)
//
// This method computes the strain-displacement matrix.
// -------------------------------------------------------------------------
//
// virtual void  NMatrix(int nn, double *shpfunc, cMatrix &N)
//
//  nn      - number of nodes                                         (in)
//  shpfunc - shape functions                                         (in)
//  N       - displacement interpolation matrix                      (out)
//
// This method computes the displacement interpolation matrix.
// -------------------------------------------------------------------------
//
// virtual void  NMatrix(int nn, double *shpfunc, sNodeCoord *drvshp, cMatrix &N)
//
//  nn      - number of nodes                                         (in)
//  shpfunc - shape functions                                         (in)
//  drvshp  - derivatives of the shape functions w.r.t. xyz           (in)
//  N       - displacement interpolation matrix                      (out)
//
// This method computes the displacement interpolation matrix.
// -------------------------------------------------------------------------
//
// virtual void  ExpandStress(cVector &sig, cVector &expsig)
//
//  sig     - stress vector with variable size                        (in)
//  expsig  - stress vector with fixed (6) size                      (out)
//
// This method expands a given stress vector, which size depends on the
// used analysis model, to a 3D stress vector with 6 components.
// ------------------------------------------------------------------------
//
// virtual void  CalcDev(double p, cVector &sig, cVector &sdev)
//  
//  p       - hidrostatic stress (p = I1/3)                           (in)
//  sig     - stress vector                                           (in)
//  sdev    - deviatoric stress vector                               (out)
// 
// This method returns the deviatoric stress vector given the stress vector
// and the hydrostatic stress.
// ------------------------------------------------------------------------
//
// virtual void  CalcTen(double p, cVector &sdev, cVector &sig)
//
//  p       - hydrostatic stress (p = I1/3)                           (in)
//  sdev    - deviatoric stress vector                                (in)
//  sig     - stress vector                                          (out) 
//
// This method evaluates the resultant stress vector given the deviatoric
// stress and the hydrostatic stress.
// ------------------------------------------------------------------------
//
// virtual void  CalcFOID(cMatrix &FO)
//
//  FO      - fourth-order identity tensor                           (out)
//
// This method assembles the fourth-order symmetric identity tensor stored
// in matrix form.
// ------------------------------------------------------------------------
//
// virtual void  CalcSOID(cVector &SO)
//
//  SO      - second-order identity tensor                           (out)
//
// This method evaluates the second-order identity tensor stored in vctor
// form.
// ------------------------------------------------------------------------
//
// virtual double  CalcI1(cVector &sig)
//
//  sig     - stress vector                                           (in)
//
// This method returns the value of the first stress invariant (I1),  
// corresponding to the trace of the stress tensor.
// ------------------------------------------------------------------------
//
// virtual double  CalcJ2(cVector &sig) 
//
//  sig     - stress vector                                           (in)
//
// This method returns the value of the second invariant of the deviatoric
// stresses (J2).
// ------------------------------------------------------------------------


#ifndef _ANMODEL_H
#define _ANMODEL_H

#include "node.h"

// ------------------------------------------------------------------------
// Forward declarations:
//
class cVector;
class cMatrix;
class cMecModel;
class cHeatModel;

// -------------------------------------------------------------------------
// Analysis model types:
//
typedef enum
{
  PLANE_STRESS,
  PLANE_STRAIN,
  AXISYMMETRIC,
  SOLID,
  THICK_PLATE,
  HSDT_PLATE,
  SHALLOW_SHELL,
  DONNELL_SHELL,
  SHELL,
  INTERFACE_2D,
  INTERFACE_3D,
  PLANE_HEAT_TRANSFER
} eAnmType;

// ------------------------------------------------------------------------
// Definition of the Analysis Model class:
//
/*  TODO apagar!
class cAnModel
{
 protected:
  eAnmType Type;

 public:
  static cAnModel *CreateModel(eAnmType);

                    cAnModel(void);
  virtual          ~cAnModel(void);
          eAnmType  GetType(void) { return Type; }
  virtual int       GetNumDofNode(void) = 0;
  virtual void      GetActDir(int *) = 0;
  virtual void      GetStrLabels(int *) = 0;
  virtual double    VolCoeff(double, int, double *, sNodeCoord *) = 0;
  virtual void      PrincStress(cVector &, cVector &) = 0;
  virtual double    CalcI1(cVector &) { return 0.0; }
  virtual double    CalcJ2(cVector &) { return 0.0; }
  virtual void      GradI1(cVector &, cVector &) { }
  virtual void      GradJ2(cVector &, cVector &) { }
  virtual void      QMatrix(double*,cMatrix&) = 0;
  virtual void      QMatrixOrtho(double*,cMatrix&) = 0;
  virtual void      DMatrix(double *, cMatrix &) = 0;
  virtual void      TMatrix(double, cMatrix &) { }
  virtual void      TMatrix(cVector &, cVector &, cVector &, cMatrix &) { }
  virtual void      BMatrix(int, double *, double *,
                            sNodeCoord *, sNodeCoord *, cMatrix &) = 0;
  virtual void      BlMatrix(int, double *, double *, cVector &,
                             sNodeCoord *, sNodeCoord *, cMatrix &) = 0;
  virtual void      BnlMatrix(int, double *, double *,
                              sNodeCoord *, sNodeCoord *, cMatrix &) = 0;
  virtual void      SMatrix(cVector &, cMatrix &) = 0;
  virtual void      NMatrix(int, double *, cMatrix &) = 0;
  virtual void      ExpandStress(cVector &, cVector &) { };
};
*/

// ------------------------------------------------------------------------
// Definition of the Analysis Model class:
//
class cAnModel
{
 protected:
  eAnmType Type;

 public:
  static cAnModel *CreateModel(eAnmType);

                      cAnModel(void);
  virtual            ~cAnModel(void);
          eAnmType    GetType(void) { return Type; }
	  cMecModel*  GetMecMod(void);
	  cHeatModel* GetHeatMod(void);
  virtual int         GetNumDofNode(void) = 0;
  virtual void        GetActDir(int *) = 0;
  virtual int         GetDimBMatrix(void) = 0;
  virtual void        GetStrLabels(int *) = 0;
  virtual int         GetDimNMatrix(void)  { return GetNumDofNode( ); }
  virtual int         GetNumNodStrCmp(void)  { return GetDimBMatrix( ); }
  virtual void        GetNodStrLabels(int*l) { return GetStrLabels(l); }
  virtual double      VolCoeff(double,int,double *,sNodeCoord *) = 0;
  virtual void        BMatrix(int, double *, double *,
                              sNodeCoord *, sNodeCoord *, cMatrix &) = 0;
  virtual void        BMatrix(double, int, double *, double *,
                              sNodeCoord *, sNodeCoord *, sNodeDrv *, cMatrix &) { }
};


// ------------------------------------------------------------------------
// Definition of the Mechanical Analysis Model class:
//
class cMecModel : public cAnModel
{
 public:
                    cMecModel(void);
  virtual          ~cMecModel(void);
  virtual int       GetDimQMatrix(void) = 0;
  virtual int       GetDimBMatrix(void) = 0;
  virtual int       GetDimBnlMatrix(void);
  virtual void      PrincStress(cVector &, cVector &) = 0;
  virtual double    CalcI1(cVector &) { return 0.0; }
  virtual double    CalcJ2(cVector &) { return 0.0; }
  virtual void      GradI1(cVector &, cVector &) { }
  virtual void      GradJ2(cVector &, cVector &) { }
  virtual void      QMatrix(double *, cMatrix &) = 0;
  virtual void      QMatrixOrtho(double *, cMatrix &) = 0;
  virtual void      DMatrix(double *, cMatrix &) = 0;
  virtual void      TMatrix(double, cMatrix &) { }
  virtual void      TMatrix(cVector &, cVector &, cVector &, cMatrix &) { }
  virtual void      BlMatrix(int, double *, double *, cVector &,
                             sNodeCoord *, sNodeCoord *, cMatrix &) = 0;
  virtual void      BnlMatrix(int, double *, double *,
                              sNodeCoord *, sNodeCoord *, cMatrix &) = 0;
  virtual void      SMatrix(cVector &, cMatrix &) = 0;
  virtual void      NMatrix(int, double *, cMatrix &);
  virtual void      NMatrix(int, double *, sNodeCoord *, cMatrix &);
  virtual void      ExpandStress(cVector &, cVector &) { }
  virtual void      CalcDev(double, cVector &, cVector &) { }
  virtual void      CalcTen(double, cVector &, cVector &) { }
  virtual void      CalcFOID(cMatrix &) { }
  virtual void      CalcSOID(cVector &) { }
  virtual void      AlphaVec(double *, cVector &);
  virtual void      AlphaVecOrtho(double *, cVector &);
};


// ------------------------------------------------------------------------
// Definition of the Thermal Analysis Model class:
//
class cHeatModel : public cAnModel
{
 public:
                    cHeatModel(void);
  virtual          ~cHeatModel(void);
  virtual void      CondMatrix(double*,cMatrix &) = 0;
  virtual void      CondMatrixOrtho(double*,cMatrix &) = 0;
};


// ------------------------------------------------------------------------
// Definition of the Plane Stress class:
//
class cPlaneStress : public cMecModel
{
 public:
         cPlaneStress(void);
        ~cPlaneStress(void);
  int    GetNumDofNode(void) { return 2; }
  int    GetDimQMatrix(void) { return 3; }
  int    GetDimBMatrix(void) { return 3; }
  int    GetDimBnlMatrix(void) { return 4; }

  void   GetActDir(int *);
  void   GetStrLabels(int *);
  void   PrincStress(cVector &, cVector &);
  double CalcI1(cVector &);
  double CalcJ2(cVector &);
  void   GradI1(cVector &, cVector &);
  void   GradJ2(cVector &, cVector &);
  double VolCoeff(double t, int, double *, sNodeCoord *) { return t; }
  void   QMatrix(double *, cMatrix &);
  void   QMatrixOrtho(double *, cMatrix &);
  void   DMatrix(double *, cMatrix &);
  void   TMatrix(double, cMatrix &);
  void   BMatrix(int, double *, double *, sNodeCoord *,
                 sNodeCoord *, cMatrix &);
  void   BlMatrix(int, double *, double *, cVector &,
                  sNodeCoord *, sNodeCoord *, cMatrix &);
  void   BnlMatrix(int, double *, double *,
                   sNodeCoord *, sNodeCoord *, cMatrix &);
  void   SMatrix(cVector &, cMatrix &);
  void   NMatrix(int, double *, cMatrix &);
  void   ExpandStress(cVector &, cVector &);
  void   AlphaVec(double *, cVector &);
  void   AlphaVecOrtho(double *, cVector &);
};


// ------------------------------------------------------------------------
// Definition of the Plane Strain class:
//
class cPlaneStrain : public cMecModel
{
 public:
         cPlaneStrain(void);
        ~cPlaneStrain(void);
  int    GetNumDofNode(void) { return 2; }
  int    GetDimQMatrix(void) { return 4; }
  int    GetDimBMatrix(void) { return 4; }
  int    GetDimBnlMatrix(void) { return 5; }
  void   GetActDir(int *);
  void   GetStrLabels(int *);
  void   PrincStress(cVector &, cVector &);
  double CalcI1(cVector &);
  double CalcJ2(cVector &);
  void   GradI1(cVector &, cVector &);
  void   GradJ2(cVector &, cVector &);
  double VolCoeff(double, int, double *, sNodeCoord *) { return 1.0; }
  void   QMatrix(double *, cMatrix &);
  void   QMatrixOrtho(double *, cMatrix &);
  void   DMatrix(double *, cMatrix &);
  void   BMatrix(int, double *, double *, sNodeCoord *,
                 sNodeCoord *, cMatrix &);
  void   BlMatrix(int, double *, double *, cVector &,
                  sNodeCoord *, sNodeCoord *, cMatrix &);
  void   BnlMatrix(int, double *, double *,
                   sNodeCoord *, sNodeCoord *, cMatrix &);
  void   SMatrix(cVector &, cMatrix &);
  void   NMatrix(int, double *, cMatrix &);
  void   ExpandStress(cVector &, cVector &);
  void   CalcDev(double, cVector &, cVector &);
  void   CalcTen(double, cVector &, cVector &);
  void   CalcFOID(cMatrix &);
  void   CalcSOID(cVector &);
};


// ------------------------------------------------------------------------
// Definition of the Axisymmetric solid class:
//
class cAxisymmetric : public cMecModel
{
 public:
         cAxisymmetric(void);
        ~cAxisymmetric(void);
  int    GetNumDofNode(void) { return 2; }
  int    GetDimQMatrix(void) { return 4; }
  int    GetDimBMatrix(void) { return 4; }
  int    GetDimBnlMatrix(void) { return 5; }
  void   GetActDir(int *);
  void   GetStrLabels(int *);
  void   PrincStress(cVector &, cVector &);
  double CalcI1(cVector &);
  double CalcJ2(cVector &);
  void   GradI1(cVector &, cVector &);
  void   GradJ2(cVector &, cVector &);
  double VolCoeff(double, int, double *, sNodeCoord *);
  void   QMatrix(double *, cMatrix &);
  void   QMatrixOrtho(double *, cMatrix &);
  void   DMatrix(double *, cMatrix &);
  void   BMatrix(int, double *, double *, sNodeCoord *,
                 sNodeCoord *, cMatrix &);
  void   BlMatrix(int, double *, double *, cVector &,
                  sNodeCoord *, sNodeCoord *, cMatrix &);
  void   BnlMatrix(int, double *, double *,
                   sNodeCoord *, sNodeCoord *, cMatrix &);
  void   SMatrix(cVector &, cMatrix &);
  void   NMatrix(int, double *, cMatrix &);
  void   ExpandStress(cVector &, cVector &);
  void   CalcDev(double, cVector &, cVector &);
  void   CalcTen(double, cVector &, cVector &);
  void   CalcFOID(cMatrix &);
  void   CalcSOID(cVector &);
};


// ------------------------------------------------------------------------
// Definition of the Solid class:
//
class cSolid : public cMecModel
{
 public:
         cSolid(void);
        ~cSolid(void);
  int    GetNumDofNode(void) { return 3; }
  int    GetDimQMatrix(void) { return 6; }
  int    GetDimBMatrix(void) { return 6; }
  int    GetDimBnlMatrix(void) { return 9; }
  void   GetActDir(int *);
  void   GetStrLabels(int *);
  void   PrincStress(cVector &, cVector &);
  double CalcI1(cVector &);
  double CalcJ2(cVector &);
  void   GradI1(cVector &, cVector &);
  void   GradJ2(cVector &, cVector &);
  double VolCoeff(double, int, double *, sNodeCoord *) { return 1.0; }
  void   QMatrix(double *, cMatrix &);
  void   QMatrixOrtho(double *, cMatrix &);
  void   DMatrix(double *, cMatrix &);
  void   TMatrix(cVector &, cVector &, cVector &, cMatrix &);
  void   BMatrix(int, double *, double *, sNodeCoord *,
                 sNodeCoord *, cMatrix &);
  void   BlMatrix(int, double *, double *, cVector &,
                  sNodeCoord *, sNodeCoord *, cMatrix &);
  void   BnlMatrix(int, double *, double *,
                   sNodeCoord *, sNodeCoord *, cMatrix &);
  void   SMatrix(cVector &, cMatrix &);
  void   NMatrix(int, double *, cMatrix &);
  void   ExpandStress(cVector &, cVector &);
  void   CalcDev(double, cVector &, cVector &);
  void   CalcTen(double, cVector &, cVector &);
  void   CalcFOID(cMatrix &);
  void   CalcSOID(cVector &);
};


// ------------------------------------------------------------------------
// Definition of the Thick Plate (Reissner-Mindlin) class:
//
class cThickPlate : public cMecModel
{
 public:
         cThickPlate(void);
        ~cThickPlate(void);
  int    GetNumDofNode(void) { return 3; }
  int    GetDimQMatrix(void) { return 5; }
  int    GetDimBMatrix(void) { return 5; }
  void   GetActDir(int *);
  void   GetStrLabels(int *);
  void   PrincStress(cVector &, cVector &);
  double VolCoeff(double, int, double *, sNodeCoord *) { return 1.0; }
  void   QMatrix(double *, cMatrix &);
  void   QMatrixOrtho(double *, cMatrix &);
  void   DMatrix(double *, cMatrix &);
  void   TMatrix(double, cMatrix &);
  void   BMatrix(int, double *, double *, sNodeCoord *,
                 sNodeCoord *, cMatrix &);
  void   BlMatrix(int, double *, double *, cVector &,
                  sNodeCoord *, sNodeCoord *, cMatrix &);
  void   BnlMatrix(int, double *, double *,
                   sNodeCoord *, sNodeCoord *, cMatrix &);
  void   SMatrix(cVector &, cMatrix &);
  void   NMatrix(int, double *, cMatrix &);
  void   ExpandStress(cVector &, cVector &);
};


// ------------------------------------------------------------------------
// Definition of the HSDT Plate class:
//
class cHSDTPlate : public cMecModel
{
 public:
         cHSDTPlate(void);
        ~cHSDTPlate(void);
  int    GetNumDofNode(void) { return 5; }
  int    GetDimNMatrix(void) { return 7; }
  int    GetDimQMatrix(void) { return 5; }
  int    GetDimBMatrix(void) { return 11; }
  int    GetDimBnlMatrix(void) { return 2; }
  void   GetActDir(int *);
  void   GetStrLabels(int *);
  void   PrincStress(cVector &, cVector &);
  double VolCoeff(double, int, double *, sNodeCoord *) { return 1.0; }
  void   QMatrix(double *, cMatrix &);
  void   QMatrixOrtho(double *, cMatrix &);
  void   DMatrix(double *, cMatrix &);
  void   TMatrix(double, cMatrix &);
  void   BMatrix(int, double *, double *, sNodeCoord *,
                 sNodeCoord *, cMatrix &) { }
  void   BMatrix(double, int, double *, double *, sNodeCoord *,
                 sNodeCoord *, sNodeDrv *, cMatrix &);
  void   BlMatrix(int, double *, double *, cVector &,
                 sNodeCoord *, sNodeCoord *, cMatrix &);
  void   BnlMatrix(int, double *, double *,
                 sNodeCoord *, sNodeCoord *, cMatrix &);
  void   SMatrix(cVector &, cMatrix &);
  void   NMatrix(int, double *, sNodeCoord *, cMatrix &);
  void   ExpandStress(cVector &, cVector &);
};


// ------------------------------------------------------------------------
// Definition of the Shallow Shell class:
//
class cShallowShell : public cMecModel
{
 public:
         cShallowShell(void);
        ~cShallowShell(void);
  int    GetNumDofNode(void) { return 5; }
  int    GetDimQMatrix(void) { return 5; }
  int    GetDimBMatrix(void) { return 8; }
  int    GetDimBnlMatrix(void) { return 2; }
  void   GetActDir(int *);
  void   GetStrLabels(int *);
  void   PrincStress(cVector &, cVector &);
  double VolCoeff(double, int, double *, sNodeCoord *) { return 1.0; }
  void   QMatrix(double *, cMatrix &);
  void   QMatrixOrtho(double *, cMatrix &);
  void   DMatrix(double *, cMatrix &);
  void   TMatrix(double, cMatrix &);
  void   BMatrix(int, double *, double *, sNodeCoord *,
                 sNodeCoord *, cMatrix &);
  void   BlMatrix(int, double *, double *, cVector &,
                 sNodeCoord *, sNodeCoord *, cMatrix &);
  void   BnlMatrix(int, double *, double *,
                 sNodeCoord *, sNodeCoord *, cMatrix &);
  void   SMatrix(cVector &, cMatrix &);
  void   NMatrix(int, double *, cMatrix &);
  void   ExpandStress(cVector &, cVector &);
  void   AlphaVec(double *, cVector &);
  void   AlphaVecOrtho(double *, cVector &);
};


// ------------------------------------------------------------------------
// Definition of the Donnell Shell class:
//
class cDonnellShell : public cMecModel
{
 public:
         cDonnellShell(void);
        ~cDonnellShell(void);
  int    GetNumDofNode(void) { return 5; }
  int    GetDimQMatrix(void) { return 5; }
  int    GetDimBMatrix(void) { return 8; }
  int    GetDimBnlMatrix(void) { return 2; }
  void   GetActDir(int *);
  void   GetStrLabels(int *);
  void   PrincStress(cVector &, cVector &);
  double VolCoeff(double, int, double *, sNodeCoord *) { return 1.0; }
  void   QMatrix(double *, cMatrix &);
  void   QMatrixOrtho(double *, cMatrix &);
  void   DMatrix(double *, cMatrix &);
  void   TMatrix(double, cMatrix &);
  void   BMatrix(int, double *, double *, sNodeCoord *,
                 sNodeCoord *, cMatrix &);
  void   BlMatrix(int, double *, double *, cVector &,
                 sNodeCoord *, sNodeCoord *, cMatrix &);
  void   BnlMatrix(int, double *, double *,
                 sNodeCoord *, sNodeCoord *, cMatrix &);
  void   SMatrix(cVector &, cMatrix &);
  void   NMatrix(int, double *, cMatrix &);
  void   ExpandStress(cVector &, cVector &);
};

// ------------------------------------------------------------------------
// Definition of the Shell class:
//
class cShell : public cMecModel
{
 public:
         cShell(void);
        ~cShell(void);
  int    GetNumDofNode(void) { return 5; }
  int    GetDimQMatrix(void) { return 6; }
  int    GetDimBMatrix(void) { return 6; }
  int    GetDimBnlMatrix(void) { return 9; }
  int    GetNumNodStrCmp(void) { return 8; }

  void   GetActDir(int *);
  void   GetStrLabels(int *);
  void   GetNodStrLabels(int *);
  void   PrincStress(cVector &, cVector &);
  double VolCoeff(double, int, double *, sNodeCoord *) { return 1.0; }
  void   QMatrix(double *, cMatrix &);
  void   QMatrixOrtho(double *, cMatrix &);
  void   DMatrix(double *, cMatrix &);
  void   TMatrix(double, cMatrix &);
  void   TMatrix(cVector &, cVector &, cVector &, cMatrix &);
  void   BMatrix(int, double *, double *, sNodeCoord *,
                 sNodeCoord *, cMatrix &);
  void   BlMatrix(int, double *, double *, cVector &,
                 sNodeCoord *, sNodeCoord *, cMatrix &);
  void   BnlMatrix(int, double *, double *,
                 sNodeCoord *, sNodeCoord *, cMatrix &);
  void   SMatrix(cVector &, cMatrix &);
  void   NMatrix(int, double *, cMatrix &);
  void   ExpandStress(cVector &, cVector &);
};


// ------------------------------------------------------------------------
// Definition of the Interface2D class:
//
class cInterface2D : public cMecModel
{
 public:
         cInterface2D(void);
        ~cInterface2D(void);
  int    GetNumDofNode(void) { return 2; }
  int    GetDimQMatrix(void) { return 2; }
  int    GetDimBMatrix(void) { return 2; }
  void   GetActDir(int *);
  void   GetStrLabels(int *);
  void   PrincStress(cVector &, cVector &) { }
  double VolCoeff(double t, int, double *, sNodeCoord *) { return t; }
  void   QMatrix(double *, cMatrix &) { }
  void   QMatrixOrtho(double *, cMatrix &) { }
  void   DMatrix(double *, cMatrix &) { }
  void   BMatrix(int, double *, double *, sNodeCoord *, sNodeCoord *, cMatrix &);
  void   BlMatrix(int, double *, double *, cVector &,
                  sNodeCoord *, sNodeCoord *, cMatrix &) { }
  void   BnlMatrix(int, double *, double *,
                   sNodeCoord *, sNodeCoord *, cMatrix &) { }
  void   SMatrix(cVector &, cMatrix &) { }
  void   NMatrix(int, double *, cMatrix &) { }
};


// ------------------------------------------------------------------------
// Definition of the Interface3D class:
//
class cInterface3D : public cMecModel
{
 public:
         cInterface3D(void);
        ~cInterface3D(void);
  int    GetNumDofNode(void) { return 3; }
  int    GetDimQMatrix(void) { return 3; }
  int    GetDimBMatrix(void) { return 3; }
  void   GetActDir(int *);
  void   GetStrLabels(int *);
  void   PrincStress(cVector &, cVector &) { }
  double VolCoeff(double, int, double *, sNodeCoord *) { return 1.0; }
  void   QMatrix(double *, cMatrix &) { }
  void   QMatrixOrtho(double *, cMatrix &) { }
  void   DMatrix(double *, cMatrix &) { }
  void   BMatrix(int, double *, double *, sNodeCoord *, sNodeCoord *, cMatrix &);
  void   BlMatrix(int, double *, double *, cVector &,
                  sNodeCoord *, sNodeCoord *, cMatrix &) { }
  void   BnlMatrix(int, double *, double *,
                   sNodeCoord *, sNodeCoord *, cMatrix &) { }
  void   SMatrix(cVector &, cMatrix &) { }
  void   NMatrix(int, double *, cMatrix &) { }
};


// ------------------------------------------------------------------------
// Definition of the Plane Heat Transfer class:
//
class cPlaneHeatTransfer : public cHeatModel
{
 public:
         cPlaneHeatTransfer(void);
        ~cPlaneHeatTransfer(void);
  int    GetNumDofNode(void) { return 1; }
  int    GetDimBMatrix(void) { return 2; }

  void   GetActDir(int *);
  void   GetStrLabels(int *);
  double VolCoeff(double t, int, double *, sNodeCoord *) { return t; }
  void   CondMatrix(double*,cMatrix &);
  void   CondMatrixOrtho(double*,cMatrix &);
  void   BMatrix(int, double *, double *, sNodeCoord *,
                 sNodeCoord *, cMatrix &);
  //void   NMatrix(int, double *, cMatrix &); TODO
};

#endif
