// -------------------------------------------------------------------------
// secanalysis.h - definition of the SectionAnalysis class.
// -------------------------------------------------------------------------
//
// The SecAnalysis class carry-out tasks inherent to element cross sections.
// It creates the constitutive models for each integration point and,
// depending on the section type, also creates through-thickness integration
// points and performs the integration of the section stresses and tangent
// stiffness.
//
// SecAnalysis
// |-- SecAnHomogeneous => Plane Stress, Plane Strain, Axisymmetric and Solid
// |-- SecAnHomOrtho (*)
// |-- SecAnHomPlate
// |-- SecAnHomShell
// |-- SecAnFGM2D => Plane Stress, Plane Strain and Axisymmetric
// |-- SecAnFGM3D => Solid
// |-- SecAnFGMShell
// |-- SecAnHSDTPlate
// |-- SecAnHSDTFGMPlate
// |-- SecAnGenShell
// |-- SecAnLaminatedMembr (*) => Plane Stress
// |-- SecAnLaminatedPlate (*)
// |-- SecAnLaminatedShell
// |-- SecAnLaminatedSolid
// |-- SecAnPlFrameRectFiber
// |-- SecAnPlFrameRectGauss
// |-- SecAnPlFrameRCRect
//
// Not implemented yet (*)
// -------------------------------------------------------------------------
// Static public methods:
// -------------------------------------------------------------------------
//
// static cSecAnalysis *Create(cElement *elm, cIntPoint *ipt)
//
// This method creates a section analysis object based on the section type
// associated with the parent element and returns a pointer to the created
// section. It also stores pointers to the element and integration point
// associated with the section.
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
// virtual void  CMatrix(cMatrix &C)
//
//  C      -  tangent constitutive matrix                            (out)
//
// This method computes the tangent constitutive matrix based on the
// current stress state and stress history (depending on the constitutive
// model).
// -------------------------------------------------------------------------
//
// virtual void  MbMatrix(cMatrix &Mb)
//
//  Mb     -  section mass matrix                                    (out)
//
// This method computes the cross-sectional mass matrix according to the
// section geometry and material density.
// -------------------------------------------------------------------------

#ifndef _SECANALYSIS_H
#define _SECANALYSIS_H

#include "vec.h"

// -------------------------------------------------------------------------
// Forward declarations:
//
class cIntPoint;
class cElement;
class cConstModel;
class cSection;
class cAnModel;
class cIntPoint;

// -------------------------------------------------------------------------
// Definition of the Section Analysis class:
//
class cSecAnalysis
{
 protected:
  cElement  *Elm;
  cSection  *Sec;

 public:
  static        cSecAnalysis *Create(cElement *, cIntPoint *);

                cSecAnalysis(cElement *, cIntPoint *);
  virtual      ~cSecAnalysis(void);
  virtual int   Stress(cVector &, cVector &) = 0;
  virtual int   Stress(double, cVector &, cVector &, cVector &);
  virtual void  CMatrix(cMatrix &) = 0;
  virtual void  CondMatrix(cMatrix &){ }
  virtual void  UpdateState(void) { }
  virtual void  MbMatrix(cMatrix &);
};


// -------------------------------------------------------------------------
// Definition of the Homogeneous Section Analysis class:
//
class cSecAnHomogeneous : public cSecAnalysis
{
 protected:
  cConstModel *Cmod;

 public:
                cSecAnHomogeneous(cElement *, cIntPoint *);
  virtual      ~cSecAnHomogeneous(void);
          int   Stress(cVector &, cVector &);
          int   Stress(double, cVector &, cVector &, cVector &);
          void  CMatrix(cMatrix &);
          void  CondMatrix(cMatrix&);
          void  UpdateState(void);
          void  MbMatrix(cMatrix &);
};


// -------------------------------------------------------------------------
// Definition of the Homogeneous Plate Section Analysis class:
//
class cSecAnHomPlate : public cSecAnalysis
{
 protected:
  cConstModel *Cmod;

 public:
                cSecAnHomPlate(cElement *, cIntPoint *);
  virtual      ~cSecAnHomPlate(void);
          int   Stress(cVector &, cVector &);
          void  CMatrix(cMatrix &);
          void  UpdateState(void);
          void  MbMatrix(cMatrix &);
};


// -------------------------------------------------------------------------
// Definition of the Homogeneous Shell Section Analysis class:
//
class cSecAnHomShell : public cSecAnalysis
{
 protected:
  cConstModel *Cmod;

 public:
                cSecAnHomShell(cElement *, cIntPoint *);
  virtual      ~cSecAnHomShell(void);
          int   Stress(cVector &, cVector &);
          int   Stress(double, cVector &, cVector &, cVector &);
          void  CMatrix(cMatrix &);
          void  UpdateState(void);
          void  MbMatrix(cMatrix &);
};


// -------------------------------------------------------------------------
// Definition of the Homogeneous degenerated Shell Section Analysis class:
//
class cSecAnHomDegShell : public cSecAnalysis
{
 protected:
  cConstModel *Cmod;
  cIntPoint   *Ipt;

 public:
                cSecAnHomDegShell(cElement *, cIntPoint *);
  virtual      ~cSecAnHomDegShell(void);
          int   Stress(cVector &, cVector &);
          void  CMatrix(cMatrix &);
          void  UpdateState(void);
          void  MbMatrix(cMatrix &);
};


// -------------------------------------------------------------------------
// Definition of the FGM 2D Section Analysis class:
//
class cSecAnFGM2D : public cSecAnalysis
{
 protected:
  cConstModel *Cmod;

 public:
                cSecAnFGM2D(cElement *, cIntPoint *);
  virtual      ~cSecAnFGM2D(void);
          int   Stress(cVector &, cVector &);
          void  CMatrix(cMatrix &);
          void  UpdateState(void);
};


// -------------------------------------------------------------------------
// Definition of the FGM 3D Section Analysis class:
//
class cSecAnFGM3D : public cSecAnalysis
{
 protected:
  cConstModel *Cmod;

 public:
                cSecAnFGM3D(cElement *, cIntPoint *);
  virtual      ~cSecAnFGM3D(void);
          int   Stress(cVector &, cVector &);
          void  CMatrix(cMatrix &);
          void  UpdateState(void);
};


// -------------------------------------------------------------------------
// Definition of the FGM Shell Section Analysis class:
//
class cSecAnFGMShell : public cSecAnalysis
{
 protected:
  cAnModel     *Anm;
  cConstModel **Cmod;
  int           NumSecPnt;  // Number of section integration points
  cIntPoint    *IntSecPnt;  // Section integration points

          void  CalcCz(double, cMatrix &, cMatrix &);

 public:
                cSecAnFGMShell(cElement *, cIntPoint *);
  virtual      ~cSecAnFGMShell(void);
          int   Stress(cVector &, cVector &);
          int   Stress(double, cVector &, cVector &, cVector &);
          void  CMatrix(cMatrix &);
          void  UpdateState(void);
          void  MbMatrix(cMatrix &);
};


// -------------------------------------------------------------------------
// Definition of the HSDT Plate Section Analysis class:
//
class cSecAnHSDTPlate : public cSecAnalysis
{
 protected:
  cAnModel     *Anm;
  cConstModel **Cmod;
  int           NumSecPnt;  // Number of section integration points
  cIntPoint    *IntSecPnt;  // Section integration points

          void  CalcCz(double , double, cMatrix &, cMatrix &);

 public:
                cSecAnHSDTPlate(cElement *, cIntPoint *);
  virtual      ~cSecAnHSDTPlate(void);
          int   Stress(cVector &, cVector &);
          void  CMatrix(cMatrix &);
          void  UpdateState(void);
          void  MbMatrix(cMatrix &);
};


// -------------------------------------------------------------------------
// Definition of the HSDT FGM Plate Section Analysis class:
//
class cSecAnHSDTFGMPlate : public cSecAnalysis
{
 protected:
  cAnModel     *Anm;
  cConstModel **Cmod;
  int           NumSecPnt;  // Number of section integration points
  cIntPoint    *IntSecPnt;  // Section integration points

          void  CalcCz(double, double, cMatrix &, cMatrix &);

 public:
                cSecAnHSDTFGMPlate(cElement *, cIntPoint *);
  virtual      ~cSecAnHSDTFGMPlate(void);
          int   Stress(cVector &, cVector &);
          void  CMatrix(cMatrix &);
          void  UpdateState(void);
          void  MbMatrix(cMatrix &);
};


// -------------------------------------------------------------------------
// Definition of the General Shell Section Analysis class:
//
class cSecAnGenShell : public cSecAnalysis
{
 public:
                cSecAnGenShell(cElement *, cIntPoint *);
  virtual      ~cSecAnGenShell(void);
          int   Stress(cVector &, cVector &);
          void  CMatrix(cMatrix &);
          void  MbMatrix(cMatrix &);
};


// -------------------------------------------------------------------------
// Definition of the Laminated Shell Section Analysis class:
//
class cSecAnLaminatedShell : public cSecAnalysis
{
 protected:
  cAnModel     *Anm;
  cConstModel **Cmod;

          void  CalcCz(double, cMatrix &, cMatrix &);

 public:
                cSecAnLaminatedShell(cElement *, cIntPoint *);
  virtual      ~cSecAnLaminatedShell(void);
          int   Stress(cVector &, cVector &);
          int   Stress(double, cVector &, cVector &, cVector &);
          void  CMatrix(cMatrix &);
          void  UpdateState(void);
          void  MbMatrix(cMatrix &);
};


// -------------------------------------------------------------------------
// Definition of the Laminated Solid Section Analysis class:
//
class cSecAnLaminatedSolid : public cSecAnalysis
{
 protected:
  cAnModel    *Anm;
  cIntPoint   *Ipt;
  sLamina     *Ply;
  cConstModel *Cmod;

 public:
                cSecAnLaminatedSolid(cElement *, cIntPoint *, sLamina *);
  virtual      ~cSecAnLaminatedSolid(void);
          void  LaminaSys(cVector &, cVector &, cVector &);
          int   Stress(cVector &, cVector &);
          void  CMatrix(cMatrix &);
          void  UpdateState(void);
};


// -------------------------------------------------------------------------
// Definition of the Plane Frame Rect. Section Analysis (Fiber Method) class:
//
class cSecAnPlFrameRectFiber : public cSecAnalysis
{
 protected:
  int           NumFiber;  // Number of fibers
  cConstModel **Cmod;      // Constitutive model of each fiber

 public:
                cSecAnPlFrameRectFiber(cElement *, cIntPoint *);
  virtual      ~cSecAnPlFrameRectFiber(void);
          int   Stress(cVector &, cVector &);
          void  CMatrix(cMatrix &);
          void  UpdateState(void);
};


// -------------------------------------------------------------------------
// Definition of the Plane Frame Rect. Section Analysis (Gauss Method) class:
//
class cSecAnPlFrameRectGauss : public cSecAnalysis
{
 protected:
  int           NumSecPnt;  // Number of integration points
  cIntPoint    *IntSecPnt;  // Integration points
  cConstModel **Cmod;       // Constitutive model of each fiber

 public:
                cSecAnPlFrameRectGauss(cElement *, cIntPoint *);
  virtual      ~cSecAnPlFrameRectGauss(void);
          int   Stress(cVector &, cVector &);
          void  CMatrix(cMatrix &);
          void  UpdateState(void);
};


// -------------------------------------------------------------------------
// Definition of the Plane Frame RC Rectangular Section Analysis class:
//
class cSecAnPlFrameRCRect : public cSecAnalysis
{
 protected:
  int           NumFiber;  // Number of fibers
  cConstModel **ModFiber;  // Constitutive model of each fiber
  cConstModel **ModRebar;  // Constitutive model of each rebar
  cConstModel **ModVoid;   // Constitutive model of each rebar void

 public:
                cSecAnPlFrameRCRect(cElement *, cIntPoint *);
  virtual      ~cSecAnPlFrameRCRect(void);
          int   Stress(cVector &, cVector &);
          void  CMatrix(cMatrix &);
          void  UpdateState(void);
};

#endif
