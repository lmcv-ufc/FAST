// ------------------------------------------------------------------------
// elmplfrm.h - file containing the definition of plane frame classes.
// ------------------------------------------------------------------------
//
// The cPlFrame class implements the plane frame element based on the
// Euler-Bernoulli theory of beams. The  element has 3 dofs (u,v,rz) per 
// node and three stress components: normal force (FORCE_X), shear force
// (FORCE_Y), and bending moment (MOMENT_Z). The class hierarchy is:
//
// PlFrame
// |-- PlFrameNI
// |-- PlFrameCR
// |   |-- PlFrameCR1 
// |   |   |-- PlFrameCR1NI
// |   |-- PlFrameCR2 
// |   |   |-- PlFrameCR2NI
// |-- PlFrameTL1
// |   |-- PlFrameTL2
//
// ------------------------------------------------------------------------
//
// PlFrame corresponds to the linear plane frame element. It uses the
// standard local and global systems with 6 dof of the Matrix Structural
// Analysis (Direct Stiffness Method). The stiffness, mass and geometric
// geometric stiffness matrices are computed using analytical expressions 
// (e.g [1]) in the local system and transformed to the global system
// using the standard 6x6 transformation matrix [T].
// ------------------------------------------------------------------------
//
// PlFrameNI is the Numerically Integrated version of PlFrame. This element
// should be used for materially nonlinear analysis of plane frames with 
// small displacements.
// ------------------------------------------------------------------------
//
// PlFrameCR defines the general behavior of plane frame elements based on
// the co-rotational formulation for large displacement analysis [2]. It
// uses two different coordinate systems:
//   1) A local (natural or corotational) system with only 3 dofs (1 axial
//      displacement and 2 local rotations) which moves with the element.
//      Therefore, there are no rigid body motions in this system.
//   2) A fixed global system with 6 dofs parallel to the global dofs.
//
// The transformation between these two coordinate systems is performed
// using a displacement-dependent rectangular (3x6) [T] matrix. This trans-
// formation depends only on the local and global dofs and is handled by the
// base PlFrameCR class. On the other hand, different strain definitions 
// can be used in the local system leading to different expressions for the
// internal force and tangent stiffness matrices. Thus, these aspects are 
// handled by the PlFrameCR subclasses.
// ------------------------------------------------------------------------
//
// References:
//
// [1] Cook et al, 2002.
// [2] ????
// ------------------------------------------------------------------------

#ifndef _ELMPLFRM_H
#define _ELMPLFRM_H

#include "element.h"
#include "section.h"

// ------------------------------------------------------------------------
// Forward declarations:
//
class cIntPoint;
class cSecAnalysis;

// ------------------------------------------------------------------------
// Definition of the plane frame element class:
//
class cPlFrame : public cElement
{
 protected:
          void  CalcConfig(sNodeCoord *, double *, double *, double *);
          void  GetTrnMat(double, double, cMatrix &);
          void  LocStiffMat(double, cMatrix &);
          void  EqvForces(cVector &);
 public:
  static  bool  ValidSection(eSecType);

                cPlFrame(int, cSection *);
  virtual      ~cPlFrame(void);
          int   GetNumDofNode(void) { return 3; }
          int   GetNumStrCmp(void)  { return 3; }
          void  Read(void);
          void  GetActDir(int *);
          void  GetStrLabels(int *);
  virtual int   IntForce(cVector &);
  virtual void  StiffMat(cMatrix &);
  virtual void  MassMat(cMatrix &);
  virtual void  GeomStiff(cMatrix &);
  virtual void  NodalStress(cMatrix &);
  virtual void  IntPntStress(cMatrix &, cMatrix &);
};

// -------------------------------------------------------------------------
// Definition of the corotational plane frame element class:
//
class cPlFrameNI : public cPlFrame
{
 protected:
  int            NumIntPnt;  // Number of int. points in x-axis
  cIntPoint     *IntPnt;     // Array of int. points in x-axis
  cSecAnalysis **SecAn;      // Section for integration

          void  BMatrix(double, double, cMatrix &);
  
 public:
                cPlFrameNI(int, cSection *);
  virtual      ~cPlFrameNI(void);
          int   GetNumIntPnt(void) { return 0; /*NumIntPnt;*/ }
          int   IntForce(cVector &);
          void  StiffMat(cMatrix &);
          void  UpdateState(void);
          void  IntPntStress(cMatrix &, cMatrix &);
};

// -------------------------------------------------------------------------
// Definition of the corotational plane frame element class:
//
class cPlFrameCR : public cPlFrame
{
 protected:
  double RigRot;           // Rigid body rotation
  
          double CalcLength(sNodeCoord *);
          double GetCurrRot(sNodeCoord *, cVector &);
          void   UpdateGeom(sNodeCoord *, cVector &);
          void   GetTrnMat(double, double, double, cMatrix &);
          void   GetTrnVec(double, double, cVector &, cVector &);
          void   EqvForces(cVector &);
  
 public:
                 cPlFrameCR(int, cSection *);
  virtual       ~cPlFrameCR(void);
          int    IntForce(cVector &);
          void   StiffMat(cMatrix &);
          void   NodalStress(cMatrix &);
          void   UpdateState(void);
  virtual void   UpdateSection(void) { }
  virtual void   LocIntForce(cVector &) = 0;
  virtual void   LocStiffMat(cMatrix &) = 0;
};

// -------------------------------------------------------------------------
// Definition of the standard corotation plane frame element class:
//
class cPlFrameCR1 : public cPlFrameCR
{
 public:
                cPlFrameCR1(int, cSection *);
  virtual      ~cPlFrameCR1(void);
          void  LocIntForce(cVector &);
          void  LocStiffMat(cMatrix &);
};

// -------------------------------------------------------------------------
// Definition of the numerically integrated standard corotation plane frame 
// element class:
//
class cPlFrameCR1NI : public cPlFrameCR1
{
 protected:
  int            NumIntPnt;  // Number of int. points in x-axis
  cIntPoint     *IntPnt;     // Array of int. points in x-axis
  cSecAnalysis **SecAn;      // Section for integration

          void  BMatrix(double, double, cMatrix &);
  
 public:
                cPlFrameCR1NI(int, cSection *);
  virtual      ~cPlFrameCR1NI(void);
          int   GetNumIntPnt(void) { return 0; /*NumIntPnt;*/ }
          void  LocIntForce(cVector &);
          void  LocStiffMat(cMatrix &);
          void  UpdateSection(void);
};

// -------------------------------------------------------------------------
// Definition of the corotational plane frame element with shallow arch
// strain class:
//
class cPlFrameCR2 : public cPlFrameCR
{
 protected:
          void  BimMatrix(double, double, double, cMatrix &);
  
 public:
                cPlFrameCR2(int, cSection *);
  virtual      ~cPlFrameCR2(void);
          void  LocIntForce(cVector &);
          void  LocStiffMat(cMatrix &);
};

// -------------------------------------------------------------------------
// Definition of the numerically integrated corotation plane frame element 
// with shallow arch strain class:
//
class cPlFrameCR2NI : public cPlFrameCR2
{
 protected:
  int            NumIntPnt;  // Number of int. points in x-axis
  cIntPoint     *IntPnt;     // Array of int. points in x-axis
  cSecAnalysis **SecAn;      // Section for integration

          void  BbarMatrix(double, double, double, double, cMatrix &);
  
 public:
                cPlFrameCR2NI(int, cSection *);
  virtual      ~cPlFrameCR2NI(void);
          int   GetNumIntPnt(void) { return 0; /*NumIntPnt;*/ }
          void  LocIntForce(cVector &);
          void  LocStiffMat(cMatrix &);
          void  UpdateSection(void);
};

// -------------------------------------------------------------------------
// Definition of the Total Lagrangian plane frame TL1 element class:
//
class cPlFrameTL1 : public cPlFrame
{
 protected:
          void  AMatrix(double, cMatrix &);
          void  B0mMatrix(double, cMatrix &);
          void  Blm1Matrix(double, cVector &, cMatrix &);
          void  Blm2Matrix(double, cVector &, cMatrix &);
          void  KbeMatrix(double, cMatrix &);
  
 public:
                cPlFrameTL1(int, cSection *);
  virtual      ~cPlFrameTL1(void);
  virtual int   IntForce(cVector &);
  virtual void  StiffMat(cMatrix &);
};

// -------------------------------------------------------------------------
// Definition of the Total Lagrangian plane frame TL2 (shallow arc)
// element class:
//
class cPlFrameTL2 : public cPlFrameTL1
{
 public:
                  cPlFrameTL2(int, cSection *);
   virtual       ~cPlFrameTL2(void);
   virtual int   IntForce(cVector &);
   virtual void  StiffMat(cMatrix &);
};

#endif
