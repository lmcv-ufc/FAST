// ------------------------------------------------------------------------
// elmtrs3d.h - file containing the definition of 3D truss classes.
// ------------------------------------------------------------------------
//
// The cTruss3D class implements the linear 3D truss element and is the
// base class for the geometrically nonlinear elements. The 3D truss
// element has three degrees of freedom (u,v,w) per node, and only one
// stress component, the normal force (FORCE_X).
//
// ------------------------------------------------------------------------
//
// The cTruss3DTL class implements the geometrically nonlinear 3D truss
// element based on the total Lagrangian approach. Therefore, the element
// response is the normal force corresponding to the Piola-Kirchoff II
// computed from the Green-Lagrange strains.
//
// ------------------------------------------------------------------------
//
// The cTruss3DCR class implements the geometrically nonlinear 3D truss
// element based on the corotational approach. Therefore, the element
// response is the normal force corresponding to the stresses associated
// with the engineering strains.
//
// ------------------------------------------------------------------------

#ifndef _ELMTRS3D_H
#define _ELMTRS3D_H

#include "element.h"

// ------------------------------------------------------------------------
// Forward declarations:
//
class cSecBar;
class cConstModel;

// ------------------------------------------------------------------------
// Definition of the 3D truss element class:
//
class cTruss3D : public cElement
{
 protected:
  cSecBar     *Section;
  cConstModel *ConstMod;

          void  B0Matrix(sNodeCoord *, double *, cVector &);

 public:
                cTruss3D(int);
  virtual      ~cTruss3D(void);
          int   GetNumDofNode(void) { return 3; }
          int   GetNumStrCmp(void)  { return 1; }
          void  GetStrLabels(int *);
          void  GetActDir(int *);
  virtual void  Read(void);
  virtual int   IntForce(cVector &);
  virtual void  StiffMat(cMatrix &);
  virtual void  MassMat(cMatrix &);
  virtual void  GeomStiff(cMatrix &);
  virtual void  NodalStress(cMatrix &);
  virtual void  UpdateState(void);
};


// ------------------------------------------------------------------------
// Definition of the Total Lagrangian 3D truss element class:
//
class cTruss3DTL : public cTruss3D
{
 protected:
  void  BlMatrix(sNodeCoord *, cVector &, double *, cVector &);

 public:
        cTruss3DTL(int);
       ~cTruss3DTL(void);
  int   IntForce(cVector &);
  void  StiffMat(cMatrix &);
  void  NodalStress(cMatrix &);
};


// ------------------------------------------------------------------------
// Definition of the corotational 3D truss element class:
//
class cTruss3DCR : public cTruss3D
{
 protected:
  void  RVector(sNodeCoord *, cVector &, double *, double *, cVector &);

 public:
        cTruss3DCR(int);
       ~cTruss3DCR(void);
  int   IntForce(cVector &);
  void  StiffMat(cMatrix &);
  void  NodalStress(cMatrix &);
};

#endif
