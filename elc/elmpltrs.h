// ------------------------------------------------------------------------
// elmpltrs.h - file containing the definition of plane truss classes.
// ------------------------------------------------------------------------
//
// The cPlTruss class implements the linear plane truss element and is the
// base class for the geometrically nonlinear elements. The plane truss
// element has two degrees of freedom (u,v) per node, and only one stress
// component, the normal force (FORCE_X).
//
// ------------------------------------------------------------------------
//
// The cPlTrussTL class implements the geometrically nonlinear plane truss
// element based on the total Lagrangian approach. Therefore, the element
// response is the normal force corresponding to the Piola-Kirchoff II
// computed from the Green-Lagrange strains.
//
// ------------------------------------------------------------------------
//
// The cPlTrussCR class implements the geometrically nonlinear plane truss
// element based on the corotational approach. Therefore, the element
// response is the normal force corresponding to the stresses associated
// with the engineering strains.
//
// ------------------------------------------------------------------------

#ifndef _ELMPLTRS_H
#define _ELMPLTRS_H

#include "element.h"

// ------------------------------------------------------------------------
// Forward declarations:
//
class cSecBar;
class cConstModel;

// ------------------------------------------------------------------------
// Definition of the plane truss element class:
//
class cPlTruss : public cElement
{
 protected:
  cSecBar     *Section;
  cConstModel *ConstMod;

          void  B0Matrix(sNodeCoord *, double *, cVector &);

 public:
                cPlTruss(int);
  virtual      ~cPlTruss(void);
          int   GetNumDofNode(void) { return 2; }
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
// Definition of the total Lagrangian plane truss element class:
//
class cPlTrussTL : public cPlTruss
{
 protected:
  void  BlMatrix(sNodeCoord *, cVector &, double *, cVector &);

 public:
        cPlTrussTL(int);
       ~cPlTrussTL(void);
  int   IntForce(cVector &);
  void  StiffMat(cMatrix &);
  void  NodalStress(cMatrix &);
};


// ------------------------------------------------------------------------
// Definition of the corotational plane truss element class:
//
class cPlTrussCR : public cPlTruss
{
 protected:
  void  RVector(sNodeCoord *, cVector &, double *, double *, cVector &);

 public:
        cPlTrussCR(int);
       ~cPlTrussCR(void);
  int   IntForce(cVector &);
  void  StiffMat(cMatrix &);
  void  NodalStress(cMatrix &);
};

#endif
