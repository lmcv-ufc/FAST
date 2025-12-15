// -------------------------------------------------------------------------
// elmfrm3dcr.h - file containing the definition of co-rotational 3D frame
//                element classes.
// -------------------------------------------------------------------------
//
// The cFrame3DEICR class implements the Element Independent Co-Rotational
// (EICR) approach for nonlinear analysis of 3D beams with large displace-
// ments and rotations. This is a abstract class derived from cFrame3D and
// has 6 dofs (u, v, w, rx, ry, rz). Different elements can be derived from 
// this class provided that the approapriate methods for computation of 
// local stiffness and internal forces are defined. The class hierarchy is:
// 
// Frame3DEICR
// |-- Frame3DCR
// |-- Frame3DCRTL
//
// The cFrame3DCR class implements the local 3D linear frame element
// based on the Classical Beam Theory (Euler-Bernoulli).
//
// The cFrame3DCRTL class implements the local 3D frame element based 
// on the Classical Beam Theory (Euler-Bernoulli) with Green-Lagrange
// (aka "shallow arch") strains.
//
// -------------------------------------------------------------------------
// [1] Mororo, L. A. T. "Analise nao linear geometrica de vigas laminadas de 
//     parede fina", Dissertacao de Mestrado, Universidade Federal do Ceara,
//     2013.
// -------------------------------------------------------------------------

#ifndef _ELMFRM3DCR_H
#define _ELMFRM3DCR_H

#include "element.h"
#include "elmfrm3d.h"
#include "secbar.h"

// -------------------------------------------------------------------------
// Forward declarations:
//
class cVector;
class cMatrix;

// -------------------------------------------------------------------------
// Definition of the EICR 3D frame element class:
//
class cFrame3DEICR : public cFrame3D
{
 protected:
  cVector delta_d;   // Vector of ??
  cMatrix R_1;       // Tensor of curvatures attached at node 1.
  cMatrix R_2;       // Tensor of curvatures attached at node 2.
  
          void    CalcL0(sNodeCoord *, double *);
          void    CalcL(sNodeCoord *, cVector &, double *);
          void    Spurrier(cMatrix &, cVector &);
          void    T0Mat(sNodeCoord *, cMatrix &);
          void    TrMat(sNodeCoord *, cVector &, cMatrix &);
          void    SMat(cVector &, cMatrix &);
          void    H1Mat(cVector &, cVector &, cMatrix &);
          void    H2Mat(cVector &, cVector &, cVector &, cVector &, cMatrix &);
          void    RMat(cVector &, cMatrix &);
          void    AlfaMat(cVector &, cMatrix &);
          void    OmegaMat(cVector &, cVector &, cMatrix &);
          void    HMat(cVector &, cVector &, cMatrix &);
          void    OMat(cVector &, cVector &, cVector &, cVector &, cMatrix &);
          void    YMat(cVector &, cVector &, double, double, cMatrix &);
          void    PMat(cMatrix &, cMatrix &, double, double *, cMatrix &, 
                       cMatrix &, cMatrix &);
          void    GetTrnMat(cMatrix &, cMatrix &);
          void    CorNodalDisp(double, double, cMatrix &, cMatrix &, cVector &);
          void    EqvForces(cVector &) { }
  
 public: 
                  cFrame3DEICR(int, cSection *);
  virtual        ~cFrame3DEICR(void);
          int     IntForce(cVector &);
          void    StiffMat(cMatrix &);
          void    MassMat(cMatrix &){ }
          void    GeomStiff(cMatrix &){ }
          void    NodalStress(cMatrix &);
          void    UpdateItera(cVector &);
  virtual void    LocIntForce(cVector &) = 0;
  virtual void    LocStiffMat(cMatrix &) = 0;
};


// -------------------------------------------------------------------------
// Definition of the co-rotational 3D frame element based on Classical 
// Beam Theory (Euler-Bernoulli) class:
//
class cFrame3DCR : public cFrame3DEICR
{
 protected:
          void    ConvLocStiff(double, cMatrix &);
          void    CoupLocStiff(double, cMatrix &);

 public:
                  cFrame3DCR(int, cSection *);
  virtual        ~cFrame3DCR(void);
          void    LocIntForce(cVector &);
          void    LocStiffMat(cMatrix &);
};


// ------------------------------------------------------------------------
// Definition of the Corotational 3D frame element based on CB Theory
// with local Total Lagrangian strains class:
//
class cFrame3DCRTL : public cFrame3DEICR
{
 protected:
          void    AMat(double, double, double, cMatrix &);
          void    AvMat(double, cMatrix &);
          void    AwMat(double, cMatrix &);
          void    AthetaMat(double, double, double, cMatrix &);
          void    B0Mat(double, cMatrix &);
          void    BLvMat(double, cVector &, cMatrix &);
          void    BLwMat(double, cVector &, cMatrix &);
          void    BLthetaMat(double, double, double, cVector &, cMatrix &);
          void    BLMat(double, double, double, cVector &, cMatrix &);
          void    ConvStiffMat(double, cVector &, cMatrix &);
          void    CoupStiffMat(double, cVector &, cMatrix &);
          void    ConvIntForce(double, cVector &, cVector &);
          void    CoupIntForce(double, cVector &, cVector &);
          double  ConvNormal(double, cVector &);
          double  CoupNormal(double, cVector &);
  
 public:
                  cFrame3DCRTL(int, cSection *);
  virtual        ~cFrame3DCRTL(void);
          void    LocIntForce(cVector &);
          void    LocStiffMat(cMatrix &);
};

#endif
