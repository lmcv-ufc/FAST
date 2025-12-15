// -------------------------------------------------------------------------
// elmfrm3d.h - file containing the definition of 3D frame element classes.
// -------------------------------------------------------------------------
//
// The cFrame3D class implements the linear 3D frame element based on the
// Euler theory of beams. The  element has 6 dofs (u, v, w, rx, ry, rz)
// per node and six stress components: normal force (FORCE_X), two shear 
// force (FORCE_Y and FORCE_Z), two bending moment (MOMENT_Y and MOMENT_Z),
// and torsion moment (MOMENT_X).
//
// The cFrame3DTL class implements the nonlinear 3D frame element based on
// the Total Lagrangian (TL) approach. This element should be used only in
// the geometrically nonlinear analysis of frames with moderate rotations.
// -------------------------------------------------------------------------
//
// References:
//
// [1] W. McGuire, R. H. Gallgher, R. D. Ziemian, "Matrix Structural
//     Analysis", John Wiley & Sons, Inc., 2 Ed, 2000.
// -------------------------------------------------------------------------

#ifndef _ELMFRM3D_H
#define _ELMFRM3D_H

#include "element.h"
#include "secbar.h"

// -------------------------------------------------------------------------
// Definition of the three-dimensional frame element class:
//
class cFrame3D : public cElement
{
 protected:
  cSecBar *Section;   // Element cross-section

          double  CalcLength(sNodeCoord *);
          void    GetTrnMat(cMatrix &);
          void    LocStiffMat(double, cMatrix &);
          void    EqvForces(cVector &);

 public:
  static  bool    ValidSection(eSecType);

                  cFrame3D(int, cSection *);
  virtual        ~cFrame3D(void);
          int     GetNumDofNode(void) { return 6; }
          int     GetNumStrCmp(void)  { return 6; }
          void    Read(void);
          void    GetActDir(int *);
          void    GetStrLabels(int *);
          int     IntForce(cVector &);
          void    StiffMat(cMatrix &);
          void    MassMat(cMatrix &){ }
          void    GeomStiff(cMatrix &){ }
          void    NodalStress(cMatrix &);
          void    CalcRotMat(cMatrix &);
};


// -------------------------------------------------------------------------
// Definition of the Total Lagrangian three-dimensional frame element class:
//
class cFrame3DTL : public cFrame3D
{
  protected:
           void   LocIntForce(double, double, double, cVector &, cVector &);
           void   AMat(double, double, double, cMatrix &);
           double NormalForce(double, double, double, cVector &);
           void   LocStiffMat(double, double, double, cVector &, cMatrix &);
  
  public:
                  cFrame3DTL(int, cSection *);
   virtual       ~cFrame3DTL(void);
           int    IntForce(cVector &);
           void   StiffMat(cMatrix &);
};

#endif
