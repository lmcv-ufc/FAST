// ------------------------------------------------------------------------------
// elmgrid.h - file containing the definition of grid classes.
// ------------------------------------------------------------------------------
//
// The cGrid class implements the linear grid element
// based on the Euler theory of beams. The  element has 3 dofs (w,rx,ry)
// per node and three stress components: one shear force
// (FORCE_Z), bending moment (MOMENT_Y), and
// torsion moment (MOMENT_X).
//
// ------------------------------------------------------------------------------

#ifndef _ELMGRID_H
#define _ELMGRID_H

#include "element.h"
#include "secbar.h"

// ------------------------------------------------------------------------------
// Definition of the grid element class:
//
class cGrid : public cElement
{
 protected:
  cSecBar *Section;   // Element cross-section

          double  CalcLength(sNodeCoord *);
          void    GetTrnMat(cMatrix &);
          void    LocStiffMat(double, cMatrix &);
          void    EqvForces(cVector &);

 public:
  static  bool    ValidSection(eSecType);

                  cGrid(int, cSection *);
  virtual        ~cGrid(void);
          int     GetNumDofNode(void) { return 3; }
          int     GetNumStrCmp(void)  { return 3; }
          void    Read(void);
          void    GetActDir(int *);
          void    GetStrLabels(int *);
          int     IntForce(cVector &);
          void    StiffMat(cMatrix &);
          void    MassMat(cMatrix &){ }
          void    GeomStiff(cMatrix &){ }
          void    NodalStress(cMatrix &);
        //  void    CalcRotMat(cMatrix &);   -> não está utilizando nesta
	                                   //  implementação!
};

#endif
