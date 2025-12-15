// -------------------------------------------------------------------------
// ctrlnlstab.h - definition of he control class for nonlinear stability
// analysis.
// -------------------------------------------------------------------------
// The implementation is based on the papers:
//
// [1] P. Wriggers, W. Wagner & C. Miehe, "A Quadratically Convergent
//     Procedure for the Calculation of Stability Points in Finite Element
//     Analysis", CMAME, vol 70, pp 329-347, 1988.
//
// [2] P. Wriggers & J. C. Simo, "A General Procedure for Direct
//     Computation of Turning and Bifurcation Points", IJNME, vol 30,
//     pp 155-176, 1990.
// -------------------------------------------------------------------------

#ifndef	_CTRLNLSTAB_H
#define	_CTRLNLSTAB_H

#include "ctrl.h"

// -------------------------------------------------------------------------
// Forward declarations:
//
class cVector;
class cSysMatrix;

// -------------------------------------------------------------------------
// Definition of nonlinear stability analysis control class:
//
class cCtrlNlStab : public cControl
{
 protected:
  static double RelPert;

          void  DirecDeriv (cVector &, cVector &, cVector &, cVector &,
                            cVector &, cVector &, cVector &, cVector &);
          void  PrintMode  (cVector &, cVector &);

 public:
  static  void  ReadData   (void);
                cCtrlNlStab(void);
  virtual      ~cCtrlNlStab(void) { }
  virtual void  Solver     (void);
  virtual void  Solver     (cSysMatrix *, cVector &, cVector &, cVector &,
                            cVector &, cVector &, cVector &, cVector &, 
                            cVector &, int *);
};

#endif
