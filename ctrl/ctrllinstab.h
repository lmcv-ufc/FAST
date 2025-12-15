// ------------------------------------------------------------------------
// ctrllinstab.h - definition of the control class for linear stability
// analysis (eigenvalue buckling).
// ------------------------------------------------------------------------
// References:
// [1] K. J. Bathe, "Finite Element Procedures", Prentice Hall, 1996.
// ------------------------------------------------------------------------

#ifndef _CTRLLINSTAB_H
#define _CTRLLINSTAB_H

#include "ctrl.h"

// ------------------------------------------------------------------------
// Definition of linear stability analysis control class:
//
class cCtrlLinStab : public cControl
{
 protected:
          void   PrintMode(int, double, cVector &);

 public:
  static  void   ReadEigenAlg(void);


                 cCtrlLinStab(void);
                ~cCtrlLinStab(void);
  virtual void   Solver(void);
};

#endif
