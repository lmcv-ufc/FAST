// -------------------------------------------------------------------------
// ctrlqstc.h - definition of the control class for linear Quasi-Static
//              problems.
// -------------------------------------------------------------------------
//
// This algorithm should be used in the analisys of linear time-dependent
// problems when the inertia effects can be neglected. A typical application
// is in the analysis of structures with linear viscoelastic behavior.
//
// -------------------------------------------------------------------------

#ifndef _CTRLQSTC_H
#define _CTRLQSTC_H

#include "ctrl.h"

// -------------------------------------------------------------------------
// Definition of linear quasi-static control class:
//
class cCtrlQuasiStatic : public cControl
{
 public:
        cCtrlQuasiStatic(void);
       ~cCtrlQuasiStatic(void);
  void  Solver(void);
};

#endif
