// ------------------------------------------------------------------------
// ctrllins.h - definition of the control class for linear static problems.
// ------------------------------------------------------------------------
//
// The cLinStat class implements the steps required to perform the static
// analysis of linear  structures.
//
// ------------------------------------------------------------------------

#ifndef _CTRLLINS_H
#define _CTRLLINS_H

#include "ctrl.h"

// ------------------------------------------------------------------------
// Definition of linear static control class:
//
class cLinStat : public cControl
{
 protected:
  void  SolverMicro(void);  // Temporary - Evandro

 public:
        cLinStat(void);
       ~cLinStat(void);
  void  Solver(void);
};

#endif
