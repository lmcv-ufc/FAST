// ------------------------------------------------------------------------
// ctrlnr.h - file containing the definition of the cCtrlNR class.
// ------------------------------------------------------------------------
//
// The cCtrlNR class implements the analysis of nonlinear static structures
// by the Load Control Method with Newton-Raphson iterations. Both
// proportional and non-proportional (time-dependent) loads can be
// considered.
//
// ------------------------------------------------------------------------

#ifndef _CTRLNR_H
#define _CTRLNR_H

#include "ctrl.h"

// ------------------------------------------------------------------------
// Definition of the Newton-Raphson control class:
//
class cCtrlNR : public cControl
{
 protected:
  void  SolverMicro(void);  // Temporary - Evandro

 public:
        cCtrlNR(void);
       ~cCtrlNR(void);
  void  Solver(void);
};

#endif
