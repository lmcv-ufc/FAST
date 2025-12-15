// ------------------------------------------------------------------------
// ctrlsec.h - file containing the definition of the cCtrlSecant class.
// ------------------------------------------------------------------------
//
// The cCtrlSecant class implements the analysis of nonlinear static
// structures by the secant method. Both proportional and non-proportional
// (time-dependent) loads can be considered.
//
// ------------------------------------------------------------------------

#ifndef _CTRLSEC_H
#define _CTRLSEC_H

#include "ctrl.h"

// ------------------------------------------------------------------------
// Definition of the Secant control class:
//
class cCtrlSecant : public cControl
{
 public:
        cCtrlSecant(void);
       ~cCtrlSecant(void);
  void  Solver(void);
};

#endif
