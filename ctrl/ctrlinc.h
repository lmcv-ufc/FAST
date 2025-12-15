// ------------------------------------------------------------------------
// ctrlinc.h - file containing the definition of the cCtrlIncremental
//             class.
// ------------------------------------------------------------------------
//
// The cCtrlIncremental class implements the analysis of nonlinear static
// structures by the pure incremental (Euler) method. Both proportional and
// non-proportional (time-dependent) loads can be considered.
//
// ------------------------------------------------------------------------

#ifndef _CTRLINC_H
#define _CTRLINC_H

#include "ctrl.h"

// ------------------------------------------------------------------------
// Definition of the Incremental control class:
//
class cCtrlIncremental : public cControl
{
 public:
        cCtrlIncremental(void);
       ~cCtrlIncremental(void);
  void  Solver(void);
};

#endif
