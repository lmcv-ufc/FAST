// -------------------------------------------------------------------------
// ctrldsp.h - file containing the definition of the cCtrlDispl class.
// -------------------------------------------------------------------------
//
// The cCtrlDispl class implements the analysis of nonlinear static
// structures by the Displacement Control Method with Newton-Raphson
// iterations.
//
// -------------------------------------------------------------------------

#ifndef _CTRLDSP_H
#define _CTRLDSP_H

#include "ctrl.h"
#include "ctrlpath.h"

// ------------------------------------------------------------------------
// Definition of the Newton-Raphson control class:
//
class cCtrlDispl : public cCtrlPath
{
 private:
  static int CtrlNode;       // Controlled node
  static int CtrlDir;        // Controlled direction (1 to 6)
  static int CtrlDof;        // Controlled degree of freedom

  static void    GetCtrlDof(void);

 public:
  static void    ReadData(void);

                 cCtrlDispl(void);
                ~cCtrlDispl(void);
         double  IncLoadFac (cVector &,cVector &, cVector &);
};

#endif
