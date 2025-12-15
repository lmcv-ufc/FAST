// ------------------------------------------------------------------------
// ctrlheat.h - definition of the control class for heat transfer problems.
// ------------------------------------------------------------------------
//
// The cHeatTransf class implements the steps required to perform the heat 
// transfer linear analysis.
//
// ------------------------------------------------------------------------

#ifndef _CTRLHEATCOND_H
#define _CTRLHEATCOND_H

#include "ctrl.h"

// ------------------------------------------------------------------------
// Definition of heat conduction control class:
//
class cHeatTransf : public cControl
{
 public:
        cHeatTransf(void);
       ~cHeatTransf(void);
  void  GlbHeatVector(cVector&);
  void  ConductivityMatrix(cSysMatrix*);
  void  AssignTemp(cVector&);
  void  PrintResult(void);
  void  PrintIntPntFlux(void);
  void  PrintTempErrNorm(void);
  void  PrintEqvMshTemp(void);
  void  Solver(void);
};

#endif
