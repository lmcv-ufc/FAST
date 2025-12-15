// ------------------------------------------------------------------------
// ctrlmodal.h - definition of the control class for modal analysis with
//               computation of vibration modes and natural frequecies.
// ------------------------------------------------------------------------
// References:
// [1] K. J. Bathe, "Finite Element Procedures", Prentice Hall, 1996.
// ------------------------------------------------------------------------

#ifndef _CTRLMODAL_H
#define _CTRLMODAL_H

#include "ctrl.h"

// ------------------------------------------------------------------------
// Definition of modal analysis control class:
//
class cCtrlModal : public cControl
{
 protected:
          void  PrintMode(int, double, cVector&);

 public:
                cCtrlModal(void);
               ~cCtrlModal(void);
  virtual void  Solver(void);
};

#endif
