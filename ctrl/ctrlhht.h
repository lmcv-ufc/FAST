// ------------------------------------------------------------------------
// ctrlhht.h - definition of the control classes for dynamic analysis 
//             using the HHT algorithm.
// ------------------------------------------------------------------------
// References:
// [1] K. J. Bathe, "Finite Element Procedures", Prentice Hall, 1996.
// ------------------------------------------------------------------------

#ifndef _CTRLHHT_H
#define _CTRLHHT_H

#include "ctrlnewmark.h"

// ------------------------------------------------------------------------
// Definition of linear dynamic (HHT algorithm) control class:
//
class cCtrlHHT : public cCtrlNewmark
{
 protected:
  static double Alpha;     // HHT parameter

 public:
  static  void  ReadHHTParam(void);

                cCtrlHHT(void);
               ~cCtrlHHT(void);
  virtual void  Solver(void);
};


// -------------------------------------------------------------------------
// Definition of nonlinear dynamic (HHT algorithm) control class:
//
class cCtrlNonlinHHT : public cCtrlHHT
{
 public:
           cCtrlNonlinHHT(void);
          ~cCtrlNonlinHHT(void);
  void     Solver(void);
};

#endif
