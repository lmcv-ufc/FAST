// ------------------------------------------------------------------------
// ctrlnewmark.h - definition of the control class for dynamic analysis
// using the Newmark´s algorithm.
// ------------------------------------------------------------------------
// There are 3 different versions of the algorithm:
//  CtrlNewmark - for linear problems based on total displacements. 
//  CtrlIncNwmk - for linear problems based on incremental displacements. 
//  cCtrlNonlinNewmark - for nonlinear problems.
//
// Obs:
//  - The 3 versions give the same results for linear problems.
//  - cCtrlNonlinNewmark can be used for linear dynamic analysis with
//    viscoelastic materials.
//  - The damping matrix is evaluated using the Rayleigh damping model
//    [C] = alpha*[M] + beta*[K].
//  - The nonlinear algorithm considers that the mass matrix [M] is constant
//    and the stiffness matrix [K] and damping matrix [C] are updated at
//    each Newton-Raphson equilibrium iteration.
// ------------------------------------------------------------------------
// References:
// [1] Cook et al., "Concepts and Applications of Finite Element Analysis",
//     4th ed.,  Element Procedures", John Wiley & Sons, 2002.
// [2] K. J. Bathe, "Finite Element Procedures", Prentice Hall, 1996.
// ------------------------------------------------------------------------

#ifndef _CTRLNEWMARK_H
#define _CTRLNEWMARK_H

#include "ctrl.h"

// ------------------------------------------------------------------------
// Definition of linear dynamic Newmark control class:
//
class cCtrlNewmark : public cControl
{
 protected:
  static double Gamma;   // Newmark´s parameter
  static double Beta;    // Newmark´s parameter
  static bool   Damping; // Flag for Rayleigh damping
  static double Omega1;  // Rayleigh damping 1st natural frequency
  static double Omega2;  // Rayleigh damping 2nd natural frequency
  static double Xi1;     // Rayleigh damping 1st critical ratio
  static double Xi2;     // Rayleigh damping 2nd critical ratio

 protected:
          void  DampingMatrix(cSysMatrix *, cSysMatrix *, cSysMatrix *);
	  void  PrintVelAcc(cVector &, cVector &); 

 public:
  static  void  ReadNewmarkParam(void);
  static  void  ReadRayleighDamp(void);

                cCtrlNewmark(void);
               ~cCtrlNewmark(void);
  virtual void  Solver(void);
};


// -------------------------------------------------------------------------
// Definition of linear dynamic Newmark (incremental) control class:
//
class cCtrlIncNmk : public cCtrlNewmark
{
 public:
        cCtrlIncNmk(void);
       ~cCtrlIncNmk(void);
  void  Solver(void);
};


// -------------------------------------------------------------------------
// Definition of nonlinear dynamic Newmark control class:
//
class cCtrlNonlinNewmark : public cCtrlNewmark
{
 public:
           cCtrlNonlinNewmark(void);
          ~cCtrlNonlinNewmark(void);
  void     Solver(void);
};
#endif
