// ------------------------------------------------------------------------
// ctrlgalpha.h - definition of the control classes for dynamic analysis 
// using the Generalized-Alpha implicit method.
// ------------------------------------------------------------------------
// References:
// [1] Chung, J. Hulbert, G. M. "Improved Numerical Dissipation: The 
//     Generalized-Alpha Method", Journal of Applied Mechanics, 1993.
// [2] Erlicher, S. Bonaventura, L. Bursi, O. S. "The analysis of the 
//     Generalized-a method for non-linear dynamic problems", Computational 
//     Mechanics, 2002.
// [3] Kuhl, D. Crisfield, M. A. "Energy-Conserving and Decaying Algorithms 
//     in Non-Linear Structural Dynamics", Int. J. Numer. Meth. Eng., 1999.
// ------------------------------------------------------------------------
// The implementation of the Generalized-Alpha method is based on the 
// solution of the equilibrium equations:
//   M*a + C*v + g(u) = f(t)
// at the discrete time [1]:
//   t(n+1-Alphaf) = (1 - Alphaf)*t(n+1) +  Alphaf*t(n),
// where the displacements (u), velocities (v), and accelerations (a) are 
// approximated as:
//   u(n+1-Alphaf) = (1 - Alphaf)*u(n+1) +  Alphaf*u(n),
//   v(n+1-Alphaf) = (1 - Alphaf)*v(n+1) +  Alphaf*v(n),
//   a(n+1-Alpham) = (1 - Alpham)*a(n+1) +  Alpham*a(n),
// and the external and internal forces are approximated using the 
// generalized trapezoidal rule [2]:
//   f(n+1-Alphaf) = (1 - Alphaf)*f(n+1) +  Alphaf*f(n) 
//   g(n+1-Alphaf) = (1 - Alphaf)*g(n+1) +  Alphaf*g(n),
// where,
//   n = increment number.
// The displacements and velocities at t(n+1) are approximated as in the
// Newmark method with [1]:
//   gamma = 1/2 - Alpham + Alphaf               Eq. (17)
//   beta  = 1/4*(1.0 - Alpham + Alphaf)^2       Eq. (19)
// These parameters ensure that the method is second-order accurate and the
// high frequency dissipation is maximized [1].
// ------------------------------------------------------------------------

#ifndef _CTRLGENALPHA_H
#define _CTRLGENALPHA_H

#include "ctrlnewmark.h"

// ------------------------------------------------------------------------
// Definition of linear dynamic Generalized-Alpha algorithm control class:
//
class cCtrlGenAlpha : public cCtrlNewmark
{
 protected:
  static double Alpham;     // Alpha-m parameter
  static double Alphaf;     // Alpha-f parameter

 public:
  static  void  ReadGenAlphaParam(void);

                cCtrlGenAlpha(void);
               ~cCtrlGenAlpha(void);
  virtual void  Solver(void);
};


// -------------------------------------------------------------------------
// Definition of nonlinear dynamic Generalized-Alpha algorithm control class:
//
class cCtrlNonlinGenAlpha : public cCtrlGenAlpha
{
 public:
         cCtrlNonlinGenAlpha(void);
        ~cCtrlNonlinGenAlpha(void);
  void   Solver(void);
};

#endif
