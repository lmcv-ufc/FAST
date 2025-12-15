// -------------------------------------------------------------------------
// ctrlpath.h - file containing the definition of the cCtrlPath class.
// -------------------------------------------------------------------------
//
// The cCtrlPath class implements the analysis of nonlinear static
// structures using path-following methods, such as displacement control
// and arc-length, using with Newton-Raphson iterations.
//
// Path
// |-- DCM      (Displacement Control)
// |-- INIORTAL (Initially orthogonal arc-length method - Riks)
// |-- UPDORTAL (Updated orthogonal arc-length method - Ramm)
// |-- CILAL    (Cylindrical arc-length - Crisfield)
// |-- CONLINAL (Linearized arc-length)
//
// -------------------------------------------------------------------------
// Pure virtual methods:
// -------------------------------------------------------------------------
//
// virtual double IncLoadFac (cVector &du1 ,cVector &du2, cVector &Du)
//
//  du1  -  [K]{du1} = {q}                                         (in)
//  du2  -  [K]{du2} = {r}                                         (in)
//  Du   -  Total Displacement Increment (step)                    (in)
//
// -------------------------------------------------------------------------

#ifndef _CTRLPATH_H
#define _CTRLPATH_H

#include "ctrl.h"

// -------------------------------------------------------------------------
// Path-Following Algorithm types:
//
typedef enum
{
  DCM,          // Displacement control method
  INIORTAL,     // Initially orthogonal arc-length (Riks)
  UPDORTAL,     // Updated orthogonal arc-length (Ramm)
  CILAL,        // Cylindrical (quadratic) arc-length (Crisfield)
  CONLINAL      // Cons. linearized arc-length
} ePathType;

// -------------------------------------------------------------------------
// Definition of cCtrlPath class:
//
class cCtrlPath : public cControl
{
 protected:
  static  ePathType PathType;   // Path-following algorithm
  static  int       CritPoint;  // Evaluate critical points along the path
  static  int       BranchSwt;  // Perform branch-switching
  static  double    SwtFactor;  // Factor for mode injection (branch-switching)
          double    CSP0,CSP1;  // Current stiffness parameters

  virtual double    IncLoadFac(cVector &,cVector &, cVector &) = 0;
          int       StabSolver(cSysMatrix *, cVector &, cVector &, 
                               cVector &, cVector &, cVector &, cVector &, 
                               cVector &, cVector &, cVector &, cVector &);

 public:
  static  void      ReadStabData(void);
                    cCtrlPath(void);
                   ~cCtrlPath(void);
  virtual void      Solver(void);
};

#endif
