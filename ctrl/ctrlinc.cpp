// -------------------------------------------------------------------------
// ctrlinc.cpp - implementation of the control class for nonlinear static
//               problems using the pure incremental (Euler) Method.
// -------------------------------------------------------------------------
// Created:      11-Feb-2009     Evandro Parente Junior
//
// Modified:
// -------------------------------------------------------------------------

#include <iostream>
#include <iomanip>

using namespace std;

#include "ctrlinc.h"
#include "vec.h"
#include "sysmat.h"
#include "gblvar.h"

// -------------------------------------------------------------------------
// Public methods:
//

// =========================== cCtrlIncremental ============================

cCtrlIncremental :: cCtrlIncremental(void)
{
  Type = INCREMENTAL;
}

// ========================== ~cCtrlIncremental ============================

cCtrlIncremental :: ~cCtrlIncremental(void)
{
}

// ================================ Solver =================================

void cCtrlIncremental :: Solver(void)
{
  // Create local arrays.

  cVector q(NumEq);            // Reference load vector
  cVector f(NumEq);            // External load vector
  cVector f0(NumEq);           // Previous load vector
  cVector df(NumEq);           // Incremental load vector
  cVector u(NumEq);            // Displacement vector
  cVector du(NumEq);           // Displacement increment

  // Create the tangent stiffness matrix.

  cSysMatrix *K = cSysMatrix :: CreateMatrix(SysMatType, NumEq, Profile);

  // Compute the reference load vector.

  if (PropLoad)
    ExtLoadVector(q);

  // Initialize some variables.

  TotFactor = 0.0;
  u.Zero( );
  f0.Zero( );
  IntForceVector(f0);   // To include the effects of initial stresses

  // Incremental loop.

  for (CurrStep = 1; CurrStep <= MaxStep; CurrStep++)
  {
    if (Feedback)
    {
      cout << "\n\tIncremental step " << setw(2) << CurrStep; 
      cout << " .....................\n\n";
    }

    // Compute the new external load factor.

    TotFactor += StepFactor;
    
    if (PropLoad)
      f = TotFactor*q;
    else
      ExtLoadVector(TotFactor, f);

    // Evaluate the incremental load vector.

    df = f - f0;

    // Compute the current the stiffness matrix.

    StiffnessMatrix(K);

    // Compute the displacement increment (Solve [K]{du} = {df}).

    K->Solve(df, du);

    // Update the nodal displacements.

    u += du;
    AssignDispl(u);

    // Print the computed results.

    if (CurrStep%PrintStep == 0)
    {
      if (Feedback) cout << "\n\tPrinting computed results ..............." << endl;
      PrintResult( );
    }

    // Update element variables and the previous load vector.

    UpdateState( );
    f0 = f;
  }

  // Adjust the current step.

  CurrStep--;
}

// ======================================================= End of file =====
