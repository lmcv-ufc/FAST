// -------------------------------------------------------------------------
// ctrlsec.cpp - implementation of the control class for nonlinear static
//               problems by the secant method (fixed point iterations).
// -------------------------------------------------------------------------
// Created:      19-Feb-2009     Evandro Parente Junior
//
// Modified:     08-Apr-2015     Elias Saraiva Barroso
//               Input/output operations using c++ streams.
// -------------------------------------------------------------------------

#include <iostream>
#include <iomanip>

using namespace std;

#include "ctrlsec.h"
#include "vec.h"
#include "sysmat.h"
#include "gblvar.h"

// -------------------------------------------------------------------------
// Public methods:
//

// ============================= cCtrlSecant ===============================

cCtrlSecant :: cCtrlSecant(void)
{
  Type = SECANT;
}

// ============================= ~cCtrlSecant ==============================

cCtrlSecant :: ~cCtrlSecant(void)
{
}

// ================================ Solver =================================

void cCtrlSecant :: Solver(void)
{
  // Create local arrays.

  cVector q(NumEq);            // Reference load vector
  cVector f(NumEq);            // External load vector
  cVector f0(NumEq);           // Initial load vector
  cVector g(NumEq);            // Internal force vector
  cVector u(NumEq);            // Current displacement vector
  cVector u0(NumEq);           // Previous displacement vector
  cVector du(NumEq);           // Displacement increment

  // Create the tangent stiffness matrix.

  cSysMatrix *K = cSysMatrix :: CreateMatrix(SysMatType, NumEq, Profile);

  // Compute the reference load vector.

  if (PropLoad)
    ExtLoadVector(q);

  // Initialize some variables.

  TotFactor = 0.0;
  u.Zero( );
  u0.Zero( );
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

    TotFactor += StepFactor;
    if (PropLoad)
      f = TotFactor*q;
    else
      ExtLoadVector(TotFactor, f);

    f -= f0;              // Evaluate the effective external force

    // Iterative loop.

    int conv = 0;
    for (CurrIter = 1; CurrIter <= MaxIter; CurrIter++)
    {
      // Compute the current the stiffness matrix.

      StiffnessMatrix(K);

      // Compute the new displacement vector (Solve [K]{u} = {f}).

      K->Solve(f, u);

      // Update the nodal displacements and compute the internal force.

      AssignDispl(u);
      if (!IntForceVector(g)) break;   // Only to compute the new stresses

      // Compute the displacement increment and check convergence.

      du = u - u0;
      u0 = u;
      conv = Convergence(u.Length( ), du.Length( ));
      if (conv) break;
    }

    // Print the computed results.

    if (conv)
    {
      if (CurrStep%PrintStep == 0)
      {
        if (Feedback) cout << "\n\tPrinting computed results ..............." << endl;
        PrintResult( );
      }
    }
    else
    {
      cout << "Convergence not achived!!!\n";
      break;
    }

    // Update element variables.

    UpdateState( );
  }

  // Adjust the current step.

  CurrStep--;
}

// ======================================================= End of file =====
