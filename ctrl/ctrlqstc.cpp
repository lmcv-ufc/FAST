// -------------------------------------------------------------------------
// ctrlqstc.cpp - implementation of the control class for linear
//                quasi-static problems.
// -------------------------------------------------------------------------
// Created:      13-Jul-2005     Evandro Parente Junior
//
// Modified:     08-Apr-2015     Elias Saraiva Barroso
//               Input/output operations using c++ streams.
// -------------------------------------------------------------------------

#include <iostream>
#include <iomanip>

using namespace std;

#include "ctrlqstc.h"
#include "vec.h"
#include "sysmat.h"
#include "gblvar.h"

// -------------------------------------------------------------------------
// Public methods:
//

// =========================== cCtrlQuasiStatic ============================

cCtrlQuasiStatic :: cCtrlQuasiStatic(void)
{
  Type = QUASISTATIC;
}

// ========================== ~cCtrlQuasiStatic ============================

cCtrlQuasiStatic :: ~cCtrlQuasiStatic(void)
{
}

// ================================ Solver =================================

void cCtrlQuasiStatic :: Solver(void)
{
  // Create the local vectors.

  cVector f(NumEq);            // External load vector
  cVector g(NumEq);            // Internal force vector
  cVector r(NumEq);            // Residual vector
  cVector u(NumEq);            // Displacement vector
  cVector du(NumEq);           // Displacement increment

  // Initialize some variables.

  u.Zero( );
  g.Zero( );
  TotTime = 0.0;

  // Compute the stiffness matrix.

  cSysMatrix *K = new cSymSkylMatrix(NumEq, Profile);
  StiffnessMatrix(K);

  // Incremental loop.

  for (CurrStep = 1; CurrStep <= MaxStep; CurrStep++)
  {
    // Update the total time.

    TotTime += TimeStep;
    if (Feedback)
    {
      cout << "\n\tStep = " << setw(3) << CurrStep << " Time = ";
      cout << scientific << setprecision(3) << TotTime << " ............\n";
    }

    // Compute the external and internal forces.

    ExtLoadVector(TotTime, f);
    IntForceVector(g);

    // Compute the displacement increment (Solve [K]{du} = {r}).

    r = f - g;
    K->Solve(r, du);

    // Update the nodal displacements.

    u += du;
    AssignDispl(u);

    // Compute the new residual (TEMPORARIO P/ ATUALIZAR AS TENSOES) !!!

    IntForceVector(g);
    r = f - g;
    cout << "|r| = " << scientific << r.Length() << "\n";
    cout << resetiosflags(ios::scientific);

    // Print the computed results.

    TotFactor = TotTime;   // Temporario: so´ para impressao !!!
    if (CurrStep%PrintStep == 0)
    {
      if (Feedback) cout << "\n\tPrinting computed results ..............." << endl;
      PrintResult( );
    }

    // Update element variables.

    UpdateState( );
  }

  // Adjust the current step.

  CurrStep--;
}

// ======================================================= End of file =====
