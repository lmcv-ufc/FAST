// -------------------------------------------------------------------------
// ctrllins.cpp - implementation of the control class for lineat static
//               problems.
// -------------------------------------------------------------------------
// Created:      16-Nov-2000     Evandro Parente Junior
//
// Modified:     03-Oct-2001     Evandro Parente Junior
//               Use of vector/matrix classes.
//
// Modified:     23-Oct-2024     Elias Saraiva Barroso
//               Implementation of MPCs.
// -------------------------------------------------------------------------

#include <iostream>
#include <iomanip>

using namespace std;

#include "ctrllins.h"
#include "vec.h"
#include "sysmat.h"
#include "gblvar.h"
#include "node.h"
#include "element.h"
#include "section.h"

// -------------------------------------------------------------------------
// Public methods:
//

// =============================== cLinStat ================================

cLinStat :: cLinStat(void)
{
  Type = LINEAR_STATIC;
}

// =============================== ~cLinStat ===============================

cLinStat :: ~cLinStat(void)
{
}

// ================================ Solver =================================

void cLinStat :: Solver(void)
{
  // Temporary - Evandro (09-Oct-2025).

  if (MicroModel)
  {
    SolverMicro( );
    return;
  }

  if (Feedback) cout << "\n\tComputing the external load vector ......" << endl;

  cVector f(NumEq);
  ExtLoadVector(f);

  if (Feedback) cout << "\n\tAssembling the stiffness matrix ........." << endl;

  cSysMatrix *K = cSysMatrix :: CreateMatrix(SysMatType, NumEq, Profile);
  K->SetParam(MaxIter,Tol);
  StiffnessMatrix(K);
 
  if (Feedback) cout << "\n\tSolving system equations ................" << endl;

  cVector u(NumEq);
  K->Solve(f, u);
  
  // Enforce multi-point constraints: [P]{u} = {b}.

  int nc = cNode :: GetNumMPC( );
  if (nc > 0)
  {
    // Get constraint equations.	  

    cMatrix P(nc, NumEq);
    cVector b(nc); 
    cNode :: GetMPCData(P, b);

    // Solve the equilibrium equations considering the MPCs.

    cVector lm(nc);             // Lagrange multipliers
    cVector ul(NumEq);          // Nodal displacements due to contraints
    cVector fl(NumEq);          // Nodal forces due to contraints
    cMatrix A(nc, nc);          // Auxiliary matrix
    ul.Zero( );
    fl.Zero( );
    K->SolveMPC(P, b, u, A, lm, ul, fl);
    u -= ul;
  }

  if (Feedback) cout << "\n\tPrinting computed results ..............." << endl;

  AssignDispl(u);
  PrintResult( );
}

// -------------------------------------------------------------------------
// Protected methods:
//

// ============================= SolverMicro ===============================

void cLinStat :: SolverMicro(void)
{
  if (Feedback)
  {
    cout << "\n-----------------------------------------------------\n";
    cout << "Micromechanical analysis - Lagrange Multiplier Method\n";
    cout << "-----------------------------------------------------\n";
    cout << "MicroModel = " << MicroModel << "D\n";
    cout << "Regularization Factor = " << scientific << MicroRegFac << endl;
    cout << "{StrainM} = " << scientific;
    MacroStrain.Print( );
  }

  // Get the PBC constraint equations [P]{u} = {b}.	  

  int nc = MicroNumPBC( );
  cVector b(nc); 
#if 0
  cMatrix P(nc, NumEq);
  MicroGetPBC(P, b);
#else
  cSprsMatrix P(nc);
  MicroGetPBC(P, b);
#endif
  if (Feedback)
    cout << "Number of PBC constraints = " << nc << endl;

  if (Feedback) cout << "\n\tComputing the external load vector ......" << endl;

  cVector f(NumEq);
  ExtLoadVector(f);

  if (Feedback) cout << "\n\tAssembling the stiffness matrix ........." << endl;

  cSysMatrix *K = cSysMatrix :: CreateMatrix(SysMatType, NumEq, Profile);
  K->SetParam(MaxIter, Tol);

  // Compute the stiffness matrix and add the regularization factor.

  StiffnessMatrix(K);
  int idx = K->GetSmlPvt( );
  double kmin = K->Get(idx, idx);
  double kreg = kmin*MicroRegFac;
  //cout << "kmin = " << scientific << kmin << "  ";
  //cout << "kreg = " << kreg << endl;
  for (int i = 0; i < NumEq; i++) K->Add(i, i, kreg);
 
  if (Feedback) cout << "\n\tSolving system equations ................" << endl;

  cVector u(NumEq);
  K->Solve(f, u);
  
  // Solve the equilibrium equations considering the PBCs.

  if (Feedback) cout << "\n\tComputing Lagrange multipliers .........." << endl;
  cVector lm(nc);             // Lagrange multipliers
  cVector ul(NumEq);          // Nodal displacements due to contraints
  cVector fl(NumEq);          // Nodal forces due to contraints
  cMatrix A(nc, nc);          // Auxiliary matrix
  ul.Zero( );
  fl.Zero( );
  int code = K->SolveMPC(P, b, u, A, lm, ul, fl);
  if (!code)
  {
    cout << "Error in the micro model solution using the Lagrage Multiplier Method!\n";
    return;
  } 

  // Compute the nodal displacements considering the PBCs.

  u -= ul;

  // Compute the macro stresses using the Lagrange multipliers.

  if (Feedback) cout << "\n\tPerforming homogenization ..............." << endl;
  CalcMacroStrLag(lm, A);

  if (Feedback) cout << "\n\tPrinting computed results ..............." << endl;

  AssignDispl(u);
  PrintResult( );
}

// ======================================================= End of file =====
