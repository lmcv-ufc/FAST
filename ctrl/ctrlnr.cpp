// -------------------------------------------------------------------------
// ctrlnr.cpp - implementation of the control class for nonlinear static
//              problems (Load Control with Newton-Rapson iterations).
// -------------------------------------------------------------------------
// Created:      26-Nov-2000     Evandro Parente Junior
//
// Modified:     03-Oct-2001     Evandro Parente Junior
//               Use of vector/matrix classes.
//
// Modified:     08-Apr-2015     Elias Saraiva Barroso
//               Input/output operations using c++ streams.
//
// Modified:     28-Jul-2025     Evandro Parente Junior
//               Implementation of MPCs.
// -------------------------------------------------------------------------

#include <iostream>
#include <iomanip>

using namespace std;

#include "ctrlnr.h"
#include "vec.h"
#include "sysmat.h"
#include "gblvar.h"
#include "element.h"

// -------------------------------------------------------------------------
// Public methods:
//

// =============================== cCtrlNR =================================

cCtrlNR :: cCtrlNR(void)
{
  Type = LOAD_CONTROL;
}

// =============================== ~cCtrlNR ================================

cCtrlNR :: ~cCtrlNR(void)
{
}

// ================================ Solver =================================

void cCtrlNR :: Solver(void)
{
  // Temporary - Evandro (10-Oct-2025).

  if (MicroModel)
  {
    SolverMicro( );
    return;
  }

  // Create local arrays.

  cVector q(NumEq);            // Reference load vector
  cVector f(NumEq);            // External load vector
  cVector g(NumEq);            // Internal force vector
  cVector r(NumEq);            // Residual vector
  cVector u(NumEq);            // Displacement vector
  cVector du(NumEq);           // Displacement increment

  // Create the tangent stiffness matrix.

  cSysMatrix *K = cSysMatrix :: CreateMatrix(SysMatType, NumEq, Profile);

  // Compute the reference load vector.

  if (PropLoad)
    ExtLoadVector(q);

  // Get MPC data.

  cVector b, bq, h, lm, dlm; 
  cVector dul;         // Nodal displacements due to contraints
  cVector fl;          // Nodal forces due to contraints
  cMatrix P, A;
  int nc = cNode :: GetNumMPC( );
//  cout << "nc = " << nc << endl;
//  cout << scientific << setprecision(2);
  if (nc > 0)
  {
    b.Resize(nc);
    bq.Resize(nc);
    h.Resize(nc);
    lm.Resize(nc);
    dlm.Resize(nc);
    dul.Resize(NumEq);
    fl.Resize(NumEq);
    A.Resize(nc, nc);
    P.Resize(nc, NumEq);
    cNode :: GetMPCData(P, bq);
//    cout << "[P]\n";
//    P.Print( );
//    cout << "{bq} = ";
//    bq.Print( );
    lm.Zero( );
    fl.Zero( );
    dlm.Zero( );
    dul.Zero( );
  }

  // Initialize some variables.

  TotFactor = 0.0;
  u.Zero( );
  g.Zero( );
  IntForceVector(g);   // To include the effects of initial stresses

  // Incremental loop.

  int totiter = 0;
  for (CurrStep = 1; CurrStep <= MaxStep; CurrStep++)
  {
    if (Feedback)
    {
      cout <<"\n\tIncremental step " << right << setw(2) << CurrStep;
      cout << " .....................\n" << endl;
    }

    // Evaluate the external load vector.

    TotFactor += StepFactor;
    if (PropLoad)
      f = TotFactor*q;
    else
      ExtLoadVector(TotFactor, f);
    double flen = f.Length( );

    if (nc > 0)
      b = TotFactor*bq;

    // Evaluate the residual vector.

    r = f - g;
    cout << scientific << setprecision(2);
//    cout << "{r} = ";
//    r.Print( );
    if (nc > 0)
    { 
      fl = t(P)*lm;
//      cout << "{fl} = ";
//      fl.Print( );
      r -= fl;
    } 

    // Iterative loop.

    int conv = 0;
    for (CurrIter = 1; CurrIter <= MaxIter; CurrIter++)
    {
      // Compute the current the stiffness matrix.

      StiffnessMatrix(K);
      totiter++;

      // Compute the displacement increment (Solve [K]{du} = {r}).

//      cout << "{r} = ";
//      r.Print( );
      K->Solve(r, du);
//      cout << "{du} = ";
//      du.Print( );

      // Enforce the MPCs.

      if (nc > 0)
      {
        cVector a(nc);
        a = P*u;
        h = b - a;
//        cout << "{b} = ";
//        b.Print( );
//        cout << "{h} = ";
//        h.Print( );
        K->SolveMPC(P, h, du, A, dlm, dul, fl);
//        cout << "{dul} = ";
//        dul.Print( );
//        cout << "{dlm} = ";
//        dlm.Print( );
        du -= dul;
        lm += dlm;
//        cout << "{lm} = ";
//        lm.Print( );
      }

      // Update internal variabels of each element.

      if (cElement::GetLargeRot3D( ))
      {
        UpdateItera(du);
      }

      // Update the nodal displacements and compute the internal force.

      u += du;
//      cout << "{u} = ";
//      u.Print( );
      AssignDispl(u);
      if (!IntForceVector(g)) break;

      // Compute the new residual.

      r = f - g;
      if (nc > 0)
      { 
        fl = t(P)*lm;
//        cout << "{fl} = ";
//        fl.Print( );
        r -= fl;
      } 

      // Check convergence.

//      cout << "flen = " << flen << "  ";
//      cout << "glen = " << g.Length( ) << "  ";
//      cout << "rlen = " << r.Length( ) << "\n";
      if (flen > Tol)
        conv = Convergence(f, r);
      else
        conv = Convergence(g, r);
      if (conv) break;
    }

    // Print the computed results.

    if (conv)
    {
      // Compute the macro stresses using the Lagrange multipliers.

      if (MicroModel > 0 && nc > 0)
        CalcMacroStrLag(lm, A);

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

  // Print the iteration count.
  
  if (Feedback)
  {
    cout << "\n";
    cout << "TotIter = " << totiter << endl;
    cout << "TotStep = " << CurrStep << endl;
    cout << "AvgIter = " << fixed << setprecision(2) << double(totiter)/double(CurrStep) << endl;
  }
}

// -------------------------------------------------------------------------
// Protected methods:
//

// ============================= SolverMicro ===============================

void cCtrlNR :: SolverMicro(void)
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

  // Create local arrays.

  cVector q(NumEq);            // Reference load vector
  cVector f(NumEq);            // External load vector
  cVector g(NumEq);            // Internal force vector
  cVector r(NumEq);            // Residual vector
  cVector u(NumEq);            // Displacement vector
  cVector du(NumEq);           // Displacement increment

  // Create the tangent stiffness matrix.

  cSysMatrix *K = cSysMatrix :: CreateMatrix(SysMatType, NumEq, Profile);

  // Compute the reference load vector.

  if (PropLoad)
    ExtLoadVector(q);

  // Get the PBC constraint equations [P]{u} = {bq}.	 

  int nc = MicroNumPBC( );
//  cMatrix P(nc, NumEq);       // Constraint matrix
  cSprsMatrix P(nc);          // Constraint matrix
  cVector bq(nc);             // Reference displacement difference
  cVector b(nc);              // {b} = TotFactor*{bq}
  MicroGetPBC(P, bq);
  if (Feedback)
    cout << "Number of PBC constraints = " << nc << endl;

  // Create PBC solver data.

  cVector lm(nc);             // Lagrange multipliers
  cVector dlm(nc);            // Iterative correction of {lm}
  cVector ul(NumEq);          // Nodal displacements due to contraints
  cVector dul(NumEq);         // Iterative correction of {ul}
  cVector fl(NumEq);          // Nodal forces due to contraints
  cVector h(nc);              // Constraint residual {h} = {b} - [P]{u}
  cVector a(nc);              // Auxiliary vector {a} = [P]{u}
  cMatrix A(nc, nc);          // Auxiliary matrix

  // Initialization.

  TotFactor = 0.0;
  u.Zero( );
  g.Zero( );
  lm.Zero( );
  dlm.Zero( );
  dul.Zero( );
  IntForceVector(g);   // To include the effects of initial stresses

  // Incremental loop.

  int totiter = 0;
  for (CurrStep = 1; CurrStep <= MaxStep; CurrStep++)
  {
    if (Feedback)
    {
      cout <<"\n\tIncremental step " << right << setw(2) << CurrStep;
      cout << " .....................\n" << endl;
    }

    // Evaluate the external load vector.

    TotFactor += StepFactor;
    if (PropLoad)
      f = TotFactor*q;
    else
      ExtLoadVector(TotFactor, f);
    double flen = f.Length( );

    // Evaluate the current displacement difference (PBC).

    b = TotFactor*bq;

    // Evaluate the residual vector {r} = {f} - {g} - {fl}.

    r  = f - g;
    P.MultTVect(lm, fl); //fl = t(P)*lm;
    r -= fl;

    // Iterative loop.

    int conv = 0;
    for (CurrIter = 1; CurrIter <= MaxIter; CurrIter++)
    {
      // Compute the stiffness matrix and add the regularization factor.

      StiffnessMatrix(K);
      int idx = K->GetSmlPvt( );
      double kmin = K->Get(idx, idx);
      double kreg = kmin*MicroRegFac;
      //cout << "kmin = " << scientific << kmin << "  ";
      //cout << "kreg = " << kreg << endl;
      for (int i = 0; i < NumEq; i++) K->Add(i, i, kreg);
      totiter++;

      // Compute the displacement increment (Solve [K]{du} = {r}).

      K->Solve(r, du);

      // Enforce the PBCs.

      P.MultVect(u, a); //a = P*u;
      h = b - a;
      int code = K->SolveMPC(P, h, du, A, dlm, dul, fl);
      if (!code)
      {
        cout << "Error in the micro model solution using the Lagrage Multiplier Method!\n";
        return;
      } 
      du -= dul;
      lm += dlm;

      // Update internal variabels of each element.

      if (cElement::GetLargeRot3D( ))
      {
        UpdateItera(du);
      }

      // Update the nodal displacements and compute the internal force.

      u += du;
      AssignDispl(u);
      if (!IntForceVector(g)) break;

      // Compute the new residual {r} = {f} - {g} - {fl}.

      r  = f - g;
      P.MultTVect(lm, fl); //fl = t(P)*lm;
      r -= fl;

      // Check convergence.

      if (flen > Tol)
        conv = Convergence(f, r);
      else
        conv = Convergence(g, r);
      if (conv) break;
    }

    // Print the computed results.

    if (conv)
    {
      // Compute the macro stresses using the Lagrange multipliers.

      CalcMacroStrLag(lm, A);

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

  // Print the iteration count.
  
  if (Feedback)
  {
    cout << "\n";
    cout << "TotIter = " << totiter << endl;
    cout << "TotStep = " << CurrStep << endl;
    cout << "AvgIter = " << fixed << setprecision(2) << double(totiter)/double(CurrStep) << endl;
  }
}

// ======================================================= End of file =====
