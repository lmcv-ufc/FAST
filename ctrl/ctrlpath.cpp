// -------------------------------------------------------------------------
// ctrlpath.cpp - implementation of the control class for nonlinear static
//                problems (Path-Following Methods)
// -------------------------------------------------------------------------
// Created:      22-Out-2011     Iuri Barcelos Rocha
//
// Modified:
// -------------------------------------------------------------------------

#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <math.h>

using namespace std;

#include "ctrlpath.h"
#include "ctrlnlstab.h"
#include "node.h"
#include "vec.h"
#include "sysmat.h"
#include "gblvar.h"
#include "element.h"

// -------------------------------------------------------------------------
// Static variables
//
ePathType cCtrlPath :: PathType;
int       cCtrlPath :: CritPoint = 0;
int       cCtrlPath :: BranchSwt = 0;
double    cCtrlPath :: SwtFactor = 1.0e-3;

// -------------------------------------------------------------------------
// Public methods:
//

// ============================ ReadStabData ===============================

void cCtrlPath :: ReadStabData(void)
{
  if (!(in >> CritPoint) || !(in >> BranchSwt) || !(in >> SwtFactor))
  {
    cout << "Error in the input of the stability data!\n";
    exit(0);
  }
}

// ============================== cCtrlPath ================================

cCtrlPath :: cCtrlPath (void)
{
  Type = PATH_FOLLOWING;
}

// ============================= ~cCtrlPath ================================

cCtrlPath :: ~cCtrlPath (void)
{

}
// =============================== Solver ==================================

void cCtrlPath :: Solver(void)
{
  // Create local vectors.

  cVector q(NumEq);            // Reference load vector
  cVector f(NumEq);            // External load vector
  cVector g(NumEq);            // Internal force vector
  cVector r(NumEq);            // Residual vector
  cVector u(NumEq);            // Displacement vector
  cVector du1(NumEq);          // [K]{du1} = {q}
  cVector du2(NumEq);          // [K]{du2} = {r}
  cVector du(NumEq);           // Displacement increment (iteration)
  cVector Du(NumEq);           // Displacement increment (step)
  cVector DuLast(NumEq);       // Last step displacement increment
  cVector fs(NumEq);         // Strain loads (d{g}/dlf)
  cVector qt(NumEq);         // {qt} = {q} + {fs}

  // Create the tangent stiffness matrix.

  cSysMatrix *K = cSysMatrix :: CreateMatrix(SysMatType, NumEq, Profile);

  // Compute the reference load vector.

  ExtLoadVector(q);
  double lenq = q.Length( );

  // Initialize some variables.

  u.Zero( );
  r.Zero( );
  TotFactor = 0.0;
  int nnp0 = 0;
  int nnp1 = 0;
  CSP0 = 1.0;
  CSP1 = 1.0;
  double sp0, sp;

  // Incremental loop.

  int totiter = 0;
  for (CurrStep = 1; CurrStep <= MaxStep; CurrStep++)
  {
    if (Feedback)
    {
      cout << "\n\tIncremental step " << right << setw(2) << CurrStep;
      cout << " .....................\n\n";
    }

    // Iterative loop.

    int conv = 0;
    double dlf, Dlf;
    for (CurrIter = 1; CurrIter <= MaxIter; CurrIter++)
    {
      // Compute the current stiffness matrix.

      StiffnessMatrix(K);
      totiter++;

      // Include strain loads (ex: thermal effects).

      qt = q;
      if (ExistStrainLoads( ))
      {
        if (StrainLoadVector(g, fs)) qt += fs;
      }

      // Compute displacements increments.

      K->Solve(qt, du1);  // [K]{du1} = {q}
      K->Solve(r, du2);   // [K]{du2} = {r}

      // Compute the incremental load factor. Do not change the order of these two blocks!

      dlf = IncLoadFac(du1, du2, Du);
      if (CurrIter == 1)
      {
        DuLast = Du;
        Du.Zero();
        Dlf = 0.0;
      }

      // Compute the incremental displacements.

      du = dlf*du1 - du2;

      // Update internal variables of each element witd large rotation.

      if (cElement::GetLargeRot3D( ))
      {
        UpdateItera(du);
      }

      // Update the load factor.

      Dlf += dlf;
      TotFactor += dlf;

      // Update the nodal displacements.

      u += du;
      Du += du;

      // Check for restarts.

      if (CurrIter == 0)
      {
        u -= Du;
        TotFactor -= Dlf;
        Du = DuLast;
      }

      // Update the internal and external force vectors and assign displacements.

      f = TotFactor*q;
      AssignDispl(u);
      if (!IntForceVector(g)) break;

      // Compute the new residual and check convergence.

      r = g - f;
      if (CurrIter != 0)
      {
        conv = Convergence(lenq, r.Length( ));
        if (conv) break;
      }
    }

    // Print the computed results.

    if (conv)
    {
      // Compute the current stiffness parameter.

      double d2 = du1*du1;
      if (sqrt(d2) > 1.0e-12)
        sp = (qt*du1)/d2;
      else
        sp = 1.0;
      if (CurrStep == 1)
      {
        sp0 = sp;
      }
      else
      {
        CSP0 = CSP1;
        CSP1 = sp/sp0;
      }

      // Stability analysis.

      if (CritPoint)
      {
        // Compute the number of negative pivots of [K] at the current point.

        nnp0 = nnp1;         // Store old value
        cVector v(NumEq);    // Create the eigenvector
        StiffnessMatrix(K);  // Assembly the current stiffness matrix
        K->Solve(du, v);     // Only to decompose [K]
        nnp1 = K->GetNgtPvt( );
        cout << "CSP0 = " << CSP0 << "  CSP1 = " << CSP1 << "  nnp0 = " << nnp0 << "  nnp1 = " << nnp1 << "\n";

        // Compute the critical point and perform branch-switching.

        if (nnp1 != nnp0)
        {
          int code = StabSolver(K, u, Du, du, du1, du2, q, r, f, g, v);
          if (code == 2) // Bifurcation point
            nnp1 = K->GetNgtPvt( );
        }
      }

      // Print results

      //if (Feedback) cout << "\n\tPrinting computed results ..............." << endl;
      PrintResult( );
    }
    else
    {
      cout << "Convergence not achieved!!!\n";
      break;
    }

    // Update step variables.

    PrevIter = CurrIter;
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

// =============================== StabSolver ==============================
//
// This function returs 1 the primary will be traced and 2 if a branch-switch
// to the secondary path was performed.
//
int cCtrlPath :: StabSolver(cSysMatrix *K, cVector &u, cVector &Du,
                            cVector &du, cVector &du1, cVector &du2,
                            cVector &q, cVector &r, cVector &f, cVector &g,
                            cVector &v)
{
  // Store the current equilibrium state.

  double  tfac = TotFactor;
  cVector unxt = u;

  // Evaluate and print the critical point.

  int code;
  cCtrlNlStab *stab = new cCtrlNlStab( );
  stab->Solver(K, u, du1, du2, q, r, f, g, v, &code);

  // Perform branch-switching for bifurcation points.

  if (BranchSwt == 1 && code == 2)
  {
    if (Feedback) cout << "\n\tBranch-switching ........................\n\n";
    BranchSwt = 0; // Perform only one branch-switch
		  
    // Set default number of target iterations.

    if (GetTargetIter( ) == 1)
    {
      SetTargetIter(3);
    }

    // Create auxiliary vectors.

    cVector fs(NumEq);         // Strain loads (d{g}/dlf)
    cVector qt(NumEq);         // {qt} = {q} + {fs}*/

    // Get the controlling dof.

    int idx = 0;
    for (int k = 1; k < NumEq; k++)
    {
      if (fabs(v[k]) > v[idx])
      {
        idx = k;
      }
    }

    // Predictor => dlf = 0.0 and {du} = fac*{v}.

    double fac = SwtFactor/v[idx];
    du1= fac*v;
    Du = du1;
    u += du1;
    AssignDispl(u);
    IntForceVector(g);
    r = g - f;

    // Perform equilibrium iterations.

    int conv = 0;
    for (CurrIter = 1; CurrIter <= MaxIter; CurrIter++)
    {
      StiffnessMatrix(K);

      // Include strain loads (ex: thermal effects)
      qt = q;
      if (ExistStrainLoads( ))
      {
        if (StrainLoadVector(g, fs)) qt += fs;
      }

      K->Solve(qt, du1);  // [K]{du1} = {q}
      K->Solve(r, du2);  // [K]{du2} = {r}

      double dlf = (Du*du2)/(Du*du1);
      TotFactor += dlf;
      f  = TotFactor*q;
      du = dlf*du1 - du2;
      u  += du;
      Du += du;

      AssignDispl(u);
      IntForceVector(g);
      r = g - f;
      conv = Convergence(q.Length( ), r.Length( ));
      if (conv) return(2);
    }
  }

  // Return to the last computed point in the primary path in case of
  // limit points or lack of convergence in the branch-switching.

  TotFactor = tfac;
  u = unxt;
  AssignDispl(u);

  return(1);
}

// ======================================================= End of file =====
