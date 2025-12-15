// -------------------------------------------------------------------------
// ctrllinstab.cpp - implementation of the control class for linearized
//                   buckling analysis.
// -------------------------------------------------------------------------
//
// Created:      04-Dec-2011     Iuri Barcelos Rocha
//
// Modified:
// -------------------------------------------------------------------------

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <algorithm>

using namespace std;

#include "ctrl.h"
#include "ctrllinstab.h"
#include "node.h"
#include "vec.h"
#include "sysmat.h"
#include "eig.h"
#include "gbldef.h"
#include "gblvar.h"
#include "shpiga.h"


// -------------------------------------------------------------------------
// Public methods:
//

// ============================= cCtrlLinStab ==============================

cCtrlLinStab :: cCtrlLinStab(void)
{
  Type = LINEAR_STABILITY;
}

// ============================ ~cCtrlLinStab ==============================

cCtrlLinStab :: ~cCtrlLinStab(void)
{
}

// ================================ Solver =================================

void cCtrlLinStab :: Solver(void)
{
  if (Feedback) cout << "\n\tLinear analysis ........................\n";

  // Get eigen solver type.

  eEigAlgType type = static_cast<eEigAlgType> (EigAlgType);

  // Using full profile in the case of Generalized Jacobi method.

  int* profile = new int[NumEq];
  if (type == GENERALIZED_JACOBI)
    for (int i = 0; i < NumEq; ++i) profile[i] = 0.0;
  else
    for (int i = 0; i < NumEq; ++i) profile[i] = Profile[i];

  // Get the reference load vector.

  TotFactor = 1.0;
  cVector f(NumEq);
  ExtLoadVector(f);

  // Get the stiffness matrix.

  cSysMatrix *K = cSysMatrix :: CreateMatrix(SysMatType, NumEq, profile);
  StiffnessMatrix(K);

  // Evaluate [K]{u} = {f}

  cVector u(NumEq);
  K->Solve(f, u);

  // Assign displacements and print results

  AssignDispl(u);
  PrintResult( );

  // Get the geometric stiffness matrix.

  if (Feedback) cout << "\n\tGeneralized eigenproblem ................" << endl;

  cSysMatrix *G = cSysMatrix :: CreateMatrix(SysMatType, NumEq, profile);
  GeomStiffMatrix(G);
  G->AddMat(-2.0, G);  // G = -1*G

  // Create eigenproblem solver object.

  cGenEigenProb *EigSolver = cGenEigenProb :: CreateEigenProbAlg(type);

  // Get the stiffness matrix.

  if (type == GENERALIZED_JACOBI)
    StiffnessMatrix(K);

  // Solve the eigenproblem: ([K] + lbd [G]){v} = {0}.

  EigSolver->setNumPair(NumModes);
  int conv = EigSolver->Solver(K, G, MaxIter, Tol);

  // Print the computed results.

  if (conv)
  {
    if (Feedback) cout << "\n\tPrinting computed results ..............." << endl;
    
    sEigPair *pair;
    int nmodes = min(EigSolver->getNumPair( ), NumModes);
    out << "\n%RESULT.CASE.STEP.BUCKLING.MODES\n";
    out << nmodes << "\n"; 

    for (int i = 0; i < nmodes; i++)
    {
      pair = EigSolver->getPair(i);  
      PrintMode(i, pair->lbd, *(pair->EigVec));
    }

    if (PrintEqvMsh)
    {
      outem << "\n%RESULT.CASE.STEP.BUCKLING.MODES\n";
      outem << nmodes << "\n";

      for (int i = 0; i < nmodes; i++)
      {
        pair = EigSolver->getPair(i);
        cShapeIGA :: PrintMode(i, pair->lbd, *(pair->EigVec),outem,"BUCKLING.FACTOR","Buckling mode");
      }
    }
  }
  else
  {
    cout << "\nConvergence not achieved in stability analysis !!!" << endl;
  }

  // Release memory.

  delete K;
  delete G;
  delete []profile;
  delete EigSolver;  
}

// -------------------------------------------------------------------------
// Protected methods:
//

// =============================== PrintMode ===============================

void cCtrlLinStab :: PrintMode(int mode, double lbd, cVector &v)
{
  if (Feedback && mode < 10)
  {	  
    cout << "Mode " << setw(2) << (mode + 1) << " lbd = "; 
    cout << scientific << setprecision(OutPrec) << lbd << "\n"; 
  }

  out << "\n%RESULT.CASE.STEP.BUCKLING.FACTOR\n";
  out << mode + 1 << "  ";
  out << showpos << right << scientific << setw(14) << setprecision(OutPrec) << lbd << "\n"; 

  int nn = cNode :: GetNumNode( );
  out << "\n%RESULT.CASE.STEP.NODAL.EIGENVECTOR";
  out << noshowpos << left << "\n" << nn << "  'Buckling mode'\n";
  for (int i = 0; i < nn; i++)
  {
    cNode *node = cNode :: GetNode(i);
    out << setw(4) << node->GetLabel( ) << " " << showpos << setprecision(OutPrec) << scientific;
    for (int j = 0; j < 6; j++)
    {
      int dof  = node->GetDof(j);
      double val = 0.0;
      if (dof > 0) val = v[dof-1];
      out << " " << val;
    }
    out << "\n" << noshowpos;
  }
  out << "\n";
}

// ======================================================= End of file =====
