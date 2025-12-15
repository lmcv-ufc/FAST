// -------------------------------------------------------------------------
// ctrlmodal.cpp - implementaton of the control class for modal analysis.
// -------------------------------------------------------------------------
//
// Created:      29-Jun-2011     Rafael Fernandes da Silva
//
// Modified:     07-Oct-2017     Elias Saraiva Barroso
//               Implementation of cGenEigenProb class.
// -------------------------------------------------------------------------

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <algorithm>

using namespace std;

#include "ctrl.h"
#include "ctrlmodal.h"
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

// ============================== cCtrlModal ===============================

cCtrlModal :: cCtrlModal(void)
{
  Type = MODAL;
}

// ============================== ~cCtrlModal ==============================

cCtrlModal :: ~cCtrlModal(void)
{
}

// ================================ Solver =================================

void cCtrlModal :: Solver(void)
{
  if (Feedback) cout << "\n\tModal analysis ........................\n";

  // Get eigen solver type.

  eEigAlgType type = static_cast<eEigAlgType> (EigAlgType);

  // Using full profile in the case of Generalized Jacobi method.
 
  int* profile = new int[NumEq];
  if (type == GENERALIZED_JACOBI)
    for(int i = 0; i < NumEq; ++i) profile[i] = 0.0;
  else
    for(int i = 0; i < NumEq; ++i) profile[i] = Profile[i];

  // Compute the mass matrix.

  cSysMatrix *M = new cSymSkylMatrix(NumEq, profile);
  MassMatrix(M);

  // Compute the stiffness matrix.

  cSysMatrix *K = new cSymSkylMatrix(NumEq, profile);
  StiffnessMatrix(K);

  // Create eigenproblem solver object.

  cGenEigenProb *EigSolver = cGenEigenProb :: CreateEigenProbAlg(type);

  // Solve the eigenproblem: [K]{v} = lbd*[M]{v}.

  EigSolver->setNumPair(NumModes);
  int flag = EigSolver->Solver(K, M, MaxIter, Tol);

  // Print the computed results.

  if (flag)
  {
    if (Feedback) cout << "\n\tPrinting computed results ..............." << endl;
    
    sEigPair *pair;
    int nmodes = min(EigSolver->getNumPair( ), NumModes);
    out << "\n%RESULT.CASE.STEP.VIBRATION.MODES\n";
    out << nmodes << "\n";

    for (int i = 0; i < nmodes; i++)
    {
      pair = EigSolver->getPair(i);  
      double w = sqrt(pair->lbd);         // Evaluate the natural frequency
      PrintMode(i, w, *(pair->EigVec));
    }

    if (PrintEqvMsh)
    {
      outem << "\n%RESULT.CASE.STEP.VIBRATION.MODES\n";
      outem << nmodes << "\n";

      for (int i = 0; i < nmodes; i++)
      {
        pair = EigSolver->getPair(i);
        cShapeIGA :: PrintMode(i, pair->lbd, *(pair->EigVec),outem,"NATURAL.FREQUENCY","Vibration mode");
      }
    }

  }
  else
  {
    cout << "\nConvergence not achieved in stability analysis !!!" << endl;
  }

  // Release memory.

  delete K;
  delete M;
  delete []profile;
  delete EigSolver;  
}

// =============================== PrintMode ===============================

void cCtrlModal :: PrintMode(int mode, double freq, cVector &v)
{
  if (Feedback && mode < 10)
    cout << "Mode " << setw(2) << (mode + 1) << " w = " << freq << "\n";

  out << "\n%RESULT.CASE.STEP.NATURAL.FREQUENCY\n";
  out << mode + 1 << "  ";
  out << showpos << right <<scientific << setw(14) << setprecision(OutPrec) << freq << "\n";

  int nn = cNode :: GetNumNode( );
  out << "\n%RESULT.CASE.STEP.NODAL.EIGENVECTOR";
  out << noshowpos << left << "\n" << nn << "  'Vibration mode'\n";

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
