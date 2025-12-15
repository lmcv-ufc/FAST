// -------------------------------------------------------------------------
// ctrlnlstab.cpp - implementation of the control class for nonlinear
// stability analysis.
// -------------------------------------------------------------------------
// Created:      01-Oct-2017     Evandro Parente Junior
//               Based on stabses.cpp from FEMOOP.
//
// Modified:
// -------------------------------------------------------------------------

#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <math.h>

using namespace std;

#include "ctrlnlstab.h"
#include "element.h"
#include "node.h"
#include "vec.h"
#include "sysmat.h"
#include "gbldef.h"
#include "gblvar.h"

// -------------------------------------------------------------------------
// Static variables : declaration and initialization
//
double cCtrlNlStab :: RelPert  = 1.0e-05;

// -------------------------------------------------------------------------
// Public methods:

// ============================== ReadData =================================

void cCtrlNlStab :: ReadData(void)
{
  if (!(in >> RelPert))
  {
    cout << "Error in the input of the number of relative perturbation!\n";
    exit(0);
  }
}

// ============================= cCtrlNlStab ===============================

cCtrlNlStab :: cCtrlNlStab(void)
{
  Type = NONLIN_STABILITY;
}

// =============================== Solver ==================================

void cCtrlNlStab :: Solver(void)
{
  // Create the global stiffness matrix.

  cSysMatrix *K = cSysMatrix :: CreateMatrix(SysMatType, NumEq, Profile);

  // Get the reference load vector and initialize the load factor.

  cVector q(NumEq);
  ExtLoadVector(q);
  TotFactor = 0.0;

  // Create local vectors and call the solver.

  int code;
  cVector u  (NumEq);    // Total displacements
  cVector du1(NumEq);    // Displacement increment
  cVector du2(NumEq);    // Displacement increment
  cVector r  (NumEq);    // Residual (unbalanced forces)
  cVector f  (NumEq);    // External load
  cVector g  (NumEq);    // Internal force
  cVector v  (NumEq);    // Eigenvector
  u.Zero( );
  r.Zero( );
  Solver(K, u, du1, du2, q, r, f, g, v, &code);

  // Release memory.

  delete K;
}

// =============================== Solver ==================================

void cCtrlNlStab :: Solver(cSysMatrix *K, cVector &u, cVector &du1,
                           cVector &du2, cVector &q, cVector &r, cVector &f,
                           cVector &g, cVector &v, int *code)
{
  if (Feedback) cout << "\n\tNonlinear stability analysis ...........\n\n";

  // Create local vectors.

  cVector dv1(NumEq);    // Eigenvector increment
  cVector dv2(NumEq);    // Eigenvector increment
  cVector w  (NumEq);    // Auxiliary vector
  cVector h1 (NumEq);    // Directional derivative
  cVector h2 (NumEq);    // Directional derivative
  cVector fs(NumEq);     // Strain loads (d{g}/dlf)
  cVector qt(NumEq);     // {qt} = {q} + {fs}


  // Initialize some variables.

  *code = 0;
  double lenq = q.Length( );
  if (lenq < Tol) lenq = 1.0;
  f.Zero( );
  r.Zero( );

  // Iterative process.

  double err = 1.0;
  for (CurrIter = 1; CurrIter <= MaxIter; CurrIter++)
  {
    // Compute the current stiffness matrix.

    StiffnessMatrix(K);

    // Include strain loads (ex: thermal effects)
    qt = q;
    if (ExistStrainLoads( ))
     {
      if (StrainLoadVector(g, fs)) qt += fs;
     }
      //SUBSTITUIR q POR qt

    // Evaluate [K] {du1} = {qt} and [K] {du2} = {r}.

    K->Solve(qt, du1);
    K->Solve(r, du2);

    // Evaluate the initial eigenvector (two steps of inverse iteration).

    if (CurrIter == 1)
    {
      v = 1.0;
      v /= v.Length( );
      for (int k = 0; k < 3; k++)
      {
        K->Solve(v, v);
        v /= v.Length( );
      }
    }

    // Evaluate the "directional derivative vectors" => {h1} and {h2}.

    DirecDeriv(u, v, du1, du2, qt, r, h1, h2);

    // Evaluate [K] {dv1} = {h1} and [K] {dv2} = {h2}.

    K->Solve(h1, dv1);
    K->Solve(h2, dv2);

    // Evaluate the load increment and update load factor.

    double lenv = v.Length( );
    double dot1 = v*dv1;
    double dot2 = v*dv2;
    double dlf  = (dot2 - lenv)/dot1;

    // Evaluate the incremental vectors.

    TotFactor += dlf;
    u += dlf*du1 - du2;
    v  = dv2 + (-dlf)*dv1;

    // Update the internal and external force vectors and assign displacements.

    f = TotFactor*q;
    AssignDispl(u);
    if (!IntForceVector(g)) break;

    // Evaluate the unbalanced force vector => {r} = {g} - {f}.

    r = g - f;
    err = r.Length( )/lenq;   // ||r||/||q||

    // Evaluate an estimate of the eigenvalue.

    K->Solve(v, w);
    double mu = v*v/(w*v);

    // Check convergence.

    if (Feedback)
    {
      cout << "Iter: "   << right << setw(3) << CurrIter;
      cout << "  EigV: " << scientific << setw(11) << setprecision(4) << mu;
      cout << "  Err: "  << scientific << setprecision(4) << err;
      cout << "  LF: "   << fixed << showpos << TotFactor << "\n" << noshowpos;
    }
    if (err < Tol && fabs(mu) < Tol && CurrIter > 1) break;
  }

  // Print results.

  if (err < Tol)
  {
    // Classify the critical point.

    double cosqv = (q*v)/lenq;
    *code = 1;
    if (fabs(cosqv) <= max(Tol, 1.0e-3)) *code = 2;
    if (Feedback)
    {
      cout << "{f}.{v}/|f||v| = " << scientific << cosqv;
      if (*code == 1)
        cout << " => limit point\n";
      else
        cout << " => bifurcation point\n";
    }

    // Normalize the critical mode.

    int imax = 0;
    for (int i = 1; i < NumEq; i++)
      if (fabs(v[i]) >  fabs(v[imax])) imax = i;
    v /= v[imax];

    if (Feedback) cout << "\n\tPrinting computed results ...............\n";
    PrintResult( );
    PrintMode(u, v);
  }
  else
  {
    if (Feedback) cout << "Convergence not achived!!!\n";
  }
}


// -------------------------------------------------------------------------
// Protected methods:
//

// ============================== PrintMode ================================

void cCtrlNlStab :: PrintMode(cVector &u, cVector &v)
{
  AssignDispl(v);

  out << "%RESULT.CASE.STEP.NODAL.EIGENVECTOR\n";
  int nnode = cNode :: GetNumNode( );
  out << nnode << "  'buckling_mode'\n";
  for (int i = 0; i < nnode; i++)
  {
    // Get the current node.

    cNode *node = cNode :: GetNode(i);
    out << left << setw(4) << fixed << node->GetLabel( ) << " ";

    // Get the nodal displacements.

    out << scientific << right << setprecision(5);
    for (int j = 0; j < 6; j++) out << setw(13) << node->GetDispl(j) << " ";
    out << "\n";
  }
  out << "\n";

  AssignDispl(u);
}

// ============================= DirecDeriv ================================

void cCtrlNlStab :: DirecDeriv(cVector &u, cVector &v, cVector &du1, cVector &du2,
                               cVector &p, cVector &q, cVector &h1, cVector &h2)
{
  // Initialize the directional derivative vectors.

  h1.Zero( );
  h2.Zero( );

  // Choose the pertubation value.

  double ratio = u.Length( )/v.Length( );
  double delta = RelPert*MAX(ratio, 1.0);

  // Evaluate forward displacement vector and perturb the displacement field.

  cVector up(NumEq);
  up = u + delta*v;
  AssignDispl(up);

  // Alloc the element stiffness matrix (max size).

  int ndof = cElement :: GetMaxDof( );
  cMatrix Kelm(ndof, ndof);
  cVector a(ndof);
  cVector b(ndof);

  // Compute and add the contribution of each element to {h1} and {h2}.

  int nelm = cElement :: GetNumElm( );
  for (int i = 0; i < nelm; i++)
  {
    // Get the current element.

    cElement *elm = cElement :: GetElm(i);

    // Resize arrays if necessary.

    ndof = elm->GetNumElmDof( );
    if (ndof != Kelm.NRow( ))
    {
      Kelm.Resize(ndof, ndof);
      a.Resize(ndof);
      b.Resize(ndof);
    }

    // Evaluate the element stiffness matrix.

    elm->StiffMat(Kelm);

    // Add [Kelm]{du1} to {h1}.

    elm->GlobToElm(du1, a);
    b = Kelm*a;
    elm->AddGlobVec(b, h1);

    // Add [Kelm]{du2} to {h2}.

    elm->GlobToElm(du2, a);
    b = Kelm*a;
    elm->AddGlobVec(b, h2);
  }

  // Evaluate the derivative vectors by forward differences.

  for (int i = 0; i < NumEq; i++)
  {
    h1[i] = (h1[i] - p[i])/delta;
    h2[i] = (h2[i] - q[i])/delta;
  }

  // Return to the current displacement field.

  AssignDispl(u);

  // Include the thermal effects
  if (ExistStrainLoads( ))
  {

    double lf = TotFactor;
    TotFactor = lf + delta;
    cVector hlf(NumEq);
    hlf.Zero( );


   for (int i = 0; i < nelm; i++)
  {
    // Get the current element.

    cElement *elm = cElement :: GetElm(i);

    // Resize arrays if necessary.

    ndof = elm->GetNumElmDof( );
    if (ndof != Kelm.NRow( ))
    {
      Kelm.Resize(ndof, ndof);
      a.Resize(ndof);
      b.Resize(ndof);
    }

    // Evaluate the element stiffness matrix.

    elm->StiffMat(Kelm);

    elm->GlobToElm(v, a);
    b = Kelm*a;
    elm->AddGlobVec(b, hlf);
  }

  TotFactor = lf - delta;


  for (int i = 0; i < nelm; i++)
  {
    // Get the current element.
    cElement *elm = cElement :: GetElm(i);

    // Resize arrays if necessary.

    ndof = elm->GetNumElmDof( );
    if (ndof != Kelm.NRow( ))
    {
      Kelm.Resize(ndof, ndof);
      a.Resize(ndof);
      b.Resize(ndof);
    }

    // Evaluate the element stiffness matrix.

    elm->StiffMat(Kelm);

    elm->GlobToElm(v, a);
    b = Kelm*a;
    b *= -1.0;
    elm->AddGlobVec(b, hlf);

  }

  hlf /= (2.0*delta);
  h1 += hlf;
  TotFactor = lf;
  }

}

// ============================================================= End of file
