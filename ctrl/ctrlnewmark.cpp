// -------------------------------------------------------------------------
// ctrlnewmark.cpp - implementation of the control classes for dynamic
// analysis using the Newmark's algorithm.
// -------------------------------------------------------------------------
// Created:      14-Mar-2020     Evandro Parente Junior
//               Based on old ctrlnwmk, ctrlincnmk.cpp, ctrlnlnmk.cpp files.
//
// Modified:
// -------------------------------------------------------------------------

#include <iostream>
#include <iomanip>
#include <stdlib.h>

using namespace std;

#include "ctrlnewmark.h"
#include "node.h"
#include "vec.h"
#include "sysmat.h"
#include "gbldef.h"
#include "gblvar.h"

// -------------------------------------------------------------------------
// Static variables:
//
double cCtrlNewmark :: Gamma = 0.50; // Constant average acceleration method
double cCtrlNewmark :: Beta  = 0.25; // Constant average acceleration method
bool   cCtrlNewmark :: Damping = false;
double cCtrlNewmark :: Omega1  = 0.0;      
double cCtrlNewmark :: Omega2  = 1.0;
double cCtrlNewmark :: Xi1     = 0.0;       
double cCtrlNewmark :: Xi2     = 0.0;       

// -------------------------------------------------------------------------
// Public methods:
//

// =========================== ReadNewmarkParam ============================

void cCtrlNewmark :: ReadNewmarkParam(void)
{
  if (!(in >> Gamma) || !(in >> Beta))
  {
    cout << "Error in the input of Newmark parameters!\n";
    exit(0);
  }
}

// ============================ ReadRayleighDamp ===========================

void cCtrlNewmark :: ReadRayleighDamp(void)
{
  if (!(in >> Omega1) || !(in >> Omega2) || !(in >> Xi1) || !(in >> Xi2))
  {
    cout << "Error in the input of Rayleigh Damping parameters!\n";
    exit(0);
  }
  Damping = true;
  cout << "Omega1 = " << Omega1  << "  Omega2 = " << Omega2 << endl; 
  cout << "Xi1 = " << Xi1  << "  Xi2 = " << Xi2 << endl; 
}

// ============================= cCtrlNewmark ==============================

cCtrlNewmark :: cCtrlNewmark(void)
{
  Type = NEWMARK;
}

// ============================= ~cCtrlNewmark =============================

cCtrlNewmark :: ~cCtrlNewmark(void)
{
}

// ================================ Solver =================================

void cCtrlNewmark :: Solver(void)
{
  if (Feedback)
  {	  
    cout << "\n\tLinear Dynamic Analysis - Newmark Method\n";
    cout << "\tGamma = " << Gamma << endl;
    cout << "\tBeta  = " << Beta  << endl;
  }

  // Create local vectors.

  cVector f(NumEq);            // External load vector
  cVector u(NumEq);            // Displacement vector
  cVector du(NumEq);           // Incremental displacements
  cVector v(NumEq);            // Velocity vector
  cVector a(NumEq);            // Acceleration vector
  cVector w(NumEq);            // Auxiliary vector
  cVector z(NumEq);            // Auxiliary vector

  // Initialize algorithm variables.

  f.Zero( );
  u.Zero( );
  du.Zero( );
  v.Zero( );
  a.Zero( );
  w.Zero( );
  TotTime = 0.0;

  // Calculate integration constants.

  double a0 = 1.0/(Beta*TimeStep*TimeStep);
  double a1 = Gamma/(Beta*TimeStep);
  double a2 = 1.0/(Beta*TimeStep);
  double a3 = 0.5/Beta - 1.0;
  double a4 = Gamma/Beta - 1.0;
  double a5 = (0.5*Gamma/Beta - 1.0)*TimeStep;
  double a6 = (1.0 - Gamma)*TimeStep;
  double a7 = Gamma*TimeStep;
#if 1
  printf("a0 = %f\n", a0);
  printf("a1 = %f\n", a1);
  printf("a2 = %f\n", a2);
  printf("a3 = %f\n", a3);
  printf("a4 = %f\n", a4);
  printf("a5 = %f\n", a5);
  printf("a6 = %f\n", a6);
  printf("a7 = %f\n", a7);
#endif

  // Compute the mass matrix.

  cSysMatrix *M = new cSymSkylMatrix(NumEq, Profile);
  MassMatrix(M);
//  cout << "[M] = ";
//  M->Print( );

  // Compute the acceleration vector at t = 0.0 (Solve [M]{a} = {f}).

  ExtLoadVector(TotTime, f);
  cSysMatrix *K = new cSymSkylMatrix(NumEq, Profile);
  K->Zero( );
  K->AddMat(1.0, M);       // [K] = [M]
  K->Solve(f, a);

  // Compute the stiffness and damping matrices.

  StiffnessMatrix(K);
//  cout << "[K] = ";
//  K->Print( );
  cSysMatrix *C = 0;
  if (Damping)
  {	  
    C = new cSymSkylMatrix(NumEq, Profile);
    DampingMatrix(M, K, C);
//  cout << "[C] = ";
//  C->Print( );
  }

  // Compute [K]eff = [K] + a0*[M] + a1*[C].

  K->AddMat(a0, M);
  if (Damping) K->AddMat(a1, C);
//  cout << "[K]eff = ";
//  K->Print( );

  // Incremental loop.

  double umax = 0.0;
  double vmax = 0.0;
  double amax = a.NormInf( );
  for (CurrStep = 1; CurrStep <= MaxStep; CurrStep++)
  {
    // Update the total time and the external force vector.

    TotTime += TimeStep;
    TotFactor = TotTime;   // For printing in output file
    if (Feedback)
    {
      cout << "\n\tStep = " << setw(3)  << CurrStep;
      cout << " Time = " << scientific << setprecision(3); 
      cout << TotTime << " ............." << endl << resetiosflags(ios::scientific);
    }
    ExtLoadVector(TotTime, f);
//    cout << "{f} = " << scientific << setprecision(4) << showpos;
//    f.Print( );

    // Add the inertia terms: {f} += [M]*(a0*{u} + a2*{v} + a3*{a}).

    w  = a0*u;
    w += a2*v;
    w += a3*a;                        // {w} = a0*{u} + a2*{v} + a3*{a}
    M->MultVect(w.Val( ), z.Val( ));  // {z} = [M]{w}
    f += z;

    // Add the damping effect: {f} += [C]*(a1*{u} + a4*{v} + a5*{a}).

    if (Damping)
    {
      w  = a1*u;
      w += a4*v;
      w += a5*a;                        // {w} = a1*{u} + a4*{v} + a5*{a}
      C->MultVect(w.Val( ), z.Val( ));  // {z} = [C]{w}
      f += z;
    }

    // Compute the displacement vector (Solve [K]{u} = {f}).

    w = u;           // Store the old displacements
    K->Solve(f, u);
    du = u - w;      // Compute the displacement increments

    // Compute the new acceleration vector.

    w = a;           // Store the old acceleration
    a *= -a3;
    a +=  a0*du;
    a += -a2*v;

    // Compute the new velocity vector.

    v += a6*w;       // Term depeding on the old acceleration
    v += a7*a;       // Term depeding on the new acceleration

    // Update the nodal displacements for stress computations.

    AssignDispl(u);

    // Print the computed results.

    umax = MAX(umax, u.NormInf( ));
    vmax = MAX(vmax, v.NormInf( ));
    amax = MAX(amax, a.NormInf( ));
    if (CurrStep%PrintStep == 0)
    {
      if (Feedback) cout << "\n\tPrinting computed results ..............." << endl;
      PrintResult( );
      PrintVelAcc(v, a);
    }

    // Update element variables.

    UpdateState( );
  }

  if (Feedback)
  {
    cout << "umax = " << scientific << setprecision(3) << umax << endl;
    cout << "vmax = " << scientific << setprecision(3) << vmax << endl;
    cout << "amax = " << scientific << setprecision(3) << amax << endl;
  }

  // Adjust the current step.

  CurrStep--;
}

// -------------------------------------------------------------------------
// Protected methods:
//

// ============================ DampingMatrix  =============================

void cCtrlNewmark :: DampingMatrix(cSysMatrix *M, cSysMatrix *K, cSysMatrix *C)
{
  if (!Damping) return;

  // Compute the Rayleigh damping parameters.

  double alphar = 2*Omega1*Omega2*(Xi1*Omega2-Xi2*Omega1)/(Omega2*Omega2-Omega1*Omega1);
  double betar  = 2*(Xi2*Omega2-Xi1*Omega1)/(Omega2*Omega2-Omega1*Omega1);
//  cout << scientific;
//  cout << "alphar = " << alphar  << "  betar = " << betar << endl; 

  // Compute the damping matrix [C] = alphar*[M] + betar*[K].

  C->Zero( );
  C->AddMat(alphar, M);   // [C] = alphar*[M]
  C->AddMat(betar, K);    // [C] = alphar*[M] + betar*[K]
}

// ============================= PrintVelAcc  ==============================

void cCtrlNewmark ::PrintVelAcc(cVector &v, cVector &a)
{
  int nnode = cNode :: GetNumNode( );

  // Print the nodal velocities.

  out << "%RESULT.CASE.STEP.VELOCITY\n";
  out << nnode << "  'velocity'\n";
  for (int i = 0; i < nnode; i++)
  {
    // Get the current node.

    cNode *node = cNode :: GetNode(i);
    out << left << setw(5) << node->GetLabel( ) << " ";

    // Get the nodal velocities.

    out << showpos << scientific << setprecision(OutPrec);
    for (int j = 0; j < 6; j++)
    {
      int dof  = node->GetDof(j);
      double val = 0.0;
      if (dof > 0) val = v[dof-1];
      out << val << "  ";
    }
    out << "\n" << noshowpos;
  }
  out << "\n";

  // Print the nodal accelerations.

  out << "%RESULT.CASE.STEP.ACCELERATION\n";
  out << nnode << "  'acceleration'\n";
  for (int i = 0; i < nnode; i++)
  {
    // Get the current node.

    cNode *node = cNode :: GetNode(i);
    out << left << setw(5) << node->GetLabel( ) << " ";

    // Get the nodal accelerations.

    out << showpos << scientific << setprecision(OutPrec);
    for (int j = 0; j < 6; j++)
    {
      int dof  = node->GetDof(j);
      double val = 0.0;
      if (dof > 0) val = a[dof-1];
      out << val << "  ";
    }
    out << "\n" << noshowpos;
  }
  out << "\n";
}

// -------------------------------------------------------------------------
// Class cCtrlIncNmk:
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Public methods:
//

// ============================= cCtrlIncNmk ===============================

cCtrlIncNmk :: cCtrlIncNmk(void)
{
  Type = NEWMARK_INCR;
}

// ============================= ~cCtrlIncNmk ==============================

cCtrlIncNmk :: ~cCtrlIncNmk(void)
{
}

// ================================ Solver =================================

void cCtrlIncNmk :: Solver(void)
{
  if (Feedback)
  {
    cout << "\n\tIncremental Dynamic Analysis - Newmark Method\n";
    cout << "\tGamma = " << Gamma << endl;
    cout << "\tBeta  = " << Beta  << endl;
  }

  // Create local vectors.

  cVector f(NumEq);            // External load vector
  cVector g(NumEq);            // Internal force vector
  cVector r(NumEq);            // Residual vector
  cVector u(NumEq);            // Displacement vector
  cVector du(NumEq);           // incremental displacements
  cVector v(NumEq);            // Velocity vector
  cVector a(NumEq);            // Acceleration vector
  cVector w(NumEq);            // Auxiliary vector
  cVector z(NumEq);            // Auxiliary vector

  // Initialize algorithm variables.

  f.Zero( );
  g.Zero( );
  r.Zero( );
  u.Zero( );
  v.Zero( );
  a.Zero( );
  w.Zero( );
  z.Zero( );
  du.Zero( );
  TotTime = 0.0;

  // Calculate integration constants.

  double a0 = 1.0/(Beta*TimeStep*TimeStep);
  double a1 = Gamma/(Beta*TimeStep);
  double a2 = 1.0/(Beta*TimeStep);
  double a3 = 0.5/Beta - 1.0;
  double a4 = Gamma/Beta - 1.0;
  double a5 = (0.5*Gamma/Beta - 1.0)*TimeStep;
  double a6 = (1.0 - Gamma)*TimeStep;
  double a7 = Gamma*TimeStep;
#if 1
  printf("a0 = %f\n", a0);
  printf("a1 = %f\n", a1);
  printf("a2 = %f\n", a2);
  printf("a3 = %f\n", a3);
  printf("a4 = %f\n", a4);
  printf("a5 = %f\n", a5);
  printf("a6 = %f\n", a6);
  printf("a7 = %f\n", a7);
#endif

  // Compute the mass matrix.

  cSysMatrix *M = new cSymSkylMatrix(NumEq, Profile);
  MassMatrix(M);
//  cout << "[M] = ";
//  M->Print( );

  // Compute the acceleration vector at t = 0.0 (Solve [M]{a} = {r}).

  ExtLoadVector(TotTime, f);
  IntForceVector(g);
  r = f - g;
  cSysMatrix *K = new cSymSkylMatrix(NumEq, Profile);
  K->Zero( );
  K->AddMat(1.0, M);       // [K] = [M]
  K->Solve(r, a);

  // Compute the stiffness and damping matrices.

  StiffnessMatrix(K);
//  cout << "[K] = ";
//  K->Print( );
  cSysMatrix *C = 0;
  if (Damping)
  {
    C = new cSymSkylMatrix(NumEq, Profile);
    DampingMatrix(M, K, C);
//  cout << "[C] = ";
//  C->Print( );
  }

  // Compute [K]eff = [K] + a0*[M] + a1*[C].

  K->AddMat(a0, M);
  if (Damping) K->AddMat(a1, C);
//  cout << "[K]eff = ";
//  K->Print( );

  // Incremental loop.

  double umax = 0.0;
  double vmax = 0.0;
  double amax = a.NormInf( );
  for (CurrStep = 1; CurrStep <= MaxStep; CurrStep++)
  {
    // Update the total time and compute the external forces.

    TotTime += TimeStep;
    TotFactor = TotTime;   // For printing in output file
    if (Feedback)
    {
      cout << "\n\tStep = " << setw(3)  << CurrStep;
      cout << " Time = " << scientific << setprecision(3);
      cout << TotTime << " ............." << endl << resetiosflags(ios::scientific);
    }
    ExtLoadVector(TotTime, f);
//    cout << "{f} = " << scientific << setprecision(4) << showpos;
//    f.Print( );

    // Compute the internal and residual force vectors.

    IntForceVector(g);
    r = f - g;

    // Add the inertia terms: {r} += [M]*(a2*{v} + a3*{a}).

    w  = a2*v;
    w += a3*a;                        // {w} = a2*{v} + a3*{a}
    M->MultVect(w.Val( ), z.Val( ));  // {z} = [M]{w}
    r += z;

    // Add the damping effect: {r} += [C]*(a4*{v} + a5*{a}).

    if (Damping)
    {
      w  = a4*v;
      w += a5*a;                        // {w} = a4*{v} + a5*{a}
      C->MultVect(w.Val( ), z.Val( ));  // {z} = [C]{w}
      r += z;
    }

    // Compute the incremental displacements (Solve [K]{du} = {r}).

    K->Solve(r, du);
    u += du;

    // Compute the new acceleration vector.

    w  = a;          // Store the old acceleration
    a *= -a3;
    a +=  a0*du;
    a += -a2*v;

    // Compute the new velocity vector.

    v += a6*w;       // Term depeding on the old acceleration
    v += a7*a;       // Term depeding on the new acceleration

    // Update the nodal displacements for stress computations.

    AssignDispl(u);

    // Print the computed results.

    umax = MAX(umax, u.NormInf( ));
    vmax = MAX(vmax, v.NormInf( ));
    amax = MAX(amax, a.NormInf( ));
    if (CurrStep%PrintStep == 0)
    {
      if (Feedback) cout << "\n\tPrinting computed results ..............." << endl;
      PrintResult( );
      PrintVelAcc(v, a);
    }

    // Update element variables.

    UpdateState( );
  }

  if (Feedback)
  {
    cout << "umax = " << scientific << setprecision(3) << umax << endl;
    cout << "vmax = " << scientific << setprecision(3) << vmax << endl;
    cout << "amax = " << scientific << setprecision(3) << amax << endl;
  }

  // Adjust the current step.

  CurrStep--;
}


// -------------------------------------------------------------------------
// Class cCtrlIncNmk:
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Public methods:
//

// ========================== cCtrlNonlinNewmark ===========================

cCtrlNonlinNewmark :: cCtrlNonlinNewmark(void)
{
  Type = NEWMARK_NL;
}

// ========================== ~cCtrlNonlinNewmark ==========================

cCtrlNonlinNewmark :: ~cCtrlNonlinNewmark(void)
{
}

// ================================ Solver =================================

void cCtrlNonlinNewmark :: Solver(void)
{
  if (Feedback)
  {	  
    cout << "\n\tNonlinear Dynamic Analysis - Newmark Method\n";
    cout << "\tGamma = " << Gamma << endl;
    cout << "\tBeta  = " << Beta  << endl;
  }

  // Create the local vectors.

  cVector f(NumEq);            // External load vector
  cVector g(NumEq);            // Internal force vector
  cVector r(NumEq);            // Residual vector
  cVector u(NumEq);            // Displacement vector
  cVector du(NumEq);           // Displacement increments
  cVector Du(NumEq);           // Displacement increments each time increment
  cVector v(NumEq);            // Velocity vector
  cVector a(NumEq);            // Acceleration vector
  cVector w(NumEq);            // Auxiliary vector
  cVector z(NumEq);            // Auxiliary vector

  // Initialize local variables.

  f.Zero( );
  g.Zero( );
  r.Zero( );
  u.Zero( );
  v.Zero( );
  a.Zero( );
  w.Zero( );
  z.Zero( );
  du.Zero( );
  Du.Zero();
  TotTime = 0.0;

  // Calculate integration constants.

  double a0 = 1.0/(Beta*TimeStep*TimeStep);
  double a1 = Gamma/(Beta*TimeStep);
  double a2 = 1.0/(Beta*TimeStep);
  double a3 = 0.5/Beta - 1.0;
  double a4 = Gamma/Beta - 1.0;
  double a5 = (0.5*Gamma/Beta - 1.0)*TimeStep;
  double a6 = (1.0 - Gamma)*TimeStep;
  double a7 = Gamma*TimeStep;
#if 1
  printf("a0 = %f\n", a0);
  printf("a1 = %f\n", a1);
  printf("a2 = %f\n", a2);
  printf("a3 = %f\n", a3);
  printf("a4 = %f\n", a4);
  printf("a5 = %f\n", a5);
  printf("a6 = %f\n", a6);
  printf("a7 = %f\n", a7);
#endif

  // Create the stiffness matrix.

  cSysMatrix *K = new cSymSkylMatrix(NumEq, Profile);

  // Create the damping matrix.

  cSysMatrix *C = 0;
  if (Damping) C = new cSymSkylMatrix(NumEq, Profile);

  // Compute the mass matrix.

  cSysMatrix *M = new cSymSkylMatrix(NumEq, Profile);
  MassMatrix(M);

  // Compute the acceleration vector at t = 0.0 (Solve [M]{a} = {r}).

  ExtLoadVector(TotTime, f);
  IntForceVector(g);
  r = f - g;
  K->Zero( );
  K->AddMat(1.0, M);       // [K] = [M]
  K->Solve(r, a);

  // Incremental loop.

  double umax = 0.0;
  double vmax = 0.0;
  double amax = a.NormInf( );
  for (CurrStep = 1; CurrStep <= MaxStep; CurrStep++)
  {
    // Update the total time and initialize displacement increment {Du}.

    TotTime += TimeStep;
    TotFactor = TotTime;   // For printing in output file
    if (Feedback)
      printf("\n\tIncremental step %-4d ...................\n\n", CurrStep);

    Du.Zero( ); // The displacement increment in each time increment is zero

    // Compute the current external load vector.

    ExtLoadVector(TotTime, f);

    // Iterative loop.

    int conv = 0;
    for (CurrIter = 1; CurrIter <= MaxIter; CurrIter ++)
    {
      // Compute the effective stiffiness: [Keff] = [K] + a0*[M] + a1*[C].

      StiffnessMatrix(K);
      if (Damping)
      {
        DampingMatrix(M, K, C);  // Use the current [K]
        K->AddMat(a1, C);
      }
      K->AddMat(a0, M);

      // Compute the effective load vector: {r} = {f} - {g} + [M]{w}.

      r = f - g;
      w  = -a0*Du;                       
      w +=  a2*v;
      w +=  a3*a;                       // {w} = -a0*{Du} + a2*{v} + a3*{a}
      M->MultVect(w.Val( ), z.Val( ));  // {z} = [M]{w}
      r += z;

      // Add the damping effect: [C]*(-a1*{Du} + a4*{v} + a5*{a}).

      if (Damping)
      {
        w  = -a1*Du;
        w +=  a4*v;
        w +=  a5*a;                       // {w} = -a1*{Du} + a4*{v} + a5*{a}
        C->MultVect(w.Val( ), z.Val( ));  // {z} = [C]{w}
        r += z;
      }    

      // Compute the iterative correction {du} and update the displacements.

      K->Solve(r, du);
      u  += du;        // Total displacements
      Du += du;        // Incremental displacements

      // Evaluate the internal force vector.

      AssignDispl(u);
      if (!IntForceVector(g)) break;

      // Check convergence.

      conv = Convergence(f.Length(), r.Length(), TotTime);
      if (conv) break;
    }

    // Update the accelerations and velocities after convergence.

    if (conv)
    {
      // Compute the new accelerations: {a}n+1 = -a3*{a}n + a0*{Du} - a2*{v}n.

      w  = a;         // Store the old acceleration: {a}n
      a *= -a3;
      a +=  a0*Du;
      a += -a2*v;

      // Compute the new velocities:  {v}n+1 = {v}n + a6*{a}n + a7*{v}n+1.

      v += a6*w;      // Term depending on the old acceleration
      v += a7*a;      // Term depending on the new acceleration

      // Print the computed results.

      umax = MAX(umax, u.NormInf( ));
      vmax = MAX(vmax, v.NormInf( ));
      amax = MAX(amax, a.NormInf( ));
      if (CurrStep%PrintStep == 0)
      {
        if (Feedback) cout << "\n\tPrinting computed results ..............." << endl;
        PrintResult( );
        PrintVelAcc(v, a);
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

  if (Feedback)
  {
    cout << "umax = " << scientific << setprecision(3) << umax << endl;
    cout << "vmax = " << scientific << setprecision(3) << vmax << endl;
    cout << "amax = " << scientific << setprecision(3) << amax << endl;
  }

  // Adjust the current step.

  CurrStep--;
}

// ======================================================= End of file =====
