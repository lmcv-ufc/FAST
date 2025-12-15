// -------------------------------------------------------------------------
// ctrlgenalpha.cpp - implementaton of the control class for dynamic
// analysis using the Generalized-Alpha implicit algorithm.
// -------------------------------------------------------------------------
// Created:      01-Dez-2017     Bergson da Silva Matias
//
// Modified:     29-Feb-2020     Bergson Matias/Evandro Parente
//               Consideration of damping and new integration constants.
// -------------------------------------------------------------------------

#include <cstdlib>
#include <iostream>
#include <iomanip>

using namespace std;

#include "ctrlgenalpha.h"
#include "vec.h"
#include "sysmat.h"
#include "gbldef.h"
#include "gblvar.h"

// -------------------------------------------------------------------------
// Static variables:
//
double cCtrlGenAlpha :: Alpham =  0.0; // Default: no numerical damping
double cCtrlGenAlpha :: Alphaf =  0.0; // Default: no numerical damping

// -------------------------------------------------------------------------
// Public methods:
//

// ========================== ReadGenAlphaParam ============================

void cCtrlGenAlpha :: ReadGenAlphaParam(void)
{
  if (!(in >> Alpham) || !(in >> Alphaf))
  {
    cout << "Error in the input of Generalized-Alpha parameters!\n";
    exit(0);
  }
 
  if (Alpham < 0.0 || Alphaf < 0.0)
  {
    cout << "Error: Generalized-Alpha parameters must be non-negative!\n";
    exit(0);
  }
}

// ============================ cCtrlGenAlpha ==============================

cCtrlGenAlpha :: cCtrlGenAlpha(void)
{
  Type = GEN_ALPHA;
}

// ============================ ~cCtrlGenAlpha =============================

cCtrlGenAlpha :: ~cCtrlGenAlpha(void)
{
}

// ================================ Solver =================================

void cCtrlGenAlpha :: Solver(void)
{
  // Create the local vectors.

  cVector f(NumEq);            // External load vector
  cVector g(NumEq);            // Internal force vector
  cVector r(NumEq);            // Effective load vector
  cVector f_t(NumEq);          // External load vector at time t
  cVector g_t(NumEq);          // Internal force vector at time t
  cVector u(NumEq);            // Displacement vector
  cVector du(NumEq);           // Incremental displacements
  cVector v(NumEq);            // Velocity vector
  cVector a(NumEq);            // Acceleration vector
  cVector w(NumEq);            // Auxiliary vector
  cVector z(NumEq);            // Auxiliary vector

  // Initialize local  variables.

  f.Zero( );
  g.Zero( );
  r.Zero ( );
  f_t.Zero( );
  g_t.Zero( );
  u.Zero( );
  du.Zero( );
  v.Zero( );
  a.Zero( );
  w.Zero( );
  TotTime = 0.0;

  // Set gamma and beta values.

  double gamma = 0.5 - Alpham + Alphaf;
  double beta  = 0.25*(1.0 - Alpham + Alphaf)*(1.0 - Alpham + Alphaf);
  if (Feedback)
  {
    cout << "\n\tLinear Dynamic Analysis - Generalized-Alpha Method\n";
    cout << "\tAlpham = " << Alpham << endl;
    cout << "\tAlphaf = " << Alphaf << endl;
    cout << "\tGamma  = " << gamma << endl;
    cout << "\tBeta   = " << beta  << endl;
  }

  // Calculate integration constants.

  double a0 = 1.0/(beta*TimeStep*TimeStep);
  double a1 = gamma/(beta*TimeStep);
  double a2 = 1.0/(beta*TimeStep);
  double a3 = 0.5/beta - 1.0;
  double a4 = gamma/beta - 1.0;
  double a5 = (0.5*gamma/beta - 1.0)*TimeStep;
  double a6 = (1.0 - gamma)*TimeStep;
  double a7 = gamma*TimeStep;
  double b1 = (1.0 - Alphaf)*a1;
  double b4 = (1.0 - Alphaf)*a4 - Alphaf;
  double b5 = (1.0 - Alphaf)*a5;
  double c0 = (1.0 - Alpham)*a0;
  double c2 = (1.0 - Alpham)*a2;
  double c3 = (1.0 - Alpham)*a3 - Alpham;
#if 1
  printf("a0 = %f\n", a0);
  printf("a1 = %f\n", a1);
  printf("a2 = %f\n", a2);
  printf("a3 = %f\n", a3);
  printf("a4 = %f\n", a4);
  printf("a5 = %f\n", a5);
  printf("a6 = %f\n", a6);
  printf("a7 = %f\n", a7);
  printf("b1 = %f\n", b1);
  printf("b4 = %f\n", b4);
  printf("b5 = %f\n", b5);
  printf("c0 = %f\n", c0);
  printf("c2 = %f\n", c2);
  printf("c3 = %f\n", c3);
#endif

  // Compute the mass matrix.

  cSysMatrix *M = new cSymSkylMatrix(NumEq, Profile);
  MassMatrix(M);
  //  printf("[M] = ");
  //  M->Print( );

  // Compute the acceleration vector at t = 0.0 (Solve [M]{a} = {f} - {g}).

  ExtLoadVector(TotTime, f);  // f(t = 0)
  f_t = f;
  IntForceVector(g);
  g_t = g;
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

  // Compute [K]eff = (1 - Alphaf)*[K] + c0*[M] + b1*[C].

  K->AddMat(-Alphaf, K);
  K->AddMat(c0, M);
  if (Damping) K->AddMat(b1, C);
  //  printf("[K]eff = ");
  //  K->Print( );

  // Incremental loop.

  double umax = 0.0;
  double vmax = 0.0;
  double amax = a.NormInf( );
  for (CurrStep = 1; CurrStep <= MaxStep; CurrStep++)
  {
    // Update the total time and the external force vector.

    TotTime += TimeStep;       // t = t + Dt
    TotFactor = TotTime;       // For printing in output file
    if (Feedback)
    printf("\n\tStep = %-4d Time = %0.4e ...........\n", CurrStep, TotTime);
    ExtLoadVector(TotTime, f); // f(t + Dt)
    // cout << "{f} = " << scientific << setprecision(4) << showpos;
    // f.Print( );

    // Compute the effective force vector {r}.

    r = (1.0 - Alphaf)*f;
    r +=  Alphaf*f_t;
    r += -Alphaf*g_t;

    // Add the inertia effect: {r} += [M]*(c0*{u} + c2*{v} + c3*{a}).

    w  = c0*u;
    w += c2*v;
    w += c3*a;                       // {w} = c0*{u} + c2*{v} + c3*{a}
    M->MultVect(w.Val( ), z.Val( )); // {z} = [M]{w}
    r += z;

    // Add the damping effect: {r} += [C]*(b1*{u} + b4*{v} + b5*{a}).

    if (Damping)
    {
      w  = b1*u;
      w += b4*v;
      w += b5*a;                        // {w} = b1*{u} + b4*{v} + b5*{a}
      C->MultVect(w.Val( ), z.Val( ));  // {z} = [C]{w}
      r += z;
    }
    // cout << "{r} = " << scientific << setprecision(4) << showpos;
    // r.Print( );

    // Compute the displacement vector (Solve [K]{u} = {r}).

    w = u;           // Store the old displacements
    K->Solve(r, u);
    du = u - w;      // Compute the displacement increments
    // cout << "{u} = " << scientific << setprecision(4);
    // u.Print( );

    // Compute the new acceleration vector.

    w = a;           // Store the old acceleration
    a *= -a3;
    a +=  a0*du;
    a += -a2*v;
    // cout << "{a} = " << scientific << setprecision(4);
    // a.Print( );

    // Compute the new velocity vector.

    v += a6*w;       // Term depending on the old acceleration
    v += a7*a;       // Term depending on the new acceleration
    // cout << "{v} = " << scientific << setprecision(4);
    // v.Print( );

    // Update the displacements and compute the internal forces.

    AssignDispl(u);
    IntForceVector(g);
//    cout << "{g} = " << scientific << setprecision(4);
//    g.Print( );

    // Store the forces at the end of the current step.

    f_t = f;
    g_t = g;

    // Print the computed results.

    umax = MAX(umax, u.NormInf( ));
    vmax = MAX(vmax, v.NormInf( ));
    amax = MAX(amax, a.NormInf( ));
    if (CurrStep%PrintStep == 0)
    {
      if (Feedback) printf("\n\tPrinting computed results ...............\n");
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
// Class cCtrlNonlinGenAlpha:
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Public methods:
//

// ========================== cCtrlNonlinGenAlpha ==========================

cCtrlNonlinGenAlpha :: cCtrlNonlinGenAlpha(void)
{
  Type = GEN_ALPHA_NL;
}

// ========================== ~cCtrlNonlinGenAlpha =========================

cCtrlNonlinGenAlpha :: ~cCtrlNonlinGenAlpha(void)
{
}

// ================================ Solver =================================

void cCtrlNonlinGenAlpha :: Solver(void)
{
  // Create local vectors.

  cVector f(NumEq);            // External load vector
  cVector g(NumEq);            // Internal force vector
  cVector f_t(NumEq);          // External load vector at time t
  cVector g_t(NumEq);          // Internal force vector at time t
  cVector r(NumEq);            // Residual vector
  cVector u(NumEq);            // Displacement vector
  cVector du(NumEq);           // Displacement increments
  cVector Du(NumEq);           // Displacement increments each time increment
  cVector v(NumEq);            // Velocity vector
  cVector a(NumEq);            // Acceleration vector
  cVector w(NumEq);            // Auxiliary vector
  cVector z(NumEq);            // Auxiliary vector

  // Initialize algorithm variables.

  f.Zero( );
  g.Zero( );
  f_t.Zero( );
  g_t.Zero( );
  r.Zero( );
  u.Zero( );
  v.Zero( );
  a.Zero( );
  w.Zero( );
  z.Zero( );
  du.Zero( );
  Du.Zero( );
  TotTime = 0.0;

  // Set gamma and beta values.

  double gamma = 0.5 - Alpham + Alphaf;
  double beta  = 0.25*(1.0 - Alpham + Alphaf)*(1.0 - Alpham + Alphaf);
  if (Feedback)
  {
    cout << "\n\tNonlinear Dynamic Analysis - Generalized-Alpha Method\n";
    cout << "\tAlpham = " << Alpham << endl;
    cout << "\tAlphaf = " << Alphaf << endl;
    cout << "\tGamma  = " << gamma << endl;
    cout << "\tBeta   = " << beta  << endl;
  }

  // Calculate integration constants.

  double a0 = 1.0/(beta*TimeStep*TimeStep);
  double a1 = gamma/(beta*TimeStep);
  double a2 = 1.0/(beta*TimeStep);
  double a3 = 0.5/beta - 1.0;
  double a4 = gamma/beta - 1.0;
  double a5 = (0.5*gamma/beta - 1.0)*TimeStep;
  double a6 = (1.0 - gamma)*TimeStep;
  double a7 = gamma*TimeStep;
  double b1 = (1.0 - Alphaf)*a1;
  double b4 = (1.0 - Alphaf)*a4 - Alphaf;
  double b5 = (1.0 - Alphaf)*a5;
  double c0 = (1.0 - Alpham)*a0;
  double c2 = (1.0 - Alpham)*a2;
  double c3 = (1.0 - Alpham)*a3 - Alpham;
#if 1
  printf("a0 = %f\n", a0);
  printf("a1 = %f\n", a1);
  printf("a2 = %f\n", a2);
  printf("a3 = %f\n", a3);
  printf("a4 = %f\n", a4);
  printf("a5 = %f\n", a5);
  printf("a6 = %f\n", a6);
  printf("a7 = %f\n", a7);
  printf("b1 = %f\n", b1);
  printf("b4 = %f\n", b4);
  printf("b5 = %f\n", b5);
  printf("c0 = %f\n", c0);
  printf("c2 = %f\n", c2);
  printf("c3 = %f\n", c3);
#endif

  // Compute the stiffness matrix.

  cSysMatrix *K = new cSymSkylMatrix(NumEq, Profile);

  // Create the damping matrix.

  cSysMatrix *C = 0;
  if (Damping) C = new cSymSkylMatrix(NumEq, Profile);

  // Compute the mass matrix.

  cSysMatrix *M = new cSymSkylMatrix(NumEq, Profile);
  MassMatrix(M);

  // Compute the acceleration vector at t = 0.0 (Solve [M]{a} = {r}).

  ExtLoadVector(TotTime, f);   // f(t = 0)
  f_t = f;
  IntForceVector(g);
  g_t = g;
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

    Du.Zero( ); // The displacement increment before each time increment is zero

    // Compute the current external load vector.

    ExtLoadVector(TotTime, f);

    // Iterative loop.

    int conv = 0;
    for (CurrIter = 1; CurrIter <= MaxIter; CurrIter ++)
    {
      // Compute [K]eff = (1 - Alphaf)*[K] + c0*[M] + b1*[C].

      StiffnessMatrix(K);
      if (Damping) DampingMatrix(M, K, C);  // Use the current [K]
      K->AddMat(-Alphaf, K);
      K->AddMat(c0, M);
      if (Damping) K->AddMat(b1, C);
//      cout << "[K]eff = " << scientific << setprecision(3) << showpos;
//      K->Print( );

      // Compute {r} = (1 - Alphaf)*(f - g) + Alpha*(f_t - g_t).

      r  = f - g;
      r *= (1.0 - Alphaf);
      w  = f_t - g_t;
      r += Alphaf*w;

      // Add the inertia effect: {r} += [M]*(-c0*{Du} + c2*{v} + c3*{a}).

      w = -c0*Du;
      w += c2*v;
      w += c3*a;                       // {w} = -c0*{Du} + c2*{v} + c3*{a}
      M->MultVect(w.Val( ), z.Val( )); // {z} = [M]{w}
      r += z;

      // Add the damping effect: {r} += [C]*(-b1*{Du} + b4*{v} + b5*{a}).

      if (Damping)
      {
        w = -b1*Du;
        w += b4*v;
        w += b5*a;                        // {w} = -b1*{Du} + b4*{v} + b5*{a}
        C->MultVect(w.Val( ), z.Val( ));  // {z} = [C]{w}
        r += z;
      }
//      cout << "{r} = " << scientific << setprecision(3) << showpos;
//      r.Print( );

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
      // Store the force vectors at the end of the current time step.	    

      g_t = g;
      f_t = f;

      // Compute the new accelerations: {a}n+1 = -a3{a}n + a0{Du} - a2{v}n.

      w  = a;         // Store the old acceleration: {a}n
      a *= -a3;
      a +=  a0*Du;
      a += -a2*v;

      // Compute the new velocities:  {v}n+1 = {v}n + a6{a}n + a7{a}n+1.

      v += a6*w;      // Term depending on the old acceleration
      v += a7*a;      // Term depending on the new acceleration

      // Print the computed results.

      umax = MAX(umax, u.NormInf( ));
      vmax = MAX(vmax, v.NormInf( ));
      amax = MAX(amax, a.NormInf( ));
      if (CurrStep%PrintStep == 0)
      {
        if (Feedback) printf("\n\tPrinting computed results ...............\n");
        PrintResult( );
        PrintVelAcc(v, a);
      }
    }
    else
    {
      printf("Convergence not achieved!!!\n");
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
