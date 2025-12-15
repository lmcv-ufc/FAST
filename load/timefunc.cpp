// -------------------------------------------------------------------------
// timefunc.cpp - implementation of the time function class.
// -------------------------------------------------------------------------
// Created:      30-Jun-2005     Evandro Parente Junior
//
// Modified:     08-Apr-2015     Elias Saraiva Barroso
//               Input/output operations using c++ streams.
// -------------------------------------------------------------------------

#include <math.h>
#include <iostream>
#include <stdlib.h>

using namespace std;

#include "timefunc.h"
#include "gblvar.h"
#include "gbldef.h"

// -------------------------------------------------------------------------
// Static variables:
//
cTimeFunc* cTimeFunc :: Head = 0;
cTimeFunc* cTimeFunc :: Tail = 0;
cTimeFunc* cTimeFunc :: Curr = 0;

// -------------------------------------------------------------------------
// Public methods:
//

// =========================== ReadTimeFuncHarmonic ========================

void cTimeFunc :: ReadTimeFuncHarmonic(void)
{
  cFuncHarmonic *func = new cFuncHarmonic( );
  Curr = func;
  func->Read( );
}

// =========================== ReadTimeFuncHalfSine ========================

void cTimeFunc :: ReadTimeFuncHalfSine(void)
{
  cFuncHalfSine *func = new cFuncHalfSine( );
  Curr = func;
  func->Read( );
}

// ============================ ReadTimeFuncTable ==========================

void cTimeFunc :: ReadTimeFuncTable(void)
{
  cFuncTable *func = new cFuncTable( );
  Curr = func;
  func->Read( );
}

// ============================ ReadTimeFuncConst ==========================

void cTimeFunc :: ReadTimeFuncConst(void)
{
  cFuncConst *func = new cFuncConst( );
  Curr = func;
  func->Read( );
}

// ================================ Destroy ================================

void cTimeFunc :: Destroy(void)
{
  cTimeFunc *func = Head;

  while (func)
  {
    cTimeFunc *f = func->Next;  // Save the next element
    delete func;                // Delete the current element
    func = f;
  }
}

// ============================== cTimeFunc ================================

cTimeFunc :: cTimeFunc(void)
{
  this->Next = 0;
  if (Head == 0)    // First element
    Head = this;
  else              // Other elements
    Tail->Next = this;
  Tail = this;
}

// ============================= ~cTimeFunc ================================

cTimeFunc :: ~cTimeFunc(void)
{
}

// -------------------------------------------------------------------------
// Class cFuncConst:
// -------------------------------------------------------------------------


// ============================== cFuncConst ===============================

cFuncConst :: cFuncConst(void)
{
}

// ============================== cFuncConst ===============================

cFuncConst :: cFuncConst(double f)
{
  Val = f;
}

// ============================= ~cFuncConst ===============================

cFuncConst :: ~cFuncConst(void)
{
}

// ================================= Read ==================================

void cFuncConst :: Read(void)
{
  if (!(in >> Val))
  {
    cout << "Error in the input of constant time function !\n";
    exit(0);
  }
}


// -------------------------------------------------------------------------
// Class cFuncHarmonic:
// -------------------------------------------------------------------------

// ============================= cFuncHarmonic =============================

cFuncHarmonic :: cFuncHarmonic(void)
{
}

// ============================= cFuncHarmonic =============================

cFuncHarmonic :: cFuncHarmonic(double a, double p, double ph, double t0,
                               double t1)
{
  Amplit = a;
  Period = p;
  Phase  = ph;
  T0 = t0;
  T1 = t1;
}

// ============================ ~cFuncHarmonic ==============================

cFuncHarmonic :: ~cFuncHarmonic(void)
{
}

// ================================= Read ==================================

void cFuncHarmonic :: Read(void)
{
  if (!(in >> Amplit) || !(in >> Period) || !(in >> Phase) ||
      !(in >> T0) || !(in >> T1))
  {
    cout << "Error in the input of the harmonic time function !\n";
    exit(0);
  }
}

// ================================ GetVal =================================

double cFuncHarmonic :: GetVal(double t)
{
  // Check for time outside the valid interval.

  if (t < T0 || t > T1) return(0.0);

  // Compute the time function value.

  double f = Amplit*sin(2.0*PI*t/Period + Phase);
  return(f);
}


// -------------------------------------------------------------------------
// Class cFuncHalfSine:
// -------------------------------------------------------------------------

// ============================= cFuncHalfSine =============================

cFuncHalfSine :: cFuncHalfSine(void)
{
}

// ============================= cFuncHalfSine =============================

cFuncHalfSine :: cFuncHalfSine(double a, double p, double r)
{
  Amplit = a;
  Period = p;
  Rest   = r;
}

// ============================ ~cFuncHarmonic ==============================

cFuncHalfSine :: ~cFuncHalfSine(void)
{
}

// ================================= Read ==================================

void cFuncHalfSine :: Read(void)
{
  if (!(in >> Amplit) || !(in >> Period) || !(in >> Rest))
  {
    cout << "Error in the input of the half-sine time function !\n";
    exit(0);
  }
}

// ================================ GetVal =================================

double cFuncHalfSine :: GetVal(double t)
{
  // Check for time outside the valid interval.

  if (t < 0.0) return(0.0);

  // Reduce the given time to the function interval.

  double tt = Period/2.0 + Rest;
//  while (t > tt) t -= tt;
  int n = int(t/tt);
  t -= n*tt;

  // Compute the time function value.

  double f = 0.0;
  if (t < Period/2.0) f = Amplit*sin(2.0*PI*t/Period);
  return(f);
}

// -------------------------------------------------------------------------
// Class cFuncTable:
// -------------------------------------------------------------------------

// ============================== cFuncTable ===============================

cFuncTable :: cFuncTable(void)
{
}

// ============================== cFuncTable ===============================

cFuncTable :: cFuncTable(int n, double *t, double *f)
{
  if (n < 2 || t[0] != 0.0)
  {
    cout << "Error in the creation of table time function !\n";
    exit(0);
  }

  NumPnt = n;
  Pnt = new sTablePnt[NumPnt];
  for (int i = 0; i < NumPnt; i++)
  {
     Pnt[i].t = t[i];
     Pnt[i].f = f[i];
  }
}

// ============================= ~cFuncTable ===============================

cFuncTable :: ~cFuncTable(void)
{
}

// ================================= Read ==================================

void cFuncTable :: Read(void)
{
  // Alloc the table points.

  if (!(in >> NumPnt) || NumPnt < 2)
  {
    cout << "Error in the input of table time function !\n";
    exit(0);
  }
  Pnt = new sTablePnt[NumPnt];

  // Read points.

  for (int i = 0; i < NumPnt; i++)
  {
    if (!(in >> Pnt[i].t) || !(in >> Pnt[i].f))
    {
      cout << "Error in the input of function points!\n";
      exit(0);
    }
  }

  // Check first point.

  if (Pnt[0].t != 0.0)
  {
    cout << "Error: Initial time should be zero!\n";
    exit(0);
  }
}

// ================================ GetVal =================================

double cFuncTable :: GetVal(double t)
{
  // Find interval.

  int i;
  for (i = 1; i < NumPnt; i++) if (t <= Pnt[i].t) break;
  if (i == NumPnt) i--; // Extrapolation

  // Interpolate (or extrapolate) the function value.

  double dt = Pnt[i].t - Pnt[i-1].t;
  double df = Pnt[i].f - Pnt[i-1].f;
  double dfdt = 0.0;
  if (dt != 0.0) dfdt = df/dt;
  double f = Pnt[i-1].f + dfdt*(t - Pnt[i-1].t);

  return(f);
}

// ======================================================= End of file =====
