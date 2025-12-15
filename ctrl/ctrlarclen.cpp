// -------------------------------------------------------------------------
// ctrlarclen.cpp - implementation of the control class for nonlinear static
//                  problems (arc-Length with Newton-Raphson iterations).
// -------------------------------------------------------------------------
// Created:      22-Out-2011     Iuri Barcelos Rocha
//
// Modified:     08-Apr-2015     Elias Saraiva Barroso
//               Input/output operations using c++ streams.
// -------------------------------------------------------------------------

#include <stdlib.h>
#include <math.h>
#include <iostream>

using namespace std;

#include "ctrlarclen.h"
#include "node.h"
#include "vec.h"
#include "sysmat.h"
#include "gblvar.h"

// -------------------------------------------------------------------------
// cCtrlArcLength base class:
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Static variables
//

double cCtrlArcLength :: LoadCoeff = 0.0;
bool   cCtrlArcLength :: Restart   = false;

// -------------------------------------------------------------------------
// Public methods:
//

// =============================== ReadData ================================

void cCtrlArcLength :: ReadData(void)
{
  if (!(in >> LoadCoeff))
  {
    cout << "Error in the input of the arc-length method control data!\n";
    exit(0);
  }
}

// ============================ cCtrlArcLength =============================

cCtrlArcLength :: cCtrlArcLength(void)
{
}

// =========================== ~cCtrlArcLength =============================

cCtrlArcLength :: ~cCtrlArcLength(void)
{
}

// -------------------------------------------------------------------------
// Protected methods:
//

// ============================= SetArcLength ==============================

void cCtrlArcLength :: SetArcLength(cVector &du1)
{
  ArcLen = StepFactor*(du1.Length( ) + LoadCoeff);
}

// ============================ FirstIncLoadFac ============================

double cCtrlArcLength :: FirstIncLoadFac(cVector &du1, cVector &Du)
{
  if ((LoadCoeff == 1.0) && ((fabs(CSP1) < 0.3) || (Du.Length() > Tol)))
  {
    cout << "ArcLen = " << ArcLen << "\n";
    ArcLen = Du.Length( );
    LoadCoeff = 0.0;
    cout << "ArcLen = " << ArcLen << "\n";
  }
  double dlf = ArcLen/(du1.Length( ) + LoadCoeff);

  if (du1*Du < 0.0)
    return(-dlf);
  else
    return(dlf);
}

// -------------------------------------------------------------------------
// cIniOrtArcLength class:
// -------------------------------------------------------------------------

// =========================== cIniOrtArcLength ============================

cIniOrtArcLength :: cIniOrtArcLength(void)
{
  PathType = INIORTAL;
}

// ========================== ~cIniOrtArcLength ============================

cIniOrtArcLength :: ~cIniOrtArcLength(void)
{
}

// =============================== IncLoadFac ==============================

double cIniOrtArcLength :: IncLoadFac(cVector &du1, cVector &du2,
                                      cVector &Du)
{
  double dlf;
  static cVector Du0;

  if (CurrStep == 1 && CurrIter == 1)
  {
    SetArcLength(du1);
    dlf = StepFactor;
    Du0.Resize(NumEq);
    Du0 = dlf*du1;
  }
  else if (CurrIter == 1)
  {
    if (TargetIter)
      ArcLen *= sqrt(double(TargetIter)/double(PrevIter));
    dlf = FirstIncLoadFac(du1, Du);
    Du0 = dlf*du1;
  }
  else
  {
    dlf = (Du0*du2)/(Du0*du1 + LoadCoeff*Dlf);
  }

  // Update load increment in current step

  if (CurrIter == 1) Dlf = 0.0;
  Dlf += dlf;

  return(dlf);
}


// -------------------------------------------------------------------------
// cUpOrtArcLength class:
// -------------------------------------------------------------------------

// =========================== cUpOrtArcLength =============================

cUpOrtArcLength :: cUpOrtArcLength(void)
{
  PathType = UPDORTAL;
}

// ========================== ~cUpOrtArcLength =============================

cUpOrtArcLength :: ~cUpOrtArcLength(void)
{
}

// =============================== IncLoadFac ==============================

double cUpOrtArcLength :: IncLoadFac(cVector &du1, cVector &du2,
                                     cVector &Du)
{
  double dlf;

  if (CurrStep == 1 && CurrIter == 1)
  {
    SetArcLength(du1);
    dlf = StepFactor;
  }
  else if (CurrIter == 1)
  {
    if (TargetIter)
      ArcLen *= sqrt(double(TargetIter)/double(PrevIter));
    dlf = FirstIncLoadFac(du1, Du);
  }
  else
  {
    dlf = (Du*du2)/(Du*du1 + LoadCoeff*Dlf);
  }

  // Update load increment in current step

  if (CurrIter == 1) Dlf = 0.0;
  Dlf += dlf;

  return(dlf);
}



// -------------------------------------------------------------------------
// cConLinArcLength class:
// -------------------------------------------------------------------------

// =========================== cConLinArcLength ============================

cConLinArcLength :: cConLinArcLength(void)
{
  PathType = CONLINAL;
}

// ========================== ~cConLinArcLength ============================

cConLinArcLength :: ~cConLinArcLength(void)
{
}

// =============================== IncLoadFac ==============================

double cConLinArcLength :: IncLoadFac(cVector &du1, cVector &du2,
                                      cVector &Du)
{
  double dlf;

  if (CurrStep == 1 && CurrIter == 1)
  {
    SetArcLength(du1);
    dlf = StepFactor;
  }
  else if (CurrIter == 1)
  {
    if (TargetIter)
      ArcLen *= sqrt(double(TargetIter)/double(PrevIter));
    dlf = FirstIncLoadFac(du1, Du);
  }
  else
  {
    double r0 = Du*Du + LoadCoeff*Dlf*Dlf - ArcLen*ArcLen;
    dlf = (Du*du2 - 0.5*r0)/(Du*du1 + LoadCoeff*Dlf);
  }

  // Update load increment in current step

  if (CurrIter == 1) Dlf = 0.0;
  Dlf += dlf;

  return(dlf);
}


// -------------------------------------------------------------------------
// cCilArcLength class:
// -------------------------------------------------------------------------


// =========================== cCilArcLength ===============================

cCilArcLength :: cCilArcLength(void)
{
  PathType = CILAL;
}

// ========================== ~cCilArcLength ===============================

cCilArcLength :: ~cCilArcLength(void)
{
}

// ============================ IncLoadFac =================================

double cCilArcLength :: IncLoadFac(cVector &du1, cVector &du2, cVector &Du)
{
  double dlf;

  if (CurrStep == 1 && CurrIter == 1 && !Restart)
  {
    SetArcLength(du1);
    dlf = StepFactor;
  }
  else if (CurrIter == 1)
  {
    if (Restart)
    {
      ArcLen /= 2.0;
      Restart = false;
    }
    else if (TargetIter)
      ArcLen *= sqrt(double(TargetIter)/double(PrevIter));
    dlf = FirstIncLoadFac(du1, Du);
  }
  else
  {
    cVector v(NumEq);
    v = Du - du2;

    double a = du1*du1 + LoadCoeff;
    double b = 2.0*(v*du1) + 2.0*LoadCoeff*Dlf;
    double c = v*v - ArcLen*ArcLen + LoadCoeff*Dlf*Dlf;
    double d = b*b - 4.0*a*c;

    if (d <= 0.0)
    {
      CurrIter = 0;  // Restart iterative process
      Restart = true;
      cout << "\nComplex roots in cilindrical arc length. Restarting with halved arc length...\n";
      return(0.0);
    }

    double r1 = (-b + sqrt(d))/(2.0*a);
    double r2 = (-b - sqrt(d))/(2.0*a);

    // Choose the root

    dlf = r1;
    double e = Du*v;
    double f = Du*du1;
    double t1 = e + f*r1;
    double t2 = e + f*r2;
    if (t2 > t1) dlf = r2;
  }

  // Update load increment in current step

  if (CurrIter == 1) Dlf = 0.0;
  Dlf += dlf;

  return(dlf);
}

// ======================================================= End of file =====
