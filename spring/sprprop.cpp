// -------------------------------------------------------------------------
// sprprop.cpp - implementation of the spring property class.
// -------------------------------------------------------------------------
// Created:      14-Ago-2011     Evandro Parente Junior
//
// Modified:     08-Apr-2015     Elias Saraiva Barroso
//               Input/output operations using c++ streams.
// -------------------------------------------------------------------------

#include <iostream>
#include <math.h>
#include <stdlib.h>

using namespace std;

#include "sprprop.h"
#include "utl.h"
#include "gblvar.h"
#include "gbldef.h"

// -------------------------------------------------------------------------
// Static variables:
//
int           cSpringProp :: NumSprProp = 0;
cSpringProp** cSpringProp :: VecSprProp = 0;

// -------------------------------------------------------------------------
// Public methods:
//

// ============================== ReadSpringProp ===========================

void cSpringProp :: ReadSpringProp(void)
{
  // Read the number of spring properties.

  if (!(in >> NumSprProp) || (NumSprProp < 1))
  {
    cout << "Error in the input of the number of spring properties!\n";
    exit(0);
  }

  // Alloc the array of spring properties.

  VecSprProp = new cSpringProp*[NumSprProp];

  // Read each spring property.

  int id;
  char type[256];
  for (int i = 0; i < NumSprProp; i++)
  {
    if (!(in >> id) || !(in >> type))
    {
      cout << "Error in the input of the spring property " << i+1 << "!\n";
      exit(0);
    }

    if ((string(type) == "Linear") || (string(type) ==  "linear") ||
        (string(type) == "LINEAR"))
    {
      VecSprProp[i] = new cSpringPropLinear(id);
    }
    else
    {
      VecSprProp[i] = new cSpringPropNonlin(id);
    }
    VecSprProp[i]->Read( );
  }
}

// ================================ Destroy ================================

void cSpringProp :: Destroy(void)
{
  // Destroy each property.

  for (int i = 0; i < NumSprProp; i++) delete VecSprProp[i];

  // Release the array of properties.

  delete []VecSprProp;
}

// ================================ Destroy ================================

cSpringProp* cSpringProp :: GetSpringProp(int label)
{
  // Check the default position.

  if (VecSprProp[label-1]->Label == label) return(VecSprProp[label-1]);

  // Search the given label in the whole vector.

  for (int i = 0; i < NumSprProp; i++)
    if (VecSprProp[i]->Label == label) return(VecSprProp[i]);

  // Not found.

  return(0);
}

// ============================== cSpringProp ==============================

cSpringProp :: cSpringProp(void)
{
}

// ============================= ~cSpringProp ==============================

cSpringProp :: ~cSpringProp(void)
{
}


// -------------------------------------------------------------------------
// Class cSpringPropLinear:
// -------------------------------------------------------------------------


// =========================== cSpringPropLinear ===========================

cSpringPropLinear :: cSpringPropLinear(int label) : cSpringProp( )
{
  Label = label;
}

// ========================== ~cSpringPropLinear ===========================

cSpringPropLinear :: ~cSpringPropLinear(void)
{
}

// ================================= Read ==================================

void cSpringPropLinear :: Read(void)
{
  if (!(in >> K))
  {
    cout << "Error in the input of spring stiffness !\n";
    exit(0);
  }
}

// =============================== GetForce ================================

double cSpringPropLinear :: GetForce(double u)
{
  // Compute the elastic force (Hooke's law).

  double f = K*u;
  return(f);
}


// -------------------------------------------------------------------------
// Class cSpringPropNonlin:
// -------------------------------------------------------------------------

// =========================== cSpringPropNonlin ===========================

cSpringPropNonlin :: cSpringPropNonlin(int label) : cSpringProp( )
{
  Label = label;
}

// ========================== ~cSpringPropNonlin ===========================

cSpringPropNonlin :: ~cSpringPropNonlin(void)
{
  delete[]Displ;
  delete[]Force;
}

// ================================= Read ==================================

void cSpringPropNonlin :: Read(void)
{
  // Alloc the table points.

  if (!(in >> NumPnt) || NumPnt < 2)
  {
    cout << "Error in the input of nonlinear spring data !\n";
    exit(0);
  }
  Displ = new double[NumPnt];
  Force = new double[NumPnt];

  // Read points.

  if (!(in >> Displ[0]) || !(in >> Force[0]))
  {
    cout << "Error in the input of nonlinear spring data points!\n";
    exit(0);
  }
  for (int i = 1; i < NumPnt; i++)
  {
    if (!(in >> Displ[i]) || !(in >> Force[i]))
    {
      cout << "Error in the input of nonlinear spring data points!\n";
      exit(0);
    }
    if (Displ[i] < Displ[i-1])
    {
      cout << "Error in the input of nonlinear spring data points!\n";
      cout << "Points not in the ascending order!\n";
      exit(0);
    }
  }
}

// =============================== GetForce ================================

double cSpringPropNonlin :: GetForce(double u)
{
  // Find interval.

  int i;
  for (i = 1; i < NumPnt; i++) if (u <= Displ[i]) break;
  if (i == NumPnt) i--; // Extrapolation

  // Interpolate (or extrapolate) the function value.

  double du = Displ[i] - Displ[i-1];
  double df = Force[i] - Force[i-1];
  double Kt = 0.0;
  if (du != 0.0) Kt = df/du;
  double f = Force[i-1] + Kt*(u - Displ[i-1]);

  return(f);
}

// =============================== GetStiff ================================

double cSpringPropNonlin :: GetStiff(double u)
{
  // Find interval.

  int i;
  for (i = 1; i < NumPnt; i++) if (u <= Displ[i]) break;
  if (i == NumPnt) i--; // Extrapolation

  // Compute the tangent stiffness.

  double du = Displ[i] - Displ[i-1];
  double df = Force[i] - Force[i-1];
  double Kt = 0.0;
  if (du != 0.0) Kt = df/du;

  return(Kt);
}

// ======================================================= End of file =====
