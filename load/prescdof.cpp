// -------------------------------------------------------------------------
// prescdof.cpp - implementation of the Prescribed Dof class.
// -------------------------------------------------------------------------
// Created:      23-Sep-2015     Evandro Parente Junior
//
// Modified:
// -------------------------------------------------------------------------

#include <iostream>
#include <math.h>
#include <stdlib.h>

using namespace std;

#include "prescdof.h"
#include "timefunc.h"
#include "node.h"
#include "vec.h"
#include "gblvar.h"
#include "gbldef.h"

// -------------------------------------------------------------------------
// Static variables:
//
int        cPrescDof :: NumPrescDof = 0;
cPrescDof* cPrescDof :: Head = 0;
cPrescDof* cPrescDof :: Tail = 0;

// -------------------------------------------------------------------------
// Public methods:
//

// ============================== ReadlPrescDofs ===========================

void cPrescDof :: ReadPrescDofs(void)
{
  // Read the number of prescribed dofs.

  int ndof;
  if (!(in >> ndof) || (ndof < 1))
  {
    cout << "Error in the input of the number of prescribed dofs!\n";
    exit(0);
  }

  // Read each dof.

  int id,dir;
  double val;
  for (int i = 0; i < ndof; i++)
  {
    if (!(in >> id) || !(in >> dir) || !(in >> val) || (dir < 1) || (dir > 8))
    {
      cout << "Error in the input of the prescribed dof " << i+1 << "!\n";
      exit(0);
    }

    cNode *node = cNode :: FindNode(id);
    if (!node)
    {
      cout << "Error in the input of the prescribed dof " << i+1 <<
              "(invalid node)!\n";
      exit(0);
    }

    // Create and store the prescribed dof.

    new cPrescDof(node, dir, val);
  }
}

// ============================== ReadlPrescTemps ==========================

void cPrescDof :: ReadPrescTemps(void)
{
  // Read the number of prescribed temperatures.

  int ndof;
  if (!(in >> ndof) || (ndof < 1))
  {
    cout << "Error in the input of the number of prescribed temperatures!\n";
    exit(0);
  }

  // Read each dof.

  int id;
  double val;
  for (int i = 0; i < ndof; i++)
  {
    if (!(in >> id) || !(in >> val))
    {
      cout << "Error in the input of the prescribed dof " << i+1 << "!\n";
      exit(0);
    }

    cNode *node = cNode :: FindNode(id);
    if (!node)
    {
      cout << "Error in the input of the prescribed temperature " << i+1 <<
              "(invalid node)!\n";
      exit(0);
    }

    // Create and store the prescribed dof.

    new cPrescDof(node, 8, val);
  }
}

// ============================= GetPrescDofs ==============================

void cPrescDof :: GetPrescDofs(int *prdof)
{
  int i = 0;
  for (cPrescDof *dof = Head; dof != 0; dof = dof->Next)
  {
    int dir = dof->Dir - 1;
    int idx = dof->Node->GetDof(dir);
    if (!idx)
    {
      cout << "Invalid prescribed dof at node ";
      cout << dof->Node->GetLabel( ) << "\n";
      exit(0);
    }
    prdof[i++] = idx;
  }
}

// ============================= GetPrescVals ==============================

void cPrescDof :: GetPrescVals(double t, cVector &prval)
{
  prval.Zero( );
  int i = 0;
  for (cPrescDof *dof = Head; dof != 0; dof = dof->Next)
  {
    prval[i++] = dof->Evaluate(t);
  }
}

// ============================ CurrNodalVals ==============================

void cPrescDof :: CurrNodalVals(cVector &crval)
{
  int i = 0;
  for (cPrescDof *dof = Head; dof != 0; dof = dof->Next)
  {
    int dir = dof->Dir - 1;
    crval[i++] = dof->Node->GetDofVal(dir);
  }
}

// ================================ Destroy ================================

void cPrescDof :: Destroy(void)
{
  cPrescDof *dof = Head;

  while (dof)
  {
    cPrescDof *d = dof->Next;  // Save the next element
    delete dof;                 // Delete the current element
    dof = d;
  }
}

// ============================== cPrescDof ================================

cPrescDof :: cPrescDof(cNode *node, int dir, double val)
{
  // Add to the dof list.

  this->Next = 0;
  if (Head == 0)    // First element
    Head = this;
  else              // Other elements
    Tail->Next = this;
  Tail = this;
  NumPrescDof++;

  // Store prescribed dof data.

  Node = node;
  Dir  = dir;
  DofVal = val;

  // Activate support flag.

  Node->SetSupp( );

  // Assign the current time function or use the default f(t) = 1.

  Func = cTimeFunc :: GetCurr( );
  if (!Func) Func = new cFuncConst(1.0);
}

// ============================= ~cPrescDof ================================

cPrescDof :: ~cPrescDof(void)
{
}

// =============================== Evaluate ================================

double cPrescDof :: Evaluate(double t)
{
  double tf = Func->GetVal(t);
  return(DofVal*tf);
}

// ======================================================= End of file =====
