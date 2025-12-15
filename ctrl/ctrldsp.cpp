// -------------------------------------------------------------------------
// ctrldsp.cpp - implementation of the control class for nonlinear static
//               problems (Displacement Control with Newton-Raphson
//               iterations).
// -------------------------------------------------------------------------
// Created:      06-Sep-2001     Evandro Parente Junior
//
// Modified:     03-Oct-2001     Evandro Parente Junior
//               Use of vector/matrix classes.
//
//               22-Oct-2011     Iuri Barcelos Rocha
//               Displacement control is now a subclass of cCtrlPath.
// -------------------------------------------------------------------------

#include <iostream>
#include <stdlib.h>

using namespace std;

#include "ctrldsp.h"
#include "node.h"
#include "vec.h"
#include "sysmat.h"
#include "gblvar.h"

// -------------------------------------------------------------------------
// Static variables:
//
int cCtrlDispl :: CtrlNode = 0;
int cCtrlDispl :: CtrlDir  = 0;
int cCtrlDispl :: CtrlDof  = 0;

// -------------------------------------------------------------------------
// Private methods:
//
// =============================== GetCtrlDof ==============================

void cCtrlDispl :: GetCtrlDof (void)
{
  // Get the controlled dof.

  cNode *node = cNode :: FindNode(CtrlNode);
  if (!node)
  {
    cout << "Invalid node for displacement control!\n";
    exit(0);
  }
  CtrlDof = node->GetDof(CtrlDir-1) - 1;
}

// -------------------------------------------------------------------------
// Public methods:
//

// =============================== ReadData ================================

void cCtrlDispl :: ReadData(void)
{
  if (!(in >> CtrlNode) || !(in >> CtrlDir))
  {
    cout << "Error in the input of the displacement control data!\n";
    exit(0);
  }

  if (CtrlDir < 1 || CtrlDir > 6)
  {
    cout << "Invalid controlled dof!\n";
    exit(0);
  }
}

// ============================== cCtrlDispl ===============================

cCtrlDispl :: cCtrlDispl(void)
{
  PathType = DCM;
}

// ============================== ~cCtrlDispl ==============================

cCtrlDispl :: ~cCtrlDispl(void)
{
}

// =============================== IncLoadFac ==============================

double cCtrlDispl :: IncLoadFac(cVector &du1, cVector &du2, cVector &Du)
{
  double dlf;

  if (CurrStep == 1 && CurrIter == 1) GetCtrlDof();

  if (CurrIter == 1)
        dlf = StepFactor/du1[CtrlDof];
  else
        dlf = du2[CtrlDof]/du1[CtrlDof];

  return(dlf);
}

// ======================================================= End of file =====
