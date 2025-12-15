// -------------------------------------------------------------------------
// elminterf.cpp - implementation of interface parametric elements.
// -------------------------------------------------------------------------
// Created:      16-Dez-2013     Edson Moreira Dantas Junior
//
// Modified:     07-May-2014     Evandro Parente Junior
//               Local-global transformation.
// -------------------------------------------------------------------------

#include <stdio.h>
#include <stdlib.h>
#include <iostream>

#include "elminterf.h"
#include "anmodel.h"
#include "node.h"
#include "shape.h"
#include "material.h"
#include "section.h"
#include "secanalysis.h"
#include "vec.h"
#include "mat.h"
#include "gblvar.h"
#include "gbldef.h"

// -------------------------------------------------------------------------
// Public methods:
//

// =============================== cElmInterf ===============================

cElmInterf :: cElmInterf(int id, eShpType tshp, eAnmType tanm):cElmParam(id,tshp, tanm)
{
  Type  = INTERFACE;
}

// ============================= ~cElmInterf ================================

cElmInterf :: ~cElmInterf(void)
{
}

// =============================== IntForce ================================

int cElmInterf :: IntForce(cVector &g)
{
  // Get the nodal coordinates.

  int nn = Shape->GetNumNode( );
  sNodeCoord *coord = new sNodeCoord[nn];
  Shape->NodalCoord(coord);

  // Get the nodal displacements.

  int ndof = nn*Model->GetNumDofNode( );
  cVector u(ndof);
  NodalDispl(u);

  // Create local arrays.

  int nsc = Model->GetDimBMatrix( );
  cMatrix B(nsc, ndof);
  cMatrix R(nsc, nsc);
  cVector el(nsc);   // Local strains
  cVector eg(nsc);   // Global strains
  cVector sl(nsc);   // Local stresses
  cVector sg(nsc);   // Global stresses

  // Get memory for shape functions.

  double *shpfnc = new double[nn];
  double *mapfnc = new double[nn];

  // Loop over integration points

  g.Zero( );
  int strok = 1;
  for (int i = 0; i < NumIntPnt; i++)
  {
    // Get integration point data.

    sNatCoord p = IntPnt[i].GetCoord( );
    double wgt  = IntPnt[i].GetWeight( );

    // Compute the global strains {eg} = [B]{u}.

    Shape->ShpFunc(p, shpfnc);
    Shape->MapFunc(p, mapfnc);
    Model->BMatrix(nn, shpfnc, mapfnc, coord, 0, B);
//    cout << "[B]\n";
//    B.Print( );
    eg = B*u;
//    cout << "{eg}\n";
//    eg.Print("%e ");

    // Compute the local stresses.

    double detJ;
    Shape->RMatrix(p, coord, &detJ, R);
    el = R*eg;

    // Compute the local stresses.

    if (!SecAn[i]->Stress(el, sl))
    {
      strok = 0;
      break;
    }

/*
    if (Label == 170)
    {
      cout << "[R]\n";
      R.Print( );
      cout << "{el}\n";
      el.Print("%e ");
      cout << "{sl}\n";
      sl.Print("%e ");
    }
*/

    // Compute the global stresses.

    sg = t(R)*sl;
//    cout << "{sg}\n";
//    sg.Print("%e ");

    // Compute [g] += coeff*t(B)*{sg}.

    double volc  = Model->VolCoeff(Section->GetThickness( ), nn, mapfnc, coord);
    double coeff = wgt*detJ*volc;
    MultTAcc(coeff, B, sg, g);
  }

  // Release memory.

  delete []coord;
  delete []shpfnc;
  delete []mapfnc;

  return(strok);
}

// =============================== StiffMat ================================

void cElmInterf :: StiffMat(cMatrix &K)
{
  // Get the nodal coordinates.

  int nn = Shape->GetNumNode( );
  sNodeCoord *coord = new sNodeCoord[nn];
  Shape->NodalCoord(coord);

  // Create local arrays.

  int ndof = nn*Model->GetNumDofNode( );
  int nsc = Model->GetDimBMatrix( );
  cMatrix Cl(nsc, nsc);
  cMatrix Cg(nsc, nsc);
  cMatrix R(nsc, nsc);
  cMatrix B(nsc, ndof);

  // Get memory for shape functions and derivatives.

  double *shpfnc = new double[nn];
  double *mapfnc = new double[nn];

  // Loop over integration points

  K.Zero( );
  for (int i = 0; i < NumIntPnt; i++)
  {
    // Get integration point data.

    sNatCoord p = IntPnt[i].GetCoord( );
    double wgt  = IntPnt[i].GetWeight( );

    // Compute [R] and [B] matrices.

    double detJ;
    Shape->ShpFunc(p, shpfnc);
    Shape->RMatrix(p, coord, &detJ, R);
    Model->BMatrix(nn, shpfnc, mapfnc, coord, 0, B);
//    cout << "[R]\n";
//    R.Print( );
//    cout << "[B]\n";
//    B.Print( );

    // Comput the global [C] matrix.

    SecAn[i]->CMatrix(Cl);        
/*
    if (Label == 170)
    {
      cout << "[Cl]\n";
      Cl.Print( );
    }
*/
    Cg.Zero( );                 // Do not change this line
    MatTripBtCB(R, Cl, 1.0, Cg);
//    cout << "[Cg]\n";
//    Cg.Print( );

    // Compute [K] += coeff*[B]t[Cg][B].

    double volc  = Model->VolCoeff(Section->GetThickness( ), nn, mapfnc, coord);
    double coeff = wgt*detJ*volc;
//    cout << "|J| = %f  volc = %f  coeff = %f\n", detJ, volc, coeff;
    MatTripBtCB(B, Cg, coeff, K);
  }
//  cout << "[K]\n";
//  K.Print( );

  // Release memory.

  delete []coord;
  delete []shpfnc;
  delete []mapfnc;
}

// =============================== MassMat =================================

void cElmInterf :: MassMat(cMatrix &M)
{
  M.Zero( );
}

// ============================== GeomStiff ================================

void cElmInterf :: GeomStiff(cMatrix &Kg)
{
  Kg.Zero( );
}

// ============================= NodalStress ===============================

void cElmInterf :: NodalStress(cMatrix &S)
{
  S.Zero( );
}

// ============================= IntPntStress ==============================

void cElmInterf :: IntPntStress(cMatrix &Strain, cMatrix &Stress)
{
  // Get the nodal coordinates.

  int nn = Shape->GetNumNode( );
  sNodeCoord *coord = new sNodeCoord[nn];
  Shape->NodalCoord(coord);

  // Get the nodal displacements.

  int ndof = nn*Model->GetNumDofNode( );
  cVector u(ndof);
  NodalDispl(u);

  // Create local arrays.

  int nsc = Model->GetDimBMatrix( );
  cMatrix B(nsc, ndof);
  cMatrix R(nsc, nsc);
  cVector el(nsc);   // Local strains
  cVector eg(nsc);   // Global strains
  cVector sl(nsc);   // Local stresses

  // Get memory for shape functions.

  double *shpfnc = new double[nn];
  double *mapfnc = new double[nn];

  // Loop over integration points.

  Strain.Zero( );
  Stress.Zero( );
  for (int i = 0; i < NumIntPnt; i++)
  {
    // Get integration point data.

    sNatCoord p = IntPnt[i].GetCoord( );

    // Compute the global strains {eg} = [B]{u}.

    Shape->ShpFunc(p, shpfnc);
    Shape->MapFunc(p, mapfnc);
    Model->BMatrix(nn, shpfnc, mapfnc, coord, 0, B);
    eg = B*u;

    // Compute the local stresses.

    double detJ;
    Shape->RMatrix(p, coord, &detJ, R);
    el = R*eg;

    // Compute the local stresses.

    SecAn[i]->Stress(el, sl);

    // Return the computed values and store the current stresses.

    for (int j = 0; j < nsc; j++)
    {
      Strain[i][j] = el[j];
      Stress[i][j] = StrIpt[i][j] = sl[j];
    }
  }

  // Release memory.

  delete []coord;
  delete []shpfnc;
  delete []mapfnc;
}

// ======================================================= End of file =====
