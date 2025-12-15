// -------------------------------------------------------------------------
// ctrllins.cpp - implementation of the control class for lineat static
//               problems.
// -------------------------------------------------------------------------
// Created:      11-Dec-2019     Elias Saraiva Barroso
//
// -------------------------------------------------------------------------

#include <iostream>
#include <iomanip>

using namespace std;

#include "ctrlheat.h"
#include "vec.h"
#include "sysmat.h"
#include "gblvar.h"
#include "element.h"
#include "shape.h"
#include "shpiga.h"
#include "prescdof.h"
#include "heat.h"
#include "field.h"
#include "intpoint.h"

// -------------------------------------------------------------------------
// Public methods:
//

// =============================== cHeatTransf =============================

cHeatTransf :: cHeatTransf(void)
{
  Type = HEAT_TRANSFER;
}

// =============================== ~cHeatTransf ============================

cHeatTransf :: ~cHeatTransf(void)
{
}

// ============================= GlbHeatVector =============================

void cHeatTransf :: GlbHeatVector(cVector &f)
{
  // Initialize the heat flux vector.

  f.Zero( );

  // Compute and add the contribution of each element.

  int maxdof = cHeat :: GetMaxDof( );
  cVector felm(maxdof);
  for (cHeat *heat = cHeat::GetHead( ); heat != 0; heat = heat->GetNext( ))
  {
    heat->ExtHeat(1.0, felm);
    heat->AddGlobVec(felm, f);
  }

  // Modify the heat flux vector to consider the prescribed dofs.

  int npd = cPrescDof :: GetNumPrescDofs( );
  if (npd > 0)
  {
    // Get the prescribed dofs and values.

    int* prdof = new int[npd];
    cVector prval(npd);
    cPrescDof :: GetPrescDofs(prdof);
    cPrescDof :: GetPrescVals(1.0, prval);

    // Add the penalty force.

    for (int i = 0; i < npd; i++)
      f[prdof[i]-1] += PENALTY_STIFFNESS*prval[i];

    delete[] prdof;
  }
}

// ========================== ConductivityMatrix ===========================

void cHeatTransf :: ConductivityMatrix(cSysMatrix *K)
{
  int i;

  // Initialize the global stiffness matrix.

  K->Zero( );

  // Alloc the element conductivity matrix (max size).

  int maxdof = cElement :: GetMaxDof( );
  cMatrix Kelm(maxdof, maxdof);

  // Compute and add the contribution of each element.

  int nelm = cElement :: GetNumElm( );
  cElement* elm = 0;

  #pragma omp parallel for firstprivate(elm,Kelm)
  for (i = 0; i < nelm; i++)
  {
    // Get the current element.
    elm = cElement :: GetElm(i);

    // Evaluate the element conductivity matrix.
    elm->ConducMat(Kelm);

    // Add the element matrix to the global matrix.
    #pragma omp critical
    elm->AddGlobMat(Kelm, K);
  }

  // Modify the conductivity matrix to consider the prescribed dofs.

  int npd = cPrescDof :: GetNumPrescDofs( );
  if (npd > 0)
  {
    // Get array of prescribed dofs.

    int* prdof = new int[npd];
    cPrescDof :: GetPrescDofs(prdof);

    // Add the penalty conduction to the matrix diagonal.

    for (i = 0; i < npd; i++)
      K->Add(prdof[i]-1, prdof[i]-1, PENALTY_STIFFNESS);

    delete[] prdof;
  }
}

// ============================== AssignTemp ===============================

void cHeatTransf :: AssignTemp(cVector &t)
{
  int nnode = cNode :: GetNumNode( );
  for (int i = 0; i < nnode; i++)
  {
    // Get the current node.
    cNode *node = cNode :: GetNode(i);

    // Get the nodal temperature from the global vector.
    int dof = node->GetDof(7);
    if (dof)
      node->SetTemp(t[dof-1]);
    else
      node->SetTemp(0.0);
  }
}

// ================================ Solver =================================

void cHeatTransf :: Solver(void)
{
  if (Feedback) cout << "\n\tComputing the global heat vector ......" << endl;

  cVector f(NumEq);
  GlbHeatVector(f);

  if (Feedback) cout << "\n\tAssembling the conductivity matrix ........." << endl;

  cSysMatrix *K = cSysMatrix :: CreateMatrix(SysMatType, NumEq, Profile);
  K->SetParam(MaxIter,Tol);
  ConductivityMatrix(K);

  if (Feedback) cout << "\n\tSolving system equations ................" << endl;

  cVector t(NumEq);
  K->Solve(f, t);

  if (Feedback) cout << "\n\tPrinting computed results ..............." << endl;

  AssignTemp(t);
  PrintResult( );
}

// ============================== PrintResult ==============================

void cHeatTransf :: PrintResult(void)
{
  // Begin printing.

  Printing = true;

  // Print the result header.

  if (++PrintedSteps == 1)
  {
    out << "%RESULT\n";
    out << "1\n";
    out << "1  'result1'\n\n";
    
    out << "%RESULT.CASE" << endl; // Do not change endl to \n (use of tellp)
    CasePos[0] = out.tellp( );        // Store current position
    for (int i = 0; i <= 1.1*MaxStep; i++)
      out << "             \n";
    
    if (PrintEqvMsh)
    {
      outem << "%RESULT\n";
      outem << "1\n";
      outem << "1  'result1'\n\n";
    
      outem << "%RESULT.CASE" << endl; // Do not change endl to \n (use of tellp)
      CasePos[1] = outem.tellp( );        // Store current position
      for (int i = 0; i <= 1.1*MaxStep; i++)
        outem << "             \n";
    }
  }

  // Print the step header.

  out << "%RESULT.CASE.STEP\n";
  out << 1 << "\n\n";

  out << "%RESULT.CASE.STEP.FACTOR\n";
  out << scientific << setprecision(OutPrec) << TotFactor << "\n\n";

  if (PrintEqvMsh)
  {
    outem << "%RESULT.CASE.STEP\n";
    outem << 1 << "\n\n";

    outem << "%RESULT.CASE.STEP.FACTOR\n";
    outem << scientific << setprecision(OutPrec) << TotFactor << "\n\n";
  }

  // Print the nodal temperature.

  out << "%RESULT.CASE.STEP.NODAL.SCALAR\n";
  out << "1\n\'TEMPERATURE\'\n\n";

  out << "%RESULT.CASE.STEP.NODAL.SCALAR.DATA\n";
  int nnode = cNode :: GetNumNode( );
  out << nnode << "\n";
  for (int i = 0; i < nnode; i++)
  {
    // Get the current node.

    cNode *node = cNode :: GetNode(i);
    out << left << setw(5) << node->GetLabel( ) << " ";

    // Get the nodal displacements.

    out << showpos << scientific << setprecision(OutPrec);
    out << node->GetTemp( ) << "\n";
    out << resetiosflags(ios::showpos | ios::scientific);
  }
  out << "\n";

  // Print temperature L2 norm.
  if (cField :: GetInstance( )->hasTempField( ))
    PrintTempErrNorm( );

  // Print fluxes.
  PrintIntPntFlux( );

  // Print IGA equivament nodal temperature and fluxes.
  
  if (PrintEqvMsh)
  {
    PrintEqvMshTemp( );
    //PrintEqvMshNodalFluxes( );
  }

  // End printing.
  Printing = false;
}

// =========================== PrintIntPntFlux =============================

void cHeatTransf :: PrintIntPntFlux(void)
{
  int i,j,k;

  // Check if there are integration points.

  int maxp = cElement :: GetMaxIntPnt( );
  if (maxp == 0) return;

  // Print the integration point scalar labels.

  out << "%RESULT.CASE.STEP.ELEMENT.GAUSS.SCALAR\n";
  int  ngsl = cElement :: GetNumSclLab( );
  int *vgsl = new int[ngsl];
  cElement :: GetVecSclLab(vgsl);
  out << ngsl << "\n";
  for (i = 0; i < ngsl; i++) out << "'" << VEC_SCL_LAB[vgsl[i]] << "'  ";
  out << "\n\n";

  // Print the integration point scalar data.

  out << "%RESULT.CASE.STEP.ELEMENT.GAUSS.SCALAR.DATA\n";
  int nelm = cElement :: GetNumElm( );
  out << nelm << "\n";

  int *vel = new int[ngsl];
  int *idx = new int[ngsl];
  cVector prn(ngsl);
  cMatrix fluxes(maxp,ngsl);
  for (i = 0; i < nelm; i++)
  {
    // Get the current element.

    cElement *elm = cElement :: GetElm(i);
    out << elm->GetLabel( ) << "\n";
    int npt = elm->GetNumIntPnt( );
    out << npt << "\n";
    if (npt == 0) continue;

    // Get element stress labels and the respective indices.

    int nel = elm->GetNumStrCmp( );
    elm->GetStrLabels(vel);
    GetPrintIdx(nel, vel, ngsl, vgsl, idx);

    // Get integration point stresses.

    elm->IntPntFlux(fluxes);

    // Print stresses in the correct order.

    out << showpos << scientific << setprecision(OutPrec);

    for (j = 0; j < npt; j++)
    {
      prn.Zero( );
      for (k = 0; k < nel; k++) prn[idx[k]] = fluxes[j][k];
      for (k = 0; k < ngsl; k++) out << prn[k] << " ";
      out << "\n";
    }

    out << resetiosflags(ios::scientific | ios::showpos);
  }
  out << "\n";

  // Release memory.

  delete []vgsl;
  delete []vel;
  delete []idx;
}

// ============================ PrintErrNorm ===============================
/*
void cHeatTransf :: PrintTempErrNorm( )
{
  // Auxiliary variables.
  double norm = 0.0; // L2 error norm.
  double vol  = 0.0; // model volum.
  int i,j;           // Counters.
  sNodeCoord cart;   // Cartesian coordinates of integration points.

  // Check if there are integration points.

  int maxp = cElement :: GetMaxIntPnt( );
  if (maxp == 0) return;

  // Alloc memory.
  int maxnode       = cElement :: GetMaxNode( );
  int  ngsl         = cElement :: GetNumSclLab( );
  double *N         = new double[maxnode];
  sNodeCoord *coord = new sNodeCoord[maxnode];

  // Evaluate temperature L2 error norm in each element.
  cVector tvec(maxp);
  int nelm = cElement :: GetNumElm( );
  for (i = 0; i < nelm; i++)
  {
    // Get the current element.
    cElement *elm = cElement :: GetElm(i);

    // Get integration point temperature.
    elm->IntPntTemp(tvec);

    // Evaluate error function.
    sNatCoord pnt;
    double    w, temp;
    int       nip = elm->GetNumIntPnt( );
    cVector error(nip);
    error.Zero( );
    for(j = 0; j < nip; ++j)
    {
      // Get analytical solution.
      elm->GetIntPnt(j,pnt,w);
      elm->GetShape( )->Evaluate(pnt,N,coord,cart);
      temp = cField :: GetInstance( )->Temp(cart.x,cart.y,cart.z);

      error[j] = pow(tvec[j]-temp,2.0);
    }

    // Integrate error function.
    norm += elm->IntFunc(error,vol);
  }
  norm = sqrt(norm);

cout << scientific << setprecision(10) << endl;
  cout << "Temperature L2 norm: " << norm << " " << norm/vol << ", vol " << vol  << endl;
  out << "\n";
}
*/


void cHeatTransf :: PrintTempErrNorm( )
{
  int i,j,k;
  int ord            = ErrIntOrd;
  int NumIntPnt      = 0;
  int maxnode        = cElement :: GetMaxNode( );
  double norm        = 0.0;
  double vol         = 0.0;
  double *N          = new double[maxnode];
  eTopType ttype;
  cIntPoint *IntPnt  = 0;
  sNodeCoord *coord  = new sNodeCoord[maxnode];

   // Loop over each element.
  int nelm = cElement :: GetNumElm( );
  int ngsl = cElement :: GetNumSclLab( );
  cMatrix stress;
  cMatrix strain;
  for (i = 0; i < nelm; i++)
  {
    // Get the current element.
    cElement *elm = cElement :: GetElm(i);

    // Create quadrature points.
    if (i == 0 || ttype != elm->GetShape( )->GetTopologyType( ))
    {
      // Release memory of old integration points.
      if (IntPnt)  delete [] IntPnt;

      int order[3]      = {ord,ord,ord};
      ttype             = elm->GetShape( )->GetTopologyType( );
      eQuadType qtype   = static_cast<eQuadType>(elm->GetIntType( ));
      eTopType  stype   = elm->GetShape( )->GetTopologyType( );
                 IntPnt = cIntPoint::CreateIntPoints(stype,qtype,order,&NumIntPnt);
    }

    // Evaluate error function.
    cVector error(NumIntPnt);
    error.Zero( );
    cVector t(maxnode);
    for(j = 0; j < NumIntPnt; ++j)
    {
      // Evaluate numerical temperature result.
      t.Zero( );
      elm->NodalTemp(t);
      elm->GetShape( )->ShpFunc(IntPnt[j].GetCoord( ),N);
      double tp = 0;
      for(int i = 0; i < elm->GetShape( )->GetNumNode( ); ++i)
        tp += t[i]*N[i];

      // Get analytical solution.
      sNodeCoord pos;
      elm->GetShape( )->Evaluate(IntPnt[j].GetCoord( ),N,coord,pos);
      double temp = cField :: GetInstance( )->Temp(pos.x,pos.y,pos.z);

      // Compute error function.
      error[j] += pow((temp - tp),2.0);
    }

    // Integrate error function.
    norm += elm->IntFunc(NumIntPnt,IntPnt,error,vol);
  }
  norm = sqrt(norm);

  cout << scientific << setprecision(6);
  cout << "volume: " << vol << endl;
  cout << "L2 Temperature Norm: " << norm << endl; 
  out  << "%L2.TEMPERATURE.NORM\n" << norm << "\n\n";

  // Release memory.
  delete [] IntPnt;
}

void cHeatTransf :: PrintEqvMshTemp( )
{
  cShapeIGA :: PrintEqvMshTemp(outem); 
}

// ============================ PrintEqvMshNodalFluxes =====================

/* TODO: Must be done!
void cControl :: PrintEqvMshNodalFluxes( )
{
  int i,j,k,e;

  // Print the integration point scalar labels.

  outem << "%RESULT.CASE.STEP.ELEMENT.NODAL.SCALAR\n";
  int  maxen  = 20; // Eqv Elem with most nodes is brick20.
       maxen *= max(cShapeIGA :: GetEqvMshDisc(0)+1,1); // Multiplied by user
       maxen *= max(cShapeIGA :: GetEqvMshDisc(1)+1,1); // defined discretization.
       maxen *= max(cShapeIGA :: GetEqvMshDisc(2)+1,1);
  int  ngsl  = cElement :: GetNumSclLab( );
  int *vgsl  = new int[ngsl];
  cElement :: GetVecSclLab(vgsl);
  outem << ngsl << "\n";
  for (i = 0; i < ngsl; i++) outem << "'" << VEC_SCL_LAB[vgsl[i]] << "'  ";
  outem << "\n\n";

  int GblEqvElmID = 0;

  // Print the element nodal scalar data.
  const ElmStdVec &elmvec = cShapeIGA :: GetParentElmVec( );

  outem << "%RESULT.CASE.STEP.ELEMENT.NODAL.SCALAR.DATA\n";
  int nelm = elmvec.size( );
  outem << cShapeIGA :: GetNumEqvElm( ) << "\n";

  int *vel = new int[ngsl];
  int *idx = new int[ngsl];
  cVector prn(ngsl);
  cMatrix stress(maxen,ngsl);
  cMatrix strain(maxen,ngsl);
  for (i = 0; i < nelm; i++)
  {
    // Get element stress labels and the respective indices.

    int nel = elmvec[i]->GetNumStrCmp( );
    elmvec[i]->GetStrLabels(vel);
    GetPrintIdx(nel, vel, ngsl, vgsl, idx);

    // Get equivalent element type.
    eTopType t = elmvec[i]->GetShape( )->GetTopologyType( );

    // Get
    sEqvElmQuad *Quad = &(cShapeIGA :: GetEqvElmPnt(t));
    
    // Get all equivalent element stresses.  

    elmvec[i]->PntStress(&(Quad->coords[0]),Quad->NumPnt,strain,stress);

    // Print stresses in the correct order.

    int ne  = Quad->NumElm;
    for(e = 0; e < ne; ++e)
    {
      int nen = Quad->NumElmPnt;

      outem << ++GblEqvElmID << "\n";
      outem << showpos << scientific << setprecision(OutPrec);
      for (j = 0; j < nen; j++)
      {
        prn.Zero( );
        for (k = 0; k < nel; k++) prn[idx[k]] = stress[j][k];
        for (k = 0; k < ngsl; k++) outem << prn[k] << " ";
        outem << "\n";
      }
      outem << noshowpos;

      outem << resetiosflags(ios::scientific | ios::showpos);
    }
  }
  outem << "\n";

  // Release memory.

  delete []vgsl;
  delete []vel;
  delete []idx;
}
*/

// ======================================================= End of file =====
