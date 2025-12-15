// -------------------------------------------------------------------------
// elmparamC1.cpp - implementation of parametric elements with C1 continuity
// -------------------------------------------------------------------------
// Created:      14-Nov-2023     Renan Melo Barros
//               Separated from elmparam.cpp
// 
// Modified:
// -------------------------------------------------------------------------

#include <cmath>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <vector>

using namespace std;

#include "elmparamC1.h"
#include "anmodel.h"
#include "node.h"
#include "shape.h"
#include "shpbsp.h"
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

// ============================ ReadElmBSP =================================

void cElmParamC1 :: ReadElmBSP (string name, eShpType shp, eAnmType anm, bool tl)
{
  if (tl) name.append(" Total Lagrangian");

  // Read the number of elements elements.

  int n;
  if (!(in >> n) || (n < 1) || (n > cElement :: GetNumElm( )))
  {
    cout << "Error in the input of the number of " << name <<" elements!\n";
    exit(0);
  }

  // Create and read each element.

  int id;
  cElmParamC1 *elm;
  for (int i = 0; i < n; i++)
  {
    if (!(in >> id) || (id < 1) || (id > cElement :: GetNumElm( )))
    {
      cout << "Error in the input of " << name << " elements!\n";
      exit(0);
    }
    elm = (tl) ? new cElmParamC1(id,shp,anm) : new cElmParamC1(id,shp,anm);
    elm->Read( );

    // TODO: Descobrir uma maneira melhor de garantir isto para malhas com
    // diferentes quadraturas.
    if (elm->GetNumIntPnt( ) > MaxIntPnt) MaxIntPnt = elm->GetNumIntPnt( );

    // Register element into BSP patch\shape map.
    dynamic_cast<cShapeBsp*> (elm->GetShape( ))->PatElmMapAdd(id-1);
  }
  elm->UpdGlbVars( );  // Upd. global vars depending on the element type
}

// ============================== cElmParamC1 ==============================

cElmParamC1 :: cElmParamC1(int i, eShpType s, eAnmType a) : cElmParam(i,s,a)
{
}

// =============================== IntForce ================================

int cElmParamC1 :: IntForce(cVector &g)
{
  // Get the nodal coordinates.

  int nn = Shape->GetNumNode( );
  sNodeCoord *coord = new sNodeCoord[nn];
  Shape->NodalCoord(coord);

  // Get the nodal displacements.

  int ndof = nn*Model->GetNumDofNode( );
  cVector u(ndof);
  NodalDispl(u);

  // Create local matrices.

  int nsc = Model->GetDimBMatrix( );
  cMatrix B(nsc, ndof);
  cVector e(nsc);    // Strain
  cVector s(nsc);    // Stress

  // Get memory for shape functions and derivatives.

  double *shpfnc = new double[nn];
  double *mapfnc = new double[nn];
  sNodeCoord *shpdrv  = new sNodeCoord[nn];
  sNodeDrv *shpdrv2  = new sNodeDrv[nn];

  // Loop over integration points

  g.Zero( );
  int strok = 1;
  for (int i = 0; i < NumIntPnt; i++)
  {
    // Get integration point data.

    sNatCoord p = IntPnt[i].GetCoord( );
    double wgt  = IntPnt[i].GetWeight( );

    // Compute element strains: {e} = [B]{u}.

    double detJ;
    Shape->ShpFunc(p, shpfnc);
    Shape->MapFunc(p, mapfnc);
    Shape->DrvShpXYZ(p, coord, &detJ, shpdrv, shpdrv2);
    Model->BMatrix(Section->GetThickness( ), nn, shpfnc, mapfnc, coord, shpdrv, shpdrv2, B);
    e = B*u;

    // Compute the element stresses.

    if (!SecAn[i]->Stress(e, s))
    {
      strok = 0;
      break;
    }

    // Compute [g] += coeff*[B]t{s}.

    double volc  = Model->VolCoeff(Section->GetThickness( ), nn, mapfnc, coord);
    double coeff = wgt*detJ*volc;
    MultTAcc(coeff, B, s, g);
  }

  // Release memory.

  delete []coord;
  delete []shpfnc;
  delete []mapfnc;
  delete []shpdrv;
  delete []shpdrv2;

  return(strok);
}

// =============================== StiffMat ================================

void cElmParamC1 :: StiffMat(cMatrix &K)
{
  // Get the nodal coordinates.
  
  int nn = Shape->GetNumNode( );
  sNodeCoord *coord = new sNodeCoord[nn];
  Shape->NodalCoord(coord);

  // Create local arrays.

  int ndof = nn*Model->GetNumDofNode( );
  int nsc = Model->GetDimBMatrix( );
  cMatrix C(nsc, nsc);
  cMatrix B(nsc, ndof);

  // Get memory for shape functions and derivatives.

  double *shpfnc     = new double[nn];
  double *mapfnc     = new double[nn];
  sNodeCoord *shpdrv = new sNodeCoord[nn];
  sNodeDrv *shpdrv2  = new sNodeDrv[nn];

  // Loop over integration points

  K.Zero( );
  for (int i = 0; i < NumIntPnt; i++)
  {
    // Get integration point data.

    sNatCoord p = IntPnt[i].GetCoord( );
    double wgt  = IntPnt[i].GetWeight( );

    // Compute [B] matrix.

    double detJ;
    Shape->ShpFunc(p, shpfnc);
    Shape->MapFunc(p, mapfnc);
    Shape->DrvShpXYZ(p, coord, &detJ, shpdrv, shpdrv2);
    Model->GetMecMod( )->BMatrix(Section->GetThickness( ), nn, shpfnc, mapfnc, coord, shpdrv, shpdrv2, B);

    // Get [C] matrix.

    SecAn[i]->CMatrix(C);

    // Compute [K] += coeff*[B]t[C][B].

    double volc  = Model->VolCoeff(Section->GetThickness( ), nn, mapfnc, coord);
    double coeff = wgt*detJ*volc;
    MatTripBtCB(B, C, coeff, K);
  }
  
  // Release memory.

  delete []coord;
  delete []shpfnc;
  delete []mapfnc;
  delete []shpdrv;
  delete []shpdrv2;
}

// =============================== MassMat ================================

void cElmParamC1 :: MassMat(cMatrix &M)
{
  // Get the nodal coordinates.

  int nn = Shape->GetNumNode( );
  sNodeCoord *coord = new sNodeCoord[nn];
  Shape->NodalCoord(coord);

  // Create local arrays.

  int dimN  = Model->GetDimNMatrix( );
  int ndofn = Model->GetNumDofNode( );
  int ndof  = nn*ndofn;
  cMatrix Mb(dimN, dimN);
  cMatrix N(dimN, ndof);
//  cout << "dimN = " << dimN << "  ";
//  cout << "ndofn = " << ndofn << "  ";
//  cout << "ndof = " << ndof << "  ";

  // Get memory for shape functions and derivatives.

  int     nshp   = Shape->GetNumShp( );
  int     nmap   = Shape->GetNumMap( );
  double *shpfnc = new double[nshp];
  double *mapfnc = new double[nmap];
  sNodeCoord *shpdrv = new sNodeCoord[nmap];

  // Create the integration points.

  int NumMassIntPnt;
  eShpType type = Shape->GetType( );
  cIntPoint *MassIntPnt = cIntPoint::CreateIntPoints(Shape->GetTopologyType( ), VecIntOrd[OrdIdx].type,VecIntOrd[OrdIdx].M, &NumMassIntPnt);
//  cout << "NumMassIntPnt = " << NumMassIntPnt << endl;

  // Loop over integration points.

  M.Zero( );
  for (int i = 0; i < NumMassIntPnt; i++)
  {
    // Get integration point data.

    sNatCoord p = MassIntPnt[i].GetCoord( );
    double wgt  = MassIntPnt[i].GetWeight( );

    // Compute [N] matrix.

    double detJ;
    Shape->ShpFunc(p, shpfnc);
    Shape->MapFunc(p, mapfnc);
    Shape->DrvShpXYZ(p, coord, &detJ, shpdrv);
    Model->GetMecMod( )->NMatrix(nn, shpfnc, shpdrv, N);

    // Evaluate [Mb] matrix.

    SecAn[i]->MbMatrix(Mb);

    // Compute [M] += coeff*[N]t[Mb][N].

    double volc  = Model->VolCoeff(Section->GetThickness( ), nn, mapfnc, coord);
    double coeff = wgt*detJ*volc;
    MatTripBtCB(N, Mb, coeff, M);
  }

  // Release memory.

  delete []coord;
  delete []shpfnc;
  delete []mapfnc;
  delete []shpdrv;
  delete []MassIntPnt;
}

// ============================== GeomStiff ================================

void cElmParamC1 :: GeomStiff(cMatrix &G)
{
  // Get the nodal coordinates.

  int nn = Shape->GetNumNode( );
  sNodeCoord *coord = new sNodeCoord[nn];
  Shape->NodalCoord(coord);

  // Create local arrays.

  int ndof = nn*Model->GetNumDofNode( );
  int nsc = Model->GetDimBMatrix( );
  int bnldim = Model->GetMecMod( )->GetDimBnlMatrix();
  cMatrix C(nsc, nsc);
  cMatrix B(nsc, ndof);
  cMatrix Bnl(bnldim, ndof);
  cMatrix S(bnldim,bnldim);
  cVector e(nsc);
  cVector s(nsc);

  // Get the nodal displacements.

  cVector u(ndof);
  NodalDispl(u);

  // Get memory for shape functions and derivatives.

//  double *shpfnc = new double[nn];
//  double *mapfnc = new double[nn];
//  sNodeCoord *shpdrv  = new sNodeCoord[nn];
//  sNodeDrv *shpdrv2  = new sNodeDrv[nn];
  int     nshp   = Shape->GetNumShp( );
  int     nmap   = Shape->GetNumMap( );
  double *shpfnc = new double[nshp];
  double *mapfnc = new double[nmap];
  sNodeCoord *shpdrv = new sNodeCoord[nshp];
  sNodeDrv *shpdrv2  = new sNodeDrv[nshp];

  // Loop over integration points

  G.Zero( );
  for (int i = 0; i < NumIntPnt; i++)
  {
    // Get integration point data.

    sNatCoord p = IntPnt[i].GetCoord( );
    double wgt  = IntPnt[i].GetWeight( );

    // Compute [B] and [Bnl] matrices.

    double detJ;
    Shape->ShpFunc(p, shpfnc);
    Shape->MapFunc(p, mapfnc);
    Shape->DrvShpXYZ(p, coord, &detJ, shpdrv, shpdrv2);
    Model->BMatrix(Section->GetThickness( ), nn, shpfnc, mapfnc, coord, shpdrv, shpdrv2, B);
    Model->GetMecMod( )->BnlMatrix(nn, shpfnc, mapfnc, coord, shpdrv, Bnl);

    // Get [C] matrix.

    SecAn[i]->CMatrix(C);

    // Compute strains and stresses.

    e = B*u;
    if (!SecAn[i]->Stress(e, s)) break;

    // Evaluate the S matrix.

    Model->GetMecMod( )->SMatrix(s, S);

    // Compute [G] += [Bnl]t[S][Bnl].

    double volc  = Model->VolCoeff(Section->GetThickness( ), nn, mapfnc, coord);
    double coeff = wgt*detJ*volc;
    MatTripBtCB(Bnl, S, coeff, G);
  }

  // Release memory.

  delete []coord;
  delete []shpfnc;
  delete []mapfnc;
  delete []shpdrv;
  delete []shpdrv2;
}

// ============================= IntPntStress ==============================

void cElmParamC1 :: IntPntStress(cMatrix &Strain, cMatrix &Stress)
{
  // Get the nodal coordinates.

  int nn = Shape->GetNumNode( );
  sNodeCoord *coord = new sNodeCoord[nn];
  Shape->NodalCoord(coord);

  // Get the nodal displacements.

  int ndof = nn*Model->GetNumDofNode( );
  cVector u(ndof);
  NodalDispl(u);

  // Create local matrices.

  int nsc = Model->GetDimBMatrix( );
  cMatrix B(nsc, ndof);
  cVector e(nsc);    // Strain
  cVector s(nsc);    // Stress

  // Get memory for shape functions and derivatives.

  double *shpfnc = new double[nn];
  double *mapfnc = new double[nn];
  sNodeCoord *shpdrv  = new sNodeCoord[nn];
  sNodeDrv *shpdrv2  = new sNodeDrv[nn];

  // Loop over integration points

  Strain.Zero( );
  Stress.Zero( );
  for (int i = 0; i < NumIntPnt; i++)
  {
    // Get integration point data.

    sNatCoord p = IntPnt[i].GetCoord( );
 
    // Compute element strains: {e} = [B]{u}.

    double detJ;
    Shape->ShpFunc(p, shpfnc);
    Shape->MapFunc(p, mapfnc);
    Shape->DrvShpXYZ(p, coord, &detJ, shpdrv, shpdrv2);
    Model->BMatrix(Section->GetThickness( ), nn, shpfnc, mapfnc, coord, shpdrv, shpdrv2, B);
    e = B*u;

    // Compute the element stresses.

    SecAn[i]->Stress(e, s);

    // Return the computed values and store the current stresses.

    for (int j = 0; j < nsc; j++)
    {
      Strain[i][j] = e[j];
      Stress[i][j] = StrIpt[i][j] = s[j];
    }
  }

  // Release memory.

  delete []coord;
  delete []shpfnc;
  delete []mapfnc;
  delete []shpdrv;
  delete []shpdrv2;
}

// ============================= PntStress =================================

void cElmParamC1 :: PntStress(sNatCoord *Pnts, int NumPnt, cMatrix &Strain, 
                            cMatrix &Stress, bool opt)
{
  // Get the nodal coordinates.

  int nn = Shape->GetNumNode( );
  sNodeCoord *coord = new sNodeCoord[nn];
  Shape->NodalCoord(coord);

  // Get the nodal displacements.

  int ndof = nn*Model->GetNumDofNode( );
  cVector u(ndof);
  NodalDispl(u);

  // Create local matrices.

  int nsc = Model->GetDimBMatrix( );
  cMatrix B(nsc, ndof);
  cVector e(nsc);    // Strain
  cVector s(nsc);    // Stress

  // Get memory for shape functions and derivatives.

  double *shpfnc = new double[nn];
  double *mapfnc = new double[nn];
  sNodeCoord *shpdrv = new sNodeCoord[nn];
  sNodeDrv *shpdrv2  = new sNodeDrv[nn];

  // Evaluate the closest integration point to each input point.
  static eShpType   lShp; 
  static int        lOrd = -1;
  static vector<int> VecID;

  if (!opt || (lShp != Shape->GetType( )) || (lOrd != OrdIdx))
  {
    // Compute new index vector.
    NearestIntPnt(Shape->NumDim( ),Pnts,NumPnt,IntPnt,NumIntPnt,VecID);

    // Update last input data, shape type and integration orden index.
    lShp = Shape->GetType( );
    lOrd = OrdIdx;
  }

  // Loop over input points

  Strain.Zero( );
  Stress.Zero( );
  sNatCoord p;
  for (int i = 0; i < NumPnt; i++)
  {
    // Get point data.

    p = Pnts[i];

    // Compute element strains: {e} = [B]{u}.

    double detJ;
    Shape->ShpFunc(p, shpfnc);
    Shape->MapFunc(p, mapfnc);
    Shape->DrvShpXYZ(p, coord, &detJ, shpdrv, shpdrv2);
    Model->BMatrix(Section->GetThickness( ), nn, shpfnc, mapfnc, coord, shpdrv, shpdrv2, B);
    e = B*u;

    // Compute the element stresses.

    SecAn[VecID[i]]->Stress(e, s);

    // Return the computed values and store the current stresses.

    for (int j = 0; j < nsc; j++)
    {
      Strain[i][j] = e[j];
      Stress[i][j] = s[j];
    }
  }

  // Release memory.

  delete []coord;
  delete []shpfnc;
  delete []mapfnc;
  delete []shpdrv;
  delete []shpdrv2;
}


// -------------------------------------------------------------------------
// Methods of cElmParamC1TL class:
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Public methods:
//

// ============================ ReadElmBSP =================================

void cElmParamC1TL :: ReadElmBSP (string name, eShpType shp, eAnmType anm, bool tl)
{
  if (tl) name.append(" Total Lagrangian");

  // Read the number of elements elements.

  int n;
  if (!(in >> n) || (n < 1) || (n > cElement :: GetNumElm( )))
  {
    cout << "Error in the input of the number of " << name <<" elements!\n";
    exit(0);
  }

  // Create and read each element.

  int id;
  cElmParamC1TL *elm;
  for (int i = 0; i < n; i++)
  {
    if (!(in >> id) || (id < 1) || (id > cElement :: GetNumElm( )))
    {
      cout << "Error in the input of " << name << " elements!\n";
      exit(0);
    }
    elm = (tl) ? new cElmParamC1TL(id,shp,anm) : new cElmParamC1TL(id,shp,anm);
    elm->Read( );

    // TODO: Descobrir uma maneira melhor de garantir isto para malhas com
    // diferentes quadraturas.
    if (elm->GetNumIntPnt( ) > MaxIntPnt) MaxIntPnt = elm->GetNumIntPnt( );

    // Register element into BSP patch\shape map.
    dynamic_cast<cShapeBsp*> (elm->GetShape( ))->PatElmMapAdd(id-1);
  }
  elm->UpdGlbVars( );  // Upd. global vars depending on the element type
}

// ============================== cElmParamC1 ==============================

cElmParamC1TL :: cElmParamC1TL(int i, eShpType s, eAnmType a) : cElmParamC1(i,s,a)
{
  Type = PARAMETRICTL;
}

// =============================== IntForce ================================

int cElmParamC1TL :: IntForce(cVector &g)
{
  // Get the nodal coordinates.

  int nn = Shape->GetNumNode( );
  sNodeCoord *coord = new sNodeCoord[nn];
  Shape->NodalCoord(coord);

  // Get the nodal displacements.

  int ndof = nn*Model->GetNumDofNode( );
  cVector u(ndof);
  NodalDispl(u);

  // Create local matrices.

  int nsc = Model->GetDimBMatrix( );
  cMatrix B(nsc, ndof);
  cMatrix Bl(nsc, ndof);
  cVector e(nsc);    // Strain
  cVector s(nsc);    // Stress

  // Get memory for shape functions and derivatives.

  int     nshp   = Shape->GetNumShp( );
  int     nmap   = Shape->GetNumMap( );
  double *shpfnc = new double[nshp];
  double *mapfnc = new double[nmap];
  sNodeCoord *shpdrv = new sNodeCoord[nn];
  sNodeDrv *shpdrv2  = new sNodeDrv[nshp];

  // Loop over integration points

  g.Zero( );
  int strok = 1;
  for (int i = 0; i < NumIntPnt; i++)
  {
    // Get integration point data.

    sNatCoord p = IntPnt[i].GetCoord( );
    double wgt  = IntPnt[i].GetWeight( );

    // Compute element strains.

    double detJ;
    Shape->ShpFunc(p, shpfnc);
    Shape->MapFunc(p, mapfnc);
    Shape->DrvShpXYZ(p, coord, &detJ, shpdrv, shpdrv2);
    Model->BMatrix(Section->GetThickness( ), nn, shpfnc, mapfnc, coord, shpdrv, shpdrv2, B);
    Model->GetMecMod( )->BlMatrix(nn, shpfnc, mapfnc, u, coord, shpdrv, Bl);
    GreenStrain(B, Bl, u, e);

    // Compute the element stresses.

    if (!SecAn[i]->Stress(e, s))
    {
      strok = 0;
      break;
    }

    // Calculate [B] = [B] + [Bl]

    B += Bl;

    // Compute [g] += coeff*[B]t{s}.

    double volc  = Model->VolCoeff(Section->GetThickness( ), nn, mapfnc, coord);
    double coeff = wgt*detJ*volc;
    MultTAcc(coeff, B, s, g);
  }


  // Release memory.

  delete []coord;
  delete []shpfnc;
  delete []mapfnc;
  // delete []shpdrv; // Tirar?
  delete []shpdrv2;



  return(strok);
}

// =============================== StiffMat ================================

void cElmParamC1TL :: StiffMat(cMatrix &Kt)
{

  // Get the nodal coordinates.

  int nn = Shape->GetNumNode( );
  sNodeCoord *coord = new sNodeCoord[nn];
  Shape->NodalCoord(coord);

  // Create local arrays.

  int ndof = nn*Model->GetNumDofNode( );
  int nsc = Model->GetDimBMatrix( );
  int bnldim = Model->GetMecMod( )->GetDimBnlMatrix();
  cMatrix C(nsc, nsc);
  cMatrix B(nsc, ndof);
  cMatrix Bl(nsc, ndof);
  cMatrix Bnl(bnldim, ndof);
  cMatrix S(bnldim,bnldim);
  cVector e(nsc);
  cVector s(nsc);

  // Get the nodal displacements.

  cVector u(ndof);
  NodalDispl(u);

  // Get memory for shape functions and derivatives.

  int     nshp   = Shape->GetNumShp( );
  int     nmap   = Shape->GetNumMap( );
  double *shpfnc = new double[nshp];
  double *mapfnc = new double[nmap];
  sNodeCoord *shpdrv = new sNodeCoord[nshp];
  sNodeDrv *shpdrv2  = new sNodeDrv[nshp];
  
  // Loop over integration points

  Kt.Zero( );
  for (int i = 0; i < NumIntPnt; i++)
  {
    // Get integration point data.

    sNatCoord p = IntPnt[i].GetCoord( );
    double wgt  = IntPnt[i].GetWeight( );

    // Compute [B], [Bl] and [Bnl] matrices.

    double detJ;
    Shape->ShpFunc(p, shpfnc);
    Shape->MapFunc(p, mapfnc);
    Shape->DrvShpXYZ(p, coord, &detJ, shpdrv, shpdrv2);
    Model->BMatrix(Section->GetThickness( ), nn, shpfnc, mapfnc, coord, shpdrv, shpdrv2, B);
    Model->GetMecMod( )->BlMatrix(nn, shpfnc, mapfnc, u, coord, shpdrv, Bl);
    Model->GetMecMod( )->BnlMatrix(nn, shpfnc, mapfnc, coord, shpdrv, Bnl);

    // Get [C] matrix.

    SecAn[i]->CMatrix(C);

    // Evaluate the Piola-Kirchhoff II stress vector.

    GreenStrain(B, Bl, u, e);

    cVector temp(2);

    int therm = GetStrnTemp(p, shpfnc, temp);

    if (therm)
    {
     double tref = 0;
     if (!SecAn[i]->Stress(tref, temp, e, s)) break;
    }
    else
    {
     if (!SecAn[i]->Stress(e, s)) break;
    }

    // Evaluate the S matrix.

    Model->GetMecMod( )->SMatrix(s, S);

    // Add [Bl] to [B].

    B += Bl;

    // Compute [K] += coeff*([B]t[C][B] + [Bnl]t[S][Bnl]).

    double volc  = Model->VolCoeff(Section->GetThickness( ), nn, mapfnc, coord);
    double coeff = wgt*detJ*volc;
    MatTripBtCB(B, C, coeff, Kt);
    MatTripBtCB(Bnl, S, coeff, Kt);

//    static double vol = 0;
//    vol += coeff;
//    cout << "volume = " << vol << endl;
  }

  // Release memory.

  delete []coord;
  delete []shpfnc;
  delete []mapfnc;
  delete []shpdrv;
  delete []shpdrv2;
}

// ============================= IntPntStress ==============================

void cElmParamC1TL :: IntPntStress(cMatrix &Strain, cMatrix &Stress)
{
  // Get the nodal coordinates.

  int nn = Shape->GetNumNode( );
  sNodeCoord *coord = new sNodeCoord[nn];
  Shape->NodalCoord(coord);

  // Get the nodal displacements.

  int ndof = nn*Model->GetNumDofNode( );
  cVector u(ndof);
  NodalDispl(u);

  // Create local matrices.

  int nsc = Model->GetDimBMatrix( );
  cMatrix B(nsc, ndof);
  cMatrix Bl(nsc,ndof);
  cVector e(nsc);    // Strain
  cVector s(nsc);    // Stress

  // Get memory for shape functions and derivatives.

  int     nshp   = Shape->GetNumShp( );
  int     nmap   = Shape->GetNumMap( );
  double *shpfnc = new double[nshp];
  double *mapfnc = new double[nmap];
  sNodeCoord *shpdrv = new sNodeCoord[nshp];
  sNodeDrv *shpdrv2  = new sNodeDrv[nshp];

  // Loop over integration points

  Strain.Zero( );
  Stress.Zero( );
  for (int i = 0; i < NumIntPnt; i++)
  {
    // Get integration point data.

    sNatCoord p = IntPnt[i].GetCoord( );

    // Compute element strains.

    double detJ;
    Shape->ShpFunc(p, shpfnc);
    Shape->MapFunc(p, mapfnc);
    Shape->DrvShpXYZ(p, coord, &detJ, shpdrv, shpdrv2);
    Model->BMatrix(Section->GetThickness( ), nn, shpfnc, mapfnc, coord, shpdrv, shpdrv2, B);
    Model->GetMecMod( )->BlMatrix(nn, shpfnc, mapfnc, u, coord, shpdrv, Bl);

    GreenStrain(B, Bl, u, e);

    // Compute the element stresses.

    cVector temp(2);

    int therm = GetStrnTemp(p, shpfnc, temp);

    if (therm)
    {
     double tref = 0;
     SecAn[i]->Stress(tref, temp, e, s);
    }
    else
    {
     SecAn[i]->Stress(e, s);
    }

    // Return the computed values and store the current stresses.

    for (int j = 0; j < nsc; j++)
    {
      Strain[i][j] = e[j];
      Stress[i][j] = StrIpt[i][j] = s[j];
    }
  }

  // Release memory.

  delete []coord;
  delete []shpfnc;
  delete []mapfnc;
  delete []shpdrv;
  delete []shpdrv2;
}

// ============================= PntStress =================================

void cElmParamC1TL :: PntStress(sNatCoord *Pnts, int NumPnt, cMatrix &Strain, 
                            cMatrix &Stress, bool opt)
{
  // Get the nodal coordinates.

  int nn = Shape->GetNumNode( );
  sNodeCoord *coord = new sNodeCoord[nn];
  Shape->NodalCoord(coord);

  // Get the nodal displacements.

  int ndof = nn*Model->GetNumDofNode( );
  cVector u(ndof);
  NodalDispl(u);

  // Create local matrices.

  int nsc = Model->GetDimBMatrix( );
  cMatrix B(nsc, ndof);
  cVector e(nsc);    // Strain
  cVector s(nsc);    // Stress

  // Get memory for shape functions and derivatives.

  int     nshp   = Shape->GetNumShp( );
  int     nmap   = Shape->GetNumMap( );
  double *shpfnc = new double[nshp];
  double *mapfnc = new double[nmap];
  sNodeCoord *shpdrv = new sNodeCoord[nshp];
  sNodeDrv *shpdrv2  = new sNodeDrv[nshp];

  // Evaluate the closest integration point to each input point.
  static eShpType   lShp; 
  static int        lOrd = -1;
  static vector<int> VecID;

  if (!opt || (lShp != Shape->GetType( )) || (lOrd != OrdIdx))
  {
    // Compute new index vector.
    NearestIntPnt(Shape->NumDim( ),Pnts,NumPnt,IntPnt,NumIntPnt,VecID);

    // Update last input data, shape type and integration orden index.
    lShp = Shape->GetType( );
    lOrd = OrdIdx;
  }

  // Loop over input points

  Strain.Zero( );
  Stress.Zero( );
  sNatCoord p;
  for (int i = 0; i < NumPnt; i++)
  {
    // Get point data.

    p = Pnts[i];

    // Compute element strains: {e} = [B]{u}.

    double detJ;
    Shape->ShpFunc(p, shpfnc);
    Shape->MapFunc(p, mapfnc);
    Shape->DrvShpXYZ(p, coord, &detJ, shpdrv, shpdrv2);
    Model->BMatrix(Section->GetThickness( ), nn, shpfnc, mapfnc, coord, shpdrv, shpdrv2, B);
    e = B*u;

    // Compute the element stresses.

    SecAn[VecID[i]]->Stress(e, s);

    // Return the computed values and store the current stresses.

    for (int j = 0; j < nsc; j++)
    {
      Strain[i][j] = e[j];
      Stress[i][j] = s[j];
    }
  }

  // Release memory.

  delete []coord;
  delete []shpfnc;
  delete []mapfnc;
  delete []shpdrv;
  delete []shpdrv2;
}

// -------------------------------------------------------------------------
// Protected methods:
//

// ============================ GreenStrain ================================

void cElmParamC1TL :: GreenStrain(cMatrix &B, cMatrix &Bl, cVector &u, cVector &str)
{
  // Get number of stress components from Analysis Model Class.

  int nsc = Model->GetDimBMatrix( );

  // Evaluate linear strain vector => {s1} = [B]*{d}.

  cVector str1(nsc);
  str1 = B*u;

  // Evaluate nonlinear strain vector => {s2} = 0.5*[Bl]{d}.

  cVector str2(nsc);
  str2 = Bl*u;
  str2 *= 0.5;

  // Add linear and nonlinear terms.

  str = str1 + str2;
} 

// ======================================================= End of file =====
