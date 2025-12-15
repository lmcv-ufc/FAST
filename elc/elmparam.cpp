// -------------------------------------------------------------------------
// elmparam.cpp - implementation of parametric elements.
// -------------------------------------------------------------------------
// Created:      30-Apr-2005     Evandro Parente Junior
//
// Modified:     14-May-2011     Evandro Parente Junior
//               Generalization to handle non-isoparametric elements.
//
// Modified:     15-Jul-2011     Iuri Barcelos Rocha
//               Removed element reading functions, moved to cSection class.
//
// Modified:     22-Jun-2012     Iuri Barcelos Rocha
//               Implementation of Total Lagrangian non-linear elements.
//
// Modified:     20-Sep-2012     Iuri Barcelos Rocha
//               Implementation of geometric stiffness.
//
// Modified:     08-Mar-2012     Iuri Barcelos Rocha
//               Use of the new cSection and cSecAnalysis classes.
//
// Modified:     08-Apr-2015     Elias Saraiva Barroso
//               Input/output operations using c++ streams.
//
// Modified:     19-Feb-2016     Elias Saraiva Barroso
//               Change the function GetIntOrd to return the maximum integration
//               order considered by the element.
//
// Modified:     31-Jan-2018     Elias Saraiva Barroso
//               Implementation of NodalStress( ) without extrapolation, and
//               Creation of PntStress( ) method.
//
// Modified:     04-May-2022     Renan Melo Barros and Elias Saraiva Barroso
//               Implementation of C1 parametric elements.
// -------------------------------------------------------------------------

#include <cmath>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <vector>

using namespace std;

#include "elmparam.h"
#include "elmdsh.h"
#include "ctrl.h"
#include "anmodel.h"
#include "node.h"
#include "shape.h"
#include "shpbsp.h"
#include "material.h"
#include "vec.h"
#include "mat.h"
#include "gblvar.h"
#include "gbldef.h"
#include "section.h"
#include "secanalysis.h"

// -------------------------------------------------------------------------
// Static variables:
//
cMatrix  cElmParam :: TR;
int      cElmParam :: PrvOrdIdx = -1;
eShpType cElmParam :: PrvShpType;

// -------------------------------------------------------------------------
// Template Methods:
//

// ============================ ReadElmParam ====================================

template <typename T>
void cElmParam :: ReadElm(string name, eShpType shp, eAnmType anm)
{
  // Read the number of elements elements.

  int n;
  if (!(in >> n) || (n < 1) || (n > cElement :: GetNumElm( )))
  {
    cout << "Error in the input of the number of " << name <<" elements!\n";
    exit(0);
  }

  // Create and read each element.

  int id;
  cElmParam *elm;
  for (int i = 0; i < n; i++)
  {
    if (!(in >> id) || (id < 1) || (id > cElement :: GetNumElm( )))
    {
      cout << "Error in the input of " << name << " elements!\n";
      exit(0);
    }
    elm = new T(id,shp,anm);
    elm->Read( );

    // TODO: Descobrir uma maneira melhor de garantir isto para malhas com
    // diferentes quadraturas.
    if (elm->GetNumIntPnt( ) > MaxIntPnt) MaxIntPnt = elm->GetNumIntPnt( );

    // Register element into BSP patch\shape map.
    cShapeBsp *shpbsp = dynamic_cast<cShapeBsp*> (elm->GetShape( ));
    if (shpbsp) shpbsp->PatElmMapAdd(id-1);
  }
  elm->UpdGlbVars( );  // Upd. global vars depending on the element type
}
template void cElmParam :: ReadElm<cElmParam>(string,eShpType,eAnmType);
template void cElmParam :: ReadElm<cElmParamTL>(string,eShpType,eAnmType);
template void cElmParam :: ReadElm<cElmDegShell>(string,eShpType,eAnmType);
template void cElmParam :: ReadElm<cElmDegShellTL>(string,eShpType,eAnmType);

// -------------------------------------------------------------------------
// Public methods:
//

// ============================ ReadElmFEM =================================

void cElmParam :: ReadElmFEM (string name, eShpType shp, eAnmType anm, bool tl)
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
  cElmParam *elm;
  for (int i = 0; i < n; i++)
  {
    if (!(in >> id) || (id < 1) || (id > cElement :: GetNumElm( )))
    {
      cout << "Error in the input of " << name << " elements!\n";
      exit(0);
    }
    elm = (tl) ? new cElmParamTL(id,shp,anm) : new cElmParam(id,shp,anm);
    elm->Read( );
  }
  elm->UpdGlbVars( );  // Upd. global vars depending on the element type
}

// ============================ ReadElmBSP =================================

void cElmParam :: ReadElmBSP (string name, eShpType shp, eAnmType anm, bool tl)
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
  cElmParam *elm;
  for (int i = 0; i < n; i++)
  {
    if (!(in >> id) || (id < 1) || (id > cElement :: GetNumElm( )))
    {
      cout << "Error in the input of " << name << " elements!\n";
      exit(0);
    }
    elm = (tl) ? new cElmParamTL(id,shp,anm) : new cElmParam(id,shp,anm);
    elm->Read( );

    // TODO: Descobrir uma maneira melhor de garantir isto para malhas com
    // diferentes quadraturas.
    if (elm->GetNumIntPnt( ) > MaxIntPnt) MaxIntPnt = elm->GetNumIntPnt( );

    // Register element into BSP patch\shape map.
    cShapeBsp *shpbsp = dynamic_cast<cShapeBsp*> (elm->GetShape( ));
    if (shpbsp) shpbsp->PatElmMapAdd(id-1);
  }
  elm->UpdGlbVars( );  // Upd. global vars depending on the element type
}

// =============================== cElmParam ===============================

cElmParam :: cElmParam(int id, eShpType tshp, eAnmType tanm):cElement(id)
{
  Type  = PARAMETRIC;
  Shape = cShape::CreateShape(tshp);
  Model = cAnModel::CreateModel(tanm);
}

// ============================= ~cElmParam ================================

cElmParam :: ~cElmParam(void)
{
  delete Model;
  delete []IntPnt;
  if (SecAn)
  {
    for (int i = 0; i < NumIntPnt; i++) delete SecAn[i];
    delete []SecAn;
  }
}

// ============================== UpdateState ==============================

void cElmParam :: UpdateState(void)
{
  for (int i = 0; i < NumIntPnt; i++) SecAn[i]->UpdateState( );
}

// ================================= Read ==================================

void cElmParam :: Read(void)
{
  int sec,ord;

  if (!(in >> sec) || !(in >> ord))
  {
    cout << "Error in the input of element " << Label << "!\n";
    exit(0);
  }

  Section = cSection::GetSection(sec);
  if (!Section)
  {
    cout << "Invalid element section for element " << Label << "!\n";
    exit(0);
  }

  if (ord < 1 || ord > NumIntOrd)
  {
    cout << "Invalid integration order for element " << Label << "!\n";
    exit(0);
  }
  OrdIdx = ord - 1;

  // Read element incidence.

  Shape->Read( );

  // Create SecAnalysis objects.

  eSecType sectype = Section->GetType( );
  eTopType shptype = Shape->GetTopologyType( );
  if (sectype == SEC_LAMINATED_SOLID)
  {
    double ttop, tbot;
    int nlam = Section->GetNumLam( );
    int npts = nlam*VecIntOrd[OrdIdx].K[0]*VecIntOrd[OrdIdx].K[1]*VecIntOrd[OrdIdx].K[2];
    sLamina *Lamina = Section->GetLayup( );
    IntPnt = new cIntPoint[npts];
    SecAn = new cSecAnalysis*[npts];

    int k = 0;
    cIntPoint *p;

    tbot = -1.0;
    for (int i = 0; i < nlam; i++)
    {
      ttop = tbot + 2*Lamina[i].Thk;
      double a = (ttop + tbot)/2.0;
      double b = (ttop - tbot)/2.0;

      p = cIntPoint::CreateIntPoints(shptype, VecIntOrd[OrdIdx].type,VecIntOrd[OrdIdx].K, &NumIntPnt);

      for (int j = 0; j < NumIntPnt; j++)
      {
        IntPnt[k] = p[j];
        SecAn[k] = new cSecAnLaminatedSolid(this, &IntPnt[k], &Lamina[i]);

        // Modify the integration point weight and t coordinate.

        double oldwgt = IntPnt[k].GetWeight( );
        double oldt   = IntPnt[k].GetCoord( ).t;
        IntPnt[k].SetWeight(oldwgt*b);
        IntPnt[k].SetTCoord(a + b*oldt);
        k++;
      }
      tbot = ttop;
    }

    // Set the total number of integration points.

    NumIntPnt = npts;
    StrIpt.Resize(NumIntPnt, Model->GetDimBMatrix( ));
  }
  else if (Model->GetType( ) == PLANE_HEAT_TRANSFER)
  {
    // Create the integration points.
    //
    IntPnt = cIntPoint::CreateIntPoints(shptype, VecIntOrd[OrdIdx].type,VecIntOrd[OrdIdx].K, &NumIntPnt);
    StrIpt.Resize(NumIntPnt, Model->GetDimBMatrix( ));

    SecAn = 0;
  }
  else
  {
    // Create the integration points.

    IntPnt = cIntPoint::CreateIntPoints(shptype, VecIntOrd[OrdIdx].type,VecIntOrd[OrdIdx].K, &NumIntPnt);
    StrIpt.Resize(NumIntPnt, Model->GetDimBMatrix( ));

    // Create the constitutive model for each integration point.

    SecAn = new cSecAnalysis*[NumIntPnt];
    for (int i = 0; i < NumIntPnt; i++)
      SecAn[i] = cSecAnalysis::Create(this, &IntPnt[i]);
  }
}

// ============================== ReadIniStr ===============================

void cElmParam :: ReadIniStr(void)
{
  int nsc = Model->GetDimBMatrix( );
  IniStr.Resize(nsc);

  double val;
  int nread = 0;
  for (int i = 0; i < nsc; i++)
  {
    if (in >> val)
      nread++;
    else
    {
      cout << "Error in the input of initial stresses of element " << Label << "!\n";
      exit(0);
    }
    IniStr[i] = val;
  }

  if (nread != nsc)
  {
    cout << "Error in the input of initial stresses of element " << Label << "!\n";
    exit(0);
  }
}

// =============================== GetIntOrd ===============================

int cElmParam :: GetIntOrd(void)
{
  int numdim = Shape->NumDim( ), Max = VecIntOrd[OrdIdx].K[0];

  if (numdim > 1)
    Max = max(VecIntOrd[OrdIdx].K[1],Max);

  if (numdim > 2)
    Max = max(VecIntOrd[OrdIdx].K[2],Max);

  return Max;
}

// =============================== GetIntPnt ===============================

void cElmParam :: GetIntPnt(int id, sNatCoord &pnt, double &w)
{
  pnt.r = IntPnt[id].GetCoord( ).r;
  pnt.s = IntPnt[id].GetCoord( ).s;
  pnt.t = IntPnt[id].GetCoord( ).t;
  w     = IntPnt[id].GetWeight( );
}

// =============================== IntForce ================================
/*
double cElmParam :: EnergyIntForce(cVector &g)
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
    Shape->DrvShpXYZ(p, coord, &detJ, shpdrv);
    Model->BMatrix(nn, shpfnc, mapfnc, coord, shpdrv, B);
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

  return(strok);
}
*/

// =============================== IntForce ================================

int cElmParam :: IntForce(cVector &g)
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
    Shape->DrvShpXYZ(p, coord, &detJ, shpdrv);
    Model->BMatrix(nn, shpfnc, mapfnc, coord, shpdrv, B);
    e = B*u;

    // Compute the element stresses.

    cVector temp(2);
    int therm = GetStrnTemp(p, shpfnc, temp);
    if (therm)
    {
      double tref = 0;
      if(!SecAn[i]->Stress(tref, temp, e, s))
      {
        strok = 0;
        break;
      }
    }
    else
    {
      if (!SecAn[i]->Stress(e, s))
      {
       strok = 0;
       break;
      }
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

  return(strok);
}

// =============================== StiffMat ================================

void cElmParam :: StiffMat(cMatrix &K)
{
  // Get the nodal coordinates.
  
  int nn = Shape->GetNumNode( );
  sNodeCoord *coord = new sNodeCoord[nn];
  Shape->NodalCoord(coord);

  // Create local arrays.

  int ndof = nn*Model->GetNumDofNode( );
  int nsc  = Model->GetDimBMatrix( );
  cMatrix C(nsc, nsc);
  cMatrix B(nsc, ndof);

  // Get memory for shape functions and derivatives.

  int     nshp   = Shape->GetNumShp( );
  int     nmap   = Shape->GetNumMap( );
  double *shpfnc = new double[nshp];
  double *mapfnc = new double[nmap];
  sNodeCoord *shpdrv = new sNodeCoord[nshp];

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
    Shape->DrvShpXYZ(p, coord, &detJ, shpdrv);
    Model->BMatrix(nn, shpfnc, mapfnc, coord, shpdrv, B);

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
}

// =============================== ConducMat ===============================

void cElmParam :: ConducMat(cMatrix &K)
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

  double *shpfnc = new double[nn];
  double *mapfnc = new double[nn];
  sNodeCoord *shpdrv = new sNodeCoord[nn];

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
    Shape->DrvShpXYZ(p, coord, &detJ, shpdrv);
    Model->BMatrix(nn, shpfnc, mapfnc, coord, shpdrv, B);

    // Get Conductivity matrix.
    
    double param[3];
    Section->GetMaterial( )->GetThermProp( )->GetParam(param);
    if (Section->GetMaterial( )->GetThermProp( )->GetType( ) == MAT_THERMAL_ISOTROPIC)
      Model->GetHeatMod( )->CondMatrix(param,C);
    else
      Model->GetHeatMod( )->CondMatrixOrtho(param,C);

    // Compute [K] += coeff*[B]t[C][B].

    double volc  = Model->VolCoeff(Section->GetThickness( ), nn, mapfnc, coord);
    double coeff = wgt*detJ*volc;
    MatTripBtCB(B, C, coeff, K);

    static double volume = 0;
    #pragma omp critical
    volume += coeff;
    cout << setprecision(16);
    cout << "volume " << volume << endl;
  }

  // Release memory.

  delete []coord;
  delete []shpfnc;
  delete []mapfnc;
  delete []shpdrv;
}

// =============================== MassMat ================================

void cElmParam :: MassMat(cMatrix &M)
{
  // Get the nodal coordinates.

  int nn = Shape->GetNumNode( );
  sNodeCoord *coord = new sNodeCoord[nn];
  Shape->NodalCoord(coord);

  // Create local arrays.

  int ndofn = Model->GetNumDofNode( );
  int ndof = nn*ndofn;
  cMatrix Mb(ndofn, ndofn);
  cMatrix N(ndofn, ndof);

  // Get memory for shape functions and derivatives.

  int     nshp   = Shape->GetNumShp( );
  int     nmap   = Shape->GetNumMap( );
  double *shpfnc = new double[nshp];
  double *mapfnc = new double[nmap];
  sNodeCoord *shpdrv = new sNodeCoord[nmap];

  // Create the integration points.

  //int NumMassIntPnt;
  //eShpType type = Shape->GetType( );
  //cIntPoint *MassIntPnt = cIntPoint::CreateIntPoints(Shape->GetTopologyType( ), VecIntOrd[OrdIdx].type,VecIntOrd[OrdIdx].M, &NumMassIntPnt);
  //cout << "NumMassIntPnt  = " << NumMassIntPnt << endl;

  // Loop over integration points.

  M.Zero( );
  for (int i = 0; i < NumIntPnt; i++)
  {
    // Get integration point data.

    sNatCoord p = IntPnt[i].GetCoord( );
    double wgt  = IntPnt[i].GetWeight( );

    // Compute [N] matrix.

    double detJ;
    Shape->ShpFunc(p, shpfnc);
    Shape->MapFunc(p, mapfnc);
    Shape->DrvShpXYZ(p, coord, &detJ, shpdrv);
    //Model->GetMecMod( )->NMatrix(nn, shpfnc, N);
    Shape->NMatrix(ndofn,shpfnc,N);

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
  //delete []MassIntPnt;
}

// ============================== GeomStiff ================================

void cElmParam :: GeomStiff(cMatrix &G)
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

  int     nshp   = Shape->GetNumShp( );
  int     nmap   = Shape->GetNumMap( );
  double *shpfnc = new double[nshp];
  double *mapfnc = new double[nmap];
  sNodeCoord *shpdrv = new sNodeCoord[nshp];

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
    Shape->DrvShpXYZ(p, coord, &detJ, shpdrv);
    Model->BMatrix(nn, shpfnc, mapfnc, coord, shpdrv, B);
    Model->GetMecMod( )->BnlMatrix(nn, shpfnc, mapfnc, coord, shpdrv, Bnl);

    // Get [C] matrix.

    SecAn[i]->CMatrix(C);

    // Compute strains and stresses.

    e = B*u;

    cVector temp(2);
    int therm = GetStrnTemp(p, shpfnc, temp);
    if (therm)
    {
     double tref = 0;
     if (!SecAn[i]->Stress(tref, temp, e, s)) break;
    }
    else
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
}

// ============================= IntPntTemp ================================

void cElmParam :: IntPntTemp(cVector &Temp)
{
  // Get memory for shape functions.
  int nn = Shape->GetNumNode( );
  double *shpfnc = new double[nn];

  // Loop over integration points
  Temp.Zero( );
  for (int i = 0; i < NumIntPnt; i++)
  {
    // Get integration point data.
    sNatCoord p = IntPnt[i].GetCoord( );

    // Compute temperature at current point.
    Shape->ShpFunc(p, shpfnc);
    for(int n = 0; n < nn; ++n)
      Temp[i] += shpfnc[n] * Shape->GetNode(n)->GetTemp( );
  }

  // Release memory.
  delete []shpfnc;
}

// ============================= IntPntStress ==============================

void cElmParam :: IntPntStress(cMatrix &Strain, cMatrix &Stress)
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

  // Loop over integration points.

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
    Shape->DrvShpXYZ(p, coord, &detJ, shpdrv);
    Model->BMatrix(nn, shpfnc, mapfnc, coord, shpdrv, B);
    e = B*u;

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
}

// ============================= IntPntFlux ================================

void cElmParam :: IntPntFlux(cMatrix &Flux)
{
  // Get the nodal coordinates.

  int nn = Shape->GetNumNode( );
  sNodeCoord *coord = new sNodeCoord[nn];
  Shape->NodalCoord(coord);

  // Get the nodal temperatures.

  cVector t(nn);
  NodalTemp(t);

  // Create local matrices.

  int nsc = Model->GetDimBMatrix( );
  cMatrix B(nsc, nn);
  cMatrix C(nsc, nsc);
  cVector dt(nsc);    // Derivative
  cVector fl(nsc);    // Fluxes

  // Get memory for shape functions and derivatives.

  double *shpfnc = new double[nn];
  double *mapfnc = new double[nn];
  sNodeCoord *shpdrv = new sNodeCoord[nn];

  // Loop over integration points

  Flux.Zero( );
  for (int i = 0; i < NumIntPnt; i++)
  {
    // Get integration point data.

    sNatCoord p = IntPnt[i].GetCoord( );

    // Compute element derivatives: {dt} = [B]{t}.

    double detJ;
    Shape->ShpFunc(p, shpfnc);
    Shape->MapFunc(p, mapfnc);
    Shape->DrvShpXYZ(p, coord, &detJ, shpdrv);
    Model->BMatrix(nn, shpfnc, mapfnc, coord, shpdrv, B);
    dt = B*t;

    // Get Conductivity matrix.

    double param[3];
    Section->GetMaterial( )->GetThermProp( )->GetParam(param);
    if (Section->GetMaterial( )->GetThermProp( )->GetType( ) == MAT_THERMAL_ISOTROPIC)
      Model->GetHeatMod( )->CondMatrix(param,C);
    else
      Model->GetHeatMod( )->CondMatrixOrtho(param,C);

    // Compute element fluxes.

    fl = C*dt;

    // Return the computed values.

    for (int j = 0; j < nsc; j++)
      Flux[i][j] = -fl[j];
  }

  // Release memory.

  delete []coord;
  delete []shpfnc;
  delete []mapfnc;
  delete []shpdrv;
}

// ============================== CalcStrInt ===============================

void cElmParam :: CalcStrInt(double &elmvol, cVector &strint)
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

  // Loop over integration points.

  int strok = 1;
  elmvol = 0.0;
  strint.Zero( );
  for (int i = 0; i < NumIntPnt; i++)
  {
    // Get integration point data.

    sNatCoord p = IntPnt[i].GetCoord( );
    double wgt  = IntPnt[i].GetWeight( );

    // Compute element strains: {e} = [B]{u}.

    double detJ;
    Shape->ShpFunc(p, shpfnc);
    Shape->MapFunc(p, mapfnc);
    Shape->DrvShpXYZ(p, coord, &detJ, shpdrv);
    Model->BMatrix(nn, shpfnc, mapfnc, coord, shpdrv, B);
    e = B*u;

    // Compute the element stresses.

    cVector temp(2);
    int therm = GetStrnTemp(p, shpfnc, temp);
    if (therm)
    {
      double tref = 0;
      if(!SecAn[i]->Stress(tref, temp, e, s))
      {
        strok = 0;
        break;
      }
    }
    else
    {
      if (!SecAn[i]->Stress(e, s))
      {
       strok = 0;
       break;
      }
    }

    // Compute {strint} += coeff*{s}.

    double thk   = Section->GetThickness( );
    double volc  = Model->VolCoeff(thk, nn, mapfnc, coord);
    double coeff = wgt*detJ*volc;
    strint += coeff*s;
    elmvol += coeff;
  }
//  cout << "elmvol = " << fixed << setprecision(3) << elmvol << endl;
//  cout << "{strint} = " << scientific << setprecision(5);
//  strint.Print( );

  // Release memory.

  delete []coord;
  delete []shpfnc;
  delete []mapfnc;
  delete []shpdrv;
}

// ============================= IntFunc ===================================

double cElmParam :: IntFunc(int NumPnt, cIntPoint *Pnts, cVector &Func, double &vol)
{
  // Get the nodal coordinates.
  int nn = Shape->GetNumNode( );
  sNodeCoord *coord = new sNodeCoord[nn];
  Shape->NodalCoord(coord);

  // Get memory for shape functions and derivatives.

  double *shpfnc     = new double [nn];
  double *mapfnc     = new double [nn];
  sNodeCoord *shpdrv = new sNodeCoord[nn];

  // Loop over input integration points
  sNatCoord p;
  double wgt, sum = 0.0;
  for (int i = 0; i < NumPnt; i++)
  {
    // Get integration point data.

    p   = Pnts[i].GetCoord( );
    wgt = Pnts[i].GetWeight( );

    // Compute detJ.

    double detJ;
    Shape->ShpFunc(p, shpfnc);
    Shape->MapFunc(p, mapfnc);
    Shape->DrvShpXYZ(p, coord, &detJ, shpdrv);

    // Compute sum += coeff * Func[i].

    double volc   = Model->VolCoeff(Section->GetThickness( ), nn, mapfnc, coord);
    double coeff  = wgt*detJ*volc;
    vol          += coeff;
    sum          += Func[i] * coeff; 
  }

  // Release memory.

  delete []coord;
  delete []shpfnc;
  delete []mapfnc;
  delete []shpdrv;

  return sum;
}

double cElmParam :: IntFunc(cVector &Func, double &vol)
{
  // Get the nodal coordinates.
  int nn = Shape->GetNumNode( );
  sNodeCoord *coord = new sNodeCoord[nn];
  Shape->NodalCoord(coord);

  // Get memory for shape functions and derivatives.

  double *shpfnc     = new double [nn];
  double *mapfnc     = new double [nn];
  sNodeCoord *shpdrv = new sNodeCoord[nn];

  // Loop over input integration points
  sNatCoord p;
  double sum = 0.0;
  for (int i = 0; i < NumIntPnt; i++)
  {
    // Get integration point data.

    sNatCoord p = IntPnt[i].GetCoord( );
    double wgt  = IntPnt[i].GetWeight( );

    // Compute detJ.

    double detJ;
    Shape->ShpFunc(p, shpfnc);
    Shape->MapFunc(p, mapfnc);
    Shape->DrvShpXYZ(p, coord, &detJ, shpdrv);

    // Compute sum += coeff * Func[i].

    double volc   = Model->VolCoeff(Section->GetThickness( ), nn, mapfnc, coord);
    double coeff  = wgt*detJ*volc;
    vol          += coeff;
    sum          += Func[i] * coeff;
  }

  // Release memory.

  delete []coord;
  delete []shpfnc;
  delete []mapfnc;
  delete []shpdrv;

  return sum;
}

// ============================= PntStress =================================

void cElmParam :: PntStress(sNatCoord *Pnts, int NumPnt, cMatrix &Strain, 
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
    Shape->DrvShpXYZ(p, coord, &detJ, shpdrv);
    Model->BMatrix(nn, shpfnc, mapfnc, coord, shpdrv, B);
    e = B*u;

    // Compute the element stresses.
 
    cVector temp(2);
    int therm = GetStrnTemp(p, shpfnc, temp);
    if (therm)
    {
     double tref = 0;
     SecAn[VecID[i]]->Stress(tref, temp, e, s);
    }
    else
    {
     SecAn[VecID[i]]->Stress(e, s);
    }

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
}

// ============================= NodalStress ===============================

void cElmParam :: NodalStress(cMatrix &Stress)
{
  if (cControl :: GetStrExtFlag( ))
  {
    // Compute the nodal stress => [S]node = [TR]*[S]ipt.

    CompTRMatrix( );
    Stress = TR*StrIpt;
  }
  else
  {
    // Define strain matrix with same size of stress matrix.
    static cMatrix Strain(Stress);

    if ((Strain.NRow() != Stress.NRow()) || (Strain.NCol() != Stress.NCol()))
      Strain.Resize(Stress.NRow( ),Stress.NCol( ));

    // Get nodal parametric coordinates.
    int numpnts = Shape->GetNumNode( );
    sNatCoord *NodPnt = new sNatCoord [numpnts];
    Shape->NodeNatCoord(NodPnt);

    // Evaluate stresses.
    PntStress(NodPnt,numpnts,Strain,Stress,true);

    // Release memory.
    delete []NodPnt;
  }
}

// =========================== InitialStresses =============================

void cElmParam :: InitialStresses(sNatCoord p, cVector &str)
{
  // Begin with zero initial stresses.

  str.Zero( );
  if (!IniStrFlag) return;

  // Check for given element stresses.

  if (IniStr.Dim( ) > 0)
  {
    str = IniStr;
    return;
  }

  // Evaluate the cartesian coordinates.

  int nn = Shape->GetNumNode( );
  sNodeCoord *coord = new sNodeCoord[nn];
  Shape->NodalCoord(coord);
  double *mapfnc = new double[nn];
  Shape->MapFunc(p, mapfnc);

  sNodeCoord c = { 0.0, 0.0, 0.0 };
  int i;
  for (i = 0; i < nn; i++)
  {
    c.x += mapfnc[i]*coord[i].x;
    c.y += mapfnc[i]*coord[i].y;
    c.z += mapfnc[i]*coord[i].z;
  }
  double vcoord = c.y;
  if (VertAxis == Z_AXIS)
    vcoord = c.z;
  else if (VertAxis == X_AXIS)
    vcoord = c.x;

  // Compute the geostatic stresses.

  double sv = 0.0;  // Vertical stresses.
  for (i = 1; i < NumSoil; i++)
  {
    if (vcoord > VecSoil[i].top) break;
    sv += -VecSoil[i-1].gamma*(VecSoil[i-1].top - VecSoil[i].top);
  }
  sv += -VecSoil[i-1].gamma*(VecSoil[i-1].top - vcoord);
  double sh = sv*VecSoil[i-1].k0;  // Lateral stresses.

  // Assign compute values to stress vector.

  if (Model->GetType( ) == AXISYMMETRIC)
  {
    str[0] = sh;
    str[1] = sv;
    str[2] = sh;
  }
  else if (Model->GetType( ) == SOLID)
  {
    str[0] = sh;
    str[1] = sh;
    str[2] = sv;
  }
  else
  {
    str[0] = sh;
    str[1] = sv;
  }

  // Release memory.

  delete []coord;
  delete []mapfnc;
}

// ============================== GetStrnTemp ==============================

int cElmParam :: GetStrnTemp(sNatCoord p, double *N, cVector &temp)
{
  // Get the element temperature data (if exists).

  sTempData  key = {Label, 0, 0};
  sTempData *elm = (sTempData *)bsearch(&key, ElemTemp, NumElemTemp,
                   sizeof(sTempData), sTempDataCmp);

  // Evaluate the current element temperature as the product of the
  // reference temperature by the current load factor.
  double lfac = cControl :: GetTotFactor( );
  if (elm)
  {
    temp[0] = elm->tinf*lfac;
    temp[1] = elm->tsup*lfac;
  
    //cout << "lfac - "<< lfac << "\n";
    //cout << "temp[0] - "<< temp[0] << "\n";
    //cout << "temp[1] - "<< temp[1] << "\n";
  
    return(1);
  }
  return(0);
}


// -------------------------------------------------------------------------
// Protected methods:
//

// ============================= CompTRMatrix ==============================

void cElmParam :: CompTRMatrix(int numipt, cIntPoint *ipt)
{
  // Load input data. 

  ipt    = (ipt)    ? ipt    : IntPnt;
  numipt = (numipt) ? numipt : NumIntPnt;

  int i,j;

  // Check if there is necessary to compute a new extrapolation matrix.

  if ((PrvOrdIdx == OrdIdx) && (PrvShpType == Shape->GetType( )) &&
      (TR.NRow( ) != 0)) return;
  PrvOrdIdx  = OrdIdx;
  PrvShpType = Shape->GetType( );

  // Define the matrix dimensions.

  int nn = Shape->GetNumNode( );
  TR.Resize(nn, numipt);

  // Handle special cases e.g. (order == 1).

  if (VecIntOrd[OrdIdx].K[0] == 1) // TODO: Descobrir como tratar isto considerando multiplas ordens
  {
    double val = 1.0/double(numipt); // Use the average value
    for (i = 0; i < nn; i++)
      for (j = 0; j < numipt; j++) TR[i][j] = val;
    return;
  }

  // Define problem size and create aux. arrays.

  int dim = Shape->NumDim( ) + 1;
  cVector u(4);
  cVector v(dim);
  cVector w(dim);

  // Assembly least-square matrix [A].

  int prof[] = {0, 0, 0, 0};
  cSymSkylMatrix A(dim, prof);
  A.Zero( );
  for (i = 0; i < numipt; i++)
  {
    sNatCoord p = ipt[i].GetCoord( );
    u[0] = 1.0;
    u[1] = p.r;
    u[2] = p.s;
    u[3] = p.t;
    for (int k = 0; k < dim; k++)
      for (int l = 0; l <= k; l++) A(k, l) += u[k]*u[l];
  }

  // Compute [A][Q] = [CI].

  cMatrix Q(dim, numipt);
  for (i = 0; i < numipt; i++)
  {
    sNatCoord p = ipt[i].GetCoord( );
    u[0] = 1.0;
    u[1] = p.r;
    u[2] = p.s;
    u[3] = p.t;
    for (j = 0; j < dim; j++) v[j] = u[j];
    A.Solve(v, w);

    for (j = 0; j < dim; j++) Q[j][i] = w[j];
  }

  // Assembly nodal point coordinate matrix [CN].

  sNatCoord *natcoord = new sNatCoord[nn];
  Shape->NodeNatCoord(natcoord);
  cMatrix CN(nn, dim);
  for (i = 0; i < nn; i++)
  {
    u[0] = 1.0;
    u[1] = natcoord[i].r;
    u[2] = natcoord[i].s;
    u[3] = natcoord[i].t;
    for (j = 0; j < dim; j++) CN[i][j] = u[j];
  }

  // Compute the [TR] = [CN]*inv[A]*[CI] = [CN]*[Q].

  TR = CN*Q;

  // Release memory.

  delete []natcoord;
}

// ============================= NearestIntPnt =============================

void cElmParam :: NearestIntPnt(int dim ,sNatCoord *Pnts, int NumPnt,
                                cIntPoint *IntPnt, int NumIntPnt, 
                                vector<int> &VecID)
{
  if (static_cast<int>(VecID.size( )) < NumPnt)
    VecID.resize(NumPnt);

  double dist, mindist;
  sNatCoord p, ip;
  for (int i = 0; i < NumPnt; i++)
  {
    // Set initial point equal to zero and compute minimum distance.
    p        = Pnts[i];
    ip       = IntPnt[0].GetCoord( );
    VecID[i] = 0;

    mindist = 0;
    switch (dim)
    {
      case (3):
      mindist  += pow(p.t-ip.t,2.0);
      case (2):
      mindist  += pow(p.s-ip.s,2.0);
      case (1):
      mindist  += pow(p.r-ip.r,2.0);
    }

    // Loop over others integration points.
    for (int j = 1; j < NumIntPnt; ++j)
    {
      // Get integration point data.

      ip   = IntPnt[j].GetCoord( );

      // Evaluate distance
      dist = 0;
      switch (dim)
      {
        case (3):
        dist  += pow(p.t-ip.t,2.0);
        case (2):
        dist  += pow(p.s-ip.s,2.0);
        case (1):
        dist  += pow(p.r-ip.r,2.0);
      }

      // Update the current nearest integration point index.

      if (dist < mindist)
      {
        mindist = dist;
        VecID[i] = j;
      }
    }
  }
}


// -------------------------------------------------------------------------
// Methods of cElmParamTL class:
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Public methods:
//

// ============================= cElmParamTL ===============================

cElmParamTL :: cElmParamTL(int id, eShpType tshp, eAnmType tanm)
             : cElmParam(id, tshp, tanm)
{
  Type = PARAMETRICTL;
}

// =============================== IntForce ================================

int cElmParamTL :: IntForce(cVector &g)
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
  sNodeCoord *shpdrv = new sNodeCoord[nshp];

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
    Shape->DrvShpXYZ(p, coord, &detJ, shpdrv);
    Model->BMatrix(nn, shpfnc, mapfnc, coord, shpdrv, B);
    Model->GetMecMod( )->BlMatrix(nn, shpfnc, mapfnc, u, coord, shpdrv, Bl);
    GreenStrain(B, Bl, u, e);

    // Compute the element stresses.

    cVector temp(2);
    int therm = GetStrnTemp(p, shpfnc, temp);
    if (therm)
    {
      double tref = 0;
      if(!SecAn[i]->Stress(tref, temp, e, s))
      {
        strok = 0;
        break;
      }
    }
    else
    {
      if (!SecAn[i]->Stress(e, s))
      {
        strok = 0;
        break;
      }
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
  delete []shpdrv;

  return(strok);
}

// ============================== StiffMat ================================

void cElmParamTL :: StiffMat(cMatrix &Kt)
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
    Shape->DrvShpXYZ(p, coord, &detJ, shpdrv);
    Model->BMatrix(nn, shpfnc, mapfnc, coord, shpdrv, B);
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
  }

  // Release memory.

  delete []coord;
  delete []shpfnc;
  delete []mapfnc;
  delete []shpdrv;
 }

// ============================ IntPntStress ===============================

void cElmParamTL :: IntPntStress(cMatrix &Strain, cMatrix &Stress)
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

  // Loop over integration points.

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
    Shape->DrvShpXYZ(p, coord, &detJ, shpdrv);
    Model->BMatrix(nn, shpfnc, mapfnc, coord, shpdrv, B);
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
}

// -------------------------------------------------------------------------
// Protected methods:
//

// ============================ GreenStrain ================================

void cElmParamTL :: GreenStrain(cMatrix &B, cMatrix &Bl, cVector &u, cVector &str)
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
