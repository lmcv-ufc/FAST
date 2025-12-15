// -------------------------------------------------------------------------
// elmdsh.cpp - implementation of degenerated shell element.
// -------------------------------------------------------------------------
// Created:      30-Sep-2020     Elias Saraiva Barroso
//
// -------------------------------------------------------------------------

#include <cmath>
#include <iostream>
#include <iomanip>
#include <vector>

using namespace std;

#include "elmdsh.h"
#include "ctrl.h"
#include "gblvar.h"
#include "shpshell.h"
#include "utl.h"

// -------------------------------------------------------------------------
// Static variables:
//
static bool PrintNormal = false;
static int  NumDegSh    = 0;
static int  UpdateOption = 4;   // Temporary (options = 1, 2, 3, 4, 5)

// -------------------------------------------------------------------------
// Class cElmDegShell:
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Public methods:
//

// ============================= CompDegShlAxes ============================

void cElmDegShell :: CompDegShlAxes( )
{
  // Abort if there is no shell elements.
  if (!NumDegSh) return; 

  // Auxiliary variables.
  int nn, label;
  int maxn = cElement :: GetMaxNode( )+1;
  int gnn  = cNode    :: GetNumNode( );
  sNodeCoord *NodCoord = new sNodeCoord[maxn];
  cVector    *GlbNodNorm = new cVector[gnn]; 
  cVector    *GlbNodV1 = new cVector[gnn];
  cVector    *GlbNodV2 = new cVector[gnn];
  sNodeAxes  *ElmAxes  = new sNodeAxes[maxn];

  // Initialize global nodal normal.
  for(int n = 0; n < gnn; ++n)
  {
    GlbNodNorm[n].Resize(3);
    GlbNodNorm[n].Zero( );
    
    GlbNodV1[n].Resize(3);
    GlbNodV1[n].Zero( );
    GlbNodV2[n].Resize(3);
    GlbNodV2[n].Zero( );
  }

  // Add element normal to global nodal normals vector.
  int ne = cElement :: GetNumElm( );
  //cElmDegShell *elm;
  cShapeShell *shp;
  for(int e = 0; e < ne; ++e)
    if ((shp = dynamic_cast<cShapeShell*>(cElement :: GetElm(e)->GetShape( ))))
    {
      // Get element normal vectors.
      shp->NodalCoord(NodCoord);

      nn = shp->GetNumNode( );
      for(int n = 0; n < nn; ++n)
      {
        label                  = shp->GetNode(n)->GetLabel( );
        GlbNodNorm[label-1]   += NodCoord[n].axes->V3;
      }
    }

  // Normalize nodal normals.
  for(int n = 0; n < gnn; ++n) GlbNodNorm[n].Normalize( );

  // Set element normals.
  for(int e = 0; e < ne; ++e)
    if ((shp = dynamic_cast<cShapeShell*>(cElement :: GetElm(e)->GetShape( ))))
    {
      // Get smoothed normal vectors.
      nn = shp->GetNumNode( );
      for(int n = 0; n < nn; ++n)
      {
        label         = shp->GetNode(n)->GetLabel( );
        ElmAxes[n].V3 = GlbNodNorm[label-1];
      }

      // Setup element normal.
      shp->SetNodalAxes(ElmAxes);

//      if (cElement :: GetElm(e)->GetType( ) == PARAMETRICTL) // Elias (12-Jun-2024)
        shp->LoadInitAxes( );
      
      // A partir daqui pode imprimir no .pos
    }
    
    
    
  //------------------------------------------------------------------------------------------------------------  
  // Set element normals.
  for(int e = 0; e < ne; ++e)
    if ((shp = dynamic_cast<cShapeShell*>(cElement :: GetElm(e)->GetShape( ))))
    {
    	
      // Get element normal vectors.
      shp->NodalCoord(NodCoord);
		
      nn = shp->GetNumNode( );
      for(int n = 0; n < nn; ++n)
      {
        label         = shp->GetNode(n)->GetLabel( );
        GlbNodV1[label-1]   += NodCoord[n].axes->V1;
        GlbNodV2[label-1]   += NodCoord[n].axes->V2;
      }
    }    
     
  // Normalize nodal normals.
  for(int n = 0; n < gnn; ++n) GlbNodV1[n].Normalize( );
  for(int n = 0; n < gnn; ++n) GlbNodV2[n].Normalize( );
  
  //cout << "Print normal -------------------------------------------------------------------";
  if (PrintNormal) PrintNormalVec(gnn,GlbNodV1,GlbNodV2,GlbNodNorm);
    
   //------------------------------------------------------------------------------------------------------------  

  // Release memory.
  delete []GlbNodNorm;  
  delete []NodCoord;
  delete []ElmAxes;
  delete []GlbNodV1;
  delete []GlbNodV2;
  // cout << "saiu comp normal" << endl;
}

// ============================= PrintNormal ===============================

void cElmDegShell :: PrintNormalVec(int nn, cVector *gV1, cVector *gV2, cVector *gV3)
{
  out << "%NODE.VECTOR\n" << nn << "\n";
  
  for(int n = 0; n < nn; ++n)
  {
    out << n+1           << "     ";
    out << showpos << scientific << setprecision(6);
    out << gV1[n][0] << "   ";
    out << gV1[n][1] << "   ";
    out << gV1[n][2] << "   ";
    out << gV2[n][0] << "   ";
    out << gV2[n][1] << "   ";
    out << gV2[n][2] << "   ";
    out << gV3[n][0] << "   ";
    out << gV3[n][1] << "   ";
    out << gV3[n][2] << "\n";
    out << resetiosflags(ios::scientific | ios::showpos);
  }
  out << "\n";

  //cout << "saiu do PrintNOrmal " << endl;
}
 
// ============================= cElmDegShell ==============================

cElmDegShell :: cElmDegShell(int id, eShpType tshp, eAnmType tanm)
             : cElmParam(id, tshp, tanm)
{
  NumDegSh++;
}

// ============================= ~cElmDegShell =============================

cElmDegShell :: ~cElmDegShell( )
{
  NumDegSh--;
}

// ============================= NodalStress ===============================

void cElmDegShell :: NodalStress(cMatrix &Stress)
{
  // Extract quadrature information from zeta direction.
  static vector<double> vecW;
  static vector<double> vecJ;

  vecW.clear( );
  vecJ.clear( );
  int sizeRS = 0;

  for(int i = 0; i < NumIntPnt; ++i)
  {
    if (vecJ.size( ) == 1) ++sizeRS;
    
    if (vecJ.empty( ) || fabs(vecJ.back( ) - IntPnt[i].GetCoord( ).t) > 1e-6)
    {
      vecJ.push_back(IntPnt[i].GetCoord( ).t); 
      vecW.push_back(IntPnt[i].GetWeight( )); 
    }
    else
      vecW.back( ) += IntPnt[i].GetWeight( ); 
  }
  int jsize         = vecW.size( );
  if (Shape->GetTopologyType( ) == WEDGE_TOPOLOGY)
    for (int i = 0; i < jsize; ++i) vecW[i] *= 2.0;
  else
    for (int i = 0; i < jsize; ++i) vecW[i] /= 4.0;  // Hex case.

  // Evaluate extrapolation matrix.
  if (cControl :: GetStrExtFlag( )) CompTRMatrix(sizeRS);

  // Get nodal coordinates.
  int nn  = Shape->GetNumNode( );
  sNatCoord  *NodPnt   = new sNatCoord [nn];
  sNodeCoord *NodCoord = new sNodeCoord[nn];
  Shape->NodeNatCoord(NodPnt);
  Shape->NodalCoord(NodCoord);

  // Create local matrices.
  int nsc = Model->GetDimBMatrix( );
  cMatrix strain(nn,nsc);
  cMatrix stress(nn,nsc);
  cMatrix zetaStrIpt(sizeRS,nsc);

  Stress.Zero( );
  for(int p = 0; p < jsize; ++p)
  {
    // Load current zeta natcoord.
    double zeta = vecJ[p];
    for(int n = 0; n < nn; ++n) NodPnt[n].t = zeta;
      
    // Evaluate stresses in element local system.
    if (cControl :: GetStrExtFlag( ))
    {
      // Compute the nodal stress through int. pnt. extrapolation.
      for(int i = 0; i < sizeRS; ++i)
        for(int j = 0; j < nsc; ++j)
          zetaStrIpt[i][j] = StrIpt[sizeRS*p+i][j];

      stress = TR*zetaStrIpt;
    }
    else
      PntStress(NodPnt,nn,strain,stress);

    for(int n = 0; n < nn; ++n)
    {
      // Add contribution of each node stress component on force integration.
      double sxx   = stress[n][0];
      double syy   = stress[n][1];
      double sxy   = stress[n][3];
      double sxz   = stress[n][4];
      double syz   = stress[n][5];
      double w     = vecW[p];
      double t     = *(NodCoord[n].thk);

      Stress[n][0] += t * w * sxx / 2.0;
      Stress[n][1] += t * w * syy / 2.0;
      Stress[n][2] += t * w * sxy / 2.0;
      Stress[n][3] -= t * t * w * zeta * sxx / 4.0;
      Stress[n][4] -= t * t * w * zeta * syy / 4.0;
      Stress[n][5] += t * t * w * zeta * sxy / 4.0;
      Stress[n][6] += t * w * sxz / 2.0;
      Stress[n][7] += t * w * syz / 2.0;            
    }
  }                        
                           
  // Release memory.
  delete []NodPnt;
  delete []NodCoord;
}

// -------------------------------------------------------------------------
// Class cElmDegShellTL:
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Public methods:
//

// ============================= ReadUpdateOption ==========================

void cElmDegShellTL :: ReadUpdateOption(void)
{
  // Read the update option label.

  char label[100];
  if (!Utl::ReadString(in, label))
  {
    cout << "Error in the input of the algorithm label!\n";
    exit(0);
  }

  // Choose the appropriate update option.

  if (string(label) == "MR" || string(label) == "mr")
    UpdateOption = 1;
  else if (string(label) == "LU" || string(label) == "lu")
    UpdateOption = 2;
  else if (string(label) == "LUN" || string(label) == "lun")
    UpdateOption = 3;
  else if (string(label) == "QU" || string(label) == "qu")
    UpdateOption = 4;
  else
  {
    cout << "Unknown update option: " << label << "\n";
    exit(0);
  }
}

// =========================== cElmDegShellTL ==============================

cElmDegShellTL :: cElmDegShellTL(int id, eShpType tshp, eAnmType tanm)
             : cElmParamTL(id, tshp, tanm)
{
  LargeRot3D = 1;
  NumDegSh++;
}

// ============================= ~cElmDegShellTL ===========================

cElmDegShellTL :: ~cElmDegShellTL( )
{
  NumDegSh--;
}

// ============================= NodalStress ===============================

void cElmDegShellTL :: NodalStress(cMatrix &Stress)
{
  // Extract quadrature information from zeta direction.
  static vector<double> vecW;
  static vector<double> vecJ;

  vecW.clear( );
  vecJ.clear( );
  int sizeRS = 0;

  for(int i = 0; i < NumIntPnt; ++i)
  {
    if (vecJ.size( ) == 1) ++sizeRS;
    
    if (vecJ.empty( ) || fabs(vecJ.back( ) - IntPnt[i].GetCoord( ).t) > 1e-6)
    {
      vecJ.push_back(IntPnt[i].GetCoord( ).t); 
      vecW.push_back(IntPnt[i].GetWeight( )); 
    }
    else
      vecW.back( ) += IntPnt[i].GetWeight( ); 
  }
  int jsize         = vecW.size( );
  if (Shape->GetTopologyType( ) == WEDGE_TOPOLOGY)
    for (int i = 0; i < jsize; ++i) vecW[i] *= 2.0;
  else
    for (int i = 0; i < jsize; ++i) vecW[i] /= 4.0;  // Hex case.

  // Evaluate extrapolation matrix.
  if (cControl :: GetStrExtFlag( )) CompTRMatrix(sizeRS);

  // Get nodal coordinates.
  int nn  = Shape->GetNumNode( );
  sNatCoord  *NodPnt   = new sNatCoord [nn];
  sNodeCoord *NodCoord = new sNodeCoord[nn];
  Shape->NodeNatCoord(NodPnt);
  Shape->NodalCoord(NodCoord);

  // Create local matrices.
  int nsc = Model->GetDimBMatrix( );
  cMatrix strain(nn,nsc);
  cMatrix stress(nn,nsc);
  cMatrix zetaStrIpt(sizeRS,nsc);

  Stress.Zero( );
  for(int p = 0; p < jsize; ++p)
  {
    // Load current zeta natcoord.
    double zeta = vecJ[p];
    for(int n = 0; n < nn; ++n) NodPnt[n].t = zeta;
      
    // Evaluate stresses in element local system.
    if (cControl :: GetStrExtFlag( ))
    {
      // Compute the nodal stress through int. pnt. extrapolation.
      for(int i = 0; i < sizeRS; ++i)
        for(int j = 0; j < nsc; ++j)
          zetaStrIpt[i][j] = StrIpt[sizeRS*p+i][j];

      stress = TR*zetaStrIpt;
    }
    else
      PntStress(NodPnt,nn,strain,stress);

    for(int n = 0; n < nn; ++n)
    {
      // Add contribution of each node stress component on force integration.
      double sxx   = stress[n][0];
      double syy   = stress[n][1];
      double sxy   = stress[n][3];
      double sxz   = stress[n][4];
      double syz   = stress[n][5];
      double w     = vecW[p];
      double t     = *(NodCoord[n].thk);

      Stress[n][0] += t * w * sxx / 2.0;
      Stress[n][1] += t * w * syy / 2.0;
      Stress[n][2] += t * w * sxy / 2.0;
      Stress[n][3] -= t * t * w * zeta * sxx / 4.0;
      Stress[n][4] -= t * t * w * zeta * syy / 4.0;
      Stress[n][5] += t * t * w * zeta * sxy / 4.0;
      Stress[n][6] += t * w * sxz / 2.0;
      Stress[n][7] += t * w * syz / 2.0;            
    }
  }                        
                           
  // Release memory.

  delete []NodPnt;
  delete []NodCoord;
}

// =============================== IntForce ================================

int cElmDegShellTL :: IntForce(cVector &g)
{
//  cout << "\ncElmDegShellTL :: IntForce\n";

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
  cVector e(nsc);    // Strain
  cVector s(nsc);    // Stress
  cMatrix B0(nsc, ndof);
  cMatrix Bl(nsc, ndof);
  cMatrix Bb(nsc, ndof);

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

    // Compute the element strains.

    double detJ;
    Shape->ShpFunc(p, shpfnc);
    Shape->MapFunc(p, mapfnc);
    Shape->DrvShpXYZ(p, coord, &detJ, shpdrv);    

    // Evalute the strain-displacement matrices.

    Model->BMatrix(nn, shpfnc, mapfnc, coord, shpdrv, B0);
    Model->GetMecMod( )->BlMatrix(nn, shpfnc, mapfnc, u, coord, shpdrv, Bl);   
    Bb = B0 + Bl;

    // Compute the strains and stresses.

    GreenStrain(nn, shpfnc, u, coord, shpdrv, e);
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
//    cout << "{e} = " << scientific << setprecision(6) << showpos;
//    e.Print();
//    cout << "{s} = " << scientific << setprecision(6) << showpos;
//    s.Print();
//    cout << endl;

    // Compute [g] += coeff*[Bb]t{s}.

    double volc  = Model->VolCoeff(Section->GetThickness( ), nn, mapfnc, coord);
    double coeff = wgt*detJ*volc;
    MultTAcc(coeff, Bb, s, g);
  }
//  cout << "{g} = " << scientific << setprecision(4) << showpos;
//  g.Print();

  // Release memory.

  delete []coord;
  delete []shpfnc;
  delete []mapfnc;
  delete []shpdrv;
//  cout << "End of cElmDegShellTL :: IntForce\n\n";

  return(strok);
}

// ============================== StiffMat ================================

void cElmDegShellTL :: StiffMat(cMatrix &Kt)
{
//  cout << "\ncElmDegShellTL :: StiffMat\n";
 
  // Get the nodal coordinates.

  int nn = Shape->GetNumNode( );
  sNodeCoord *coord = new sNodeCoord[nn];
  Shape->NodalCoord(coord);

  // Create local arrays.

  int ndof = nn*Model->GetNumDofNode( );
  int nsc = Model->GetDimBMatrix( );
  int bnldim = Model->GetMecMod( )->GetDimBnlMatrix();
  cVector e(nsc);
  cVector s(nsc);
  cMatrix C(nsc, nsc);
  cMatrix B0(nsc, ndof);
  cMatrix Bb(nsc, ndof);
  cMatrix Bl(nsc, ndof);
  cMatrix Bnl(bnldim, ndof);
  cMatrix S(bnldim, bnldim);
  cMatrix Kgnl(ndof, ndof);

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

    // Compute the [B] matrices.

    double detJ;
    Shape->ShpFunc(p, shpfnc);
    Shape->MapFunc(p, mapfnc);
    Shape->DrvShpXYZ(p, coord, &detJ, shpdrv);
    Model->BMatrix(nn, shpfnc, mapfnc, coord, shpdrv, B0);
    Model->GetMecMod( )->BlMatrix(nn, shpfnc, mapfnc, u, coord, shpdrv, Bl);
    Model->GetMecMod( )->BnlMatrix(nn, shpfnc, mapfnc, coord, shpdrv, Bnl);
    Bb = B0 + Bl;  // Incremental strain-displacement matrix

    // Get [C] matrix.

    SecAn[i]->CMatrix(C);

    // Evaluate the Green-Lagrande strains and Piola-Kirchhoff II stresses.

    GreenStrain(nn, shpfnc, u, coord, shpdrv, e);
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

    // Compute [K] += coeff*([Bb]t[C][Bb] + [Bnl]t[S][Bnl]).

    double volc  = Model->VolCoeff(Section->GetThickness( ), nn, mapfnc, coord);
    double coeff = wgt*detJ*volc;
    MatTripBtCB(Bb, C, coeff, Kt);
    MatTripBtCB(Bnl, S, coeff, Kt);

    // Add the nonlinear geometric stiffness matrix.

    KgnlMatrix(nn, shpfnc, coord, shpdrv, u, s, Kgnl);
    Kt += coeff*Kgnl;
   
//    static double vol = 0;
//    vol += coeff;
//    cout << "volume = " << vol << endl;
  }

  // Release memory.

  delete []coord;
  delete []shpfnc;
  delete []mapfnc;
  delete []shpdrv;
//  cout << "End of cElmDegShellTL :: StiffMat\n\n";
}

// ============================ IntPntStress ===============================

void cElmDegShellTL :: IntPntStress(cMatrix &Strain, cMatrix &Stress)
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
    Shape->DrvShpXYZ(p, coord, &detJ, shpdrv);
    Model->BMatrix(nn, shpfnc, mapfnc, coord, shpdrv, B);
    Model->GetMecMod( )->BlMatrix(nn, shpfnc, mapfnc, u, coord, shpdrv, Bl);
    GreenStrain(nn, shpfnc, u, coord, shpdrv, e);

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

// ============================== UpdateItera ==============================

void RotationMatrix(double alpha, double beta, cMatrix &Q)
{
  // Auxiliary terms.
  
  double theta = sqrt(alpha*alpha + beta*beta);
  double k1 = 1.0;
  double k2 = 0.5;
  if (theta > 1.0e-12)
  {	  
    k1 = sin(theta)/theta;
    k2 = (1.0 - cos(theta))/(theta*theta);
  }

  // Rotation matrix.
  
  Q[0][0] = 1.0 - k2*beta*beta;
  Q[0][1] = k2*alpha*beta;
  Q[0][2] = k1*beta;
  Q[1][1] = 1.0 - k2*alpha*alpha;
  Q[1][2] = -k1*alpha;
  Q[2][2] = 1.0 - k2*(alpha*alpha + beta*beta);
  Q[1][0] =  Q[0][1];
  Q[2][0] = -Q[0][2];
  Q[2][1] = -Q[1][2];
}	

void cElmDegShellTL :: UpdateItera(cVector &dug)
{
  // Get the nodal directors.

  int nn = Shape->GetNumNode( );
  sNodeCoord *coord = new sNodeCoord[nn];
  Shape->NodalCoord(coord);

  // Extract element local vector from the global one.
  
  int ndof = nn*Model->GetNumDofNode( );
  cVector due(ndof);
  GlobToElm(dug, due);
 
  // Evaluate updated normal vectors.

  static double maxlen3 = 1.0;  // Temporary
  sNodeAxes *naxes = new sNodeAxes[nn]; 
  for (int i = 0; i < nn; i++)
  {
    double thk2 = *(coord[i].thk)/2.0;
    double alpha = due[5*i+3];    // Iterative rotation
    double beta  = due[5*i+4];    // Iterative rotation
    naxes[i] = *(coord[i].axes);  // Current nodel axes

    // Update {V3} and evaluate the derivatives of {p}.

    double L;  // |V3|
    if (UpdateOption == 1)  // Moderate rotations
    {
      // Update {V3}.

      sNodeAxes iaxes = *(coord[i].iaxes);
      for (int j = 0; j < 3; j++)
      {
        naxes[i].V3[j] += iaxes.V1[j]*beta - iaxes.V2[j]*alpha;
      }
      L = naxes[i].V3.Length( );

      // Evalute the first derivatives of {p}.      

      *coord[i].pa = -thk2*coord[i].iaxes->V2; // d{p}/dalpha
      *coord[i].pb =  thk2*coord[i].iaxes->V1; // d{p}/dbeta

      // Evalute the second derivatives of {p}.      

      (*coord[i].paa).Zero( );
      (*coord[i].pab).Zero( );
      (*coord[i].pbb).Zero( );
    }
    else if (UpdateOption == 2)  // Linear update
    {
      // Update {V3}.

      for (int j = 0; j < 3; j++)
      {
        naxes[i].V3[j] += naxes[i].V1[j]*beta - naxes[i].V2[j]*alpha;
      }
      L = naxes[i].V3.Length( );
     
      // Evalute the first derivatives of {p}.      
     
      *coord[i].pa = -thk2*naxes[i].V2; // d{p}/dalpha
      *coord[i].pb =  thk2*naxes[i].V1; // d{p}/dbeta

      // Evalute the second derivatives of {p}.      

      (*coord[i].paa).Zero( );
      (*coord[i].pab).Zero( );
      (*coord[i].pbb).Zero( );
    }
    else if (UpdateOption == 3)  // Linear update with normalization
    {
      // Update {V3}.

      for (int j = 0; j < 3; j++)
      {
        naxes[i].V3[j] += naxes[i].V1[j]*beta - naxes[i].V2[j]*alpha;
      }
      L = naxes[i].V3.Length( );
      naxes[i].V3 /= L;          // Normalization
     
      // Evaluate the derivatives of {p}.
     
      double La = alpha/L;
      double Lb =  beta/L;
      double Laa = (1.0 - La*La)/L;
      double Lbb = (1.0 - Lb*Lb)/L;
      double Lab = -La*Lb/L;
      for (int j = 0; j < 3; j++)
      {
	// Evalute the first derivatives of {p}.      
     
        double V3aj = (-naxes[i].V2[j] - La*naxes[i].V3[j])/L;
        double V3bj = ( naxes[i].V1[j] - Lb*naxes[i].V3[j])/L;
        (*coord[i].pa)[j] = V3aj*thk2;
        (*coord[i].pb)[j] = V3bj*thk2;
     
	// Evalute the second derivatives of {p}.      
     
        double V3aaj = -(2.0*La*V3aj + Laa*naxes[i].V3[j])/L;
        double V3bbj = -(2.0*Lb*V3bj + Lbb*naxes[i].V3[j])/L;
        double V3abj = -(La*V3bj + Lb*V3aj + Lab*naxes[i].V3[j])/L;
        (*coord[i].paa)[j] = V3aaj*thk2; 
        (*coord[i].pab)[j] = V3abj*thk2; 
        (*coord[i].pbb)[j] = V3bbj*thk2; 
      }
    }
    else if (UpdateOption == 4)  // Quadratic update (Bathe)
    {
      double k1 = 0.5*(alpha*alpha + beta*beta);
      for (int j = 0; j < 3; j++)
      {
	// Evalute the first derivatives of {p}.      
     
        double V3aj = -alpha*naxes[i].V3[j] - naxes[i].V2[j];
        double V3bj =  -beta*naxes[i].V3[j] + naxes[i].V1[j];
        (*coord[i].pa)[j] = V3aj*thk2;
        (*coord[i].pb)[j] = V3bj*thk2;
     
	// Evalute the second derivatives of {p}.      
     
        double V3aaj = -naxes[i].V3[j];
        double V3abj = 0.0;
        double V3bbj = -naxes[i].V3[j];
        (*coord[i].paa)[j] = V3aaj*thk2; 
        (*coord[i].pab)[j] = V3abj*thk2; 
        (*coord[i].pbb)[j] = V3bbj*thk2; 
     
	// Update {V3}.
     
        naxes[i].V3[j] += beta*naxes[i].V1[j] - alpha*naxes[i].V2[j] - k1*naxes[i].V3[j];
      }
      L = naxes[i].V3.Length( );
    }
    else  // Rotation matrix
    {
      cout << "Rotation matrix update not implemented yet!\n";
      cMatrix Q(3, 3);
      RotationMatrix(alpha, beta, Q);
      cout << "[Q]\n";
      Q.Print( );

      cVector V3n(3);
      V3n = Q*naxes[i].V3;
      naxes[i].V3 = V3n;
      cout << "{V3n} = ";
      naxes[i].V3.Print( );
      cout << "\n";
      L = naxes[i].V3.Length( );

      exit(0); 
    }

    // Temporary: print the length of {V3}.

    if (L > maxlen3)
    {
      maxlen3 = L;
      //cout << "MaxLenV3 = " << fixed << setprecision(8) << maxlen3 << endl;
    }

    // Update {V1} and {V2}.

    double ex[] = {1.0, 0.0, 0.0};
    double ey[] = {0.0, 1.0, 0.0};
    cVector e1 = cVector(3, ex);
    cVector e2 = cVector(3, ey);
    CrossProd(e2, naxes[i].V3, naxes[i].V1);             // V1 = e2 x V3
    if (naxes[i].V1.Length( ) > 1.0e-6)
    {	    
      CrossProd(naxes[i].V3, naxes[i].V1, naxes[i].V2);  // V2 = V3 x V1
    }
    else
    {	    
      cout << "|V1| = " << naxes[i].V1.Length( ) << endl;
      CrossProd(naxes[i].V3, e1, naxes[i].V2);           // V2 = V3 x e1
      CrossProd(naxes[i].V2, naxes[i].V3, naxes[i].V1);  // V1 = V2 x V3 
    }
    naxes[i].V1.Normalize();	
    naxes[i].V2.Normalize();	

//    cout << "{V1}n = ";
//    naxes[i].V1.Print();
//    cout << "{V2}n = ";
//    naxes[i].V2.Print();
//    cout << "{V3}n = ";
//    naxes[i].V3.Print();
//    cout << endl;
  }
  
  // Update director vectors.
  
  Shape->UpdNodalAxes(naxes);
  
  // Release memory.
  
  delete []coord;
  delete []naxes;
}


// -------------------------------------------------------------------------
// Protected methods:
//

// ============================ GreenStrain ================================

void cElmDegShellTL :: GreenStrain(int nn, double *shpfunc, cVector &disp, 
                                   sNodeCoord *coord, sNodeCoord *drvshp,
                                   cVector &e)
{	
  // Get the through-thickness coordinate and its derivatives.

  double zeta = shpfunc[nn];
  double dzetadx = drvshp[nn].x;
  double dzetady = drvshp[nn].y;
  double dzetadz = drvshp[nn].z;

  // Compute the displacement derivatives.

  double ux = 0.0;
  double uy = 0.0;
  double uz = 0.0;
  double vx = 0.0;
  double vy = 0.0;
  double vz = 0.0;
  double wx = 0.0;
  double wy = 0.0;
  double wz = 0.0;
  for (int i = 0; i < nn; i++)
  {
    double u = disp[5*i];
    double v = disp[5*i+1];
    double w = disp[5*i+2];
    
    double thk2 = *(coord[i].thk)/2.0;
    double dHdx = drvshp[i].x*zeta + shpfunc[i]*dzetadx;
    double dHdy = drvshp[i].y*zeta + shpfunc[i]*dzetady;
    double dHdz = drvshp[i].z*zeta + shpfunc[i]*dzetadz;
    
    double px = (coord[i].axes->V3[0] - coord[i].iaxes->V3[0])*thk2;
    double py = (coord[i].axes->V3[1] - coord[i].iaxes->V3[1])*thk2;
    double pz = (coord[i].axes->V3[2] - coord[i].iaxes->V3[2])*thk2;

    ux += drvshp[i].x*u + dHdx*px;
    uy += drvshp[i].y*u + dHdy*px;
    uz += drvshp[i].z*u + dHdz*px;
                   
    vx += drvshp[i].x*v + dHdx*py;
    vy += drvshp[i].y*v + dHdy*py;
    vz += drvshp[i].z*v + dHdz*py;
                   
    wx += drvshp[i].x*w + dHdx*pz;
    wy += drvshp[i].y*w + dHdy*pz;
    wz += drvshp[i].z*w + dHdz*pz;
  }
//  cout << "{beta} = " << scientific << setprecision(6) << showpos;
//  cout << ux << " " << uy << " " << uz << " ";
//  cout << vx << " " << vy << " " << vz << " ";
//  cout << wx << " " << wy << " " << wz << "\n";

  // Evaluate the Green-Lagrange strains.
  
  e[0] = ux;
  e[1] = vy;
  e[2] = wz;
  e[3] = uy + vx; 
  e[4] = wx + uz; 
  e[5] = vz + wy; 
  
//  cout << "{str1} = " << scientific << setprecision(4) << showpos;
//  e.Print();
  
  e[0] += (ux*ux + vx*vx + wx*wx)/2.0;
  e[1] += (uy*uy + vy*vy + wy*wy)/2.0;
  e[2] += (uz*uz + vz*vz + wz*wz)/2.0;
  e[3] += (ux*uy + vx*vy + wx*wy); 
  e[4] += (ux*uz + vx*vz + wx*wz); 
  e[5] += (uy*uz + vy*vz + wy*wz);
  
//  cout << "{str2} = " << scientific << setprecision(4) << showpos;
//  e.Print();
//  cout << "\n";
}

// ============================= KgnlMatrix ================================

void cElmDegShellTL :: KgnlMatrix(int nn, double *shpfunc, sNodeCoord *coord, 
                                  sNodeCoord *drvshp, cVector &disp,
                                  cVector &s, cMatrix &Kgnl)
{
  // Get the through-thickness coordinate and its derivatives.

  double zeta = shpfunc[nn];
  double dzetadx = drvshp[nn].x;
  double dzetady = drvshp[nn].y;
  double dzetadz = drvshp[nn].z;

  // Compute the displacement derivatives.

  double ux = 0.0;
  double uy = 0.0;
  double uz = 0.0;
  double vx = 0.0;
  double vy = 0.0;
  double vz = 0.0;
  double wx = 0.0;
  double wy = 0.0;
  double wz = 0.0;
  for (int i = 0; i < nn; i++)
  {
    double u = disp[5*i];
    double v = disp[5*i+1];
    double w = disp[5*i+2];
    
    double thk2 = *(coord[i].thk)/2.0;
    double dHdx = drvshp[i].x*zeta + shpfunc[i]*dzetadx;
    double dHdy = drvshp[i].y*zeta + shpfunc[i]*dzetady;
    double dHdz = drvshp[i].z*zeta + shpfunc[i]*dzetadz;
    
    double px = (coord[i].axes->V3[0] - coord[i].iaxes->V3[0])*thk2;
    double py = (coord[i].axes->V3[1] - coord[i].iaxes->V3[1])*thk2;
    double pz = (coord[i].axes->V3[2] - coord[i].iaxes->V3[2])*thk2;

    ux += drvshp[i].x*u + dHdx*px;
    uy += drvshp[i].y*u + dHdy*px;
    uz += drvshp[i].z*u + dHdz*px;
                   
    vx += drvshp[i].x*v + dHdx*py;
    vy += drvshp[i].y*v + dHdy*py;
    vz += drvshp[i].z*v + dHdz*py;
                   
    wx += drvshp[i].x*w + dHdx*pz;
    wy += drvshp[i].y*w + dHdy*pz;
    wz += drvshp[i].z*w + dHdz*pz;
  }

  // Get stresses.

  double sxx = s[0];
  double syy = s[1];
  double szz = s[2];
  double sxy = s[3];
  double sxz = s[4];
  double syz = s[5];

  // Evaluate {As} = [A]{s}.

  double As1 = sxx*ux + sxy*uy + sxz*uz;
  double As2 = sxy*ux + syy*uy + syz*uz;
  double As3 = sxz*ux + syz*uy + szz*uz;
  double As4 = sxx*vx + sxy*vy + sxz*vz;
  double As5 = sxy*vx + syy*vy + syz*vz;
  double As6 = sxz*vx + syz*vy + szz*vz;
  double As7 = sxx*wx + sxy*wy + sxz*wz;
  double As8 = sxy*wx + syy*wy + syz*wz;
  double As9 = sxz*wx + syz*wy + szz*wz;

  // Compute the [Kgnl] matrix.

  Kgnl.Zero( );
  for (int i = 0; i < nn; i++)
  {
    int ia = 5*i + 3;  // Index of alpha nodal rotation
    int ib = 5*i + 4;  // Index of beta  nodal rotation
    
    //double thk2 = *(coord[i].thk)/2.0;
    double dHdx = drvshp[i].x*zeta + shpfunc[i]*dzetadx;
    double dHdy = drvshp[i].y*zeta + shpfunc[i]*dzetady;
    double dHdz = drvshp[i].z*zeta + shpfunc[i]*dzetadz;
    
    double paax = (*coord[i].paa)[0];
    double paay = (*coord[i].paa)[1];
    double paaz = (*coord[i].paa)[2];
    double pabx = (*coord[i].pab)[0];
    double paby = (*coord[i].pab)[1];
    double pabz = (*coord[i].pab)[2];
    double pbbx = (*coord[i].pbb)[0];
    double pbby = (*coord[i].pbb)[1];
    double pbbz = (*coord[i].pbb)[2];

    // Evaluate [Kgnl2] terms.
    
    double cx = dHdx*sxx + dHdy*sxy + dHdz*sxz;
    double cy = dHdx*sxy + dHdy*syy + dHdz*syz;
    double cz = dHdx*sxz + dHdy*syz + dHdz*szz;

    double Kaa2 = cx*paax + cy*paay + cz*paaz; 
    double Kab2 = cx*pabx + cy*paby + cz*pabz; 
    double Kbb2 = cx*pbbx + cy*pbby + cz*pbbz; 

    Kgnl[ia][ia] = Kaa2;
    Kgnl[ia][ib] = Kab2;
    Kgnl[ib][ia] = Kab2;
    Kgnl[ib][ib] = Kbb2;

    // Evaluate [Kgnl3] terms.
    
    cx = As1*dHdx + As2*dHdy + As3*dHdz;
    cy = As4*dHdx + As5*dHdy + As6*dHdz;
    cz = As7*dHdx + As8*dHdy + As9*dHdz;

    double Kaa3 = cx*paax + cy*paay + cz*paaz; 
    double Kab3 = cx*pabx + cy*paby + cz*pabz; 
    double Kbb3 = cx*pbbx + cy*pbby + cz*pbbz; 

    Kgnl[ia][ia] += Kaa3;
    Kgnl[ia][ib] += Kab3;
    Kgnl[ib][ia] += Kab3;
    Kgnl[ib][ib] += Kbb3;
  }
}

// ======================================================= End of file =====
