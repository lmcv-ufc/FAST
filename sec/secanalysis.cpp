// -------------------------------------------------------------------------
// secanalysis.cpp - implementation of the Section Analysis class.
// -------------------------------------------------------------------------
// Created:      08-Mar-2013     Iuri Barcelos Rocha
//
// Modified:     20-Jul-2013     Iuri Barcelos Rocha
//               Implementation of Laminated Solid section.
//
// Modified:     03-Aug-2014     Evandro Parente Junior
//               Implementation of Retangular Plane Frame section.
//
// Modified:     08-Feb-2015     Evandro Parente Junior
//               Plane frame integration section using Gauss method.
//
// Modified:     12-Sep-2018     Samir Parente Auad
//               Implementation of General Shell section.
//
// Modified:     06-Nov-2019     Evandro Parente Junior
//               Implementation of MbMatrix methods.
//
// Modified:     16-Sep-2022     Pedro Ygor Rodrigues Mesquita
//               Implementation of HSDT plate section.
//
// Modified:     10-Mar-2023     Evandro Parente Junior
//               Refactoration of Stress and CMatrix methods.
// -------------------------------------------------------------------------

#include <math.h>
#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <iostream>

using namespace std;

#include "section.h"
#include "secbar.h"
#include "secanalysis.h"
#include "intpoint.h"
#include "element.h"
#include "elmparam.h"
#include "vec.h"
#include "mat.h"
#include "gblvar.h"
#include "ctrl.h"

const double Weight[3] = {0.33333333333333,1.33333333333333,0.33333333333333};
const double Point[3] = {-1.00000000000000,0.00000000000000,1.00000000000000};

// -------------------------------------------------------------------------
// Public methods:
//

// ================================= Create ================================

cSecAnalysis *cSecAnalysis :: Create(cElement *elm, cIntPoint *ipt)
{
//  cout << "cSecAnalysis :: Create\n";
  cSecAnalysis *an = 0;

  cSection *sec = elm->GetSection( );
  cAnModel *anm = elm->GetAnModel( );
  switch(sec->GetType( ))
  {
    case SEC_HOM_ISO_2D:
      if (anm->GetType( ) == THICK_PLATE)
        an = new cSecAnHomPlate(elm, ipt);
      else if (anm->GetType( ) == SHALLOW_SHELL || anm->GetType( ) == DONNELL_SHELL)
        an = new cSecAnHomShell(elm, ipt);
      else
        an = new cSecAnHomogeneous(elm, ipt);
    break;

    case SEC_HOM_ISO_3D:
      if (anm->GetType( ) == SHELL)
        an = new cSecAnHomDegShell(elm, ipt);
      else
        an = new cSecAnHomogeneous(elm, ipt);
    break;

    case SEC_HOM_ORTHO_2D:
    case SEC_HOM_ORTHO_3D:
      cout << "Homogenous-orthotropic section not implemented yet!\n";
      exit(0);
    break;

//    case SEC_LAMINATED_2D:
//      if (anm->GetType( ) == THICK_PLATE)
//        an = new cSecAnLaminatedPlate(elm, ipt);
//      else
//        an = new cSecAnLaminatedMembr(elm, ipt);
//    break;

    case SEC_FGM_2D:
      an = new cSecAnFGM2D(elm, ipt);
    break;

    case SEC_FGM_3D:
      an = new cSecAnFGM3D(elm, ipt);
    break;

    case SEC_FGM_SHELL:
      an = new cSecAnFGMShell(elm, ipt);
    break;

    case SEC_HSDT_PLATE:
    an = new cSecAnHSDTPlate(elm, ipt);
     break;

    case SEC_HSDT_FGM_PLATE:
    an = new cSecAnHSDTFGMPlate(elm, ipt);
    break;

    case SEC_GENERAL_SHELL:
      an = new cSecAnGenShell(elm, ipt);
    break;

    case SEC_LAMINATED_SHELL:
      an = new cSecAnLaminatedShell(elm, ipt);
    break;

    case SEC_INTERFACE_2D:
    case SEC_INTERFACE_3D:
      an = new cSecAnHomogeneous(elm, ipt);
    break;

    case SEC_BAR_HOM_RECT:
      if (sec->IntegrationMethod( ) == SEC_GAUSS)
        an = new cSecAnPlFrameRectGauss(elm, ipt);
      else
        an = new cSecAnPlFrameRectFiber(elm, ipt);
    break;

    case SEC_BAR_RC_RECT:
      an = new cSecAnPlFrameRCRect(elm, ipt);
    break;

    default:
      cout << "Section analysis not implemented yet!\n";
      exit(0);
    break;
  }
//  cout << "End of cSecAnalysis :: Create\n";

  return(an);
}

// ============================= cSecAnalysis ==============================

cSecAnalysis :: cSecAnalysis(cElement *elm, cIntPoint *ipt)
{
  Elm = elm;
  Sec = Elm->GetSection( );
}

// ============================ ~cSecAnalysis ==============================

cSecAnalysis :: ~cSecAnalysis(void)
{
}

// =============================== MbMatrix ================================

void cSecAnalysis :: MbMatrix(cMatrix &Mb)
{
  cout << "Error: MbMatrix not implemented for this section!\n";
  exit(0);
}

// ================================ Stress =================================

int cSecAnalysis :: Stress(double tref, cVector &temp, cVector &eps, cVector &sig)
{
  cout << "Error: Stress with temperature not implemented for this section!\n";
  exit(0);
}


// -------------------------------------------------------------------------
// Class cSecHomogeneous:
// -------------------------------------------------------------------------

// =========================== cSecHomogeneous =============================

cSecAnHomogeneous :: cSecAnHomogeneous(cElement *elm, cIntPoint *ipt) :
                     cSecAnalysis(elm, ipt)
{
  cMaterial *mat = Sec->GetMaterial( );
  Cmod = cConstModel::Create(elm, mat);
}

// =========================== ~cSecHomogeneous ============================

cSecAnHomogeneous :: ~cSecAnHomogeneous(void)
{
  delete Cmod;
}

// ================================ Stress =================================

int cSecAnHomogeneous :: Stress(cVector &eps, cVector &sig)
{
  // Get the stresses from the associated Constitutive Model.

  return Cmod->Stress(eps, sig);
}

// ================================ Stress =================================

int cSecAnHomogeneous :: Stress(double tref, cVector &temp, cVector &eps, cVector &sig)
{
  // Get the stresses from the associated Constitutive Model.

  return Cmod->Stress(tref, temp, eps, sig);
}

// =============================== CMatrix =================================

void cSecAnHomogeneous :: CMatrix(cMatrix &C)
{
  // Get the constitutive matrix from the associated Constitutive Model.

  Cmod->TangMat(C);
}

// =============================== CondMatrix ==============================

void cSecAnHomogeneous :: CondMatrix(cMatrix &D)
{
  // Get the conductive matrix from the associated Constitutive Model.

  Cmod->CondMatrix(D);
}

// ============================= UpdateState ===============================

void cSecAnHomogeneous :: UpdateState(void)
{
  // Update the state of the associated Constitutive Model.

  Cmod->UpdateState( );
}

// =============================== MbMatrix ================================

void cSecAnHomogeneous :: MbMatrix(cMatrix &Mb)
{
  cMaterial *mat = Sec->GetMaterial( );
  double rho = mat->GetDensity( );
  Mb.Zero( );
  int n = Mb.NRow( );
  for (int i = 0; i < n; i++) Mb[i][i] = rho;
}


// -------------------------------------------------------------------------
// Class cSecHomPlate:
// -------------------------------------------------------------------------

// ============================= cSecHomPlate ==============================

cSecAnHomPlate :: cSecAnHomPlate(cElement *elm, cIntPoint *ipt) :
                  cSecAnalysis(elm, ipt)
{
  // Check material compatibility and create the associated constitutive model.	

  cMaterial *mat = Sec->GetMaterial( );
  if (mat->GetMecProp( )->GetType( ) != MAT_ELASTIC_ISOTROPIC)
  {
    cout << "Error: chosen material is not compatible with Homogeneous Plate Section!\n";
    exit(0);
  }	  
  Cmod = cConstModel::Create(elm, mat);
}

// ============================= ~cSecHomPlate =============================

cSecAnHomPlate :: ~cSecAnHomPlate(void)
{
  delete Cmod;
}

// ================================ Stress =================================

int cSecAnHomPlate :: Stress(cVector &eps, cVector &sig)
{
  // Evaluate the generalized stresses assuming an elastic material.

  cMatrix C(5, 5);
  CMatrix(C);
  sig = C*eps;

  return(1);
}

// =============================== CMatrix =================================

void cSecAnHomPlate :: CMatrix(cMatrix &C)
{
  // Get the constitutive matrix from the associated Constitutive Model.

  Cmod->TangMat(C);

  // "Integrate" the elastic constitutive matrix in the plate thickness.

  double h = Sec->GetThickness( );
  double d = h*h*h/12.0;

  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) C[i][j] *= d;    // Bending

  double ks = 5.0/6.0;
  for (int i = 3; i < 5; i++)
    for (int j = 3; j < 5; j++) C[i][j] *= ks*h; // Shear
}

// ============================= UpdateState ===============================

void cSecAnHomPlate :: UpdateState(void)
{
  // Update the state of the associated Constitutive Model.

  Cmod->UpdateState( );
}

// =============================== MbMatrix ================================

void cSecAnHomPlate :: MbMatrix(cMatrix &Mb)
{
  cMaterial *mat = Sec->GetMaterial( );
  double rho = mat->GetDensity( );
  double h = Sec->GetThickness( );
  double I0 = rho*h;
  double I2 = rho*h*h*h/12.0;
  Mb.Zero( );
  Mb[0][0] = I0;
  Mb[1][1] = Mb[2][2] = I2;
}


// -------------------------------------------------------------------------
// Class cSecHomShell:
// -------------------------------------------------------------------------

// ============================= cSecHomShell ==============================

cSecAnHomShell :: cSecAnHomShell(cElement *elm, cIntPoint *ipt) :
                  cSecAnalysis(elm, ipt)
{
  // Check material compatibility and create the associated constitutive model.	

  cMaterial *mat = Sec->GetMaterial( );
  if (mat->GetMecProp( )->GetType( ) != MAT_ELASTIC_ISOTROPIC)
  {
    cout << "Error: chosen material is not compatible with Homogeneous Shell Section!\n";
    exit(0);
  }	  
  Cmod = cConstModel::Create(elm, mat);
}

// ============================= ~cSecHomShell =============================

cSecAnHomShell :: ~cSecAnHomShell(void)
{
  delete Cmod;
}

// ================================ Stress =================================

int cSecAnHomShell :: Stress(cVector &eps, cVector &sig)
{
  // Evaluate the generalized stresses assuming an elastic material.

  cMatrix C(8, 8);
  CMatrix(C);
  sig = C*eps;

  return(1);
}

// ================================ Stress =================================

int cSecAnHomShell :: Stress(double tref, cVector &temp, cVector &eps, cVector &sig)
{
  // Evaluate the section elastic matrix.

  cMatrix C(8, 8);
  CMatrix(C);

  // Get thermal expansion vector.

  cVector alpha(5);
  Cmod->AlphaVec(alpha);       // Evandro: por que? Não deveria ser do AnModel?

  // Evaluate temperature variation.

  double h  = Sec->GetThickness( );
  double Tb = temp[0] - tref;  // Bottom variation
  double Tt = temp[1] - tref;  // Top variation
  double Tm = (Tb + Tt)/2.0;   // Average variation
  double gT = (Tb - Tt)/h;     // Thickness gradient

  // Evaluate the generalized thermal strains.

  cVector eth(8);
  eth.Zero( );
  eth[0] = alpha[0]*Tm;        // Mebrane strain x
  eth[1] = alpha[1]*Tm;        // Mebrane strain y
  eth[3] = alpha[0]*gT;        // Curvature x
  eth[4] = alpha[1]*gT;        // Curvature y

  // Evaluate the section generalized stresses.

  cVector eff(8);
  eff = eps - eth;             // Effective generalized strains
  sig = C*eff;

  return(1);
}

// =============================== CMatrix =================================

void cSecAnHomShell :: CMatrix(cMatrix &C)
{
  // Get the constitutive matrix from the associated Constitutive Model.

  cMatrix Q(5, 5);
  Cmod->TangMat(Q);

  // "Integrate" the elastic constitutive matrix in the shell thickness.

  C.Zero( );
  double h = Sec->GetThickness( );
  double d = h*h*h/12.0;
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
    {
      C[i][j]     = Q[i][j]*h;     // Membrane
      C[i+3][j+3] = Q[i][j]*d;     // Bending
    }

  double ks = 5.0/6.0;
  for (int i = 3; i < 5; i++)
    for (int j = 3; j < 5; j++)
      C[i+3][j+3] = Q[i][j]*ks*h;  // Shear
}

// ============================= UpdateState ===============================

void cSecAnHomShell :: UpdateState(void)
{
  // Update the state of the associated Constitutive Model.

  Cmod->UpdateState( );
}

// =============================== MbMatrix ================================

void cSecAnHomShell :: MbMatrix(cMatrix &Mb)
{
  cMaterial *mat = Sec->GetMaterial( );
  double rho = mat->GetDensity( );
  double h  = Sec->GetThickness( );
  double I0 = rho*h;
  double I2 = rho*h*h*h/12.0;
  Mb.Zero( );
  Mb[0][0] = Mb[1][1] = Mb[2][2] = I0;
  Mb[3][3] = Mb[4][4] = I2;
}


// -------------------------------------------------------------------------
// Class cSecHomDegShell:
// -------------------------------------------------------------------------

// ============================= cSecHomDegShell ===========================

cSecAnHomDegShell :: cSecAnHomDegShell(cElement *elm, cIntPoint *ipt) :
                  cSecAnalysis(elm, ipt)
{
  // Check material compatibility and create the associated constitutive model.	

  cMaterial *mat = Sec->GetMaterial( );
  if (mat->GetMecProp( )->GetType( ) != MAT_ELASTIC_ISOTROPIC)
  {
    cout << "Error: chosen material is not compatible with Homogeneous DegShell Section!\n";
    exit(0);
  }	  
  Cmod = cConstModel::Create(elm, mat);

  // Store given integration points.

  Ipt = ipt;
}

// ============================= ~cSecHomDegShell ==========================

cSecAnHomDegShell :: ~cSecAnHomDegShell(void)
{
  delete Cmod;
}

// ================================ Stress =================================

int cSecAnHomDegShell :: Stress(cVector &eps, cVector &sig)
{
  // Get Nodal coordinates.
  int nn            = Elm->GetShape( )->GetNumNode( );
  sNodeCoord *coord = new sNodeCoord[nn];
  Elm->GetShape( )->NodalCoord(coord);

  // Compute local system.
  cMatrix R(3,3);
  Elm->GetShape( )->LocalSys(Ipt->GetCoord( ),coord,R);

  // Get element local vectors.
  cVector v1(3),v2(3),v3(3);
  for(int i = 0; i < 3; ++i)
  {
    v1(i) = R[i][0];
    v2(i) = R[i][1];
    v3(i) = R[i][2];
  }

  // Evaluate the transformation matrix.
  cMatrix T(6,6);
  Elm->GetAnModel( )->GetMecMod( )->TMatrix(v1,v2,v3,T);

  // Compute the local strains.
  cVector epsl(6);
  epsl = T*eps;

  // Evaluate the local stresses using the associated Constitutive Model.
  cVector sigl(6);
  Cmod->Stress(epsl, sigl);

  // Return the local stresses.
  if (Printing)  // Local stresses
    sig = sigl;
  else
    sig = t(T)*sigl;

  // Release memory.
  delete []coord;

  // Return success.
  return(1);
}

// =============================== CMatrix =================================

void cSecAnHomDegShell :: CMatrix(cMatrix &C)
{
  // Get Nodal coordinates.
  int nn            = Elm->GetShape( )->GetNumNode( );
  sNodeCoord *coord = new sNodeCoord[nn];
  Elm->GetShape( )->NodalCoord(coord);

  // Compute local system.
  cMatrix R(3,3);
  Elm->GetShape( )->LocalSys(Ipt->GetCoord( ),coord,R);

  // Get element local vectors.
  cVector v1(3),v2(3),v3(3);
  for(int i = 0; i < 3; ++i)
  {
    v1(i) = R[i][0];
    v2(i) = R[i][1];
    v3(i) = R[i][2];
  }

  // Evaluate the transformation matrix.
  cMatrix T(6,6);
  Elm->GetAnModel( )->GetMecMod( )->TMatrix(v1,v2,v3,T);

  // Get the local constitutive matrix from the associated Constitutive Model.
  cMatrix lC(6, 6);
  Cmod->TangMat(lC);

  // Evaluate the global constitutive matrix.
  C.Zero( );
  MatTripBtCB(T, lC, 1.0, C);

  // Release memory.
  delete []coord;
}

// ============================= UpdateState ===============================

void cSecAnHomDegShell :: UpdateState(void)
{
  // Update the state of the associated Constitutive Model.

  Cmod->UpdateState( );
}

// =============================== MbMatrix ================================

void cSecAnHomDegShell :: MbMatrix(cMatrix &Mb)
{
  cMaterial *mat = Sec->GetMaterial( );
  double rho = mat->GetDensity( );
  Mb.Zero( );
  int n = Mb.NRow( );
  for (int i = 0; i < n; i++) Mb[i][i] = rho;
}


// -------------------------------------------------------------------------
// Class cSecAnFGM2D:
// -------------------------------------------------------------------------

// ============================== cSecAnFGM2D ==============================

cSecAnFGM2D :: cSecAnFGM2D(cElement *elm, cIntPoint *ipt) :
               cSecAnalysis(elm, ipt)
{
  // Evaluate the shape functions at the integration point.

  cShape *shp = elm->GetShape( );
  int nnode = shp->GetNumNode( );
  double *shpfnc = new double[nnode];
  sNatCoord p = ipt->GetCoord( );
  shp->ShpFunc(p, shpfnc);

  // Evaluate the volume fraction at the integration point.

  cSecFGM2D *sec = (cSecFGM2D *)Sec;
  double V2 = 0.0;
  for (int i = 0; i < nnode; i++)
  {
    // Search element node in the section node map.	  

    cNode *node = shp->GetNode(i);
    VFNodeMap::iterator it = sec->SecVF.find(node->GetLabel( ));
    if (it == sec->SecVF.end())  // Node not in the section map
    {
      cout << "Error: node " << node->GetLabel( ) << " not in Section " << sec->Label << "!\n";
      exit(0);
    }	  

    // Get the nodal volume fraction and interpolate to the integration point.

    double V2i = it->second;  // Nodal volume fraction
//    cout << "node = " << node->GetLabel( ) << "  v2 = " << V2i << "\n";
    V2 += shpfnc[i]*V2i;
  }
//  cout << "V2 = " << V2 << "\n";

  // Code for comparison with examples where the material gradation is defined 
  // by an analytical expression in Cartesian coordinates (x, y).

#if 0  
  // Get the nodal coordinates.

  sNodeCoord *coord = new sNodeCoord[nnode];
  shp->NodalCoord(coord);

  // Evaluate Cartesian coordinates of the integration point.

  double *mapfnc = new double[nnode];
  shp->MapFunc(p, mapfnc);
  double x = 0.0;
  double y = 0.0;
  for (int i = 0; i < nnode; i++)
  {
    x += mapfnc[i]*coord[i].x;
    y += mapfnc[i]*coord[i].y;
  }
  cout << "x = " << x << "  ";
  cout << "y = " << y << "  ";
  delete []coord;
  delete []mapfnc;

  // Evaluate the exact gradation using the problem expression.

//  double V2ex = x/2.0;
  double V2ex = exp((2.0 - y)*log(4.0)/2.0)/4.0; // CILAMCE 2022
  double dif = V2/V2ex - 1.0;
  cout << "V2ex = " <<  V2ex << "  dif = " << 100*dif << " %\n";
#endif

  // Create the associated constitutive model.

  double param[1];
  cMaterial *mat = Sec->GetMaterial( );
  param[0] = V2;
  Cmod = cConstModel::Create(elm, mat, param);

  // Release memory.

  delete []shpfnc;
}

// ============================= ~cSecAnFGM2D ==============================

cSecAnFGM2D :: ~cSecAnFGM2D(void)
{
  delete Cmod;
}

// ================================ Stress =================================

int cSecAnFGM2D :: Stress(cVector &eps, cVector &sig)
{
  // Get the stresses from the associated Constitutive Model.

  return Cmod->Stress(eps, sig);
}

// =============================== CMatrix =================================

void cSecAnFGM2D :: CMatrix(cMatrix &C)
{
  // Get the constitutive matrix from the associated Constitutive Model.

  Cmod->TangMat(C);
}

// ============================= UpdateState ===============================

void cSecAnFGM2D :: UpdateState(void)
{
  // Update the state of the associated Constitutive Model.

  Cmod->UpdateState( );
}


// -------------------------------------------------------------------------
// Class cSecAnFGM3D:
// -------------------------------------------------------------------------

// ============================== cSecAnFGM3D ==============================

cSecAnFGM3D :: cSecAnFGM3D(cElement *elm, cIntPoint *ipt) :
               cSecAnalysis(elm, ipt)
{
  // Evaluate the shape functions at the integration point.

  cShape *shp = elm->GetShape( );
  int nnode = shp->GetNumNode( );
  double *shpfnc = new double[nnode];
  sNatCoord p = ipt->GetCoord( );
  shp->ShpFunc(p, shpfnc);

  // Evaluate the volume fraction at the integration point.

  cSecFGM3D *sec = (cSecFGM3D *)Sec;
  double V2 = 0.0;
  for (int i = 0; i < nnode; i++)
  {
    // Search element node in the section node map.	  

    cNode *node = shp->GetNode(i);
    VFNodeMap::iterator it = sec->SecVF.find(node->GetLabel( ));
    if (it == sec->SecVF.end())  // Node not in the section map
    {
      cout << "Error: node " << node->GetLabel( ) << " not in Section " << sec->Label << "!\n";
      exit(0);
    }	  

    // Get the nodal volume fraction and interpolate to the integration point.

    double V2i = it->second;  // Nodal volume fraction
//    cout << "node = " << node->GetLabel( ) << "  v2 = " << V2i << "\n";
    V2 += shpfnc[i]*V2i;
  }
//  cout << "V2 = " << V2 << "\n";

  // Code for comparison with examples where the material gradation is defined 
  // by an analytical expression in Cartesian coordinates (x, y).

#if 0  
  // Get the nodal coordinates.

  sNodeCoord *coord = new sNodeCoord[nnode];
  shp->NodalCoord(coord);

  // Evaluate Cartesian coordinates of the integration point.

  double *mapfnc = new double[nnode];
  shp->MapFunc(p, mapfnc);
  double x = 0.0;
  double y = 0.0;
  double z = 0.0;
  for (int i = 0; i < nnode; i++)
  {
    x += mapfnc[i]*coord[i].x;
    y += mapfnc[i]*coord[i].y;
    z += mapfnc[i]*coord[i].z;
  }
  cout << "x = " << x << "  ";
  cout << "y = " << y << "  ";
  cout << "z = " << z << "  ";
  delete []coord;
  delete []mapfnc;

  // Evaluate the exact gradation using the problem expression.

  double V2ex = 1.0 - z/10.0; // CILAMCE 2022
  double dif = V2/V2ex - 1.0;
  cout << "V2ex = " <<  V2ex << "  dif = " << 100*dif << " %\n";
#endif

  // Create the associated constitutive model.

  double param[1];
  cMaterial *mat = Sec->GetMaterial( );
  param[0] = V2;
  Cmod = cConstModel::Create(elm, mat, param);

  // Release memory.

  delete []shpfnc;
}

// ============================= ~cSecAnFGM3D ==============================

cSecAnFGM3D :: ~cSecAnFGM3D(void)
{
}

// ================================ Stress =================================

int cSecAnFGM3D :: Stress(cVector &eps, cVector &sig)
{
  // Get the stresses from the associated Constitutive Model.

  return Cmod->Stress(eps, sig);
}

// =============================== CMatrix =================================

void cSecAnFGM3D :: CMatrix(cMatrix &C)
{
  // Get the constitutive matrix from the associated Constitutive Model.

  Cmod->TangMat(C);
}

// ============================= UpdateState ===============================

void cSecAnFGM3D :: UpdateState(void)
{
  // Update the state of the associated Constitutive Model.

  Cmod->UpdateState( );
}


// -------------------------------------------------------------------------
// Class cSecAnFGMShell:
// -------------------------------------------------------------------------

// ============================ cSecAnFGMShell =============================

cSecAnFGMShell :: cSecAnFGMShell(cElement *elm, cIntPoint *ipt) :
                  cSecAnalysis(elm, ipt)
{
  // Get section data.

  //cout << "cSecAnFGMShell :: cSecAnFGMShell \n";
  Anm = elm->GetAnModel( );
  cSecFGMShell *sec = (cSecFGMShell *)Sec;
  double expn = sec->ExpN;          // Exponent of the volume fraction
  double h = Sec->GetThickness( );  // Thickness
  cMaterial *mat = Sec->GetMaterial( );
  //cout << "expn  = " << expn  << "\n";

  // Create the section integration points and associated constitutive models.

  NumSecPnt = Sec->GetNumSecPnt( );
  IntSecPnt = cIntPoint :: CreateLinePoints(NumSecPnt, GAUSS, &NumSecPnt);
  Cmod = new cConstModel*[NumSecPnt];
  double zb = -h/2.0;
  double zt = zb + h;
  double param[1];
  for (int i = 0; i < NumSecPnt; i++)
  {
    double z = IntSecPnt[i].GetCoord( ).r*(zt - zb)/2.0;
    double a = 0.5 + z/h;
    double V2 = pow(a, expn);
//    cout << "z  = " << z  << "  ";
//    cout << "V2 = " << V2 << "\n";
    param[0] = V2;
    Cmod[i] = cConstModel::Create(elm, mat, param);
  }
  //cout << "end of cSecAnFGMShell :: cSecAnFGMShell \n";
}

// =========================== ~cSecAnFGMShell ============================

cSecAnFGMShell :: ~cSecAnFGMShell(void)
{
//  cout << "cSecAnFGMShell :: ~cSecAnFGMShell\n";
  for (int i = 0; i < NumSecPnt; i++)
    delete Cmod[i];

  delete []Cmod;
//  cout << "end of cSecAnFGMShell :: ~cSecAnFGMShell\n";
}

// ================================ Stress =================================

int cSecAnFGMShell :: Stress(cVector &eps, cVector &sig)
{
  //cout << "cSecAnFGMShell :: Stress\n";

  // Evaluate the generalized stress in the section through numeric integration.

  sig.Zero( );
  cVector epsg(5);    // Strain at a integration point
  cVector sigg(5);    // Stress at a integration point
  double ks = 5.0/6.0;
  double h  = Sec->GetThickness( );
  for (int i = 0; i < NumSecPnt; i++)
  {
    // Integration point position.

    double z = IntSecPnt[i].GetCoord( ).r*h/2.0;

    // Integration point strains.

    epsg[0] = eps[0] + z*eps[3];
    epsg[1] = eps[1] + z*eps[4];
    epsg[2] = eps[2] + z*eps[5];
    epsg[3] = eps[6];
    epsg[4] = eps[7];
    //cout << "{epsg}\n";
    //epsg.Print( );

    // Integration point stresses (material).

    if (!Cmod[i]->Stress(epsg, sigg))
      return(0);
    //cout << "{sigg}\n";
    //sigg.Print( );

    // Generalized stresses (section).

    double coeff = IntSecPnt[i].GetWeight( )*h/2.0;
    sig[0] +=    coeff*sigg[0];   // Nx
    sig[1] +=    coeff*sigg[1];   // Ny
    sig[2] +=    coeff*sigg[2];   // Nxy
    sig[3] +=  z*coeff*sigg[0];   // Mx
    sig[4] +=  z*coeff*sigg[1];   // My
    sig[5] +=  z*coeff*sigg[2];   // Mxy
    sig[6] += ks*coeff*sigg[3];   // Vx
    sig[7] += ks*coeff*sigg[4];   // Vy
  }
  //cout << "{sig}\n";
  //sig.Print( );

  //cout << "End of cSecAnFGMShell :: Stress\n";
  return(1);
}

// ================================ Stress =================================

int cSecAnFGMShell :: Stress(double tref, cVector &temp, cVector &eps, cVector &sig)
{
  //cout << "cSecAnFGMShell :: Stress (temperature)\n";

  // Evaluate the generalized stress in the section through numeric integration.

  sig.Zero( );
  cVector epsg(5);    // Strain at a integration point
  cVector sigg(5);    // Stress at a integration point
  cVector Alpha(5);   // Alpha vector at a integration point
  double ks = 5.0/6.0;
  double h = Sec->GetThickness( );
  double a = (temp[1] - temp[0])/h;
  double b = (temp[1] + temp[0])/2.0;
  for (int i = 0; i < NumSecPnt; i++)
  {
    // Integration point position.

    double z = IntSecPnt[i].GetCoord( ).r*h/2.0;

    // Integration point temperature (linear variation).

    double Tz = a*z + b;
    double dT = Tz - tref;

    // Integration point strains.

    Cmod[i]->AlphaVec(Alpha);
    epsg[0] = eps[0] + z*eps[3] - Alpha[0]*dT;
    epsg[1] = eps[1] + z*eps[4] - Alpha[1]*dT;
    epsg[2] = eps[2] + z*eps[5] - Alpha[2]*dT;
    epsg[3] = eps[6];
    epsg[4] = eps[7];

    // Integration point stresses (material).

    if (!Cmod[i]->Stress(epsg, sigg))
        return(0);

    // Generalized stresses (section).

    double coeff = IntSecPnt[i].GetWeight( )*h/2.0;
    sig[0] +=    coeff*sigg[0];   // Nx
    sig[1] +=    coeff*sigg[1];   // Ny
    sig[2] +=    coeff*sigg[2];   // Nxy
    sig[3] +=  z*coeff*sigg[0];   // Mx
    sig[4] +=  z*coeff*sigg[1];   // My
    sig[5] +=  z*coeff*sigg[2];   // Mxy
    sig[6] += ks*coeff*sigg[3];   // Vx
    sig[7] += ks*coeff*sigg[4];   // Vy
  }
  //cout << "End of cSecAnFGMShell :: Stress (temperature)\n";
  return(1);
}

// =============================== CMatrix =================================

void cSecAnFGMShell :: CMatrix(cMatrix &C)
{
  //cout << "cSecAnFGMShell :: CMatrix\n";
  // Get section data.

  int nsc = Anm->GetMecMod( )->GetDimQMatrix( );  // Number of strains
  int ngs = Anm->GetMecMod( )->GetDimBMatrix( );  // Number of gen. strains  
//  cout << "nsc = " << nsc << "  ";
//  cout << "ngs = " << ngs << "\n";
  double h = Sec->GetThickness( );

  // Evaluate the constitutive matrix through numeric integration.

  C.Zero( );
  cMatrix Q(nsc, nsc);
  cMatrix Cz(ngs, ngs);
  for (int i = 0; i < NumSecPnt; i++)
  {
    // Integration point position.

    double z = IntSecPnt[i].GetCoord( ).r*h/2.0;

    // Integration point constitutive matrix (material).

    Cmod[i]->TangMat(Q);

    // Integrand: [Q], [Q]*z, [Q]*z^2, [Qs]*ks.

    CalcCz(z, Q, Cz);

    // Section constitutive matrix.

    double coeff = IntSecPnt[i].GetWeight( )*h/2.0;
    C += coeff*Cz;
  }
  //cout << "[C]\n";
  //cout << scientific << setprecision(4);
  //C.Print( );
  //cout << "End of cSecAnFGMShell :: CMatrix\n";
}

// ============================= UpdateState ===============================

void cSecAnFGMShell :: UpdateState(void)
{
  // Update the state of each integration point.

  for (int i = 0; i < NumSecPnt; i++)
    Cmod[i]->UpdateState( );
}

// =============================== MbMatrix ================================

void cSecAnFGMShell :: MbMatrix(cMatrix &Mb)
{
//  cout << "cSecAnFGMShell :: MbMatrix\n";
  // Get section data.

  Anm = Elm->GetAnModel( );
  cSecFGMShell *sec = (cSecFGMShell *)Sec;
  double expn = sec->ExpN;          // Exponent of the volume fraction
  double h = Sec->GetThickness( );  // Thickness
  cMaterial *mat = Sec->GetMaterial( );
  cVector matrho(2);
  mat->GetDensity(matrho);
//  cout << "Des0 " << matrho[0] << endl;
//  cout << "Des1 " << matrho[1] << endl;
//  cout << "expn  = " << expn  << "\n";

  // Evaluate [Mb] matrix through integration.

  double I0 = 0.0;
  double I1 = 0.0;
  double I2 = 0.0;
  double zb = -h/2.0;
  double zt = zb + h;
  double J  = h/2.0;
  for (int i = 0; i < NumSecPnt; i++)
  {
    double w = IntSecPnt[i].GetWeight( );
    double z = IntSecPnt[i].GetCoord( ).r*(zt - zb)/2.0;
    double a = 0.5 + z/h;
    double V2 = pow(a, expn);
    double rho = matrho[0]*(1.0 - V2) + matrho[1]*V2;
    I0 += rho*w*J;
    I1 += rho*z*w*J;
    I2 += rho*z*z*w*J;
  }
  Mb.Zero( );
  Mb[0][0] = Mb[1][1] = Mb[2][2] = I0;
  Mb[3][3] = Mb[4][4] = I2;
  Mb[0][4] = Mb[4][0] = I1;
  Mb[1][3] = Mb[3][1] = -I1;
//  cout << "[Mb]\n";
//  cout << scientific << setprecision(4);
//  Mb.Print( );
//  cout << "End of cSecAnFGMShell :: MbMatrix\n";
}

// -------------------------------------------------------------------------
// Protected methods:
//

// =============================== CalcCz =================================

void cSecAnFGMShell :: CalcCz(double z, cMatrix &Q, cMatrix &Cz)
{
  Cz.Zero( );

  // Membrane and bending terms.

  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
    {
      Cz[i][j] = Q[i][j];                   // [A]
      Cz[3+i][j] = Cz[i][3+j] = z*Q[i][j];  // [B]
      Cz[3+i][3+j] = z*z*Q[i][j];           // [D]
    }

  // Transverse shear terms.

  double ks = 5.0/6.0;
  for (int i = 3; i < 5; i++)
    for (int j = 3; j < 5; j++)
      Cz[3+i][3+j] = ks*Q[i][j];            // [G]
}


// -------------------------------------------------------------------------
// Class cSecAnHSDTPlate:
// -------------------------------------------------------------------------

// ============================ cSecAnHSDTPlate =============================

cSecAnHSDTPlate :: cSecAnHSDTPlate(cElement *elm, cIntPoint *ipt) :
                  cSecAnalysis(elm, ipt)
{
  // Get section data.

  Anm = elm->GetAnModel( );
  cMaterial *mat = Sec->GetMaterial( );

  // Create the section integration points and associated constitutive models.

  NumSecPnt = 10;
  IntSecPnt = cIntPoint :: CreateLinePoints(NumSecPnt, GAUSS, &NumSecPnt);
  Cmod = new cConstModel*[NumSecPnt];
  for (int i = 0; i < NumSecPnt; i++)
  {
    Cmod[i] = cConstModel::Create(elm, mat);
  }
}

// =========================== ~cSecAnHSDTPlate ============================

cSecAnHSDTPlate :: ~cSecAnHSDTPlate(void)
{
  for (int i = 0; i < NumSecPnt; i++)
    delete Cmod[i];

  delete []Cmod;
}

// ================================ Stress =================================

int cSecAnHSDTPlate :: Stress(cVector &eps, cVector &sig)
{
  // Evaluate the generalized shell stresses assuming an elastic material.

  cMatrix C(11, 11);
  CMatrix(C);
  sig = C*eps;

  return(1);
}

// =============================== CMatrix =================================

void cSecAnHSDTPlate :: CMatrix(cMatrix &C)
{
  // Get section data.

  int nsc = Anm->GetMecMod( )->GetDimQMatrix( );  // Number of strains
  int ngs = Anm->GetMecMod( )->GetDimBMatrix( );  // Number of gen. strains  
  double h = Sec->GetThickness( );

  // Evaluate the constitutive matrix through numeric integration.

  C.Zero( );
  cMatrix Q(nsc, nsc);
  cMatrix Cz(ngs, ngs);
  for (int i = 0; i < NumSecPnt; i++)
  {
    // Integration point position.

    double z = IntSecPnt[i].GetCoord( ).r*h/2.0;

    // Integration point constitutive matrix (material).

    Cmod[i]->TangMat(Q);

    // Integrand: [Q], [Q]*z, [Q]*z^2, [Qs], ...

    CalcCz(h, z, Q, Cz);

    // Section constitutive matrix.

    double coeff = IntSecPnt[i].GetWeight( )*h/2.0;
    C += coeff*Cz;
  }
}

// ============================= UpdateState ===============================

void cSecAnHSDTPlate :: UpdateState(void)
{
  // Update the state of each integration point.

  for (int i = 0; i < NumSecPnt; i++)
    Cmod[i]->UpdateState( );
}

// =============================== MbMatrix ================================

void cSecAnHSDTPlate :: MbMatrix(cMatrix &Mb)
{
  // Get section data.

  Anm = Elm->GetAnModel( );
  cSecHSDTPlate *sec = (cSecHSDTPlate *)Sec;
  double h = Sec->GetThickness( );  // Thickness
  cMaterial *mat = Sec->GetMaterial( );
  cVector matrho(1);
  mat->GetDensity(matrho);
  double rho = matrho[0];

  // Evaluate [Mb] matrix through integration.

  double I0 = 0.0;
  double I1 = 0.0;
  double I2 = 0.0;
  double J1 = 0.0;
  double J2 = 0.0;
  double K1 = 0.0;
  double J  = h/2.0;
  for (int i = 0; i < NumSecPnt; i++)
  {
    double w = IntSecPnt[i].GetWeight( );
    double z = IntSecPnt[i].GetCoord( ).r*(h)/2.0;

    // Evaluate f(z)

    double fz = z;                 // First-Order Shear Deformation Theory
    if (sec->Theory == TSDT)       // Third-Order Shear Deformation Theory
    {
      fz  = z - 4.0*z*z*z/(3.0*h*h);
    }
    else if (sec->Theory == ESDT)  // Exponential Shear Deformation Theory
    {
      double alpha = exp(1);   // Karama et al. (2003) - ESDT1
      //double alpha = 2.85;   // Mantari el al. (2011) - ESDT2
      fz  = z*pow(alpha,-2*pow(z/h,2));
    }
    else if (sec->Theory == SSDT)  // Sinusoidal Shear Deformation Theory
    {
      fz  = (h/M_PI)*sin(M_PI*z/h);
    }
    else if (sec->Theory == SHSDT) // Soldatos Hyperbolic Shear Deformation Theory
    {
      fz  = h*sinh(z/h) - z*cosh(1.0/2.0);
    }

    I0 += rho*w*J;
    I1 += rho*z*w*J;
    I2 += rho*z*z*w*J;
    J1 += rho*fz*w*J;
    J2 += rho*fz*fz*w*J;
    K1 += rho*z*fz*w*J;
  }
  Mb.Zero( );
  Mb[0][0] = Mb[1][1] = Mb[2][2] = I0;
  Mb[0][3] = Mb[3][0] = Mb[1][4] = Mb[4][1] = -J1;
  Mb[3][3] = Mb[4][4] = J2;
  Mb[0][5] = Mb[5][0] = Mb[1][6] = Mb[6][1] = J1-I1;
  Mb[3][5] = Mb[5][3] = Mb[4][6] = Mb[6][4] = K1-J2;
  Mb[5][5] = Mb[6][6] = J2-2.0*K1+I2;
}

// -------------------------------------------------------------------------
// Protected methods:
//

// =============================== CalcCz =================================

void cSecAnHSDTPlate :: CalcCz(double h, double z, cMatrix &Q, cMatrix &Cz)
{
  Cz.Zero( );

  // Evaluate f(z), f'(z) and ks.

  cSecHSDTPlate *sec = (cSecHSDTPlate *)Sec;
  double fz  = z;                // First-Order Shear Deformation Theory
  double dfz = 1.0;
  double ks  = 5.0/6.0;
  if (sec->Theory == TSDT)       // Third-Order Shear Deformation Theory 
  {
    fz  = z - 4.0*z*z*z/(3.0*h*h);
    dfz = 1.0 - 4.0*z*z/(h*h);
    ks  = 1.0;
  }
  else if (sec->Theory == ESDT)  // Exponential Shear Deformation Theory
  {
    double alpha = exp(1);   // Karama et al. (2003) - ESDT1
    //double alpha = 2.85;   // Mantari el al. (2011) - ESDT2
    fz  = z*pow(alpha,-2*pow(z/h,2));
    dfz = (1-4*z*z/(h*h))*pow(alpha,-2*pow(z/h,2));
    ks  = 1.0;
  }
  else if (sec->Theory == SSDT)  // Sinusoidal Shear Deformation Theory
  {
    fz  = (h/M_PI)*sin(M_PI*z/h);
    dfz = cos(M_PI*z/h);
    ks  = 1.0;
  }
  else if (sec->Theory == SHSDT) // Soldatos Hyperbolic Shear Deformation Theory
  {
    fz  = h*sinh(z/h) - z*cosh(1.0/2.0);
    dfz = cosh(z/h) - cosh(1.0/2.0);
    ks  = 1.0;
  }

  // Membrane and bending terms.

  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
    {
      Cz[i][j] = Q[i][j];                               // [A]
      Cz[3+i][j] = Cz[i][3+j] = z*Q[i][j];              // [B]
      Cz[3+i][3+j] = z*z*Q[i][j];                       // [D]
      Cz[6+i][j] = Cz[i][6+j] = fz*Q[i][j];             // [E]
      Cz[3+i][6+j] = Cz[6+i][3+j] = z*fz*Q[i][j];       // [F]
      Cz[6+i][6+j] = fz*fz*Q[i][j];                     // [H]
    }

  // Transverse shear terms.

  Cz[9][9] = Cz[10][10] = ks*dfz*dfz*Q[4][4];           // [Ds]
}


// -------------------------------------------------------------------------
// Class cSecAnHSDTFGMPlate:
// -------------------------------------------------------------------------

// ============================ cSecAnHSDTFGMPlate =============================

cSecAnHSDTFGMPlate :: cSecAnHSDTFGMPlate(cElement *elm, cIntPoint *ipt) :
                  cSecAnalysis(elm, ipt)
{
  // Get section data.

  //cout << "cSecAnHSDTFGMPlate :: cSecAnHSDTFGMPlate \n";
  Anm = elm->GetAnModel( );
  cSecHSDTFGMPlate *sec = (cSecHSDTFGMPlate *)Sec;
  double expn = sec->ExpN;          // Exponent of the volume fraction
  double h = Sec->GetThickness( );  // Thickness
  cMaterial *mat = Sec->GetMaterial( );
  //cout << "expn  = " << expn  << "\n";

  // Create the section integration points and associated constitutive models.

  NumSecPnt = Sec->GetNumSecPnt( );
  IntSecPnt = cIntPoint :: CreateLinePoints(NumSecPnt, GAUSS, &NumSecPnt);
  Cmod = new cConstModel*[NumSecPnt];
  double zb = -h/2.0;
  double zt = zb + h;
  double param[1];
  for (int i = 0; i < NumSecPnt; i++)
  {
    double z = IntSecPnt[i].GetCoord( ).r*(zt - zb)/2.0;
    double a = 0.5 + z/h;
    double V2 = pow(a, expn);
    //cout << "z  = " << z  << "  ";
    //cout << "V2 = " << V2 << "\n";
    param[0] = V2;
    Cmod[i] = cConstModel::Create(elm, mat, param);
  }
  //cout << "end of cSecAnHSDTFGMPlate :: cSecAnHSDTFGMPlate \n";
}

// =========================== ~cSecAnHSDTFGMPlate ============================

cSecAnHSDTFGMPlate :: ~cSecAnHSDTFGMPlate(void)
{
  //cout << "cSecAnHSDTFGMPlate :: ~cSecAnHSDTFGMPlate\n";
  for (int i = 0; i < NumSecPnt; i++)
    delete Cmod[i];

  delete []Cmod;
  //cout << "end of cSecAnHSDTFGMPlate :: ~cSecAnHSDTFGMPlate\n";
}

// ================================ Stress =================================

int cSecAnHSDTFGMPlate :: Stress(cVector &eps, cVector &sig)
{
  // Evaluate the generalized shell stresses assuming an elastic material.

  cMatrix C(11, 11);
  CMatrix(C);
  sig = C*eps;

  return(1);
}

// =============================== CMatrix =================================

void cSecAnHSDTFGMPlate :: CMatrix(cMatrix &C)
{
  // Get section data.

  //cout << "cSecAnHSDTFGMPlate:: CMatrix\n";
  int nsc = Anm->GetMecMod( )->GetDimQMatrix( );  // Number of strains
  int ngs = Anm->GetMecMod( )->GetDimBMatrix( );  // Number of gen. strains  
  //cout << "nsc = " << nsc << "  ";
  //cout << "ngs = " << ngs << "\n";
  double h = Sec->GetThickness( );

  // Evaluate the constitutive matrix through numeric integration.

  C.Zero( );
  cMatrix Q(nsc, nsc);
  cMatrix Cz(ngs, ngs);
  for (int i = 0; i < NumSecPnt; i++)
  {
    // Integration point position.

    double z = IntSecPnt[i].GetCoord( ).r*h/2.0;

    // Integration point constitutive matrix (material).

    Cmod[i]->TangMat(Q);

    // Integrand: [Q], [Q]*z, [Q]*z^2, [Qs], ...

    CalcCz(h, z, Q, Cz);

    // Section constitutive matrix.

    double coeff = IntSecPnt[i].GetWeight( )*h/2.0;
    C += coeff*Cz;
  }
  //cout << "[C]\n";
  //cout << scientific << setprecision(4);
  //C.Print( );
  //cout << "End of cSecAnHSDTPlate :: CMatrix\n";
}

// ============================= UpdateState ===============================

void cSecAnHSDTFGMPlate :: UpdateState(void)
{
  // Update the state of each integration point.

  for (int i = 0; i < NumSecPnt; i++)
    Cmod[i]->UpdateState( );
}
 
// =============================== MbMatrix ================================

void  cSecAnHSDTFGMPlate :: MbMatrix(cMatrix &Mb)
{
  // Get section data.

  Anm = Elm->GetAnModel( );
  cSecHSDTFGMPlate *sec = (cSecHSDTFGMPlate *)Sec;
  double ExpN = sec->ExpN;             // Exponent of the volume fraction
  double h = Sec->GetThickness( );  // Thickness
  cMaterial *mat = Sec->GetMaterial( );
  cVector matrho(2);
  mat->GetDensity(matrho);

  // Evaluate [Mb] matrix through integration.

  double I0 = 0.0;
  double I1 = 0.0;
  double I2 = 0.0;
  double J1 = 0.0;
  double J2 = 0.0;
  double K1 = 0.0;
  double J  = h/2.0;
  for (int i = 0; i < NumSecPnt; i++)
  {
    double w = IntSecPnt[i].GetWeight( );
    double z = IntSecPnt[i].GetCoord( ).r*(h)/2.0;
    double a = 0.5 + z/h;
    double V2 = pow(a, ExpN);
    double rho = matrho[0]*(1.0 - V2) + matrho[1]*V2;

    // Evaluate f(z)

    double fz  = z;                // First-Order Shear Deformation Theory
    if (sec->Theory == TSDT)       // Third-Order Shear Deformation Theory
    {
      fz  = z - 4.0*z*z*z/(3.0*h*h);
    }
    else if (sec->Theory == ESDT)  // Exponential Shear Deformation Theory
    {
      double alpha = exp(1);   // Karama et al. (2003) - ESDT1
      //double alpha = 2.85;   // Mantari el al. (2011) - ESDT2
      fz  = z*pow(alpha,-2*pow(z/h,2));
    }
    else if (sec->Theory == SSDT)  // Sinusoidal Shear Deformation Theory
    {
      fz  = (h/M_PI)*sin(M_PI*z/h);
    }
    else if (sec->Theory == SHSDT) // Soldatos Hyperbolic Shear Deformation Theory
    {
      fz  = h*sinh(z/h) - z*cosh(1.0/2.0);
    }

    I0 += rho*w*J;
    I1 += rho*z*w*J;
    I2 += rho*z*z*w*J;
    J1 += rho*fz*w*J;
    J2 += rho*fz*fz*w*J;
    K1 += rho*z*fz*w*J;
  }
  Mb.Zero( );
  Mb[0][0] = Mb[1][1] = Mb[2][2] = I0;
  Mb[0][3] = Mb[3][0] = Mb[1][4] = Mb[4][1] = -J1;
  Mb[3][3] = Mb[4][4] = J2;
  Mb[0][5] = Mb[5][0] = Mb[1][6] = Mb[6][1] = J1-I1;
  Mb[3][5] = Mb[5][3] = Mb[4][6] = Mb[6][4] = K1-J2;
  Mb[5][5] = Mb[6][6] = J2-2.0*K1+I2;
}

// -------------------------------------------------------------------------
// Protected methods:
//

// =============================== CalcCz =================================

void cSecAnHSDTFGMPlate :: CalcCz(double h, double z, cMatrix &Q, cMatrix &Cz)
{
  Cz.Zero( );

  // Evaluate f(z), f'(z) and ks.

  cSecHSDTFGMPlate *sec = (cSecHSDTFGMPlate *)Sec;
  double fz  = z;                // First-Order Shear Deformation Theory
  double dfz = 1.0;
  double ks  = 5.0/6.0;
  if (sec->Theory == TSDT)       // Third-Order Shear Deformation Theory
  {
    fz  = z - 4.0*z*z*z/(3.0*h*h);
    dfz = 1.0 - 4.0*z*z/(h*h);
    ks  = 1.0;
  }
  else if (sec->Theory == ESDT)  // Exponential Shear Deformation Theory
  {
    double alpha = exp(1);   // Karama et al. (2003) - ESDT1
    //double alpha = 2.85;   // Mantari el al. (2011) - ESDT2
    fz  = z*pow(alpha,-2*pow(z/h,2));
    dfz = (1-4*z*z/(h*h))*pow(alpha,-2*pow(z/h,2));
    ks  = 1.0;
  }
  else if (sec->Theory == SSDT)  // Sinusoidal Shear Deformation Theory
  {
    fz  = (h/M_PI)*sin(M_PI*z/h);
    dfz = cos(M_PI*z/h);
    ks  = 1.0;
  }
  else if (sec->Theory == SHSDT) // Soldatos Hyperbolic Shear Deformation Theory
  {
    fz  = h*sinh(z/h) - z*cosh(1.0/2.0);
    dfz = cosh(z/h) - cosh(1.0/2.0);
    ks  = 1.0;
  }

  // Membrane and bending terms.

  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
    {
      Cz[i][j] = Q[i][j];                               // [A]
      Cz[3+i][j] = Cz[i][3+j] = z*Q[i][j];              // [B]
      Cz[3+i][3+j] = z*z*Q[i][j];                       // [D]
      Cz[6+i][j] = Cz[i][6+j] = fz*Q[i][j];             // [E]
      Cz[3+i][6+j] = Cz[6+i][3+j] = z*fz*Q[i][j];       // [F]
      Cz[6+i][6+j] = fz*fz*Q[i][j];                     // [H]
    }

  // Transverse shear terms.

  Cz[9][9] = Cz[10][10] = ks*dfz*dfz*Q[4][4];           // [Ds]
}


// -------------------------------------------------------------------------
// Class cSecGenShell:
// -------------------------------------------------------------------------

// ============================= cSecGenShell ==============================

cSecAnGenShell :: cSecAnGenShell(cElement *elm, cIntPoint *ipt) :
                  cSecAnalysis(elm, ipt)
{
}

// ============================= ~cSecGenShell =============================

cSecAnGenShell :: ~cSecAnGenShell(void)
{
}

// ================================ Stress =================================

int cSecAnGenShell :: Stress(cVector &eps, cVector &sig)
{
  // Evaluate the generalized shell stresses assuming an elastic material.

  cMatrix C(8, 8);
  CMatrix(C);
  sig = C*eps;

  return(1);
}

// =============================== CMatrix =================================

void cSecAnGenShell :: CMatrix(cMatrix &C)
{
  // Get the elastic constitutive matrix from the associated section.

  cSecGeneralShell *sec = (cSecGeneralShell *)Sec;
  C = sec->C;
//  cout << "[C]\n";
//  cout << scientific;
//  C.Print( );
}

// =============================== MbMatrix ================================

void cSecAnGenShell :: MbMatrix(cMatrix &Mb)
{
  // Get [Mb] mass matrix from the associated section.

  cSecGeneralShell *sec = (cSecGeneralShell *)Sec;
  Mb = sec->Mb;
//  cout << "[Mb]\n";
//  cout << scientific;
//  Mb.Print( );
}


// -------------------------------------------------------------------------
// Class cSecAnLaminatedShell:
// -------------------------------------------------------------------------

// ========================= cSecAnLaminatedShell ==========================

cSecAnLaminatedShell :: cSecAnLaminatedShell(cElement *elm, cIntPoint *ipt) :
                        cSecAnalysis(elm, ipt)
{
  Anm = elm->GetAnModel( );
  sLamina *Lamina = Sec->GetLayup( );
  int nlam = Sec->GetNumLam( );
  Cmod = new cConstModel*[3*nlam];

  for (int i = 0; i < nlam; i++)
    for (int j = 0; j < 3; j++)
      Cmod[3*i + j] = cConstModel::Create(elm, Lamina[i].Mat);
}

// ========================= ~cSecAnLaminatedShell =========================

cSecAnLaminatedShell :: ~cSecAnLaminatedShell(void)
{
  int nlam = Sec->GetNumLam( );
  for (int i = 0; i < (3*nlam); i++)
    delete Cmod[i];

  delete []Cmod;
}

// ================================ Stress =================================

int cSecAnLaminatedShell :: Stress(cVector &eps, cVector &sig)
{
  // Get section data.

  int nlam = Sec->GetNumLam( );
  sLamina *Lamina = Sec->GetLayup( );
  cVector a(3);
  Sec->GetOrientationVector(a);
  double beta = atan2(a[1], a[0]);

  // Evaluate the generalized stress in the section through numeric integration.

  sig.Zero( );
  cMatrix T(5,5);
  cVector epsg(5);
  cVector epsl(5);
  cVector sigl(5);
  cVector sigg(5);
  double ks = 5.0/6.0;
  double zb = -Sec->GetThickness( )/2.0;
  for (int i = 0; i < nlam; i++)
  {
    // Ply orientation.

    double theta = beta + Lamina[i].Ang*PI/180.0;
    Anm->GetMecMod( )->TMatrix(theta, T);

    // Integration in each ply.

    double zt = zb + Lamina[i].Thk;
    for (int j = 0; j < 3; j++)   // Evandro: deveria ser uma variável
    {
      // Integration point position.

      double z = (zb + zt)/2.0 + Point[j]*(zt - zb)/2.0;

      // Integration point strains in global system.

      epsg[0] = eps[0] + z*eps[3];
      epsg[1] = eps[1] + z*eps[4];
      epsg[2] = eps[2] + z*eps[5];
      epsg[3] = eps[6];
      epsg[4] = eps[7];

      // Integration point strains in local system.

      epsl = T*epsg;

      // Integration point stresses (material) in local system.

      if (!Cmod[3*i+j]->Stress(epsl, sigl))
        return(0);

      // Integration point stresses in global system.

      sigg = t(T)*sigl;

      // Generalized stresses (section).

      double coeff = Weight[j]*Lamina[i].Thk/2.0;
      sig[0] +=    coeff*sigg[0];   // Nx
      sig[1] +=    coeff*sigg[1];   // Ny
      sig[2] +=    coeff*sigg[2];   // Nxy
      sig[3] +=  z*coeff*sigg[0];   // Mx
      sig[4] +=  z*coeff*sigg[1];   // My
      sig[5] +=  z*coeff*sigg[2];   // Mxy
      sig[6] += ks*coeff*sigg[3];   // Vx
      sig[7] += ks*coeff*sigg[4];   // Vy
    }
    zb = zt;
  }
  return(1);
}

// ================================ Stress =================================

int cSecAnLaminatedShell :: Stress(double tref, cVector &temp, cVector &eps, cVector &sig)
{
  // Get section data.

  int nlam = Sec->GetNumLam( );
  sLamina *Lamina = Sec->GetLayup( );
  cVector a(3);
  Sec->GetOrientationVector(a);
  double beta = atan2(a[1], a[0]);
//  cout << "beta = " << beta << endl;

  // Evaluate temperature profile (linear variation).

  double h = Sec->GetThickness( );
  double c = (temp[1] - temp[0])/h;
  double d = (temp[1] + temp[0])/2.0;

  // Evaluate the generalized stress in the section through numeric integration.

  sig.Zero( );
  cMatrix T(5,5);
  cVector epsg(5);
  cVector epsl(5);
  cVector sigl(5);
  cVector sigg(5);
  cVector alpha(5);
  double ks = 5.0/6.0;
  double zb = -h/2.0;
  for (int i = 0; i < nlam; i++)
  {
    // Ply orientation.

    double theta = beta + Lamina[i].Ang*PI/180.0;
    Anm->GetMecMod( )->TMatrix(theta, T);

    // Integration in each ply.

    double zt = zb + Lamina[i].Thk;
    for (int j = 0; j < 3; j++)
    {
      // Integration point position.

      double z = (zb + zt)/2.0 + Point[j]*(zt - zb)/2.0;

      // Integration point strains in global system.

      epsg[0] = eps[0] + z*eps[3];
      epsg[1] = eps[1] + z*eps[4];
      epsg[2] = eps[2] + z*eps[5];
      epsg[3] = eps[6];
      epsg[4] = eps[7];

      // Integration point strains in local system.

      epsl = T*epsg;

      // Effective strains in local system.
  
      double Tz = c*z + d;
      double dT = Tz - tref;
      Cmod[3*i+j]->AlphaVec(alpha);  // Evandro: verificar
      epsl[0] -= alpha[0]*dT;
      epsl[1] -= alpha[1]*dT;
      epsl[2] -= alpha[2]*dT;

      // Integration point stresses (material) in local system.

      if (!Cmod[3*i+j]->Stress(epsl, sigl))
        return(0);

      // Integration point stresses in global system.

      sigg = t(T)*sigl;

      // Generalized stresses (section).

      double coeff = Weight[j]*Lamina[i].Thk/2.0;
      sig[0] +=    coeff*sigg[0];   // Nx
      sig[1] +=    coeff*sigg[1];   // Ny
      sig[2] +=    coeff*sigg[2];   // Nxy
      sig[3] +=  z*coeff*sigg[0];   // Mx
      sig[4] +=  z*coeff*sigg[1];   // My
      sig[5] +=  z*coeff*sigg[2];   // Mxy
      sig[6] += ks*coeff*sigg[3];   // Vx
      sig[7] += ks*coeff*sigg[4];   // Vy
    }
    zb = zt;
  }
  return(1);
}

// =============================== CMatrix =================================

void cSecAnLaminatedShell :: CMatrix(cMatrix &C)
{
//  cout << "End of cSecAnLaminatedShell :: CMatrix\n";

  // Get section data.

  int nlam = Sec->GetNumLam( );
  int nsc = Anm->GetMecMod( )->GetDimQMatrix( );  // Number of strains
  int ngs = Anm->GetMecMod( )->GetDimBMatrix( );  // Number of gen. strains
  sLamina *Lamina = Sec->GetLayup( );
  cVector a(3);
  Sec->GetOrientationVector(a);
  double beta = atan2(a[1], a[0]);
//  cout << "beta = " << beta << endl;

  // Evaluate the constitutive matrix through numeric integration.

  C.Zero( );
  cMatrix Q(nsc, nsc);
  cMatrix Qb(nsc, nsc);
  cMatrix T(nsc, nsc);
  cMatrix Cz(ngs, ngs);
  double zb = -Sec->GetThickness( )/2.0;
  for (int i = 0; i < nlam; i++)
  {
    // Ply orientation.

//    cout << "ply = " << i+1 << "  ";
//    cout << "theta = " << fixed << Lamina[i].Ang << endl;
    double theta = beta + Lamina[i].Ang*PI/180.0;
    Anm->GetMecMod( )->TMatrix(theta, T);

    // Integration in each ply.

    double zt = zb + Lamina[i].Thk;
    for (int j = 0; j < 3; j++)
    {
      // Integration point position.

      double z = (zb + zt)/2.0 + Point[j]*(zt - zb)/2.0;

      // Integration point constitutive matrix (material) in local system.

      Cmod[3*i+j]->TangMat(Q);

      // Integration point tangent matrix in global system.

      Qb.Zero( );
      MatTripBtCB(T, Q, 1.0, Qb);

      // Integrand: [Qb], [Qb]*z, [Qb]*z^2, [Qbs]*ks.

      CalcCz(z, Qb, Cz);

      // Section constitutive matrix.

      double coeff = Weight[j]*Lamina[i].Thk/2.0;
      C += coeff*Cz;
    }
    zb = zt;
  }
//  cout << "[C]\n";
//  cout << scientific << setprecision(4);
//  C.Print( );
//  cout << "End of cSecAnLaminatedShell :: CMatrix\n";
}

// ============================= UpdateState ===============================

void cSecAnLaminatedShell :: UpdateState(void)
{
  // Update the state of each integration point.

  int nlam = Sec->GetNumLam( );
  for (int i = 0; i < 3*nlam; i++)
    Cmod[i]->UpdateState( );
}

// =============================== MbMatrix ================================

void cSecAnLaminatedShell :: MbMatrix(cMatrix &Mb)
{
  // Get laminate data.

  int nlam = Sec->GetNumLam( );
  sLamina *Lamina = Sec->GetLayup( );

  // Evaluate [Mb] matrix through integration.

  double zb = -Sec->GetThickness( )/2.0;
  double I0 = 0.0;
  double I1 = 0.0;
  double I2 = 0.0;
  for (int i = 0; i < nlam; i++)
  {
    cMaterial *mat = Lamina[i].Mat;
    double rho = mat->GetDensity( );
    double zt = zb + Lamina[i].Thk;
    I0 += rho*(zt - zb);
    I1 += rho*(zt*zt - zb*zb)/2.0;
    I2 += rho*(zt*zt*zt - zb*zb*zb)/3.0;
    zb = zt;
  }
  Mb.Zero( );
  Mb[0][0] = Mb[1][1] = Mb[2][2] = I0;
  Mb[3][3] = Mb[4][4] = I2;
  Mb[0][4] = Mb[4][0] = I1;
  Mb[1][3] = Mb[3][1] = -I1;
}

// -------------------------------------------------------------------------
// Protected methods:
//

// =============================== CalcCz =================================

void cSecAnLaminatedShell :: CalcCz(double z, cMatrix &Q, cMatrix &Cz)
{
  Cz.Zero( );

  // Membrane and bending terms.

  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
    {
      Cz[i][j] = Q[i][j];                   // [A]
      Cz[3+i][j] = Cz[i][3+j] = z*Q[i][j];  // [B]
      Cz[3+i][3+j] = z*z*Q[i][j];           // [D]
    }

  // Transverse shear terms.

  double ks = 5.0/6.0;
  for (int i = 3; i < 5; i++)
    for (int j = 3; j < 5; j++)
      Cz[3+i][3+j] = ks*Q[i][j];            // [G]
}


// -------------------------------------------------------------------------
// Class cSecAnLaminatedSolid:
// -------------------------------------------------------------------------

// ========================= cSecAnLaminatedSolid ==========================

cSecAnLaminatedSolid :: cSecAnLaminatedSolid(cElement *elm, cIntPoint *ipt,
                        sLamina *ply) : cSecAnalysis(elm, ipt)
{
  Anm = elm->GetAnModel( );
  Ipt = ipt;
  Ply = ply;
  Cmod = cConstModel::Create(elm, ply->Mat);
}

// ========================= ~ccSecAnLaminatedSolid ========================

cSecAnLaminatedSolid :: ~cSecAnLaminatedSolid(void)
{
  delete Cmod;
}

// ============================ LaminaSys ==================================

void cSecAnLaminatedSolid :: LaminaSys(cVector& e1, cVector& e2, cVector& e3)
{
  // Evaluate the unit vector normal to the laminate => {e3}.

  cShape *shp = Elm->GetShape( );
  int nnode = shp->GetNumNode( );
  sNodeCoord *coord = new sNodeCoord[nnode];
  sNatCoord p = Ipt->GetCoord( );
  shp->NodalCoord(coord);
  shp->tVector(p, coord, e3);
  delete []coord;

//  cout << "{e3} = ";
//  e3.Print("%f ");
//  cout << "\n";

  // Get the laminate orientation vector (1-axis).

  cVector a(3);
  Sec->GetOrientationVector(a);
  a.Normalize( );
//  double lena = a.Length( );
//  a /= lena;
  //cout << "{a} = ";
  //a.Print("%f ");
  //cout << "\n";

  // Evaluate {e2} = {e3} x {a}.

  CrossProd(e3, a, e2);
  //cout << "{e2} = ";
  //e2.Print("%f ");
  //cout << "\n";

  // Modify {a} to be normal to {e2} and {e3}.

  CrossProd(e2, e3, a);
  //cout << "{a} = ";
  //a.Print("%f ");
  //cout << "\n";

  // Evaluate the [S] matrix.

  cMatrix S(3, 3);
  S[0][0] =  0;
  S[0][1] = -e3[2];
  S[0][2] =  e3[1];
  S[1][0] =  e3[2];
  S[1][1] =  0;
  S[1][2] = -e3[0];
  S[2][0] = -e3[1];
  S[2][1] =  e3[0];
  S[2][2] =  0;

  // Evaluate [S]^2 = [S]*[S].

  cMatrix S2(3, 3);
  S2 = S*S;

  // Initialize [R] = [I].

  cMatrix R(3, 3);
  R[0][0] =  1;
  R[0][1] =  0;
  R[0][2] =  0;
  R[1][0] =  0;
  R[1][1] =  1;
  R[1][2] =  0;
  R[2][0] =  0;
  R[2][1] =  0;
  R[2][2] =  1;

  // [R] = [I] + sin(theta)*[S] + (1 - cos(theta))*[S]*[S].

  R += sin(Ply->Ang*PI/180.0)*S;
  R += (1.0 - cos(Ply->Ang*PI/180.0))*S2;

  // Evaluate {e1} = [R]{a}

  e1 = R*a;
  //cout << "{e1} = ";
  //e1.Print("%f ");
  //cout << "\n";

  // Evaluate {e2} = {e3} x {e1}.

  CrossProd(e3, e1, e2);
/*
  cout << "ang = %f\n", Ply->Ang;
  cout << "{e1} = ";
  e1.Print("%f ");
  cout << "\n";

  cout << "{e2} = ";
  e2.Print("%f ");
  cout << "\n";

  cout << "{e3} = ";
  e3.Print("%f ");
  cout << "\n";
*/
}

// =============================== Stress ==================================

int cSecAnLaminatedSolid :: Stress(cVector &eps, cVector &sig)
{
  // Evaluate the lamina vectors.

  cVector e1(3),e2(3),e3(3);
  LaminaSys(e1, e2, e3);

  // Evaluate the transformation matrix.

  cMatrix T(6,6);
  Anm->GetMecMod( )->TMatrix(e1, e2, e3, T);

  // Compute the local strains.

  cVector epsl(6);
  epsl = T*eps;

  // Evaluate the local stresses using the associated Constitutive Model.

  cVector sigl(6);
  Cmod->Stress(epsl, sigl);

  // Return the global or local stresses.

  if (cControl :: GetLamLocStr( ) && Printing)  // Local stresses
    sig = sigl;
  else           // Global stresses
    sig = t(T)*sigl;

  return(1);
}

// =============================== CMatrix =================================

void cSecAnLaminatedSolid :: CMatrix(cMatrix &C)
{

  // Get local vectors.

  cVector e1(3),e2(3),e3(3);
  LaminaSys(e1, e2, e3);

  // Get the transformation matrix

  cMatrix T(6, 6);
  Anm->GetMecMod( )->TMatrix(e1, e2, e3, T);

  // Get the local constitutive matrix from the associated Constitutive Model

  cMatrix lC(6, 6);
  Cmod->TangMat(lC);

  // Evaluate the global constitutive matrix

  C.Zero( );
  MatTripBtCB(T, lC, 1.0, C);
}

// ============================= UpdateState ===============================

void cSecAnLaminatedSolid :: UpdateState(void)
{
  // Update the state of the associated Constitutive Model.

  Cmod->UpdateState( );
}


// -------------------------------------------------------------------------
// Class cSecAnPlFrameRectFiber:
// -------------------------------------------------------------------------

// ======================== cSecAnPlFrameRectFiber =========================

cSecAnPlFrameRectFiber :: cSecAnPlFrameRectFiber(cElement *elm, cIntPoint *ipt) :
                          cSecAnalysis(elm, ipt)
{
  NumFiber = Sec->NumIntPntY( );

  // Create the constitutive model of each fiber.

  cMaterial *mat = Sec->GetMaterial( );
  Cmod = new cConstModel*[NumFiber];
  for (int i = 0; i < NumFiber; i++)
    Cmod[i] = cConstModel::Create(elm, mat);
}

// ======================== ~cSecAnPlFrameRectFiber ========================

cSecAnPlFrameRectFiber :: ~cSecAnPlFrameRectFiber(void)
{
  for (int i = 0; i < NumFiber; i++) delete Cmod[i];
  delete []Cmod;
}

// ================================ Stress =================================

int cSecAnPlFrameRectFiber :: Stress(cVector &eps, cVector &sig)
{
  // Get section data.

  cSecBarRect *sec = (cSecBarRect *)Sec;
  double b = sec->b;
  double h = sec->h;

  // Compute the stress resultants using the fiber (layer) method.

  cVector el(1);
  cVector sl(1);
  double dy = h/NumFiber;         // Layer height
  double Al = b*dy;               // Layer area
  double yl = -h/2 + dy/2;        // Center of the bottom layer
  double N  = 0.0;
  double M  = 0.0;
  for (int i = 0; i < NumFiber; i++)
  {
    el[0] = eps[0] - yl*eps[1];               // Layer strain
    if (!Cmod[i]->Stress(el, sl)) return(0);  // Layer stress
    N  += sl[0]*Al;                           // Normal force
    M  -= sl[0]*Al*yl;                        // Bending moment
    yl += dy;                                 // Next layer
  }

  // Return the stress resultants.

  sig[0] = N;
  sig[1] = M;

  return(1);
}

// =============================== CMatrix =================================

void cSecAnPlFrameRectFiber :: CMatrix(cMatrix &C)
{
  // Get section data.

  cSecBarRect *sec = (cSecBarRect *)Sec;
  double b = sec->b;
  double h = sec->h;

  // Compute the section stiffness by the layer method.

  cMatrix Et(1, 1);
  double dy = h/NumFiber;         // Layer height
  double Al = b*dy;               // Layer area
  double yl = -h/2 + dy/2;        // Center of the bottom layer
  double EA = 0.0;
  double ES = 0.0;
  double EI = 0.0;
  for (int i = 0; i < NumFiber; i++)
  {
    Cmod[i]->TangMat(Et);         // Layer tangent modulus
    EA += Et[0][0]*Al;
    ES -= Et[0][0]*Al*yl;
    EI += Et[0][0]*Al*(yl*yl + dy*dy/12.0);
    yl += dy;                     // Next layer
  }

  // Return the constitutive matrix.

  C[0][0] = EA;
  C[0][1] = C[1][0] = ES;
  C[1][1] = EI;
}

// ============================= UpdateState ===============================

void cSecAnPlFrameRectFiber :: UpdateState(void)
{
  // Update the state of each constitutive model.

  for (int i = 0; i < NumFiber; i++)
    Cmod[i]->UpdateState( );
}


// -------------------------------------------------------------------------
// Class cSecAnPlFrameRectGauss:
// -------------------------------------------------------------------------

// ======================== cSecAnPlFrameRectGauss =========================

cSecAnPlFrameRectGauss :: cSecAnPlFrameRectGauss(cElement *elm, cIntPoint *ipt) :
                          cSecAnalysis(elm, ipt)
{
  IntSecPnt = cIntPoint :: CreateLinePoints(Sec->NumIntPntY( ), GAUSS, &NumSecPnt);

  // Create the constitutive model of each integration point.

  cMaterial *mat = Sec->GetMaterial( );
  Cmod = new cConstModel*[NumSecPnt];
  for (int i = 0; i < NumSecPnt; i++)
    Cmod[i] = cConstModel::Create(elm, mat);
}

// ======================== ~cSecAnPlFrameRectGauss ========================

cSecAnPlFrameRectGauss :: ~cSecAnPlFrameRectGauss(void)
{
  delete []IntSecPnt;
  for (int i = 0; i < NumSecPnt; i++) delete Cmod[i];
  delete []Cmod;
}

// ================================ Stress =================================

int cSecAnPlFrameRectGauss :: Stress(cVector &eps, cVector &sig)
{
  // Get section data.

  cSecBarRect *sec = (cSecBarRect *)Sec;
  double b = sec->b;
  double h = sec->h;

  // Compute the stress resultants using Gauss integration.

  cVector ep(1);
  cVector sp(1);
  double N = 0.0;
  double M = 0.0;
  for (int i = 0; i < NumSecPnt; i++)
  {
    sNatCoord p = IntSecPnt[i].GetCoord( );          // Parametric coord
    double yp =   (h/2.0)*p.r;                       // Coordinate
    double cp = b*(h/2.0)*IntSecPnt[i].GetWeight( ); // Coeficient
    ep[0] = eps[0] - yp*eps[1];                      // Strain
    if (!Cmod[i]->Stress(ep, sp)) return(0);         // Stress
    N  +=    sp[0]*cp;                               // Normal force
    M  -= yp*sp[0]*cp;                               // Bending moment
  }

  // Return the stress resultants.

  sig[0] = N;
  sig[1] = M;

  return(1);
}

// =============================== CMatrix =================================

void cSecAnPlFrameRectGauss :: CMatrix(cMatrix &C)
{
  // Get section data.

  cSecBarRect *sec = (cSecBarRect *)Sec;
  double b = sec->b;
  double h = sec->h;

  // Compute the section stiffness using Gauss integration.

  cMatrix Et(1, 1);
  double EA = 0.0;
  double ES = 0.0;
  double EI = 0.0;
  for (int i = 0; i < NumSecPnt; i++)
  {
    sNatCoord p = IntSecPnt[i].GetCoord( );          // Parametric coord
    double yp =   (h/2.0)*p.r;                       // Coordinate
    double cp = b*(h/2.0)*IntSecPnt[i].GetWeight( ); // Coeficient
    Cmod[i]->TangMat(Et);                            // Tangent modulus
    EA +=       Et[0][0]*cp;
    ES -=    yp*Et[0][0]*cp;
    EI += yp*yp*Et[0][0]*cp;
  }

  // Return the constitutive matrix.

  C[0][0] = EA;
  C[0][1] = C[1][0] = ES;
  C[1][1] = EI;
}

// ============================= UpdateState ===============================

void cSecAnPlFrameRectGauss :: UpdateState(void)
{
  // Update the state of each constitutive model.

  for (int i = 0; i < NumSecPnt; i++)
    Cmod[i]->UpdateState( );
}


// -------------------------------------------------------------------------
// Class cSecAnPlFrameRCRect:
// -------------------------------------------------------------------------

// ========================= cSecAnPlFrameRCRect ===========================

cSecAnPlFrameRCRect :: cSecAnPlFrameRCRect(cElement *elm, cIntPoint *ipt) :
                       cSecAnalysis(elm, ipt)
{
  NumFiber = Sec->NumIntPntY( );

  // Create the constitutive model of each fiber.

  cMaterial *mat = Sec->GetMaterial( );
  ModFiber = new cConstModel*[NumFiber];
  for (int i = 0; i < NumFiber; i++)
    ModFiber[i] = cConstModel::Create(elm, mat);

  // Create the constitutive model of each rebar and material void.

  cSecBarRCRect *sec = (cSecBarRCRect *)Sec;
  ModRebar = new cConstModel*[sec->NumBar];
  ModVoid  = new cConstModel*[sec->NumBar];
  for (int i = 0; i < sec->NumBar; i++)
  {
    ModRebar[i] = cConstModel::Create(elm, sec->MatBar[i]);
    ModVoid[i]  = cConstModel::Create(elm, mat);
  }
}

// ========================= ~cSecAnPlFrameRCRect ==========================

cSecAnPlFrameRCRect :: ~cSecAnPlFrameRCRect(void)
{
  for (int i = 0; i < NumFiber; i++) delete ModFiber[i];
  delete []ModFiber;

  cSecBarRCRect *sec = (cSecBarRCRect *)Sec;
  for (int i = 0; i < sec->NumBar; i++) delete ModRebar[i];
  for (int i = 0; i < sec->NumBar; i++) delete ModVoid[i];
  delete []ModRebar;
  delete []ModVoid;
}

// ================================ Stress =================================

int cSecAnPlFrameRCRect :: Stress(cVector &eps, cVector &sig)
{
  // Get section data.

  cSecBarRCRect *sec = (cSecBarRCRect *)Sec;
  double b = sec->b;
  double h = sec->h;

  // Compute the stress resultants using the fiber (layer) method.

  cVector el(1);
  cVector sl(1);
  double dy = h/NumFiber;         // Layer height
  double Al = b*dy;               // Layer area
  double yl = -h/2 + dy/2;        // Center of the bottom layer
  double N  = 0.0;
  double M  = 0.0;
  for (int i = 0; i < NumFiber; i++)
  {
    el[0] = eps[0] - yl*eps[1];                  // Layer strain
    if (!ModFiber[i]->Stress(el, sl)) return(0); // Layer stress
    N  += sl[0]*Al;                              // Normal force
    M  -= sl[0]*Al*yl;                           // Bending moment
    yl += dy;                                    // Next layer
  }

  // Add the contribution of rebars.

  cVector eb(1);
  cVector sb(1);
  cVector sv(1);
  for (int i = 0; i < sec->NumBar; i++)
  {
    eb[0] = eps[0] - sec->YBar[i]*eps[1];        // Rebar strain
    if (!ModRebar[i]->Stress(eb, sb)) return(0); // Rebar stress
    if (!ModVoid[i]->Stress(eb, sv)) return(0);  // Void stress
    sb[0] -= sv[0];
    N += sb[0]*sec->ABar[i];                     // Normal force
    M -= sb[0]*sec->ABar[i]*sec->YBar[i];        // Bending moment
  }

  // Return the stress resultants.

  sig[0] = N;
  sig[1] = M;

  return(1);
}

// =============================== CMatrix =================================

void cSecAnPlFrameRCRect :: CMatrix(cMatrix &C)
{
  // Get section data.

  cSecBarRCRect *sec = (cSecBarRCRect *)Sec;
  double b = sec->b;
  double h = sec->h;

  // Compute the section stiffness by the layer method.

  cMatrix Et(1, 1);
  double dy = h/NumFiber;         // Layer height
  double Al = b*dy;               // Layer area
  double yl = -h/2 + dy/2;        // Center of the bottom layer
  double EA = 0.0;
  double ES = 0.0;
  double EI = 0.0;
  for (int i = 0; i < NumFiber; i++)
  {
    ModFiber[i]->TangMat(Et);     // Layer tangent modulus
    EA += Et[0][0]*Al;
    ES -= Et[0][0]*Al*yl;
    EI += Et[0][0]*Al*(yl*yl + dy*dy/12.0);
    yl += dy;                     // Next layer
  }

  // Add the contribution of rebars.

  cMatrix Ev(1, 1);
  for (int i = 0; i < sec->NumBar; i++)
  {
    ModRebar[i]->TangMat(Et);     // Rebar tangent modulus
    ModVoid[i]->TangMat(Ev);      // Void tangent modulus
    Et[0][0] -= Ev[0][0];
    EA += Et[0][0]*sec->ABar[i];
    ES -= Et[0][0]*sec->ABar[i]*sec->YBar[i];
    EI += Et[0][0]*sec->ABar[i]*sec->YBar[i]*sec->YBar[i];
  }

  // Return the constitutive matrix.

  C[0][0] = EA;
  C[0][1] = C[1][0] = ES;
  C[1][1] = EI;
}

// ============================= UpdateState ===============================

void cSecAnPlFrameRCRect :: UpdateState(void)
{
  // Update the state of each constitutive model.

  for (int i = 0; i < NumFiber; i++)
    ModFiber[i]->UpdateState( );

  cSecBarRCRect *sec = (cSecBarRCRect *)Sec;
  for (int i = 0; i < sec->NumBar; i++)
    ModRebar[i]->UpdateState( );
}

// ======================================================= End of file =====
