// -------------------------------------------------------------------------
// secbar.cpp - implementation of bar cross-section class.
// -------------------------------------------------------------------------
// Created:      22-Sep-2012     Evandro Parente Junior
//
// Modified:     06-Oct-2013     Evandro Parente Junior
//               Unification with Section class.
//
// Modified:     05-Feb-2015     Evandro Parente Junior
//               Elimination of integration methods.
//
// Modified:     08-Apr-2015     Elias Saraiva Barroso
//               Input/output operations using c++ streams.
// -------------------------------------------------------------------------

#include <iostream>
#include <math.h>
#include <stdlib.h>

using namespace std;

#include "secbar.h"
#include "material.h"
#include "vec.h"
#include "mat.h"
#include "gbldef.h"
#include "gblvar.h"

// -------------------------------------------------------------------------
// Public methods:
//

// ================================ cSecBar ================================

cSecBar :: cSecBar(int id) : cSection(id)
{
  A  = 0.0;
  Ay = 0.0;
  Az = 0.0;
  Iy = 0.0;
  Iz = 0.0;
  Jt = 0.0;
  Cw = 0.0;
}

// ================================ ~cSecBar ===============================

cSecBar :: ~cSecBar(void)
{
}


// -------------------------------------------------------------------------
// cSecBarHom class:
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Public methods:
//

// =============================== ~cSecBarHom =============================

cSecBarHom :: cSecBarHom(int id) : cSecBar(id)
{
}

// =============================== ~cSecBarHom =============================

cSecBarHom :: ~cSecBarHom(void)
{
}

// ============================== PreIntegrated ============================

bool cSecBarHom :: PreIntegrated(void)
{
  if (Mat->GetMecProp( )->GetType( ) == MAT_ELASTIC_ISOTROPIC)
    return(true);

  return(false);
}

// ================================== Read =================================

void cSecBarHom :: Read(void)
{
  // Read section material.

  int mat;
  in >> mat;
  Mat = cMaterial::GetMaterial(mat);
  if (!Mat)
  {
    cout << "Error in the input of section " << Label << "!\n";
    cout << "Invalid material!\n";
    exit(0);
  }

  // Read geometric data.

  ReadGeom( );
}

// ============================== SecStiffness2D ===========================

void cSecBarHom :: SecStiffness2D(cMatrix &C)
{
  // Get material parameters.

  int np = Mat->GetMecProp( )->NumParam( );
  double* matpar = new double[np];
  Mat->GetMecProp( )->GetParam(matpar);
  double E = matpar[0];

  // Return the constitutive matrix.

  C.Zero( );
  C[0][0] = E*A;
  C[1][1] = E*Iz;

  delete [] matpar;
}

// ============================== SecStiffness3D ===========================

void cSecBarHom :: SecStiffness3D(cMatrix &C)
{
  // Get material parameters.

  int np = Mat->GetMecProp( )->NumParam( );
  double* matpar = new double[np];
  Mat->GetMecProp( )->GetParam(matpar);
  double E  = matpar[0];
  double nu = matpar[0];
  double G  = 0.5*E/(1.0 + nu);

  // Compute the constitutive matrix.

  C.Zero( );
  C[0][0] = E*A;
  C[1][1] = E*Iy;
  C[2][2] = E*Iz;
  C[3][3] = G*Jt;

  delete [] matpar;
}


// -------------------------------------------------------------------------
// cSecGen class:
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Public methods:
//

// =============================== cSecBarGen ==============================

cSecBarGen :: cSecBarGen(int id) : cSecBarHom(id)
{
  Type = SEC_BAR_GEN;
}

// =============================== ~cSecBarGen =============================

cSecBarGen :: ~cSecBarGen(void)
{
}

// ================================ ReadGeom ===============================

void cSecBarGen :: ReadGeom(void)
{
  // Read section properties.

  double Hy,Hz;
  if (!(in >> A)  || !(in >> Ay) || !(in >> Az) || !(in >> Jt) || !(in >> Iy) ||
      !(in >> Iz) || !(in >> Cw) || !(in >> Hy) || !(in >> Hz))
  {
    cout << "Error in the input of section " << Label << "!\n";
    exit(0);
  }
}


// -------------------------------------------------------------------------
// cSecRect class:
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Public methods:
//

// ============================== ~cSecBarRect =============================

cSecBarRect :: cSecBarRect(int id) : cSecBarHom(id)
{
  Type = SEC_BAR_HOM_RECT;
}

// ============================== ~cSecBarRect =============================

cSecBarRect :: ~cSecBarRect(void)
{
}

// ================================ ReadGeom ===============================

void cSecBarRect :: ReadGeom(void)
{
  if (!(in >> b) || !(in >> h))
  {
    cout << "Error in the input of section " << Label << "!\n";
    exit(0);
  }

  // Compute the geometric properties.

  A  = b*h;
  Ay = 5.0*A/6.0;
  Az = Ay;
  Iy = b*b*b*h/12.0;
  Iz = b*h*h*h/12.0;
  Jt = Iy + Iz;      // Temporario: procurar expressao mais exata
  Cw = 0.0;
}


// -------------------------------------------------------------------------
// cSecBarCirc class:
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Public methods:
//

// =============================== cSecBarCirc =============================

cSecBarCirc :: cSecBarCirc(int id) : cSecBarHom(id)
{
  Type = SEC_BAR_HOM_CIRCLE;
}

// ============================== ~cSecBarCirc =============================

cSecBarCirc :: ~cSecBarCirc(void)
{
}

// ================================ ReadGeom ===============================

void cSecBarCirc :: ReadGeom(void)
{
  if (!(in >> R))
  {
    cout << "Error in the input of section " << Label << "!\n";
    exit(0);
  }

  // Compute the geometric properties.

  A  = PI*R*R;
  Ay = 0.9*A;
  Az = Ay;
  Iy = PI*R*R*R*R/4.0;
  Iz = Iy;
  Jt = Iy + Iz;
  Cw = 0.0;
}


// -------------------------------------------------------------------------
// cSecTube class:
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Public methods:
//

// ============================== ~cSecBarTube =============================

cSecBarTube :: cSecBarTube(int id) : cSecBarHom(id)
{
  Type = SEC_BAR_HOM_TUBE;
}

// ============================== ~cSecBarTube =============================

cSecBarTube :: ~cSecBarTube(void)
{
}

// ================================ ReadGeom ===============================

void cSecBarTube :: ReadGeom(void)
{
  if (!(in >> R) || !(in >> t))
  {
    cout << "Error in the input of section " << Label << "!\n";
    exit(0);
  }

  // Compute the geometric properties.

  double r = R - t;  // Inner radius
  A  = PI*(R*R - r*r);
  Ay = A/2.0;
  Az = Ay;
  Iy = PI*(R*R*R*R - r*r*r*r)/4.0;
  Iz = Iy;
  Jt = Iy + Iz;
  Cw = 0.0;
}


// -------------------------------------------------------------------------
// cSecBarB3DCoupled class:
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Public methods:
//

// ============================ cSecBarB3DCoupled ==========================

cSecBarB3DCoupled :: cSecBarB3DCoupled(int id) : cSecBar(id)
{
  Type = SEC_BAR_B3DCOUPLED;
  C.Resize(4, 4);
}

// =========================== ~cSecBarB3DCoupled ==========================

cSecBarB3DCoupled :: ~cSecBarB3DCoupled(void)
{
}

// ================================== Read =================================

void cSecBarB3DCoupled :: Read(void)
{
  // Read section material.

  int mat;
  if (in >> mat)
    Mat = cMaterial::GetMaterial(mat);
  if (!Mat)
  {
    cout << "Error in the input of section " << Label << "!\n";
    cout << "Invalid material!\n";
    exit(0);
  }

  if (Mat->GetMecProp( )->GetType( ) != MAT_ELASTIC_ISOTROPIC)
  {
    cout << "Error in the input of section " << Label << "!\n";
    cout << "Beam coupled section can be used only with linear elastic material!\n";
    exit(0);
  }

  // Read section properties.

  double Hy,Hz;
  if (!(in >> A)  || !(in >> Ay) || !(in >> Az) || !(in >> Jt) || !(in >> Iy) ||
      !(in >> Iz) || !(in >> Cw) || !(in >> Hy) || !(in >> Hz))
  {
    cout << "Error in the input of section " << Label << "!\n";
    exit(0);
  }

  // Read section matrix.

  for (int i = 0; i < 4; i++)
    for (int j = 0; j < 4; j++) in >> C[i][j];
}

// ============================== SecStiffness2D ===========================

void cSecBarB3DCoupled :: SecStiffness2D(cMatrix &Ct)
{
  Ct.Zero( );
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++) Ct[i][j] = C[i][j];
}

// ============================== SecStiffness3D ===========================

void cSecBarB3DCoupled :: SecStiffness3D(cMatrix &Ct)
{
  Ct = C;
}


// -------------------------------------------------------------------------
// cSecRC class:
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Public methods:
//

// =============================== ~cSecBarRC ==============================

cSecBarRC :: cSecBarRC(int id) : cSecBar(id)
{
}

// =============================== ~cSecBarRC ==============================

cSecBarRC :: ~cSecBarRC(void)
{
}

// ================================== Read =================================

void cSecBarRC :: Read(void)
{
  // Read section material.

  int idmat;
  if (in >> idmat)
    Conc = cMaterial::GetMaterial(idmat);
  if (!Conc)
  {
    cout << "Error in the input of section " << Label << "!\n";
    cout << "Invalid material!\n";
    exit(0);
  }

  // Read geometric data.

  ReadGeom( );

  // Read rebar data.

  in >> NumBar;
  XBar = new double[NumBar];
  YBar = new double[NumBar];
  ABar = new double[NumBar];
  MatBar = new cMaterial*[NumBar];
  for (int i = 0; i < NumBar; i++)
  {
    in >> idmat;
    MatBar[i] = cMaterial::GetMaterial(idmat);
    if (!MatBar[i])
    {
      cout << "Error in the input of section " << Label << "!\n";
      cout << "Invalid material of rebar " << i+1 << "!\n";
      exit(0);
    }

    if (!(in >> XBar[i])  || !(in >> YBar[i]) || !(in >> ABar[i]))
    {
      cout << "Error in the input of section " << Label << "!\n";
      cout << "Invalid material of rebar " << i+1 << "!\n";
      exit(0);
    }
  }

  // Transform the coordinates to the centroid system.

  double xc,yc;
  GetCentroid(&xc, &yc);
  for (int i = 0; i < NumBar; i++)
  {
    XBar[i] -= xc;
    YBar[i] -= yc;
  }
}


// -------------------------------------------------------------------------
// cSecBarRCRect class:
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Public methods:
//

// ============================= ~cSecBarRCRect ============================

cSecBarRCRect :: cSecBarRCRect(int id) : cSecBarRC(id)
{
  Type = SEC_BAR_RC_RECT;
}

// ============================= ~cSecBarRCRect ============================

cSecBarRCRect :: ~cSecBarRCRect(void)
{
}

// ================================ ReadGeom ===============================

void cSecBarRCRect :: ReadGeom(void)
{
  if (!(in >> b)  || !(in >> h))
  {
    cout << "Error in the input of section " << Label << "!\n";
    exit(0);
  }

  // Compute the geometric properties.

  A  = b*h;
  Ay = 5.0*A/6.0;
  Az = Ay;
  Iy = b*b*b*h/12.0;
  Iz = b*h*h*h/12.0;
  Jt = Iy + Iz;      // Temporario: procurar expressao mais exata
  Cw = 0.0;
}

// =============================== GetCentroid =============================

void cSecBarRCRect :: GetCentroid(double *xc, double *yc)
{
  *xc = 0.0;
  *yc = 0.0;
}

// ======================================================= End of file =====
