// -------------------------------------------------------------------------
// material.cpp - implementation of the Material class.
// -------------------------------------------------------------------------
// Created:      14-Jul-2005     Evandro Parente Junior
//
// Modified:     06-Nov-2012     Evandro Parente Junior
//               Creation of NBR6118 material.
//
// Modified:     18-Jul-2013     Iuri Barcelos Rocha
//               Creation of FRC materials.
//
// Modified:     08-Nov-2013     Leandro Soares Moreira
//               Creation of EuroCEB material.
//
// Modified:     25-Fev-2014     Evandro Parente Junior
//               NBR6118 modified to NBRCEB.
//
// Modified:     09-May-2014     Edson Dantas/Evandro Parente
//               Creation of CZM material.
//
// Modified:     04-Aug-2014     Carlos David/Evandro Parente
//               Creation of ElastoPlastic materials.
//
// Modified:     08-Apr-2015     Elias Saraiva Barroso
//               Input/output operations using c++ streams.
//
// Modified:     02-Mar-2020     Bergson Matias/Evandro Parente
//               Creation of Mazars, MuModel and Lee-Fenves materials.
//
// Modified:     01-Feb-2023     Renan Melo Barros
//               Creation of elastoplastic material with isotropic hardening.
// -------------------------------------------------------------------------

#include <iostream>
#include <stdlib.h>

using namespace std;

#include "material.h"
#include "gblvar.h"

// -------------------------------------------------------------------------
// Static variables:
//
int        cMaterial :: NumMat  = 0;
cMaterial* cMaterial :: VecMat  = 0;

// -------------------------------------------------------------------------
// Protected methods:
//

// ============================== SetMecProp ===============================

void cMaterial :: SetMecProp(cMatMec *mecmat)
{
  MecProp = mecmat;
  MecProp->mat = this;
}

// ============================== SetThermProp =============================

void cMaterial :: SetThermProp(cMatTherm *tmat)
{
  ThermProp      = tmat;
  ThermProp->mat = this;
}

// -------------------------------------------------------------------------
// Public methods:
//

// ============================== ReadNumMat ===============================

void cMaterial :: ReadNumMat(void)
{
  // Read the number of materials.

  if (!(in >> NumMat) || (NumMat <= 0))
  {
    cout << "Error in the input of the number of materials!\n";
    exit(0);
  }

  // Alloc the array of materials.

  VecMat = new cMaterial[NumMat];

  // Set material labels.

  for(int i = 0; i < NumMat; ++i) VecMat[i].Label = i+1;
}

// =============================== ReadElastIso ============================

void cMaterial :: ReadElastIso(void)
{
  // Read the number of elastic isotropic materials.

  int n;
  if (!(in >> n))
  {
    cout << "Error in the input of the number of isotropic materials!\n";
    exit(0);
  }

  // Read each material.

  int id;
  for (int i = 0; i < n; i++)
  {
    if (!(in >> id) || id < 1 || id > NumMat)
    {
      cout << "Error in the input of elastic isotropic material " << i+1 << "!\n";
      exit(0);
    }
    VecMat[id-1].SetMecProp(new cMatElastIso( ));
    VecMat[id-1].MecProp->Read( );
  }
}

// ============================== ReadElastOrtho ===========================

void cMaterial :: ReadElastOrtho(void)
{
  // Read the number of elastic orthotropic materials.

  int n;
  if (!(in >> n))
  {
    cout << "Error in the input of the number of orthotropic materials!\n";
    exit(0);
  }

  // Read each material.

  int id;
  for (int i = 0; i < n; i++)
  {
    if (!(in >> id) || id < 1 || id > NumMat)
    {
      cout << "Error in the input of elastic orthotropic material " << i+1 << "!\n";
      exit(0);
    }
    VecMat[id-1].SetMecProp(new cMatElastOrtho( ));
    VecMat[id-1].MecProp->Read( );
  }
}

// ============================= ReadElastNonlin ===========================

void cMaterial :: ReadElastNonlin(void)
{
  // Read the number of nonlinear elastic (isotropic) materials.

  int n;
  if (!(in >> n))
  {
    cout << "Error in the input of the number of nonlinear materials!\n";
    exit(0);
  }

  // Read each material.

  int id;
  for (int i = 0; i < n; i++)
  {
    if (!(in >> id) || id < 1 || id > NumMat)
    {
      cout << "Error in the input of nonlinear elastic material " << i+1 << "!\n";
      exit(0);
    }
    VecMat[id-1].SetMecProp(new cMatElastNonlin( ));
    VecMat[id-1].MecProp->Read( );
  }
}

// ================================= ReadFGM ===============================

void cMaterial :: ReadFGM(void)
{
  // Read the number of FGMs.

  int n;
  if (!(in >> n))
  {
    cout << "Error in the input of the number of FGM materials!\n";
    exit(0);
  }

  // Read each material.

  int id;
  for (int i = 0; i < n; i++)
  {
    if (!(in >> id) || id < 1 || id > NumMat)
    {
      cout << "Error in the input of FGM material " << i+1 << "!\n";
      exit(0);
    }
    VecMat[id-1].SetMecProp(new cMatFGM( ));
    VecMat[id-1].MecProp->Read( );
  }
}

// ================================ ReadFGMTTO =============================

void cMaterial :: ReadFGMTTO(void)
{
  // Read the number of FGMs.

  int n;
  if (!(in >> n))
  {
    cout << "Error in the input of the number of FGM_TTO materials!\n";
    exit(0);
  }

  // Read each material.

  int id;
  for (int i = 0; i < n; i++)
  {
    if (!(in >> id) || id < 1 || id > NumMat)
    {
      cout << "Error in the input of FGM_TTO material " << i+1 << "!\n";
      exit(0);
    }
    VecMat[id-1].SetMecProp(new cMatFGMTTO( ));
    VecMat[id-1].MecProp->Read( );
  }
}

// =========================== ReadConcreteEuroCEB ==========================

void cMaterial :: ReadConcEuroCEB(void)
{
  // Read the number of concrete materials.

  int n;
  if (!(in >> n))
  {
    cout << "Error in the input of the number of concrete materials!\n";
    exit(0);
  }

  // Read each material.

  int id;
  for (int i = 0; i < n; i++)
  {
    if (!(in >> id) || id < 1 || id > NumMat)
    {
      cout << "Error in the input of concrete material " << i+1 << "!\n";
      exit(0);
    }
    VecMat[id-1].SetMecProp(new cMatConcEuroCEB( ));
    VecMat[id-1].MecProp->Read( );
  }
}

// ========================== ReadConcreteNBRECB ===========================

void cMaterial :: ReadConcNBRCEB(void)
{
  // Read the number of concrete materials.

  int n;
  if (!(in >> n))
  {
    cout << "Error in the input of the number of concrete materials!\n";
    exit(0);
  }

  // Read each material.

  int id;
  for (int i = 0; i < n; i++)
  {
    if (!(in >> id) || id < 1 || id > NumMat)
    {
      cout << "Error in the input of concrete material " << i+1 << "!\n";
      exit(0);
    }
    VecMat[id-1].SetMecProp(new cMatConcNBRCEB( ));
    VecMat[id-1].MecProp->Read( );
  }
}

// ============================= ReadViscoElastic ==========================

void cMaterial :: ReadViscoElastic(void)
{
  // Read the number of viscoelastic materials.

  int n;
  if (!(in >> n))
  {
    cout << "Error in the input of the number of viscoelastic materials!\n";
    exit(0);
  }

  // Read each material.

  int id;
  for (int i = 0; i < n; i++)
  {
    if (!(in >> id) || id < 1 || id > NumMat)
    {
      cout << "Error in the input of viscoelastic material " << i+1 << "!\n";
      exit(0);
    }
    VecMat[id-1].SetMecProp(new cMatViscoElastic( ));
    VecMat[id-1].MecProp->Read( );
  }
}

// ============================ ReadPlastIsoLinHard ========================

void cMaterial :: ReadPlastIsoLinHard(void)
{
  // Read the number of elastoplastic  materials.

  int n;
  if (!(in >> n))
  {
    cout << "Error in the input of the number of elastoplastic materials!\n";
    exit(0);
  }

  // Read each material.

  int id;
  for (int i = 0; i < n; i++)
  {
    if (!(in >> id) || id < 1 || id > NumMat)
    {
      cout << "Error in the input of elastoplastic material " << i+1 << "!\n";
      exit(0);
    }
    VecMat[id-1].SetMecProp(new cMatPlastIsoLinHard( ));
    VecMat[id-1].MecProp->Read( );
  }
}

// =========================== ReadPlastLinHard =========================

void cMaterial :: ReadPlastLinHard(void)
{
  // Read the number of elastoplatic materials with linear hardening.

  int n;
  if (!(in >> n))
  {
    cout << "Error in the input of the number of ElastoPlastic materials!\n";
    exit(0);
  }

  // Read each material.

  int id;
  for (int i = 0; i < n; i++)
  {
    if (!(in >> id) || id < 1 || id > NumMat)
    {
      cout << "Error in the input of number of ElastoPlastic material " << i+1 << "!\n";
      exit(0);
    }
    VecMat[id-1].SetMecProp(new cMatPlastLinHard( ));
    VecMat[id-1].MecProp->Read( );
  }
}

// ========================== ReadPlastIsoNonlinHard =======================

void cMaterial :: ReadPlastIsoNonlinHard(void)
{
  // Read the number of elastoplatic materials with nonlinear hardening.

  int n;
  if (!(in >> n))
  {
    cout << "Error in the input of the number of ElastoPlastic materials!\n";
    exit(0);
  }

  // Read each material.

  int id;
  for (int i = 0; i < n; i++)
  {
    if (!(in >> id) || id < 1 || id > NumMat)
    {
      cout << "Error in the input of number of ElastoPlastic material " << i+1 << "!\n";
      exit(0);
    }
    VecMat[id-1].SetMecProp(new cMatPlastIsoNonlinHard( ));
    VecMat[id-1].MecProp->Read( );
  }
}

// ========================== ReadPlastIsoExpHard ==========================

void cMaterial :: ReadPlastIsoExpHard(void)
{
  // Read the number of elastoplatic materials with exponential hardening.

  int n;
  if (!(in >> n))
  {
    cout << "Error in the input of the number of ElastoPlastic material!\n";
    exit(0);
  }

  // Read each material.

  int id;
  for (int i = 0; i < n; i++)
  {
    if (!(in >> id) || id < 1 || id > NumMat)
    {
      cout << "Error in the input of number of ElastoPlastic material " << i+1 << "!\n";
      exit(0);
    }
    VecMat[id-1].SetMecProp(new cMatPlastIsoExpHard( ));
    VecMat[id-1].MecProp->Read( );
  }
}

// =========================== ReadDamageMazars ============================

void cMaterial :: ReadDamageMazars(void)
{
  // Read the number of Damage Mazars.

  int n;
  if (!(in >> n))
  {
    cout << "Error in the input of the number of Damage Mazars material!\n";
    exit(0);
  }

  // Read each material.

  int id;
  for (int i = 0; i < n; i++)
  {
    if (!(in >> id) || id < 1 || id > NumMat)
    {
      cout << "Error in the input of number of Damage Mazars material " << i+1 << "!\n";
      exit(0);
    }
    VecMat[id-1].SetMecProp(new cMatDamageMazars( ));
    VecMat[id-1].MecProp->Read( );
 }
}

// =========================== ReadDamageMuModel ===========================

void cMaterial :: ReadDamageMuModel(void)
{
  // Read the number of Damage Mu-Model Mazars.
  int n;
  if (!(in >> n))
  {
    cout << "Error in the input of the number of Damage Mu-Model material!\n";
    exit(0);
  }

  // Read each material.

  int id;
  for (int i = 0; i < n; i++)
  {
    if (!(in >> id) || id < 1 || id > NumMat)
    {
      cout << "Error in the input of number of Damage Mu-Model material " << i+1 << "!\n";
      exit(0);
    }
    VecMat[id-1].SetMecProp(new cMatDamageMuModel ( ));
    VecMat[id-1].MecProp->Read( );
 }
}

// ============================= ReadLeeFenves =============================

void cMaterial :: ReadLeeFenves(void)
{
  // Read the number of Plastic Damage Lee-Fenves Material.
  int n;
  if (!(in >> n))
  {
    cout << "Error in the input of the number of Lee-Fenves material!\n";
    exit(0);
  }

  // Read each material.

  int id;
  for (int i = 0; i < n; i++)
  {
    if (!(in >> id) || id < 1 || id > NumMat)
    {
      cout << "Error in the input of number of Lee-Fenves material " << i+1 << "!\n";
      exit(0);
    }
    VecMat[id-1].SetMecProp(new cMatLeeFenves( ));
    VecMat[id-1].MecProp->Read( );
 }
}

// =========================== ReadProgFailHashin ==========================

void cMaterial :: ReadProgFailHashin(void)
{
  // Read the number of composite materials.

  int n;
  if (!(in >> n))
  {
    cout << "Error in the input of the number of FRC materials (Hashin)!\n";
    exit(0);
  }

  // Read each material.

  int id;
  for (int i = 0; i < n; i++)
  {
    if (!(in >> id) || id < 1 || id > NumMat)
    {
      cout << "Error in the input of FRC material (Hashin) " << i+1 << "!\n";
      exit(0);
    }
    VecMat[id-1].SetMecProp(new cMatProgFailHashin( ));
    VecMat[id-1].MecProp->Read( );
  }
}

// =========================== ReadProgFailTsai ============================

void cMaterial :: ReadProgFailTsai(void)
{
  // Read the number of composite materials.

  int n;
  if (!(in >> n))
  {
    cout << "Error in the input of the number of FRC materials (Tsai)!\n";
    exit(0);
  }

  // Read each material.

  int id;
  for (int i = 0; i < n; i++)
  {
    if (!(in >> id) || id < 1 || id > NumMat)
    {
      cout << "Error in the input of FRC material (Tsai) " << i+1 << "!\n";
      exit(0);
    }
    VecMat[id-1].SetMecProp(new cMatProgFailTsai( ));
    VecMat[id-1].MecProp->Read( );
  }
}

// ======================= ReadProgFailEngelstad ==========================

void cMaterial :: ReadProgFailEngelstad(void)
{
  // Read the number of composite materials.

  int n;
  if (!(in >> n))
  {
    cout << "Error in the input of the number of FRC materials (Engelstad)!\n";
    exit(0);
  }

  // Read each material.

  int id;
  for (int i = 0; i < n; i++)
  {
    if (!(in >> id) || id < 1 || id > NumMat)
    {
      cout << "NumMat: " << NumMat << " id: " << id << "\n";
      cout << "Error in the input of FRC material (Engelstad) " << i+1 << "!\n";
      exit(0);
    }
    VecMat[id-1].SetMecProp(new cMatProgFailEngelstad( ));
    VecMat[id-1].MecProp->Read( );
  }
}

// ======================== ReadCZMModeIBilinear ===========================

void cMaterial :: ReadCZMModeIBilinear(void)
{
  // Read the number of CZM materials.

  int n;
  if (!(in >> n))
  {
    cout << "Error in the input of the number of CZM materials (ModeI Bilinear)!\n";
    exit(0);
  }

  // Read each material.

  int id;
  for (int i = 0; i < n; i++)
  {
    if (!(in >> id) || id < 1 || id > NumMat)
    {
      cout << "NumMat: " << NumMat << " id: " << id << "\n";
      cout << "Error in the input of CZM material (ModeI Bilinear) " << i+1 << "!\n";
      exit(0);
    }
    VecMat[id-1].SetMecProp(new cMatCZMModeIBilinear( ));
    VecMat[id-1].MecProp->Read( );
  }
}

// ======================= ReadCZMModeIExponential =========================
void cMaterial :: ReadCZMModeIExponential(void)
{
  // Read the number of CZM materials.

  int n;
  if (!(in >> n))
  {
    cout << "Error in the input of the number of CZM materials (ModeI Exponential)!\n";
    exit(0);
  }

  // Read each material.

  int id;
  for (int i = 0; i < n; i++)
  {
    if (!(in >> id) || id < 1 || id > NumMat)
    {
      cout << "NumMat: " << NumMat << " id: " << id << "\n";
      cout << "Error in the input of CZM material (ModeI Exponential) " << i+1 << "!\n";
      exit(0);
    }
    VecMat[id-1].SetMecProp(new cMatCZMModeIExponential( ));
    VecMat[id-1].MecProp->Read( );
  }
}

// =============================== ReadThermalIso ==========================

void cMaterial :: ReadThermalIso(void)
{
  // Read the number of thermal isotropic materials.

  int n;
  if (!(in >> n))
  {
    cout << "Error in the input of the number of isotropic materials!\n";
    exit(0);
  }

  // Read each material.

  int id;
  for (int i = 0; i < n; i++)
  {
    if (!(in >> id) || id < 1 || id > NumMat)
    {
      cout << "Error in the input of thermal isotropic material " << i+1 << "!\n";
      exit(0);
    }
    VecMat[id-1].SetThermProp(new cMatThermalIso( ));
    VecMat[id-1].ThermProp->Read( );
  }
}

// ============================== ReadThermalOrtho =========================

void cMaterial :: ReadThermalOrtho(void)
{
  // Read the number of thermal orthotropic materials.

  int n;
  if (!(in >> n))
  {
    cout << "Error in the input of the number of orthotropic materials!\n";
    exit(0);
  }

  // Read each material.

  int id;
  for (int i = 0; i < n; i++)
  {
    if (!(in >> id) || id < 1 || id > NumMat)
    {
      cout << "Error in the input of thermal orthotropic material " << i+1 << "!\n";
      exit(0);
    }
    VecMat[id-1].SetThermProp(new cMatThermalOrtho( ));
    VecMat[id-1].ThermProp->Read( );
  }
}


// ================================ Destroy ================================

void cMaterial :: Destroy(void)
{
  // Release the array of materials.

  delete []VecMat;
}

// =============================== ReadDensity =============================

void cMaterial :: ReadDensity(void)
{
  int n;
  if (!(in >> n))
  {
    cout << "\n Error on reading number of densities !!!\n\n";
    exit(0);
  }
  int label;
  double dens;
  for (int i = 0; i < n; i++)
  {
    if (!(in >> label) || !(in >> dens))
    {
      cout << "\n Error on reading densities " << i+1 << "!!!\n\n";
      exit(0);
    }

    cMaterial *pcMat = GetMaterial(label);
    if (!pcMat)
    {
      cout << "\n Error on reading density " << i+1 << "!!!\n\n";
      exit(0);
    }
    pcMat->Density.Resize(1);
    pcMat->Density = dens;
 }
}

// ============================= ReadThermExpan ============================

void cMaterial :: ReadThermExpan(void)
{
  // Read the number of input thermal expansion coeficient.

  int n;
  if (!(in >> n))
  {
    cout << "\n Error on reading number of thermal expansion coeficient!!!\n\n";
    exit(0);
  }

  // Read each coeficient.

  int id;
  for (int i = 0; i < n; i++)
  {
    if (!(in >> id) || id < 1 || id > NumMat)
    {
      cout << "\n Error on reading coeficient of thermal expansion " << i+1 << "!!!\n\n";
      exit(0);
    }
    VecMat[id-1].MecProp->ReadAlpha( );
  }
}

// ============================== GetMaterial ==============================

cMaterial *cMaterial :: GetMaterial(int id)
{
  if (id > 0 && id <= NumMat) return(&VecMat[id-1]);
  return(0);
}

// ============================== cMaterial ================================

cMaterial :: cMaterial(void)
{
  Density   = 0.0;
  MecProp   = 0;
  ThermProp = 0;
}

// ============================== ~cMaterial ===============================

cMaterial :: ~cMaterial(void)
{
  delete MecProp;
  delete ThermProp;
}

// -------------------------------------------------------------------------
// Class cMatElastIso:
// -------------------------------------------------------------------------

// ============================= cMatElastIso ==============================

cMatElastIso :: cMatElastIso( )
{
  Type = MAT_ELASTIC_ISOTROPIC;
  E  = 1.0;
  Nu = 0.0;
  Alpha = 0.0;
}

// ============================= ~cMatElastIso =============================

cMatElastIso :: ~cMatElastIso(void)
{
}

// ================================== Read =================================

void cMatElastIso :: Read(void)
{
  if (!(in >> E) || !(in >> Nu))
  {
    cout << "Error in the input of material " << mat->GetLabel( ) << "!\n";
    exit(0);
  }
}

// ================================== ReadAlpha ============================

void cMatElastIso :: ReadAlpha(void)
{
  if (!(in >> Alpha))
  {
    cout << "Error in the input of thermal expansion coeficient of material " << mat->GetLabel( ) << "!\n";
    exit(0);
  }
}

// ================================ GetParam ===============================

void cMatElastIso :: GetParam(double *param)
{
  param[0] = E;
  param[1] = Nu;
  param[2] = Alpha;
}

// -------------------------------------------------------------------------
// Class cMatElastOrtho:
// -------------------------------------------------------------------------

// =============================== cMatOrtho ===============================

cMatElastOrtho :: cMatElastOrtho( )
{
  Type = MAT_ELASTIC_ORTHOTROPIC;
  E1   = 1.0;
  E2   = 1.0;
  E3   = 1.0;
  Nu12 = 0.0;
  Nu13 = 0.0;
  Nu23 = 0.0;
  G12  = 0.0;
  G13  = 0.0;
  G23  = 0.0;
  Alpha1 = 0.0;
  Alpha2 = 0.0;
}

// ============================== ~cMatOrtho ===============================

cMatElastOrtho :: ~cMatElastOrtho(void)
{
}

// ================================== Read =================================

void cMatElastOrtho :: Read(void)
{
   if (!(in >> E1)   || !(in >> E2)   || !(in >> E3)   ||
       !(in >> Nu12) || !(in >> Nu13) || !(in >> Nu23) ||
       !(in >> G12)  || !(in >> G13)  || !(in >> G23))
{
    cout << "Error in the input of material " << mat->GetLabel( ) << "!\n";
    exit(0);
  }
}

// ================================== ReadAlpha ============================

void cMatElastOrtho :: ReadAlpha(void)
{
  if (!(in >> Alpha1) || !(in >> Alpha2))
  {
    cout << "Error in the input of thermal expansion coeficients of material " << mat->GetLabel( ) << "!\n";
    exit(0);
  }
}

// ================================ GetParam ===============================

void cMatElastOrtho :: GetParam(double *param)
{
  param[0]  = E1;
  param[1]  = E2;
  param[2]  = E3;
  param[3]  = Nu12;
  param[4]  = Nu13;
  param[5]  = Nu23;
  param[6]  = G12;
  param[7]  = G13;
  param[8]  = G23;
  param[9]  = Alpha1;
  param[10] = Alpha2;
}


// -------------------------------------------------------------------------
// Class cMatElastNonlin:
// -------------------------------------------------------------------------

// ============================= cMatElastNonlin ===========================

cMatElastNonlin :: cMatElastNonlin( )
{
  Type = MAT_ELASTIC_NONLINEAR;
  Nu = 0.0;
  NumPnt = 0;
}

// ============================= ~cMatElastNonlin ==========================

cMatElastNonlin :: ~cMatElastNonlin(void)
{
  delete []Eps;
  delete []Sig;
}

// ================================== Read =================================

void cMatElastNonlin :: Read(void)
{
  if (!(in >> NumPnt) || !(in >> Nu))
  {
    cout << "Error in the input of material " << mat->GetLabel( ) << "!\n";
    exit(0);
  }

  // Read the strain x stress table.

  Eps = new double[NumPnt];
  Sig = new double[NumPnt];
  for (int i = 0; i < NumPnt; i++)
  {
    if (!(in >> Eps[i]) || !(in >> Sig[i]))
    {
      cout << "Error in the input of the data of material " << mat->GetLabel( ) << "!\n";
      exit(0);
    }
  }
}

// ================================ GetParam ===============================

void cMatElastNonlin :: GetParam(double *param)
{
  param[0] = Nu;
  for (int i = 0; i < NumPnt; i++)
  {
    param[2*i+1] = Eps[i];
    param[2*i+2] = Sig[i];
  }
}


// -------------------------------------------------------------------------
// Class cMatFGM:
// -------------------------------------------------------------------------

// =============================== cMatFGM =================================

cMatFGM :: cMatFGM( )
{
  Type = MAT_FGM;
  E1   = 1.0;
  Nu1  = 0.0;
  E2   = 1.0;
  Nu2  = 0.0;
  Alpha1 = 1.0;
  Alpha2 = 1.0;
}

// ================================ ~cMatFGM ===============================

cMatFGM :: ~cMatFGM(void)
{
}

// ================================== Read =================================

void cMatFGM :: Read(void)
{
  GetDensity( ).Resize(2);
  if (!(in >> E1) || !(in >> Nu1) || !(in >> GetDensity( )[0]) ||
      !(in >> E2) || !(in >> Nu2) || !(in >> GetDensity( )[1]) || !(in >> HomMethod))
  {
    cout << "Error in the input of material " << mat->GetLabel( ) << "!\n";
    exit(0);
  }

  if (HomMethod < 1 || HomMethod > 2)
  {
    cout << "Error: invalid FGM homogenization method!\n";
    exit(0);
  }
}

// ================================== ReadAlpha ============================

void cMatFGM :: ReadAlpha(void)
{
  if (!(in >> Alpha1) || !(in >> Alpha2))
  {
    cout << "Error in the input of thermal expansion coeficients of material " << mat->GetLabel( ) << "!\n";
    exit(0);
  }
}

// ================================ GetParam ===============================

void cMatFGM :: GetParam(double *param)
{
  param[0] = E1;
  param[1] = Nu1;
  param[2] = E2;
  param[3] = Nu2;
  param[4] = HomMethod;
  param[5] = Alpha1;
  param[6] = Alpha2;
}


// -------------------------------------------------------------------------
// Class cMatFGMTTO:
// -------------------------------------------------------------------------

// ============================== cMatFGMTTO ===============================

cMatFGMTTO :: cMatFGMTTO( )
{
  Type = MAT_FGM_TTO;
  E1   = 1.0;
  Nu1  = 0.0;
  E2   = 1.0;
  Nu2  = 0.0;
  SigYm = 1.0;
  Hm   = 1.0;
  qtto = 1.0;
}

// ============================== ~cMatFGMTTO ==============================

cMatFGMTTO :: ~cMatFGMTTO(void)
{
}

// ================================== Read =================================

void cMatFGMTTO :: Read(void)
{
  GetDensity( ).Resize(2);
  if (!(in >> E1) || !(in >> Nu1) || !(in >> GetDensity( )[0]) || !(in >> E2) || !(in >> Nu2) ||
      !(in >> GetDensity( )[1]) || !(in >> HomMethod) || !(in >> SigYm) || !(in >> Hm) || !(in >> qtto))
  {
    cout << "Error in the input of material " << mat->GetLabel( ) << endl;
    exit(0);
  }

  if (HomMethod < 1 || HomMethod > 2)
  {
    cout << "Error: invalid FGM homogenization method!\n";
    exit(0);
  }
}

// ================================ GetParam ===============================

void cMatFGMTTO :: GetParam(double *param)
{
  param[0] = E1;
  param[1] = Nu1;
  param[2] = E2;
  param[3] = Nu2;
  param[4] = HomMethod;
  param[5] = SigYm;
  param[6] = Hm;
  param[7] = qtto;
}


// -------------------------------------------------------------------------
// Class cMatConcEuroCEB:
// -------------------------------------------------------------------------

// =========================== cMatConcEuroCEB =============================

cMatConcEuroCEB :: cMatConcEuroCEB( )
{
  Type = MAT_CONCRETE_EUROCEB;
  fcm = 1.0;
  ec1 = 0.0;
  Ecm = 1.0;
  fct = 0.0;
  Eci = 1.0;
  ey  = 0.0;
  Es  = 0.0;
  Rho = 0.0;
}

// =========================== ~cMatConcEuroCEB ============================

cMatConcEuroCEB :: ~cMatConcEuroCEB(void)
{
}

// ================================== Read =================================

void cMatConcEuroCEB :: Read(void)
{
  if (!(in >> fcm) || !(in >> ec1) || !(in >> Ecm) || !(in >> fct) ||
      !(in >> Eci) || !(in >> ey) || !(in >> Es) || !(in >> Rho))
  {
    cout << "Error in the input of material " << mat->GetLabel( ) << "!\n";
    exit(0);
  }

  if (fcm < 0.0 || ec1 >= 0.0 || Ecm <= 0.0 || Eci <= 0.0 || fct < 0.0 ||
      ey  < 0.0 || Es  < 0.0  || Rho < 0.0)
  {
    cout << "Invalid parameters of CEB material " << mat->GetLabel( ) << "!\n";
    cout << "fcm = " << fcm << " ec1 = " << ec1 << " Ecm = " << Ecm << " fct = "<< fct;
    cout << "Eci = " << Eci << " ey = " << ey << " Es = " << Es << " Rho = " << Rho << "\n";
    exit(0);
  }
}

// ================================= GetParam ===============================

void cMatConcEuroCEB :: GetParam(double *param)
{
  param[0] = fcm;
  param[1] = ec1;
  param[2] = Ecm;
  param[3] = fct;
  param[4] = Eci;
  param[5] = ey;
  param[6] = Es;
  param[7] = Rho;
}


// -------------------------------------------------------------------------
// Class cMatConcNBRCEB:
// -------------------------------------------------------------------------

// =========================== cMatConcNBRCEB ==============================

cMatConcNBRCEB :: cMatConcNBRCEB( )
{
  Type = MAT_CONCRETE_NBRCEB;
  fc  = 1.0;
  fct = 0.0;
  Eci = 1.0;
  ey  = 0.0;
  Es  = 0.0;
  Rho = 0.0;
}

// =========================== ~cMatConcNBRCEB =============================

cMatConcNBRCEB :: ~cMatConcNBRCEB(void)
{
}

// ================================== Read =================================

void cMatConcNBRCEB :: Read(void)
{
  if (!(in >> fc) || !(in >> fct) || !(in >> Eci) || !(in >> ey) ||
      !(in >> Es) || !(in >> Rho))
  {
    cout << "Error in the input of material " << mat->GetLabel( ) << "!\n";
    exit(0);
  }

  if (fc < 0.0 || Eci <= 0.0 || fct < 0.0 || ey < 0.0 || Es < 0.0 || Rho < 0.0)
  {
    cout << "Invalid parameters of CEB material " << mat->GetLabel( ) << "!\n";
    cout << "fct = " << fc << " fct = " << fct;
    cout << "Eci = " << Eci << " ey = " << ey << " Es = " << Es << " Rho = ";
    cout << Rho << "\n";
    exit(0);
  }
}

// ================================ GetParam ===============================

void cMatConcNBRCEB :: GetParam(double *param)
{
  param[0] = fc;
  param[1] = fct;
  param[2] = Eci;
  param[3] = ey;
  param[4] = Es;
  param[5] = Rho;
}


// -------------------------------------------------------------------------
// Class cMatViscoElastic:
// -------------------------------------------------------------------------

// ============================= cMatViscoElastic ==========================

cMatViscoElastic :: cMatViscoElastic( )
{
  Type = MAT_VISCOELASTIC;
  Nu = 0.0;
  RelMod.Eoo = 1.0;
  RelMod.n = 0;
}

// ============================ ~cMatViscoElastic ==========================

cMatViscoElastic :: ~cMatViscoElastic(void)
{
  delete []RelMod.E;
  delete []RelMod.rho;
}

// ================================== Read =================================

void cMatViscoElastic :: Read(void)
{
  if (!(in >> RelMod.Eoo) || !(in >> Nu) || !(in >> RelMod.n))
  {
    cout << "Error in the input of material " << mat->GetLabel( ) << "!\n";
    exit(0);
  }

  // Read the Prony series.

  RelMod.E   = new double[RelMod.n];
  RelMod.rho = new double[RelMod.n];
  for (int i = 0; i < RelMod.n; i++)
  {
    if (!(in >> RelMod.E[i]) || !(in >> RelMod.rho[i]))
    {
      cout << "Error in the input of the Prony series of material " << mat->GetLabel( );
      cout << "!\n,";
      exit(0);
    }
  }
}

// ================================ GetParam ===============================

void cMatViscoElastic :: GetParam(double *param)
{
  param[0] = RelMod.Eoo;
  param[1] = Nu;
  for (int i = 0; i < RelMod.n; i++)
  {
    param[2*i+2] = RelMod.E[i];
    param[2*i+3] = RelMod.rho[i];
  }
}


// -------------------------------------------------------------------------
// Class cMatPlastIsoLinHard:
// -------------------------------------------------------------------------

// ========================== cMatPlastIsoLinHard ==========================

cMatPlastIsoLinHard :: cMatPlastIsoLinHard( )
{
  Type = MAT_PLASTIC_ISOLINHARD;
  E    = 1.0;
  Nu   = 0.0;
  SigY0= 0.0;
  H    = 0.0;
}

// ========================== ~cMatPlastIsoLinHard =========================

cMatPlastIsoLinHard :: ~cMatPlastIsoLinHard(void)
{
}

// ================================== Read =================================

void cMatPlastIsoLinHard :: Read(void)
{
  if (!(in >> E) || !(in >> Nu) || !(in >> SigY0) || !(in >> H))
  {
    cout << "Error in the input of material " << mat->GetLabel( ) << "!\n";
    exit(0);
  }
}

// ================================ GetParam ===============================

void cMatPlastIsoLinHard :: GetParam(double *param)
{
  param[0] = E;
  param[1] = Nu;
  param[2] = SigY0;
  param[3] = H;
}


// -------------------------------------------------------------------------
// Class cMatPlastLinHard:
// -------------------------------------------------------------------------

// ========================== cMatPlastLinHard ==========================

cMatPlastLinHard :: cMatPlastLinHard( )
{
  Type = MAT_PLASTIC_LINHARD;
  E    = 0.0;
  Nu   = 0.0;
  SigY = 0.0;
  K    = 0.0;
  H    = 0.0;
}

// ========================= ~cMatPlastLinHard ==========================

cMatPlastLinHard :: ~cMatPlastLinHard(void)
{
}

// ================================== Read =================================

void cMatPlastLinHard :: Read(void)
{
  if (!(in >> E) || !(in >> SigY) || !(in >> K) || !(in >> H) || !(in >> Nu))
  {
    cout << "Error in the input of material " << mat->GetLabel( ) << "!\n";
    exit(0);
  }
//  cout << "E = " << E << endl;
//  cout << "SigY = " << SigY << endl;
//  cout << "K = " << K << endl;
//  cout << "H = " << H << endl;
//  cout << "Nu = " << Nu << endl;
}

// ================================ GetParam ===============================

void cMatPlastLinHard :: GetParam(double *param)
{
  param[0] = E;
  param[1] = Nu;
  param[2] = SigY;
  param[3] = K;
  param[4] = H;
}


// -------------------------------------------------------------------------
// Class cMatPlastIsoNonlinHard:
// -------------------------------------------------------------------------

// ========================== cMatPlastIsoNonlinHard =======================

cMatPlastIsoNonlinHard :: cMatPlastIsoNonlinHard( )
{
  Type = MAT_PLASTIC_ISONONLINHARD;
  E    = 0.0;
  SigY = 0.0;
  K    = 0.0;
  ExpH = 0.0;
  Nu   = 0.0;
}

// =========================== ~cMatPlastIsoNonlinHard =====================

cMatPlastIsoNonlinHard :: ~cMatPlastIsoNonlinHard(void)
{
}

// ================================== Read =================================

void cMatPlastIsoNonlinHard :: Read(void)
{
  if (!(in >> E) || !(in >> SigY) || !(in >> K) || !(in >> ExpH) || !(in >> Nu))
  {
    cout << "Error in the input of material " << mat->GetLabel( ) << "!\n";
    exit(0);
  }
}

// ================================ GetParam ===============================

void cMatPlastIsoNonlinHard :: GetParam(double *param)
{
  param[0] = E;
  param[1] = SigY;
  param[2] = K;
  param[3] = ExpH;
  param[4] = Nu;
}


// -------------------------------------------------------------------------
// Class cMatPlastIsoExpHard:
// -------------------------------------------------------------------------

// ========================== cMatPlastIsoExpHard ==========================

cMatPlastIsoExpHard :: cMatPlastIsoExpHard( )
{
  Type = MAT_PLASTIC_ISOEXPHARD;
  E    = 0.0;
  SigY = 0.0;
  Mu   = 0.0;
  Nu   = 0.0;
}

// ============================= ~cMatUniaxialElastoPlastic =================

cMatPlastIsoExpHard :: ~cMatPlastIsoExpHard(void)
{
}

// ================================== Read =================================

void cMatPlastIsoExpHard :: Read(void)
{
  if (!(in >> E) || !(in >> Nu) || !(in >> SigY) || !(in >> Mu))
  {
    cout << "Error in the input of material " << mat->GetLabel( ) << "!\n";
    exit(0);
  }
}

// ================================ GetParam ===============================

void cMatPlastIsoExpHard :: GetParam(double *param)
{
  param[0] = E;
  param[1] = SigY;
  param[2] = Mu;
  param[3] = Nu;
}


// -------------------------------------------------------------------------
// Class cMatDamageMazars:
// -------------------------------------------------------------------------

// ============================ cMatDamageMazars ==========================

cMatDamageMazars :: cMatDamageMazars( )
{
  Type = MAT_DAMAGE_MAZARS;
  E    = 0.0;
  Ac   = 0.0;
  Bc   = 0.0;
  At   = 0.0;
  Bt   = 0.0;
  epsd0= 0.0;
  Nu   = 0.0;
}

// =========================== ~cMatDamageMazars ==========================

cMatDamageMazars :: ~cMatDamageMazars(void)
{
}

// ================================== Read =================================

void cMatDamageMazars :: Read(void)
{
  if (!(in >> E) || !(in >> Ac) || !(in >> Bc) || !(in >> At) || !(in >> Bt) ||
      !(in >> epsd0) || !(in >> Nu))
  {
    cout << "Error in the input of material " << mat->GetLabel( ) << "!\n";
    exit(0);
  }
}

// ================================ GetParam ===============================

void cMatDamageMazars :: GetParam(double *param)
{
  param[0] = E;
  param[1] = Ac;
  param[2] = Bc;
  param[3] = At;
  param[4] = Bt;
  param[5] = epsd0;
  param[6] = Nu;
}

// -------------------------------------------------------------------------
// Class cMatDamageMuModel:
// -------------------------------------------------------------------------

// ============================ cMatDamageMuModel ==========================

cMatDamageMuModel :: cMatDamageMuModel( )
{
  Type = MAT_DAMAGE_MU_MODEL;
  E    = 0.0;
  Ac   = 0.0;
  Bc   = 0.0;
  At   = 0.0;
  Bt   = 0.0;
  epsdt0= 0.0;
  epsdc0= 0.0;
  Nu   = 0.0;
}

// =========================== ~cMatDamageMuModel ==========================

cMatDamageMuModel :: ~cMatDamageMuModel(void)
{
}

// ================================== Read =================================

void cMatDamageMuModel :: Read(void)
{
  if (!(in >> E) || !(in >> Ac) || !(in >> Bc) || !(in >> At) || !(in >> Bt) ||
      !(in >> epsdt0) || !(in >> epsdc0) || !(in >> Nu))
  {
    cout << "Error in the input of material " << mat->GetLabel( ) << "!\n";
    exit(0);
  }
}

// ================================ GetParam ===============================

void cMatDamageMuModel :: GetParam(double *param)
{
  param[0] = E;
  param[1] = Ac;
  param[2] = Bc;
  param[3] = At;
  param[4] = Bt;
  param[5] = epsdt0;
  param[6] = epsdc0;
  param[7] = Nu;
}


// -------------------------------------------------------------------------
// Class cMatLeeFenves:
// -------------------------------------------------------------------------

// =============================== cMatLeeFenves ===========================

cMatLeeFenves :: cMatLeeFenves( )
{
  Type = MAT_LEE_FENVES;
  E    = 0.0;
  fb0  = 0.0;
  fc0  = 0.0;
  ft0  = 0.0;
  a_c  = 0.0;
  b_c  = 0.0;
  d_c  = 0.0;
  a_t  = 0.0;
  b_t  = 0.0;
  d_t  = 0.0;
  Gt   = 0.0;
  Gc   = 0.0;
  lt   = 0.0;
  lc   = 0.0;
  s0   = 0.0;
}

// ============================== ~cMatLeeFenves ===========================

cMatLeeFenves :: ~cMatLeeFenves(void)
{
}

// ================================== Read =================================

void cMatLeeFenves :: Read(void)
{
  if (!(in >> E) || !(in >> fb0) || !(in >> fc0) || !(in >> ft0)   ||
      !(in >> a_c) || !(in >> b_c) || !(in >> d_c) || !(in >> a_t) ||
      !(in >> b_t) || !(in >> d_t) || !(in >> Gt) || !(in >> Gc)   ||
      !(in >> lt) || !(in >> lc) || !(in >> s0))
  {
    cout << "Error in the input of material " << mat->GetLabel( ) << "!\n";
    exit(0);
  }
}

// ================================ GetParam ===============================

void cMatLeeFenves :: GetParam(double *param)
{
  param[0] = E;
  param[1] = fb0;
  param[2] = fc0;
  param[3] = ft0;
  param[4] = a_c;
  param[5] = b_c;
  param[6] = d_c;
  param[7] = a_t;
  param[8] = b_t;
  param[9] = d_t;
  param[10] = Gt;
  param[11] = Gc;
  param[12] = lt;
  param[13] = lc;
  param[14] = s0;
}


// -------------------------------------------------------------------------
// Class cMatProgFailHashin:
// -------------------------------------------------------------------------

// =========================== cMatProgFailHashin ==========================

cMatProgFailHashin :: cMatProgFailHashin( )
{
  Type = MAT_PROGFAIL_HASHIN;
  E1   = 1.0;
  E2   = 1.0;
  E3   = 1.0;
  Nu12 = 0.0;
  Nu13 = 0.0;
  Nu23 = 0.0;
  G12  = 0.0;
  G13  = 0.0;
  G23  = 0.0;
  Xt   = 0.0;
  Xc   = 0.0;
  Yt   = 0.0;
  Yc   = 0.0;
  Zt   = 0.0;
  Zc   = 0.0;
  S12  = 0.0;
  S13  = 0.0;
  S23  = 0.0;
  Alpha = 0.0;
}

// ========================== ~cMatProgFailHashin ==========================

cMatProgFailHashin :: ~cMatProgFailHashin(void)
{
}

// ================================== Read =================================

void cMatProgFailHashin :: Read(void)
{
    if (!(in >> E1)   || !(in >> E2)   || !(in >> E3)   ||
        !(in >> Nu12) || !(in >> Nu13) || !(in >> Nu23) ||
        !(in >> G12)  || !(in >> G13)  || !(in >> G23)  ||
        !(in >> Xt)   || !(in >> Xc)   || !(in >> Yt)   ||
        !(in >> Yc)   || !(in >> Zt)   || !(in >> Zc)   ||
        !(in >> S12)  || !(in >> S13)  || !(in >> S23)  ||
        !(in >> Alpha))
  {
    cout << "Error in the input of material " << mat->GetLabel( ) << "!\n";
    exit(0);
  }
}

// ================================ GetParam ===============================

void cMatProgFailHashin :: GetParam(double *param)
{
  param[0] = E1;
  param[1] = E2;
  param[2] = E3;
  param[3] = Nu12;
  param[4] = Nu13;
  param[5] = Nu23;
  param[6] = G12;
  param[7] = G13;
  param[8] = G23;
  param[9] = Xt;
  param[10] = Xc;
  param[11] = Yt;
  param[12] = Yc;
  param[13] = Zt;
  param[14] = Zc;
  param[15] = S12;
  param[16] = S13;
  param[17] = S23;
  param[18] = Alpha;
}


// -------------------------------------------------------------------------
// Class cMatProgFailTsai:
// -------------------------------------------------------------------------

// =========================== cMatProgFailTsai ==========================

cMatProgFailTsai :: cMatProgFailTsai( )
{
  Type = MAT_PROGFAIL_TSAI;
  E1   = 1.0;
  E2   = 1.0;
  E3   = 1.0;
  Nu12 = 0.0;
  Nu13 = 0.0;
  Nu23 = 0.0;
  G12  = 0.0;
  G13  = 0.0;
  G23  = 0.0;
  Xt   = 0.0;
  Xc   = 0.0;
  Yt   = 0.0;
  Yc   = 0.0;
  Zt   = 0.0;
  Zc   = 0.0;
  S12  = 0.0;
  S13  = 0.0;
  S23  = 0.0;
}

// ========================== ~cMatProgFailTsai =========================

cMatProgFailTsai :: ~cMatProgFailTsai(void)
{
}

// ================================ Read ================================

void cMatProgFailTsai :: Read(void)
{
  if   (!(in >> E1)   || !(in >> E2)   || !(in >> E3)   ||
        !(in >> Nu12) || !(in >> Nu13) || !(in >> Nu23) ||
        !(in >> G12)  || !(in >> G13)  || !(in >> G23)  ||
        !(in >> Xt)   || !(in >> Xc)   || !(in >> Yt)   ||
        !(in >> Yc)   || !(in >> Zt)   || !(in >> Zc)   ||
        !(in >> S12)  || !(in >> S13)  || !(in >> S23))
  {
    cout << "Error in the input of material " << mat->GetLabel( ) << "!\n";
    exit(0);
  }
}

// ============================= GetParam ===============================

void cMatProgFailTsai :: GetParam(double *param)
{
  param[0] = E1;
  param[1] = E2;
  param[2] = E3;
  param[3] = Nu12;
  param[4] = Nu13;
  param[5] = Nu23;
  param[6] = G12;
  param[7] = G13;
  param[8] = G23;
  param[9] = Xt;
  param[10] = Xc;
  param[11] = Yt;
  param[12] = Yc;
  param[13] = Zt;
  param[14] = Zc;
  param[15] = S12;
  param[16] = S13;
  param[17] = S23;
}


// -------------------------------------------------------------------------
// Class cMatProgFailEngelstad:
// -------------------------------------------------------------------------

// ======================= cMatProgFailEngelstad ===========================

cMatProgFailEngelstad :: cMatProgFailEngelstad( )
{
  Type = MAT_PROGFAIL_ENGELSTAD;
  E1   = 1.0;
  E2   = 1.0;
  E3   = 1.0;
  Nu12 = 0.0;
  Nu13 = 0.0;
  Nu23 = 0.0;
  G12  = 0.0;
  G13  = 0.0;
  G23  = 0.0;
  Xt   = 0.0;
  Xc   = 0.0;
  Yt   = 0.0;
  Yc   = 0.0;
  Zt   = 0.0;
  Zc   = 0.0;
  S12  = 0.0;
  S13  = 0.0;
  S23  = 0.0;
  Alpha = 1e-6;
}

// ===================== ~cMatProgFailEngelstad ===========================

cMatProgFailEngelstad :: ~cMatProgFailEngelstad(void)
{
}

// ============================== Read ====================================

void cMatProgFailEngelstad :: Read(void)
{
  if (!(in >> E1)   || !(in >> E2)   || !(in >> E3)   ||
      !(in >> Nu12) || !(in >> Nu13) || !(in >> Nu23) ||
      !(in >> G12)  || !(in >> G13)  || !(in >> G23)  ||
      !(in >> Xt)   || !(in >> Xc)   || !(in >> Yt)   ||
      !(in >> Yc)   || !(in >> Zt)   || !(in >> Zc)   ||
      !(in >> S12)  || !(in >> S13)  || !(in >> S23)  ||
      !(in >> Alpha))
  {
    cout << "Error in the input of material " << mat->GetLabel( ) << "!\n";
    exit(0);
  }
}

// ============================ GetParam ==================================

void cMatProgFailEngelstad :: GetParam(double *param)
{
  param[0] = E1;
  param[1] = E2;
  param[2] = E3;
  param[3] = Nu12;
  param[4] = Nu13;
  param[5] = Nu23;
  param[6] = G12;
  param[7] = G13;
  param[8] = G23;
  param[9] = Xt;
  param[10] = Xc;
  param[11] = Yt;
  param[12] = Yc;
  param[13] = Zt;
  param[14] = Zc;
  param[15] = S12;
  param[16] = S13;
  param[17] = S23;
  param[18] = Alpha;
}


// -------------------------------------------------------------------------
// Class cMatCZMModeIBilinear:
// -------------------------------------------------------------------------

// ========================= cMatCZMModeIBilinear ==========================

cMatCZMModeIBilinear :: cMatCZMModeIBilinear( )
{
  Type = MAT_CZM_MODEI_BILINEAR;
  Ft  = 0.0;
  GIc = 0.0;
  K0  = 0.0;
}

// ========================= cMatCZMModeIBilinear ==========================

cMatCZMModeIBilinear :: ~cMatCZMModeIBilinear(void)
{
}

// ================================== Read =================================

void cMatCZMModeIBilinear :: Read(void)
{
  if (!(in >> Ft) || !(in >> GIc) || !(in >> K0))
  {
    cout << "Error in the input of material " << mat->GetLabel( ) << "!\n";
    exit(0);
  }
}

// ================================ GetParam ===============================

void cMatCZMModeIBilinear :: GetParam(double *param)
{
  param[0] = Ft;
  param[1] = GIc;
  param[2] = K0;
}


// -------------------------------------------------------------------------
// Class cMatCZMModeIExponencial:
// -------------------------------------------------------------------------

// ======================== cMatCZMModeIExponential ========================

cMatCZMModeIExponential :: cMatCZMModeIExponential( )
{
  Type = MAT_CZM_MODEI_EXPONENTIAL;
  Ft  = 0.0;
  GIc = 0.0;
}

// ========================= ~cMatCZMModeIBilinear =========================

cMatCZMModeIExponential :: ~cMatCZMModeIExponential(void)
{
}

// ================================== Read =================================

void cMatCZMModeIExponential :: Read(void)
{
  if (!(in >> Ft) || !(in >> GIc))
  {
    cout << "Error in the input of material " << mat->GetLabel( ) << "!\n";
    exit(0);
  }
}

// ================================ GetParam ===============================

void cMatCZMModeIExponential :: GetParam(double *param)
{
  param[0] = Ft;
  param[1] = GIc;
}


// -------------------------------------------------------------------------
// Class cMatThermalIso:
// -------------------------------------------------------------------------

// ============================= cMatThermalIso ============================

cMatThermalIso :: cMatThermalIso( )
{
  Type = MAT_THERMAL_ISOTROPIC;
  k    = 1.0;
}

// ============================= ~cMatThermalIso ===========================

cMatThermalIso :: ~cMatThermalIso(void)
{
}

// ================================== Read =================================

void cMatThermalIso :: Read(void)
{
  if (!(in >> k))
  {
    cout << "Error in the input of material " << mat->GetLabel( ) << "!\n";
    exit(0);
  }
}

// ================================ GetParam ===============================

void cMatThermalIso :: GetParam(double *param)
{
  param[0] = k;
}


// -------------------------------------------------------------------------
// Class cMatThermalOrtho:
// -------------------------------------------------------------------------

// ============================= cMatThermalOrtho ==========================

cMatThermalOrtho :: cMatThermalOrtho( )
{
  Type = MAT_THERMAL_ORTHOTROPIC;
  k1   = 1.0;
  k2   = 1.0;
  k3   = 1.0;
}

// ============================== ~cMatOrtho ===============================

cMatThermalOrtho :: ~cMatThermalOrtho(void)
{
}

// ================================== Read =================================

void cMatThermalOrtho :: Read(void)
{
   if (!(in >> k1) || !(in >> k2) || !(in >> k3))
   {
     cout << "Error in the input of material " << mat->GetLabel( ) << "!\n";
     exit(0);
   }
}

// ================================ GetParam ===============================

void cMatThermalOrtho :: GetParam(double *param)
{
  param[0] = k1;
  param[1] = k2;
  param[2] = k3;
}

// ======================================================= End of file =====
