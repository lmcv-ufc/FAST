// -------------------------------------------------------------------------
// section.cpp - implementation of the Section class.
// -------------------------------------------------------------------------
// Created:      07-Mar-2013     Iuri Barcelos Rocha
//
// Modified:     06-Oct-2013     Evandro Parente Junior
//               Unification with SecBar class.
//
// Modified:     05-May-2014     Evandro Parente/Edson Dantas
//               Implementation of interface sections (2D/3D).
//
// Modified:     08-Apr-2015     Elias Saraiva Barroso
//               Input/output operations using c++ streams.
//
// Modified:     12-Sep-2018     Samir Parente Auad
//               Implementation of General Shell section.
//
// Modified:     16-Sep-2022     Pedro Ygor Rodrigues Mesquita
//               Implementation of HSDT Plate section.
// -------------------------------------------------------------------------

#include <iostream>
#include <math.h>
#include <stdlib.h>

using namespace std;

#include "section.h"
#include "secbar.h"
#include "material.h"
#include "mat.h"
#include "element.h"
#include "utl.h"
#include "gblvar.h"

// -------------------------------------------------------------------------
// Static variables:
//
int        cSection :: NumSec    = 0;
cSection** cSection :: VecSec    = 0;
int        cSection :: NumIntX   = 100;
int        cSection :: NumIntY   = 100;
eSecIntMet cSection :: IntMethod = SEC_DEFAULT;

// -------------------------------------------------------------------------
// Public methods:
//

// ============================== ReadNumSec ===============================

void cSection :: ReadNumSec(void)
{
  // Read the number of sections.

  if (!(in >> NumSec) || (NumSec <= 0))
  {
    cout << "Error in the input of the number of sections!\n";
    exit(0);
  }

  // Alloc the array of materials.

  VecSec = new cSection*[NumSec];
}

// ============================== ReadHomIso2D =============================

void cSection :: ReadHomIso2D(void)
{
  // Read the number of 2D homogeneous isotropic sections.

  int n;
  if (!(in >> n))
  {
    cout << "Error in the input of the number of 2D homogeneous sections!\n";
    exit(0);
  }

  // Read each section.

  int id;
  for (int i = 0; i < n; i++)
  {
    if (!(in >> id) || id < 1 || id > NumSec)
    {
      cout << "Error in the input of 2D homogeneous section " << i+1 << "!\n";
      exit(0);
    }
    cSecHomIso2D *sec = new cSecHomIso2D(id);
    sec->Read( );
  }
}

// ============================== ReadHomIso3D =============================

void cSection :: ReadHomIso3D(void)
{
  // Read the number of 3D homogeneous isotropic sections.

  int n;
  if (!(in >> n))
  {
    cout << "Error in the input of the number of 3D homogeneous sections!\n";
    exit(0);
  }

  // Read each section.

  int id;
  for (int i = 0; i < n; i++)
  {
    if (!(in >> id) || id < 1 || id > NumSec)
    {
      cout << "Error in the input of 3D homogeneous section " << i+1 << "!\n";
      exit(0);
    }
    cSecHomIso3D *sec = new cSecHomIso3D(id);
    sec->Read( );
  }
}

// ============================= ReadHomOrtho2D ============================

void cSection :: ReadHomOrtho2D(void)
{
  // Read the number of 2D homogeneous orthotropic sections.

  int n;
  if (!(in >> n))
  {
    cout << "Error in the input of the number of 2D homogeneous sections!\n";
    exit(0);
  }

  // Read each section.

  int id;
  for (int i = 0; i < n; i++)
  {
    if (!(in >> id) || id < 1 || id > NumSec)
    {
      cout << "Error in the input of 2D homogeneous section " << i+1 << "!\n";
      exit(0);
    }
    cSecHomOrtho2D *sec = new cSecHomOrtho2D(id);
    sec->Read( );
  }
}

// ============================= ReadHomOrtho3D ============================

void cSection :: ReadHomOrtho3D(void)
{
  // Read the number of 3D homogeneous orthotropic sections.

  int n;
  if (!(in >> n))
  {
    cout << "Error in the input of the number of 3D homogeneous sections!\n";
    exit(0);
  }

  // Read each section.

  int id;
  for (int i = 0; i < n; i++)
  {
    if (!(in >> id) || id < 1 || id > NumSec)
    {
      cout << "Error in the input of 3D homogeneous section " << i+1 << "!\n";
      exit(0);
    }
    cSecHomOrtho3D *sec = new cSecHomOrtho3D(id);
    sec->Read( );
  }
}

// ================================ ReadFGM2D ==============================

void cSection :: ReadFGM2D(void)
{
  // Read the number of FGM 2D sections.

  int n;
  if (!(in >> n))
  {
    cout << "Error in the input of the number of sections!\n";
    exit(0);
  }

  // Read the section.

  int id;
  for (int i = 0; i < n; i++)
  {
    if (!(in >> id) || id < 1 || id > NumSec)
    {
      cout << "Error in the input of FGM 2D section " << i+1 << "!\n";
      exit(0);
    }
    cSecFGM2D *sec = new cSecFGM2D(id);
    sec->Read( );
  }
}

// ================================ ReadFGM3D ==============================

void cSection :: ReadFGM3D(void)
{
  // Read the number of FGM 3D sections.

  int n;
  if (!(in >> n))
  {
    cout << "Error in the input of the number of sections!\n";
    exit(0);
  }

  // Read the section.

  int id;
  for (int i = 0; i < n; i++)
  {
    if (!(in >> id) || id < 1 || id > NumSec)
    {
      cout << "Error in the input of FGM 3D section " << i+1 << "!\n";
      exit(0);
    }
    cSecFGM3D *sec = new cSecFGM3D(id);
    sec->Read( );
  }
}

// ============================== ReadHSDTPlate =============================

void cSection :: ReadHSDTPlate(void)
{
  // Read the number of HSDT plate sections.

  int n;
  if (!(in >> n))
  {
    cout << "Error in the input of the number of sections!\n";
    exit(0);
  }

  // Read the section.

  int id;
  for (int i = 0; i < n; i++)
  {
    if (!(in >> id) || id < 1 || id > NumSec)
    {
      cout << "Error in the input of HSDT plate section " << i+1 << "!\n";
      exit(0);
    }
    cSecHSDTPlate *sec = new cSecHSDTPlate(id);
    sec->Read( );
  }
}

// ============================== ReadHSDTFGMPlate =============================

void cSection :: ReadHSDTFGMPlate(void)
{
  // Read the number of HSDT FGM plate sections.

  int n;
  if (!(in >> n))
  {
    cout << "Error in the input of the number of sections!\n";
    exit(0);
  }

  // Read the section.

  int id;
  for (int i = 0; i < n; i++)
  {
    if (!(in >> id) || id < 1 || id > NumSec)
    {
      cout << "Error in the input of HSDT FGM plate section " << i+1 << "!\n";
      exit(0);
    }
    cSecHSDTFGMPlate *sec = new cSecHSDTFGMPlate(id);
    sec->Read( );
  }
}

// ============================== ReadFGMShell =============================

void cSection :: ReadFGMShell(void)
{
  // Read the number of FGM shell sections.

  int n;
  if (!(in >> n))
  {
    cout << "Error in the input of the number of sections!\n";
    exit(0);
  }

  // Read the section.

  int id;
  for (int i = 0; i < n; i++)
  {
    if (!(in >> id) || id < 1 || id > NumSec)
    {
      cout << "Error in the input of FGM shell section " << i+1 << "!\n";
      exit(0);
    }
    cSecFGMShell *sec = new cSecFGMShell(id);
    sec->Read( );
  }
}

// ============================== ReadGeneralShell =============================

void cSection :: ReadGeneralShell(void)
{
  // Read the number of General shell sections.

  int n;
  if (!(in >> n))
  {
    cout << "Error in the input of the number of sections!\n";
    exit(0);
  }

  // Read the section.

  int id;
  for (int i = 0; i < n; i++)
  {
    if (!(in >> id) || id < 1 || id > NumSec)
    {
      cout << "Error in the input of General shell section " << i+1 << "!\n";
      exit(0);
    }
    cSecGeneralShell *sec = new cSecGeneralShell(id);
    sec->Read( );
  }
}

// ============================ ReadLaminated2D ============================

void cSection :: ReadLaminated2D(void)
{
  // Read the number of laminated 2D sections.

  int n;
  if (!(in >> n))
  {
    cout << "Error in the input of the number of laminated sections!\n";
    exit(0);
  }

  // Read each section.

  int id;
  for (int i = 0; i < n; i++)
  {
    if (!(in >> id) || id < 1 || id > NumSec)
    {
      cout << "Error in the input of laminated 2D section " << i+1 << "!\n";
      exit(0);
    }
    cSecLaminated2D *sec = new cSecLaminated2D(id);
    sec->Read( );
  }
}

// =========================== ReadLaminatedShell ==========================

void cSection :: ReadLaminatedShell(void)
{
  // Read the number of laminated shell sections.

  int n;
  if (!(in >> n))
  {
    cout << "Error in the input of the number of laminated sections!\n";
    exit(0);
  }

  // Read each section.

  int id;
  for (int i = 0; i < n; i++)
  {
    if (!(in >> id) || id < 1 || id > NumSec)
    {
      cout << "Error in the input of laminated shell section " << i+1 << "!\n";
      exit(0);
    }
    cSecLaminatedShell *sec = new cSecLaminatedShell(id);
    sec->Read( );
  }
}

// ========================== ReadLaminatedSolid ===========================

void cSection :: ReadLaminatedSolid(void)
{
  // Read the number of laminated solid sections.

  int n;
  if (!(in >> n))
  {
    cout << "Error in the input of the number of laminated solid sections!\n";
    exit(0);
  }

  // Read each section.

  int id;
  for (int i = 0; i < n; i++)
  {
    if (!(in >> id) || id < 1 || id > NumSec)
    {
      cout << "Error in the input of laminated solid section " << i+1 << "!\n";
      exit(0);
    }
    cSecLaminatedSolid *sec = new cSecLaminatedSolid(id);
    sec->Read( );
  }
}

// ============================ ReadInterface2D ============================

void cSection :: ReadInterface2D(void)
{
  // Read the number of Interface2D sections.

  int n;
  if (!(in >> n))
  {
    cout << "Error in the input of the number of Interface 2D sections!\n";
    exit(0);
  }

  // Read each section.

  int id;
  for (int i = 0; i < n; i++)
  {
    if (!(in >> id) || id < 1 || id > NumSec)
    {
      cout << "Error in the input of Interface 2D section " << i+1 << "!\n";
      exit(0);
    }
    cSecInterface2D *sec = new cSecInterface2D(id);
    sec->Read( );
  }
}

// ============================ ReadInterface3D ============================

void cSection :: ReadInterface3D(void)
{
  // Read the number of Interface3D sections.

  int n;
  if (!(in >> n))
  {
    cout << "Error in the input of the number of Interface 3D sections!\n";
    exit(0);
  }

  // Read each section.

  int id;
  for (int i = 0; i < n; i++)
  {
    if (!(in >> id) || id < 1 || id > NumSec)
    {
      cout << "Error in the input of Interface 3D section " << i+1 << "!\n";
      exit(0);
    }
    cSecInterface3D *sec = new cSecInterface3D(id);
    sec->Read( );
  }
}

// ============================== ReadSecGen ===============================

void cSection :: ReadSecGen(void)
{
  // Read the number of generic cross sections.

  int n;
  if ((!(in >> n)) || (n < 1) || (n > NumSec))
  {
    cout << "Error in the input of the number of generic sections!\n";
    exit(0);
  }

  // Read each section.

  int id;
  cSection *sec;
  for (int i = 0; i < n; i++)
  {
    if ((!(in >> id)) || (id < 1) || (id > NumSec))
    {
      cout << "Error in the input of generic bar sections!\n";
      exit(0);
    }
    sec = new cSecBarGen(id);
    sec->Read( );
  }
}

// ========================== ReadSecB3DCoupled ============================

void cSection :: ReadSecB3DCoupled(void)
{
  // Read the number of Beam 3D coupled cross-sections.

  int n;
  if ((!(in >> n)) || (n < 1) || (n > NumSec))
  {
    cout << "Error in the input of the number of B3D coupled sections!\n";
    exit(0);
  }

  // Read each section.

  int id;
  cSection *sec;
  for (int i = 0; i < n; i++)
  {
    if ((!(in >> id)) || (id < 1) || (id > NumSec))
    {
      cout << "Error in the input of B3D coupled sections!\n";
      exit(0);
    }
    sec = new cSecBarB3DCoupled(id);
    sec->Read( );
  }
}

// ============================== ReadSecRect ==============================

void cSection :: ReadSecRect(void)
{
  // Read the number of rectangular cross sections.

  int n;
  if ((!(in >> n)) || (n < 1) || (n > NumSec))
  {
    cout << "Error in the input of the number of retangular sections!\n";
    exit(0);
  }

  // Read each section.

  int id;
  cSection *sec;
  for (int i = 0; i < n; i++)
  {
    if ((!(in >> id)) || (id < 1) || (id > NumSec))
    {
      cout << "Error in the input of rectangular bar sections!\n";
      exit(0);
    }
    sec = new cSecBarRect(id);
    sec->Read( );
  }
}

// ============================== ReadSecCirc ==============================

void cSection :: ReadSecCirc(void)
{
  // Read the number of rectangular cross sections.

  int n;
  if ((!(in >> n)) || (n < 1) || (n > NumSec))
  {
    cout << "Error in the input of the number of circular cross sections!\n";
    exit(0);
  }

  // Read each section.

  int id;
  cSection *sec;
  for (int i = 0; i < n; i++)
  {
    if ((!(in >> id)) || (id < 1) || (id > NumSec))
    {
      cout << "Error in the input of circular bar sections!\n";
      exit(0);
    }
    sec = new cSecBarCirc(id);
    sec->Read( );
  }
}

// ============================== ReadSecTube ==============================

void cSection :: ReadSecTube(void)
{
  int n;
  if ((!(in >> n)) || (n < 1) || (n > NumSec))
  {
    cout << "Erro in the input of the number of tube cross section!\n";
    exit(0);
  }

  // Read each section.

  int id;
  cSection *sec;
  for (int i = 0; i < n; i++)
  {
    if ((!(in >> id)) || (id < 1) || (id > NumSec))
    {
      cout << "Error in the input of tube bar sections!\n";
      exit(0);
    }
    sec = new cSecBarTube(id);
    sec->Read( );
  }
}

// ============================= ReadSecRCRect =============================

void cSection :: ReadSecRCRect(void)
{
  // Read the number of rectangular reinforced concrete cross sections.

  int n;
  if ((!(in >> n)) || (n < 1) || (n > NumSec))
  {
    cout << "Error in the input of the number of retangular RC sections!\n";
    exit(0);
  }

  // Read each section.

  int id;
  cSection *sec;
  for (int i = 0; i < n; i++)
  {
    if ((!(in >> id)) || (id < 1) || (id > NumSec))
    {
      cout << "Error in the input of rectangular RC cross sections!\n";
      exit(0);
    }
    sec = new cSecBarRCRect(id);
    sec->Read( );
  }
}

// ============================= ReadSecInteg ==============================

void cSection :: ReadSecInteg(void)
{
  // Read the section integration method.

  char label[100];
  Utl::ReadString(in, label);

  if ((string(label) == "Default") || (string(label) == "default"))
    IntMethod = SEC_DEFAULT;
  else if ((string(label) == "Fiber") || (string(label) == "fiber"))
    IntMethod = SEC_FIBER;
  else if ((string(label) == "Gauss") || (string(label) == "gauss"))
    IntMethod = SEC_GAUSS;
  else
  {
    cout << "Invalid section integration method - using default method!\n";
    IntMethod = SEC_DEFAULT;
  }

  // Read integration parameters.

  if (IntMethod != SEC_DEFAULT)
  {
    if (!(in >> NumIntX) || !(in >> NumIntY))
    {
      cout << "Error in the input of section integration data!\n";
      exit(0);
    }
  }
}

// ================================ Destroy ================================

void cSection :: Destroy(void)
{
  // Destroy each section.

  for (int i = 0; i < NumSec; i++) delete VecSec[i];

  // Release the array of sections.

  delete []VecSec;
}

// ============================== GetSection ===============================

cSection *cSection :: GetSection(int id)
{
  if (id > 0 && id <= NumSec) return(VecSec[id-1]);
  return(0);
}

// =============================== cSection ================================

cSection :: cSection(int id)
{
  Label = id;           // Store the section label
  VecSec[id-1] = this;  // Add the section to the section vector
}

// ============================== ~cSection ================================

cSection :: ~cSection(void)
{
}

// ========================= GetOrientationVector ==========================

void cSection :: GetOrientationVector(cVector &v)
{
  v[0] = 1.0;
  v[1] = 0.0;
  v[2] = 0.0;
}


// -------------------------------------------------------------------------
// Class cSecHomIso2D:
// -------------------------------------------------------------------------

// ============================ cSecHomIso2D ===============================

cSecHomIso2D :: cSecHomIso2D(int id) : cSection(id)
{
  Type = SEC_HOM_ISO_2D;
  Thk  = 1.0;
}

// =========================== ~cSecHomIso2D ===============================

cSecHomIso2D :: ~cSecHomIso2D(void)
{
}

// ================================== Read =================================

void cSecHomIso2D :: Read(void)
{
  int mat;
  if (!(in >> mat) || !(in >> Thk))
  {
    cout << "Error in the input of section " << Label << "!\n";
    exit(0);
  }

  Mat = cMaterial::GetMaterial(mat);
  if (!Mat)
  {
    cout << "Error in the input of section " << Label << "!\n";
    cout << "Invalid material!\n";
    exit(0);
  }
}


// -------------------------------------------------------------------------
// Class cSecHomIso3D:
// -------------------------------------------------------------------------

// ============================ cSecHomIso3D ===============================

cSecHomIso3D :: cSecHomIso3D(int id) : cSection(id)
{
  Type = SEC_HOM_ISO_3D;
}

// =========================== ~cSecHomIso2D ===============================

cSecHomIso3D :: ~cSecHomIso3D(void)
{
}

// ================================ Read ===================================

void cSecHomIso3D :: Read(void)
{
  int mat;
  if (!(in >> mat))
  {
    cout << "Error in the input of section " << Label << "!\n";
    exit(0);
  }

  Mat = cMaterial::GetMaterial(mat);
  if (!Mat)
  {
    cout << "Error in the input of section " << Label << "!\n";
    cout << "Invalid material!\n";
    exit(0);
  }
}


// -------------------------------------------------------------------------
// Class cSecHomOrtho2D:
// -------------------------------------------------------------------------

// =========================== cSecHomOrtho2D ==============================

cSecHomOrtho2D :: cSecHomOrtho2D(int id) : cSection(id)
{
  Type = SEC_HOM_ORTHO_2D;
  Thk  = 1.0;
}

// ========================== ~cSecHomOrtho2D ==============================

cSecHomOrtho2D :: ~cSecHomOrtho2D(void)
{
}

// ========================= GetOrientationVector ==========================

void cSecHomOrtho2D :: GetOrientationVector(cVector &v)
{
  v[0] = cos(Ang);
  v[1] = sin(Ang);
  v[2] = 0.0;
}

// ================================ Read ===================================

void cSecHomOrtho2D :: Read(void)
{
  int mat;
  if (!(in >> mat) || !(in >> Ang)  || !(in >> Thk))
  {
    cout << "Error in the input of section " << Label << "!\n";
    exit(0);
  }

  Mat = cMaterial::GetMaterial(mat);
  if (!Mat)
  {
    cout << "Error in the input of section " << Label << "!\n";
    cout << "Invalid material!\n";
    exit(0);
  }
}


// -------------------------------------------------------------------------
// Class cSecHomOrtho3D:
// -------------------------------------------------------------------------

// =========================== cSecHomOrtho3D ==============================

cSecHomOrtho3D :: cSecHomOrtho3D(int id) : cSection(id)
{
  Type = SEC_HOM_ORTHO_3D;

  OriVec.Resize(3);
  OriVec[0] = 1.0;
  OriVec[1] = 0.0;
  OriVec[2] = 0.0;
}

// ========================== ~cSecHomOrtho3D ==============================

cSecHomOrtho3D :: ~cSecHomOrtho3D(void)
{
}

// ================================ Read ===================================

void cSecHomOrtho3D :: Read(void)
{
  int mat;
  if (!(in >> mat))
  {
    cout << "Error in the input of section " << Label << "!\n";
    exit(0);
  }

  Mat = cMaterial::GetMaterial(mat);
  if (!Mat)
  {
    cout << "Error in the input of section " << Label << "!\n";
    cout << "Invalid material!\n";
    exit(0);
  }

  // Read the orientation vector.

  double a1,a2,a3;
  if (!(in >> a1) || !(in >> a2) || !(in >> a3))
  {
    cout << "Error in the input of section orientation vector\n";
    exit(0);
  }
  OriVec[0] = a1;
  OriVec[1] = a2;
  OriVec[2] = a3;
}


// -------------------------------------------------------------------------
// Class cSecFGM2D:
// -------------------------------------------------------------------------

// ============================== cSecFGM2D ================================

cSecFGM2D :: cSecFGM2D(int id) : cSection(id)
{
  Type = SEC_FGM_2D;
  Thk  = 1.0;
}

// ============================= ~cSecFGM2D ================================

cSecFGM2D :: ~cSecFGM2D(void)
{
}

// ================================== Read =================================

void cSecFGM2D :: Read(void)
{
  int mat,nnode;
  if (!(in >> mat) || !(in >> Thk) || !(in >> nnode))
  {
    cout << "Error in the input of section " << Label << "!\n"; 
    exit(0);
  }

  Mat = cMaterial::GetMaterial(mat);
  if (!Mat)
  {
    cout << "Error in the input of section " << Label << "!\n";
    cout << "Invalid material!\n";

    exit(0);
  }

  int node;
  double v2;
  for (int i = 0; i < nnode; i++)
  {	  
    if (!(in >> node) || !(in >> v2) || !cNode::FindNode(node))
    {
      cout << "Error in the input of section " << Label << "!\n"; 
      cout << "Node: " << node << " v2 = " << v2 << "!\n"; 
      exit(0);
    }	  
    SecVF[node] = v2;
  }	  
}


// -------------------------------------------------------------------------
// Class cSecFGM3D:
// -------------------------------------------------------------------------

// ============================== cSecFGM3D ================================

cSecFGM3D :: cSecFGM3D(int id) : cSection(id)
{
  Type = SEC_FGM_3D;
}

// ============================= ~cSecFGM3D ================================

cSecFGM3D :: ~cSecFGM3D(void)
{
}

// ================================== Read =================================

void cSecFGM3D :: Read(void)
{
  int mat,nnode;
  if (!(in >> mat) || !(in >> nnode))
  {
    cout << "Error in the input of section " << Label << "!\n"; 
    exit(0);
  }

  Mat = cMaterial::GetMaterial(mat);
  if (!Mat)
  {
    cout << "Error in the input of section " << Label << "!\n";
    cout << "Invalid material!\n";

    exit(0);
  }

  int node;
  double v2;
  for (int i = 0; i < nnode; i++)
  {	  
    if (!(in >> node) || !(in >> v2) || !cNode::FindNode(node))
    {
      cout << "Error in the input of section " << Label << "!\n"; 
      cout << "Node: " << node << " v2 = " << v2 << "!\n"; 
      exit(0);
    }	  
    SecVF[node] = v2;
  }	  
}


// -------------------------------------------------------------------------
// Class cSecHSDTPlate:
// -------------------------------------------------------------------------

// ============================ cSecHSDTPlate ==============================

cSecHSDTPlate :: cSecHSDTPlate(int id) : cSection(id)
{
  Type   = SEC_HSDT_PLATE;
  Theory = FSDT;
  Thk    = 1.0;
}

// =========================== ~cSecHSDTPlate ==============================

cSecHSDTPlate :: ~cSecHSDTPlate(void)
{
}

// ================================== Read =================================

void cSecHSDTPlate :: Read(void)
{
  int mat;
  string str;
  if (!(in >> mat) || !(in >> Thk) || !(in >> str))
  {
    cout << "Error in the input of section " << Label << "!\n";
    exit(0);
  }

  if (str == "FSDT")
    Theory = FSDT;
  else if (str == "TSDT")
    Theory = TSDT;
  else if (str == "ESDT")
    Theory = ESDT;
  else if (str == "SSDT")
    Theory = SSDT;
  else if (str == "SHSDT")
    Theory = SHSDT;
  else
  {
    cout << "Error in the input of section " << Label << "!\n";
    cout << "Invalid plate theory " << str << "!\n";
    exit(0);
  }

  Mat = cMaterial::GetMaterial(mat);
  if (!Mat)
  {
    cout << "Error in the input of section " << Label << "!\n";
    cout << "Invalid material!\n";
    exit(0);
  }
}


// -------------------------------------------------------------------------
// Class cSecHSDTFGMPlate:
// -------------------------------------------------------------------------

// ========================== cSecHSDTFGMPlate =============================

cSecHSDTFGMPlate :: cSecHSDTFGMPlate(int id) : cSection(id)
{
  Type   = SEC_HSDT_FGM_PLATE;
  Theory = FSDT;  
  Thk    = 0.1;  // Thickness
  ExpN   = 1.0;  // Exponent of the volume fraction
  NumIpt = 5;    // Integration points in the thickness
}

// ========================= ~cSecHSDTFGMPlate =============================

cSecHSDTFGMPlate :: ~cSecHSDTFGMPlate(void)
{
}

// ================================== Read =================================

void cSecHSDTFGMPlate :: Read(void)
{
  int mat;
  string str;
  if (!(in >> mat) || !(in >> Thk) || !(in >> ExpN) || !(in >> NumIpt) || 
      !(in >> str))
  {
    cout << "Error in the input of section " << Label << "!\n";
    exit(0);
  }

  if (str == "FSDT")
    Theory = FSDT;
  else if (str == "TSDT")
    Theory = TSDT;
  else if (str == "ESDT")
    Theory = ESDT;
  else if (str == "SSDT")
    Theory = SSDT;
  else if (str == "SHSDT")
    Theory = SHSDT;
  else
  {
    cout << "Error in the input of section " << Label << "!\n";
    cout << "Invalid plate theory " << str << "!\n";
    exit(0);
  }

  Mat = cMaterial::GetMaterial(mat);
  if (!Mat)
  {
    cout << "Error in the input of section " << Label << "!\n";
    cout << "Invalid material!\n";

    exit(0);
  }
}


// -------------------------------------------------------------------------
// Class cSecFGMShell:
// -------------------------------------------------------------------------

// ============================ cSecFGMShell ===============================

cSecFGMShell :: cSecFGMShell(int id) : cSection(id)
{
  Type = SEC_FGM_SHELL;
  Thk    = 1.0;  // Thickness
  ExpN   = 1.0;  // Exponent of the volume fraction
  NumIpt = 3;    // Integration points in the thickness
}

// =========================== ~cSecFGMShell ===============================

cSecFGMShell :: ~cSecFGMShell(void)
{
}

// ================================== Read =================================

void cSecFGMShell :: Read(void)
{
  cout << "cSecFGMShell :: Read\n";

  int mat;
  if (!(in >> mat) || !(in >> Thk) || !(in >> ExpN) || !(in >> NumIpt))
  {
    cout << "Error in the input of section " << Label << "!\n"; 
    exit(0);
  }
  cout << "Mat    = " << mat    << "  ";
  cout << "Thk    = " << Thk    << "  ";
  cout << "ExpN   = " << ExpN   << "  ";
  cout << "NumIpt = " << NumIpt << "\n";

  Mat = cMaterial::GetMaterial(mat);
  if (!Mat)
  {
    cout << "Error in the input of section " << Label << "!\n";
    cout << "Invalid material!\n";

    exit(0);
  }
  cout << "End of cSecFGMShell :: Read\n";
}


// -------------------------------------------------------------------------
// Class cSecGeneralShell:
// -------------------------------------------------------------------------

// ========================== cSecGeneralShell =============================

cSecGeneralShell :: cSecGeneralShell(int id) : cSection(id)
{
  Type = SEC_GENERAL_SHELL;
  C.Resize(8,8);
  Mb.Resize(5,5);
}

// =========================== ~cSecFGMShell ===============================

cSecGeneralShell :: ~cSecGeneralShell(void)
{
}

// ================================== Read =================================

void cSecGeneralShell :: Read(void)
{
  cout << "cSecGeneralShell :: Read\n";
  
  for (int i = 0; i < 8; i++)
  {
    for (int j = 0; j < 8; j++)
    {
      if (!(in >> C[i][j]))
      {
        cout << "Error in the input of section " << Label << "!\n";
        exit(0);
      }
    }
  }
  cout << "[C]\n";
  C.Print( );

  for (int i = 0; i < 5; i++)
  {
    for (int j = 0; j < 5; j++)
    {
      if (!(in >> Mb[i][j]))
      {
        cout << "Error in the input of section " << Label << "!\n";
        exit(0);
      }
    }
  }
  cout << "[Mb]\n";
  Mb.Print( );

  cout << "End of cSecGeneralShell :: Read\n";
}


// -------------------------------------------------------------------------
// Class cSecLaminated2D:
// -------------------------------------------------------------------------

// =========================== cSecLaminated2D =============================

cSecLaminated2D :: cSecLaminated2D(int id) : cSection(id)
{
  Type = SEC_LAMINATED_2D;
  NumLam = 0;
  LamThk = 0.0;
  Layup  = 0;
}

// =========================== ~cSecLaminated2D ============================

cSecLaminated2D :: ~cSecLaminated2D(void)
{
  delete []Layup;
}

// ================================ Read ===================================

void cSecLaminated2D :: Read(void)
{
  // Read the section layup.

  if (!(in >> NumLam))
  {
    cout << "Error in the input of laminated shell section " << Label << "!\n";
    exit(0);
  }

  if (NumLam < 2 || (NumLam%2 != 0))
  {
    cout << "Invalid number of layers (" << NumLam << ") in section ";
    cout << Label << ".\n";
    exit(0);
  }

  Layup = new sLamina[NumLam];

  int mat;
  double thk;
  LamThk = 0.0;
  for (int i = 0; i < NumLam; i++)
  {
    if (!(in >> mat) || !(in >> thk) || !(in >> Layup[i].Ang))
    {
      cout << "Error on reading laminated shell section data\n";
      exit(0);
    }

    LamThk += thk;
    Layup[i].Thk = thk;
    Layup[i].Mat = cMaterial::GetMaterial(mat);

    if (!Layup[i].Mat)
    {
      cout << "Error in the input of section " << Label << "!\n";
      cout << "Invalid material in layer " << i+1 << "!\n";
      exit(0);
    }

    // Check the layup symmetry.

    if (i >= NumLam/2)
    {
      int j = NumLam - i;
      if (Layup[i].Mat != Layup[j].Mat || Layup[i].Thk != Layup[j].Thk ||
          Layup[i].Ang != Layup[j].Ang)
      {
        cout << "Error: asymmetric layup in section " << Label << "!\n";
        exit(0);
      }
    }
  }
}

// ============================== GetDensity ===============================

double cSecLaminated2D :: GetDensity(void)
{
  double dens = 0.0;
  for (int i = 0; i < NumLam; i++)
    dens += Layup[i].Thk*Layup[i].Mat->GetDensity( );

  return(dens/LamThk);  // Return the average density
}


// -------------------------------------------------------------------------
// Class cSecLaminatedShell:
// -------------------------------------------------------------------------

// ======================== cSecLaminatedShell =============================

cSecLaminatedShell :: cSecLaminatedShell(int id) : cSecLaminated2D(id)
{
  Type = SEC_LAMINATED_SHELL;

  OriVec.Resize(3);
  OriVec[0] = 1.0;
  OriVec[1] = 0.0;
  OriVec[2] = 0.0;
}

// ======================== ~cSecLaminatedShell ============================

cSecLaminatedShell :: ~cSecLaminatedShell(void)
{
}

// ================================ Read ===================================

void cSecLaminatedShell :: Read(void)
{
  // Read the orientation vector.

  double a1,a2,a3;
  if (!(in >> a1) || !(in >> a2) || !(in >> a3))
  {
    cout << "Error in the input of section orientation vector\n";
    exit(0);
  }
  OriVec[0] = a1;
  OriVec[1] = a2;
  OriVec[2] = a3;

  // Read the section layup.

  if (!(in >> NumLam))
  {
    cout << "Error in the input of laminated shell section " << Label << "!\n";
    exit(0);
  }

  if (NumLam < 1)
  {
    cout << "Invalid number of layers (" << NumLam << ") in section ";
    cout << Label << ".\n";
    exit(0);
  }

  Layup = new sLamina[NumLam];

  int mat;
  double thk;
  LamThk = 0.0;
  for (int i = 0; i < NumLam; i++)
  {
    if (!(in >> mat) || !(in >> thk) || !(in >> Layup[i].Ang))
    {
      cout << "Error on reading laminated shell section data\n";
      exit(0);
    }

    LamThk += thk;
    Layup[i].Thk = thk;
    Layup[i].Mat = cMaterial::GetMaterial(mat);

    if (!Layup[i].Mat)
    {
      cout << "Error in the input of section " << Label << "!\n";
      cout << "Invalid material in layer " << i+1 << "!\n";
      exit(0);
    }
  }
}


// -------------------------------------------------------------------------
// Class cSecLaminatedSolid:
// -------------------------------------------------------------------------

// ======================== cSecLaminatedSolid =============================

cSecLaminatedSolid :: cSecLaminatedSolid(int id) : cSection(id)
{
  Type = SEC_LAMINATED_SOLID;
  NumLam = 0;
  Layup  = 0;

  OriVec.Resize(3);
  OriVec[0] = 1.0;
  OriVec[1] = 0.0;
  OriVec[2] = 0.0;
}

// ======================== ~cSecLaminatedShell ============================

cSecLaminatedSolid :: ~cSecLaminatedSolid(void)
{
  delete []Layup;
}

// ================================ Read ===================================

void cSecLaminatedSolid :: Read(void)
{
  // Read the orientation vector.

  double a1,a2,a3;
  if (!(in >> a1) || !(in >> a2) || !(in >> a3))
  {
    cout << "Error in the input of laminate orientation vector\n";
    exit(0);
  }
  OriVec.Resize(3);
  OriVec[0] = a1;
  OriVec[1] = a2;
  OriVec[2] = a3;

  // Read the laminate data (layup).

  if (!(in >> NumLam))
  {
    cout << "Error in the input of laminated shell section " << Label << "!\n";
    exit(0);
  }

  if (NumLam < 1)
  {
    cout << "Invalid number of layers: " << NumLam << ".";
    exit(0);
  }

  Layup = new sLamina[NumLam];

  int mat;
  double thk;
  double thksum = 0.0;
  for (int i = 0; i < NumLam; i++)
  {
    if (!(in >> mat) || !(in >> thk) || !(in >> Layup[i].Ang))
    {
      cout << "Error on reading laminated solid section data\n";
      exit(0);
    }

    thksum += thk;
    Layup[i].Thk = thk;
    Layup[i].Mat = cMaterial::GetMaterial(mat);

    if (!Layup[i].Mat)
    {
      cout << "Error in the input of section " << Label << "!\n";
      cout << "Invalid material in layer " << i+1 << "!\n";
      exit(0);
    }
  }

  // Thickness normalization.

  for (int i = 0; i < NumLam; i++)
    Layup[i].Thk /= thksum;
}

// ============================== GetDensity ===============================

double cSecLaminatedSolid :: GetDensity(void)
{
  double dens = 0.0;
  for (int i = 0; i < NumLam; i++)
    dens += Layup[i].Thk*Layup[i].Mat->GetDensity( );

  return(dens);  // Return the average density
}


// -------------------------------------------------------------------------
// Class cSecInterface2D:
// -------------------------------------------------------------------------

// ============================ cSecInterface2D ============================

cSecInterface2D :: cSecInterface2D(int id) : cSection(id)
{
  Type = SEC_INTERFACE_2D;
}

// =========================== ~cSecInterface2D ============================

cSecInterface2D :: ~cSecInterface2D(void)
{
}

// ================================ Read ===================================

void cSecInterface2D :: Read(void)
{
  int mat;
  if (!(in >> mat) || !(in >> Thk))
  {
    cout << "Error in the input of section " << Label << "!\n";
    exit(0);
  }

  Mat = cMaterial::GetMaterial(mat);
  if (!Mat)
  {
    cout << "Error in the input of section " << Label << "!\n";
    cout << "Invalid material!\n";
    exit(0);
  }
}


// -------------------------------------------------------------------------
// Class cSecInterface3D:
// -------------------------------------------------------------------------

// ============================ cSecInterface3D ============================

cSecInterface3D :: cSecInterface3D(int id) : cSection(id)
{
  Type = SEC_INTERFACE_3D;
}

// =========================== ~cSecInterface3D ============================

cSecInterface3D :: ~cSecInterface3D(void)
{
}

// ================================ Read ===================================

void cSecInterface3D :: Read(void)
{
  int mat;
  if (!(in >> mat))
  {
    cout << "Error in the input of section " << Label << "!\n";
    exit(0);
  }

  Mat = cMaterial::GetMaterial(mat);
  if (!Mat)
  {
    cout << "Error in the input of section " << Label << "!\n";
    cout << "Invalid material!\n";
    exit(0);
  }
}

// ======================================================= End of file =====
