// -------------------------------------------------------------------------
// element.cpp - implementation of element class.
// -------------------------------------------------------------------------
// Created:      16-Nov-2000     Evandro Parente Junior
//
// Modified:     03-Oct-2001     Evandro Parente Junior
//               Use of vector/matrix classes.
//
// Modified:     10-Oct-2013    Iuri Barcelos Rocha/Evandro Parente Junior
//               Creation of Section class, implementation of TL elements.
//
// Modified:     08-Apr-2015     Elias Saraiva Barroso
//               Input/output operations using c++ streams.
//
// Modified:     06-Sep-2015     Evandro Parente Junior
//               Implementation of 3D truss elements.
//
// Modified:     13-Oct-2015     Evandro Parente Junior
//               Evaluation of support reactions.
//
// Modified:     16-Sep-2022     Pedro Ygor Rodrigues Mesquita
//               Implementation of HSDT plate elements.
//
// Modified:     09-Jun-2023     Evandro Parente Junior
//               Implementation of tetrahedral elements.
// -------------------------------------------------------------------------

#include <iostream>
#include <stdlib.h>
#include <string>

using namespace std;

#include "element.h"
#include "elmparam.h"
#include "elmparamC1.h"
#include "elmdsh.h"
#include "elminterf.h"
#include "elmpltrs.h"
#include "elmtrs3d.h"
#include "elmplfrm.h"
#include "elmgrid.h"
#include "elmfrm3d.h"
#include "elmfrm3dcr.h"
#include "intpoint.h"
#include "node.h"
#include "shape.h"
#include "shpbsp.h"
#include "shpbez.h"
#include "material.h"
#include "secbar.h"
#include "vec.h"
#include "mat.h"
#include "sysmat.h"
#include "gblvar.h"
#include "gbldef.h"
#include "patch.h"

// -------------------------------------------------------------------------
// Static variables:
//
int          cElement :: NumElm         = 0;
cElement**   cElement :: VecElm         = 0;
int          cElement :: NumIntOrd      = 0;
sIntOrd*     cElement :: VecIntOrd      = 0;
int          cElement :: MaxNode        = 0;
int          cElement :: MaxDof         = 0;
int          cElement :: MaxIntPnt      = 0;
int          cElement :: NumSclLab      = 0;
int*         cElement :: VecSclLab      = 0;
int          cElement :: NumNodSclLab   = 0;
int*         cElement :: VecNodSclLab   = 0;
int          cElement :: IniStrFlag     = 0;
eAxis        cElement :: VertAxis       = Y_AXIS;
int          cElement :: NumSoil        = 0;
sGeoStat*    cElement :: VecSoil        = 0;
eOutConv     cElement :: OutputConv     = DIRECT_STIFFNESS;
int          cElement :: NumBeamOrt     = 0;
sBeamOrt*    cElement :: VecBeamOrt     = 0;
int          cElement :: LargeRot3D     = 0;
int          cElement :: NumElemTemp    = 0;
sTempData*   cElement :: ElemTemp      = NULL;

// -------------------------------------------------------------------------
// Public methods:
//

// ============================ ReadGeoVertAxis ============================

void cElement :: ReadGeoVertAxis(void)
{
  // Read the vertical axis.

  int axis;
  if (!(in >> axis) || (axis < 1) || (axis > 3))
  {
    cout << "Error in the input of the vertical axis!\n";
    exit(0);
  }
  VertAxis = (eAxis)axis;
}

// ============================== ReadGeoStat ==============================

void cElement :: ReadGeoStat(void)
{
  // Read the number of soil layers.

  if (!(in >> NumSoil) || (NumSoil <= 0))
  {
    cout << "Error in the input of the number of soil layers!\n";
    exit(0);
  }
  VecSoil = new sGeoStat[NumSoil];

  // Read soil layers.

  double top,gam,k0;
  for (int i = 0; i < NumSoil; i++)
  {
    if (!(in >> top) || !(in >> gam) || !(in >> k0))
    {
      cout << "Error in the input of soil layers!\n";
      exit(0);
    }
    VecSoil[i].top   = top;
    VecSoil[i].gamma = gam;
    VecSoil[i].k0    = k0;
  }

  // Activate initial stress flag.

  IniStrFlag = 1;
}

// =========================== ReadInitialStress ===========================

void cElement :: ReadInitialStress(void)
{
  int nelm;
  if (!(in >> nelm))
  {
    cout << "Error in the input of initial stresses!\n";
    cout << "Error in the input of initial stresses!\n";
    exit(0);
  }

  for (int i = 0; i < nelm; i++)
  {
    int elmid;
    cElement *elem;
    if (!(in >> elmid))
    {
      cout << "Error in the input of initial stresses!\n";
      exit(0);
    }
    elem = cElement :: FindElm(elmid);
    if (!elem)
    {
      cout << "Error in the input of initial stresses!\n";
      exit(0);
    }
    elem->ReadIniStr( );
  }

  // Activate initial stress flag.

  IniStrFlag = 1;
}

// ============================ ReadOutputConv =============================

void cElement :: ReadOutputConv(void)
{
  // Read the convention for output of frame element forces.

  int conv;
  if (!(in >> conv))
  {
    cout << "Error in the input of output convention for frame forces!\n";
    exit(0);
  }

  // Store the convention.

  if (conv == 2)
    OutputConv = STATICS;
  else
    OutputConv = DIRECT_STIFFNESS;
}

// ============================== ReadIntOrd ===============================

void cElement :: ReadIntOrd(void)
{
  // Read the number of integration orders.

  if (!(in >> NumIntOrd) || NumIntOrd <= 0)
  {
    cout << "Error in the input of the number of integration orders!\n";
    exit(0);
  }
  VecIntOrd = new sIntOrd[NumIntOrd];

  // Read integration orders.

  int id,ordk[3],ordm[3];
  for (int i = 0; i < NumIntOrd; i++)
  {
    // Use only the first number as the integration order.

    if (!(in >> id)       || !(in >> ordk[0]) || !(in >> ordk[1]) ||
        !(in >> ordk[2])  || !(in >> ordm[0]) || !(in >> ordm[1]) ||
	!(in >> ordm[2]) )
    {
      cout << "Error in the input of integration orders elements!\n";
      exit(0);
    }
    VecIntOrd[i].K[0] = ordk[0];
    VecIntOrd[i].K[1] = ordk[1];
    VecIntOrd[i].K[2] = ordk[2];
    VecIntOrd[i].M[0] = ordm[0];
    VecIntOrd[i].M[1] = ordm[1];
    VecIntOrd[i].M[2] = ordm[2];

    VecIntOrd[i].type = 0; // Default is Gauss
  }
}

// ============================== ReadIntType ==============================

void cElement :: ReadIntType(void)
{
  // Read the number of integration types.

  int num, id;

  if (!(in >> num) || num  <= 0)
  {
    cout << "Error in the input of the number of integration types!\n";
    exit(0);
  }

  // Read integration types and apply it.

  for (int i = 0; i < num; i++)
  {
    // Read the integration order index

    if (!(in >> id) || id <= 0 || (id-1) > NumIntOrd)
    {
      cout << "Error in the input of integration order index!\n";
      exit(0);
    }

    // Read the quadrature type and store it in the integration order data

    eQuadType type;
    in >> type;

    VecIntOrd[id-1].type = static_cast<int>(type);
  }
}

// ============================== ReadBeamOrt ==============================

void cElement :: ReadBeamOrt(void)
{
  // Read the number of beam orientation.

  if (!(in >> NumBeamOrt) || NumBeamOrt <= 0)
  {
    cout << "Error on reading number of beam orientation sequence!\n";
    exit(0);
  }
  if (NumBeamOrt == 0) return;

  // Allocate and read the beam orientations.

  int i,orientid;
  double l,m,n;
  VecBeamOrt = new sBeamOrt[NumBeamOrt];

  for (i = 0; i < NumBeamOrt; i++)
  {
    if (!(in >> orientid) || !(in >> l) || !(in >> m) || !(in >> n))
    {
      cout << "Error on reading angle of orientation " << i+1 << "!\n";
      exit(0);
    }

    VecBeamOrt[i].idort = orientid;
    VecBeamOrt[i].e0y[0] = l;  // ly
    VecBeamOrt[i].e0y[1] = m;  // my
    VecBeamOrt[i].e0y[2] = n;  // ny
  }
}

// ============================== ReadNumElm ===============================

void cElement :: ReadNumElm(void)
{
  // Read the number of mesh elements.

  if (!(in >> NumElm) || NumElm <= 0)
  {
    cout << "Error in the input of the number of elements!\n";
    exit(0);
  }

  // Alloc the array of elements.

  VecElm = new cElement*[NumElm];
}

// ============================== ReadPlTruss ==============================

void cElement :: ReadPlTruss(void)
{
  // Read the number of plane truss elements.

  int n;
  if (!(in >> n) || (n < 1) || (n > NumElm))
  {
    cout << "Error in the input of the number of PLTRUSS elements!\n";
    exit(0);
  }

  // Create and read each element.

  int id;
  cElement *elm;
  for (int i = 0; i < n; i++)
  {
    if (!(in >> id) || (id < 1) || (id > NumElm))
    {
      cout << "Error in the input of PLTRUSS elements!\n";
      exit(0);
    }
    elm = new cPlTruss(id);
    elm->Read( );
  }
  elm->UpdGlbVars( );  // Upd. global vars depending on the element type
}

// ============================= ReadPlTrussTL =============================

void cElement :: ReadPlTrussTL(void)
{
  // Read the number of Total Lagrangian plane truss elements.

  int n;
  if (!(in >> n) || (n < 1) || (n > NumElm))
  {
    cout << "Error in the input of the number of PLTRUSSTL elements!\n";
    exit(0);
  }

  // Create and read each element.

  int id;
  cElement *elm;
  for (int i = 0; i < n; i++)
  {
    if (!(in >> id) || (id < 1) || (id > NumElm))
    {
      cout << "Error in the input of PLTRUSSTL element!\n";
      exit(0);
    }
    elm = new cPlTrussTL(id);
    elm->Read( );
  }
  elm->UpdGlbVars( );  // Upd. global vars depending on the element type
}

// ============================= ReadPlTrussCR =============================

void cElement :: ReadPlTrussCR(void)
{
  // Read the number of corotational plane truss elements.

  int n;
  if (!(in >> n) || (n < 1) || (n > NumElm))
  {
    cout << "Error in the input of the number of PLTRUSSCR elements!\n";
    exit(0);
  }

  // Create and read each element.

  int id;
  cElement *elm;
  for (int i = 0; i < n; i++)
  {
    if (!(in >> id) || (id < 1) || (id > NumElm))
    {
      cout << "Error in the input of PLTRUSSCR element!\n";
      exit(0);
    }
    elm = new cPlTrussCR(id);
    elm->Read( );
  }
  elm->UpdGlbVars( );  // Upd. global vars depending on the element type
}

// ============================== ReadTruss3D ==============================

void cElement :: ReadTruss3D(void)
{
  // Read the number of 3D truss elements.

  int n;
  if (!(in >> n) || (n < 1) || (n > NumElm))
  {
    cout << "Error in the input of the number of TRUSS3D elements!\n";
    exit(0);
  }

  // Create and read each element.

  int id;
  cElement *elm;
  for (int i = 0; i < n; i++)
  {
    if (!(in >> id) || (id < 1) || (id > NumElm))
    {
      cout << "Error in the input of TRUSS3D elements!\n";
      exit(0);
    }
    elm = new cTruss3D(id);
    elm->Read( );
  }
  elm->UpdGlbVars( );  // Upd. global vars depending on the element type
}

// ============================= ReadTruss3DTL =============================

void cElement :: ReadTruss3DTL(void)
{
  // Read the number of Total Lagrangian 3D truss elements.

  int n;
  if (!(in >> n) || (n < 1) || (n > NumElm))
  {
    cout << "Error in the input of the number of TRUSS3DTL elements!\n";
    exit(0);
  }

  // Create and read each element.

  int id;
  cElement *elm;
  for (int i = 0; i < n; i++)
  {
    if (!(in >> id) || (id < 1) || (id > NumElm))
    {
      cout << "Error in the input of TRUSS3DTL element!\n";
      exit(0);
    }
    elm = new cTruss3DTL(id);
    elm->Read( );
  }
  elm->UpdGlbVars( );  // Upd. global vars depending on the element type
}

// ============================= ReadTruss3DCR =============================

void cElement :: ReadTruss3DCR(void)
{
  // Read the number of corotational 3D truss elements.

  int n;
  if (!(in >> n) || (n < 1) || (n > NumElm))
  {
    cout << "Error in the input of the number of TRUSS3DCR elements!\n";
    exit(0);
  }

  // Create and read each element.

  int id;
  cElement *elm;
  for (int i = 0; i < n; i++)
  {
    if (!(in >> id) || (id < 1) || (id > NumElm))
    {
      cout << "Error in the input of TRUSS3DCR element!\n";
      exit(0);
    }
    elm = new cTruss3DCR(id);
    elm->Read( );
  }
  elm->UpdGlbVars( );  // Upd. global vars depending on the element type
}

// ============================== ReadPlFrame ==============================

void cElement :: ReadPlFrame(void)
{
  // Read the number of plane frame elements.

  int n;
  if (!(in >> n) || (n < 1) || (n > NumElm))
  {
    cout << "Error in the input of the number of PLFRAME elements!\n";
    exit(0);
  }

  // Create and read each element.

  int id,sid;
  cElement *elm;
  cSection *sec;
  for (int i = 0; i < n; i++)
  {
    if (!(in >> id) || !(in >> sid) || (id < 1) || (id > NumElm))
    {
      cout << "Error in the input of PLFRAME elements!\n";
      exit(0);
    }
    sec = cSection::GetSection(sid);
    if (!sec || !cPlFrame::ValidSection(sec->GetType( )))
    {
      cout << "Error in the input of element " << id << "!\nInvalid Section!\n";
      exit(0);
    }
    if (sec->PreIntegrated( ))
      elm = new cPlFrame(id, sec);
    else
      elm = new cPlFrameNI(id, sec);
    elm->Read( );
  }
  elm->UpdGlbVars( );  // Upd. global vars depending on the element type
}

// ============================= ReadPlFrameCR1 ============================

void cElement :: ReadPlFrameCR1(void)
{
  // Read the number of corotational plane frame elements.

  int n;
  if (!(in >> n) || (n < 1) || (n > NumElm))
  {
    cout << "Error in the input of the number of PLFRAMECR1 elements!\n";
    exit(0);
  }

  // Create and read each element.

  int id,sid;
  cElement *elm;
  cSection *sec;
  for (int i = 0; i < n; i++)
  {
    if (!(in >> id) || !(in >> sid) || (id < 1) || (id > NumElm))
    {
      cout << "Error in the input of PLFRAMECR1 elements!\n";
      exit(0);
    }
    sec = cSection::GetSection(sid);
    if (!sec || !cPlFrame::ValidSection(sec->GetType( )))
    {
      cout << "Error in the input of element " << id << "!\nInvalid Section!\n";
      exit(0);
    }
    if (sec->PreIntegrated( ))
      elm = new cPlFrameCR1(id, sec);
    else
      elm = new cPlFrameCR1NI(id, sec);
    elm->Read( );
  }
  elm->UpdGlbVars( );  // Upd. global vars depending on the element type
}

// ============================ ReadPlFrameCR2 =============================

void cElement :: ReadPlFrameCR2(void)
{
  // Read the number of enhanced corotational plane frame elements.

  int n;
  if (!(in >> n) || (n < 1) || (n > NumElm))
  {
    cout << "Error in the input of the number of PLFRAMECR2 elements!\n";
    exit(0);
  }

  // Create and read each element.

  int id,sid;
  cElement *elm;
  cSection *sec;
  for (int i = 0; i < n; i++)
  {
    if (!(in >> id) || !(in >> sid) || (id < 1) || (id > NumElm))
    {
      cout << "Error in the input of PLFRAMECR2 elements!\n";
      exit(0);
    }
    sec = cSection::GetSection(sid);
    if (!sec || !cPlFrame::ValidSection(sec->GetType( )))
    {
      cout << "Error in the input of element " << id << "!\nInvalid Section!\n";
      exit(0);
    }
    if (sec->PreIntegrated( ))
      elm = new cPlFrameCR2(id, sec);
    else
      elm = new cPlFrameCR2NI(id, sec);
    elm->Read( );
  }
  elm->UpdGlbVars( );  // Upd. global vars depending on the element type
}

// ============================= ReadPlFrameTL1 ============================

void cElement :: ReadPlFrameTL1(void)
{
  // Read the number of enhanced corotational plane frame elements.

  int n;
  if (!(in >> n) || (n < 1) || (n > NumElm))
  {
    cout << "Error in the input of the number of PLFRAMETL1 elements!\n";
    exit(0);
  }

  // Create and read each element.

  int id,sid;
  cElement *elm;
  cSection *sec;
  for (int i = 0; i < n; i++)
  {
    if (!(in >> id) || !(in >> sid) || (id < 1) || (id > NumElm))
    {
      cout << "Error in the input of PLFRAMETL1 elements!\n";
      exit(0);
    }
    sec = cSection::GetSection(sid);
    if (!sec || !cPlFrame::ValidSection(sec->GetType( )))
    {
      cout << "Error in the input of element " << id << "!\nInvalid Section!\n";
      exit(0);
    }
    elm = new cPlFrameTL1(id, sec);
    elm->Read( );
  }
  elm->UpdGlbVars( );  // Upd. global vars depending on the element type
}

// ============================= ReadPlFrameTL2 ============================

void cElement :: ReadPlFrameTL2(void)
{
  // Read the number of enhanced corotational plane frame elements.

  int n;
  if (!(in >> n) || (n < 1) || (n > NumElm))
  {
    cout << "Error in the input of the number of PLFRAMETL2 elements!\n";
    exit(0);
  }

  // Create and read each element.

  int id,sid;
  cElement *elm;
  cSection *sec;
  for (int i = 0; i < n; i++)
  {
    if (!(in >> id) || !(in >> sid) || (id < 1) || (id > NumElm))
    {
      cout << "Error in the input of PLFRAMETL2 elements!\n";
      exit(0);
    }
    sec = cSection::GetSection(sid);
    if (!sec || !cPlFrame::ValidSection(sec->GetType( )))
    {
      cout << "Error in the input of element " << id << "!\nInvalid Section!\n";
      exit(0);
    }
    elm = new cPlFrameTL2(id, sec);
    elm->Read( );
  }
  elm->UpdGlbVars( );  // Upd. global vars depending on the element type
}

// ============================== ReadGrid =================================

void cElement :: ReadGrid(void)
{
  // Read the number of three dimensional elements.

  int n;
  if (!(in >> n) || (n < 1) || (n > NumElm))
  {
    cout << "Error in the input of the number of FRAME3D elements!\n";
    exit(0);
  }

  // Create and read each element.

  int id,sid;
  cElement *elm;
  cSection *sec;
  for (int i = 0; i < n; i++)
  {
    if (!(in >> id) || !(in >> sid) || (id < 1) || (id > NumElm))
    {
      cout << "Error in the input of GRID elements!\n";
      exit(0);
    }
    sec = cSection::GetSection(sid);
    if (!sec || !cGrid::ValidSection(sec->GetType( )))
    {
      cout << "Error in the input of element " << id << "!\nInvalid Section!\n";
      exit(0);
    }
    if (sec->PreIntegrated( ))
      elm = new cGrid(id, sec);
    else
    {
      cout << "grid element with nonlinear material not implemented yet!\n";
      exit(0);
    }
    elm->Read( );
  }
  elm->UpdGlbVars( );  // Upd. global vars depending on the element type
}

// ============================== ReadFrame3D ==============================

void cElement :: ReadFrame3D(void)
{
  // Read the number of three dimensional elements.

  int n;
  if (!(in >> n) || (n < 1) || (n > NumElm))
  {
    cout << "Error in the input of the number of FRAME3D elements!\n";
    exit(0);
  }

  // Create and read each element.

  int id,sid;
  cElement *elm;
  cSection *sec;
  for (int i = 0; i < n; i++)
  {
    if (!(in >> id) || !(in >> sid) || (id < 1) || (id > NumElm))
    {
      cout << "Error in the input of FRAME3D elements!\n";
      exit(0);
    }
    sec = cSection::GetSection(sid);
    if (!sec || !cFrame3D::ValidSection(sec->GetType( )))
    {
      cout << "Error in the input of element " << id << "!\nInvalid Section!\n";
      exit(0);
    }
    if (sec->PreIntegrated( ))
      elm = new cFrame3D(id, sec);
    else
    {
      cout << "3D frame element with nonlinear material not implemented yet!\n";
      exit(0);
    }
    elm->Read( );
  }
  elm->UpdGlbVars( );  // Upd. global vars depending on the element type
}

// ============================= ReadFrame3DTL =============================

void cElement :: ReadFrame3DTL(void)
{
  // Read the number of three dimensional elements.

  int n;
  if (!(in >> n) || (n < 1) || (n > NumElm))
  {
    cout << "Error in the input of the number of FRAME3DTL elements!\n";
    exit(0);
  }

  // Create and read each element.

  int id,sid;
  cElement *elm;
  cSection *sec;
  for (int i = 0; i < n; i++)
  {
    if (!(in >> id) || !(in >> sid) || (id < 1) || (id > NumElm))
    {
      cout << "Error in the input of FRAME3D elements!\n";
      exit(0);
    }
    sec = cSection::GetSection(sid);
    if (!sec || !cFrame3D::ValidSection(sec->GetType( )))
    {
      cout << "Error in the input of element " << id << "!\nInvalid Section!\n";
      exit(0);
    }
    if (sec->PreIntegrated( ))
      elm = new cFrame3DTL(id, sec);
    else
    {
      cout << "3D frame element with nonlinear material not implemented yet!\n";
      exit(0);
    }
    elm->Read( );
  }
  elm->UpdGlbVars( );  // Upd. global vars depending on the element type
}

// ============================= ReadFrame3DCR =============================

void cElement :: ReadFrame3DCR(void)
{
  // Read the number of corotational three dimensional elements.

  int n;
  if (!(in >> n) || (n < 1) || (n > NumElm))
  {
    cout << "Error in the input of the number of FRAME3DCR elements!\n";
    exit(0);
  }

  // Create and read each element.

  int id, sid;
  cElement *elm;
  cSection *sec;
  for (int i = 0; i < n; i++)
  {
    if (!(in >> id) || !(in >> sid) || (id < 1) || (id > NumElm))
    {
      cout << "Error in the input of FRAME3DCRB elements!\n";
      exit(0);
    }

    sec = cSection::GetSection(sid);
    if (!sec || !cFrame3D::ValidSection(sec->GetType( )))
    {
      cout << "Error in the input of element " << id << "!\nInvalid Section!\n";
      exit(0);
    }
    if (sec->PreIntegrated( ))
      elm = new cFrame3DCR(id, sec);
    else
    {
      cout << "3D frame element with nonlinear material not implemented yet!\n";
      exit(0);
    }
    elm->Read( );
  }
  elm->UpdGlbVars( );  // Upd. global vars depending on the element type
}

// ============================ ReadFrame3DCRTL ============================

void cElement :: ReadFrame3DCRTL(void)
{
  // Read the number of corotational three dimensional elements.

  int n;
  if (!(in >> n) || (n < 1) || (n > NumElm))
  {
    cout << "Error in the input of the number of FRAME3DCRTL elements!\n";
    exit(0);
  }

  // Create and read each element.

  int id, sid;
  cElement *elm;
  cSection *sec;
  for (int i = 0; i < n; i++)
  {
    if (!(in >> id) || !(in >> sid) || (id < 1) || (id > NumElm))
    {
      cout << "Error in the input of FRAME3DCRB elements!\n";
      exit(0);
    }

    sec = cSection::GetSection(sid);
    if (!sec || !cFrame3D::ValidSection(sec->GetType( )))
    {
      cout << "Error in the input of element " << id << "!\nInvalid Section!\n";
      exit(0);
    }
    if (sec->PreIntegrated( ))
      elm = new cFrame3DCRTL(id, sec);
    else
    {
      cout << "3D frame element with nonlinear material not implemented yet!\n";
      exit(0);
    }
    elm->Read( );
  }
  elm->UpdGlbVars( );  // Upd. global vars depending on the element type
}

// ============================= ReadPlStressQ4 ============================

void cElement :: ReadPlStressQ4(void)
{
  // Read the number of Q4 elements.

  int n;
  if (!(in >> n) || (n < 1) || (n > NumElm))
  {
    cout << "Error in the input of the number of PLSTRESSQ4 elements!\n";
    exit(0);
  }

  // Create and read each element.

  int id;
  cElement *elm;
  for (int i = 0; i < n; i++)
  {
    if (!(in >> id) || (id < 1) || (id > NumElm))
    {
      cout << "Error in the input of PLSTRESSQ4 elements!\n";
      exit(0);
    }
    elm = new cElmParam(id, SHAPE_Q4, PLANE_STRESS);
    elm->Read( );
  }
  elm->UpdGlbVars( );  // Upd. global vars depending on the element type
}

// ============================ ReadPlStressQ4TL ===========================

void cElement :: ReadPlStressQ4TL(void)
{
  // Read the number of Q4 elements.

  int n;
  if (!(in >> n) || (n < 1) || (n > NumElm))
  {
    cout << "Error in the input of the number of PLSTRESSQ4 elements!\n";
    exit(0);
  }

  // Create and read each element.

  int id;
  cElement *elm;
  for (int i = 0; i < n; i++)
  {
    if (!(in >> id) || (id < 1) || (id > NumElm))
    {
      cout << "Error in the input of PLSTRESSQ4 elements!\n";
      exit(0);
    }
    elm = new cElmParamTL(id, SHAPE_Q4, PLANE_STRESS);
    elm->Read( );
  }
  elm->UpdGlbVars( );  // Upd. global vars depending on the element type
}

// ============================= ReadPlStressQ8 ============================

void cElement :: ReadPlStressQ8(void)
{
  // Read the number of Q8 elements.

  int n;
  if (!(in >> n) || (n < 1) || (n > NumElm))
  {
    cout << "Error in the input of the number of PLSTRESSQ8 elements!\n";
    exit(0);
  }

  // Create and read each element.

  int id;
  cElement *elm;
  for (int i = 0; i < n; i++)
  {
    if (!(in >> id) || (id < 1) || (id > NumElm))
    {
      cout << "Error in the input of PLSTRESSQ8 elements!\n";
      exit(0);
    }
    elm = new cElmParam(id, SHAPE_Q8, PLANE_STRESS);
    elm->Read( );
  }
  elm->UpdGlbVars( );  // Upd. global vars depending on the element type
}

// ========================== ReadPlStressQ8TL ============================

void cElement :: ReadPlStressQ8TL(void)
{
  // Read the number of Q8 elements.

  int n;
  if (!(in >> n) || (n < 1) || (n > NumElm))
  {
    cout << "Error in the input of the number of PLSTRESSQ8 elements!\n";
    exit(0);
  }

  // Create and read each element.

  int id;
  cElement *elm;
  for (int i = 0; i < n; i++)
  {
    if (!(in >> id) || (id < 1) || (id > NumElm))
    {
      cout << "Error in the input of PLSTRESSQ8 elements!\n";
      exit(0);
    }
    elm = new cElmParamTL(id, SHAPE_Q8, PLANE_STRESS);
    elm->Read( );
  }
  elm->UpdGlbVars( );  // Upd. global vars depending on the element type
}

// ============================= ReadPlStressQ9 ============================

void cElement :: ReadPlStressQ9(void)
{
  // Read the number of Q9 elements.

  int n;
  if (!(in >> n) || (n < 1) || (n > NumElm))
  {
    cout << "Error in the input of the number of PLSTRESSQ9 elements!\n";
    exit(0);
  }

  // Create and read each element.

  int id;
  cElement *elm;
  for (int i = 0; i < n; i++)
  {
    if (!(in >> id) || (id < 1) || (id > NumElm))
    {
      cout << "Error in the input of PLSTRESSQ9 elements!\n";
      exit(0);
    }
    elm = new cElmParam(id, SHAPE_Q9, PLANE_STRESS);
    elm->Read( );
  }
  elm->UpdGlbVars( );  // Upd. global vars depending on the element type
}

// =========================== ReadPlStressQ9TL ============================

void cElement :: ReadPlStressQ9TL(void)
{
  // Read the number of Q9 elements.

  int n;
  if (!(in >> n) || (n < 1) || (n > NumElm))
  {
    cout << "Error in the input of the number of PLSTRESSQ9 elements!\n";
    exit(0);
  }

  // Create and read each element.

  int id;
  cElement *elm;
  for (int i = 0; i < n; i++)
  {
    if (!(in >> id) || (id < 1) || (id > NumElm))
    {
      cout << "Error in the input of PLSTRESSQ9 elements!\n";
      exit(0);
    }
    elm = new cElmParamTL(id, SHAPE_Q9, PLANE_STRESS);
    elm->Read( );
  }
  elm->UpdGlbVars( );  // Upd. global vars depending on the element type
}

// ============================= ReadPlStressT3 ============================

void cElement :: ReadPlStressT3(void)
{
  // Read the number of T3 elements.

  int n;
  if (!(in >> n) || (n < 1) || (n > NumElm))
  {
    cout << "Error in the input of the number of PLSTRESST3 elements!\n";
    exit(0);
  }

  // Create and read each element.

  int id;
  cElement *elm;
  for (int i = 0; i < n; i++)
  {
    if (!(in >> id) || (id < 1) || (id > NumElm))
    {
      cout << "Error in the input of PLSTRESST3 elements!\n";
      exit(0);
    }
    elm = new cElmParam(id, SHAPE_T3, PLANE_STRESS);
    elm->Read( );
  }
  elm->UpdGlbVars( );  // Upd. global vars depending on the element type
}

// =========================== ReadPlStressT3TL ============================

void cElement :: ReadPlStressT3TL(void)
{
  // Read the number of T3 elements.

  int n;
  if (!(in >> n) || (n < 1) || (n > NumElm))
  {
    cout << "Error in the input of the number of PLSTRESST3 elements!\n";
    exit(0);
  }

  // Create and read each element.

  int id;
  cElement *elm;
  for (int i = 0; i < n; i++)
  {
    if (!(in >> id) || (id < 1) || (id > NumElm))
    {
      cout << "Error in the input of PLSTRESST3 elements!\n";
      exit(0);
    }
    elm = new cElmParamTL(id, SHAPE_T3, PLANE_STRESS);
    elm->Read( );
  }
  elm->UpdGlbVars( );  // Upd. global vars depending on the element type
}

// ============================= ReadPlStressT6 ============================

void cElement :: ReadPlStressT6(void)
{
  // Read the number of T6 elements.

  int n;
  if (!(in >> n) || (n < 1) || (n > NumElm))
  {
    cout << "Error in the input of the number of PLSTRESST6 elements!\n";
    exit(0);
  }

  // Create and read each element.

  int id;
  cElement *elm;
  for (int i = 0; i < n; i++)
  {
    if (!(in >> id) || (id < 1) || (id > NumElm))
    {
      cout << "Error in the input of PLSTRESST6 elements!\n";
      exit(0);
    }
    elm = new cElmParam(id, SHAPE_T6, PLANE_STRESS);
    elm->Read( );
  }
  elm->UpdGlbVars( );  // Upd. global vars depending on the element type
}

// =========================== ReadPlStressT6TL ============================

void cElement :: ReadPlStressT6TL(void)
{
  // Read the number of T6 elements.

  int n;
  if (!(in >> n) || (n < 1) || (n > NumElm))
  {
    cout << "Error in the input of the number of PLSTRESST6 elements!\n";
    exit(0);
  }

  // Create and read each element.

  int id;
  cElement *elm;
  for (int i = 0; i < n; i++)
  {
    if (!(in >> id) || (id < 1) || (id > NumElm))
    {
      cout << "Error in the input of PLSTRESST6 elements!\n";
      exit(0);
    }
    elm = new cElmParamTL(id, SHAPE_T6, PLANE_STRESS);
    elm->Read( );
  }
  elm->UpdGlbVars( );  // Upd. global vars depending on the element type
}

// =========================== ReadPlStressBezTri ==========================

void cElement :: ReadPlStressBezTri(void)
{
  cElmParam :: ReadElmFEM("Plane Stress Bezier Triangle",
		           SHAPE_BEZTRI_PLSURF, PLANE_STRESS);
}

// =========================== ReadPlStressBezTriTL ========================

void cElement :: ReadPlStressBezTriTL(void)
{
  cElmParam :: ReadElmFEM("Plane Stress Bezier Triangle",
		           SHAPE_BEZTRI_PLSURF, PLANE_STRESS,true);
}

// =========================== ReadPlStressBezSurf =========================

void cElement :: ReadPlStressBezSurf(void)
{
  cElmParam :: ReadElmFEM("Plane Stress Bezier Surface",
		          SHAPE_BEZ_PLSURF, PLANE_STRESS);
}

// =========================== ReadPlStressBezSurfTL =======================

void cElement :: ReadPlStressBezSurfTL(void)
{
  cElmParam :: ReadElmFEM("Plane Stress Bezier Surface",
		          SHAPE_BEZ_PLSURF, PLANE_STRESS,true);

}

// =========================== ReadPlStressBsp =============================

void cElement :: ReadPlStressBsp(void)
{
  cElmParam :: ReadElmBSP("Plane Stress B-Spline Surface",
		          SHAPE_BSP_PLSURF, PLANE_STRESS);
}

// =========================== ReadPlStressBspTL ===========================

void cElement :: ReadPlStressBspTL(void)
{
  cElmParam :: ReadElmBSP("Plane Stress B-Spline Surface",
		          SHAPE_BSP_PLSURF, PLANE_STRESS,true);
}

// ============================= ReadPlStrainQ4 ============================

void cElement :: ReadPlStrainQ4(void)
{
  // Read the number of Q4 elements.

  int n;
  if (!(in >> n) || (n < 1) || (n > NumElm))
  {
    cout << "Error in the input of the number of PLSTRAINQ4 elements!\n";
    exit(0);
  }

  // Create and read each element.

  int id;
  cElement *elm;
  for (int i = 0; i < n; i++)
  {
    if (!(in >> id) || (id < 1) || (id > NumElm))
    {
      cout << "Error in the input of PLSTRAINQ4 elements!\n";
      exit(0);
    }
    elm = new cElmParam(id, SHAPE_Q4, PLANE_STRAIN);
    elm->Read( );
  }
  elm->UpdGlbVars( );  // Upd. global vars depending on the element type
}

// =========================== ReadPlStrainQ4TL ============================

void cElement :: ReadPlStrainQ4TL(void)
{
  // Read the number of Q4 elements.

  int n;
  if (!(in >> n) || (n < 1) || (n > NumElm))
  {
    cout << "Error in the input of the number of PLSTRAINQ4 elements!\n";
    exit(0);
  }

  // Create and read each element.

  int id;
  cElement *elm;
  for (int i = 0; i < n; i++)
  {
    if (!(in >> id) || (id < 1) || (id > NumElm))
    {
      cout << "Error in the input of PLSTRAINQ4 elements!\n";
      exit(0);
    }
    elm = new cElmParamTL(id, SHAPE_Q4, PLANE_STRAIN);
    elm->Read( );
  }
  elm->UpdGlbVars( );  // Upd. global vars depending on the element type
}

// ============================= ReadPlStrainQ8 ============================

void cElement :: ReadPlStrainQ8(void)
{
  // Read the number of Q8 elements.

  int n;
  if (!(in >> n) || (n < 1) || (n > NumElm))
  {
    cout << "Error in the input of the number of PLSTRAINQ8 elements!\n";
    exit(0);
  }

  // Create and read each element.

  int id;
  cElement *elm;
  for (int i = 0; i < n; i++)
  {
    if (!(in >> id) || (id < 1) || (id > NumElm))
    {
      cout << "Error in the input of PLSTRAINQ8 elements!\n";
      exit(0);
    }
    elm = new cElmParam(id, SHAPE_Q8, PLANE_STRAIN);
    elm->Read( );
  }
  elm->UpdGlbVars( );  // Upd. global vars depending on the element type
}

// =========================== ReadPlStrainQ8TL ============================

void cElement :: ReadPlStrainQ8TL(void)
{
  // Read the number of Q8 elements.

  int n;
  if (!(in >> n) || (n < 1) || (n > NumElm))
  {
    cout << "Error in the input of the number of PLSTRAINQ8 elements!\n";
    exit(0);
  }

  // Create and read each element.

  int id;
  cElement *elm;
  for (int i = 0; i < n; i++)
  {
    if (!(in >> id) || (id < 1) || (id > NumElm))
    {
      cout << "Error in the input of PLSTRAINQ8 elements!\n";
      exit(0);
    }
    elm = new cElmParamTL(id, SHAPE_Q8, PLANE_STRAIN);
    elm->Read( );
  }
  elm->UpdGlbVars( );  // Upd. global vars depending on the element type
}

// ============================= ReadPlStrainQ9 ============================

void cElement :: ReadPlStrainQ9(void)
{
  // Read the number of Q9 elements.

  int n;
  if (!(in >> n) || (n < 1) || (n > NumElm))
  {
    cout << "Error in the input of the number of PLSTRAINQ9 elements!\n";
    exit(0);
  }

  // Create and read each element.

  int id;
  cElement *elm;
  for (int i = 0; i < n; i++)
  {
    if (!(in >> id) || (id < 1) || (id > NumElm))
    {
      cout << "Error in the input of PLSTRAINQ9 elements!\n";
      exit(0);
    }
    elm = new cElmParam(id, SHAPE_Q9, PLANE_STRAIN);
    elm->Read( );
  }
  elm->UpdGlbVars( );  // Upd. global vars depending on the element type
}

// =========================== ReadPlStrainQ9TL ============================

void cElement :: ReadPlStrainQ9TL(void)
{
  // Read the number of Q9 elements.

  int n;
  if (!(in >> n) || (n < 1) || (n > NumElm))
  {
    cout << "Error in the input of the number of PLSTRAINQ9 elements!\n";
    exit(0);
  }

  // Create and read each element.

  int id;
  cElement *elm;
  for (int i = 0; i < n; i++)
  {
    if (!(in >> id) || (id < 1) || (id > NumElm))
    {
      cout << "Error in the input of PLSTRAINQ9 elements!\n";
      exit(0);
    }
    elm = new cElmParamTL(id, SHAPE_Q9, PLANE_STRAIN);
    elm->Read( );
  }
  elm->UpdGlbVars( );  // Upd. global vars depending on the element type
}

// ============================= ReadPlStrainT3 ============================

void cElement :: ReadPlStrainT3(void)
{
  // Read the number of T3 elements.

  int n;
  if (!(in >> n) || (n < 1) || (n > NumElm))
  {
    cout << "Error in the input of the number of PLSTRAINT3 elements!\n";
    exit(0);
  }

  // Create and read each element.

  int id;
  cElement *elm;
  for (int i = 0; i < n; i++)
  {
    if (!(in >> id) || (id < 1) || (id > NumElm))
    {
      cout << "Error in the input of PLSTRAINT3 elements!\n";
      exit(0);
    }
    elm = new cElmParam(id, SHAPE_T3, PLANE_STRAIN);
    elm->Read( );
  }
  elm->UpdGlbVars( );  // Upd. global vars depending on the element type
}

// =========================== ReadPlStrainT3TL ============================

void cElement :: ReadPlStrainT3TL(void)
{
  // Read the number of T3 elements.

  int n;
  if (!(in >> n) || (n < 1) || (n > NumElm))
  {
    cout << "Error in the input of the number of PLSTRAINT3 elements!\n";
    exit(0);
  }

  // Create and read each element.

  int id;
  cElement *elm;
  for (int i = 0; i < n; i++)
  {
    if (!(in >> id) || (id < 1) || (id > NumElm))
    {
      cout << "Error in the input of PLSTRAINT3 elements!\n";
      exit(0);
    }
    elm = new cElmParamTL(id, SHAPE_T3, PLANE_STRAIN);
    elm->Read( );
  }
  elm->UpdGlbVars( );  // Upd. global vars depending on the element type
}

// ============================= ReadPlStrainT6 ============================

void cElement :: ReadPlStrainT6(void)
{
  // Read the number of T6 elements.

  int n;
  if (!(in >> n) || (n < 1) || (n > NumElm))
  {
    cout << "Error in the input of the number of PLSTRAINT6 elements!\n";
    exit(0);
  }

  // Create and read each element.

  int id;
  cElement *elm;
  for (int i = 0; i < n; i++)
  {
    if (!(in >> id) || (id < 1) || (id > NumElm))
    {
      cout << "Error in the input of PLSTRAINT6 elements!\n";
      exit(0);
    }
    elm = new cElmParam(id, SHAPE_T6, PLANE_STRAIN);
    elm->Read( );
  }
  elm->UpdGlbVars( );  // Upd. global vars depending on the element type
}

// =========================== ReadPlStrainT6TL ============================

void cElement :: ReadPlStrainT6TL(void)
{
  // Read the number of T6 elements.

  int n;
  if (!(in >> n) || (n < 1) || (n > NumElm))
  {
    cout << "Error in the input of the number of PLSTRAINT6 elements!\n";
    exit(0);
  }

  // Create and read each element.

  int id;
  cElement *elm;
  for (int i = 0; i < n; i++)
  {
    if (!(in >> id) || (id < 1) || (id > NumElm))
    {
      cout << "Error in the input of PLSTRAINT6 elements!\n";
      exit(0);
    }
    elm = new cElmParamTL(id, SHAPE_T6, PLANE_STRAIN);
    elm->Read( );
  }
  elm->UpdGlbVars( );  // Upd. global vars depending on the element type
}

// =========================== ReadPlStrainBezTri ==========================

void cElement :: ReadPlStrainBezTri(void)
{
  cElmParam :: ReadElmFEM("Plane Strain Bezier Triangle",
		           SHAPE_BEZTRI_PLSURF, PLANE_STRAIN);
}

// =========================== ReadPlStrainBezTriTL ========================

void cElement :: ReadPlStrainBezTriTL(void)
{
  cElmParam :: ReadElmFEM("Plane Strain Bezier Triangle",
		           SHAPE_BEZTRI_PLSURF, PLANE_STRAIN,true);
}

// =========================== ReadPlStrainBezSurf =========================

void cElement :: ReadPlStrainBezSurf(void)
{
  cElmParam :: ReadElmFEM("Plane Strain Bezier Surface",
		          SHAPE_BEZ_PLSURF, PLANE_STRAIN);
}

// =========================== ReadPlStrainBezSurfTL =======================

void cElement :: ReadPlStrainBezSurfTL(void)
{
  cElmParam :: ReadElmFEM("Plane Strain Bezier Surface",
		          SHAPE_BEZ_PLSURF, PLANE_STRAIN,true);

}

// =========================== ReadPlStrainBsp =============================

void cElement :: ReadPlStrainBsp(void)
{
  cElmParam :: ReadElmBSP("Plane Strain B-Spline Surface",
		          SHAPE_BSP_PLSURF, PLANE_STRAIN);
}

// =========================== ReadPlStrainBspTL ===========================

void cElement :: ReadPlStrainBspTL(void)
{
  cElmParam :: ReadElmBSP("Plane Strain B-Spline Surface",
		          SHAPE_BSP_PLSURF, PLANE_STRAIN,true);
}

// ============================= ReadAxiSymL6Inf ===========================

void cElement :: ReadAxiSymL6Inf(void)
{
  // Read the number of L6INF elements.

  int n;
  if (!(in >> n) || (n < 1) || (n > NumElm))
  {
    cout << "Error in the input of the number of AXISYML6INF elements!\n";
    exit(0);
  }

  // Create and read each element.

  int id;
  cElement *elm;
  for (int i = 0; i < n; i++)
  {
    if (!(in >> id) || (id < 1) || (id > NumElm))
    {
      cout << "Error in the input of AXISYML6INF element " << id << "!\n";
      exit(0);
    }
    elm = new cElmParam(id, SHAPE_L6_INF, AXISYMMETRIC);
    elm->Read( );
  }
  elm->UpdGlbVars( );  // Upd. global vars depending on the element type
}

// ============================== ReadAxiSymQ4 =============================

void cElement :: ReadAxiSymQ4(void)
{
  // Read the number of Q4 elements.

  int n;
  if (!(in >> n) || (n < 1) || (n > NumElm))
  {
    cout << "Error in the input of the number of AXISYMQ4 elements!\n";
    exit(0);
  }

  // Create and read each element.

  int id;
  cElement *elm;
  for (int i = 0; i < n; i++)
  {
    if (!(in >> id) || (id < 1) || (id > NumElm))
    {
      cout << "Error in the input of AXISYMQ4 elements!\n";
      exit(0);
    }
    elm = new cElmParam(id, SHAPE_Q4, AXISYMMETRIC);
    elm->Read( );
  }
  elm->UpdGlbVars( );  // Upd. global vars depending on the element type
}

// =========================== ReadAxiSymQ4TL =============================

void cElement :: ReadAxiSymQ4TL(void)
{
  // Read the number of Q4 elements.

  int n;
  if (!(in >> n) || (n < 1) || (n > NumElm))
  {
    cout << "Error in the input of the number of AXISYMQ4 elements!\n";
    exit(0);
  }

  // Create and read each element.

  int id;
  cElement *elm;
  for (int i = 0; i < n; i++)
  {
    if (!(in >> id) || (id < 1) || (id > NumElm))
    {
      cout << "Error in the input of AXISYMQ4 elements!\n";
      exit(0);
    }
    elm = new cElmParamTL(id, SHAPE_Q4, AXISYMMETRIC);
    elm->Read( );
  }
  elm->UpdGlbVars( );  // Upd. global vars depending on the element type
}

// ============================== ReadAxiSymQ8 =============================

void cElement :: ReadAxiSymQ8(void)
{
  // Read the number of Q8 elements.

  int n;
  if (!(in >> n) || (n < 1) || (n > NumElm))
  {
    cout << "Error in the input of the number of AXISYMQ8 elements!\n";
    exit(0);
  }

  // Create and read each element.

  int id;
  cElement *elm;
  for (int i = 0; i < n; i++)
  {
    if (!(in >> id) || (id < 1) || (id > NumElm))
    {
      cout << "Error in the input of AXISYMQ8 elements!\n";
      exit(0);
    }
    elm = new cElmParam(id, SHAPE_Q8, AXISYMMETRIC);
    elm->Read( );
  }
  elm->UpdGlbVars( );  // Upd. global vars depending on the element type
}

// ============================ ReadAxiSymQ8TL =============================

void cElement :: ReadAxiSymQ8TL(void)
{
  // Read the number of Q8 elements.

  int n;
  if (!(in >> n) || (n < 1) || (n > NumElm))
  {
    cout << "Error in the input of the number of AXISYMQ8 elements!\n";
    exit(0);
  }

  // Create and read each element.

  int id;
  cElement *elm;
  for (int i = 0; i < n; i++)
  {
    if (!(in >> id) || (id < 1) || (id > NumElm))
    {
      cout << "Error in the input of AXISYMQ8 elements!\n";
      exit(0);
    }
    elm = new cElmParamTL(id, SHAPE_Q8, AXISYMMETRIC);
    elm->Read( );
  }
  elm->UpdGlbVars( );  // Upd. global vars depending on the element type
}

// ============================== ReadAxiSymQ9 =============================

void cElement :: ReadAxiSymQ9(void)
{
  // Read the number of Q9 elements.

  int n;
  if (!(in >> n) || (n < 1) || (n > NumElm))
  {
    cout << "Error in the input of the number of AXISYMQ9 elements!\n";
    exit(0);
  }

  // Create and read each element.

  int id;
  cElement *elm;
  for (int i = 0; i < n; i++)
  {
    if (!(in >> id) || (id < 1) || (id > NumElm))
    {
      cout << "Error in the input of AXISYMQ9 elements!\n";
      exit(0);
    }
    elm = new cElmParam(id, SHAPE_Q9, AXISYMMETRIC);
    elm->Read( );
  }
  elm->UpdGlbVars( );  // Upd. global vars depending on the element type
}

// ============================ ReadAxiSymQ9TL =============================

void cElement :: ReadAxiSymQ9TL(void)
{
  // Read the number of Q9 elements.

  int n;
  if (!(in >> n) || (n < 1) || (n > NumElm))
  {
    cout << "Error in the input of the number of AXISYMQ9 elements!\n";
    exit(0);
  }

  // Create and read each element.

  int id;
  cElement *elm;
  for (int i = 0; i < n; i++)
  {
    if (!(in >> id) || (id < 1) || (id > NumElm))
    {
      cout << "Error in the input of AXISYMQ9 elements!\n";
      exit(0);
    }
    elm = new cElmParamTL(id, SHAPE_Q9, AXISYMMETRIC);
    elm->Read( );
  }
  elm->UpdGlbVars( );  // Upd. global vars depending on the element type
}

// ============================== ReadAxiSymT3 =============================

void cElement :: ReadAxiSymT3(void)
{
  // Read the number of T3 elements.

  int n;
  if (!(in >> n) || (n < 1) || (n > NumElm))
  {
    cout << "Error in the input of the number of AXISYMT3 elements!\n";
    exit(0);
  }

  // Create and read each element.

  int id;
  cElement *elm;
  for (int i = 0; i < n; i++)
  {
    if (!(in >> id) || (id < 1) || (id > NumElm))
    {
      cout << "Error in the input of AXISYMT3 elements!\n";
      exit(0);
    }
    elm = new cElmParam(id, SHAPE_T3, AXISYMMETRIC);
    elm->Read( );
  }
  elm->UpdGlbVars( );  // Upd. global vars depending on the element type
}

// ============================= ReadAxiSymT3TL =============================

void cElement :: ReadAxiSymT3TL(void)
{
  // Read the number of T3 elements.

  int n;
  if (!(in >> n) || (n < 1) || (n > NumElm))
  {
    cout << "Error in the input of the number of AXISYMT3 elements!\n";
    exit(0);
  }

  // Create and read each element.

  int id;
  cElement *elm;
  for (int i = 0; i < n; i++)
  {
    if (!(in >> id) || (id < 1) || (id > NumElm))
    {
      cout << "Error in the input of AXISYMT3 elements!\n";
      exit(0);
    }
    elm = new cElmParamTL(id, SHAPE_T3, AXISYMMETRIC);
    elm->Read( );
  }
  elm->UpdGlbVars( );  // Upd. global vars depending on the element type
}

// ============================== ReadAxiSymT6 =============================

void cElement :: ReadAxiSymT6(void)
{
  // Read the number of T6 elements.

  int n;
  if (!(in >> n) || (n < 1) || (n > NumElm))
  {
    cout << "Error in the input of the number of AXISYMT6 elements!\n";
    exit(0);
  }

  // Create and read each element.

  int id;
  cElement *elm;
  for (int i = 0; i < n; i++)
  {
    if (!(in >> id) || (id < 1) || (id > NumElm))
    {
      cout << "Error in the input of AXISYMT6 elements!\n";
      exit(0);
    }
    elm = new cElmParam(id, SHAPE_T6, AXISYMMETRIC);
    elm->Read( );
  }
  elm->UpdGlbVars( );  // Upd. global vars depending on the element type
}

// ============================ ReadAxiSymT6TL =============================

void cElement :: ReadAxiSymT6TL(void)
{
  // Read the number of T6 elements.

  int n;
  if (!(in >> n) || (n < 1) || (n > NumElm))
  {
    cout << "Error in the input of the number of AXISYMT6 elements!\n";
    exit(0);
  }

  // Create and read each element.

  int id;
  cElement *elm;
  for (int i = 0; i < n; i++)
  {
    if (!(in >> id) || (id < 1) || (id > NumElm))
    {
      cout << "Error in the input of AXISYMT6 elements!\n";
      exit(0);
    }
    elm = new cElmParamTL(id, SHAPE_T6, AXISYMMETRIC);
    elm->Read( );
  }
  elm->UpdGlbVars( );  // Upd. global vars depending on the element type
}

// ============================= ReadPlThermT3 =============================

void cElement :: ReadPlThermT3(void)
{
  cElmParam :: ReadElmFEM("Plane Thermal T3",
		           SHAPE_T3, PLANE_HEAT_TRANSFER);
}

// ============================= ReadPlThermT6 =============================

void cElement :: ReadPlThermT6(void)
{
  cElmParam :: ReadElmFEM("Plane Thermal T6",
		           SHAPE_T6, PLANE_HEAT_TRANSFER);
}

// ============================= ReadPlThermQ4 =============================

void cElement :: ReadPlThermQ4(void)
{
  cElmParam :: ReadElmFEM("Plane Thermal Q4",
		           SHAPE_Q4, PLANE_HEAT_TRANSFER);
}

// ============================= ReadPlThermQ8 =============================

void cElement :: ReadPlThermQ8(void)
{
  cElmParam :: ReadElmFEM("Plane Thermal Q8",
		           SHAPE_Q8, PLANE_HEAT_TRANSFER);
}

// =========================== ReadPlThermBezTri ===========================

void cElement :: ReadPlThermBezTri(void)
{
  cElmParam :: ReadElmFEM("Plane Thermal Bezier Triangle",
		           SHAPE_BEZTRI_PLSURF, PLANE_HEAT_TRANSFER);
}

// =========================== ReadPlThermBezSurf ==========================

void cElement :: ReadPlThermBezSurf(void)
{
  cElmParam :: ReadElmFEM("Plane Thermal Bezier Surface",
		           SHAPE_BEZ_PLSURF, PLANE_HEAT_TRANSFER);
}

// =============================== ReadBrick8 ==============================

void cElement :: ReadBrick8(void)
{
  // Read the number of BRICK8 elements.

  int n;
  if (!(in >> n) || (n < 1) || (n > NumElm))
  {
    cout << "Error in the input of the number of BRICK8 elements!\n";
    exit(0);
  }

  // Create and read each element.

  int id;
  cElement *elm;
  for (int i = 0; i < n; i++)
  {
    if (!(in >> id) || (id < 1) || (id > NumElm))
    {
      cout << "Error in the input of BRICK8 elements!\n";
      exit(0);
    }
    elm = new cElmParam(id, SHAPE_BRICK8, SOLID);
    elm->Read( );
  }
  elm->UpdGlbVars( );  // Upd. global vars depending on the element type
}

// ============================= ReadBrick8TL ==============================

void cElement :: ReadBrick8TL(void)
{
  // Read the number of BRICK8TL elements.

  int n;
  if (!(in >> n) || (n < 1) || (n > NumElm))
  {
    cout << "Error in the input of the number of BRICK8TL elements!\n";
    exit(0);
  }

  // Create and read each element.

  int id;
  cElement *elm;
  for (int i = 0; i < n; i++)
  {
    if (!(in >> id) || (id < 1) || (id > NumElm))
    {
      cout << "Error in the input of BRICK8TL elements!\n";
      exit(0);
    }
    elm = new cElmParamTL(id, SHAPE_BRICK8, SOLID);
    elm->Read( );
  }
  elm->UpdGlbVars( );  // Upd. global vars depending on the element type
}

// =============================== ReadBrick20 =============================

void cElement :: ReadBrick20(void)
{
  // Read the number of BRICK20 elements.

  int n;
  if (!(in >> n) || (n < 1) || (n > NumElm))
  {
    cout << "Error in the input of the number of BRICK20 elements!\n";
    exit(0);
  }

  // Create and read each element.

  int id;
  cElement *elm;
  for (int i = 0; i < n; i++)
  {
    if (!(in >> id) || (id < 1) || (id > NumElm))
    {
      cout << "Error in the input of BRICK20 elements!\n";
      exit(0);
    }
    elm = new cElmParam(id, SHAPE_BRICK20, SOLID);
    elm->Read( );
  }
  elm->UpdGlbVars( );  // Upd. global vars depending on the element type
}

// ============================= ReadBrick20TL ==============================

void cElement :: ReadBrick20TL(void)
{
  // Read the number of BRICK20TL elements.

  int n;
  if (!(in >> n) || (n < 1) || (n > NumElm))
  {
    cout << "Error in the input of the number of BRICK20TL elements!\n";
    exit(0);
  }

  // Create and read each element.

  int id;
  cElement *elm;
  for (int i = 0; i < n; i++)
  {
    if (!(in >> id) || (id < 1) || (id > NumElm))
    {
      cout << "Error in the input of BRICK20TL elements!\n";
      exit(0);
    }
    elm = new cElmParamTL(id, SHAPE_BRICK20, SOLID);
    elm->Read( );
  }
  elm->UpdGlbVars( );  // Upd. global vars depending on the element type
}

// ================================ ReadTet4 ===============================

void cElement :: ReadTet4(void)
{
  // Read the number of TET4 elements.

  int n;
  if (!(in >> n) || (n < 1) || (n > NumElm))
  {
    cout << "Error in the input of the number of TET4 elements!\n";
    exit(0);
  }

  // Create and read each element.

  int id;
  cElement *elm;
  for (int i = 0; i < n; i++)
  {
    if (!(in >> id) || (id < 1) || (id > NumElm))
    {
      cout << "Error in the input of TET4 elements!\n";
      exit(0);
    }
    elm = new cElmParam(id, SHAPE_TET4, SOLID);
    elm->Read( );
  }
  elm->UpdGlbVars( );  // Upd. global vars depending on the element type
}

// =============================== ReadTet4TL ==============================

void cElement :: ReadTet4TL(void)
{
  // Read the number of TET4TL elements.

  int n;
  if (!(in >> n) || (n < 1) || (n > NumElm))
  {
    cout << "Error in the input of the number of TET4TL elements!\n";
    exit(0);
  }

  // Create and read each element.

  int id;
  cElement *elm;
  for (int i = 0; i < n; i++)
  {
    if (!(in >> id) || (id < 1) || (id > NumElm))
    {
      cout << "Error in the input of TET4TL elements!\n";
      exit(0);
    }
    elm = new cElmParamTL(id, SHAPE_TET4, SOLID);
    elm->Read( );
  }
  elm->UpdGlbVars( );  // Upd. global vars depending on the element type
}

// ================================ ReadTet10 ==============================

void cElement :: ReadTet10(void)
{
  // Read the number of TET10 elements.

  int n;
  if (!(in >> n) || (n < 1) || (n > NumElm))
  {
    cout << "Error in the input of the number of TET10 elements!\n";
    exit(0);
  }

  // Create and read each element.

  int id;
  cElement *elm;
  for (int i = 0; i < n; i++)
  {
    if (!(in >> id) || (id < 1) || (id > NumElm))
    {
      cout << "Error in the input of TET10 elements!\n";
      exit(0);
    }
    elm = new cElmParam(id, SHAPE_TET10, SOLID);
    elm->Read( );
  }
  elm->UpdGlbVars( );  // Upd. global vars depending on the element type
}

// =============================== ReadTet10TL =============================

void cElement :: ReadTet10TL(void)
{
  // Read the number of TET10TL elements.

  int n;
  if (!(in >> n) || (n < 1) || (n > NumElm))
  {
    cout << "Error in the input of the number of TET10TL elements!\n";
    exit(0);
  }

  // Create and read each element.

  int id;
  cElement *elm;
  for (int i = 0; i < n; i++)
  {
    if (!(in >> id) || (id < 1) || (id > NumElm))
    {
      cout << "Error in the input of TET10TL elements!\n";
      exit(0);
    }
    elm = new cElmParamTL(id, SHAPE_TET10, SOLID);
    elm->Read( );
  }
  elm->UpdGlbVars( );  // Upd. global vars depending on the element type
}

// ============================= ReadSolidBsp ==============================

void cElement :: ReadSolidBsp(void)
{
  cElmParam :: ReadElmBSP("Solid B-Spline",SHAPE_BSP_SOLID,SOLID);
}

// ============================= ReadSolidBspTL ============================

void cElement :: ReadSolidBspTL(void)
{
  cElmParam :: ReadElmBSP("Solid B-Spline",SHAPE_BSP_SOLID,SOLID,true);
}

// ============================== ReadThickPltT3 ===========================

void cElement :: ReadThickPltT3(void)
{
  cElmParam :: ReadElmFEM("Thick Plate T3",SHAPE_T3, THICK_PLATE);
}

// ============================== ReadThickPltT6 ===========================

void cElement :: ReadThickPltT6(void)
{
  cElmParam :: ReadElmFEM("Thick Plate T6",SHAPE_T6, THICK_PLATE);
}

// ============================== ReadThickPltQ4 ===========================

void cElement :: ReadThickPltQ4(void)
{
  // Read the number of THICKPLTQ4 elements.

  int n;
  if (!(in >> n) || (n < 1) || (n > NumElm))
  {
    cout << "Error in the input of the number of THICKPLTQ4 elements!\n";
    exit(0);
  }

  // Create and read each element.

  int id;
  cElement *elm;
  for (int i = 0; i < n; i++)
  {
    if (!(in >> id) || (id < 1) || (id > NumElm))
    {
      cout << "Error in the input of THICKPLTQ4 elements!\n";
      exit(0);
    }
    elm = new cElmParam(id, SHAPE_Q4, THICK_PLATE);
    elm->Read( );
  }
  elm->UpdGlbVars( );  // Upd. global vars depending on the element type
}

// ============================== ReadThickPltQ8 ===========================

void cElement :: ReadThickPltQ8(void)
{
  // Read the number of THICKPLTQ8 elements.

  int n;
  if (!(in >> n) || (n < 1) || (n > NumElm))
  {
    cout << "Error in the input of the number of THICKPLTQ8 elements!\n";
    exit(0);
  }

  // Create and read each element.

  int id;
  cElement *elm;
  for (int i = 0; i < n; i++)
  {
    if (!(in >> id) || (id < 1) || (id > NumElm))
    {
      cout << "Error in the input of THICKPLTQ8 elements!\n";
      exit(0);
    }
    elm = new cElmParam(id, SHAPE_Q8, THICK_PLATE);
    elm->Read( );
  }
  elm->UpdGlbVars( );  // Upd. global vars depending on the element type
}

// ============================= ReadShallowShellQ8 ========================

void cElement :: ReadShallowShellQ8(void)
{
  // Read the number of SHALLOWSHELLQ8 elements.

  int n;
  if (!(in >> n) || (n < 1) || (n > NumElm))
  {
    cout << "Error in the input of the number of SHALLOWSHELLQ8 elements!\n";
    exit(0);
  }

  // Create and read each element.

  int id;
  cElement *elm;
  for (int i = 0; i < n; i++)
  {
    if (!(in >> id) || (id < 1) || (id > NumElm))
    {
      cout << "Error in the input of SHALLOWSHELLQ8 elements!\n";
      exit(0);
    }
    elm = new cElmParam(id, SHAPE_Q8, SHALLOW_SHELL);
    elm->Read( );
  }
  elm->UpdGlbVars( );  // Upd. global vars depending on the element type
}

// ========================= ReadShallowShellQ8TL ==========================

void cElement :: ReadShallowShellQ8TL(void)
{
  // Read the number of SHALLOWSHELLQ8TL elements.

  int n;
  if (!(in >> n) || (n < 1) || (n > NumElm))
  {
    cout << "Error in the input of the number of SHALLOWSHELLQ8TL elements!\n";
    exit(0);
  }

  // Create and read each element.

  int id;
  cElement *elm;
  for (int i = 0; i < n; i++)
  {
    if (!(in >> id) || (id < 1) || (id > NumElm))
    {
      cout << "Error in the input of SHALLOWSHELLQ8TL elements!\n";
      exit(0);
    }
    elm = new cElmParamTL(id, SHAPE_Q8, SHALLOW_SHELL);
    elm->Read( );
  }
  elm->UpdGlbVars( );  // Upd. global vars depending on the element type
}

// ============================= ReadDonnellShellQ8 ========================

void cElement :: ReadDonnellShellQ8(void)
{
  // Read the number of DonnellShellQ8 elements.

  int n;
  if (!(in >> n) || (n < 1) || (n > NumElm))
  {
    cout << "Error in the input of the number of DONNELLSHELLQ8 elements!\n";
    exit(0);
  }

  // Create and read each element.

  int id;
  cElement *elm;
  for (int i = 0; i < n; i++)
  {
    if (!(in >> id) || (id < 1) || (id > NumElm))
    {
      cout << "Error in the input of DONNELLSHELLQ8 elements!\n";
      exit(0);
    }
    elm = new cElmParam(id, SHAPE_Q8, DONNELL_SHELL);
    elm->Read( );
  }
  elm->UpdGlbVars( );  // Upd. global vars depending on the element type
}

// ========================= ReadDonnellShellQ8TL ==========================

void cElement :: ReadDonnellShellQ8TL(void)
{
  // Read the number of DONNELLSHELLQ8TL elements.

  int n;
  if (!(in >> n) || (n < 1) || (n > NumElm))
  {
    cout << "Error in the input of the number of DONNELLSHELLQ8TL elements!\n";
    exit(0);
  }

  // Create and read each element.

  int id;
  cElement *elm;
  for (int i = 0; i < n; i++)
  {
    if (!(in >> id) || (id < 1) || (id > NumElm))
    {
      cout << "Error in the input of DONNELLSHELLQ8TL elements!\n";
      exit(0);
    }
    elm = new cElmParamTL(id, SHAPE_Q8, DONNELL_SHELL);
    elm->Read( );
  }
  elm->UpdGlbVars( );  // Upd. global vars depending on the element type
}

// ========================== ReadDonnellShellBsp ==========================

void cElement :: ReadDonnellShellBsp(void)
{
  cElmParam :: ReadElmBSP("Donnell Shell B-Spline Surface",
		          SHAPE_BSP_SURF, DONNELL_SHELL);
}

// ========================= ReadDonnellShellBspTL =========================

void cElement :: ReadDonnellShellBspTL(void)
{
  cElmParam :: ReadElmBSP("Donnell Shell B-Spline Surface",
		          SHAPE_BSP_SURF, DONNELL_SHELL,true);
}


// ========================== ReadShallowShellT6 ===========================

void cElement :: ReadShallowShellT6(void)
{
  cElmParam :: ReadElmFEM("Shallow Shell T6",
		          SHAPE_T6, SHALLOW_SHELL);
}

// ========================== ReadShallowShellT6TL =========================

void cElement :: ReadShallowShellT6TL(void)
{
  cElmParam :: ReadElmFEM("Shallow Shell T6",
		          SHAPE_T6, SHALLOW_SHELL,true);
}

// ========================== ReadShallowShellBezTrian =====================

void cElement :: ReadShallowShellBezTrian(void)
{
  cElmParam :: ReadElmFEM("Shallow Shell Bezier Triangle",
		          SHAPE_BEZTRI_SURF, SHALLOW_SHELL);
}

// ========================== ReadShallowShellBezTrianTL ===================

void cElement :: ReadShallowShellBezTrianTL(void)
{
  cElmParam :: ReadElmFEM("Shallow Shell Bezier Triangle",
		          SHAPE_BEZTRI_SURF,SHALLOW_SHELL,true);
}

// ========================== ReadShallowShellBezSurf ======================

void cElement :: ReadShallowShellBezSurf(void)
{
  cElmParam :: ReadElmFEM("Shallow Shell Bezier Surface",
		          SHAPE_BEZ_SURF, SHALLOW_SHELL);
}

// ========================== ReadShallowShellBezSurfTL ====================

void cElement :: ReadShallowShellBezSurfTL(void)
{
  cElmParam :: ReadElmFEM("Shallow Shell Bezier Surface",
		          SHAPE_BEZ_SURF, SHALLOW_SHELL,true);
}

// ============================ ReadHSDTPlateBsp ===========================

void cElement :: ReadHSDTPlateBsp(void)
{
  cElmParamC1 :: ReadElmBSP("HSDT Plate B-Spline Surface",
                  SHAPE_BSP_PLSURF, HSDT_PLATE);
}

// ========================== ReadHSDTPlateBspTL ===========================

void cElement :: ReadHSDTPlateBspTL(void)
{
  cElmParamC1TL :: ReadElmBSP("HSDT Plate B-Spline Surface",
                  SHAPE_BSP_PLSURF, HSDT_PLATE, true);
}

// ========================== ReadShallowShellBsp ==========================

void cElement :: ReadShallowShellBsp(void)
{
  cElmParam :: ReadElmBSP("Shallow Shell B-Spline Surface",
		          SHAPE_BSP_SURF, SHALLOW_SHELL);
}

// ========================= ReadShallowShellBspTL =========================

void cElement :: ReadShallowShellBspTL(void)
{
  cElmParam :: ReadElmBSP("Shallow Shell B-Spline Surface",
		          SHAPE_BSP_SURF, SHALLOW_SHELL,true);
}

// ============================= ReadShellT3 ===============================

void cElement :: ReadShellT3(void)
{
  cElmParam :: ReadElm<cElmDegShell>("Shell T3", SHAPE_T3_SHELL, SHELL);
}

// ============================= ReadShellT6 ===============================

void cElement :: ReadShellT6(void)
{
  cElmParam :: ReadElm<cElmDegShell>("Shell T6", SHAPE_T6_SHELL, SHELL);
}

// ============================= ReadShellQ8 ===============================

void cElement :: ReadShellQ8(void)
{
  cElmParam :: ReadElm<cElmDegShell>("Shell Q8", SHAPE_Q8_SHELL, SHELL);
}

// ============================= ReadShellQ8TL =============================

void cElement :: ReadShellQ8TL(void)
{
  cElmParam :: ReadElm<cElmDegShellTL>("Total Lagrangian Shell Q8", SHAPE_Q8_SHELL, SHELL);
}

// ========================== ReadShellQ9 ==================================

void cElement :: ReadShellQ9(void)
{
  cElmParam :: ReadElm<cElmDegShell>("Shell Q9", SHAPE_Q9_SHELL, SHELL);
}

// ======================== ReadShellBezTrian ==============================

void cElement :: ReadShellBezTrian(void)
{
  cElmParam :: ReadElm<cElmDegShell>("Shell Bezier Triangle", SHAPE_BEZTRI_SHELL, SHELL);
}

// ======================== ReadShellBezTrianTL ============================

void cElement :: ReadShellBezTrianTL(void)
{
  cElmParam :: ReadElm<cElmDegShellTL>("Total Lagrangian Shell Bezier Triangle", SHAPE_BEZTRI_SHELL, SHELL);
}

// ========================== ReadShellBez =================================

void cElement :: ReadShellBez(void)
{
  cElmParam :: ReadElm<cElmDegShell>("Shell Bez", SHAPE_BEZ_SHELL, SHELL);
}

// ========================== ReadShellBezTL ===============================

void cElement :: ReadShellBezTL(void)
{
  cElmParam :: ReadElm<cElmDegShellTL>("Total Lagrangian Shell Bez", SHAPE_BEZ_SHELL, SHELL);
}

// ========================== ReadShellBsp =================================

void cElement :: ReadShellBsp(void)
{
  cElmParam :: ReadElm<cElmDegShell>("Shell BSpline", SHAPE_BSP_SHELL, SHELL);
}

// ========================== ReadShellBspTL ===============================

void cElement :: ReadShellBspTL(void)
{
  cElmParam :: ReadElm<cElmDegShellTL>("Total Lagrangian Shell BSpline", SHAPE_BSP_SHELL, SHELL);
}

// ============================= ReadInterfL2 ==============================

void cElement :: ReadInterfL2(void)
{
  // Read the number of INTERFACE.L2 elements.

  int n;
  if (!(in >> n) || (n < 1) || (n > NumElm))
  {
    cout << "Error in the input of the number of INTERFACE.L2 elements!\n";
    exit(0);
  }

  // Create and read each element.

  int id;
  cElement *elm;
  for (int i = 0; i < n; i++)
  {
    if (!(in >> id) || (id < 1) || (id > NumElm))
    {
      cout << "Error in the input of INTERFACE.L2 elements!\n";
      exit(0);
    }
    elm = new cElmInterf(id, SHAPE_INTERF_L2, INTERFACE_2D);
    elm->Read( );
  }
  elm->UpdGlbVars( );  // Upd. global vars depending on the element type
}

// ============================= ReadInterfL3 ==============================

void cElement :: ReadInterfL3(void)
{
  // Read the number of INTERFACE.L3 elements.

  int n;
  if (!(in >> n) || (n < 1) || (n > NumElm))
  {
    cout << "Error in the input of the number of INTERFACE.L3 elements!\n";
    exit(0);
  }

  // Create and read each element.

  int id;
  cElement *elm;
  for (int i = 0; i < n; i++)
  {
    if (!(in >> id) || (id < 1) || (id > NumElm))
    {
      cout << "Error in the input of INTERFACE.L3 elements!\n";
      exit(0);
    }
    elm = new cElmInterf(id, SHAPE_INTERF_L3, INTERFACE_2D);
    elm->Read( );
  }
  elm->UpdGlbVars( );  // Upd. global vars depending on the element type
}

// ============================= ReadInterfQ4 ==============================

void cElement :: ReadInterfQ4(void)
{
  // Read the number of INTERFACE.Q4 elements.

  int n;
  if (!(in >> n) || (n < 1) || (n > NumElm))
  {
    cout << "Error in the input of the number of INTERFACE.Q4 elements!\n";
    exit(0);
  }

  // Create and read each element.

  int id;
  cElement *elm;
  for (int i = 0; i < n; i++)
  {
    if (!(in >> id) || (id < 1) || (id > NumElm))
    {
      cout << "Error in the input of INTERFACE.Q4 elements!\n";
      exit(0);
    }
    elm = new cElmInterf(id, SHAPE_INTERF_Q4, INTERFACE_3D);
    elm->Read( );
  }
  elm->UpdGlbVars( );  // Upd. global vars depending on the element type
}

// ============================= ReadInterfQ8 ==============================

void cElement :: ReadInterfQ8(void)
{
  // Read the number of INTERFACE.Q8 elements.

  int n;
  if (!(in >> n) || (n < 1) || (n > NumElm))
  {
    cout << "Error in the input of the number of INTERFACE.Q8 elements!\n";
    exit(0);
  }

  // Create and read each element.

  int id;
  cElement *elm;
  for (int i = 0; i < n; i++)
  {
    if (!(in >> id) || (id < 1) || (id > NumElm))
    {
      cout << "Error in the input of INTERFACE.Q8 elements!\n";
      exit(0);
    }
    elm = new cElmInterf(id, SHAPE_INTERF_Q8, INTERFACE_3D);
    elm->Read( );
  }
  elm->UpdGlbVars( );  // Upd. global vars depending on the element type
}

// ========================== ReadInterfLineBsp ============================

void cElement :: ReadInterfLineBsp( )
{
/*
  // Read the number of INTERFACE.LINE.IGA elements.

  int n;
  if (!(in >> n) || (n < 1) || (n > NumElm))
  {
    cout << "Error in the input of the number of INTERFACE.LINE.BSP  elements!\n";
    exit(0);
  }

  // Create and read each element.

  int id;
  cElement *elm;
  for (int i = 0; i < n; i++)
  {
    if (!(in >> id) || (id < 1) || (id > NumElm))
    {

  int id;
  cElement *elm;
  for (int i = 0; i < n; i++)
  {
    if (!(in >> id) || (id < 1) || (id > NumElm))
    {
      cout << "Error in the input of SHALLOWSHELLQ8TL elements!\n";
      exit(0);
    }
    elm = new cElmParamTL(id, SHAPE_Q8, SHALLOW_SHELL);
    elm->Read( );
  }
  elm->UpdGlbVars( );  // Upd. global vars depending on the element type
  */
}


// ========================== ReadElemTemp ==============================
void cElement :: ReadElemTemp(void)
{
  // Read the number of element temperatures

  if (!(in >> NumElemTemp) || (NumElemTemp < 1))
  {
    cout << "Error in the input of the number of element temperatures!\n";
    exit(0);
  }

  if (NumElemTemp <= 0) return;

  ElemTemp = new sTempData[NumElemTemp];

  int label;
  double ti,ts;
  //cElement *elm;

  for (int i = 0; i < NumElemTemp; i++)
  {
    if (!(in >> label) || !(in >> ti) || !(in >> ts))
    {
      cout << "\n Error on reading element temperatures " << i+1 << "!\n";
      exit(0);
    }

    ElemTemp[i].label = label;
    ElemTemp[i].tinf = ti;
    ElemTemp[i].tsup = ts;

    //cout << "label = " << ElemTemp[i].label << " - " << "tinf = " << ElemTemp[i].tinf << " - " << "tsup = " << ElemTemp[i].tsup << " \n";

    //sTempData *temp;
    //temp = elm->GetTemp();

   //  if (!temp)
   // {
   //   cout << "\n Error on reading element temperatures %d !!!\n\n" << i+1 << "!\n";
   //   exit(0);
   // }

    }

  // Sort the tempetature data according to the element label.

  qsort(ElemTemp, NumElemTemp, sizeof(sTempData), sTempDataCmp);

} // End of ReadElemTemp


// ================================ GetTemp ================================
sTempData *cElement :: GetTemp(void)
{
  for(int i = 0; i < NumElemTemp;++i)
 {
  if (ElemTemp[i].label == Label) return &(ElemTemp[i]);
 }
  return (0);
} // End of GetTemp


// =========================== ExistThermalLoad ============================

int cElement :: ExistThermalLoad(void)
{
  if (NumElemTemp > 0) return(1);

  return(0);

} // End of ExistThermalLoad


// ================================ Destroy ================================

void cElement :: Destroy(void)
{
  // Destroy each element.

  for (int i = 0; i < NumElm; i++) delete VecElm[i];

  // Release the array of elements.

  delete []VecElm;
  VecElm = 0;

  // Release the array of scalar labels.

  delete []VecSclLab;

  // Release the array of beam orientation.

  delete []VecBeamOrt;
}

// ================================ FindElm ================================

cElement *cElement :: FindElm(int label)
{
  if (label < 1 || label > NumElm) return(0);

  return(VecElm[label-1]);
}

// ============================= GetVecSclLab ==============================

void cElement :: GetVecSclLab(int *label)
{
  for (int i = 0; i < NumSclLab; i++) label[i] = VecSclLab[i];
}

// ============================= GetVecNodSclLab ===========================

void cElement :: GetVecNodSclLab(int *label)
{
  for (int i = 0; i < NumNodSclLab; i++) label[i] = VecNodSclLab[i];
}

// =============================== cElement ================================

cElement :: cElement(int id)
{
  Label = id;           // Store the element label
  VecElm[id-1] = this;  // Add the element to the element vector
}

// ============================== ~cElement ================================

cElement :: ~cElement(void)
{
  delete Shape;
}

// ============================ GetNumElmDof ===============================

int cElement :: GetNumElmDof(void)
{
  return Shape->GetNumNode( )*GetNumDofNode( );
}

// ============================== GetElmDofs ===============================

void cElement :: GetElmDofs(int *dof)
{
  // Get element active directions.

  int dir[8] = {0,0,0,0,0,0,0,0};
  GetActDir(dir);

  // Get element equations from connected nodes.

  int k = 0;
  int nn = Shape->GetNumNode( );
  for (int i = 0; i < nn; i++)
  {
    cNode *node = Shape->GetNode(i);
    for (int j = 0; j < 8; j++) if (dir[j])
    {
      dof[k++] = node->GetDof(j);
    }
  }
}

// ============================== ActivateDof ==============================

void cElement :: ActivateDof(void)
{
  // Get element active directions.

  int dir[8] = {0,0,0,0,0,0,0,0};
  GetActDir(dir);

  // Activate the dofs of the element nodes.

  int nn = Shape->GetNumNode( );
  for (int i = 0; i < nn; i++)
  {
    cNode *node = Shape->GetNode(i);
    for (int j = 0; j < 8; j++)
    {
      int dof = node->GetDof(j) || dir[j];
      node->SetDof(j, dof);
    }
  }
}

// ============================== AddGlobVec ===============================

void cElement :: AddGlobVec(cVector &elmvec, cVector &glbvec)
{
  // Get element dofs.

  int ndof = GetNumElmDof( );
  int *dof = new int[ndof];
  GetElmDofs(dof);

  // Add element vector to the global vector.

  for (int i = 0; i < ndof; i++) if (dof[i])
  {
    glbvec[dof[i]-1] += elmvec[i];
  }

  // Release the dof array.

  delete []dof;
}

// ============================== AddGlobMat ===============================

void cElement :: AddGlobMat(cMatrix &elmmat, cSysMatrix *glbmat)
{
  // Get element dofs.

  int ndof = GetNumElmDof( );
  int *dof = new int[ndof];
  GetElmDofs(dof);

  // Add element matrix to the global matrix.

  for (int i = 0; i < ndof; i++) if (dof[i])
  {
    for (int j = 0; j < ndof; j++) if (dof[j])
    {
      // Add if element is below/at the main diagonal or the matrix is
      // not symmetric.

      if ((dof[i] >= dof[j]) || (glbmat->Symmetric( ) == 0))
      {
        glbmat->Add(dof[i]-1, dof[j]-1, elmmat[i][j]);
      }
    }
  }

  // Release the dof array.

  delete []dof;
}

// =============================== GlobToElm ===============================

void cElement :: GlobToElm(cVector &glbvec, cVector &elmvec)
{
  // Get element dofs.

  int ndof = GetNumElmDof( );
  int *dof = new int[ndof];
  GetElmDofs(dof);

  // Get element data from the global vector.

  for (int i = 0; i < ndof; i++)
  {
    if (dof[i])
      elmvec[i] = glbvec[dof[i]-1];
    else
      elmvec[i] = 0.0;
  }

  // Release the dof array.

  delete []dof;
}

// ============================== NodalDispl ===============================

void cElement :: NodalDispl(cVector &u)
{
  // Get element active directions.

  int dir[8] = {0,0,0,0,0,0,0,0};
  GetActDir(dir);

  // Get the displacements from connected nodes.

  int k = 0;
  int nn = Shape->GetNumNode( );
  for (int i = 0; i < nn; i++)
  {
    cNode *node = Shape->GetNode(i);
    for (int j = 0; j < 6; j++) if (dir[j])
    {
      u[k++] = node->GetDispl(j);
    }
  }
}

// ============================== NodalTemp ================================

void cElement :: NodalTemp(cVector &t)
{
  // Get the displacements from connected nodes.

  int k = 0;
  int nn = Shape->GetNumNode( );
  for (int i = 0; i < nn; i++)
    t[i] = Shape->GetNode(i)->GetTemp( );
}

// ============================== NodalDu ===================================

void cElement :: NodalDu(cVector &du, cVector &duelm)
{
  // Get Element active direction.

  int dir[8] = {0,0,0,0,0,0,0,0};
  GetActDir(dir);

  int k = 0;
  int nn = Shape->GetNumNode( );

  for (int i = 0; i < nn; i++)
  {
    cNode *node = Shape->GetNode(i);
    for (int j = 0; j < 8; j++) if (dir[j])
    {
      int dof = node->GetDof(j);
      if (dof)
        duelm[k++] = du[dof -1];
      else
        duelm[k++] = 0.0;
    }
  }
}

// ================================ HasSupp ================================

bool cElement :: HasSupp(void)
{
  int nn = Shape->GetNumNode( );
  for (int i = 0; i < nn; i++)
  {
    cNode *node = Shape->GetNode(i);
    if (node->HasSupp( )) return(true);
  }

  return(false);
}

// ============================= AddReactions ==============================

void cElement :: AddReactions(cVector &gelm, cMatrix &R)
{
  // Get element active directions.

  int dir[8] = {0,0,0,0,0,0,0,0};
  GetActDir(dir);

  // Add internal forces to nodal reactions.

  int k = 0;
  int nn = Shape->GetNumNode( );
  for (int i = 0; i < nn; i++)
  {
    cNode *node = Shape->GetNode(i);
    int idx = node->GetLabel( ) - 1;
    for (int j = 0; j < 8; j++) if (dir[j])
    {
      R(idx, j) += gelm[k++];
    }
  }
}

// -------------------------------------------------------------------------
// Protected methods:
//

// ============================== UpdGlbVars ===============================
//
// This function updates some static (global) variables, such as MaxNode,
// MaxDof, and VecSclLab, according to the element data.
//
void cElement :: UpdGlbVars(void)
{
  // Update the max number of nodes.

  int nn = max(Shape->GetNumShp( ),Shape->GetNumShp( ));
  if (nn > MaxNode) MaxNode = nn;

  // Update the max number of integration points.

  int np = GetNumIntPnt( );
  if (np > MaxIntPnt) MaxIntPnt = np;

  // Update the max number of dofs.

  int ndof = GetNumElmDof( );
  if (ndof > MaxDof) MaxDof = ndof;

  // Create the vector of scalar labels.

  if (NumSclLab == 0)
  {
    VecSclLab = new int[MAX_SCL_LAB];
  }

  // Get the vector of scalar labels of the current element.

  int nel = GetNumStrCmp( );
  int elmlab[MAX_SCL_LAB];
  GetStrLabels(elmlab);

  // Add the scalar labels of the current element to the global vector.

  for (int i = 0; i < nel; i++)
  {
    int found = 0;
    for (int l = 0; l < NumSclLab; l++) // Search the element label
    {
      if (elmlab[i] == VecSclLab[l])
      {
        found = 1;
        break;
      }
    }
    if (!found)  // New label => add to the global vector
    {
      VecSclLab[NumSclLab++] = elmlab[i];
    }
  }

  // Create the vector of element nodal scalar labels.

  if (NumNodSclLab == 0)
  {
    VecNodSclLab = new int[MAX_SCL_LAB];
  }

  // Get the vector of scalar labels of the current element.

  int neln = GetNumNodStrCmp( );
  int elmnodlab[MAX_SCL_LAB];
  GetNodStrLabels(elmnodlab);

  // Add the scalar labels of the current element to the global vector.

  for (int i = 0; i < neln; i++)
  {
    int found = 0;
    for (int l = 0; l < NumNodSclLab; l++) // Search the element node label
    {
      if (elmnodlab[i] == VecNodSclLab[l])
      {
        found = 1;
        break;
      }
    }
    if (!found)  // New label => add to the global vector
    {
      VecNodSclLab[NumNodSclLab++] = elmnodlab[i];
    }
  }
}

// ============================== MassMat ===============================

void cElement :: MassMat(cMatrix &M)
{
  // Default case: element without mass matrix.

  M.Zero( );
}

// ======================================================= End of file =====

