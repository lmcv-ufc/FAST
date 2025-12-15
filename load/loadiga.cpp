// -------------------------------------------------------------------------
// loadiga.cpp - implementation of load iga class.
// -------------------------------------------------------------------------
// Created:      20-Feb-2015     Elias Saraiva Barroso
//
// Modified:
//
// -------------------------------------------------------------------------

#include <iostream>
#include <sstream>
#include <stdlib.h>

using namespace std;

#include "loadiga.h"
#include "load.h"
#include "shpiga.h"
#include "shpbsp.h"
#include "patch.h"
#include "element.h"
#include "gblvar.h"
#include "utl.h"

// -------------------------------------------------------------------------
// Class cPatConcLoad:
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Public methods:
//

// =============================== cPatConcLoad ===============================

cPatConcLoad :: cPatConcLoad(void)
{
}

// ============================== ~cPatConcLoad ===============================

cPatConcLoad :: ~cPatConcLoad(void)
{
}

// ================================= Read ==================================

void cPatConcLoad :: Read(void)
{
  //  Auxiliary data.
  int elmid;
  int patid;
  double SpanLim[3][2];

  // Get loaded patch.
  cPatch *pat = 0;
  if (in >> patid)
    pat = cPatch :: FindPatch(patid);
  if (!pat)
  {
    cout << "Error in the input of patch concentrade forces - patch " << patid << " not found!\n";
    exit(0);
  }

  // Get parametric coord
  int dim = pat->getNumParVar( );
  double par[3];
  for(int i = 0; i < dim; i++)
    if (!(in >> par[i]))
    {
      cout << "Error in the input of patch concentrate force!\n";
      exit(0);
    }

  // Check if the given parametric point within patch bounds
  double inf, sup;
  for(int i = 0; i < dim; i++)
  {
    // Get parametric variable limits in dimension i.
    pat->getParVarLimits(i,inf,sup);

    if (par[i] < inf || par[i] > sup)
    {
      cout << "Error in the input of patch concentrate force!\n";
      cout << "Parametric coordenates should be defined in the patch bounds!\n";
      exit(0);
    }
  }

  // Get element label.
  elmid = pat->getSpanID(par[0],par[1],par[2]);

  // Test if shape object exists.
  const PatElmStdMap &Map = cShapeIGA :: GetPatElmMap( );
  bool skip = false;

  if (Map.count(pat) == 0)
    skip = true;
  else if (Map.find(pat)->second[elmid] == 0)
    skip = true;

  if (skip)
  {
    cout << "Cannot find corresponding element in given patch." << endl;
    exit(0);
  }

  // Get element and shape.
  Elem  = Map.find(pat)->second[elmid];
  Shape = Elem->GetShape( ); 

  // Get shape span limits.
  pat->getSpanLimits(elmid,SpanLim);

  // Compute parametric coordinates in element reference [-1,1].
  cout << "elmid: " << Elem->GetLabel( ) << endl;
  for(int pvar = 0; pvar < pat->getNumParVar( ); ++pvar)
  {
    ParVal[pvar]  = par[pvar] - SpanLim[pvar][0];
    ParVal[pvar] *= 2.0/(SpanLim[pvar][1] - SpanLim[pvar][0]);
    ParVal[pvar] -= 1.0;
  }

  // Read load data.
  if (!(in >> LoadVal[0]) || !(in >> LoadVal[1]) || !(in >> LoadVal[2]))
  {
    cout << "Error in the input of parametric concentrate force!\n";
    exit(0);
  }
}

// -------------------------------------------------------------------------
// Class cEdgeUnifLoadIGA:
// -------------------------------------------------------------------------

// ============================ cEdgeUnifLoadIGA ===========================

cEdgeUnifLoadIGA :: cEdgeUnifLoadIGA(void) : cEdgeUnifLoad( )
{
}

// =========================== ~cEdgeUnifLoadIGA ===========================

cEdgeUnifLoadIGA :: ~cEdgeUnifLoadIGA(void)
{
}

// ================================= Read ==================================

void cEdgeUnifLoadIGA :: Read(void)
{
  // Get loaded element.

  int elmid;
  if (in >> elmid)
    Elem = cElement :: FindElm(elmid);
  if (!Elem)
  {
    cout << "Error in the input of uniform line loads !\n";
    exit(0);
  }

  //cShapeBsp *elmshp = dynamic_cast<cShapeBsp*>(Elem->GetShape( ));
  cShapeBsp *elmshp = dynamic_cast<cShapeBsp*>(Elem->GetShape( )->GetBaseShp( )); 

  if (!elmshp)
  {
    cout << "IGA loads can only be applied on isogeometric B-Spline elements!\n";
    exit(0);
  }

  // Read the face and edge where the force is  applied.
  eFaceType ftype;
  eEdgeType etype;

  in >> ftype;
  in >> etype;

  facetype = (int) ftype;
  edgetype = (int) etype;

  // Get Edge patch.
  cPatch *pat;
  eShpType type;
  int span[1];

  pat = elmshp->GetPatEdge(facetype,edgetype,type,span);

  Shape = cShape :: CreateShape(type);
  cShapeBsp *shpiga = dynamic_cast<cShapeBsp*>(Shape);

  if (!shpiga)
  {
    cout << "Error in the function GetEdgePat, it is returning non-iga shape type!\n";
    exit(0);
  }

  shpiga->SetTopology(pat,span);

  // Read load data.
  if (!(in >> Local) || !(in >> LoadVal[0]) || !(in >> LoadVal[1]) || !(in >> LoadVal[2]))
  {
    cout << "Error in the input of edge uniform loads!\n";
    exit(0);
  }
}

// -------------------------------------------------------------------------
// Class cEdgeUnifMomentIGA:
// -------------------------------------------------------------------------

// ============================ cEdgeUnitMomentIGA =========================

cEdgeUnifMomentIGA :: cEdgeUnifMomentIGA(void) : cEdgeUnifMoment( )
{
}

// =========================== ~cEdgeUnifMomentIGA =========================

cEdgeUnifMomentIGA :: ~cEdgeUnifMomentIGA(void)
{
}

// ================================= Read ==================================

void cEdgeUnifMomentIGA :: Read(void)
{
  // Get loaded element.

  int elmid;
  if (in >> elmid)
    Elem = cElement :: FindElm(elmid);
  if (!Elem)
  {
    cout << "Error in the input of uniform line loads !\n";
    exit(0);
  }

  //cShapeBsp *elmshp = dynamic_cast<cShapeBsp*>(Elem->GetShape( ));
  cShapeBsp *elmshp = dynamic_cast<cShapeBsp*>(Elem->GetShape( )->GetBaseShp( )); 

  if (!elmshp)
  {
    cout << "IGA loads can only be applied on isogeometric B-Spline elements!\n";
    exit(0);
  }

  // Read the face and edge where the force is  applied.
  eFaceType ftype;
  eEdgeType etype;

  in >> ftype;
  in >> etype;

  facetype = (int) ftype;
  edgetype = (int) etype;

  // Get Edge patch.
  cPatch *pat;
  eShpType type;
  int span[1];

  pat = elmshp->GetPatEdge(facetype,edgetype,type,span);

  Shape = cShape :: CreateShape(type);
  cShapeBsp *shpiga = dynamic_cast<cShapeBsp*>(Shape);

  if (!shpiga)
  {
    cout << "Error in the function GetEdgePat, it is returning non-iga shape type!\n";
    exit(0);
  }

  shpiga->SetTopology(pat,span);

  // Read load data.
  if (!(in >> Local) || !(in >> LoadVal[0]) || !(in >> LoadVal[1]) || !(in >> LoadVal[2]))
  {
    cout << "Error in the input of edge uniform loads!\n";
    exit(0);
  }
}

// -------------------------------------------------------------------------
// Class cFaceUnifLoadIGA:
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Public methods:
//

// ============================ cFaceUnifLoadIGA ===========================

cFaceUnifLoadIGA :: cFaceUnifLoadIGA(void)
{
}

// =========================== ~cFaceUnifLoadIGA ===========================

cFaceUnifLoadIGA :: ~cFaceUnifLoadIGA(void)
{
}

// ================================= Read ==================================

void cFaceUnifLoadIGA :: Read(void)
{
  // Get loaded element.
  int elmid;
  if (in >> elmid)
    Elem = cElement :: FindElm(elmid);
  if (!Elem)
  {
    cout << "Error in the input of uniform line loads !\n";
    exit(0);
  }

  // Check if the element is isogeometric

  cShapeBsp *elmshp = dynamic_cast<cShapeBsp*>(Elem->GetShape( ));
  if (!elmshp)
  {
    cout << "IGA loads can only be applied on IGA elements!\n";
    exit(0);
  }

  // Read the face where the force is applied

  eFaceType ftype;
  in >> ftype;
  facetype = (int) ftype;

  // Get Edge patch

  cPatch *pat;
  eShpType type;
  int span[2];
  pat = elmshp->GetPatFace(facetype,type,span);
  Shape = cShape :: CreateShape(type);
  cShapeBsp *shpiga = dynamic_cast<cShapeBsp*>(Shape);
  if (!shpiga)
  {
    cout << "Error in the function GetPatFace, it is returning non-iga shape type!\n";
    exit(0);
  }
  shpiga->SetTopology(pat,span);

  // Read load data.

  if (!(in >> Local) || !(in >> LoadVal[0]) || !(in >> LoadVal[1]) || !(in >> LoadVal[2]))
  {
    cout << "Error in the input of area uniform loads!\n";
    exit(0);
  }
}

// -------------------------------------------------------------------------
// Class cEdgeGenLoadIGA:
// -------------------------------------------------------------------------

// ============================ cEdgeGenLoadIGA ==========================

cEdgeGenLoadIGA :: cEdgeGenLoadIGA(void) : cEdgeGenLoad( )
{
}

// =========================== ~cEdgeUnifLoadIGA ===========================

cEdgeGenLoadIGA :: ~cEdgeGenLoadIGA(void)
{
}

// ================================= Read ==================================

void cEdgeGenLoadIGA :: Read(void)
{
  // Get loaded element.

  int elmid;
  if (in >> elmid)
    Elem = cElement :: FindElm(elmid);
  if (!Elem)
  {
    cout << "Error in the input of uniform line loads !\n";
    exit(0);
  }

  cShapeBsp *elmshp = dynamic_cast<cShapeBsp*>(Elem->GetShape( ));

  if (!elmshp)
  {
    cout << "IGA loads can only be applied on isogeometric B-Spline elements!\n";
    exit(0);
  }

  // Read the face and edge where the force is  applied.
  eFaceType ftype;
  eEdgeType etype;

  in >> ftype;
  in >> etype;

  facetype = (int) ftype;
  edgetype = (int) etype;

  // Get Edge patch.
  cPatch *pat;
  eShpType type;
  int span[1];

  pat = elmshp->GetPatEdge(facetype,edgetype,type,span);

  Shape = cShape :: CreateShape(type);
  cShapeBsp *shpiga = dynamic_cast<cShapeBsp*>(Shape);

  if (!shpiga)
  {
    cout << "Error in the function GetEdgePat, it is returning non-iga shape type!\n";
    exit(0);
  }

  shpiga->SetTopology(pat,span);

  // Get input expression label.
  int exprid;
  if (!(in >> exprid) || (exprid <= 0) || (exprid > Utl :: GetNumExp( )))
  {
    cout << "Invalid expression label: " << exprid << "!\n";
    exit(0);
  }

  // Setup variables table.
  expr.AddVar("x",x);
  expr.AddVar("y",y);
  expr.AddVar("z",z);
  expr.AddVar("fx",fx);
  expr.AddVar("fy",fy);
  expr.AddVar("fz",fz);

  // Compile input expression.
  expr.register_symbol_table( );
  int c = expr.compile(Utl :: GetExp(exprid));
  if (!c)
  {
    cout << "Error in compilation of the expression " << exprid << "!" << endl;
    exit(0);
  }
}

// ======================================================= End of file =====
