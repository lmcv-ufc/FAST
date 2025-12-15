// -------------------------------------------------------------------------
// heat.cpp - implementation of heat class.
// -------------------------------------------------------------------------
// Created:      11-Dec-2019     Elias Saraiva Barroso
//
// -------------------------------------------------------------------------

#include <iostream>
#include <sstream>
#include <math.h>
#include <stdlib.h>

using namespace std;

#include "heat.h"
#include "timefunc.h"
#include "node.h"
#include "element.h"
#include "shape.h"
#include "intpoint.h"
#include "vec.h"
#include "mat.h"
#include "gblvar.h"
#include "gbldef.h"
#include "utl.h"
#include "field.h"

// -------------------------------------------------------------------------
// Static variables:
//
cHeat* cHeat :: Head = 0;
cHeat* cHeat :: Tail = 0;

// -------------------------------------------------------------------------
// Public methods:
//

// ============================= ReadNodalSources ==========================

void cHeat :: ReadNodalSources(void)
{
  int nsource;

  // Read the number of nodal sources.

  if (!(in >> nsource) || (nsource < 1))
  {
    cout << "Error in the input of the number of nodal sources!\n";
    exit(0);
  }

  // Read each source.

  int id;
  double val;
  for (int i = 0; i < nsource; i++)
  {
    if (!(in >> id) || !(in >> val))
    {
      cout << "Error in the input of the nodal source " << i+1 << "!\n";
      exit(0);
    }

    cNode *node = cNode :: FindNode(id);
    if (!node)
    {
      cout << "Error in the input of the nodal source " << i+1 << "(invalid node)!\n";
      exit(0);
    }

    new cNodalSource(node, val); // Create the nodal source
  }
}

// ============================= ReadBodySources ===========================

void cHeat :: ReadBodySources(void)
{
  int nsource;

  // Read the number of body sources.

  if (!(in >> nsource) || (nsource < 1))
  {
    cout << "Error in the input of the number of body forces!\n";
    exit(0);
  }

  // Read each source.

  for (int i = 0; i < nsource; i++)
  {
    cBodySource *source = new cBodySource( );
    source->Read( );
  }
}

// ============================ ReadBodyGenSources ==========================

void cHeat :: ReadBodyGenSources(void)
{
  int nsource;

  // Read the number of body sources.

  if (!(in >> nsource) || (nsource < 1))
  {
    cout << "Error in the input of the number of body forces!\n";
    exit(0);
  }

  // Read each source.

  for (int i = 0; i < nsource; i++)
  {
    cBodyGenSource *source = new cBodyGenSource( );
    source->Read( );
  }
}

// ============================= ReadFieldBodySources ======================

void cHeat :: ReadFieldBodySources(void)
{
  int nsource;

  // Read the number of body sources.

  if (!(in >> nsource) || (nsource < 1))
  {
    cout << "Error in the input of the number of body forces!\n";
    exit(0);
  }

  // Read each source.

  for (int i = 0; i < nsource; i++)
  {
    cBodySrcField *source = new cBodySrcField( );
    source->Read( );
  }
}

// ============================ ReadLineUnifFluxes =========================

void cHeat :: ReadLineUnifFluxes(void)
{
  int nflux;

  // Read the number of line fluxes.

  if (!(in >> nflux) || (nflux < 1))
  {
    cout << "Error in the input of the number of uniform line fluxes!\n";
    exit(0);
  }

  // Read each load.

  for (int i = 0; i < nflux; i++)
  {
    cEdgeUnifFlux *flux = new cEdgeUnifFlux( );
    flux->Read( );
  }
}

// ============================ ReadLineFieldFluxes =========================

void cHeat :: ReadLineFieldFluxes(void)
{
  int nflux;

  // Read the number of line fluxes.

  if (!(in >> nflux) || (nflux < 1))
  {
    cout << "Error in the input of the number of uniform line fluxes!\n";
    exit(0);
  }

  // Read each load.

  for (int i = 0; i < nflux; i++)
  {
    cEdgeFieldFlux *flux = new cEdgeFieldFlux( );
    flux->Read( );
  }
}

// ============================ ReadLineGenFluxes ==========================

void cHeat :: ReadLineGenFluxes(void)
{
  int nflux;

  // Read the number of line fluxes.

  if (!(in >> nflux) || (nflux < 1))
  {
    cout << "Error in the input of the number of uniform line fluxes!\n";
    exit(0);
  }

  // Read each load.

  for (int i = 0; i < nflux; i++)
  {
    cEdgeGenFlux *flux = new cEdgeGenFlux( );
    flux->Read( );
  }
}

// ============================ ReadFaceUnifLoads ==========================
/*
void cLoad :: ReadFaceUnifLoads(void)
{
  int nload;

  // Read the number of face loads.

  if (!(in >> nload) || (nload < 1))
  {
    cout << "Error in the input of the number of uniform face loads!\n";
    exit(0);
  }

  // Read each load.

  for (int i = 0; i < nload; i++)
  {
    cFaceUnifLoad *load = new cFaceUnifLoad( );
    load->Read( );
  }
}
*/

// ============================== EqvHeat ==================================

void cHeat :: EqvHeat(cElement *elm, double t, cVector &f)
{
  int n = f.Dim( );
  cVector fe(n);

  // Compute the total equivalent heat applied in the given element.

  f.Zero( );
  for (cHeat *heat = Head; heat != 0; heat = heat->Next)
  {
    if (heat->Type != NODAL_SOURCE && ((cElemHeat*)heat)->Elem == elm)
    {
      heat->ExtHeat(t, fe);
      f += fe;
    }
  }
}

// ================================ Destroy ================================

void cHeat :: Destroy(void)
{
  cHeat *heat = Head;

  while (heat)
  {
    cHeat *h = heat->Next;  // Save the next element
    delete heat;            // Delete the current element
    heat = h;
  }
}

// ================================ cHeat ==================================

cHeat :: cHeat(void)
{
  this->Next = 0;
  if (Head == 0)    // First element
    Head = this;
  else                  // Other elements
    Tail->Next = this;
  Tail = this;

  // Assign the current time function or use the default f(t) = 1.

  Func = cTimeFunc :: GetCurr( );
  if (!Func) Func = new cFuncConst(1.0);
}

// =============================== ~cHeat ==================================

cHeat :: ~cHeat(void)
{
}

// -------------------------------------------------------------------------
// Class cNodalSource:
// -------------------------------------------------------------------------


// ============================== cNodalSource =============================

cNodalSource :: cNodalSource(cNode *node, double val)
{
  Type = NODAL_SOURCE;
  Node = node;

  Value = val;
}

// ============================= ~cNodalSource =============================

cNodalSource:: ~cNodalSource(void)
{
}

// =============================== ExtHeat =================================

void cNodalSource :: ExtHeat(double t, cVector &f)
{
  double tf = Func->GetVal(t);
  f[0] = tf*Value;
}

// ============================== AddGlobVec ===============================

void cNodalSource :: AddGlobVec(cVector &felm, cVector &fglob)
{
  int dof = Node->GetDof(7);
  if (dof) fglob[dof-1] += felm[0];
}

// -------------------------------------------------------------------------
// Class cElemHeat:
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Public methods:
//

// =============================== cElemHeat ===============================

cElemHeat :: cElemHeat(void)
{
}

// ============================== ~cElemHeat ===============================

cElemHeat :: ~cElemHeat(void)
{
}

// =============================== ExtHeat =================================

void cElemHeat :: ExtHeat(double t, cVector &f)
{
  int i,j,k;
  double heatval;

  // Create the integration points.

  int nipt;
  int type = Elem->GetIntType( );
  int ord[] = {Elem->GetIntOrd( ), Elem->GetIntOrd( ), Elem->GetIntOrd( )};
  cIntPoint *intpnt = cIntPoint::CreateIntPoints(Shape->GetTopologyType(), type, ord, &nipt);

  // Get the nodal coordinates.

  int nnode = Shape->GetNumNode( );
  sNodeCoord *coord = new sNodeCoord[nnode];
  Shape->NodalCoord(coord);

  // Create n vector.
  cVector n(nnode);

  // Get memory for shape functions and derivatives.

  double *shpfnc = new double[nnode];
  double *mapfnc = new double[nnode];
  sNodeCoord *shpdrv = new sNodeCoord[nnode];

  // Loop over integration points

  f.Zero( );
  for (i = 0; i < nipt; i++)
  {
    // Get integration point data.

    sNatCoord p = intpnt[i].GetCoord( );
    double wgt  = intpnt[i].GetWeight( );

    // Get the external heat value.

    heatval = GetVal(p, coord);

    // Compute {n} vector.

    Shape->ShpFunc(p, shpfnc);
    Shape->MapFunc(p, mapfnc);
    n.Zero( );
    for (k = 0; k < nnode; k++) n[k] = shpfnc[k];

    // Compute {f} += coeff * heatval * {n}.

    double detJ;
    Shape->DrvShpXYZ(p, coord, &detJ, shpdrv);
    double volc  = Elem->VolCoeff(1.0, nnode, mapfnc, coord);
    double coeff = wgt*detJ*volc*heatval;
    for (k = 0; k < nnode; k++) f[k] += coeff * n[k];
  }

  // Multiply by the time function value.

  double tf = Func->GetVal(t);
  f *= tf;

  // Release memory.

  delete []coord;
  delete []intpnt;
  delete []shpfnc;
  delete []mapfnc;
  delete []shpdrv;
}

// ============================== AddGlobVec ===============================

void cElemHeat :: AddGlobVec(cVector &felm, cVector &fglob)
{
  // Get the associated dofs.

  int ndof = GetNumDofs( );
  int *dof = new int[ndof];
  GetDofs(dof);

  // Add local vector to the correct positions of the global one.

  for (int i = 0; i < ndof; i++)
    if (dof[i]) fglob[dof[i]-1] += felm[i];

  delete []dof;
}

// -------------------------------------------------------------------------
// Protected methods:
//

// =============================== GetNumDofs ==============================

int cElemHeat :: GetNumDofs(void)
{
  return Shape->GetNumNode( )*Elem->GetNumDofNode( );
}

// ================================ GetDofs ================================

void cElemHeat :: GetDofs(int *dof)
{
  // Get element active directions.

  int dir[8] = {0,0,0,0,0,0,0,0};
  Elem->GetActDir(dir);

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

// -------------------------------------------------------------------------
// Class cBodySource:
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Public methods:
//

// =============================== cBodySource ==============================

cBodySource :: cBodySource(void)
{
  Type  = BODY_SOURCE;
  Elem  = 0;
  Shape = 0;
}

// ============================== ~cBodySource ==============================

cBodySource :: ~cBodySource(void)
{
}

// ================================= Read ==================================

void cBodySource :: Read(void)
{
  // Get element with the heat source.

  int elmid;
  if (in >> elmid)
    Elem = cElement :: FindElm(elmid);
  if (!Elem)
  {
    cout << "Error in the input of body sources - element " << elmid << " not found!\n";
    exit(0);
  }
  Shape = Elem->GetShape( );

  // Read heat source data.

  if (!(in >> Value))
  {
    cout << "Error in the input of body sources!\n";
    exit(0);
  }
}

// ================================ GetVal =================================

double cBodySource :: GetVal(sNatCoord p, sNodeCoord *coord)
{
  return Value;
}

// -------------------------------------------------------------------------
// Class cEdgeUnifFlux:
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Public methods:
//

// ============================= cEdgeUnifFlux =============================

cEdgeUnifFlux :: cEdgeUnifFlux(void)
{
  Type  = LINE_FLUX;
  Elem  = 0;
  Shape = 0;
}

// ============================ ~cEdgeUnifFlux =============================

cEdgeUnifFlux :: ~cEdgeUnifFlux(void)
{
}

// ================================= Read ==================================

void cEdgeUnifFlux :: Read(void)
{
  // Get loaded element.

  int elmid;
  if (in >> elmid)
    Elem = cElement :: FindElm(elmid);
  if (!Elem)
  {
    cout << "Error in the input of uniform line fluxes!\n";
    cout << "Invalid element identify: " << elmid << "!\n";
    exit(0);
  }
  cShape *elmshp = Elem->GetShape( );

  // Get loaded edge.

  int corners[2];
  int nnode;
  eShpType type;
  cNode **conn = new cNode*[elmshp->GetNumNode( )];
  int code = 0;
  if ((in >> corners[0]) && (in >> corners[1]))
    code = elmshp->GetEdge(corners, &type, &nnode, conn);
  if (!code)
  {
    cout << "Error in the input of uniform line fluxes! \n";
    cout << "Invalid element corners input! \n";
    exit(0);
  }
  Shape = cShape :: CreateShape(type);
  Shape->SetNodes(conn, nnode);

  // Read load data.

  if (!(in >> FluxVal))
  {
    cout << "Error in the input of edge uniform fluxes!\n";
    exit(0);
  }

  delete [] conn;
}

// ================================ GetVal =================================

double cEdgeUnifFlux :: GetVal(sNatCoord p, sNodeCoord *coord)
{
  return -FluxVal;
}

// -------------------------------------------------------------------------
// Class cBodySrcField:
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Public methods:
//

// =============================== cBodySrcField ===========================

cBodySrcField :: cBodySrcField(void)
{
  Type  = BODY_SOURCE;
  Elem  = 0;
  Shape = 0;
}

// ============================== ~cBodySrcField ===========================

cBodySrcField :: ~cBodySrcField(void)
{
}

// ================================= Read ==================================

void cBodySrcField :: Read(void)
{
  // Get element with the heat source.

  int elmid;
  if (in >> elmid)
    Elem = cElement :: FindElm(elmid);
  if (!Elem)
  {
    cout << "Error in the input of body sources - element " << elmid << " not found!\n";
    exit(0);
  }
  Shape = Elem->GetShape( );
}

// ================================ GetVal =================================

double cBodySrcField :: GetVal(sNatCoord p, sNodeCoord *coord)
{
  // Evaluate cartesian coordinates.
  int nnode = Shape->GetNumNode( );
  double *shpfnc = new double[nnode];
  Shape->ShpFunc(p, shpfnc);

  double x = 0.0, y = 0.0, z = 0.0;
  for(int n = 0; n < nnode; ++n)
  {
    x += shpfnc[n] * coord[n].x;   
    y += shpfnc[n] * coord[n].y;   
    z += shpfnc[n] * coord[n].z;   
  }

  return cField :: GetInstance( )->Source(x,y,z);
}

// -------------------------------------------------------------------------
// Class cEdgeGenFlux:
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Public methods:
//

// ============================= cEdgeGenFlux ==============================

cEdgeGenFlux :: cEdgeGenFlux(void)
{
  Type  = LINE_FLUX;
  Elem  = 0;
  Shape = 0;
}

// ============================ ~cEdgeGenFlux ==============================

cEdgeGenFlux :: ~cEdgeGenFlux(void)
{
}

// ================================= Read ==================================

void cEdgeGenFlux :: Read(void)
{
  // Get loaded element.

  int elmid;
  if (in >> elmid)
    Elem = cElement :: FindElm(elmid);
  if (!Elem)
  {
    cout << "Error in the input of uniform general fluxes !\n";
    cout << "Invalid element identify: " << elmid << "!\n";
    exit(0);
  }
  cShape *elmshp = Elem->GetShape( );

  // Get loaded edge.

  int corners[2];
  int nnode;
  eShpType type;
  cNode **conn = new cNode*[elmshp->GetNumNode( )];
  int code = 0;
  if ((in >> corners[0]) && (in >> corners[1]))
    code = elmshp->GetEdge(corners, &type, &nnode, conn);
  if (!code)
  {
    cout << "Error in the input of uniform general fluxes! \n";
    cout << "Invalid element corners input! \n";
    exit(0);
  }
  Shape = cShape :: CreateShape(type);
  Shape->SetNodes(conn, nnode);

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
  expr.AddVar("f",f);

  // Compile input expression.
  expr.register_symbol_table( );
  int c = expr.compile(Utl :: GetExp(exprid));
  if (!c)
  {
    cout << "Error in compilation of the expression " << exprid << "!" << endl;
    cout << Utl :: GetExp(exprid) << endl;
    exit(0);
  }

  delete [] conn;
}

// ================================ GetVal =================================

double cEdgeGenFlux :: GetVal(sNatCoord p, sNodeCoord *coord)
{
  // Evaluate cartesian coordinates.
  int nnode = Shape->GetNumNode( );
  double *shpfnc = new double[nnode];
  Shape->ShpFunc(p, shpfnc);

  x  = y  = z  = 0.0;
  f  = 0.0;
  for(int n = 0; n < nnode; ++n)
  {
    x += shpfnc[n] * coord[n].x;   
    y += shpfnc[n] * coord[n].y;   
    z += shpfnc[n] * coord[n].z;   
  }

  // Evaluate expression.
  expr.evaluate( );

  // Release memory.
  delete [] shpfnc;

  return f;
}

// -------------------------------------------------------------------------
// Class cBodyGenSource:
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Public methods:
//

// ============================= cBodyGenSource ============================

cBodyGenSource :: cBodyGenSource(void)
{
  Type  = BODY_SOURCE;
  Elem  = 0;
  Shape = 0;
}

// ============================ ~cBodyGenSource ============================

cBodyGenSource :: ~cBodyGenSource(void)
{
}

// ================================= Read ==================================

void cBodyGenSource :: Read(void)
{
  // Get element with the heat source.

  int elmid;
  if (in >> elmid)
    Elem = cElement :: FindElm(elmid);
  if (!Elem)
  {
    cout << "Error in the input of body sources - element " << elmid << " not found!\n";
    exit(0);
  }
  Shape = Elem->GetShape( );

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
  expr.AddVar("f",f);

  // Compile input expression.
  expr.register_symbol_table( );
  int c = expr.compile(Utl :: GetExp(exprid));
  if (!c)
  {
    cout << "Error in compilation of the expression " << exprid << "!" << endl;
    cout << Utl :: GetExp(exprid) << endl;
    exit(0);
  }
}

// ================================ GetVal =================================

double cBodyGenSource :: GetVal(sNatCoord p, sNodeCoord *coord)
{
  // Evaluate cartesian coordinates.
  int nnode = Shape->GetNumNode( );
  double *shpfnc = new double[nnode];
  Shape->ShpFunc(p, shpfnc);

  x  = y  = z  = 0.0;
  f  = 0.0;
  for(int n = 0; n < nnode; ++n)
  {
    x += shpfnc[n] * coord[n].x;   
    y += shpfnc[n] * coord[n].y;   
    z += shpfnc[n] * coord[n].z;   
  }

  // Evaluate expression.
  expr.evaluate( );

  // Release memory.
  delete [] shpfnc;

  return f;
}

// -------------------------------------------------------------------------
// Class cEdgeFieldFlux:
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Public methods:
//

// ============================= cEdgeFieldFlux ============================

cEdgeFieldFlux :: cEdgeFieldFlux(void)
{
  Type  = LINE_FLUX;
  Elem  = 0;
  Shape = 0;
}

// ============================ ~cEdgeUnifFlux =============================

cEdgeFieldFlux :: ~cEdgeFieldFlux(void)
{
}

// ================================= Read ==================================

void cEdgeFieldFlux :: Read(void)
{
  // Get loaded element.

  int elmid;
  if (in >> elmid)
    Elem = cElement :: FindElm(elmid);
  if (!Elem)
  {
    cout << "Error in the input of uniform line fluxes!\n";
    cout << "Invalid element identify: " << elmid << "!\n";
    exit(0);
  }
  cShape *elmshp = Elem->GetShape( );

  // Get loaded edge.

  int corners[2];
  int nnode;
  eShpType type;
  cNode **conn = new cNode*[elmshp->GetNumNode( )];
  int code = 0;
  if ((in >> corners[0]) && (in >> corners[1]))
    code = elmshp->GetEdge(corners, &type, &nnode, conn);
  if (!code)
  {
    cout << "Error in the input of uniform line fluxes! \n";
    cout << "Invalid element corners input! \n";
    exit(0);
  }
  Shape = cShape :: CreateShape(type);
  Shape->SetNodes(conn, nnode);

  delete [] conn;
}

// ================================ GetVal =================================

double cEdgeFieldFlux :: GetVal(sNatCoord p, sNodeCoord *coord)
{
  // Evaluate cartesian coordinates.
  int nnode = Shape->GetNumNode( );
  double *shpfnc = new double[nnode];
  Shape->ShpFunc(p, shpfnc);

  double x = 0.0, y = 0.0, z = 0.0;
  for(int n = 0; n < nnode; ++n)
  {
    x += shpfnc[n] * coord[n].x;   
    y += shpfnc[n] * coord[n].y;   
    z += shpfnc[n] * coord[n].z;   
  }

  // Get flux vector.
  double flux[3];
  cField :: GetInstance( )->Flux(x,y,z,flux);

  cout << "func" << flux[0] << " " << flux[1] << " " << flux[2] << endl;

  // Evaluate rotation matrix.
  cMatrix R(3,3);
  Shape->LocalSys(p,coord,R);

  // Compute normal vector.
  cVector nx(3), n(3);
  nx.Zero( );
  nx[1] = 1.0;
  n     = R * nx;
  n.Normalize( );
  n    *= -1.0;

  cout << "n " << n[0] << " " << n[1] << " " << n[2] << endl;

  // Release memory.
  delete [] shpfnc;
  
  // Compute Forces.
  cout << "x " << x <<" y " << y << " ";
  cout << "flux " << -(flux[0] * n[0] + flux[1] * n[1] + flux[2] * n[2]) << endl;
  return (flux[0] * n[0] + flux[1] * n[1] + flux[2] * n[2]);
}

// -------------------------------------------------------------------------
// Class cFaceUnifFlux:
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Public methods:
//

/*
// ============================= cFaceUnifLoad =============================

cFaceUnifLoad :: cFaceUnifLoad(void)
{
  Type  = FACE_LOAD;
  Elem  = 0;
  Shape = 0;
  Local = 1;
  LoadVal.Resize(3);
  LoadVal.Zero( );
}

// ============================ ~cFaceUnifLoad =============================

cFaceUnifLoad :: ~cFaceUnifLoad(void)
{
}

// ================================= Read ==================================

void cFaceUnifLoad :: Read(void)
{
  // Get loaded element.

  int elmid;
  if (in >> elmid)
    Elem = cElement :: FindElm(elmid);
  if (!Elem)
  {
    cout << "Error in the input of uniform area loads !\n";
    exit(0);
  }
  cShape *elmshp = Elem->GetShape( );

  // Get loaded face.

  int corners[3];
  int nnode;
  eShpType type;
  cNode **conn = new cNode*[elmshp->GetNumNode( )];
  int code = 0;
  if ((in >> corners[0]) && (in >> corners[1]) && (in >> corners[2]))
    code = elmshp->GetFace(corners, &type, &nnode, conn);
  if (!code)
  {
    cout << "Error in the input of uniform face (area) loads !\n";
    exit(0);
  }
  Shape = cShape :: CreateShape(type);
  Shape->SetNodes(conn);

  // Read load data.

  if (!(in >> Local) || !(in >> LoadVal[0]) || !(in >> LoadVal[1]) || !(in >> LoadVal[2]))
  {
    cout << "Error in the input of face uniform loads!\n";
    exit(0);
  }
}

// ================================ GetVal =================================

void cFaceUnifLoad :: GetVal(sNatCoord p, sNodeCoord *coord, double *val)
{
  if (Local)
  {
    cMatrix R(3, 3);
    Shape->LocalSys(p, coord, R);  // Get rotation matrix
    cVector qr(3);
    qr = R*LoadVal;                // Rotate to global system
    for (int i = 0; i < 3; i++) val[i] = qr[i];
  }
  else
  {
    for (int i = 0; i < 3; i++) val[i] = LoadVal[i];
  }
  for (int i = 3; i < 6; i++) val[i] = 0.0;
  
  return;
}
*/

// ======================================================= End of file =====
