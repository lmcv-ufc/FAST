// -------------------------------------------------------------------------
// shpiga.cpp - implementation of the isogeometric shape class.
// -------------------------------------------------------------------------
// Created:      25-May-2015     Elias Saraiva Barroso
//
// Modified:     
//               
// -------------------------------------------------------------------------

using namespace std;

#include "shpbsp.h"
#include "shpiga.h"
#include "gblvar.h"
#include "element.h"
#include "bspbasis.h"
#include "bernbasis.h"
#include "ctrlpnt.h"
#include "vec.h"
#include "subpat.h"

// -------------------------------------------------------------------------
// Class cShapeBsp:
// -------------------------------------------------------------------------

// =========================== cShapeBsp ===================================

cShapeBsp :: cShapeBsp( )
{
  SpanInd[0] = 0;
  SpanInd[1] = 0;
  SpanInd[2] = 0;
}

// =========================== ~cShapeBsp ==================================

cShapeBsp :: ~cShapeBsp(void)
{
  if (PatElmMap[pat])
  {  
    delete [] PatElmMap[pat];
    PatElmMap[pat] = 0;
  }
}

// ============================== PatElmMapAdd =============================

void cShapeBsp :: PatElmMapAdd(int elmid)
{
  // Insert the owner of this shape in patch/element map.
  int id = pat->getSpanID(SpanInd[0],SpanInd[1],SpanInd[2]);
  GetPatElmMap( ).find(pat)->second[id] = cElement :: GetElm(elmid);
}

// ================================= Read ==================================

void cShapeBsp :: ReadIGA(void)
{
  // Read the associated patch to this element.
  int idpat;

  if (!(in >> idpat))
  {
    cout << "Error in the read of patch label!\n";
    exit(0);
  }

  cPatch *patch;
  patch = cPatch :: FindPatch(idpat);

  if (!patch)
  {
    cout << "Invalid patch label: " << idpat << "\n";
    exit(0);
  }

  pat = dynamic_cast<cBspPat*> (patch);

  if (!pat)
  {
    cout << "Error in the read of isogeometric B-Spline element!\n";  
    cout << "Input patch must be B-Spline!\n";
    exit(0);
  }

  // Read the element knot span and load bézier 
  // extraction matrices in each direction.
  for(int i = 0; i < pat->getNumParVar( ); i++)
  {
    if (!(in >> SpanInd[i]) && !(SpanInd[i] > 0))
    {
      cout << "Error in the input of IGA element knot span!\n";
      exit(0);
    }
    SpanInd[i] -= 1;

    // Load Bézier extraction matrix.
    if (pat->getExtMat(i,0) == 0)
      pat->LoadExtMat( );
  }

  // Alloc memory for this patch in patch/shape map.
  if (PatElmMap.count(pat) == 0)
  {
    int nelm = pat->getNumSpan( );
    cElement** Evec = new cElement* [nelm];

    for(int e = 0; e < nelm; ++e)
      Evec[e] = 0; 

    PatElmMap[pat] = Evec;
  }

  // Initialize variables related to the patch subclass.
  Init( );
}

void cShapeBsp :: CompBasisIndex(const int &pvar, int &nvar, int* &basind)
{
  // Get the number of shape basis function.
  nvar = pat->getBasis(pvar).getDegree( ) + 1;

  // Get local basis function index.
  basind = new int [nvar];
  pat->getBasInd(pvar,SpanInd[pvar], basind);
}


// ============================== Init =====================================
/*
void cShapeIGA :: Init(void)
{
  // Auxiliary data.
  int    nid;         // Node label.
  cNode *node;        // Node pointer.
  int    dim;

  dim     = pat->getNumParVar( );
  NumNode = 1;
  for(int pvar = 0; pvar < dim; ++pvar)
  {
    // Get the number of shape basis function.
    NumBas[pvar] = pat->getBasis(pvar).getDegree( ) + 1;
    
    // Get local basis function index.
    BasInd[pvar] = new int [NumBas[pvar]];
    pat->GetBasInd(pvar,SpanInd[pvar], BasInd[pvar]);
    
    // Get control points of the element.
    NumNode *= NumBas[pvar];
  }

  ElmNode = new cNode*[NumNode];

  int i,j,k;
  for(int i = 0; i < NumBas[0]; i++)
  {
    // Evaluate control points tensor indices.
    i = j = k = 0;
    switch dim:
    {
      case 3: k = a / (NumBas[0] * NumBas[1]);
      case 2: j = (a - k * NumBas[0] * NumBas[1] ) / NumBas[0];
      case 1: i = a - k * NumBas[0] * NumBas[1] - j * NumBas[0];
    }

    // Get control point vector id.
    if (dim == 1)
      cpid = pat->GetCtrlPnt(BasInd[0][i]);
    else if (dim == 2)
      cpid = pat2d->GetCPInd(BasInd[0][j],BasInd[1][i]);
     
    cpid 

    // Get node label of current control point.
    if (dim == 1)
      nid = pat->GetCtrlPnt(BasInd[0][i])->getLabel( );
    else if (dim == 2)
      cpid = pat2d->GetCPInd(BasInd[0][j],BasInd[1][i]);
      nid = pat->GetCtrlPnt(BasInd[0][i],BasInd[][j)->getLabel( );


    // Get node label and correspondent node pointer.
    cout << "Bas " <<  BasInd[0][i] << endl;
    cout << "cpref " << pat->GetCtrlPnt(BasInd[0][i]) << endl;
    nid = pat->GetCtrlPnt(BasInd[0][i])->getLabel( );
    node = cNode :: FindNode(nid);

    if (!node)
    {
      cout << "Error in the read of control points!" << endl;
      cout << "Invalid node label." << endl;
      exit(0);
    }

    ElmNode[i] = node;

    
    if (i 
  }
}
*/

// ================================= SetTopology ===========================

void cShapeBsp :: SetTopology(cPatch *Pat, int span[])
{
  // Store the patch pointer.

  cBspPat* Patch = dynamic_cast<cBspPat*> (Pat);

  if (!Patch)
  {
    cout << "Invalid dynamic_cast in cShapeIGA :: SetTopology" << endl;
    exit(0);
  }

  pat = Patch;

  // Set the element span index according with patch dimension. 
  
  for(int i = 0; i < pat->getNumParVar( ); i++)
    SpanInd[i] = span[i];

  // Load patch extraction matrix.
  pat->LoadExtMat( );

  // Alloc memory for this patch in patch/shape map.
  if (PatElmMap.count(pat) == 0)
  {
    int nelm = 1;
    for(int i = 0; i < pat->getNumParVar( ); i++)
      nelm *= pat->getBasis(i).getKnotVec( ).getNumSpan( );
    PatElmMap[pat] = new cElement* [nelm];
  }

  // Initialize variables related to the patch subclass.

  Init( );
}

// ================================ GetEdge ================================
/*
int cShapeIGA :: GetEdge(int *corner, eShpType *type, int *nnode,
                        cNode **conn)
{
  cout << "IGA elements cannot be loaded by FEM loads! (Use IGA Loads)\n";
  return 0;
}*/


// ============================= EvalBspBas ================================

void cShapeBsp :: EvalBspBas(int bind, double t, cVector &Nr)
{
  // Map parameter coordinate from [-1,1] to [0,1].
  t = (t + 1.0)/2.0;

  // Compute the univariate Bernstein basis functions in the parent domain.
  int deg      = pat->getBasis(bind).getDegree( );
  double *Bern = new double [deg+1];
  nBernsteinBasis :: CompPols(deg,t,Bern);

  // Compute the univariate B-Spline functions w.r.t. the parent domain.
  double **C = pat->getExtMat(bind,SpanInd[bind]);
  
  for(int i = 0; i < deg+1; i++)
    for(int j = 0; j < deg+1; j++)
      Nr[i] += C[i][j] * Bern[j];

  // Release memory
  delete [] Bern;
}

// =========================== EvalBspBasDerv ==============================

void cShapeBsp :: EvalBspBasDerv(int bind, double t, cVector &dNr)
{
  // Map parameter coordinate from [-1,1] to [0,1].
  t = (t + 1.0)/2.0;

  // Compute the univariate Bernstein derivatives in the parent domain.
  int deg          = pat->getBasis(bind).getDegree( );
  double *BernDerv = new double [deg+1];
  nBernsteinBasis :: CompDervPols(deg,t,BernDerv);

  // Apply jacobian transformation to derivatives.
  for(int i = 0; i < deg+1; ++i)
    BernDerv[i] /= 2.0;

  // Compute the univariate B-Spline derivatives w.r.t. the parent domain.
  double **C = pat->getExtMat(bind,SpanInd[bind]);

  for(int i = 0; i < deg+1; i++)
    for(int j = 0; j < deg+1; j++)
      dNr[i] += C[i][j] * BernDerv[j];

  // Release memory
  delete [] BernDerv;
}

// =========================== GetPatSpanID ================================
/* TODO: REMOVE ME
int cShapeBsp :: GetPatSpanID( )
{
  int id = 0;
  switch (pat->getNumParVar( ))
  {
    case 3:
      id += SpanInd[2] * pat->getBasis(0).GetKnotVec( ).GetNumSpan( )
                       * pat->getBasis(1).GetKnotVec( ).GetNumSpan( );

    case 2:
      id += SpanInd[1] * pat->getBasis(0).GetKnotVec( ).GetNumSpan( );

    case 1:
      id += SpanInd[0];
  }

  return id;
}
*/
// =========================== MapStrPoint ================================

/*
bool cShapeIGA :: MapStrPoint(sPatPnt patpnt, sNatCoord &p)
{
  // Check patch label

  if (pat != cPatch :: GetPatch(patpnt.patid - 1))
    return false;

  // Get Span limits and check if the given point not belongs there

  double par[3][2];
  for(int i = 0; i < pat->getNumParVar( ); i++)
  {
    pat->GetSpanVal(i,SpanInd[i],par[i][0],par[i][1]); 

    if (!(patpnt.par[i] >= par[i][0] && patpnt.par[i] <= par[i][1]))
      return false;
  }

  // Map the patch point to element reference 

  double IntPntCoord[3] = {0.0,0.0,0.0};
  for(int i = 0; i < pat->getNumParVar( ); i++)
    IntPntCoord[i] = 2*(patpnt.par[i]-par[i][0])/(par[i][1]-par[i][0]) - 1;

  p.r = IntPntCoord[0];
  p.s = IntPntCoord[1];
  p.t = IntPntCoord[2];

  return true;
}
*/

// -------------------------------------------------------------------------
// Class cShapePlBspCurve:
// -------------------------------------------------------------------------

// =========================== cShapePlBspCurve ============================

cShapePlBspCurve :: cShapePlBspCurve(void)
{
  Type    = SHAPE_BSP_PLCURV;
  pat     = 0;
}

// ========================== ~cShapePlBspCurve ============================

cShapePlBspCurve :: ~cShapePlBspCurve(void)
{
}
 
// ============================== Init =====================================

void cShapePlBspCurve :: Init(void)
{
  // Auxiliary data.
  int    nid;         // Node label.
  int   *basind;      // Basis index 
  cNode *node;        // Node pointer.

  // Compute the shape function basis index w.r.t parent patch domain.
  CompBasisIndex(0,NumBas[0],basind);

  // Get control points of the element.
  NumNode = NumBas[0];
  ElmNode = new cNode*[NumNode];

  for(int i = 0; i < NumBas[0]; i++)
  {
    // Get node label and correspondent node pointer.
    nid = pat->getCtrlPnt(basind[i])->getLabel( );
    node = cNode :: FindNode(nid);

    ElmNode[i] = node;
  }

  // Release memory.
  delete [] basind;
}

// ============================== NodeNatCoord =============================

void cShapePlBspCurve :: NodeNatCoord(sNatCoord *c)
{
  c[0].r = -1.0;
  c[1].r =  0.0;
  c[2].r =  1.0;
}

// ================================= ShpFunc ===============================

void cShapePlBspCurve :: ShpFunc(sNatCoord p, double *N)
{
  double r = p.r;

  cVector Nr(NumNode);
  Nr.Zero( );

  // Compute the univariate B-Spline basis functions w.r.t. the parent domain.

  EvalBspBas(0,r,Nr);

  // Compute B-spline Shape function the the case of BSPLINE_1D patch.

  int a = 0;
  if (!pat->isRational( ))
  {
    for(int i = 0; i < NumBas[0]; i++)
      N[a++] = Nr[i];

    return;
  }

  // Compute the numerators and denominator for the NURBS functions.

  a = 0;
  double w = 0;

  for(int i = 0; i < NumBas[0]; i++)
  {
    N[a] = Nr[i] * ElmNode[i]->GetW( );
    w += N[a];
    a++;
  }

  // Divide by the denominators to complete the computation of the
  // functions w.r.t. the parent domain.

  for(int i = 0; i < a; i++)
    N[i] /= w; 
}

// =============================== DrvShpRST ===============================

void cShapePlBspCurve :: DrvShpRST(sNatCoord p, sNatCoord *dN)
{
  double r = p.r;

  cVector N(NumNode);
  cVector dNr(NumNode);
  cVector Nr(NumNode);

  N.Zero( );
  Nr.Zero( );
  dNr.Zero( );

  // Compute the univariate B-Spline basis functions and derivatives w.r.t. 
  // the parent domain.

  EvalBspBas(0,r,Nr);
  EvalBspBasDerv(0,r,dNr);

  // Compute B-spline Shape function in the case of BSPLINE_1D patch.

  if (!pat->isRational( ))
  {
    for(int i = 0; i < NumNode; i++)
    {
      N[i] = Nr[i];
      dN[i].r = dNr[i];
    }

    return;
  }

  // Compute the numerators and denominator for the NURBS functions.

  double w = 0, dwr = 0;
  for(int i = 0; i < NumNode; i++)
  {
    N[i] = Nr[i] * ElmNode[i]->GetW( );
    w += N[i];

    dN[i].r = dNr[i] * ElmNode[i]->GetW( );
    dwr += dN[i].r;
  }

  // Divide by the denominators to complete the computation of the
  // functions w.r.t. the parent domain.

  for(int i = 0; i < NumNode; i++)
  {
    N[i] /= w; 
    dN[i].r = (dN[i].r - N[i]*dwr) / w;
  }
}

// -------------------------------------------------------------------------
// Class cShapeBspCurve:
// -------------------------------------------------------------------------

// ============================= cShapeBspCurve ============================

cShapeBspCurve :: cShapeBspCurve(void)
{
  Type = SHAPE_BSP_CURV;
}

// ============================ ~cShapeBspCurve =============================

cShapeBspCurve :: ~cShapeBspCurve(void)
{
}

// -------------------------------------------------------------------------
// Class cShapePlBspSurf:
// -------------------------------------------------------------------------

// =========================== cShapePlBspSurf =============================

cShapePlBspSurf :: cShapePlBspSurf(void)
{
  TopType = QUADRILATERAL_TOPOLOGY;
  Type    = SHAPE_BSP_PLSURF;
  pat     = 0;
  ElmNode = 0;
}

// ========================== ~cShapePlBspSurf =============================

cShapePlBspSurf :: ~cShapePlBspSurf(void)
{
}
 
// ============================== Init =====================================

void cShapePlBspSurf :: Init(void)
{
  // Auxiliary data.
  int    nid;         // Node label.
  int   *basind[2];   // Basis index. 
  cNode *node;        // Node pointer.
  int    cpid;        // Control point label.

  // Compute the shape function basis index w.r.t parent patch domain.
  CompBasisIndex(0,NumBas[0],basind[0]);
  CompBasisIndex(1,NumBas[1],basind[1]);

  //cout << "Element span " << SpanInd[0] << "," << SpanInd[1] << endl;

  // Get control points of the element.

  NumNode = NumBas[0] * NumBas[1];
  ElmNode = new cNode*[NumNode];


//  static int count = 0;
//  cout << "elem " << ++count << endl;

  int    counter = 0;
  for(int j = 0; j < NumBas[1]; j++)
    for(int i = 0; i < NumBas[0]; i++, counter++)
    {
      cpid = pat->getCPInd(basind[0][i],basind[1][j]);
      nid  = pat->getCtrlPnt(cpid)->getLabel( );
      node = cNode :: FindNode(nid);

      ElmNode[counter] = node;
    }
}

// ============================== NodeNatCoord =============================

void cShapePlBspSurf :: NodeNatCoord(sNatCoord *c)
{
  double r,s;
  int    counter = 0;

  for(int j = 0; j < NumBas[1]; j++)
    for(int i = 0; i < NumBas[0]; i++)
    {
      r = -1.0+(2.0/(NumBas[0]-1))*i;
      s = -1.0+(2.0/(NumBas[1]-1))*j;

      c[counter].r = r;
      c[counter].s = s;
      counter++;
    } 
}

// ================================= ShpFunc ===============================

void cShapePlBspSurf :: ShpFunc(sNatCoord p, double *N)
{
  double r = p.r;
  double s = p.s;

  cVector Nr(NumBas[0]), Ns(NumBas[1]);

  Nr.Zero( );
  Ns.Zero( );

  // Compute the univariate B-Spline basis functions w.r.t. the parent domain.

  EvalBspBas(0,r,Nr);
  EvalBspBas(1,s,Ns);

  // Compute B-spline Shape function in the case of BSPLINE_2D patch.

  int a = 0;
  if (!pat->isRational( ))
  {
    for(int i = 0; i < NumBas[1]; i++)
      for(int j = 0; j < NumBas[0]; j++)
        N[a++] = Nr[j] * Ns[i];

     return;
  }

  // Compute the numerators and denominator for the tensor product
  // NURBS functions and derivatives.

  a = 0;
  double w = 0;

  for(int i = 0; i < NumBas[1]; i++)
    for(int j = 0; j < NumBas[0]; j++)
    {
      N[a] = Nr[j] * Ns[i] * ElmNode[a]->GetW( );
      w += N[a];

      a++;
    }

  // Divide by the denominators to complete the computation of the
  // functions and derivatives w.r.t. the parent domain.

  for(int i = 0; i < a; i++)
    N[i] /= w; 
}

// =============================== DrvShpRST ===============================

void cShapePlBspSurf :: DrvShpRST(sNatCoord p, sNatCoord *dN)
{ 
  double r = p.r;
  double s = p.s;

  cVector N(NumNode);
  cVector dNr(NumBas[0]), dNs(NumBas[1]);
  cVector Nr(NumBas[0]), Ns(NumBas[1]);

  N.Zero( );
  Nr.Zero( );
  Ns.Zero( );
  dNr.Zero( );
  dNs.Zero( );

  // Compute the univariate B-Spline basis functions and derivatives w.r.t. 
  // the parent domain.
  EvalBspBas(0,r,Nr);
  EvalBspBas(1,s,Ns);
  EvalBspBasDerv(0,r,dNr);
  EvalBspBasDerv(1,s,dNs);

  // Compute B-spline Shape function in the case of BSPLINE_2D patch.

  int a = 0;
  if (!pat->isRational( ))
  {
    for(int i = 0; i < NumBas[1]; i++)
      for(int j = 0; j < NumBas[0]; j++)
      {
        N[a] = Nr[j]*Ns[i];
    
        dN[a].r = dNr[j]*Ns[i];
        dN[a].s = dNs[i]*Nr[j];

	a++;
      }

     return;
  }

  // Compute the numerators and denominator for the tensor product
  // NURBS functions and derivatives.

  a = 0;
  double w = 0, dws = 0, dwr = 0;
  
  for(int i = 0; i < NumBas[1]; i++)
    for(int j = 0; j < NumBas[0]; j++)
    {

      N[a] = Nr[j] * Ns[i] * ElmNode[a]->GetW( );
      w += N[a];

      dN[a].r = dNr[j] *  Ns[i] * ElmNode[a]->GetW( );
      dN[a].s =  Nr[j] * dNs[i] * ElmNode[a]->GetW( );

      dwr += dN[a].r;
      dws += dN[a].s;
      a++;
    }

  // Divide by the denominators to complete the computation of the
  // functions and derivatives w.r.t. the parent domain.

  for(int i = 0; i < a; i++)
  {
    N[i] /= w; 
    dN[i].s = (dN[i].s - N[i]*dws) / w;
    dN[i].r = (dN[i].r - N[i]*dwr) / w;
  }
}

void cShapePlBspSurf :: DrvShpRST(sNatCoord p, sNatCoord *dN, sNatDrv *ddN)
{ 
  int i, j, a;
  sNatCoord pd;
  double tol = 1e-6;

  // Evaluate ,rr components.
  pd   = p;
  pd.r = p.r-tol/2.0; 
  DrvShpRST(pd,dN);

  for(i = 0, a = 0; i < NumBas[1]; i++)
    for(j = 0; j < NumBas[0]; ++j, ++a)
      ddN[a].rr = dN[a].r;

  pd   = p;
  pd.r = p.r+tol/2.0; 
  DrvShpRST(pd,dN);

  for(i = 0, a = 0; i < NumBas[1]; i++)
    for(j = 0; j < NumBas[0]; ++j, ++a)
      ddN[a].rr = (dN[a].r - ddN[a].rr)/tol;

  // Evaluate ,rs and ,ss components.
  pd   = p;
  pd.s = p.s-tol/2.0; 
  DrvShpRST(pd,dN);

  for(i = 0, a = 0; i < NumBas[1]; i++)
    for(j = 0; j < NumBas[0]; ++j, ++a)
    {
      ddN[a].rs = dN[a].r;
      ddN[a].ss = dN[a].s;
    }

  pd   = p;
  pd.s = p.s+tol/2.0; 
  DrvShpRST(pd,dN);

  for(i = 0, a = 0; i < NumBas[1]; i++)
    for(j = 0; j < NumBas[0]; ++j, ++a)
    {
      ddN[a].rs = (dN[a].r - ddN[a].rs)/tol;
      ddN[a].ss = (dN[a].s - ddN[a].ss)/tol;
    }
}

// ========================= GetPatEdge ====================================

cPatch* cShapePlBspSurf :: GetPatEdge(int f,int e, eShpType &t, int s[])
{ 
  // Evaluate a unique score for this sub-patch.
  int score  = f + e * 10 + pat->getLabel( ) * 100;

  // Create a new isopatch if it does not exist.
  if (IsoPatMap.count(score) == 0)
    IsoPatMap[score] = nGeomAlg :: nCurvAlg :: CreateSubPatch(pat,f,e);

  // Set shape type.
  t = SHAPE_BSP_PLCURV; 

  // inds[f][e] = activate parametric direction index.
  //
  //               Edges: {T,B,L,R}                //  Faces:
  static int inds[][4] = {{0,0,1,1}, {0,0,1,1},    //  Front and Back.
                          {2,2,1,1}, {2,2,1,1},    //  Left and Right.
                          {0,0,2,2}, {0,0,2,2}};   //  Top and Below.

  // Set span.
  s[0] = SpanInd[inds[f][e]];

  return IsoPatMap[score];
}

// ========================= GetPatFace ====================================

cPatch* cShapePlBspSurf :: GetPatFace(int f, eShpType &t, int s[])
{ 
  // Set shape type.
  t = Type;

  // Set iga shape knot span
  s[0] = SpanInd[0];
  s[1] = SpanInd[1];

  return pat;
}

// -------------------------------------------------------------------------
// Class cShapeBspSurf:
//

// ============================ cShapeBspSurf ==============================

cShapeBspSurf :: cShapeBspSurf(void)
{
  Type = SHAPE_BSP_SURF;
}

// =========================== ~cShapeBspSurf ==============================

cShapeBspSurf :: ~cShapeBspSurf(void)
{
}

// ================================ GetEdge ================================
/*
int cShapeSurfIGA :: GetEdge(int *corner, eShpType *type, int *nnode,
                            cNode **conn)
{
  static int NumEdge = 4;
  static tEdge Edge[] = { {SHAPE_L3_3D, 3, {0, 1, 2}},
                          {SHAPE_L3_3D, 3, {2, 3, 4}},
                          {SHAPE_L3_3D, 3, {4, 5, 6}},
                          {SHAPE_L3_3D, 3, {6, 7, 0}}};

  return(GetShapeEdge(NumEdge, Edge, corner, type, nnode, conn));
}
*/

// -------------------------------------------------------------------------
// Class cShapeBspSolid:
// -------------------------------------------------------------------------

// =========================== cShapeBspSolid =============================

cShapeBspSolid :: cShapeBspSolid(void)
{
  TopType = HEXAHEDRAL_TOPOLOGY;
  Type    = SHAPE_BSP_SOLID;
  pat     = 0;
  ElmNode = 0;
}

// ========================== ~cShapeBspSolid =============================

cShapeBspSolid :: ~cShapeBspSolid(void)
{
}

// ============================== Init =====================================

void cShapeBspSolid :: Init(void)
{
  // Auxiliary data.
  int    nid;         // Node label.
  int   *basind[3];   // Basis index 
  cNode *node;        // Node pointer.
  int    cpid;        // Control point label.

  // Compute the shape function basis index w.r.t parent patch domain.
  CompBasisIndex(0,NumBas[0],basind[0]);
  CompBasisIndex(1,NumBas[1],basind[1]);
  CompBasisIndex(2,NumBas[2],basind[2]);

  // Get control points of the element.

  NumNode = NumBas[0] * NumBas[1] * NumBas[2];
  ElmNode = new cNode*[NumNode];

  int    counter = 0;
  for(int k = 0; k < NumBas[2]; k++)
    for(int j = 0; j < NumBas[1]; j++)
      for(int i = 0; i < NumBas[0]; i++, counter++)
      {
        cpid = pat->getCPInd(basind[0][i],basind[1][j],basind[2][k]);
        nid  = pat->getCtrlPnt(cpid)->getLabel( );
        node = cNode :: FindNode(nid);
    
        ElmNode[counter] = node;
      }
}

// ============================== NodeNatCoord =============================

void cShapeBspSolid :: NodeNatCoord(sNatCoord *c)
{
  double r,s,t;
  int counter = 0;

  for(int i = 0; i < NumBas[2]; i++)
    for(int j = 0; j < NumBas[1]; j++)
      for(int k = 0; k < NumBas[0]; k++)
      {
        r = -1.0+(2.0/(NumBas[0]-1))*k;
        s = -1.0+(2.0/(NumBas[1]-1))*j;
        t = -1.0+(2.0/(NumBas[2]-1))*i;

        c[counter].r = r;
        c[counter].s = s;
        c[counter].t = t;
        counter++;
      } 
}

// ================================= ShpFunc ===============================

void cShapeBspSolid :: ShpFunc(sNatCoord p, double *N)
{
  double r = p.r;
  double s = p.s;
  double t = p.t;
  
  cVector Nr(NumBas[0]), Ns(NumBas[1]), Nt(NumBas[2]);

  Nr.Zero( );
  Ns.Zero( );
  Nt.Zero( );

  // Compute the univariate B-Spline basis functions w.r.t. the parent domain.

  EvalBspBas(0,r,Nr);
  EvalBspBas(1,s,Ns);
  EvalBspBas(2,t,Nt);
  
  // Compute B-spline Shape function in the case of BSPLINE_3D patch.

  int a = 0;

  if (!pat->isRational( ))
  {
    for(int i = 0; i < NumBas[2]; i++)
      for(int j = 0; j < NumBas[1]; j++)
        for(int k = 0; k < NumBas[0]; k++, a++)
          N[a] = Nr[k]*Ns[j]*Nt[i];

     return;
  }

  // Compute the numerators and denominator for the tensor product
  // NURBS functions and derivatives.

  a = 0;
  double w = 0;

  for(int i = 0; i < NumBas[2]; i++)
    for(int j = 0; j < NumBas[1]; j++)
      for(int k = 0; k < NumBas[0]; k++)
      {
        N[a] = Nr[k] * Ns[j] * Nt[i] * ElmNode[a]->GetW( );
        w += N[a];
      
        a++;
      }

  // Divide by the denominators to complete the computation of the
  // functions and derivatives w.r.t. the parent domain.

  for(int i = 0; i < a; i++)
    N[i] /= w;
}

// =============================== DrvShpRST ===============================

void cShapeBspSolid :: DrvShpRST(sNatCoord p, sNatCoord *dN)
{
  double r = p.r;
  double s = p.s;
  double t = p.t;

  cVector N(NumNode);
  cVector dNr(NumBas[0]), dNs(NumBas[1]), dNt(NumBas[2]);
  cVector Nr(NumBas[0]), Ns(NumBas[1]), Nt(NumBas[2]);

  N.Zero( );
  Nr.Zero( );
  Ns.Zero( );
  Nt.Zero( );
  dNr.Zero( );
  dNs.Zero( );
  dNt.Zero( );

  // Compute the univariate B-Spline basis functions and derivatives w.r.t. 
  // the parent domain.

  EvalBspBas(0,r,Nr);
  EvalBspBas(1,s,Ns);
  EvalBspBas(2,t,Nt);
  EvalBspBasDerv(0,r,dNr);
  EvalBspBasDerv(1,s,dNs);
  EvalBspBasDerv(2,t,dNt);

  // Compute B-spline Shape function in the case of BSPLINE_3D patch.

  int a = 0;

  if (!pat->isRational( ))
  {
    for(int i = 0; i < NumBas[2]; i++)
      for(int j = 0; j < NumBas[1]; j++)
        for(int k = 0; k < NumBas[0]; k++)
        {
          N[a] = Nr[k]*Ns[j]*Nt[i];
      
          dN[a].r = dNr[k] *  Ns[j] *  Nt[i];
          dN[a].s =  Nr[k] * dNs[j] *  Nt[i];
          dN[a].t =  Nr[k] *  Ns[j] * dNt[i];
        }

     return;
  }

  // Compute the numerators and denominator for the tensor product
  // NURBS functions and derivatives.

  a = 0;
  double w = 0, dws = 0, dwr = 0, dwt = 0;

  for(int i = 0; i < NumBas[2]; i++)
    for(int j = 0; j < NumBas[1]; j++)
      for(int k = 0; k < NumBas[0]; k++)
      {
        N[a] = Nr[k] * Ns[j] * Nt[i] * ElmNode[a]->GetW( );
        w += N[a];
     
        dN[a].r = dNr[k] *  Ns[j] *  Nt[i] * ElmNode[a]->GetW( );
        dN[a].s =  Nr[k] * dNs[j] *  Nt[i] * ElmNode[a]->GetW( );
        dN[a].t =  Nr[k] *  Ns[j] * dNt[i] * ElmNode[a]->GetW( );
        
	dwr += dN[a].r;
	dws += dN[a].s;
	dwt += dN[a].t;

        a++;
      }

  // Divide by the denominators to complete the computation of the
  // functions and derivatives w.r.t. the parent domain.

  for(int i = 0; i < a; i++)
  {
    N[i] /= w; 
    dN[i].r = (dN[i].r - N[i]*dwr) / w;
    dN[i].s = (dN[i].s - N[i]*dws) / w;
    dN[i].t = (dN[i].t - N[i]*dwt) / w;
  }
}

// ========================= GetPatFace ====================================


cPatch* cShapeBspSolid :: GetPatFace(int f, eShpType &t, int s[])
{
  int score  = f + pat->getLabel( ) * 100000;

  // Create a new isopatch if it does not exist.
  if (IsoPatMap.count(score) == 0)
    IsoPatMap[score] = nGeomAlg :: nSurfAlg :: CreateSubPatch(pat,f);

  // Set shape type.
  t = SHAPE_BSP_SURF; 

  // inds[f][2] = activate parametric direction indexes.
  //
  //                                    //  Faces:
  static int inds[][2] = {{0,1},{0,1},  //  Front and Back.
                          {1,2},{1,2},  //  Left and Right.
                          {0,2},{0,2}}; //  Top and Below.
  // Set span.
  s[0] = SpanInd[inds[f][0]];
  s[1] = SpanInd[inds[f][1]];

  return IsoPatMap[score];
}

// ======================================================= End of file =====
