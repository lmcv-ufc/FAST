// -------------------------------------------------------------------------
// shpiga.cpp - implementation of the isogeometric shape class.
// -------------------------------------------------------------------------
// Created:      25-May-2017     Elias Saraiva Barroso
//
// Modified:     
//               
// -------------------------------------------------------------------------

#include <iostream>
#include <cstdlib>

using namespace std;

#include "shpbez.h"
#include "bernbasis.h"
#include "gblvar.h"

// -------------------------------------------------------------------------
// Class cShapePlBezCurve:
// -------------------------------------------------------------------------

// =========================== cShapePlBezCurve ============================

cShapePlBezCurve :: cShapePlBezCurve(void)
{
  Type     = SHAPE_BEZ_PLCURV;
  ElmNode  = 0;
  Degree   = 0;
  Rational = true;
}

// ========================== ~cShapePlBezCurve ============================

cShapePlBezCurve :: ~cShapePlBezCurve(void)
{
}

// ================================= Read ==================================

void cShapePlBezCurve :: Read(void)
{
  // Read Rational flag and degree.
  if (!(in >> Rational) || !(in >> Degree))
  {
    cout << "Error in the input of bézier curve shape!\n";
    exit(0);
  }

  // Setting FEM node input data.
  NumNode = Degree + 1;
  ElmNode = new cNode* [NumNode];
  
  // Read each control point.
  cShape :: Read( );
}

// =============================== SetNodes ================================

void  cShapePlBezCurve :: SetNodes(cNode **node, int size)
{
  if (size == 0)
    size = NumNode;

  Degree  = size - 1;
  NumNode = size;
  ElmNode = new cNode* [size];
  for (int i = 0; i < size; i++)
  {
    ElmNode[i] = node[i];
    if (!ElmNode[i])
    {
      cout << "Error: invalied node in SetNodes function!\n";
      exit(0);
    }
  }
}

// ============================== NodeNatCoord =============================

void cShapePlBezCurve :: NodeNatCoord(sNatCoord *c)
{
/*
  if (cControl :: GetIgaEqvMesh( ))
  {
    c[0].r =  0.0; c[0].s =  0.0;
    c[1].r =  0.5; c[1].s =  0.0;
    c[2].r =  1.0; c[2].s =  0.0;
    c[3].r =  0.5; c[3].s =  0.5;
    c[4].r =  0.0; c[4].s =  1.0;
    c[5].r =  0.0; c[5].s =  0.5;
  }
  else
  {

  // Precisa implementar isso! 

    double r,s;
    int counter = 0;

    for(int i = 0; i < NumBas[1]; i++)
     for(int j = 0; j < NumBas[0]; j++)
     {
       r = -1.0+(2.0/(NumBas[0]-1))*j;
       s = -1.0+(2.0/(NumBas[1]-1))*i;

       c[counter].r = r;
       c[counter].s = s;
       counter++;
     } 
  }
*/  
}

// ================================= ShpFunc ===============================

void cShapePlBezCurve :: ShpFunc(sNatCoord p, double *N)
{
  // Map parameter coordinate from [-1,1] to [0,1].
  double r = (p.r + 1)/2.0;

  // Compute the univariate Bernstein basis functions w.r.t.
  double *Basis = new double [NumNode];
  nBernsteinBasis :: CompPols(Degree,r,Basis);

  // Compute non-rational shape functions.
  if (!Rational)
  {
    for(int i = 0; i < NumNode; i++)
        N[i] = Basis[i];

     return;
  }

  // Compute weight function for rational basis.
  double w = 0;
  for(int i = 0; i < NumNode; i++)
    w += Basis[i] * ElmNode[i]->GetW( );

  // Divide by the weight function to complete the computation.
  for(int i = 0; i < NumNode; i++)
    N[i] = Basis[i] * ElmNode[i]->GetW( ) / w;

  // Release Memory.
  delete [] Basis;
}

// =============================== DrvShpRST ===============================

void cShapePlBezCurve :: DrvShpRST(sNatCoord p, sNatCoord *dN)
{ 
  // Map parameter coordinate from [-1,1] to [0,1].
  double r = (p.r + 1.0)/2.0;

  double *N   = new double [NumNode];
  double *dNr = new double [NumNode];

  // Compute the univariate Bernstein basis functions and derivatives w.r.t. 
  // the element domain.

  nBernsteinBasis :: CompPols(Degree,r,N);
  nBernsteinBasis :: CompDervPols(Degree,r,dNr);

  // Apply jacobian transformation to derivatives.
  for(int i = 0; i < NumNode; ++i)
    dNr[i] /= 2.0;

  // Return Bézier basis in the case of non-rational patch.

  if (!Rational)
  {
    for(int i = 0; i < NumNode; i++)
      dN[i].r = dNr[i];

    return;
  }

  // Compute the numerators and denominator for the rational
  // Bézier functions and derivatives.

  double w = 0, dwr = 0;
  
  for(int i = 0; i < NumNode; ++i)
  {
    N[i] *= ElmNode[i]->GetW( ); 
    w += N[i];

    dN[i].r = dNr[i] * ElmNode[i]->GetW( );
    dwr += dN[i].r;
  }

  // Divide by the denominators to complete the computation of the
  // functions and derivatives w.r.t. the parent domain.

  for(int i = 0; i < NumNode; i++)
  {
    N[i] /= w; 
    dN[i].r = (dN[i].r - N[i]*dwr) / w;
  }

  // Release Memory.
  delete [] N;
  delete [] dNr;
}

// ================================ GetEdge ================================

int cShapePlBezCurve :: GetEdge(int *corner, eShpType *type, int *nnode, cNode **conn)
{
  // TODO: Must be implemented!
  return 0;
}

// -------------------------------------------------------------------------
// Class cShapeBezCurve:
// -------------------------------------------------------------------------

// ============================= cShapeBezCurve ============================

cShapeBezCurve :: cShapeBezCurve(void)
{
  Type = SHAPE_BEZ_CURV;
}

// ============================ ~cShapeBezCurve =============================

cShapeBezCurve :: ~cShapeBezCurve(void)
{
}

// -------------------------------------------------------------------------
// Class cShapePlBezTrianSurf:
// -------------------------------------------------------------------------

// =========================== cShapePlSurfIGA =============================

cShapePlBezTrianSurf :: cShapePlBezTrianSurf(void)
{
  Type = SHAPE_BEZTRI_PLSURF;
  TopType  = TRIANGULAR_TOPOLOGY;
  ElmNode = 0;
}

// ========================== ~cShapePlSurfIGA =============================

cShapePlBezTrianSurf :: ~cShapePlBezTrianSurf(void)
{
}

// ================================= Read ==================================

void cShapePlBezTrianSurf :: Read(void)
{

  // Read Rational flag and degree.
  if (!(in >> Rational) || !(in >> Degree))
  {
    cout << "Error in the input of bézier curve shape!\n";
    exit(0);
  }

  // Setting FEM node input data.
  NumNode = (Degree + 1) * (Degree + 2) / 2;
  ElmNode = new cNode* [NumNode];
  
  // Read each control point.
  cShape :: Read( );
}

// ============================== NodeNatCoord =============================

void cShapePlBezTrianSurf :: NodeNatCoord(sNatCoord *c)
{
  int i,j;
  int p = Degree;
  
  for(int id = 0; id < NumNode; id++)
  {
    // Evaluate berycentric indices
  
    if (id == 0)
    {
      i = p; 
      j = 0;
    }
    else
    {
      if (j == 0)
      {
        --i;
        j = p - i;
      }
      else
        --j;
    }
  
    // Evaluate parametric values.
    
    c[id].r = ((double) i)/p;
    c[id].s = ((double) j)/p;
  }
}

// ================================= ShpFunc ===============================

void cShapePlBezTrianSurf :: ShpFunc(sNatCoord p, double *N)
{
  double r = p.r;
  double s = p.s;

  // Compute the bivariates Bernstein basis functions w.r.t. the parent domain.

  double *Basis = new double [NumNode];
  nBernsteinBasis :: CompTrianPols(Degree,r,s,Basis);

  // Compute non-rationals shape functions.

  if (!Rational)
  {
    for(int i = 0; i < NumNode; i++)
      N[i] = Basis[i];

     return;
  }

  // Compute the numerators and denominator for the bézier triangle 
  // functions and derivatives.

  double w = 0;
  for(int i = 0; i < NumNode; i++)
    w += Basis[i] * ElmNode[i]->GetW( );

  // Divide by the denominators to complete the computation of the
  // functions and derivatives w.r.t. the parent domain.

  for(int i = 0; i < NumNode; i++)
    N[i] = Basis[i] * ElmNode[i]->GetW( ) / w;

  // Release memory.
  delete [] Basis;
}

// =============================== DrvShpRST ===============================

void cShapePlBezTrianSurf :: DrvShpRST(sNatCoord p, sNatCoord *dN)
{ 
  double r = p.r;
  double s = p.s;

  double *N   = new double [NumNode];
  double *dNr = new double [NumNode];
  double *dNs = new double [NumNode];

  // Compute the univariate B-Spline basis functions and derivatives w.r.t. 
  // the parent domain.
  nBernsteinBasis :: CompTrianPols(Degree,r,s,N);
  nBernsteinBasis :: CompTrianDervPols(Degree,r,s,dNr,dNs);

  // Compute B-spline Shape function in the case of BSPLINE_2D patch.
  if (!Rational)
  {
    for(int i = 0; i < NumNode; i++)
    {
      dN[i].r = dNr[i];
      dN[i].s = dNs[i];
    }
    return;
  }

  // Compute the numerators and denominator for the rational
  // Bézier functions and derivatives.
  double w = 0, dws = 0, dwr = 0, wi;
  
  for(int i = 0; i < NumNode; ++i)
  {
    wi = ElmNode[i]->GetW( );
    N[i] *= wi; 
    w += N[i];

    dN[i].r = dNr[i] * wi;
    dN[i].s = dNs[i] * wi;

    dwr += dN[i].r;
    dws += dN[i].s;
  }

  // Divide by the denominators to complete the computation of the
  // functions and derivatives w.r.t. the parent domain.

  for(int i = 0; i < NumNode; i++)
  {
    N[i] /= w; 
    dN[i].r = (dN[i].r - N[i]*dwr) / w;
    dN[i].s = (dN[i].s - N[i]*dws) / w;
  }
}

// ================================ GetEdge ================================

int cShapePlBezTrianSurf :: GetEdge(int *corner, eShpType *type, int *nnode, cNode **conn)
{
  // Barycentric index data

  int deg = Degree;
  int id1 = 0;                  // id p00
  int id2 = NumNode - 1 - deg;  // id 0p0
  int id3 = NumNode - 1;        // id 00p

  // Find valid vector index of corner control points

  int vecid[2] = {-1,-1};

  if (corner[0] == ElmNode[id1]->GetLabel( ))      // first corner
    vecid[0] = id1;
  else if (corner[0] == ElmNode[id2]->GetLabel( ))
    vecid[0] = id2;
  else if (corner[0] == ElmNode[id3]->GetLabel( ))
    vecid[0] = id3;

  if (corner[1] == ElmNode[id1]->GetLabel( ))      // second corner
    vecid[1] = id1;
  else if (corner[1] == ElmNode[id2]->GetLabel( ))
    vecid[1] = id2;
  else if (corner[1] == ElmNode[id3]->GetLabel( ))
    vecid[1] = id3;

  // Abort if do not find 

  if (vecid[0] == -1 || vecid[1] == -1)
  {
    cout << "Invalid control point label in Area Load in Bézier Triangle Shape!\n";
    return 0;
    //exit(0);
  }

  // Get barycentric indexes

  int cp1_i = 0, cp1_j = 0, cp2_i = 0, cp2_j = 0;

  if (vecid[0] == id1)
    cp1_i = deg;
  else if (vecid[0] == id2)
    cp1_j = deg;

  if (vecid[1] == id1)
    cp2_i = deg;
  else if (vecid[1] == id2)
    cp2_j = deg;


  // Store load shape nodes

  int i, j, indice;

  for(int k = 0; k < deg+1; ++k)
  {
    if (k==0)
    {
      i = cp1_i;
      j = cp1_j;
    }
    else // Modify indexes
    {
      if (cp1_i > 0)
        --i;
      if (cp1_j > 0)
        --j;
      if (cp2_i > 0)
        ++i;
      if (cp2_j > 0)
        ++j;
    }

    indice = nBernsteinBasis :: MapTrianBaryToVecID(i,j,deg);
    conn[k] = ElmNode[indice]; 
  }

  *type  = SHAPE_BEZ_PLCURV;
  *nnode = deg+1;  

  return 1;
}

// -------------------------------------------------------------------------
// Class cShapeBezTrianSurf:
//

// ============================ cShapeBezSurf ==============================

cShapeBezTrianSurf :: cShapeBezTrianSurf(void)
{
  Type = SHAPE_BSP_SURF;
}

// =========================== ~cShapeBezTrianSurf =========================

cShapeBezTrianSurf :: ~cShapeBezTrianSurf(void)
{
}

// -------------------------------------------------------------------------
// Class cShapePlBezSurf:
// -------------------------------------------------------------------------

// =========================== cShapePlBezSurfIGA ==========================

cShapePlBezSurf :: cShapePlBezSurf(void)
{
  Type    = SHAPE_BEZ_PLSURF;
  TopType = QUADRILATERAL_TOPOLOGY;
  ElmNode = 0;
}

// ========================== ~cShapePlSurfIGA =============================

cShapePlBezSurf :: ~cShapePlBezSurf(void)
{
}

// ================================= Read ==================================

void cShapePlBezSurf :: Read(void)
{
  // Read Rational flag and degree.
  if (!(in >> Rational) || !(in >> Degree[0]) || !(in >> Degree[1]))
  {
    cout << "Error in the input of bézier curve shape!\n";
    exit(0);
  }

  // Setting FEM node input data.
  NumNode = (Degree[0] + 1) * (Degree[1] + 1);
  ElmNode = new cNode* [NumNode];
  
  // Read each control point.
  cShape :: Read( );
}

// ============================== NodeNatCoord =============================

void cShapePlBezSurf :: NodeNatCoord(sNatCoord *c)
{
  double r,s;
  int counter = 0;
  
  int dr = Degree[0];
  int ds = Degree[1];
     
  for(int j = 0; j < ds+1; j++)
    for(int i = 0; i < dr+1; i++)
    {
      r = -1.0 + (2.0/dr) * i;
      s = -1.0 + (2.0/ds) * j;

      c[counter].r = r;
      c[counter].s = s;
      counter++;
    } 
}

// ================================= ShpFunc ===============================

void cShapePlBezSurf :: ShpFunc(sNatCoord p, double *N)
{
  // Map parameter coordinate from [-1,1] to [0,1].
  double r = (p.r + 1.0)/2.0;
  double s = (p.s + 1.0)/2.0;

  int nbr = Degree[0] + 1;
  int nbs = Degree[1] + 1;

  double *Nr = new double [nbr];
  double *Ns = new double [nbs];
 
  // Compute the univariate B-Spline basis functions w.r.t. the parent domain.
  nBernsteinBasis :: CompPols(nbr-1,r,Nr);
  nBernsteinBasis :: CompPols(nbs-1,s,Ns);

  // Compute the numerators and denominator for the tensor product
  // Rational Beźier functions.
  int a = 0;
  double w = 0;

  for(int j = 0; j < nbs; j++)
    for(int i = 0; i < nbr; i++)
    {
      N[a] = Nr[i] * Ns[j] * ElmNode[a]->GetW( );
      w += N[a];

      a++;
    }

  // Divide by the denominators to complete the computation of the
  // functions and derivatives w.r.t. the parent domain.

  for(int i = 0; i < a; i++)
    N[i] /= w; 

  // Release memory.
  delete [] Nr;
  delete [] Ns;
}

// =============================== DrvShpRST ===============================

void cShapePlBezSurf :: DrvShpRST(sNatCoord p, sNatCoord *dN)
{ 
  // Map parameter coordinate from [-1,1] to [0,1].
  double r = (p.r + 1.0)/2.0;
  double s = (p.s + 1.0)/2.0;

  int nbr = Degree[0] + 1;
  int nbs = Degree[1] + 1;

  double *N   = new double [nbr*nbs];
  double *dNr = new double [nbr];
  double *dNs = new double [nbs];
  double *Nr  = new double [nbr];
  double *Ns  = new double [nbs];

  // Compute the univariate B-Spline basis functions and derivatives w.r.t. 
  // the parent domain.

  nBernsteinBasis :: CompPols(nbr-1,r,Nr);
  nBernsteinBasis :: CompPols(nbs-1,s,Ns);
  nBernsteinBasis :: CompDervPols(nbr-1,r,dNr);
  nBernsteinBasis :: CompDervPols(nbs-1,s,dNs);

  // Apply jacobian transformation to derivatives.
  for(int i = 0; i < nbr; ++i)
  {
    dNr[i] /= 2.0;
    dNs[i] /= 2.0;
  }

  // Compute B-spline Shape function in the case of BSPLINE_2D patch.
  int a = 0;
  if (!Rational)
  {
    for(int j = 0; j < nbs; j++)
      for(int i = 0; i < nbr; i++)
      {
        N[a] = Nr[i] * Ns[j];
    
        dN[a].s = dNs[j] * Nr[i];
        dN[a].r = dNr[i] * Ns[j];
      }

     return;
  }

  // Compute the numerators and denominator for the tensor product
  // rational Bézier functions and derivatives.

  a = 0;
  double w = 0, dws = 0, dwr = 0;
  
  for(int j = 0; j < nbs; j++)
    for(int i = 0; i < nbr; i++)
    {

      N[a] = Nr[i] * Ns[j] * ElmNode[a]->GetW( );
      w += N[a];

      dN[a].r = dNr[i] *  Ns[j] * ElmNode[a]->GetW( );
      dN[a].s =  Nr[i] * dNs[j] * ElmNode[a]->GetW( );

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

// ================================ GetEdge ================================

int cShapePlBezSurf :: GetEdge(int *corner, eShpType *type, int *nnode, cNode **conn)
{
  // Get patch data.
  int degR = Degree[0];
  int degS = Degree[1];
  int size = (degR + 1) * (degS + 1);

  // Corner indices.
  int id_bl = ElmNode[0            ]->GetLabel( ); 
  int id_br = ElmNode[degR         ]->GetLabel( ); 
  int id_tl = ElmNode[size -1 -degR]->GetLabel( ); 
  int id_tr = ElmNode[size-1       ]->GetLabel( ); 

  int stride;
  int first;
  bool reverse = false;

  // Case Bottom.
  if ((corner[0] == id_bl && corner[1] == id_br) || 
      (corner[0] == id_br && corner[1] == id_bl)) 
  {
    *nnode   = degR + 1;
     first   = 0; 
     stride  = 0;
  }

  // Case Top.
  else if ((corner[0] == id_tl && corner[1] == id_tr) || 
           (corner[0] == id_tr && corner[1] == id_tl))
  {
     *nnode = degR + 1;
     first  = size - 1 - degR; 
     stride = 0;
     reverse = true;
  }
  
  // Case Left.
  else if ((corner[0] == id_bl && corner[1] == id_tl) || 
           (corner[0] == id_tl && corner[1] == id_bl))
  {
     *nnode  = degS + 1;
     first   = 0; 
     stride  = degR;
     reverse = true;
  }
  
  // Case Right.
  else if ((corner[0] == id_br && corner[1] == id_tr) || 
           (corner[0] == id_tr && corner[1] == id_br))
  {
     *nnode = degS + 1;
     first  = degR;
     stride = degR;
  }
  else
    return 0;
  // Get control points reference.

  if (reverse)
    for(int i = 0; i < *nnode; ++i)
      conn[*nnode - 1 -i] = ElmNode[first + (stride + 1) * i];
  else
    for(int i = 0; i < *nnode; ++i)
      conn[i] = ElmNode[first + (stride + 1) * i];

  // Set shape type.
  *type  = SHAPE_BEZ_PLCURV;

  return 1;
}

// -------------------------------------------------------------------------
// Class cShapeBezSurf:
//

// ============================ cShapeBezSurf ==============================

cShapeBezSurf :: cShapeBezSurf(void)
{
  Type = SHAPE_BSP_SURF;
}

// =========================== ~cShapeBezSurf ==============================

cShapeBezSurf :: ~cShapeBezSurf(void)
{
}

// ======================================================= End of file =====
