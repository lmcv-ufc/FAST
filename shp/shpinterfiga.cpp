// -------------------------------------------------------------------------
// shpinterfiga.cpp - implementation of the Shape class.
// -------------------------------------------------------------------------
// Created:      16-Dez-2013     Edson Moreira Dantas Junior
//
// Modified:     06-May-2014     Evandro Parente Junior
//               Creation of two-dimensional interface shapes.
// -------------------------------------------------------------------------

#include <math.h>
#include <stdio.h>
#include <iostream>
#include <cstdlib>
using namespace std;

//#include "patch.h"
#include "shpinterfiga.h"
#include "bspbasis.h"
#include "ctrlpnt.h"
#include "ctrl.h"
#include "vec.h"
#include "gblvar.h"

/*
static sNatCoord _dMrst[50];    // Mapping func derivatives/parametric coords

// =============================== RMatrix =================================

void cShapeInterf2D :: RMatrix(sNatCoord pt, sNodeCoord *coord,
                               double *detJ, cMatrix &R)
{
  int i;

  // Compute dMrs.

  DrvMapRST(pt, _dMrst);

  // Evaluate the tangent vector {t} = d{x}/dr.

  double tx = 0.0;
  double ty = 0.0;
  for (i = 0; i < NumNode; i++)
  {
    tx += _dMrst[i].r*coord[i].x/2.0;
    ty += _dMrst[i].r*coord[i].y/2.0;
  }

  // Normalize the tangent vector.

  double len = sqrt(tx*tx + ty*ty);
  if (len != 0.0)
  {
    tx /= len;

ty /= len;
  }
  *detJ = len;

  // Assembly the rotation matrix (vectors in rows).

  R[0][0] =  tx;  R[0][1] = ty;
  R[1][0] = -ty;  R[1][1] = tx;
}

*/
// -------------------------------------------------------------------------
// Class cShapeInterflineBsp:
// -------------------------------------------------------------------------

// ============================ cShapeInterflineBsp =============================

cShapeInterflineBsp :: cShapeInterflineBsp(void)
{
  Type = SHAPE_BSP_INTERF_LINE;
  pat = 0;
  ElmNode = 0;
}

// =========================== ~cShapeInterflineBsp =============================

cShapeInterflineBsp :: ~cShapeInterflineBsp(void)
{
}

// ============================== Read ======================================

void cShapeInterflineBsp :: Read(void)
{
  // Read active parametric direction.
  in >> ActParDir;

  // Read patch data.
  ReadIGA(  );
}

// ============================== Init =====================================

void cShapeInterflineBsp :: Init(void)
{
  // Check if correct patch type was stored.
  cBspSurface* pat2d = dynamic_cast<cBspSurface*>(pat);

  if (!pat2d)
  {
    cout << "Error in the input of patch used by IGA plane surface element!\n";
    exit(0);
  }

  // Compute the shape function basis index w.r.t parent patch domain.
  int* BasInd[2];
  NumBas[0] = pat->getBasis(0).getDegree( ) + 1;
  NumBas[1] = pat->getBasis(1).getDegree( ) + 1;

  BasInd[0] = new int [NumBas[0]];
  BasInd[1] = new int [NumBas[1]];

  // Get local basis function index.
  if (ActParDir == 0)
  {
    pat->getBasInd(0,SpanInd[0],   BasInd[0]);
    pat->getBasInd(1,SpanInd[1]-1, BasInd[1]);
  }
  else
  {
    pat->getBasInd(0,SpanInd[0]-1, BasInd[0]);
    pat->getBasInd(1,SpanInd[1],   BasInd[1]);
  }

  // Get control points of the element.
  NumNode = 2 * NumBas[ActParDir];
  ElmNode = new cNode* [NumNode];

  int counter = 0;
  int patid;
  int nodid;
  int id;
  cNode  *node;

  if (ActParDir == 0) // Point in
  {
    // First line;
    for(int i = 0; i < NumBas[0]; i++)
    {
      // Get patch control point.
      patid = pat2d->getCPInd(BasInd[0][i],BasInd[1][NumBas[1]-1]);

      if (!patid)
      {
        cout << "Error in the input of isogeometric interface element!\n";
        cout << "Invalid patch knot span.\n";
	exit(0);
      }

      // Get patch find node.
      nodid = pat->getCtrlPnt(patid)->getLabel( );
      node  = cNode :: FindNode(nodid);

      if (!node)
      {
        cout << "Error in the input of isogeometric interface element!\n";
        cout << "Invalid patch knot span.\n";
	exit(0);
      }

      ElmNode[counter++] = cNode :: FindNode(nodid);

    }

    // Second line;
    pat->getBasInd(1,SpanInd[1], BasInd[1]);

    for(int i = 0; i < NumBas[0]; i++)
    {
      id    = pat2d->getCPInd(BasInd[0][i],BasInd[1][0]);
      nodid = pat->getCtrlPnt(id)->getLabel( );
      ElmNode[counter++] = cNode :: FindNode(nodid);
    }
  }
  else
  {
    // First line;
    for(int i = 0; i < NumBas[1]; i++)
    {
      id = pat2d->getCPInd(BasInd[0][NumBas[1]-1],BasInd[1][i]);
      nodid = pat->getCtrlPnt(id)->getLabel( );
      ElmNode[counter++] = cNode :: FindNode(nodid);
    }

    // Second line;
    pat->getBasInd(0,SpanInd[0],BasInd[0]);

    for(int i = 0; i < NumBas[1]; i++)
    {
      id = pat2d->getCPInd(BasInd[0][0],BasInd[1][i]);
    nodid = pat->getCtrlPnt(id)->getLabel( );
      ElmNode[counter++] = cNode :: FindNode(nodid);
    }
  }
}

// ============================= NodeNatCoord ==============================

void cShapeInterflineBsp :: NodeNatCoord(sNatCoord *c)
{
  if (cControl :: GetIgaEqvMesh( ))
  {
    c[0].r = c[3].r = -1.0;
    c[1].r = c[4].r =  0.0;
    c[2].r = c[5].r =  1.0;
  }
  else
  {
    //double r;
    int counter = 0;

     for(int i = 0; i < NumBas[ActParDir]; i++)
     {
       // ???
       //r = -1.0+(2.0/(NumBas[ActParDir]-1))*i;

       c[counter].r = c[counter+3].r;
       counter++;
     }
  }
}


// ================================= ShpFunc ===============================

void cShapeInterflineBsp :: ShpFunc(sNatCoord p, double *N)
{

//NÃO TENHO CERTEZA SE ISSO ESTÁ CORRETO,
//VERIFICAR ISSO PARA O CASO DO ELEMENTO DE INTERFACE ACIMA.
//FUNÇÕES DE FORMA DO ELEMENTO DE INTERFACE ENTENDER O FUNCIONAMENTO ANTES DE MODIFICAR PARA

  double r = p.r;
  //double s = p.s;

  cVector Nr(NumBas[ActParDir]);//;, Ns(NumBas[1]);
  Nr.Zero( );

  // Compute the univariate B-Spline basis functions w.r.t. the parent domain.
  EvalBspBas(ActParDir,r,Nr);

  // Compute B-spline Shape function in the case of BSPLINE_2D patch.
  int a = 0;
  if (!pat->isRational( ))
  {
    for(int i = 0; i < NumBas[ActParDir]; i++)
      N[a] = N[a+NumBas[ActParDir]] = Nr[i];
     return;
  }

  // Compute the numerators and denominator for the tensor product
  // NURBS functions and derivatives.

  a = 0;
  double w = 0;

  for(int i = 0; i < NumBas[0]; i++)
    {
      N[a] = Nr[i] * ElmNode[a]->GetW( );
      w += N[a];
      a++;
    }

  // Divide by the denominators to complete the computation of the
  // functions and derivatives w.r.t. the parent domain.

  for(int i = 0; i < a; i++)
  {
    N[i] /= w;
    N[i+NumBas[ActParDir]] = N[i];
  }


}

// =============================== DrvShpRST ===============================

void cShapeInterflineBsp :: DrvShpRST(sNatCoord p, sNatCoord *dN)
{
  double r = p.r;

  cVector N(NumNode);
  cVector dNr(NumBas[0]);
  cVector Nr(NumBas[0]);

  N.Zero( );
  Nr.Zero( );
  dNr.Zero( );

  // Compute the univariate B-Spline basis functions and derivatives w.r.t.
  // the parent domain.

  EvalBspBas(ActParDir,r,Nr);
  EvalBspBasDerv(ActParDir,r,dNr);

  // Compute B-spline Shape function in the case of BSPLINE_1D patch.

  if (!pat->isRational( ))
  {
    for(int i = 0; i < NumBas[ActParDir]; i++)
    {
      N[i]    = N[i+NumBas[ActParDir]]    = Nr[i];
      dN[i].r = dN[i+NumBas[ActParDir]].r = dNr[i];
    }

    return;
  }

  // Compute the numerators and denominator for the NURBS functions.

  double w = 0, dwr = 0;

  for(int i = 0; i < NumBas[ActParDir]; i++)
  {
    N[i] = Nr[i] * ElmNode[i]->GetW( );
    w += N[i];

    dN[i].r = dNr[i] * ElmNode[i]->GetW( );
    dwr += dN[i].r;
  }

  // Divide by the denominators to complete the computation of the
  // functions w.r.t. the parent domain.

  for(int i = 0; i < NumBas[ActParDir]; i++)
  {
    N[i] /= w;
    dN[i].r = (dN[i].r - N[i]*dwr) / w;
    dN[i+NumBas[ActParDir]].r = dN[i].r;
   // cout << i << " " << dN[i].r << " ";
    //cout << (i+NumBas[ActParDir]) << " " << dN[i+NumBas[ActParDir]].r;
  }
//cout << endl;


}

/*
// -------------------------------------------------------------------------
// Class cShapeInterfL3:
// -------------------------------------------------------------------------

// ============================= cShapeInterfL3 ============================

cShapeInterfL3 :: cShapeInterfL3(void)
{
  Type = SHAPE_INTERF_L3;
  NumNode = 6;
  ElmNode = new cNode*[NumNode];
}

// ============================ ~cShapeInterfL3 ============================

cShapeInterfL3 :: ~cShapeInterfL3(void)
{
}

// ============================== NodeNatCoord =============================

void cShapeInterfL3 :: NodeNatCoord(sNatCoord *c)
{
  c[0].r = c[3].r = -1.0;
  c[1].r = c[4].r =  0.0;
  c[2].r = c[5].r =  1.0;
}

// ================================= ShpFunc ===============================

void cShapeInterfL3 :: ShpFunc(sNatCoord p, double *N)
{
  double r = p.r;

  N[0] = N[3] = 0.5*r*(r - 1.0);
  N[1] = N[4] = 1.0 - r*r;
  N[2] = N[5] = 0.5*r*(r + 1.0);
}

// =============================== DrvShpRST ===============================

void cShapeInterfL3 :: DrvShpRST(sNatCoord p, sNatCoord *dN)
{
  double r = p.r;

  dN[0].r = dN[3].r =  0.5*(2.0*r - 1.0);
  dN[1].r = dN[4].r = -2.0*r;
  dN[2].r = dN[5].r =  0.5*(2.0*r + 1.0);
}


// -------------------------------------------------------------------------
// Class cShapeInterf3D:
// -------------------------------------------------------------------------

// ============================= cShapeInterf3D ============================

cShapeInterf3D :: cShapeInterf3D(void)
{
}

// ============================ ~cShapeInterf3D ============================

cShapeInterf3D :: ~cShapeInterf3D(void)
{
}

// =============================== RMatrix =================================

void cShapeInterf3D :: RMatrix(sNatCoord pt, sNodeCoord *coord,
                               double *detJ, cMatrix &R)
{
  int i;

  // Compute dMrs.

  DrvMapRST(pt, _dMrst);

  // Evaluate the tangent vectors: {u} = d{x}/dr and {v} = d{x}/ds.

  cVector u(3),v(3);
  u.Zero( );
  v.Zero( );
  for (i = 0; i < NumNode; i++)
  {
    u[0] += _dMrst[i].r*coord[i].x/2.0;
    u[1] += _dMrst[i].r*coord[i].y/2.0;
    u[2] += _dMrst[i].r*coord[i].z/2.0;

    v[0] += _dMrst[i].s*coord[i].x/2.0;
    v[1] += _dMrst[i].s*coord[i].y/2.0;
    v[2] += _dMrst[i].s*coord[i].z/2.0;
  }

  // Evaluate the normal vector {w} = {u} x {v}.

  cVector w(3);
  CrossProd(u, v, w);
  *detJ = w.Length( );

  // Ensure that the 3 vectors are orthogonal to each other.

  CrossProd(w, u, v); // {v} = {w} x {u}

  // Assembly the rotation matrix (vectors in rows).

  u.Normalize( );
  v.Normalize( );
  w.Normalize( );
  for (i = 0; i < 3; i++)
  {
    R[0][i] = u[i];
    R[1][i] = v[i];
    R[2][i] = w[i];
  }
}


// -------------------------------------------------------------------------
// Class cShapeInterfQ4:
// -------------------------------------------------------------------------

// ============================ cShapeInterfQ4 =============================

cShapeInterfQ4 :: cShapeInterfQ4(void)
{
  Type = SHAPE_INTERF_Q4;
  NumNode = 8;
  ElmNode = new cNode*[NumNode];
}

// =========================== ~cShapeInterfQ4 =============================

cShapeInterfQ4  :: ~cShapeInterfQ4(void)
{
}

// ============================== NodeNatCoord =============================

void cShapeInterfQ4 :: NodeNatCoord(sNatCoord *c)
{
  c[0].r = c[4].r = -1.0;  c[0].s = c[4].s = -1.0;
  c[1].r = c[5].r =  1.0;  c[1].s = c[5].s = -1.0;
  c[2].r = c[6].r =  1.0;  c[2].s = c[6].s =  1.0;
  c[3].r = c[7].r = -1.0;  c[3].s = c[7].s =  1.0;
}

// ================================= ShpFunc ===============================

void cShapeInterfQ4 :: ShpFunc(sNatCoord p, double *N)
{
  double r = p.r;
  double s = p.s;

  N[0] = N[4] = 0.25*(1.0 - r)*(1.0 - s);
  N[1] = N[5] = 0.25*(1.0 + r)*(1.0 - s);
  N[2] = N[6] = 0.25*(1.0 + r)*(1.0 + s);
  N[3] = N[7] = 0.25*(1.0 - r)*(1.0 + s);
}

// =============================== DrvShpRST ===============================

void cShapeInterfQ4 :: DrvShpRST(sNatCoord p, sNatCoord *dN)
{
  double r = p.r;
  double s = p.s;

  dN[0].r = dN[4].r = -0.25*(1.0 - s);
  dN[1].r = dN[5].r =  0.25*(1.0 - s);
  dN[2].r = dN[6].r =  0.25*(1.0 + s);
  dN[3].r = dN[7].r = -0.25*(1.0 + s);

  dN[0].s = dN[4].s = -0.25*(1.0 - r);
  dN[1].s = dN[5].s = -0.25*(1.0 + r);
  dN[2].s = dN[6].s =  0.25*(1.0 + r);
  dN[3].s = dN[7].s =  0.25*(1.0 - r);
}


// -------------------------------------------------------------------------
// Class cShapeInterfQ8:
// -------------------------------------------------------------------------

// ============================ cShapeInterfQ8 =============================

cShapeInterfQ8 :: cShapeInterfQ8(void)
{
  Type = SHAPE_INTERF_Q8;
  NumNode = 16;
  ElmNode = new cNode*[NumNode];
}

// =========================== ~cShapeInterfQ8 =============================

cShapeInterfQ8  :: ~cShapeInterfQ8(void)
{
}

// ============================== NodeNatCoord =============================

void cShapeInterfQ8 :: NodeNatCoord(sNatCoord *c)
{
  c[0].r = c[ 8].r = -1.0;  c[0].s = c[ 8].s = -1.0;
  c[1].r = c[ 9].r =  0.0;  c[1].s = c[ 9].s = -1.0;
  c[2].r = c[10].r =  1.0;  c[2].s = c[10].s = -1.0;
  c[3].r = c[11].r =  1.0;  c[3].s = c[11].s =  0.0;
  c[4].r = c[12].r =  1.0;  c[4].s = c[12].s =  1.0;
  c[5].r = c[13].r =  0.0;  c[5].s = c[13].s =  1.0;
  c[6].r = c[14].r = -1.0;  c[6].s = c[14].s =  1.0;
  c[7].r = c[15].r = -1.0;  c[7].s = c[15].s =  0.0;
}

// ================================= ShpFunc ===============================

void cShapeInterfQ8 :: ShpFunc(sNatCoord p, double *N)
{
  double r = p.r;
  double s = p.s;

  N[0] = 0.25*(1.0 - r)*(1.0 - s) - 0.25*(1.0 - r)*(1.0 - s*s) -
         0.25*(1.0 - r*r)*(1.0 - s);
  N[1] = 0.50*(1.0 - r*r)*(1.0 - s);
  N[2] = 0.25*(1.0 + r)*(1.0 - s) - 0.25*(1.0 - r*r)*(1.0 - s) -
         0.25*(1.0 - s*s)*(1.0 + r);
  N[3] = 0.50*(1.0 - s*s)*(1.0 + r);
  N[4] = 0.25*(1.0 + r)*(1.0 + s) - 0.25*(1.0 - s*s)*(1.0 + r) -
         0.25*(1.0 + s)*(1.0 - r*r);
  N[5] = 0.50*(1.0 - r*r)*(1.0 + s);
  N[6] = 0.25*(1.0 - r)*(1.0 + s) - 0.25*(1.0 - r*r)*(1.0 + s) -
         0.25*(1.0 - s*s)*(1.0 - r);
  N[7] = 0.50*(1.0 - s*s)*(1.0 - r);

  N[ 8] = N[0];
  N[ 9] = N[1];
  N[10] = N[2];
  N[11] = N[3];
  N[12] = N[4];
  N[13] = N[5];
  N[14] = N[6];
  N[15] = N[7];
}

// =============================== DrvShpRST ===============================

void cShapeInterfQ8 :: DrvShpRST(sNatCoord p, sNatCoord *dN)
{
  double r = p.r;
  double s = p.s;

  dN[0].r = dN[ 8].r = (2.0*r - 2.0*r*s - s*s + s)/4.0;
  dN[1].r = dN[ 9].r = r*s - r;
  dN[2].r = dN[10].r = (2.0*r - 2.0*r*s + s*s - s)/4.0;
  dN[3].r = dN[11].r = (1.0 - s*s)/2.0;
  dN[4].r = dN[12].r = (2.0*r + 2.0*r*s + s*s + s)/4.0;
  dN[5].r = dN[13].r = -r - r*s;
  dN[6].r = dN[14].r = (2.0*r + 2.0*r*s - s*s - s)/4.0;
  dN[7].r = dN[15].r = (s*s - 1.0)/2.0;

  dN[0].s = dN[ 8].s = (2.0*s - r*r - 2.0*r*s + r)/4.0;
  dN[1].s = dN[ 9].s = (r*r - 1.0)/2.0;
  dN[2].s = dN[10].s = (2.0*s - r*r + 2.0*r*s - r)/ 4.0;
  dN[3].s = dN[11].s = -s - r*s;
  dN[4].s = dN[12].s = (2.0*s + r*r + 2.0*r*s + r)/4.0;
  dN[5].s = dN[13].s = (1.0 - r*r)/2.0;
  dN[6].s = dN[14].s = (2.0*s + r*r - 2.0*r*s - r)/4.0;
  dN[7].s = dN[15].s = r*s - s;
}

// ======================================================= End of file =====*/
