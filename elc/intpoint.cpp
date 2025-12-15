// -------------------------------------------------------------------------
// intpoint.cpp - implementation of the Integration Point class.
// -------------------------------------------------------------------------
// Created:      05-May-2005     Evandro Parente Junior
//
// Modified:     10-Nov-2012     Evandro Parente Junior
//               Implementation of Lobatto integration.
//
// Modified:     08-Feb-2015     Evandro Parente Junior
//               Gauss integration for n = 5 to 15.
//
// Modified:     08-Apr-2015     Elias Saraiva Barroso
//               Input/output operations using c++ streams.
//
// Modified:     18-Sep-2019     Elias Saraiva Barroso
//               Added Dunavant Triangular quarature data (p 1-20).
//
// Modified:     15-Jun-2023     Evandro Parente Junior
//               Added integration rules for tetrahedra.
// -------------------------------------------------------------------------

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <vector>

using namespace std;

#include "intpoint.h"
#include "utl.h"

static eQuadType TrianQuad = DUNAVANT;

// -------------------------------------------------------------------------
// Auxiliary types:
//
struct sIntData
{
  const double *r;
  const double *s;
  const double *t;
  const double *w;
  const int     np;

  sIntData(const double *_r,const double *_s,const double *_t,const double *_w,const int _np) : np(_np)
  {
    r = _r;
    s = _s;
    t = _t;
    w = _w;
  }
};

// -------------------------------------------------------------------------
// Quadrature data functions:
//
static const vector<sIntData>& getGaussData( );
static const vector<sIntData>& getLobattoData( );
static const vector<sIntData>& getNewtonCotesData( );
static const vector<sIntData>& getDunavantData( );
static const vector<sIntData>& getLynessData( );
static const vector<sIntData>& getTetData( );

// -------------------------------------------------------------------------
// Auxiliary methods:
//

// =========================== GetIntPntDataSize ===========================

int GetIntPntDataSize(eQuadType type)
{
  int size = 0;

  switch (type)
  {
    case GAUSS:
      size = getGaussData( ).size( );
    break;

    case LOBATTO:
      size = getLobattoData( ).size( );
    break;

    case NEWTON_COTES:
      size = getNewtonCotesData( ).size( );
    break;

    case DUNAVANT:
      size = getDunavantData( ).size( );
    break;

    case LYNESS:
      size = getLynessData( ).size( );
    break;

    case TETRAHEDRA:
      size = getTetData( ).size( );
    break;
  }

  return size;
}

// =========================== GetIntPntData ===============================

const vector<sIntData>* GetIntPntData(eQuadType type)
{
  const vector<sIntData> *Data = 0;

  switch (type)
  {
    case GAUSS:
      Data = &getGaussData( );
    break;

    case LOBATTO:
      Data = &getLobattoData( );
    break;
    
    case NEWTON_COTES:
      Data = &getNewtonCotesData( );
    break;

    case DUNAVANT:
      Data = &getDunavantData( );
    break;

    case LYNESS:
      Data = &getLynessData( );
    break;

    case TETRAHEDRA:
      Data = &getTetData( );
    break;
  }

  return Data;
}

// ======================== operator>> (eQuadType) =========================

istream& operator>> (istream &in, eQuadType &type)
{
  // Read the patch type label.

  char label[100];
  if (!Utl::ReadString(in, label))
    Utl::Exit("Error in the input of the quadrature type label");

  // Set the appropriate patch type.

  if (string(label) == "Gauss")
    type = GAUSS;
  else if (string(label) == "Lobatto")
    type = LOBATTO;
  else if (string(label) == "Newton-Cotes")
    type = NEWTON_COTES;
  else if (string(label) == "Dunavant")
    type = DUNAVANT;
  else if (string(label) == "Lyness")
    type = LYNESS;
  else if (string(label) == "Lyness")
    type = TETRAHEDRA;
  else
    Utl::Exit("Unknown quadrature type: " + string(label));

  return in;
}

// -------------------------------------------------------------------------
// Public methods:
//
  
// =========================== SetTrianQuad ================================

void cIntPoint :: SetTrianQuad(eQuadType t)
{
  TrianQuad = t; 
}

// =========================== CreateIntPoints =============================

cIntPoint *cIntPoint :: CreateIntPoints(eTopType type, int quadtype, int order[], int *n)
{
//  cout << "cIntPoint :: CreateIntPoints\n";
  eQuadType qtype = static_cast<eQuadType>(quadtype);

  // Correction to valid quadrature. Evandro: Não entendi.
  if (type != TRIANGULAR_TOPOLOGY)
    if (qtype == DUNAVANT || qtype == LYNESS)
      qtype = GAUSS;

//  cout << "toptype = " << type << "  quadtype = " << quadtype << endl;
  switch(type)
  {
    case LINE_TOPOLOGY:
     return(CreateLinePoints(order[0], qtype, n));

 /* Evandro: por que tirou isso?
    case SHAPE_INTERF_L2:
    case SHAPE_INTERF_L3:
    case SHAPE_BSP_INTERF_LINE:
     return(CreateLinePoints(order[0], NEWTON_COTES, n));
     */
     //return(CreateNCLinePoints(order, n));

    case TRIANGULAR_TOPOLOGY:
     return(CreateTriaPoints(order[0], TrianQuad, n));

    case QUADRILATERAL_TOPOLOGY:
     return(CreateQuadPoints(order, qtype, n));

    case WEDGE_TOPOLOGY:
     return(CreateWedgePoints(order, qtype, n));

 /* Evandro: por que tirou isso?
    case SHAPE_INTERF_Q4:
    case SHAPE_INTERF_Q8:
     return(CreateQuadPoints(order, NEWTON_COTES, n));
     */
     //return(CreateNCQuadPoints(order, n));

    case HEXAHEDRAL_TOPOLOGY:
     return(CreateBrickPoints(order, qtype, n));

    case TETRAHEDRAL_TOPOLOGY:
     // Evandro: essa questão de tipo/topologia está confusa (Elias)!
     return(CreateTetPoints(order[0], TETRAHEDRA, n));
  }

  cout << "\nUnknown shape type in CreateIntPoints!\n";
  exit(0);
  return(NULL);
}

// =========================== CreateLinePoints ============================

cIntPoint *cIntPoint :: CreateLinePoints(int order, eQuadType qtype, int *n)
{
  int DataSize = GetIntPntDataSize(qtype);

  if (order < 1 || order > DataSize)
  {
    cout << "\nInvalid cube integration order (" << order << ")!\n";
    if (order > DataSize)
      cout << "\nThe maximum order supported is " << DataSize << endl;
    exit(0);
  }

  const vector<sIntData> *Data = GetIntPntData(qtype);
  *n = order;
  cIntPoint *pts = new cIntPoint[*n];
  int idr = order -1;
  for (int i = 0; i < order; i++)
  {
    pts[i].pnt.r = (*Data)[idr].r[i];
    pts[i].wgt   = (*Data)[idr].w[i];
  }

  return(pts);
}

// =========================== CreateQuadPoints ============================

cIntPoint *cIntPoint :: CreateQuadPoints(int order[], eQuadType qtype, int *n)
{
  int DataSize = GetIntPntDataSize(qtype);

  for(int i = 0; i < 2; i++)
  {
    if (order[i] < 1 || order[i] > DataSize)
    {
      cout << "\nInvalid cube integration order (" << order << ")!\n";
      if (order[i] > DataSize)
        cout << "\nThe maximum order supported is " << DataSize << endl;
      exit(0);
    }
  }

  const vector<sIntData> *Data = GetIntPntData(qtype);
  *n = order[0]*order[1];
  cIntPoint *pts = new cIntPoint[*n];
  int k = 0;
  int idr = order[0] - 1;
  int ids = order[1] - 1;

  for (int i = 0; i < order[1]; i++)
  {
    for (int j = 0; j < order[0]; j++)
    {
      pts[k].pnt.r = (*Data)[idr].r[j];
      pts[k].pnt.s = (*Data)[ids].r[i];
      pts[k].wgt   = (*Data)[idr].w[j] * (*Data)[ids].w[i];

      k++;
    }
  }

  return(pts);
}

// =========================== CreateTriaPoints ============================

cIntPoint *cIntPoint :: CreateTriaPoints(int order,eQuadType qtype, int *n)
{
  int DataSize = GetIntPntDataSize(qtype);

  if (order < 1 || order > DataSize)
  {
    cout << "\nInvalid cube integration order (" << order << ")!\n";
    if (order > DataSize)
      cout << "\nThe maximum order supported is " << DataSize << endl;
    exit(0);
  }

  const vector<sIntData> *Data = GetIntPntData(qtype);
  *n = (*Data)[order-1].np;
  cIntPoint *pts = new cIntPoint[*n];
  int id = order - 1;

  for (int i = 0; i < *n; i++)
  {
      pts[i].pnt.r = (*Data)[id].r[i];
      pts[i].pnt.s = (*Data)[id].s[i];
      pts[i].wgt   = (*Data)[id].w[i];
  }

  return(pts);

/*
  cIntPoint *pts;


  if (order == 1)
  {
    *n = 1;
    pts = new cIntPoint[*n];
    pts[0].pnt.r = 0.333333333333333;
    pts[0].pnt.s = 0.333333333333333;
    pts[0].wgt   = 0.500000000000000;
  }
  else if (order == 3)
  {
    *n = 3;
    pts = new cIntPoint[*n];
    pts[0].pnt.r = 0.166666666666666;
    pts[0].pnt.s = 0.166666666666666;
    pts[0].wgt   = 0.166666666666666;

    pts[1].pnt.r = 0.166666666666666;
    pts[1].pnt.s = 0.666666666666666;
    pts[1].wgt   = 0.166666666666666;

    pts[2].pnt.r = 0.666666666666666;
    pts[2].pnt.s = 0.166666666666666;
    pts[2].wgt   = 0.166666666666666;
  }
  else if (order == 4)
  {
//    *n = 4;
//    pts = new cIntPoint[*n];
//    pts[0].pnt.r = 0.333333333333333;
//    pts[0].pnt.s = 0.333333333333333;
//    pts[0].wgt  = -0.562500000000000;
//
//    pts[1].pnt.r = 0.200000000000000;
//    pts[1].pnt.s = 0.200000000000000;
//    pts[1].wgt   = 0.520833333333333;
//
//    pts[2].pnt.r = 0.600000000000000;
//    pts[2].pnt.s = 0.200000000000000;
//    pts[2].wgt   = 0.520833333333333;
//
//    pts[3].pnt.r = 0.200000000000000;
//    pts[3].pnt.s = 0.600000000000000;
//    pts[3].wgt   = 0.520833333333333; 
    *n = 4;
    pts = new cIntPoint[*n];
    pts[0].pnt.r = 0.333333333333333;
    pts[0].pnt.s = 0.333333333333333;
    pts[0].wgt  = -0.281250000000000;

    pts[1].pnt.r = 0.200000000000000;
    pts[1].pnt.s = 0.200000000000000;
    pts[1].wgt   = 0.260416666666666;

    pts[2].pnt.r = 0.600000000000000;
    pts[2].pnt.s = 0.200000000000000;
    pts[2].wgt   = 0.260416666666666;

    pts[3].pnt.r = 0.200000000000000;
    pts[3].pnt.s = 0.600000000000000;
    pts[3].wgt   = 0.260416666666666;
  }
  else if (order == 7)
  {
    *n = 7;
    pts = new cIntPoint[*n];
    pts[0].pnt.r = 0.333333333333333;	
    pts[0].pnt.s = 0.333333333333333;	 
    pts[0].wgt   = 0.112500000000000;	

    pts[1].pnt.r = 0.059715871700000;
    pts[1].pnt.s = 0.470142064100000;
    pts[1].wgt   = 0.066197076400000;

    pts[2].pnt.r = 0.470142064100000;
    pts[2].pnt.s = 0.059715871700000;
    pts[2].wgt   = 0.066197076400000;

    pts[3].pnt.r = 0.470142064100000;
    pts[3].pnt.s = 0.470142064100000;
    pts[3].wgt   = 0.066197076400000;

    pts[4].pnt.r = 0.797426985300000;
    pts[4].pnt.s = 0.101286507300000;
    pts[4].wgt   = 0.062969590300000;

    pts[5].pnt.r = 0.101286507300000;
    pts[5].pnt.s = 0.797426985300000;
    pts[5].wgt   = 0.062969590300000;

    pts[6].pnt.r = 0.101286507300000;
    pts[6].pnt.s = 0.101286507300000;
    pts[6].wgt   = 0.062969590300000;
  }
  else if (order == 12)
  {
    double w[12] = { 
     0.050844906370206816921,
     0.050844906370206816921,
     0.050844906370206816921,
     0.11678627572637936603,
     0.11678627572637936603,
     0.11678627572637936603,
     0.082851075618373575194,
     0.082851075618373575194,
     0.082851075618373575194,
     0.082851075618373575194,
     0.082851075618373575194,
     0.082851075618373575194 };
   double xy[2*12] = { 
    0.87382197101699554332,  0.063089014491502228340, 
    0.063089014491502228340,  0.87382197101699554332, 
    0.063089014491502228340,  0.063089014491502228340, 
    0.50142650965817915742,  0.24928674517091042129, 
    0.24928674517091042129,  0.50142650965817915742, 
    0.24928674517091042129,  0.24928674517091042129, 
    0.053145049844816947353,  0.31035245103378440542, 
    0.31035245103378440542,  0.053145049844816947353, 
    0.053145049844816947353,  0.63650249912139864723, 
    0.31035245103378440542,  0.63650249912139864723, 
    0.63650249912139864723,  0.053145049844816947353, 
    0.63650249912139864723,  0.31035245103378440542 };


    *n = 12;
    pts = new cIntPoint[*n];
    for(int p = 0; p < 12; ++p)
    {
      pts[p].pnt.r = xy[p*2];
      pts[p].pnt.s = xy[p*2+1];
      pts[p].wgt   = w[p];
    }


  }
  else if (order == 20)
  {
    // Quadrature from High degree efficient symmetrical guassian quadrature rules for the triangle. International Journal for numerical methods in engineering. D. A. DUNAVANT (1985).
    double w[20] = {0.033057055541624, 0.000867019185663, 0.011660052716448, 0.022876936356421, 0.030448982673938, 0.030624891725355, 0.024368057676800, 0.015997432032024, 0.007698301815602, -0.000632060497488, 0.001751134301193, 0.016465839189576, 0.004839033540485, 0.025804906534650, 0.008471091054441, 0.018354914106280, 0.000704404677908, 0.010112684927462, 0.003573909385950};

    double r[20] = {0.333333333333333, -0.001900928704400, 0.023574084130543, 0.089726636099435, 0.196007481363421, 0.488214180481157, 0.647023488009788, 0.791658289326483, 0.893862072318140, 0.916762569607942, 0.976836157186356, 0.048741583664839, 0.006314115948605, 0.134316520547348, 0.013973893962392, 0.075549132909764, -0.008368153208227, 0.026686063258714, 0.010547719294141};
    double s[20] = {0.333333333333333, 0.500950464352200, 0.488212957934729, 0.455136681950283, 0.401996259318289, 0.255892909759421, 0.176488255995106, 0.104170855336758, 0.053068963840930, 0.041618715196029, 0.011581921406822, 0.344855770229001, 0.377843269594854, 0.306635479062357, 0.249419362774742, 0.212775724802802, 0.146965436053239, 0.137726978828923, 0.059696109149007};

    *n = 20;
    pts = new cIntPoint[*n];
    for(int p = 0; p < 20; ++p)
    {
      pts[p].pnt.r = r[p];
      pts[p].pnt.s = s[p];
      pts[p].wgt   = w[p];
    } 

  }
  else if (order == 28)
  {
    double r[28] = {0.333333333333333, 0.948021718143423, 0.025989140928288, 0.025989140928288, 0.811424994704155, 0.094287502647923, 0.094287502647923, 0.010726449965571, 0.494636775017215, 0.494636775017215, 0.585313234770972, 0.207343382614514, 0.207343382614514, 0.122184388599019, 0.438907805700491, 0.438907805700491, 0.677937654882590, 0.677937654882590, 0.044841677589131, 0.044841677589131, 0.277220667528279, 0.277220667528279, 0.858870281282636, 0.858870281282636, 0.000000000000000, 0.000000000000000, 0.141129718717364, 0.141129718717364};
    double s[28] = {0.333333333333333, 0.025989140928288, 0.948021718143423, 0.025989140928288, 0.094287502647923, 0.811424994704155, 0.094287502647923, 0.494636775017215, 0.010726449965571, 0.494636775017215, 0.207343382614514, 0.585313234770972, 0.207343382614514, 0.438907805700491, 0.122184388599019, 0.438907805700491, 0.044841677589131, 0.277220667528279, 0.677937654882590, 0.277220667528279, 0.677937654882590, 0.044841677589131, 0.000000000000000, 0.141129718717364, 0.858870281282636, 0.141129718717364, 0.858870281282636, 0.000000000000000};

    double w[28] = {0.043988650581111, 0.004372155776868, 0.004372155776868, 0.004372155776868, 0.019040785996968, 0.019040785996968, 0.019040785996968, 0.009427724028066, 0.009427724028066, 0.009427724028066, 0.036079848772370, 0.036079848772370, 0.036079848772370, 0.034664569352769, 0.034664569352769, 0.034664569352769, 0.020528157714644, 0.020528157714644, 0.020528157714644, 0.020528157714644, 0.020528157714644, 0.020528157714644, 0.003681191891650, 0.003681191891650, 0.003681191891650, 0.003681191891650, 0.003681191891650, 0.003681191891650};
	  
    *n = 28;
    pts = new cIntPoint[*n];
    for(int p = 0; p < 28; ++p)
    {
      pts[p].pnt.r = r[p];
      pts[p].pnt.s = s[p];
      pts[p].wgt   = w[p];
    }
  }
  else if (order == 37)
  {
    double r[37] = {0.3333, 0.9503, 0.0249, 0.0249, 0.1716, 0.4142, 0.4142, 0.5394, 0.2303, 0.2303, 0.7722, 0.1139, 0.1139, 0.0091, 0.4955, 0.4955, 0.0623, 0.4689, 0.4689, 0.0221, 0.0221, 0.8513, 0.8513, 0.1266, 0.1266, 0.0186, 0.0186, 0.6894, 0.6894, 0.2919, 0.2919, 0.0965, 0.0965, 0.6359, 0.6359, 0.2676, 0.2676};
    double s[37] = {0.3333, 0.0249, 0.9503, 0.0249, 0.4142, 0.1716, 0.4142, 0.2303, 0.5394, 0.2303, 0.1139, 0.7722, 0.1139, 0.4955, 0.0091, 0.4955, 0.4689, 0.0623, 0.4689, 0.8513, 0.1266, 0.0221, 0.1266, 0.0221, 0.8513, 0.6894, 0.2919, 0.0186, 0.2919, 0.0186, 0.6894, 0.6359, 0.2676, 0.0965, 0.2676, 0.0965, 0.6359};
    double w[37] = {0.0259, 0.004, 0.004, 0.004, 0.0234, 0.0234, 0.0234, 0.0233, 0.0233, 0.0233, 0.0155, 0.0155, 0.0155, 0.0054, 0.0054, 0.0054, 0.0161, 0.0161, 0.0161, 0.0077, 0.0077, 0.0077, 0.0077, 0.0077, 0.0077, 0.0089, 0.0089, 0.0089, 0.0089, 0.0089, 0.0089, 0.0185, 0.0185, 0.0185, 0.0185, 0.0185, 0.0185};

    *n = 37;
    pts = new cIntPoint[*n];
    for(int p = 0; p < 37; ++p)
    {
      pts[p].pnt.r = r[p];
      pts[p].pnt.s = s[p];
      pts[p].wgt   = w[p];
    }
  }

  else
  {
    cout << "\nInvalid triangular integration order (" << order << ")!\n";
    exit(0);
  }

  return(pts);
  */
}

// ========================== CreateWedgePoints ============================

cIntPoint *cIntPoint :: CreateWedgePoints(int order[], eQuadType type, int *n)
{
  // Get integration data in t direction.
  int DataSizeT = GetIntPntDataSize(type);

  if (order[2] < 1 || order[2] > DataSizeT)
  {
    cout << "\nInvalid cube integration order (" << order[2] << ")!\n";
    if (order[2] > DataSizeT)
      cout << "\nThe maximum order supported is " << DataSizeT << endl;
    exit(0);
  }
  const vector<sIntData> *DataT = GetIntPntData(type);

  // Get triangular data.
  int DataSize = GetIntPntDataSize(TrianQuad);

  if (order[0] < 1 || order[0] > DataSize)
  {
    cout << "\nInvalid cube integration order (" << order[0] << ")!\n";
    if (order[0] > DataSize)
      cout << "\nThe maximum order supported is " << DataSize << endl;
    exit(0);
  }

  const vector<sIntData> *Data = GetIntPntData(TrianQuad);

  int tqs        = (*Data)[order[0]-1].np;        // Triangle quadrature size.
  *n             = tqs * order[2];
  cIntPoint *pts = new cIntPoint[*n];

  // Get surface quadrature for each point in t direction.
  int idrs = order[0] - 1;
  int idt  = order[2] - 1;
  for (int i = 0; i < order[2]; i++)
    for (int j = 0; j < tqs; j++)
    {
        pts[tqs*i+j].pnt.r = (*Data)[idrs].r[j];
        pts[tqs*i+j].pnt.s = (*Data)[idrs].s[j];
        pts[tqs*i+j].pnt.t = (*DataT)[idt].r[i];
        pts[tqs*i+j].wgt   = (*Data)[idrs].w[j] * (*DataT)[idt].w[i];
    }

  return(pts);
}

// ========================== CreateBrickPoints ============================

cIntPoint *cIntPoint :: CreateBrickPoints(int order[], eQuadType type, int *n)
{
  int DataSize = GetIntPntDataSize(type);

  for(int i = 0; i < 3; i++)
  {
    if (order[i] < 1 || order[i] > DataSize)
    {
      cout << "\nInvalid cube integration order (" << order << ")!\n";
      if (order[i] > DataSize)
        cout << "\nThe maximum order supported is " << DataSize << endl;
      exit(0);
    }
  }

  const vector<sIntData> *Data = GetIntPntData(type);
  *n = order[0]*order[1]*order[2];
  cIntPoint *pts = new cIntPoint[*n];
  int k = 0;
  int idr = order[0] - 1;
  int ids = order[1] - 1;
  int idt = order[2] - 1;
  for (int i = 0; i < order[2]; i++)
  {
    for (int j = 0; j < order[1]; j++)
    {
      for (int l = 0; l < order[0]; l++)
      {
        pts[k].pnt.r = (*Data)[idr].r[l];
        pts[k].pnt.s = (*Data)[ids].r[j];
        pts[k].pnt.t = (*Data)[idt].r[i];
        pts[k].wgt   = (*Data)[idr].w[l] * (*Data)[ids].w[j] * (*Data)[idt].w[i];
        k++;
      }
    }
  }

  return(pts);
}

// ============================ CreateTetPoints ============================

cIntPoint *cIntPoint :: CreateTetPoints(int order, eQuadType qtype, int *n)
{
//  cout << "cIntPoint :: CreateTetPoints\n";
  int DataSize = GetIntPntDataSize(qtype);

  if (order < 1 || order > DataSize)
  {
    cout << "\nInvalid tetrahedra integration order (" << order << ")!\n";
    if (order > DataSize)
      cout << "\nThe maximum order supported is " << DataSize << endl;
    exit(0);
  }

  const vector<sIntData> *Data = GetIntPntData(qtype);
  *n = (*Data)[order-1].np;
  cIntPoint *pts = new cIntPoint[*n];
  int id = order - 1;

  for (int i = 0; i < *n; i++)
  {
     pts[i].pnt.r = (*Data)[id].r[i];
     pts[i].pnt.s = (*Data)[id].s[i];
     pts[i].pnt.t = (*Data)[id].t[i];
     pts[i].wgt   = (*Data)[id].w[i];
  }

  return(pts);
}

/*
// ======================= CreateLobattoLinePoints =========================

cIntPoint *cIntPoint :: CreateLobattoLinePoints(int order, int *n)
{
  if (order < 1 || order > 4)
  {
    cout << "\nInvalid Lobatto integration order (" << order << ")!\n";
    exit(0);
  }

  *n = order;
  cIntPoint *pts = new cIntPoint[*n];
  int col = order - 1;
  for (int i = 0; i < order; i++)
  {
     pts[i].pnt.r = LobattoPoint[i][col];
     pts[i].wgt   = LobattoWeight[i][col];
  }

  return(pts);
}

// ========================== CreateNCLinePoints ===========================

cIntPoint *cIntPoint :: CreateNCLinePoints(int order, int *n)
{
  if (order < 1 || order > 12)
  {
    cout << "\nInvalid Newton-Cotes integration order (" << order << ")!\n";
    exit(0);
  }

  *n = order;
  cIntPoint *pts = new cIntPoint[*n];
  int col = order - 1;
  for (int i = 0; i < order; i++)
  {
     pts[i].pnt.r = NewtonCotesPoint[i][col];
     pts[i].wgt   = NewtonCotesWeight[i][col];
  }

  return(pts);
}

// ========================== CreateNCQuadPoints ===========================

cIntPoint *cIntPoint :: CreateNCQuadPoints(int order, int *n)
{
  if (order < 1 || order > 12)
  {
    cout << "\nInvalid Newton-Cotes integration order (" << order << ")!\n";
    exit(0);
  }

  *n = order*order;
  cIntPoint *pts = new cIntPoint[*n];
  int k = 0;
  int col = order - 1;
  for (int i = 0; i < order; i++)
  {
    for (int j = 0; j < order; j++)
    {
      pts[k].pnt.r = NewtonCotesPoint[j][col];
      pts[k].pnt.s = NewtonCotesPoint[i][col];
      pts[k].wgt   = NewtonCotesWeight[i][col]*NewtonCotesWeight[j][col];
      k++;
    }
  }

  return(pts);
}
*/
// =============================== cIntPoint ===============================

cIntPoint :: cIntPoint(void)
{
}

// ============================== ~cIntPoint ===============================

cIntPoint :: ~cIntPoint(void)
{
}

// -------------------------------------------------------------------------
// Quadrature data functions:
//

// ============================== getGaussData ===============================
//

const vector<sIntData>& getGaussData( )
{
  static vector<sIntData> data;

  if (data.empty( )) 
  {
    // -------------------------------------------------------------------------
    // Gauss data:
    //
    
    // n = 1.
    
    static const double gaur1[] =
    { 0.000000000000000};
    static const double gauw1[] =
    { 2.000000000000000};
    
    // n = 2.
    
    static const double gaur2[] =
    {-0.577350269189626, 0.577350269189626};
    static const double gauw2[] =
    { 1.000000000000000, 1.000000000000000};
    
    // n = 3.
    
    static const double gaur3[] =
    {-0.774596669241483, 0.000000000000000, 0.774596669241483};
    static const double gauw3[] =
    { 0.555555555555556, 0.888888888888889, 0.555555555555556};
    
    // n = 4.
    
    static const double gaur4[] =
    {-0.861136311594053,-0.339981043584856, 0.339981043584856, 0.861136311594053};
    static const double gauw4[] =
    { 0.347854845137454, 0.652145154862546, 0.652145154862546, 0.347854845137454};
    
    // n = 5.
    
    static const double gaur5[] =
    {-0.906179845938664,-0.538469310105683, 0.000000000000000, 0.538469310105683,
      0.906179845938664};
    static const double gauw5[] =
    { 0.236926885056189, 0.478628670499367, 0.568888888888889, 0.478628670499367,
      0.236926885056189};
    
    // n = 6.
    
    static const double gaur6[] =
    {-0.932469514203152,-0.661209386466264,-0.238619186083197, 0.238619186083197,
      0.661209386466264, 0.932469514203152};
    static const double gauw6[] =
    { 0.171324492379170, 0.360761573048139, 0.467913934572691, 0.467913934572691,
      0.360761573048139, 0.171324492379170};
    
    // n = 7.
    
    static const double gaur7[] =
    {-0.949107912342758,-0.741531185599394,-0.405845151377397, 0.000000000000000,
      0.405845151377397, 0.741531185599394, 0.949107912342758};
    static const double gauw7[] =
    { 0.129484966168870, 0.279705391489277, 0.381830050505119, 0.417959183673469,
      0.381830050505119, 0.279705391489277, 0.129484966168870};
    
    // n = 8.
    
    static const double gaur8[] =
    {-0.960289856497536,-0.796666477413627,-0.525532409916329,-0.183434642495650,
      0.183434642495650, 0.525532409916329, 0.796666477413627, 0.960289856497536};
    static const double gauw8[] =
    { 0.101228536290376, 0.222381034453375, 0.313706645877887, 0.362683783378362,
      0.362683783378362, 0.313706645877887, 0.222381034453375, 0.101228536290376};
    
    // n = 9.
    
    static const double gaur9[] =
    {-0.968160239507626,-0.836031107326636,-0.613371432700590,-0.324253423403809,
      0.000000000000000, 0.324253423403809, 0.613371432700590, 0.836031107326636,
      0.968160239507626};
    static const double gauw9[] =
    { 0.081274388361574, 0.180648160694857, 0.260610696402935, 0.312347077040003,
      0.330239355001260, 0.312347077040003, 0.260610696402935, 0.180648160694857,
      0.081274388361574};
    
    // n = 10.
    
    static const double gaur10[] =
    {-0.973906528517172,-0.865063366688985,-0.679409568299024,-0.433395394129247,
     -0.148874338981631, 0.148874338981631, 0.433395394129247, 0.679409568299024,
      0.865063366688985, 0.973906528517172};
    static const double gauw10[] =
    { 0.066671344308688, 0.149451349150581, 0.219086362515982, 0.269266719309996,
      0.295524224714753, 0.295524224714753, 0.269266719309996, 0.219086362515982,
      0.149451349150581, 0.066671344308688};
    
    // n = 11.
    
    static const double gaur11[] =
    {-0.978228658146057,-0.887062599768095,-0.730152005574049,-0.519096129206812,
     -0.269543155952345, 0.000000000000000, 0.269543155952345, 0.519096129206812,
      0.730152005574049, 0.887062599768095, 0.978228658146057};
    static const double gauw11[] =
    { 0.055668567116174, 0.125580369464905, 0.186290210927734, 0.233193764591991,
      0.262804544510247, 0.272925086777901, 0.262804544510247, 0.233193764591991,
      0.186290210927734, 0.125580369464905, 0.055668567116174};
    
    // n = 12.
    
    static const double gaur12[] =
    {-0.981560634246719,-0.904117256370475,-0.769902674194305,-0.587317954286617,
     -0.367831498998180,-0.125233408511469, 0.125233408511469, 0.367831498998180,
      0.587317954286617, 0.769902674194305, 0.904117256370475, 0.981560634246719};
    static const double gauw12[] =
    { 0.047175336386512, 0.106939325995318, 0.160078328543346, 0.203167426723066,
      0.233492536538355, 0.249147045813403, 0.249147045813403, 0.233492536538355,
      0.203167426723066, 0.160078328543346, 0.106939325995318, 0.047175336386512};
    
    // n = 13.
    
    static const double gaur13[] =
    {-0.984183054718588,-0.917598399222978,-0.801578090733310,-0.642349339440340,
     -0.448492751036447,-0.230458315955135, 0.000000000000000, 0.230458315955135,
      0.448492751036447, 0.642349339440340, 0.801578090733310, 0.917598399222978,
      0.984183054718588};
    static const double gauw13[] =
    { 0.040484004765316, 0.092121499837728, 0.138873510219787, 0.178145980761946,
      0.207816047536889, 0.226283180262897, 0.232551553230874, 0.226283180262897,
      0.207816047536889, 0.178145980761946, 0.138873510219787, 0.092121499837728,
      0.040484004765316};
    
    // n = 14.
    
    static const double gaur14[] =
    {-0.986283808696812,-0.928434883663574,-0.827201315069765,-0.687292904811685,
     -0.515248636358154,-0.319112368927890,-0.108054948707344, 0.108054948707344,
      0.319112368927890, 0.515248636358154, 0.687292904811685, 0.827201315069765,
      0.928434883663574, 0.986283808696812};
    static const double gauw14[] =
    { 0.035119460331752, 0.080158087159760, 0.121518570687903, 0.157203167158193,
      0.185538397477938, 0.205198463721296, 0.215263853463158, 0.215263853463158,
      0.205198463721296, 0.185538397477938, 0.157203167158193, 0.121518570687903,
      0.080158087159760, 0.035119460331752};
    
    // n = 15.
    
    static const double gaur15[] =
    {-0.987992518020485,-0.937273392400706,-0.848206583410427,-0.724417731360170,
     -0.570972172608539,-0.394151347077563,-0.201194093997434, 0.000000000000000,
      0.201194093997434, 0.394151347077563, 0.570972172608539, 0.724417731360170,
      0.848206583410427, 0.937273392400706, 0.987992518020485};
    static const double gauw15[] =
    { 0.030753241996117, 0.070366047488108, 0.107159220467172, 0.139570677926154,
      0.166269205816994, 0.186161000015562, 0.198431485327112, 0.202578241925561,
      0.198431485327112, 0.186161000015562, 0.166269205816994, 0.139570677926154,
      0.107159220467172, 0.070366047488108, 0.030753241996117};
    
    // Array of Gauss data.
    
    data.push_back(sIntData(gaur1 ,0,0, gauw1 ,1 ));
    data.push_back(sIntData(gaur2 ,0,0, gauw2 ,2 ));
    data.push_back(sIntData(gaur3 ,0,0, gauw3 ,3 ));
    data.push_back(sIntData(gaur4 ,0,0, gauw4 ,4 ));
    data.push_back(sIntData(gaur5 ,0,0, gauw5 ,5 ));
    data.push_back(sIntData(gaur6 ,0,0, gauw6 ,6 ));
    data.push_back(sIntData(gaur7 ,0,0, gauw7 ,7 ));
    data.push_back(sIntData(gaur8 ,0,0, gauw8 ,8 ));
    data.push_back(sIntData(gaur9 ,0,0, gauw9 ,9 ));
    data.push_back(sIntData(gaur10,0,0, gauw10,10));
    data.push_back(sIntData(gaur11,0,0, gauw11,11));
    data.push_back(sIntData(gaur12,0,0, gauw12,12));
    data.push_back(sIntData(gaur13,0,0, gauw13,13));
    data.push_back(sIntData(gaur14,0,0, gauw14,14));
    data.push_back(sIntData(gaur15,0,0, gauw15,15));
  }

  return data;
}

// ============================== getLobattoData ===========================
//

const vector<sIntData>& getLobattoData( )
{
  static vector<sIntData> data;

  if (data.empty( )) 
  {
    // -------------------------------------------------------------------------
    // Lobatto data:
    //
    
    // n = 1
    
    static const double lobr1[] =
    {+0.000000000000000};
    static const double lobw1[] =
    {+2.000000000000000};
    
    // n = 2
    
    static const double lobr2[] =
    {-1.000000000000000,+1.000000000000000};
    static const double lobw2[] =
    {+1.000000000000000,+1.000000000000000};
    
    // n = 3
    
    static const double lobr3[] =
    {-1.000000000000000,+0.000000000000000,+1.000000000000000};
    static const double lobw3[] =
    {+0.333333333333333,+1.333333333333333,+0.333333333333333};
    
    // n = 4
    
    static const double lobr4[] =
    {-1.000000000000000,-0.447213595499958,+0.447213595499958,+1.000000000000000};
    static const double lobw4[] =
    {+0.166666666666667,+0.833333333333333,+0.833333333333333,+0.166666666666667};
    
    // n = 5
    
    static const double lobr5[] =
    {-1.000000000000000,-0.654653670707977,+0.000000000000000,+0.654653670707977,
     +1.000000000000000};
    static const double lobw5[] =
    {+0.100000000000000,+0.544444444444444,+0.711111111111111,+0.544444444444444,
     +0.100000000000000};
    
    // n = 6
    
    static const double lobr6[] =
    {-1.000000000000000,-0.765055323929465,-0.285231516480645,+0.285231516480645,
     +0.765055323929465,+1.000000000000000};
    static const double lobw6[] =
    {+0.066666666666667,+0.378474956297847,+0.554858377035486,+0.554858377035486,
     +0.378474956297847,+0.066666666666667};
    
    // n = 7
    
    static const double lobr7[] =
    {-1.000000000000000,-0.830223896278567,-0.468848793470714,+0.000000000000000,
     +0.468848793470714,+0.830223896278567,+1.000000000000000};
    static const double lobw7[] =
    {+0.047619047619048,+0.276826047361566,+0.431745381209863,+0.487619047619048,
     +0.431745381209863,+0.276826047361566,+0.047619047619048};
    
    // n = 8
    
    static const double lobr8[] =
    {-1.000000000000000,-0.871740148509607,-0.591700181433142,-0.209299217902479,
     +0.209299217902479,+0.591700181433142,+0.871740148509607,+1.000000000000000};
    const double lobw8[] =
    {+0.035714285714286,+0.210704227143506,+0.341122692483504,+0.412458794658704,
     +0.412458794658704,+0.341122692483504,+0.210704227143506,+0.035714285714286};
    
    // n = 9
    
    static const double lobr9[] =
    {-1.000000000000000,-0.899757995411460,-0.677186279510738,-0.363117463826178,
     +0.000000000000000,+0.363117463826178,+0.677186279510738,+0.899757995411460,
     +1.000000000000000};
    static const double lobw9[] =
    {+0.027777777777778,+0.165495361560806,+0.274538712500162,+0.346428510973046,
     +0.371519274376417,+0.346428510973046,+0.274538712500162,+0.165495361560806,
     +0.027777777777778};
    
    // n = 10
    
    static const double lobr10[] =
    {-1.000000000000000,-0.919533908166459,-0.738773865105505,-0.477924949810444,
     -0.165278957666387,+0.165278957666387,+0.477924949810444,+0.738773865105505,
     +0.919533908166459,+1.000000000000000};
    static const double lobw10[] =
    {+0.022222222222222,+0.133305990851070,+0.224889342063126,+0.292042683679684,
     +0.327539761183897,+0.327539761183897,+0.292042683679684,+0.224889342063126,
     +0.133305990851070,+0.022222222222222};
    
    // n = 11
    
    static const double lobr11[] =
    {-1.000000000000000,-0.934001430408059,-0.784483473663144,-0.565235326996205,
     -0.295758135586939,+0.000000000000000,+0.295758135586939,+0.565235326996205,
     +0.784483473663144,+0.934001430408059,+1.000000000000000};
    static const double lobw11[] =
    {+0.018181818181818,+0.109612273266995,+0.187169881780305,+0.248048104264028,
     +0.286879124779008,+0.300217595455691,+0.286879124779008,+0.248048104264028,
     +0.187169881780305,+0.109612273266995,+0.018181818181818};
    
    // n = 12
    
    static const double lobr12[] =
    {-1.000000000000000,-0.944899272222882,-0.819279321644007,-0.632876153031870,
     -0.399530940965349,-0.136552932854928,+0.136552932854928,+0.399530940965349,
     +0.632876153031870,+0.819279321644007,+0.944899272222882,+1.000000000000000};
    static const double lobw12[] =
    {+0.015151515151515,+0.091684517413196,+0.157974705564370,+0.212508417761021,
     +0.251275603199201,+0.271405240910696,+0.271405240910696,+0.251275603199201,
     +0.212508417761021,+0.157974705564370,+0.091684517413196,+0.015151515151515};
    
    // n = 13
    
    static const double lobr13[] =
    {-1.000000000000000,-0.953309846642164,-0.846347564651872,-0.686188469081757,
     -0.482909821091336,-0.249286930106240,+0.000000000000000,+0.249286930106240,
     +0.482909821091336,+0.686188469081757,+0.846347564651872,+0.953309846642164,
     +1.000000000000000};
    static const double lobw13[] =
    {+0.012820512820513,+0.077801686746819,+0.134981926689608,+0.183646865203550,
     +0.220767793566110,+0.244015790306676,+0.251930849333447,+0.244015790306676,
     +0.220767793566110,+0.183646865203550,+0.134981926689608,+0.077801686746819,
     +0.012820512820513};
    
    // n = 14
    
    static const double lobr14[] =
    {-1.000000000000000,-0.959935045267261,-0.867801053830347,-0.728868599091326,
     -0.550639402928647,-0.342724013342713,-0.116331868883704,+0.116331868883704,
     +0.342724013342713,+0.550639402928647,+0.728868599091326,+0.867801053830347,
     +0.959935045267261,+1.000000000000000};
    static const double lobw14[] =
    {+0.010989010989011,+0.066837284497681,+0.116586655898712,+0.160021851762952,
     +0.194826149373416,+0.219126253009771,+0.231612794468457,+0.231612794468457,
     +0.219126253009771,+0.194826149373416,+0.160021851762952,+0.116586655898712,
     +0.066837284497681,+0.010989010989011};
    
    // n = 15
    
    static const double lobr15[] =
    {-1.000000000000000,-0.965245926503839,-0.885082044222976,-0.763519689951815,
     -0.606253205469846,-0.420638054713672,-0.215353955363794,+0.000000000000000,
     +0.215353955363794,+0.420638054713672,+0.606253205469846,+0.763519689951815,
     +0.885082044222976,+0.965245926503839,+1.000000000000000};
    static const double lobw15[] =
    {+0.009523809523810,+0.058029893028601,+0.101660070325718,+0.140511699802428,
     +0.172789647253601,+0.196987235964613,+0.211973585926821,+0.217048116348816,
     +0.211973585926821,+0.196987235964613,+0.172789647253601,+0.140511699802428,
     +0.101660070325718,+0.058029893028601,+0.009523809523810};
    
    // n = 16
    
    static const double lobr16[] =
    {-1.000000000000000,-0.969568046270218,-0.899200533093472,-0.792008291861815,
     -0.652388702882493,-0.486059421887138,-0.299830468900763,-0.101326273521949,
     +0.101326273521949,+0.299830468900763,+0.486059421887138,+0.652388702882493,
     +0.792008291861815,+0.899200533093472,+0.969568046270218,+1.000000000000000};
    static const double lobw16[] =
    {+0.008333333333333,+0.050850361005920,+0.089393697325931,+0.124255382132514,
     +0.154026980807164,+0.177491913391704,+0.193690023825204,+0.201958308178230,
     +0.201958308178230,+0.193690023825204,+0.177491913391704,+0.154026980807164,
     +0.124255382132514,+0.089393697325931,+0.050850361005920,+0.008333333333333};
    
    // n = 17
    
    static const double lobr17[] =
    {-1.000000000000000,-0.973132176631418,-0.910879995915574,-0.815696251221770,
     -0.691028980627685,-0.541385399330102,-0.372174433565477,-0.189511973518317,
     +0.000000000000000,+0.189511973518317,+0.372174433565477,+0.541385399330102,
     +0.691028980627685,+0.815696251221770,+0.910879995915574,+0.973132176631418,
     +1.000000000000000};
    static const double lobw17[] =
    {+0.007352941176471,+0.044921940543254,+0.079198270503687,+0.110592909007028,
     +0.137987746201927,+0.160394661997622,+0.177004253515658,+0.187216339677619,
     +0.190661874753469,+0.187216339677619,+0.177004253515658,+0.160394661997622,
     +0.137987746201927,+0.110592909007028,+0.079198270503687,+0.044921940543254,
     +0.007352941176471};
    
    // n = 18
    
    static const double lobr18[] =
    {-1.000000000000000,-0.976105557412199,-0.920649185347534,-0.835593535218090,
     -0.723679329283243,-0.588504834318662,-0.434415036912124,-0.266362652878281,
     -0.089749093484652,+0.089749093484652,+0.266362652878281,+0.434415036912124,
     +0.588504834318662,+0.723679329283243,+0.835593535218090,+0.920649185347534,
     +0.976105557412199,+1.000000000000000};
    static const double lobw18[] =
    {+0.006535947712418,+0.039970628810914,+0.070637166885634,+0.099016271717503,
     +0.124210533132967,+0.145411961573802,+0.161939517237602,+0.173262109489456,
     +0.179015863439703,+0.179015863439703,+0.173262109489456,+0.161939517237602,
     +0.145411961573802,+0.124210533132967,+0.099016271717503,+0.070637166885634,
     +0.039970628810914,+0.006535947712418};
    
    // n = 19
    
    static const double lobr19[] =
    {-1.000000000000000,-0.978611766222080,-0.928901528152586,-0.852460577796646,
     -0.751494202552613,-0.628908137265221,-0.488229285680713,-0.333504847824499,
     -0.169186023409282,+0.000000000000000,+0.169186023409282,+0.333504847824499,
     +0.488229285680713,+0.628908137265221,+0.751494202552613,+0.852460577796646,
     +0.928901528152586,+0.978611766222080,+1.000000000000000};
    static const double lobw19[] =
    {+0.005847953216374,+0.035793365186176,+0.063381891762630,+0.089131757099207,
     +0.112315341477305,+0.132267280448751,+0.148413942595939,+0.160290924044061,
     +0.167556584527143,+0.170001919284827,+0.167556584527143,+0.160290924044061,
     +0.148413942595939,+0.132267280448751,+0.112315341477305,+0.089131757099207,
     +0.063381891762630,+0.035793365186176,+0.005847953216374};
    
    // n = 20
    
    static const double lobr20[] =
    {-1.000000000000000,-0.980743704893914,-0.935934498812665,-0.866877978089950,
     -0.775368260952056,-0.663776402290311,-0.534992864031886,-0.392353183713909,
     -0.239551705922986,-0.080545937238822,+0.080545937238822,+0.239551705922986,
     +0.392353183713909,+0.534992864031886,+0.663776402290311,+0.775368260952056,
     +0.866877978089950,+0.935934498812665,+0.980743704893914,+1.000000000000000};
    static const double lobw20[] =
    {+0.005263157894737,+0.032237123188489,+0.057181802127567,+0.080631763996120,
     +0.101991499699451,+0.120709227628675,+0.136300482358724,+0.148361554070917,
     +0.156580102647475,+0.160743286387846,+0.160743286387846,+0.156580102647475,
     +0.148361554070917,+0.136300482358724,+0.120709227628675,+0.101991499699451,
     +0.080631763996120,+0.057181802127567,+0.032237123188489,+0.005263157894737};
    
    // Array of Lobatto data.
    
    data.push_back(sIntData(lobr1 ,0,0, lobw1 ,1 ));
    data.push_back(sIntData(lobr2 ,0,0, lobw2 ,2 ));
    data.push_back(sIntData(lobr3 ,0,0, lobw3 ,3 ));
    data.push_back(sIntData(lobr4 ,0,0, lobw4 ,4 ));
    data.push_back(sIntData(lobr5 ,0,0, lobw5 ,5 ));
    data.push_back(sIntData(lobr6 ,0,0, lobw6 ,6 ));
    data.push_back(sIntData(lobr7 ,0,0, lobw7 ,7 ));
    data.push_back(sIntData(lobr8 ,0,0, lobw8 ,8 ));
    data.push_back(sIntData(lobr9 ,0,0, lobw9 ,9 ));
    data.push_back(sIntData(lobr10,0,0, lobw10,10));
    data.push_back(sIntData(lobr11,0,0, lobw11,11));
    data.push_back(sIntData(lobr12,0,0, lobw12,12));
    data.push_back(sIntData(lobr13,0,0, lobw13,13));
    data.push_back(sIntData(lobr14,0,0, lobw14,14));
    data.push_back(sIntData(lobr15,0,0, lobw15,15));
    data.push_back(sIntData(lobr16,0,0, lobw16,16));
    data.push_back(sIntData(lobr17,0,0, lobw17,17));
    data.push_back(sIntData(lobr18,0,0, lobw18,18));
    data.push_back(sIntData(lobr19,0,0, lobw19,19));
    data.push_back(sIntData(lobr20,0,0, lobw20,20));
  }

  return data;
}

// ============================== getNewtonCotesData =======================

const vector<sIntData>& getNewtonCotesData( )
{
  static vector<sIntData> data;

  if (data.empty( )) 
  {
    // -------------------------------------------------------------------------
    // Newton-Cotes data:
    //
    
    // n = 1
    
    static const double newcotr1[] =
    {+0.000000000000000};
    static const double newcotw1[] =
    {+2.000000000000000};
    
    // n = 2
    
    static const double newcotr2[] =
    {-1.000000000000000,+1.000000000000000};
    static const double newcotw2[] =
    {+1.000000000000000,+1.000000000000000};
    
    // n = 3
    
    static const double newcotr3[] =
    {-1.000000000000000,+0.000000000000000,+1.000000000000000};
    static const double newcotw3[] =
    {+0.333333333333333,+1.333333333333333,+0.333333333333333};
    
    // n = 4
    
    static const double newcotr4[] =
    {-1.000000000000000,-0.333333333333333,+0.333333333333333,+1.000000000000000};
    static const double newcotw4[] =
    {+0.250000000000000,+0.750000000000000,+0.750000000000000,+0.250000000000000};
    
    // n = 5
    
    static const double newcotr5[] =
    {-1.000000000000000,-0.500000000000000,+0.000000000000000,+0.500000000000000,
     +1.000000000000000};
    static const double newcotw5[] =
    {+0.155555555555556,+0.711111111111111,+0.266666666666667,+0.711111111111111,
     +0.155555555555556};
    
    // n = 6
    
    static const double newcotr6[] =
    {-1.000000000000000,-0.600000000000000,-0.200000000000000,+0.200000000000000,
     +0.600000000000000,+1.000000000000000};
    static const double newcotw6[] =
    {+0.131944444444444,+0.520833333333333,+0.347222222222222,+0.347222222222222,
     +0.520833333333333,+0.131944444444444};
    
    // n = 7
    
    static const double newcotr7[] =
    {-1.000000000000000,-0.666666666666667,-0.333333333333333,+0.000000000000000,
     +0.333333333333333,+0.666666666666667,+1.000000000000000};
    static const double newcotw7[] =
    {+0.097619047619048,+0.514285714285714,+0.064285714285714,+0.647619047619048,
     +0.064285714285714,+0.514285714285714,+0.097619047619048};
    
    // n = 8
    
    static const double newcotr8[] =
    {-1.000000000000000,-0.714285714285714,-0.428571428571429,-0.142857142857143,
     +0.142857142857143,+0.428571428571429,+0.714285714285714,+1.000000000000000};
    static const double newcotw8[] =
    {+0.086921296296296,+0.414004629629630,+0.153125000000000,+0.345949074074074,
     +0.345949074074074,+0.153125000000000,+0.414004629629630,+0.086921296296296};
    
    // n = 9
    
    static const double newcotr9[] =
    {-1.000000000000000,-0.750000000000000,-0.500000000000000,-0.250000000000000,
     +0.000000000000000,+0.250000000000000,+0.500000000000000,+0.750000000000000,
     +1.000000000000000};
    static const double newcotw9[] =
    {+0.069770723104056,+0.415379188712522,-0.065467372134039,+0.740458553791887,
     -0.320282186948854,+0.740458553791887,-0.065467372134039,+0.415379188712522,
     +0.069770723104056};
    
    // n = 10
    
    static const double newcotr10[] =
    {-1.000000000000000,-0.777777777777778,-0.555555555555556,-0.333333333333333,
     -0.111111111111111,+0.111111111111111,+0.333333333333333,+0.555555555555556,
     +0.777777777777778,+1.000000000000000};
    static const double newcotw10[] =
    {+0.063772321428571,+0.351361607142857,+0.024107142857143,+0.431785714285714,
     +0.128973214285714,+0.128973214285714,+0.431785714285714,+0.024107142857143,
     +0.351361607142857,+0.063772321428571};
    
    // n = 11
    
    static const double newcotr11[] =
    {-1.000000000000000,-0.800000000000000,-0.600000000000000,-0.400000000000000,
     -0.200000000000000,+0.000000000000000,+0.200000000000000,+0.400000000000000,
     +0.600000000000000,+0.800000000000000,+1.000000000000000};
    static const double newcotw11[] =
    {+0.053668296723852,+0.355071882849661,-0.162087141253808,+0.909892576559243,
     -0.870310245310245,+1.427529260862594,-0.870310245310245,+0.909892576559243,
     -0.162087141253808,+0.355071882849661,+0.053668296723852};
    
    // n = 12
    
    static const double newcotr12[] =
    {-1.000000000000000,-0.818181818181818,-0.636363636363636,-0.454545454545455,
     -0.272727272727273,-0.090909090909091,+0.090909090909091,+0.272727272727273,
     +0.454545454545455,+0.636363636363636,+0.818181818181818,+1.000000000000000};
    static const double newcotw12[] =
    {+0.049866461823927,+0.309710717041446,-0.074338463587596,+0.579316509589947,
     -0.220356178350970,+0.355800953483245,+0.355800953483245,-0.220356178350970,
     +0.579316509589947,-0.074338463587596,+0.309710717041446,+0.049866461823927};
    
    // n = 13
    
    static const double newcotr13[] =
    {-1.000000000000000,-0.833333333333333,-0.666666666666667,-0.500000000000000,
     -0.333333333333333,-0.166666666666667,+0.000000000000000,+0.166666666666667,
     +0.333333333333333,+0.500000000000000,+0.666666666666667,+0.833333333333333,
     +1.000000000000000};
    static const double newcotw13[] =
    {+0.043278974993261,+0.314072213500785,-0.240643927501070,+1.132997795854939,
     -1.633011274439846,+2.775519337805052,-2.784426240426241,+2.775519337805052,
     -1.633011274439846,+1.132997795854939,-0.240643927501070,+0.314072213500785,
     +0.043278974993261};
    
    // n = 14
    
    static const double newcotr14[] =
    {-1.000000000000000,-0.846153846153846,-0.692307692307692,-0.538461538461538,
     -0.384615384615385,-0.230769230769231,-0.076923076923077,+0.076923076923077,
     +0.230769230769231,+0.384615384615385,+0.538461538461538,+0.692307692307692,
     +0.846153846153846,+1.000000000000000};
    static const double newcotw14[] =
    {+0.040669438210247,+0.279752170531571,-0.155423740576828,+0.775792308487766,
     -0.753847632664235,+1.027352359112311,-0.214294903100831,-0.214294903100831,
     +1.027352359112311,-0.753847632664235,+0.775792308487766,-0.155423740576828,
     +0.279752170531571,+0.040669438210247};
    
    // n = 15
    
    static const double newcotr15[] =
    {-1.000000000000000,-0.857142857142857,-0.714285714285714,-0.571428571428571,
     -0.428571428571429,-0.285714285714286,-0.142857142857143,+0.000000000000000,
     +0.142857142857143,+0.285714285714286,+0.428571428571429,+0.571428571428571,
     +0.714285714285714,+0.857142857142857,+1.000000000000000};
    static const double newcotw15[] =
    {+0.036068942431597,+0.284175589385466,-0.308050694104706,+1.399497820880537,
     -2.647995211293051,+5.048155508871559,-6.715728979011386,+7.807754045679972,
     -6.715728979011386,+5.048155508871559,-2.647995211293051,+1.399497820880537,
     -0.308050694104706,+0.284175589385466,+0.036068942431597};
    
    // n = 16
    
    static const double newcotr16[] =
    {-1.000000000000000,-0.866666666666667,-0.733333333333333,-0.600000000000000,
     -0.466666666666667,-0.333333333333333,-0.200000000000000,-0.066666666666667,
     +0.066666666666667,+0.200000000000000,+0.333333333333333,+0.466666666666667,
     +0.600000000000000,+0.733333333333333,+0.866666666666667,+1.000000000000000};
    static const double newcotw16[] =
    {+0.034174599543252,+0.257014757354810,-0.225445810119010,+1.014085416420447,
     -1.512586229688280,+2.382720699013598,-1.936010422992496,+0.986046990467680,
     +0.986046990467680,-1.936010422992496,+2.382720699013598,-1.512586229688280,
     +1.014085416420447,-0.225445810119010,+0.257014757354810,+0.034174599543252};
    
    // n = 17
    
    static const double newcotr17[] =
    {-1.000000000000000,-0.875000000000000,-0.750000000000000,-0.625000000000000,
     -0.500000000000000,-0.375000000000000,-0.250000000000000,-0.125000000000000,
     +0.000000000000000,+0.125000000000000,+0.250000000000000,+0.375000000000000,
     +0.500000000000000,+0.625000000000000,+0.750000000000000,+0.875000000000000,
     +1.000000000000000};
    static const double newcotw17[] =
    {+0.030797894233299,+0.261282382880280,-0.367952893298676,+1.703737977809009,
     -3.950148071778393,+8.552529993440295,-13.934614237197881,+19.180342211078734,
     -20.951950514333333,+19.180342211078734,-13.934614237197881,+8.552529993440295,
     -3.950148071778393,+1.703737977809009,-0.367952893298676,+0.261282382880280,
     +0.030797894233299};
    
    // n = 18
    
    static const double newcotr18[] =
    {-1.000000000000000,-0.882352941176471,-0.764705882352941,-0.647058823529412,
     -0.529411764705882,-0.411764705882353,-0.294117647058824,-0.176470588235294,
     -0.058823529411765,+0.058823529411765,+0.176470588235294,+0.294117647058824,
     +0.411764705882353,+0.529411764705882,+0.647058823529412,+0.764705882352941,
     +0.882352941176471,+1.000000000000000};
    static const double newcotw18[] =
    {+0.029364429446790,+0.239072385160517,-0.287843192311834,+1.289734802610926,
     -2.532547749581263,+4.702695904581749,-5.791330845017044,+5.347500024845654,
     -1.996645759735495,-1.996645759735495,+5.347500024845654,-5.791330845017044,
     +4.702695904581749,-2.532547749581263,+1.289734802610926,-0.287843192311834,
     +0.239072385160517,+0.029364429446790};
    
    // n = 19
    
    static const double newcotr19[] =
    {-1.000000000000000,-0.888888888888889,-0.777777777777778,-0.666666666666667,
     -0.555555555555556,-0.444444444444444,-0.333333333333333,-0.222222222222222,
     -0.111111111111111,+0.000000000000000,+0.111111111111111,+0.222222222222222,
     +0.333333333333333,+0.444444444444444,+0.555555555555556,+0.666666666666667,
     +0.777777777777778,+0.888888888888889,+1.000000000000000};
    static const double newcotw19[] =
    {+0.026790824664820,+0.243108208883743,-0.422476206213465,+2.042174237602923,
     -5.571479168174973,+13.634004454324977,-26.122288374274994,+41.953237533490707,
     -55.115367445968609,+60.664591871329741,-55.115367445968609,+41.953237533490707,
     -26.122288374274994,+13.634004454324977,-5.571479168174973,+2.042174237602923,
     -0.422476206213465,+0.243108208883743,+0.026790824664820};
    
    // n = 20
    
    static const double newcotr20[] =
    {-1.000000000000000,-0.894736842105263,-0.789473684210526,-0.684210526315789,
     -0.578947368421053,-0.473684210526316,-0.368421052631579,-0.263157894736842,
     -0.157894736842105,-0.052631578947368,+0.052631578947368,+0.157894736842105,
     +0.263157894736842,+0.368421052631579,+0.473684210526316,+0.578947368421053,
     +0.684210526315789,+0.789473684210526,+0.894736842105263,+1.000000000000000};
    static const double newcotw20[] =
    {+0.025670822345560,+0.224489685952519,-0.344678900990309,+1.599697436697807,
     -3.846673091095298,+8.306599334472983,-13.139430424771119,+16.333513604742677,
     -13.792641220001199,+5.633452752646377,+5.633452752646377,-13.792641220001199,
     +16.333513604742677,-13.139430424771119,+8.306599334472983,-3.846673091095298,
     +1.599697436697807,-0.344678900990309,+0.224489685952519,+0.025670822345560};
    
    // n = 21
    
    static const double newcotr21[] =
    {-1.000000000000000,-0.900000000000000,-0.800000000000000,-0.700000000000000,
     -0.600000000000000,-0.500000000000000,-0.400000000000000,-0.300000000000000,
     -0.200000000000000,-0.100000000000000,+0.000000000000000,+0.100000000000000,
     +0.200000000000000,+0.300000000000000,+0.400000000000000,+0.500000000000000,
     +0.600000000000000,+0.700000000000000,+0.800000000000000,+0.900000000000000,
     +1.000000000000000};
    static const double newcotw21[] =
    {+0.023650546498063,+0.228275435289214,-0.472956741022854,+2.412373786963751,
     -7.542063453430661,+20.673596439879603,-45.417631687959023,+83.656114844387105,
     -128.150558980308006,+165.594566944945711,-180.010734270485784,+165.594566944945711,
     -128.150558980308006,+83.656114844387105,-45.417631687959023,+20.673596439879603,
     -7.542063453430661,+2.412373786963751,-0.472956741022854,+0.228275435289214,
     +0.023650546498063};
    
    // Array of Newton-Cotes data.
    
    data.push_back(sIntData(newcotr1 ,0,0, newcotw1 ,1 ));
    data.push_back(sIntData(newcotr2 ,0,0, newcotw2 ,2 ));
    data.push_back(sIntData(newcotr3 ,0,0, newcotw3 ,3 ));
    data.push_back(sIntData(newcotr4 ,0,0, newcotw4 ,4 ));
    data.push_back(sIntData(newcotr5 ,0,0, newcotw5 ,5 ));
    data.push_back(sIntData(newcotr6 ,0,0, newcotw6 ,6 ));
    data.push_back(sIntData(newcotr7 ,0,0, newcotw7 ,7 ));
    data.push_back(sIntData(newcotr8 ,0,0, newcotw8 ,8 ));
    data.push_back(sIntData(newcotr9 ,0,0, newcotw9 ,9 ));
    data.push_back(sIntData(newcotr10,0,0, newcotw10,10));
    data.push_back(sIntData(newcotr11,0,0, newcotw11,11));
    data.push_back(sIntData(newcotr12,0,0, newcotw12,12));
    data.push_back(sIntData(newcotr13,0,0, newcotw13,13));
    data.push_back(sIntData(newcotr14,0,0, newcotw14,14));
    data.push_back(sIntData(newcotr15,0,0, newcotw15,15));
    data.push_back(sIntData(newcotr16,0,0, newcotw16,16));
    data.push_back(sIntData(newcotr17,0,0, newcotw17,17));
    data.push_back(sIntData(newcotr18,0,0, newcotw18,18));
    data.push_back(sIntData(newcotr19,0,0, newcotw19,19));
    data.push_back(sIntData(newcotr20,0,0, newcotw20,20));
    data.push_back(sIntData(newcotr21,0,0, newcotw21,21));
  }

  return data;
}

// ============================== getDunavantData ==========================

const vector<sIntData>& getDunavantData( )
{
  static vector<sIntData> data;

  if (data.empty( )) 
  {
    // -------------------------------------------------------------------------
    // Dunavant data:
    //
 
    // n = 1
    
    static const double trianr1[] =
    {+0.333333333333333};
    static const double trians1[] =
    {+0.333333333333333};
    static const double trianw1[] =
    {+0.500000000000000};
    
    // n = 2
    
    static const double trianr2[] =
    {+0.666666666666667,+0.166666666666667,+0.166666666666667};
    static const double trians2[] =
    {+0.166666666666667,+0.166666666666667,+0.666666666666667};
    static const double trianw2[] =
    {+0.166666666666666,+0.166666666666666,+0.166666666666666};
    
    // n = 3
    
    static const double trianr3[] =
    {+0.333333333333333,+0.600000000000000,+0.200000000000000,+0.200000000000000};
    static const double trians3[] =
    {+0.333333333333333,+0.200000000000000,+0.200000000000000,+0.600000000000000};
    static const double trianw3[] =
    {-0.281250000000000,+0.260416666666667,+0.260416666666667,+0.260416666666667};
    
    // n = 4
    
    static const double trianr4[] =
    {+0.108103018168070,+0.445948490915965,+0.445948490915965,+0.816847572980459,
     +0.091576213509771,+0.091576213509771};
    static const double trians4[] =
    {+0.445948490915965,+0.445948490915965,+0.108103018168070,+0.091576213509771,
     +0.091576213509771,+0.816847572980459};
    static const double trianw4[] =
    {+0.111690794839005,+0.111690794839005,+0.111690794839005,+0.054975871827661,
     +0.054975871827661,+0.054975871827661};
    
    // n = 5
    
    static const double trianr5[] =
    {+0.333333333333333,+0.059715871789770,+0.470142064105115,+0.470142064105115,
     +0.797426985353087,+0.101286507323456,+0.101286507323456};
    static const double trians5[] =
    {+0.333333333333333,+0.470142064105115,+0.470142064105115,+0.059715871789770,
     +0.101286507323456,+0.101286507323456,+0.797426985353087};
    static const double trianw5[] =
    {+0.112500000000000,+0.066197076394253,+0.066197076394253,+0.066197076394253,
     +0.062969590272414,+0.062969590272414,+0.062969590272414};
    
    // n = 6
    
    static const double trianr6[] =
    {+0.501426509658179,+0.249286745170910,+0.249286745170910,+0.873821971016996,
     +0.063089014491502,+0.063089014491502,+0.053145049844817,+0.310352451033784,
     +0.636502499121399,+0.310352451033784,+0.636502499121399,+0.053145049844817};
    static const double trians6[] =
    {+0.249286745170910,+0.249286745170910,+0.501426509658179,+0.063089014491502,
     +0.063089014491502,+0.873821971016996,+0.310352451033784,+0.636502499121399,
     +0.053145049844817,+0.053145049844817,+0.310352451033784,+0.636502499121399};
    static const double trianw6[] =
    {+0.058393137863189,+0.058393137863189,+0.058393137863189,+0.025422453185103,
     +0.025422453185103,+0.025422453185103,+0.041425537809187,+0.041425537809187,
     +0.041425537809187,+0.041425537809187,+0.041425537809187,+0.041425537809187};
    
    // n = 7
    
    static const double trianr7[] =
    {+0.333333333333333,+0.479308067841920,+0.260345966079040,+0.260345966079040,
     +0.869739794195568,+0.065130102902216,+0.065130102902216,+0.048690315425316,
     +0.312865496004874,+0.638444188569810,+0.312865496004874,+0.638444188569810,
     +0.048690315425316};
    static const double trians7[] =
    {+0.333333333333333,+0.260345966079040,+0.260345966079040,+0.479308067841920,
     +0.065130102902216,+0.065130102902216,+0.869739794195568,+0.312865496004874,
     +0.638444188569810,+0.048690315425316,+0.048690315425316,+0.312865496004874,
     +0.638444188569810};
    static const double trianw7[] =
    {-0.074785022233841,+0.087807628716604,+0.087807628716604,+0.087807628716604,
     +0.026673617804419,+0.026673617804419,+0.026673617804419,+0.038556880445128,
     +0.038556880445128,+0.038556880445128,+0.038556880445128,+0.038556880445128,
     +0.038556880445128};
    
    // n = 8
    
    static const double trianr8[] =
    {+0.333333333333333,+0.081414823414554,+0.459292588292723,+0.459292588292723,
     +0.658861384496480,+0.170569307751760,+0.170569307751760,+0.898905543365938,
     +0.050547228317031,+0.050547228317031,+0.008394777409958,+0.263112829634638,
     +0.728492392955404,+0.263112829634638,+0.728492392955404,+0.008394777409958};
    static const double trians8[] =
    {+0.333333333333333,+0.459292588292723,+0.459292588292723,+0.081414823414554,
     +0.170569307751760,+0.170569307751760,+0.658861384496480,+0.050547228317031,
     +0.050547228317031,+0.898905543365938,+0.263112829634638,+0.728492392955404,
     +0.008394777409958,+0.008394777409958,+0.263112829634638,+0.728492392955404};
    static const double trianw8[] =
    {+0.072157803838894,+0.047545817133642,+0.047545817133642,+0.047545817133642,
     +0.051608685267359,+0.051608685267359,+0.051608685267359,+0.016229248811599,
     +0.016229248811599,+0.016229248811599,+0.013615157087217,+0.013615157087217,
     +0.013615157087217,+0.013615157087217,+0.013615157087217,+0.013615157087217};
    
    // n = 9
    
    static const double trianr9[] =
    {+0.333333333333333,+0.020634961602525,+0.489682519198738,+0.489682519198738,
     +0.125820817014127,+0.437089591492937,+0.437089591492937,+0.623592928761935,
     +0.188203535619033,+0.188203535619033,+0.910540973211095,+0.044729513394453,
     +0.044729513394453,+0.036838412054736,+0.221962989160766,+0.741198598784498,
     +0.221962989160766,+0.741198598784498,+0.036838412054736};
    static const double trians9[] =
    {+0.333333333333333,+0.489682519198738,+0.489682519198738,+0.020634961602525,
     +0.437089591492937,+0.437089591492937,+0.125820817014127,+0.188203535619033,
     +0.188203535619033,+0.623592928761935,+0.044729513394453,+0.044729513394453,
     +0.910540973211095,+0.221962989160766,+0.741198598784498,+0.036838412054736,
     +0.036838412054736,+0.221962989160766,+0.741198598784498};
    static const double trianw9[] =
    {+0.048567898141400,+0.015667350113570,+0.015667350113570,+0.015667350113570,
     +0.038913770502387,+0.038913770502387,+0.038913770502387,+0.039823869463605,
     +0.039823869463605,+0.039823869463605,+0.012788837829349,+0.012788837829349,
     +0.012788837829349,+0.021641769688645,+0.021641769688645,+0.021641769688645,
     +0.021641769688645,+0.021641769688645,+0.021641769688645};
    
    // n = 10
    
    static const double trianr10[] =
    {+0.333333333333333,+0.028844733232685,+0.485577633383657,+0.485577633383657,
     +0.781036849029926,+0.109481575485037,+0.109481575485037,+0.141707219414880,
     +0.307939838764121,+0.550352941820999,+0.307939838764121,+0.550352941820999,
     +0.141707219414880,+0.025003534762686,+0.246672560639903,+0.728323904597411,
     +0.246672560639903,+0.728323904597411,+0.025003534762686,+0.009540815400299,
     +0.066803251012200,+0.923655933587500,+0.066803251012200,+0.923655933587500,
     +0.009540815400299};
    static const double trians10[] =
    {+0.333333333333333,+0.485577633383657,+0.485577633383657,+0.028844733232685,
     +0.109481575485037,+0.109481575485037,+0.781036849029926,+0.307939838764121,
     +0.550352941820999,+0.141707219414880,+0.141707219414880,+0.307939838764121,
     +0.550352941820999,+0.246672560639903,+0.728323904597411,+0.025003534762686,
     +0.025003534762686,+0.246672560639903,+0.728323904597411,+0.066803251012200,
     +0.923655933587500,+0.009540815400299,+0.009540815400299,+0.066803251012200,
     +0.923655933587500};
    static const double trianw10[] =
    {+0.045408995191377,+0.018362978878233,+0.018362978878233,+0.018362978878233,
     +0.022660529717764,+0.022660529717764,+0.022660529717764,+0.036378958422710,
     +0.036378958422710,+0.036378958422710,+0.036378958422710,+0.036378958422710,
     +0.036378958422710,+0.014163621265528,+0.014163621265528,+0.014163621265528,
     +0.014163621265528,+0.014163621265528,+0.014163621265528,+0.004710833481867,
     +0.004710833481867,+0.004710833481867,+0.004710833481867,+0.004710833481867,
     +0.004710833481867};
    
    // n = 11
    
    static const double trianr11[] =
    {-0.069222096541517,+0.534611048270758,+0.534611048270758,+0.202061394068290,
     +0.398969302965855,+0.398969302965855,+0.593380199137435,+0.203309900431282,
     +0.203309900431282,+0.761298175434837,+0.119350912282581,+0.119350912282581,
     +0.935270103777448,+0.032364948111276,+0.032364948111276,+0.050178138310495,
     +0.356620648261293,+0.593201213428213,+0.356620648261293,+0.593201213428213,
     +0.050178138310495,+0.021022016536166,+0.171488980304042,+0.807489003159792,
     +0.171488980304042,+0.807489003159792,+0.021022016536166};
    static const double trians11[] =
    {+0.534611048270758,+0.534611048270758,-0.069222096541517,+0.398969302965855,
     +0.398969302965855,+0.202061394068290,+0.203309900431282,+0.203309900431282,
     +0.593380199137435,+0.119350912282581,+0.119350912282581,+0.761298175434837,
     +0.032364948111276,+0.032364948111276,+0.935270103777448,+0.356620648261293,
     +0.593201213428213,+0.050178138310495,+0.050178138310495,+0.356620648261293,
     +0.593201213428213,+0.171488980304042,+0.807489003159792,+0.021022016536166,
     +0.021022016536166,+0.171488980304042,+0.807489003159792};
    static const double trianw11[] =
    {+0.000463503164481,+0.000463503164481,+0.000463503164481,+0.038574767457406,
     +0.038574767457406,+0.038574767457406,+0.029661488690387,+0.029661488690387,
     +0.029661488690387,+0.018092270251709,+0.018092270251709,+0.018092270251709,
     +0.006829865501339,+0.006829865501339,+0.006829865501339,+0.026168555981102,
     +0.026168555981102,+0.026168555981102,+0.026168555981102,+0.026168555981102,
     +0.026168555981102,+0.010353829819571,+0.010353829819571,+0.010353829819571,
     +0.010353829819571,+0.010353829819571,+0.010353829819571};
    
    // n = 12
    
    static const double trianr12[] =
    {+0.023565220452390,+0.488217389773805,+0.488217389773805,+0.120551215411079,
     +0.439724392294460,+0.439724392294460,+0.457579229975768,+0.271210385012116,
     +0.271210385012116,+0.744847708916828,+0.127576145541586,+0.127576145541586,
     +0.957365299093579,+0.021317350453210,+0.021317350453210,+0.115343494534698,
     +0.275713269685514,+0.608943235779788,+0.275713269685514,+0.608943235779788,
     +0.115343494534698,+0.022838332222257,+0.281325580989940,+0.695836086787803,
     +0.281325580989940,+0.695836086787803,+0.022838332222257,+0.025734050548330,
     +0.116251915907597,+0.858014033544073,+0.116251915907597,+0.858014033544073,
     +0.025734050548330};
    static const double trians12[] =
    {+0.488217389773805,+0.488217389773805,+0.023565220452390,+0.439724392294460,
     +0.439724392294460,+0.120551215411079,+0.271210385012116,+0.271210385012116,
     +0.457579229975768,+0.127576145541586,+0.127576145541586,+0.744847708916828,
     +0.021317350453210,+0.021317350453210,+0.957365299093579,+0.275713269685514,
     +0.608943235779788,+0.115343494534698,+0.115343494534698,+0.275713269685514,
     +0.608943235779788,+0.281325580989940,+0.695836086787803,+0.022838332222257,
     +0.022838332222257,+0.281325580989940,+0.695836086787803,+0.116251915907597,
     +0.858014033544073,+0.025734050548330,+0.025734050548330,+0.116251915907597,
     +0.858014033544073};
    static const double trianw12[] =
    {+0.012865533220227,+0.012865533220227,+0.012865533220227,+0.021846272269019,
     +0.021846272269019,+0.021846272269019,+0.031429112108943,+0.031429112108943,
     +0.031429112108943,+0.017398056465355,+0.017398056465355,+0.017398056465355,
     +0.003083130525780,+0.003083130525780,+0.003083130525780,+0.020185778883191,
     +0.020185778883191,+0.020185778883191,+0.020185778883191,+0.020185778883191,
     +0.020185778883191,+0.011178386601152,+0.011178386601152,+0.011178386601152,
     +0.011178386601152,+0.011178386601152,+0.011178386601152,+0.008658115554329,
     +0.008658115554329,+0.008658115554329,+0.008658115554329,+0.008658115554329,
     +0.008658115554329};
    
    // n = 13
    
    static const double trianr13[] =
    {+0.333333333333333,+0.009903630120591,+0.495048184939705,+0.495048184939705,
     +0.062566729780852,+0.468716635109574,+0.468716635109574,+0.170957326397447,
     +0.414521336801277,+0.414521336801277,+0.541200855914337,+0.229399572042831,
     +0.229399572042831,+0.771151009607340,+0.114424495196330,+0.114424495196330,
     +0.950377217273082,+0.024811391363459,+0.024811391363459,+0.094853828379579,
     +0.268794997058761,+0.636351174561660,+0.268794997058761,+0.636351174561660,
     +0.094853828379579,+0.018100773278807,+0.291730066734288,+0.690169159986905,
     +0.291730066734288,+0.690169159986905,+0.018100773278807,+0.022233076674090,
     +0.126357385491669,+0.851409537834241,+0.126357385491669,+0.851409537834241,
     +0.022233076674090};
    static const double trians13[] =
    {+0.333333333333333,+0.495048184939705,+0.495048184939705,+0.009903630120591,
     +0.468716635109574,+0.468716635109574,+0.062566729780852,+0.414521336801277,
     +0.414521336801277,+0.170957326397447,+0.229399572042831,+0.229399572042831,
     +0.541200855914337,+0.114424495196330,+0.114424495196330,+0.771151009607340,
     +0.024811391363459,+0.024811391363459,+0.950377217273082,+0.268794997058761,
     +0.636351174561660,+0.094853828379579,+0.094853828379579,+0.268794997058761,
     +0.636351174561660,+0.291730066734288,+0.690169159986905,+0.018100773278807,
     +0.018100773278807,+0.291730066734288,+0.690169159986905,+0.126357385491669,
     +0.851409537834241,+0.022233076674090,+0.022233076674090,+0.126357385491669,
     +0.851409537834241};
    static const double trianw13[] =
    {+0.026260461700401,+0.005640072604665,+0.005640072604665,+0.005640072604665,
     +0.015711759181227,+0.015711759181227,+0.015711759181227,+0.023536251252097,
     +0.023536251252097,+0.023536251252097,+0.023681793268178,+0.023681793268178,
     +0.023681793268178,+0.015583764522897,+0.015583764522897,+0.015583764522897,
     +0.003987885732537,+0.003987885732537,+0.003987885732537,+0.018424201364366,
     +0.018424201364366,+0.018424201364366,+0.018424201364366,+0.018424201364366,
     +0.018424201364366,+0.008700731651911,+0.008700731651911,+0.008700731651911,
     +0.008700731651911,+0.008700731651911,+0.008700731651911,+0.007760893419522,
     +0.007760893419522,+0.007760893419522,+0.007760893419522,+0.007760893419522,
     +0.007760893419522};
    
    // n = 14
    
    static const double trianr14[] =
    {+0.022072179275643,+0.488963910362179,+0.488963910362179,+0.164710561319092,
     +0.417644719340454,+0.417644719340454,+0.453044943382323,+0.273477528308839,
     +0.273477528308839,+0.645588935174913,+0.177205532412543,+0.177205532412543,
     +0.876400233818255,+0.061799883090873,+0.061799883090873,+0.961218077502598,
     +0.019390961248701,+0.019390961248701,+0.057124757403648,+0.172266687821356,
     +0.770608554774996,+0.172266687821356,+0.770608554774996,+0.057124757403648,
     +0.092916249356972,+0.336861459796345,+0.570222290846683,+0.336861459796345,
     +0.570222290846683,+0.092916249356972,+0.014646950055654,+0.298372882136258,
     +0.686980167808088,+0.298372882136258,+0.686980167808088,+0.014646950055654,
     +0.001268330932872,+0.118974497696957,+0.879757171370171,+0.118974497696957,
     +0.879757171370171,+0.001268330932872};
    static const double trians14[] =
    {+0.488963910362179,+0.488963910362179,+0.022072179275643,+0.417644719340454,
     +0.417644719340454,+0.164710561319092,+0.273477528308839,+0.273477528308839,
     +0.453044943382323,+0.177205532412543,+0.177205532412543,+0.645588935174913,
     +0.061799883090873,+0.061799883090873,+0.876400233818255,+0.019390961248701,
     +0.019390961248701,+0.961218077502598,+0.172266687821356,+0.770608554774996,
     +0.057124757403648,+0.057124757403648,+0.172266687821356,+0.770608554774996,
     +0.336861459796345,+0.570222290846683,+0.092916249356972,+0.092916249356972,
     +0.336861459796345,+0.570222290846683,+0.298372882136258,+0.686980167808088,
     +0.014646950055654,+0.014646950055654,+0.298372882136258,+0.686980167808088,
     +0.118974497696957,+0.879757171370171,+0.001268330932872,+0.001268330932872,
     +0.118974497696957,+0.879757171370171};
    static const double trianw14[] =
    {+0.010941790684715,+0.010941790684715,+0.010941790684715,+0.016394176772063,
     +0.016394176772063,+0.016394176772063,+0.025887052253646,+0.025887052253646,
     +0.025887052253646,+0.021081294368497,+0.021081294368497,+0.021081294368497,
     +0.007216849834889,+0.007216849834889,+0.007216849834889,+0.002461701801200,
     +0.002461701801200,+0.002461701801200,+0.012332876606282,+0.012332876606282,
     +0.012332876606282,+0.012332876606282,+0.012332876606282,+0.012332876606282,
     +0.019285755393531,+0.019285755393531,+0.019285755393531,+0.019285755393531,
     +0.019285755393531,+0.019285755393531,+0.007218154056767,+0.007218154056767,
     +0.007218154056767,+0.007218154056767,+0.007218154056767,+0.007218154056767,
     +0.002505114419250,+0.002505114419250,+0.002505114419250,+0.002505114419250,
     +0.002505114419250,+0.002505114419250};
    
    // n = 15
    
    static const double trianr15[] =
    {-0.013945833716486,+0.506972916858243,+0.506972916858243,+0.137187291433955,
     +0.431406354283023,+0.431406354283023,+0.444612710305711,+0.277693644847144,
     +0.277693644847144,+0.747070217917492,+0.126464891041254,+0.126464891041254,
     +0.858383228050628,+0.070808385974686,+0.070808385974686,+0.962069659517853,
     +0.018965170241073,+0.018965170241073,+0.133734161966621,+0.261311371140087,
     +0.604954466893291,+0.261311371140087,+0.604954466893291,+0.133734161966621,
     +0.036366677396917,+0.388046767090269,+0.575586555512814,+0.388046767090269,
     +0.575586555512814,+0.036366677396917,-0.010174883126571,+0.285712220049916,
     +0.724462663076655,+0.285712220049916,+0.724462663076655,-0.010174883126571,
     +0.036843869875878,+0.215599664072284,+0.747556466051838,+0.215599664072284,
     +0.747556466051838,+0.036843869875878,+0.012459809331199,+0.103575616576386,
     +0.883964574092416,+0.103575616576386,+0.883964574092416,+0.012459809331199};
    static const double trians15[] =
    {+0.506972916858243,+0.506972916858243,-0.013945833716486,+0.431406354283023,
     +0.431406354283023,+0.137187291433955,+0.277693644847144,+0.277693644847144,
     +0.444612710305711,+0.126464891041254,+0.126464891041254,+0.747070217917492,
     +0.070808385974686,+0.070808385974686,+0.858383228050628,+0.018965170241073,
     +0.018965170241073,+0.962069659517853,+0.261311371140087,+0.604954466893291,
     +0.133734161966621,+0.133734161966621,+0.261311371140087,+0.604954466893291,
     +0.388046767090269,+0.575586555512814,+0.036366677396917,+0.036366677396917,
     +0.388046767090269,+0.575586555512814,+0.285712220049916,+0.724462663076655,
     -0.010174883126571,-0.010174883126571,+0.285712220049916,+0.724462663076655,
     +0.215599664072284,+0.747556466051838,+0.036843869875878,+0.036843869875878,
     +0.215599664072284,+0.747556466051838,+0.103575616576386,+0.883964574092416,
     +0.012459809331199,+0.012459809331199,+0.103575616576386,+0.883964574092416};
    static const double trianw15[] =
    {+0.000958437821424,+0.000958437821424,+0.000958437821424,+0.022124513635572,
     +0.022124513635572,+0.022124513635572,+0.025593274359426,+0.025593274359426,
     +0.025593274359426,+0.011843867935344,+0.011843867935344,+0.011843867935344,
     +0.006644887845011,+0.006644887845011,+0.006644887845011,+0.002374458304096,
     +0.002374458304096,+0.002374458304096,+0.019275036299796,+0.019275036299796,
     +0.019275036299796,+0.019275036299796,+0.019275036299796,+0.019275036299796,
     +0.013607907160312,+0.013607907160312,+0.013607907160312,+0.013607907160312,
     +0.013607907160312,+0.013607907160312,+0.001091038683399,+0.001091038683399,
     +0.001091038683399,+0.001091038683399,+0.001091038683399,+0.001091038683399,
     +0.010752659923865,+0.010752659923865,+0.010752659923865,+0.010752659923865,
     +0.010752659923865,+0.010752659923865,+0.003836971315524,+0.003836971315524,
     +0.003836971315524,+0.003836971315524,+0.003836971315524,+0.003836971315524};
    
    // n = 16
    
    static const double trianr16[] =
    {+0.333333333333333,+0.005238916103123,+0.497380541948438,+0.497380541948438,
     +0.173061122901295,+0.413469438549352,+0.413469438549352,+0.059082801866017,
     +0.470458599066991,+0.470458599066991,+0.518892500060958,+0.240553749969521,
     +0.240553749969521,+0.704068411554854,+0.147965794222573,+0.147965794222573,
     +0.849069624685052,+0.075465187657474,+0.075465187657474,+0.966807194753950,
     +0.016596402623025,+0.016596402623025,+0.103575692245252,+0.296555596579887,
     +0.599868711174861,+0.296555596579887,+0.599868711174861,+0.103575692245252,
     +0.020083411655416,+0.337723063403079,+0.642193524941505,+0.337723063403079,
     +0.642193524941505,+0.020083411655416,-0.004341002614139,+0.204748281642812,
     +0.799592720971327,+0.204748281642812,+0.799592720971327,-0.004341002614139,
     +0.041941786468010,+0.189358492130623,+0.768699721401368,+0.189358492130623,
     +0.768699721401368,+0.041941786468010,+0.014317320230681,+0.085283615682657,
     +0.900399064086661,+0.085283615682657,+0.900399064086661,+0.014317320230681};
    static const double trians16[] =
    {+0.333333333333333,+0.497380541948438,+0.497380541948438,+0.005238916103123,
     +0.413469438549352,+0.413469438549352,+0.173061122901295,+0.470458599066991,
     +0.470458599066991,+0.059082801866017,+0.240553749969521,+0.240553749969521,
     +0.518892500060958,+0.147965794222573,+0.147965794222573,+0.704068411554854,
     +0.075465187657474,+0.075465187657474,+0.849069624685052,+0.016596402623025,
     +0.016596402623025,+0.966807194753950,+0.296555596579887,+0.599868711174861,
     +0.103575692245252,+0.103575692245252,+0.296555596579887,+0.599868711174861,
     +0.337723063403079,+0.642193524941505,+0.020083411655416,+0.020083411655416,
     +0.337723063403079,+0.642193524941505,+0.204748281642812,+0.799592720971327,
     -0.004341002614139,-0.004341002614139,+0.204748281642812,+0.799592720971327,
     +0.189358492130623,+0.768699721401368,+0.041941786468010,+0.041941786468010,
     +0.189358492130623,+0.768699721401368,+0.085283615682657,+0.900399064086661,
     +0.014317320230681,+0.014317320230681,+0.085283615682657,+0.900399064086661};
    static const double trianw16[] =
    {+0.023437848713821,+0.003202939289292,+0.003202939289292,+0.003202939289292,
     +0.020855148369693,+0.020855148369693,+0.020855148369693,+0.013445742125032,
     +0.013445742125032,+0.013445742125032,+0.021066261380825,+0.021066261380825,
     +0.021066261380825,+0.015000133421386,+0.015000133421386,+0.015000133421386,
     +0.007100049462512,+0.007100049462512,+0.007100049462512,+0.001791231175637,
     +0.001791231175637,+0.001791231175637,+0.016386573730313,+0.016386573730313,
     +0.016386573730313,+0.016386573730313,+0.016386573730313,+0.016386573730313,
     +0.007649153124220,+0.007649153124220,+0.007649153124220,+0.007649153124220,
     +0.007649153124220,+0.007649153124220,+0.001193122096420,+0.001193122096420,
     +0.001193122096420,+0.001193122096420,+0.001193122096420,+0.001193122096420,
     +0.009542396377949,+0.009542396377949,+0.009542396377949,+0.009542396377949,
     +0.009542396377949,+0.009542396377949,+0.003425027273271,+0.003425027273271,
     +0.003425027273271,+0.003425027273271,+0.003425027273271,+0.003425027273271};
    
    // n = 17
    
    static const double trianr17[] =
    {+0.333333333333333,+0.005658918886452,+0.497170540556774,+0.497170540556774,
     +0.035647354750751,+0.482176322624625,+0.482176322624625,+0.099520061958437,
     +0.450239969020782,+0.450239969020782,+0.199467521245206,+0.400266239377397,
     +0.400266239377397,+0.495717464058095,+0.252141267970953,+0.252141267970953,
     +0.675905990683077,+0.162047004658461,+0.162047004658461,+0.848248235478508,
     +0.075875882260746,+0.075875882260746,+0.968690546064356,+0.015654726967822,
     +0.015654726967822,+0.010186928826919,+0.334319867363658,+0.655493203809423,
     +0.334319867363658,+0.655493203809423,+0.010186928826919,+0.135440871671036,
     +0.292221537796944,+0.572337590532020,+0.292221537796944,+0.572337590532020,
     +0.135440871671036,+0.054423924290583,+0.319574885423190,+0.626001190286228,
     +0.319574885423190,+0.626001190286228,+0.054423924290583,+0.012868560833637,
     +0.190704224192292,+0.796427214974071,+0.190704224192292,+0.796427214974071,
     +0.012868560833637,+0.067165782413524,+0.180483211648746,+0.752351005937729,
     +0.180483211648746,+0.752351005937729,+0.067165782413524,+0.014663182224828,
     +0.080711313679564,+0.904625504095608,+0.080711313679564,+0.904625504095608,
     +0.014663182224828};
    static const double trians17[] =
    {+0.333333333333333,+0.497170540556774,+0.497170540556774,+0.005658918886452,
     +0.482176322624625,+0.482176322624625,+0.035647354750751,+0.450239969020782,
     +0.450239969020782,+0.099520061958437,+0.400266239377397,+0.400266239377397,
     +0.199467521245206,+0.252141267970953,+0.252141267970953,+0.495717464058095,
     +0.162047004658461,+0.162047004658461,+0.675905990683077,+0.075875882260746,
     +0.075875882260746,+0.848248235478508,+0.015654726967822,+0.015654726967822,
     +0.968690546064356,+0.334319867363658,+0.655493203809423,+0.010186928826919,
     +0.010186928826919,+0.334319867363658,+0.655493203809423,+0.292221537796944,
     +0.572337590532020,+0.135440871671036,+0.135440871671036,+0.292221537796944,
     +0.572337590532020,+0.319574885423190,+0.626001190286228,+0.054423924290583,
     +0.054423924290583,+0.319574885423190,+0.626001190286228,+0.190704224192292,
     +0.796427214974071,+0.012868560833637,+0.012868560833637,+0.190704224192292,
     +0.796427214974071,+0.180483211648746,+0.752351005937729,+0.067165782413524,
     +0.067165782413524,+0.180483211648746,+0.752351005937729,+0.080711313679564,
     +0.904625504095608,+0.014663182224828,+0.014663182224828,+0.080711313679564,
     +0.904625504095608};
    static const double trianw17[] =
    {+0.016718599645402,+0.002546707720254,+0.002546707720254,+0.002546707720254,
     +0.007335432263819,+0.007335432263819,+0.007335432263819,+0.012175439176836,
     +0.012175439176836,+0.012175439176836,+0.015553775434484,+0.015553775434484,
     +0.015553775434484,+0.015628555609310,+0.015628555609310,+0.015628555609310,
     +0.012407827169833,+0.012407827169833,+0.012407827169833,+0.007028036535279,
     +0.007028036535279,+0.007028036535279,+0.001597338086889,+0.001597338086889,
     +0.001597338086889,+0.004059827659497,+0.004059827659497,+0.004059827659497,
     +0.004059827659497,+0.004059827659497,+0.004059827659497,+0.013402871141582,
     +0.013402871141582,+0.013402871141582,+0.013402871141582,+0.013402871141582,
     +0.013402871141582,+0.009229996605411,+0.009229996605411,+0.009229996605411,
     +0.009229996605411,+0.009229996605411,+0.009229996605411,+0.004238434267164,
     +0.004238434267164,+0.004238434267164,+0.004238434267164,+0.004238434267164,
     +0.004238434267164,+0.009146398385012,+0.009146398385012,+0.009146398385012,
     +0.009146398385012,+0.009146398385012,+0.009146398385012,+0.003332816002083,
     +0.003332816002083,+0.003332816002083,+0.003332816002083,+0.003332816002083,
     +0.003332816002083};
    
    // n = 18
    
    static const double trianr18[] =
    {+0.333333333333333,+0.013310382738157,+0.493344808630921,+0.493344808630921,
     +0.061578811516086,+0.469210594241957,+0.469210594241957,+0.127437208225989,
     +0.436281395887006,+0.436281395887006,+0.210307658653168,+0.394846170673416,
     +0.394846170673416,+0.500410862393686,+0.249794568803157,+0.249794568803157,
     +0.677135612512315,+0.161432193743843,+0.161432193743843,+0.846803545029257,
     +0.076598227485371,+0.076598227485371,+0.951495121293100,+0.024252439353450,
     +0.024252439353450,+0.913707265566071,+0.043146367216965,+0.043146367216965,
     +0.008430536202420,+0.358911494940944,+0.632657968856636,+0.358911494940944,
     +0.632657968856636,+0.008430536202420,+0.131186551737188,+0.294402476751957,
     +0.574410971510855,+0.294402476751957,+0.574410971510855,+0.131186551737188,
     +0.050203151565675,+0.325017801641814,+0.624779046792512,+0.325017801641814,
     +0.624779046792512,+0.050203151565675,+0.066329263810916,+0.184737559666046,
     +0.748933176523037,+0.184737559666046,+0.748933176523037,+0.066329263810916,
     +0.011996194566236,+0.218796800013321,+0.769207005420443,+0.218796800013321,
     +0.769207005420443,+0.011996194566236,+0.014858100590125,+0.101179597136408,
     +0.883962302273467,+0.101179597136408,+0.883962302273467,+0.014858100590125,
     -0.035222015287949,+0.020874755282586,+1.014347260005363,+0.020874755282586,
     +1.014347260005363,-0.035222015287949};
    static const double trians18[] =
    {+0.333333333333333,+0.493344808630921,+0.493344808630921,+0.013310382738157,
     +0.469210594241957,+0.469210594241957,+0.061578811516086,+0.436281395887006,
     +0.436281395887006,+0.127437208225989,+0.394846170673416,+0.394846170673416,
     +0.210307658653168,+0.249794568803157,+0.249794568803157,+0.500410862393686,
     +0.161432193743843,+0.161432193743843,+0.677135612512315,+0.076598227485371,
     +0.076598227485371,+0.846803545029257,+0.024252439353450,+0.024252439353450,
     +0.951495121293100,+0.043146367216965,+0.043146367216965,+0.913707265566071,
     +0.358911494940944,+0.632657968856636,+0.008430536202420,+0.008430536202420,
     +0.358911494940944,+0.632657968856636,+0.294402476751957,+0.574410971510855,
     +0.131186551737188,+0.131186551737188,+0.294402476751957,+0.574410971510855,
     +0.325017801641814,+0.624779046792512,+0.050203151565675,+0.050203151565675,
     +0.325017801641814,+0.624779046792512,+0.184737559666046,+0.748933176523037,
     +0.066329263810916,+0.066329263810916,+0.184737559666046,+0.748933176523037,
     +0.218796800013321,+0.769207005420443,+0.011996194566236,+0.011996194566236,
     +0.218796800013321,+0.769207005420443,+0.101179597136408,+0.883962302273467,
     +0.014858100590125,+0.014858100590125,+0.101179597136408,+0.883962302273467,
     +0.020874755282586,+1.014347260005363,-0.035222015287949,-0.035222015287949,
     +0.020874755282586,+1.014347260005363};
    static const double trianw18[] =
    {+0.015404969968823,+0.004536218339702,+0.004536218339702,+0.004536218339702,
     +0.009380658469797,+0.009380658469797,+0.009380658469797,+0.009720548992738,
     +0.009720548992738,+0.009720548992738,+0.013876974305405,+0.013876974305405,
     +0.013876974305405,+0.016128112675729,+0.016128112675729,+0.016128112675729,
     +0.012537016308461,+0.012537016308461,+0.012537016308461,+0.007635963985916,
     +0.007635963985916,+0.007635963985916,+0.003396961011482,+0.003396961011482,
     +0.003396961011482,-0.001111549364960,-0.001111549364960,-0.001111549364960,
     +0.003165957038203,+0.003165957038203,+0.003165957038203,+0.003165957038203,
     +0.003165957038203,+0.003165957038203,+0.013628769024569,+0.013628769024569,
     +0.013628769024569,+0.013628769024569,+0.013628769024569,+0.013628769024569,
     +0.008838392824732,+0.008838392824732,+0.008838392824732,+0.008838392824732,
     +0.008838392824732,+0.008838392824732,+0.009189742319035,+0.009189742319035,
     +0.009189742319035,+0.009189742319035,+0.009189742319035,+0.009189742319035,
     +0.004052366404096,+0.004052366404096,+0.004052366404096,+0.004052366404096,
     +0.004052366404096,+0.004052366404096,+0.003817064535363,+0.003817064535363,
     +0.003817064535363,+0.003817064535363,+0.003817064535363,+0.003817064535363,
     +0.000023093830397,+0.000023093830397,+0.000023093830397,+0.000023093830397,
     +0.000023093830397,+0.000023093830397};
    
    // n = 19
    
    static const double trianr19[] =
    {+0.333333333333333,+0.020780025853987,+0.489609987073006,+0.489609987073006,
     +0.090926214604215,+0.454536892697893,+0.454536892697893,+0.197166638701138,
     +0.401416680649431,+0.401416680649431,+0.488896691193805,+0.255551654403098,
     +0.255551654403098,+0.645844115695741,+0.177077942152130,+0.177077942152130,
     +0.779877893544096,+0.110061053227952,+0.110061053227952,+0.888942751496321,
     +0.055528624251840,+0.055528624251840,+0.974756272445543,+0.012621863777229,
     +0.012621863777229,+0.003611417848412,+0.395754787356943,+0.600633794794645,
     +0.395754787356943,+0.600633794794645,+0.003611417848412,+0.134466754530780,
     +0.307929983880436,+0.557603261588784,+0.307929983880436,+0.557603261588784,
     +0.134466754530780,+0.014446025776115,+0.264566948406520,+0.720987025817365,
     +0.264566948406520,+0.720987025817365,+0.014446025776115,+0.046933578838178,
     +0.358539352205951,+0.594527068955871,+0.358539352205951,+0.594527068955871,
     +0.046933578838178,+0.002861120350567,+0.157807405968595,+0.839331473680839,
     +0.157807405968595,+0.839331473680839,+0.002861120350567,+0.223861424097916,
     +0.075050596975911,+0.701087978926173,+0.075050596975911,+0.701087978926173,
     +0.223861424097916,+0.034647074816760,+0.142421601113383,+0.822931324069857,
     +0.142421601113383,+0.822931324069857,+0.034647074816760,+0.010161119296278,
     +0.065494628082938,+0.924344252620784,+0.065494628082938,+0.924344252620784,
     +0.010161119296278};
    static const double trians19[] =
    {+0.333333333333333,+0.489609987073006,+0.489609987073006,+0.020780025853987,
     +0.454536892697893,+0.454536892697893,+0.090926214604215,+0.401416680649431,
     +0.401416680649431,+0.197166638701138,+0.255551654403098,+0.255551654403098,
     +0.488896691193805,+0.177077942152130,+0.177077942152130,+0.645844115695741,
     +0.110061053227952,+0.110061053227952,+0.779877893544096,+0.055528624251840,
     +0.055528624251840,+0.888942751496321,+0.012621863777229,+0.012621863777229,
     +0.974756272445543,+0.395754787356943,+0.600633794794645,+0.003611417848412,
     +0.003611417848412,+0.395754787356943,+0.600633794794645,+0.307929983880436,
     +0.557603261588784,+0.134466754530780,+0.134466754530780,+0.307929983880436,
     +0.557603261588784,+0.264566948406520,+0.720987025817365,+0.014446025776115,
     +0.014446025776115,+0.264566948406520,+0.720987025817365,+0.358539352205951,
     +0.594527068955871,+0.046933578838178,+0.046933578838178,+0.358539352205951,
     +0.594527068955871,+0.157807405968595,+0.839331473680839,+0.002861120350567,
     +0.002861120350567,+0.157807405968595,+0.839331473680839,+0.075050596975911,
     +0.701087978926173,+0.223861424097916,+0.223861424097916,+0.075050596975911,
     +0.701087978926173,+0.142421601113383,+0.822931324069857,+0.034647074816760,
     +0.034647074816760,+0.142421601113383,+0.822931324069857,+0.065494628082938,
     +0.924344252620784,+0.010161119296278,+0.010161119296278,+0.065494628082938,
     +0.924344252620784};
    static const double trianw19[] =
    {+0.016453165694459,+0.005165365945636,+0.005165365945636,+0.005165365945636,
     +0.011193623631508,+0.011193623631508,+0.011193623631508,+0.015133062934734,
     +0.015133062934734,+0.015133062934734,+0.015245483901099,+0.015245483901099,
     +0.015245483901099,+0.012079606370821,+0.012079606370821,+0.012079606370821,
     +0.008025401793400,+0.008025401793400,+0.008025401793400,+0.004042290130892,
     +0.004042290130892,+0.004042290130892,+0.001039681013742,+0.001039681013742,
     +0.001039681013742,+0.001942438452491,+0.001942438452491,+0.001942438452491,
     +0.001942438452491,+0.001942438452491,+0.001942438452491,+0.012787080306011,
     +0.012787080306011,+0.012787080306011,+0.012787080306011,+0.012787080306011,
     +0.012787080306011,+0.004440451786669,+0.004440451786669,+0.004440451786669,
     +0.004440451786669,+0.004440451786669,+0.004440451786669,+0.008062273380866,
     +0.008062273380866,+0.008062273380866,+0.008062273380866,+0.008062273380866,
     +0.008062273380866,+0.001245970908745,+0.001245970908745,+0.001245970908745,
     +0.001245970908745,+0.001245970908745,+0.001245970908745,+0.009121420059476,
     +0.009121420059476,+0.009121420059476,+0.009121420059476,+0.009121420059476,
     +0.009121420059476,+0.005129281868099,+0.005129281868099,+0.005129281868099,
     +0.005129281868099,+0.005129281868099,+0.005129281868099,+0.001899964427651,
     +0.001899964427651,+0.001899964427651,+0.001899964427651,+0.001899964427651,
     +0.001899964427651};
    
    // n = 20
    
    static const double trianr20[] =
    {+0.333333333333333,-0.001900928704400,+0.500950464352200,+0.500950464352200,
     +0.023574084130543,+0.488212957934729,+0.488212957934729,+0.089726636099435,
     +0.455136681950283,+0.455136681950283,+0.196007481363421,+0.401996259318289,
     +0.401996259318289,+0.488214180481157,+0.255892909759421,+0.255892909759421,
     +0.647023488009788,+0.176488255995106,+0.176488255995106,+0.791658289326483,
     +0.104170855336758,+0.104170855336758,+0.893862072318140,+0.053068963840930,
     +0.053068963840930,+0.916762569607942,+0.041618715196029,+0.041618715196029,
     +0.976836157186356,+0.011581921406822,+0.011581921406822,+0.048741583664839,
     +0.344855770229001,+0.606402646106160,+0.344855770229001,+0.606402646106160,
     +0.048741583664839,+0.006314115948605,+0.377843269594854,+0.615842614456541,
     +0.377843269594854,+0.615842614456541,+0.006314115948605,+0.134316520547348,
     +0.306635479062357,+0.559048000390295,+0.306635479062357,+0.559048000390295,
     +0.134316520547348,+0.013973893962392,+0.249419362774742,+0.736606743262866,
     +0.249419362774742,+0.736606743262866,+0.013973893962392,+0.075549132909764,
     +0.212775724802802,+0.711675142287434,+0.212775724802802,+0.711675142287434,
     +0.075549132909764,-0.008368153208227,+0.146965436053239,+0.861402717154987,
     +0.146965436053239,+0.861402717154987,-0.008368153208227,+0.026686063258714,
     +0.137726978828923,+0.835586957912363,+0.137726978828923,+0.835586957912363,
     +0.026686063258714,+0.010547719294141,+0.059696109149007,+0.929756171556853,
     +0.059696109149007,+0.929756171556853,+0.010547719294141};
    static const double trians20[] =
    {+0.333333333333333,+0.500950464352200,+0.500950464352200,-0.001900928704400,
     +0.488212957934729,+0.488212957934729,+0.023574084130543,+0.455136681950283,
     +0.455136681950283,+0.089726636099435,+0.401996259318289,+0.401996259318289,
     +0.196007481363421,+0.255892909759421,+0.255892909759421,+0.488214180481157,
     +0.176488255995106,+0.176488255995106,+0.647023488009788,+0.104170855336758,
     +0.104170855336758,+0.791658289326483,+0.053068963840930,+0.053068963840930,
     +0.893862072318140,+0.041618715196029,+0.041618715196029,+0.916762569607942,
     +0.011581921406822,+0.011581921406822,+0.976836157186356,+0.344855770229001,
     +0.606402646106160,+0.048741583664839,+0.048741583664839,+0.344855770229001,
     +0.606402646106160,+0.377843269594854,+0.615842614456541,+0.006314115948605,
     +0.006314115948605,+0.377843269594854,+0.615842614456541,+0.306635479062357,
     +0.559048000390295,+0.134316520547348,+0.134316520547348,+0.306635479062357,
     +0.559048000390295,+0.249419362774742,+0.736606743262866,+0.013973893962392,
     +0.013973893962392,+0.249419362774742,+0.736606743262866,+0.212775724802802,
     +0.711675142287434,+0.075549132909764,+0.075549132909764,+0.212775724802802,
     +0.711675142287434,+0.146965436053239,+0.861402717154987,-0.008368153208227,
     -0.008368153208227,+0.146965436053239,+0.861402717154987,+0.137726978828923,
     +0.835586957912363,+0.026686063258714,+0.026686063258714,+0.137726978828923,
     +0.835586957912363,+0.059696109149007,+0.929756171556853,+0.010547719294141,
     +0.010547719294141,+0.059696109149007,+0.929756171556853};
    static const double trianw20[] =
    {+0.016528527770812,+0.000433509592831,+0.000433509592831,+0.000433509592831,
     +0.005830026358224,+0.005830026358224,+0.005830026358224,+0.011438468178211,
     +0.011438468178211,+0.011438468178211,+0.015224491336969,+0.015224491336969,
     +0.015224491336969,+0.015312445862677,+0.015312445862677,+0.015312445862677,
     +0.012184028838400,+0.012184028838400,+0.012184028838400,+0.007998716016012,
     +0.007998716016012,+0.007998716016012,+0.003849150907801,+0.003849150907801,
     +0.003849150907801,-0.000316030248744,-0.000316030248744,-0.000316030248744,
     +0.000875567150596,+0.000875567150596,+0.000875567150596,+0.008232919594788,
     +0.008232919594788,+0.008232919594788,+0.008232919594788,+0.008232919594788,
     +0.008232919594788,+0.002419516770243,+0.002419516770243,+0.002419516770243,
     +0.002419516770243,+0.002419516770243,+0.002419516770243,+0.012902453267325,
     +0.012902453267325,+0.012902453267325,+0.012902453267325,+0.012902453267325,
     +0.012902453267325,+0.004235545527221,+0.004235545527221,+0.004235545527221,
     +0.004235545527221,+0.004235545527221,+0.004235545527221,+0.009177457053140,
     +0.009177457053140,+0.009177457053140,+0.009177457053140,+0.009177457053140,
     +0.009177457053140,+0.000352202338954,+0.000352202338954,+0.000352202338954,
     +0.000352202338954,+0.000352202338954,+0.000352202338954,+0.005056342463731,
     +0.005056342463731,+0.005056342463731,+0.005056342463731,+0.005056342463731,
     +0.005056342463731,+0.001786954692975,+0.001786954692975,+0.001786954692975,
     +0.001786954692975,+0.001786954692975,+0.001786954692975};
    
    // Array of Dunavant data.
    
    data.push_back(sIntData(trianr1,trians1,0,trianw1,1));
    data.push_back(sIntData(trianr2,trians2,0,trianw2,3));
    data.push_back(sIntData(trianr3,trians3,0,trianw3,4));
    data.push_back(sIntData(trianr4,trians4,0,trianw4,6));
    data.push_back(sIntData(trianr5,trians5,0,trianw5,7));
    data.push_back(sIntData(trianr6,trians6,0,trianw6,12));
    data.push_back(sIntData(trianr7,trians7,0,trianw7,13));
    data.push_back(sIntData(trianr8,trians8,0,trianw8,16));
    data.push_back(sIntData(trianr9,trians9,0,trianw9,19));
    data.push_back(sIntData(trianr10,trians10,0,trianw10,25));
    data.push_back(sIntData(trianr11,trians11,0,trianw11,27));
    data.push_back(sIntData(trianr12,trians12,0,trianw12,33));
    data.push_back(sIntData(trianr13,trians13,0,trianw13,37));
    data.push_back(sIntData(trianr14,trians14,0,trianw14,42));
    data.push_back(sIntData(trianr15,trians15,0,trianw15,48));
    data.push_back(sIntData(trianr16,trians16,0,trianw16,52));
    data.push_back(sIntData(trianr17,trians17,0,trianw17,61));
    data.push_back(sIntData(trianr18,trians18,0,trianw18,70));
    data.push_back(sIntData(trianr19,trians19,0,trianw19,73));
    data.push_back(sIntData(trianr20,trians20,0,trianw20,79));
  }

  return data;
}

// ============================== getLynessData ============================

const vector<sIntData>& getLynessData( )
{
  static vector<sIntData> data;

  if (data.empty( )) 
  {
    // -------------------------------------------------------------------------
    // Dunavant data:
    //
    //
    //
    // n = 1, rule = 0
    
    static const double lynr1[] =
    {0.3333333333000000};
    
    static const double lyns1[] =
    {0.3333333333000000};
    
    static const double lynw1[] =
    {0.5000000000000000};
    
    
    // n = 2, rule = 1
    
    static const double lynr2[] =
    {0.0000000000000000,0.5000000000000000,0.5000000000000000};
    
    static const double lyns2[] =
    {0.5000000000000000,0.5000000000000000,0.0000000000000000};
    
    static const double lynw2[] =
    {0.1666666666666667,0.1666666666666667,0.1666666666666667};
    
    
    // n = 3, rule = 3
   /* 
    static const double lynr3[] =
    {0.3333333333000000,0.6000000000000000,0.2000000000000000,0.2000000000000000};
    
    static const double lyns3[] =
    {0.3333333333000000,0.2000000000000000,0.2000000000000000,0.6000000000000000};
    
    static const double lynw3[] =
    {-0.2812500000000000,0.2604166666666667,0.2604166666666667,0.2604166666666667};
    */
    
    
    // n = 3, rule = 4
    
    static const double lynr3[] =
    {0.3333333333333333,1.0000000000000000,0.0000000000000000,0.0000000000000000,
     0.0000000000000000,0.5000000000000000,0.5000000000000000};
    
    static const double lyns3[] =
    {0.3333333333333333,0.0000000000000000,0.0000000000000000,1.0000000000000000,
     0.5000000000000000,0.5000000000000000,0.0000000000000000};
    
    static const double lynw3[] =
    {0.2250000000000000,0.0250000000000000,0.0250000000000000,0.0250000000000000,
     0.0666666666666667,0.0666666666666667,0.0666666666666667};
    
    
    // n = 4, rule = 5
    
    /*
    static const double lynr4[] =
    {0.8168475729804585,0.0915762135097707,0.0915762135097708,0.1081030181680702,
     0.4459484909159649,0.4459484909159649};
    
    static const double lyns4[] =
    {0.0915762135097707,0.0915762135097708,0.8168475729804585,0.4459484909159649,
     0.4459484909159649,0.1081030181680702};
    
    static const double lynw4[] =
    {0.0549758718276609,0.0549758718276609,0.0549758718276609,0.1116907948390057,
     0.1116907948390057,0.1116907948390057};
     */
    
    
    // n = 4, rule = 6
    
    static const double lynr4[] =
    {0.3333333333333333,1.0000000000000000,0.0000000000000000,0.0000000000000000,
     0.0000000000000000,0.7886751345948128,0.2113248654051872,0.7886751345948128,
     0.2113248654051872,0.0000000000000000};
    
    static const double lyns4[] =
    {0.3333333333333333,0.0000000000000000,0.0000000000000000,1.0000000000000000,
     0.7886751345948128,0.2113248654051872,0.0000000000000000,0.0000000000000000,
     0.7886751345948128,0.2113248654051872};
    
    static const double lynw4[] =
    {0.2250000000000000,-0.0083333333333333,-0.0083333333333333,-0.0083333333333333,
     0.0500000000000000,0.0500000000000000,0.0500000000000000,0.0500000000000000,
     0.0500000000000000,0.0500000000000000};
    
    
    // n = 4, rule = 7
    /*
    static const double lynr4[] =
    {1.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,
     0.5000000000000000,0.5000000000000000,0.6228390306071100,0.1885804846964451,
     0.1885804846964451};
    
    static const double lyns4[] =
    {0.0000000000000000,0.0000000000000000,1.0000000000000000,0.5000000000000000,
     0.5000000000000000,0.0000000000000000,0.1885804846964451,0.1885804846964451,
     0.6228390306071100};
    
    static const double lynw4[] =
    {0.0102700676729667,0.0102700676729667,0.0102700676729667,0.0309877494341336,
     0.0309877494341336,0.0309877494341336,0.1254088495595664,0.1254088495595664,
     0.1254088495595664};
     */
    
    
    // n = 5, rule = 8
    
    /*
    static const double lynr5[] =
    {0.3333333333333333,0.7974269853530872,0.1012865073234563,0.1012865073234563,
     0.0597158717897698,0.4701420641051151,0.4701420641051151};
    
    static const double lyns5[] =
    {0.3333333333333333,0.1012865073234563,0.1012865073234563,0.7974269853530872,
     0.4701420641051151,0.4701420641051151,0.0597158717897698};
    
    static const double lynw5[] =
    {0.1125000000000000,0.0629695902724136,0.0629695902724136,0.0629695902724136,
     0.0661970763942531,0.0661970763942531,0.0661970763942531};
    */
    
    // n = 5, rule = 9
    
    static const double lynr5[] =
    {0.3333333333333333,1.0000000000000000,0.0000000000000000,0.0000000000000000,
     0.0000000000000000,0.5000000000000000,0.5000000000000000,0.7142857142857143,
     0.1428571428571428,0.1428571428571428};
    
    static const double lyns5[] =
    {0.3333333333333333,0.0000000000000000,0.0000000000000000,1.0000000000000000,
     0.5000000000000000,0.5000000000000000,0.0000000000000000,0.1428571428571428,
     0.1428571428571428,0.7142857142857143};
    
    static const double lynw5[] =
    {0.1265625000000000,0.0055555555555556,0.0055555555555556,0.0055555555555556,
     0.0355555555555556,0.0355555555555556,0.0355555555555556,0.0833680555555556,
     0.0833680555555556,0.0833680555555556};
    
    
    // n = 6, rule = 10
   /* 
    static const double lynr6[] =
    {0.5014265096581342,0.2492867451709329,0.2492867451709329,0.8738219710169965,
     0.0630890144915018,0.0630890144915017,0.6365024991213939,0.0531450498448322,
     0.3103524510337740,0.0531450498448322,0.3103524510337740,0.6365024991213939};
    
    static const double lyns6[] =
    {0.2492867451709329,0.2492867451709329,0.5014265096581342,0.0630890144915018,
     0.0630890144915017,0.8738219710169965,0.0531450498448322,0.3103524510337740,
     0.6365024991213939,0.6365024991213939,0.0531450498448322,0.3103524510337740};
    
    static const double lynw6[] =
    {0.0583931378631704,0.0583931378631704,0.0583931378631704,0.0254224531851027,
     0.0254224531851027,0.0254224531851027,0.0414255378091965,0.0414255378091965,
     0.0414255378091965,0.0414255378091965,0.0414255378091965,0.0414255378091965};
    */
    // n = 6, rule = 12
    
    static const double lynr6[] =
    {0.3333333333333333,0.0523383720926975,0.4738308139536513,0.4738308139536513,
     0.6557646607383649,0.1721176696308175,0.1721176696308176,0.0000000000000000,
     0.8653073540834570,0.1346926459165430,0.8653073540834570,0.1346926459165430,
     0.0000000000000000};
    
    static const double lyns6[] =
    {0.3333333333333333,0.4738308139536513,0.4738308139536513,0.0523383720926975,
     0.1721176696308175,0.1721176696308176,0.6557646607383649,0.8653073540834570,
     0.1346926459165430,0.0000000000000000,0.0000000000000000,0.8653073540834570,
     0.1346926459165430};
    
    static const double lynw6[] =
    {0.0763544833941762,0.0490679340394460,0.0490679340394460,0.0490679340394460,
     0.0647842146403128,0.0647842146403128,0.0647842146403128,0.0136815117610912,
     0.0136815117610912,0.0136815117610912,0.0136815117610912,0.0136815117610912,
     0.0136815117610912};
    
    
    // n = 7, rule = 13
    /*
    static const double lynr7[] =
    {0.3333333333333333,0.4793080678419067,0.2603459660790466,0.2603459660790467,
     0.8697397941955675,0.0651301029022162,0.0651301029022163,0.6384441885698096,
     0.0486903154253176,0.3128654960048729,0.0486903154253176,0.3128654960048729,
     0.6384441885698096};
    
    static const double lyns7[] =
    {0.3333333333333333,0.2603459660790466,0.2603459660790467,0.4793080678419067,
     0.0651301029022162,0.0651301029022163,0.8697397941955675,0.0486903154253176,
     0.3128654960048729,0.6384441885698096,0.6384441885698096,0.0486903154253176,
     0.3128654960048729};
    
    static const double lynw7[] =
    {-0.0747850222338747,0.0878076287166138,0.0878076287166138,0.0878076287166138,
     0.0266736178044194,0.0266736178044194,0.0266736178044194,0.0385568804451292,
     0.0385568804451292,0.0385568804451292,0.0385568804451292,0.0385568804451292,
     0.0385568804451292};
    */
    
    
    // n = 7, rule = 14
    
    static const double lynr7[] =
    {0.3333333333333333,1.0000000000000000,0.0000000000000000,0.0000000000000000,
     0.6901278795524791,0.1549360602237604,0.1549360602237605,0.0616985077123759,
     0.4691507461438120,0.4691507461438120,0.0000000000000000,0.8392991722729236,
     0.1607008277270764,0.8392991722729236,0.1607008277270764,0.0000000000000000};
    
    static const double lyns7[] =
    {0.3333333333333333,0.0000000000000000,0.0000000000000000,1.0000000000000000,
     0.1549360602237604,0.1549360602237605,0.6901278795524791,0.4691507461438120,
     0.4691507461438120,0.0616985077123759,0.8392991722729236,0.1607008277270764,
     0.0000000000000000,0.0000000000000000,0.8392991722729236,0.1607008277270764};
    
    static const double lynw7[] =
    {0.0881563078002626,0.0020181692212805,0.0020181692212805,0.0020181692212805,
     0.0583260292949516,0.0583260292949516,0.0583260292949516,0.0532519959070870,
     0.0532519959070870,0.0532519959070870,0.0118425181549634,0.0118425181549634,
     0.0118425181549634,0.0118425181549634,0.0118425181549634,0.0118425181549634};
    
    
    // n = 8, rule = 15
    
    /*
    static const double lynr8[] =
    {0.3333333333333333,0.0814148234145541,0.4592925882927229,0.4592925882927230,
     0.8989055433659379,0.0505472283170310,0.0505472283170311,0.6588613844964797,
     0.1705693077517601,0.1705693077517602,0.0083947774099572,0.7284923929554041,
     0.2631128296346387,0.7284923929554041,0.2631128296346387,0.0083947774099572};
    
    static const double lyns8[] =
    {0.3333333333333333,0.4592925882927229,0.4592925882927230,0.0814148234145541,
     0.0505472283170310,0.0505472283170311,0.8989055433659379,0.1705693077517601,
     0.1705693077517602,0.6588613844964797,0.7284923929554041,0.2631128296346387,
     0.0083947774099572,0.0083947774099572,0.7284923929554041,0.2631128296346387};
    
    static const double lynw8[] =
    {0.0721578038388931,0.0475458171336425,0.0475458171336425,0.0475458171336425,
     0.0162292488115991,0.0162292488115991,0.0162292488115991,0.0516086852673592,
     0.0516086852673592,0.0516086852673592,0.0136151570872174,0.0136151570872174,
     0.0136151570872174,0.0136151570872174,0.0136151570872174,0.0136151570872174};
     */
    
    
    // n = 8, rule = 16
    /*
    static const double lynr8[] =
    {1.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,
     0.5000000000000000,0.5000000000000000,0.0086372116488837,0.4956813941755582,
     0.4956813941755582,0.8193444849714693,0.0903277575142653,0.0903277575142654,
     0.5316905005853895,0.2341547497073052,0.2341547497073053,0.0000000000000000,
     0.7236067977499790,0.2763932022500210,0.7236067977499790,0.2763932022500210,
     0.0000000000000000};
    
    static const double lyns8[] =
    {0.0000000000000000,0.0000000000000000,1.0000000000000000,0.5000000000000000,
     0.5000000000000000,0.0000000000000000,0.4956813941755582,0.4956813941755582,
     0.0086372116488837,0.0903277575142653,0.0903277575142654,0.8193444849714693,
     0.2341547497073052,0.2341547497073053,0.5316905005853895,0.7236067977499790,
     0.2763932022500210,0.0000000000000000,0.0000000000000000,0.7236067977499790,
     0.2763932022500210};
    
    static const double lynw8[] =
    {0.0020121232254880,0.0020121232254880,0.0020121232254880,-0.1415263313191909,
     -0.1415263313191909,-0.1415263313191909,0.1737279114818890,0.1737279114818890,
     0.1737279114818890,0.0324538298568710,0.0324538298568710,0.0324538298568710,
     0.0751975461200220,0.0751975461200220,0.0751975461200220,0.0124007936507937,
     0.0124007936507937,0.0124007936507937,0.0124007936507937,0.0124007936507937,
     0.0124007936507937};
    */ 
    
    
    // n = 8, rule = 17
    
    static const double lynr8[] =
    {0.3333333333333333,0.0466691212356951,0.4766654393821525,0.4766654393821524,
     0.9324563118910393,0.0337718440544803,0.0337718440544804,0.4593042216691921,
     0.2703478891654040,0.2703478891654040,0.0514643354866615,0.7458294907672514,
     0.2027061737460871,0.7458294907672514,0.2027061737460871,0.0514643354866615};
    
    static const double lyns8[] =
    {0.3333333333333333,0.4766654393821525,0.4766654393821524,0.0466691212356951,
     0.0337718440544803,0.0337718440544804,0.9324563118910393,0.2703478891654040,
     0.2703478891654040,0.4593042216691921,0.7458294907672514,0.2027061737460871,
     0.0514643354866615,0.0514643354866615,0.7458294907672514,0.2027061737460871};
    
    static const double lynw8[] =
    {-0.1417091925556979,0.0349534809663262,0.0349534809663262,0.0349534809663262,
     0.0085454563358004,0.0085454563358004,0.0085454563358004,0.1094149411652251,
     0.1094149411652251,0.1094149411652251,0.0304945928589405,0.0304945928589405,
     0.0304945928589405,0.0304945928589405,0.0304945928589405,0.0304945928589405};
    
    
    // n = 9, rule = 18
    
    /*
    static const double lynr9[] =
    {0.3333333333333333,0.0206349616025259,0.4896825191987370,0.4896825191987370,
     0.1258208170141290,0.4370895914929355,0.4370895914929354,0.6235929287619356,
     0.1882035356190322,0.1882035356190322,0.9105409732110941,0.0447295133944530,
     0.0447295133944530,0.0368384120547363,0.7411985987844980,0.2219629891607657,
     0.7411985987844980,0.2219629891607657,0.0368384120547363};
    
    static const double lyns9[] =
    {0.3333333333333333,0.4896825191987370,0.4896825191987370,0.0206349616025259,
     0.4370895914929355,0.4370895914929354,0.1258208170141290,0.1882035356190322,
     0.1882035356190322,0.6235929287619356,0.0447295133944530,0.0447295133944530,
     0.9105409732110941,0.7411985987844980,0.2219629891607657,0.0368384120547363,
     0.0368384120547363,0.7411985987844980,0.2219629891607657};
    
    static const double lynw9[] =
    {0.0485678981413981,0.0156673501135699,0.0156673501135699,0.0156673501135699,
     0.0389137705023877,0.0389137705023877,0.0389137705023877,0.0398238694636045,
     0.0398238694636045,0.0398238694636045,0.0127888378293491,0.0127888378293491,
     0.0127888378293491,0.0216417696886447,0.0216417696886447,0.0216417696886447,
     0.0216417696886447,0.0216417696886447,0.0216417696886447};
     */
    
    
    // n = 9, rule = 19
    
    static const double lynr9[] =
    {0.3333333333333333,1.0000000000000000,0.0000000000000000,0.0000000000000000,
     0.0000000000000000,0.5000000000000000,0.5000000000000000,0.1004413236259677,
     0.4497793381870162,0.4497793381870161,0.9061051136018193,0.0469474431990903,
     0.0469474431990903,0.6162561745251021,0.1918719127374489,0.1918719127374490,
     0.0368384120547363,0.7411985987844980,0.2219629891607657,0.7411985987844980,
     0.2219629891607657,0.0368384120547363};
    
    static const double lyns9[] =
    {0.3333333333333333,0.0000000000000000,0.0000000000000000,1.0000000000000000,
     0.5000000000000000,0.5000000000000000,0.0000000000000000,0.4497793381870162,
     0.4497793381870161,0.1004413236259677,0.0469474431990903,0.0469474431990903,
     0.9061051136018193,0.1918719127374489,0.1918719127374490,0.6162561745251021,
     0.7411985987844980,0.2219629891607657,0.0368384120547363,0.0368384120547363,
     0.7411985987844980,0.2219629891607657};
    
    static const double lynw9[] =
    {0.0566812422299596,0.0001770956316411,0.0001770956316411,0.0001770956316411,
     0.0080056858564321,0.0080056858564321,0.0080056858564321,0.0420707167722883,
     0.0420707167722883,0.0420707167722883,0.0130320906191451,0.0130320906191451,
     0.0130320906191451,0.0412037909998841,0.0412037909998841,0.0412037909998841,
     0.0216417696886447,0.0216417696886447,0.0216417696886447,0.0216417696886447,
     0.0216417696886447,0.0216417696886447};
    
    // n = 11, rule = 21
    
    static const double lynr11[] =
    {0.3333333333333333,0.9480217181434233,0.0259891409282883,0.0259891409282884,
     0.8114249947041546,0.0942875026479227,0.0942875026479227,0.0107264499655706,
     0.4946367750172147,0.4946367750172147,0.5853132347709715,0.2073433826145142,
     0.2073433826145143,0.1221843885990187,0.4389078057004907,0.4389078057004906,
     0.0000000000000000,0.8588702812826364,0.1411297187173636,0.8588702812826364,
     0.1411297187173636,0.0000000000000000,0.0448416775891306,0.6779376548825902,
     0.2772206675282792,0.6779376548825902,0.2772206675282792,0.0448416775891306};
    
    static const double lyns11[] =
    {0.3333333333333333,0.0259891409282883,0.0259891409282884,0.9480217181434233,
     0.0942875026479227,0.0942875026479227,0.8114249947041546,0.4946367750172147,
     0.4946367750172147,0.0107264499655706,0.2073433826145142,0.2073433826145143,
     0.5853132347709715,0.4389078057004907,0.4389078057004906,0.1221843885990187,
     0.8588702812826364,0.1411297187173636,0.0000000000000000,0.0000000000000000,
     0.8588702812826364,0.1411297187173636,0.6779376548825902,0.2772206675282792,
     0.0448416775891306,0.0448416775891306,0.6779376548825902,0.2772206675282792};
    
    static const double lynw11[] =
    {0.0439886505811109,0.0043721557768681,0.0043721557768681,0.0043721557768681,
     0.0190407859969677,0.0190407859969677,0.0190407859969677,0.0094277240280656,
     0.0094277240280656,0.0094277240280656,0.0360798487723705,0.0360798487723705,
     0.0360798487723705,0.0346645693527686,0.0346645693527686,0.0346645693527686,
     0.0036811918916503,0.0036811918916503,0.0036811918916503,0.0036811918916503,
     0.0036811918916503,0.0036811918916503,0.0205281577146443,0.0205281577146443,
     0.0205281577146443,0.0205281577146443,0.0205281577146443,0.0205281577146443};

    // Array of Lyness data.
    
    data.push_back(sIntData(lynr1,lyns1,0,lynw1,1));
    data.push_back(sIntData(lynr2,lyns2,0,lynw2,3));
    data.push_back(sIntData(lynr3,lyns3,0,lynw3,7));
    data.push_back(sIntData(lynr4,lyns4,0,lynw4,10));
    data.push_back(sIntData(lynr5,lyns5,0,lynw5,10));
    data.push_back(sIntData(lynr6,lyns6,0,lynw6,13));
    data.push_back(sIntData(lynr7,lyns7,0,lynw7,16));
    data.push_back(sIntData(lynr8,lyns8,0,lynw8,16));
    data.push_back(sIntData(lynr9,lyns9,0,lynw9,22));
    data.push_back(sIntData(lynr11,lyns11,0,lynw11,28));
    data.push_back(sIntData(lynr11,lyns11,0,lynw11,28));
  }

  return data;
}

// =============================== getTetData ==============================

const vector<sIntData>& getTetData( )
{
  static vector<sIntData> data;

  if (data.empty( )) 
  {
    static const double d11 = 1.0000000000000000;
    static const double d12 = 0.5000000000000000;
    static const double d14 = 0.2500000000000000;
    static const double d16 = 0.1666666666666667;

    // Degree = 1 (npoints = 1).
    
    static const double tet_r1[] = {d14};
    static const double tet_s1[] = {d14};
    static const double tet_t1[] = {d14};
    static const double tet_w1[] = {d16};
    
    // Degree = 2 (npoints = 4).
    
    static const double a = (5.0 + 3.0*sqrt(5.0))/20.0;
    static const double b = (5.0 - 1.0*sqrt(5.0))/20.0;
    static const double w2 = 1.0/24.0;
    static const double tet_r2[] = {a, b, b, b};
    static const double tet_s2[] = {b, b, b, a};
    static const double tet_t2[] = {b, b, a, b};
    static const double tet_w2[] = {w2, w2, w2, w2};

    // Degree = 3 (npoints = 5).
    
    static const double w3a = -4.0/30.0;
    static const double w3b =  9.0/120.0;
    static const double tet_r3[] = {d14, d12, d16, d16, d16};
    static const double tet_s3[] = {d14, d16, d16, d16, d12};
    static const double tet_t3[] = {d14, d16, d16, d12, d16};
    static const double tet_w3[] = {w3a, w3b, w3b, w3b, w3b};

    // Array of Tetrahedra data.
    
    data.push_back(sIntData(tet_r1, tet_s1, tet_t1, tet_w1, 1));
    data.push_back(sIntData(tet_r2, tet_s2, tet_t2, tet_w2, 4));
    data.push_back(sIntData(tet_r3, tet_s3, tet_t3, tet_w3, 5));
  }

  return data;
}


// Old Data


// Newton-Cotes.

const double NewtonCotesPoint[12][12] = {
  {.000000000000000,-1.000000000000000,-1.000000000000000,-1.000000000000000, -1.000000000000000, -1.000000000000000,-1.000000000000000,-1.000000000000000,-1.000000000000000,-1.000000000000000, -1.000000000000000, -1.0000000000000000 },
  {.000000000000000, 1.000000000000000, 0.000000000000000,-0.333333333333333, -0.500000000000000, -0.600000000000000,-0.666666666666666,-0.714285714285714,-0.750000000000000,-0.777777777777777, -0.800000000000000, -0.8181818181818181 },
  {.000000000000000, 0.000000000000000, 1.000000000000000, 0.333333333333333,  0.000000000000000, -0.200000000000000,-0.333333333333333,-0.428571428571429,-0.500000000000000,-0.555555555555555, -0.600000000000000, -0.6363636363636363 },
  {.000000000000000, 0.000000000000000, 0.000000000000000, 1.000000000000000,  0.500000000000000,  0.200000000000000, 0.000000000000000,-0.142857142857143,-0.250000000000000,-0.333333333333333, -0.400000000000000, -0.4545454545454545 },
  {.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000,  1.000000000000000,  0.600000000000000, 0.333333333333333, 0.142857142857143, 0.000000000000000,-0.111111111111111, -0.200000000000000, -0.2727272727272727 },
  {.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000,  0.000000000000000,  1.000000000000000, 0.666666666666666, 0.428571428571429, 0.250000000000000, 0.111111111111111, 0.000000000000000,  -0.0909090909090909 },
  {.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000,  0.000000000000000,  0.000000000000000, 1.000000000000000, 0.714285714285714, 0.750000000000000, 0.333333333333333,  0.200000000000000,  0.0909090909090909 },
  {0.00000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000,  0.000000000000000,  0.000000000000000, 0.000000000000000, 1.000000000000000, 1.000000000000000, 0.555555555555555,  0.400000000000000,  0.2727272727272727 },
  {0.00000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000,  0.000000000000000,  0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000, 0.777777777777777,  0.600000000000000,  0.4545454545454545 },
  {0.00000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000,  0.000000000000000,  0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000, 1.000000000000000,  0.800000000000000,  0.6363636363636363 },
  {.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000,  0.000000000000000,  0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000,  1.000000000000000,  0.8181818181818181 },
  {.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000,  0.000000000000000,  0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000,  0.000000000000000,  1.0000000000000000 },
 };

const double NewtonCotesWeight[12][12] = {
  {4.00000000000000, 1.000000000000000, 0.333333333333333, 0.250000000000000, 0.155555555555555,  0.131944444444444,   0.097619047619048, 0.086921296296296, 0.069770723104057, 0.063772321428572, 0.053668296723852,  0.0471753363864754 },
  {0.00000000000000, 1.000000000000000, 1.333333333333333, 0.750000000000000, 0.711111111111111,  0.520833333333333,   0.514285714285714, 0.414004629629630, 0.415379188712522, 0.351361607142857, 0.355071882849661,  0.1069393259953637 },
  {0.00000000000000, 0.000000000000000, 0.333333333333333, 0.750000000000000, 0.266666666666666,  0.347222222222222,   0.064285714285714, 0.153125000000000,-0.065467372134039, 0.024107142857143,-0.162087141253808,  0.1600783285433586 },
  {0.00000000000000, 0.000000000000000, 0.000000000000000, 0.250000000000000, 0.711111111111111,  0.347222222222222,   0.647619047619048, 0.345949074074074, 0.740458553791887, 0.431785714285714, 0.909892576559243,  0.2031674267230672 },
  {.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000, 0.155555555555555,  0.520833333333333,   0.064285714285714, 0.345949074074074,-0.320282186948854, 0.128973214285714,-0.870310245310245,  0.2334925365383534 },
  {.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000,  0.131944444444444,   0.514285714285714, 0.153125000000000, 0.740458553791887, 0.128973214285714, 1.427529260862590,  0.2491470458134027 },
  {.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000,  0.000000000000000,   0.097619047619048, 0.414004629629630,-0.065467372134039, 0.431785714285714,-0.870310245310245,  0.2491470458134027 },
  {0.00000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000,  0.000000000000000,   0.000000000000000, 0.086921296296296, 0.415379188712522, 0.024107142857143, 0.909892576559243,  0.2334925365383534 },
  {0.00000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000,  0.000000000000000,   0.000000000000000, 0.000000000000000, 0.069770723104057, 0.351361607142857,-0.162087141253808,  0.2031674267230672 },
  {0.00000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000,  0.000000000000000,   0.000000000000000, 0.000000000000000, 0.000000000000000, 0.063772321428572, 0.355071882849661,  0.1600783285433586 },
  {.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000,  0.000000000000000,   0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000, 0.053668296723852,  0.1069393259953637 },
  {.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000,  0.000000000000000,   0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000,  0.0471753363864754 },
};
// ======================================================= End of file =====
