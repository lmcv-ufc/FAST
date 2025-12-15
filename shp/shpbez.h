// -------------------------------------------------------------------------
// shapeiga.h - file containing the definition of the Shape class.
// -------------------------------------------------------------------------

#ifndef _SHPBEZ_H
#define _SHPBEZ_H

#include "shpline.h"
#include "shpplane.h"
#include "shpsurf.h"

using namespace std;

// ------------------------------------------------------------------------
// Definition of the plane bezier curve shape class:
//
class cShapePlBezCurve : public virtual cShapeLine2D
{
 protected:
  int  Degree;
  bool Rational;

 public:
                cShapePlBezCurve(void);
  virtual      ~cShapePlBezCurve(void);

          void  Read(void);
          void  SetNodes(cNode **, int size);
          void  NodeNatCoord(sNatCoord *);
          void  ShpFunc(sNatCoord, double *);
          void  DrvShpRST(sNatCoord, sNatCoord *);
          int   GetEdge(int *, eShpType *, int *, cNode **);
};

// -------------------------------------------------------------------------
// Three-dimensional isogeometric B-Spline curve Shape class:
//
class cShapeBezCurve : public virtual cShapePlBezCurve,
                       public virtual cShapeLine3D
{
 public:
           cShapeBezCurve(void);
  virtual ~cShapeBezCurve(void);
};

// ------------------------------------------------------------------------
// Definition of the Plane Triangular Bézier Surface IGA shape class:
//
class cShapePlBezTrianSurf : public virtual cShapePlane
{
 protected:
  int  Degree;
  bool Rational;

 public:
                cShapePlBezTrianSurf(void);
  virtual      ~cShapePlBezTrianSurf(void);

          void  Read(void);
          //void  SetNodes(cNode **, int size);
          void  NodeNatCoord(sNatCoord *);
          void  ShpFunc(sNatCoord, double *);
          void  DrvShpRST(sNatCoord, sNatCoord *);
          int   GetEdge(int *, eShpType *, int *, cNode **);
};

// -------------------------------------------------------------------------
// Definition of the Bezier Triangle Surface shape class:
//
class cShapeBezTrianSurf : public virtual cShapeSurf,
                           public virtual cShapePlBezTrianSurf
{
 public:
              cShapeBezTrianSurf(void);
  virtual    ~cShapeBezTrianSurf(void);
};

// ------------------------------------------------------------------------
// Definition of the Plane Tensor product Bézier Surface IGA shape class:
//
class cShapePlBezSurf : public virtual cShapePlane
{
 protected:
  int  Degree[2];
  bool Rational;

 public:
                cShapePlBezSurf(void);
  virtual      ~cShapePlBezSurf(void);

          void  Read(void);
          //void  SetNodes(cNode **, int size);
          void  NodeNatCoord(sNatCoord *);
          void  ShpFunc(sNatCoord, double *);
          void  DrvShpRST(sNatCoord, sNatCoord *);
          int   GetEdge(int *, eShpType *, int *, cNode **);
};

// -------------------------------------------------------------------------
// Definition of the Bezier Surface shape class:
//
class cShapeBezSurf : public virtual cShapeSurf,
                      public virtual cShapePlBezSurf
{
 public:
              cShapeBezSurf(void);
  virtual    ~cShapeBezSurf(void);
};

#endif
