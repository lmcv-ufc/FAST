// -------------------------------------------------------------------------
// shpbsp.h - file containing the definition of BSpline Shape class.
// -------------------------------------------------------------------------

#ifndef _SHPBSP_H
#define _SHPBSP_H

#include "shpplane.h"
#include "shpsurf.h"
#include "shpsolid.h"
#include "shpiga.h"

// ------------------------------------------------------------------------
// Foward declarations:
//
class cBspPat;

// ------------------------------------------------------------------------
// Isogeometric B-Spline shape class:
//
class cShapeBsp : public cShapeIGA
{
 protected:
  cBspPat      *pat;
  int           NumBas[3];
  int           SpanInd[3];

 public:
                cShapeBsp(void);
  virtual      ~cShapeBsp(void);


  void          PatElmMapAdd(int);
  cPatch*       GetPatch(void) {return pat;}
  int           GetPatSpanID(void);
  void          ReadIGA(void);
  void          SetTopology(cPatch*, int *);
  void          EvalBspBas(int,double,cVector&); 
  void          EvalBspBasDerv(int,double,cVector&);
  void          CompBasisIndex(const int&,int&,int*&);
};

// ------------------------------------------------------------------------
// Two-dimensional isogeometric curve shape class:
//
class cShapePlBspCurve : public virtual cShapeLine2D,
                         public virtual cShapeBsp
{
 public:
                cShapePlBspCurve(void);
  virtual      ~cShapePlBspCurve(void);
          void  Read(void) { ReadIGA( ); }
          void  Init(void);
  virtual void  NodeNatCoord(sNatCoord *);
  virtual void  ShpFunc(sNatCoord, double *);
  virtual void  DrvShpRST(sNatCoord, sNatCoord *);
          int   GetEdge(int *, eShpType *, int *, cNode **) {return 0;}
          //cPatch* GetPatEdge(int,int,eShpType&,int[]);
};

// -------------------------------------------------------------------------
// Three-dimensional isogeometric B-Spline curve Shape class:
//
class cShapeBspCurve : public virtual cShapePlBspCurve,
                       public virtual cShapeLine3D
{
 public:
           cShapeBspCurve(void);
  virtual ~cShapeBspCurve(void);
};


// ------------------------------------------------------------------------
// Definition of the Plane Surface IGA shape class:
//
class cShapePlBspSurf : public virtual cShapePlane,
                        public virtual cShapeBsp
{
 public:
                cShapePlBspSurf(void);
  virtual      ~cShapePlBspSurf(void);

          void  Init(void);
          void  Read(void) { ReadIGA( ); }
          void  NodeNatCoord(sNatCoord *);
          void  ShpFunc(sNatCoord, double *);
          void  DrvShpRST(sNatCoord, sNatCoord *);
          void  DrvShpRST(sNatCoord, sNatCoord *, sNatDrv *);
          cPatch* GetPatEdge(int,int,eShpType&,int[]);
          cPatch* GetPatFace(int,eShpType&,int[]);
          int   GetEdge(int*,eShpType*,int*,cNode**) {return 0;}
};


// -------------------------------------------------------------------------
// Definition of the Surface IGA shape class:
//
class cShapeBspSurf : public virtual cShapeSurf,
                      public virtual cShapePlBspSurf
{
 public:
              cShapeBspSurf(void);
  virtual    ~cShapeBspSurf(void);
};

// ------------------------------------------------------------------------
// Definition of the Solid IGA Shape class:
//
class cShapeBspSolid : public cShapeSolid,
                       public virtual cShapeBsp
{
 public:
         cShapeBspSolid(void);
        ~cShapeBspSolid(void);
  void   Init(void);
  void  Read(void) { ReadIGA( ); }
  void   NodeNatCoord(sNatCoord *);
  void   ShpFunc(sNatCoord, double *);
  void   DrvShpRST(sNatCoord, sNatCoord *);
//         cPatch* GetPatEdge(int,int,eShpType&,int[]) {return 0;}
         cPatch* GetPatFace(int,eShpType&,int[]);
  int   GetEdge(int*,eShpType*,int*,cNode**) {return 0;}
  int   GetFace(int*,eShpType*,int*,cNode**) {return 0;}

};

#endif
