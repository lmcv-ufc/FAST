// ------------------------------------------------------------------------
// shpplane.h - file containing the definition of the Plane Shape class.
// ------------------------------------------------------------------------
//
// The numbering schemes of cShapePlane subclasses are given below:
//
// ------------------------------------------------------------------------
//
//                     SHAPE_T3               SHAPE_T6
//  natural
//  coord.                  3 o                   5 o
//                           /|                    /|
//   s                      / |                   / |
//   ^                     /  |                6 o  o 4
//   |                    /   |                 /   |
//   |                   /    |                /    |
//   -----> r         1 o-----o 2           1 o--o--o 3
//                                               2
// ------------------------------------------------------------------------
//
//                      SHAPE_Q4                SHAPE_Q8
//
//   natural                                        6
//   coord.        4 o-----------o 3        7 o-----o-----o 5
//                   |           |            |           |
//    s              |           |            |           |
//    ^              |           |          8 o           o 4
//    |              |           |            |           |
//    |              |           |            |           |
//    -----> r     1 o-----------o 2        1 o-----o-----o 3
//                                                  2
//
//                      SHAPE_Q9
//
//                         6
//                 7 o-----o-----o 5
//                   |           |
//                   |     9     |
//                 8 o     o     o 4
//                   |           |
//                   |           |
//                 1 o-----o-----o 3
//                         2
//
// ------------------------------------------------------------------------

#ifndef _SHPPLANE_H
#define _SHPPLANE_H

#include "shape.h"

// ------------------------------------------------------------------------
// Definition of the Plane Shape class:
//
class cShapePlane : public cShape
{
 public:
                 cShapePlane(void);
  virtual       ~cShapePlane(void);
          int    NumDim(void) { return 2; }
  virtual void   NodeNatCoord(sNatCoord *) = 0;
  virtual void   ShpFunc(sNatCoord, double *) = 0;
  virtual void   DrvShpRST(sNatCoord, sNatCoord *) = 0;
  virtual void   DrvShpRST(sNatCoord, sNatCoord *, sNatDrv *) { }
  virtual void   DrvShpXYZ(sNatCoord, sNodeCoord *, double *,
                            sNodeCoord *);
  virtual void   DrvShpXYZ(sNatCoord, sNodeCoord *, double *,
                           sNodeCoord *, sNodeDrv *);
  virtual int    GetEdge(int *, eShpType *, int *, cNode **) = 0;
  virtual int    GetFace(int *, eShpType *, int *, cNode **);
};


// ------------------------------------------------------------------------
// Definition of the T3 shape class:
//
class cShapeT3 : public virtual cShapePlane
{
 public:
                cShapeT3(void);
  virtual      ~cShapeT3(void);
          void  NodeNatCoord(sNatCoord *);
          void  ShpFunc(sNatCoord, double *);
          void  DrvShpRST(sNatCoord, sNatCoord *);
          int   GetEdge(int *, eShpType *, int *, cNode **);
};


// ------------------------------------------------------------------------
// Definition of the T6 shape class:
//
class cShapeT6 : public virtual cShapePlane
{
 public:
                 cShapeT6(void);
  virtual       ~cShapeT6(void);
          void   NodeNatCoord(sNatCoord *);
          void   ShpFunc(sNatCoord, double *);
          void   DrvShpRST(sNatCoord, sNatCoord *);
          int    GetEdge(int *, eShpType *, int *, cNode **);
};


// ------------------------------------------------------------------------
// Definition of the Q4 shape class:
//
class cShapeQ4 : public virtual cShapePlane
{
 public:
                cShapeQ4(void);
  virtual      ~cShapeQ4(void);
          void  NodeNatCoord(sNatCoord *);
          void  ShpFunc(sNatCoord, double *);
          void  DrvShpRST(sNatCoord, sNatCoord *);
          int   GetEdge(int *, eShpType *, int *, cNode **);
};


// ------------------------------------------------------------------------
// Definition of the Q8 shape class:
//
class cShapeQ8 : public virtual cShapePlane
{
 public:
                cShapeQ8(void);
  virtual      ~cShapeQ8(void);
          void  NodeNatCoord(sNatCoord *);
          void  ShpFunc(sNatCoord, double *);
          void  DrvShpRST(sNatCoord, sNatCoord *);
          int   GetEdge(int *, eShpType *, int *, cNode **);
};


// ------------------------------------------------------------------------
// Definition of the Q9 shape class:
//
class cShapeQ9 : public virtual cShapePlane
{
 public:
                cShapeQ9(void);
  virtual      ~cShapeQ9(void);
          void  NodeNatCoord(sNatCoord *);
          void  ShpFunc(sNatCoord, double *);
          void  DrvShpRST(sNatCoord, sNatCoord *);
          int   GetEdge(int *, eShpType *, int *, cNode **);
};

#endif
