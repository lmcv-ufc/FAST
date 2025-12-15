// ------------------------------------------------------------------------
// shpinterf.h - file containing the definition of interface element
// shape class.
// ------------------------------------------------------------------------
// The numbering scheme of two-dimensional interface element shapes are:
// ------------------------------------------------------------------------
//
//                    s
//                    ^
//                    |    (natural coord. system)
//                    |
//                    -----> r
//
//            1 o-----------o 2
//
//            3 o-----------o 4
//                                 INTERF_L2
// -------------------------------------------------------------------------
//
//                    s
//                    ^
//                    |    (natural coord. system)
//                    |
//                    -----> r
//
//                    2
//            1 o-----o-----o 3
//
//            4 o-----o-----o 6
//                    5            INTERF_L3
//
// -------------------------------------------------------------------------
// The numbering scheme of three-dimensional interface element shapes are:
// ------------------------------------------------------------------------
//
//         s
//        ^
//        |    natural coord
//        |
//         -----> r
//       /
//      /
//     t
//                               SHAPE_INTERF_Q4
//
//         8 o---------------o 7                8 o---------------o 7
//          .|              .|                   .|              .|
//         . |             . |                  . |             . |
//        .  |            .  |                 .  |            .  |
//       .   |           .   |                .   |           .   |
//    4 .    |        3 .    |             4 .    |        3 .    |
//     o---------------o     |              o---------------o     |
//     |               |     |              |     |         |     |
//     |               |     o 6            |  5  o---------|-----o 6
//     |               |    .               |    .          |    .
//     |               |   .                |   .           |   .
//     |               |  .                 |  .            |  .
//     |               | .                  | .             | .
//     |               |.                   |.              |.
//     o---------------o                    o---------------o
//    1                 2                  1                 2
//
// ------------------------------------------------------------------------
//
//         s
//        ^
//        |    natural coord
//        |
//         -----> r
//       /
//      /
//     t
//                                SHAPE_INTERF_Q8
//
//                    14                                   14
//         15 o-------o-------o 13              15 o-------o-------o 13
//           .|              .|                   .|              .|
//          . |             . |                  . |             . |
//         .  |16          .  |                 .  |16          .  |
//        .   o           .   o 12             .   o           .   o 12
//       .    |  6       .    |               .    |  6       .    |
//    7 o-------o-------o5    |            7 o-------o-------o5    |
//      |               |     |              |     |         |     |
//      |               |     o 11           |    9o------10-|-----o 11
//      |               |    .               |    .          |    .
//    8 o               o4  .              8 o   .           o4  .
//      |               |  .                 |  .            |  .   
//      |               | .                  | .             | .
//      |               |.                   |.              |.
//    1 o-------o-------o 3                1 o-------o-------o 3
//              2                                    2
// ------------------------------------------------------------------------

#ifndef _SHPINTERF_H
#define _SHPINTERF_H

#include "shape.h"

// ------------------------------------------------------------------------
// Definition of the 2D Interface Shape class:
//
class cShapeInterf2D : public cShape
{
 public:
                 cShapeInterf2D(void);
  virtual       ~cShapeInterf2D(void);
          int    NumDim(void) { return 1; }
  virtual void   NodeNatCoord(sNatCoord *) = 0;
  virtual void   ShpFunc(sNatCoord, double *) = 0;
  virtual void   DrvShpRST(sNatCoord, sNatCoord *) = 0;
  virtual void   DrvShpXYZ(sNatCoord, sNodeCoord *, double *,
                               sNodeCoord *) { }
  virtual void   RMatrix(sNatCoord, sNodeCoord *, double *, cMatrix &);
  virtual int    GetEdge(int *, eShpType *, int *, cNode **) { return 0; }
  virtual int    GetFace(int *, eShpType *, int *, cNode **) { return 0; }
};

// ------------------------------------------------------------------------
// Definition of the L2 Interface shape class:
//
//
class cShapeInterfL2 : public virtual cShapeInterf2D
{
 public:
                cShapeInterfL2(void);
  virtual      ~cShapeInterfL2(void);
          void  NodeNatCoord(sNatCoord *);
          void  ShpFunc(sNatCoord, double *);
          void  DrvShpRST(sNatCoord, sNatCoord *);
};

// ------------------------------------------------------------------------
// Definition of the L3 Interface shape class:
//
//
class cShapeInterfL3 : public virtual cShapeInterf2D
{
 public:
                cShapeInterfL3(void);
  virtual      ~cShapeInterfL3(void);
          void  NodeNatCoord(sNatCoord *);
          void  ShpFunc(sNatCoord, double *);
          void  DrvShpRST(sNatCoord, sNatCoord *);
};


// ------------------------------------------------------------------------
// Definition of the 3D Interface Shape class:
//
class cShapeInterf3D : public cShape
{
 public:
                 cShapeInterf3D(void);
  virtual       ~cShapeInterf3D(void);
          int    NumDim(void) { return 2; }
  virtual void   NodeNatCoord(sNatCoord *) = 0;
  virtual void   ShpFunc(sNatCoord, double *) = 0;
  virtual void   DrvShpRST(sNatCoord, sNatCoord *) = 0;
  virtual void   DrvShpXYZ(sNatCoord, sNodeCoord *, double *,
                               sNodeCoord *) { }
  virtual void   RMatrix(sNatCoord, sNodeCoord *, double *, cMatrix &);
  virtual int    GetEdge(int *, eShpType *, int *, cNode **) { return 0; }
  virtual int    GetFace(int *, eShpType *, int *, cNode **) { return 0; }
};

// ------------------------------------------------------------------------
// Definition of the Q4 Interface shape class:
//
//
class cShapeInterfQ4 : public virtual cShapeInterf3D
{
 public:
                cShapeInterfQ4(void);
  virtual      ~cShapeInterfQ4(void);
          void  NodeNatCoord(sNatCoord *);
          void  ShpFunc(sNatCoord, double *);
          void  DrvShpRST(sNatCoord, sNatCoord *);
};

// ------------------------------------------------------------------------
// Definition of the Q8 Interface shape class:
//
//
class cShapeInterfQ8 : public virtual cShapeInterf3D
{
 public:
                cShapeInterfQ8(void);
  virtual      ~cShapeInterfQ8(void);
          void  NodeNatCoord(sNatCoord *);
          void  ShpFunc(sNatCoord, double *);
          void  DrvShpRST(sNatCoord, sNatCoord *);
};

#endif
