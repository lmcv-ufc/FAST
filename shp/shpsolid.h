// ------------------------------------------------------------------------
// shpsolid.h - file containing the definition of the Solid Shape class.
// ------------------------------------------------------------------------
//
// The numbering schemes of cShapeSolid subclasses are given below:
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
//                               SHAPE_BRICK8
//
//         8 o---------------o 7                8 o---------------o 7
//          /               /|                   /|              /|
//         /               / |                  / |             / |
//        /               /  |                 /  |            /  |
//       /               /   |                /   |           /   |
//    4 /             3 /    |             4 /    |        3 /    |
//     o---------------o     |              o---------------o     |
//     |               |     |              |     |         |     |
//     |               |     o              |  5  o---------|-----o 6
//     |               |    /               |    /          |    /
//     |               |   /                |   /           |   /
//     |               |  /                 |  /            |  /
//     |               | /                  | /             | /
//     |               |/                   |/              |/
//     o---------------o                    o---------------o
//    1                 2                  1                 2
//
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
//                                SHAPE_BRICK20
//
//                    18                                   18
//         19 o-------o-------o 17              19 o-------o-------o 17
//           /               /|                   /|              /|
//          /               / |                  / |             / |
//      12 o            11 o  |              12 o  |20       11 o  |
//        /               /   o 16             /   o           /   o 16
//       /       6       /    |               /    |  6       /    |
//    7 o-------o-------o5    |            7 o-------o-------o5    |
//      |               |     |              |     |         |     |
//      |               |     o 15           |   13o------14-|-----o 15
//      |               |    /               |    /          |    /
//    8 o               o4  /              8 o   /           o4  /
//      |               |  o 10              |  o 9          |  o 10
//      |               | /                  | /             | /
//      |               |/                   |/              |/
//    1 o-------o-------o 3                1 o-------o-------o 3
//              2                                    2
// ------------------------------------------------------------------------
//
//                             SHAPE_TET4
//
//      The nodes are numbered according to the following rule:
//      1.  first, counterclockwise, the nodes at the face with
//          natural coordinate L4 = 0,
//      2.  and finally the node at the vertex with L4 = 1.
//
//
//                                 o 4 (L4 = 1)
//                                /|\
//                               / | \
//                              /  |  \
//                             /   |   \
//                            /    |    \
//                           /     |     \
//                          /      |      \
//                         /       |       \
//                        /        |        \
//                       /         |         \
//                      /          |          \
//          (L1 = 1) 1 o . . . . . | . . . . . o 2 (L2 = 1)
//                      \_         |         _/
//                        \_       |       _/
//                          \_     |     _/
//                            \_   |   _/
//                              \_ | _/
//                                \|/
//                                 o 3 (L3 = 1)
//
// ------------------------------------------------------------------------
//                             SHAPE_TET10
//
//        The nodes are numbered according to the following rule:
//      1.  first, counterclockwise, the nodes at the face with
//          natural coordinate L4 = 0,
//      2.  then, also counterclockwise, the nodes at the middle
//          of the edges that go from the node with L4 = 1 to the
//          the face with L4 = 0,
//      3.  and finally the node at the vertex with L4 = 1.
//
//                                 o 10 (L4 = 1)
//                                /|\
//                               / | \
//                              /  |  \
//                             /   |   \
//                            /    |    \
//                         7 o     |     o 8
//                          /      |      \
//                         /       |       \
//                        /        o 9      \
//                       /         |         \
//                      /          |          \
//          (L1 = 1) 1 o . . . . . |o 2. . . . o 3 (L2 = 1)
//                      \_         |         _/
//                        \_       |       _/
//                        6 o_     |     _o 4
//                            \_   |   _/
//                              \_ | _/
//                                \|/
//                                 o 5 (L3 = 1)
//
// ------------------------------------------------------------------------

#ifndef _SHPSOLID_H
#define _SHPSOLID_H

#include "shape.h"

// ------------------------------------------------------------------------
// Definition of the Solid Shape class:
//
class cShapeSolid : public cShape
{
 public:
                 cShapeSolid(void);
  virtual       ~cShapeSolid(void);
          int    NumDim(void) { return 3; }
  virtual void   NodeNatCoord(sNatCoord *) = 0;
  virtual void   ShpFunc(sNatCoord, double *) = 0;
  virtual void   DrvShpRST(sNatCoord, sNatCoord *) = 0;
          void   DrvShpXYZ(sNatCoord, sNodeCoord *, double *,
                            sNodeCoord *);
  virtual void   tVector(sNatCoord, sNodeCoord *, cVector &);
  virtual int    GetEdge(int *, eShpType *, int *, cNode **) = 0;
  virtual int    GetFace(int *, eShpType *, int *, cNode **) { return 0; }
};


// ------------------------------------------------------------------------
// Definition of the Brick8 shape class:
//
class cShapeBrick8 : public cShapeSolid
{
 public:
         cShapeBrick8(void);
        ~cShapeBrick8(void);
  void   NodeNatCoord(sNatCoord *);
  void   ShpFunc(sNatCoord, double *);
  void   DrvShpRST(sNatCoord, sNatCoord *);
  int    GetEdge(int *, eShpType *, int *, cNode **);
  int    GetFace(int *, eShpType *, int *, cNode **);
};


// ------------------------------------------------------------------------
// Definition of the Brick20 shape class:
//
class cShapeBrick20 : public cShapeSolid
{
 public:
         cShapeBrick20(void);
        ~cShapeBrick20(void);
  void   NodeNatCoord(sNatCoord *);
  void   ShpFunc(sNatCoord, double *);
  void   DrvShpRST(sNatCoord, sNatCoord *);
  int    GetEdge(int *, eShpType *, int *, cNode **);
  int    GetFace(int *, eShpType *, int *, cNode **);
};


// ------------------------------------------------------------------------
// Definition of the Tet4 shape class:
//
class cShapeTet4 : public cShapeSolid
{
 public:
         cShapeTet4(void);
        ~cShapeTet4(void);
  void   NodeNatCoord(sNatCoord *);
  void   ShpFunc(sNatCoord, double *);
  void   DrvShpRST(sNatCoord, sNatCoord *);
  int    GetEdge(int *, eShpType *, int *, cNode **);
  int    GetFace(int *, eShpType *, int *, cNode **);
};


// ------------------------------------------------------------------------
// Definition of the Tet10 shape class:
//
class cShapeTet10 : public cShapeSolid
{
 public:
         cShapeTet10(void);
        ~cShapeTet10(void);
  void   NodeNatCoord(sNatCoord *);
  void   ShpFunc(sNatCoord, double *);
  void   DrvShpRST(sNatCoord, sNatCoord *);
  int    GetEdge(int *, eShpType *, int *, cNode **);
  int    GetFace(int *, eShpType *, int *, cNode **);
};


#endif
