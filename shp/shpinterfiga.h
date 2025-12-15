// ------------------------------------------------------------------------
// shpinterfiga.h - file containing the definition of interface IGA element
// shape class.
// ------------------------------------------------------------------------
// The numbering scheme of two-dimensional interface element shapes are:
// ------------------------------------------------------------------------
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
//                    5          SHAPE_IGA_INTERF_LINE
//
// -------------------------------------------------------------------------
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

#ifndef _SHPINTERFIGA_H
#define _SHPINTERFIGA_H

#include "patch.h"
#include "shpinterf.h"
#include "shpbsp.h"
#include "shpline.h"
#include "shpplane.h"
#include "shpsurf.h"
//#include "shpsolid.h"

using namespace std;

#include <iosfwd>


/*
// ------------------------------------------------------------------------
// Definition of the 2D Interface Shape class:

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


class cShapeInterfL2 : public virtual cShapeInterf2D
{
 public:
                cShapeInterfL2(void);
  virtual      ~cShapeInterfL2(void);
          void  NodeNatCoord(sNatCoord *);
          void  ShpFunc(sNatCoord, double *);
          void  DrvShpRST(sNatCoord, sNatCoord *);
};
*/
// ------------------------------------------------------------------------
// Definition of the IGA Line Interface shape class:


class cShapeInterflineBsp : public cShapeInterf2D,
                            public virtual cShapeBsp
			    // public virtual cShapePlCurveIGA
{
 protected:
  int ActParDir;    // Active parametric direction.


 public:
                cShapeInterflineBsp(void);
  virtual      ~cShapeInterflineBsp(void);

          void  Init(void);
          void  Read(void);
          void  NodeNatCoord(sNatCoord *);
          void  ShpFunc(sNatCoord, double *);
          void  DrvShpRST(sNatCoord, sNatCoord *);
       //   int   GetNumNodeNodStr(void);
       //  cPatch* GetPatEdge(int, int, eShpType&, int[]);
       //  cPatch* GetPatFace(int, eShpType&, int[]);
       //   int   GetEdge(int *, eShpType *, int *, cNode **) {return 0;}
};

// ------------------------------------------------------------------------
// Definition of the L3 Interface shape class:
//
//
//class cShapeInterfL3 : public virtual cShapeInterf2D
//{
// public:
//                cShapeInterfL3(void);
//  virtual      ~cShapeInterfL3(void);
//          void  NodeNatCoord(sNatCoord *);
//          void  ShpFunc(sNatCoord, double *);
//          void  DrvShpRST(sNatCoord, sNatCoord *);
//};

#endif
