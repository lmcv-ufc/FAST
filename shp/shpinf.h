// -------------------------------------------------------------------------
// shpinf.h - this file contains methods for the 2D mapped infinite
//             elements shape.
// -------------------------------------------------------------------------
// Remarks:
// -------------------------------------------------------------------------
// The numbering scheme is:
//
//                    s
//                    ^
//                    |    (natural coord. system)
//                    |
//                    -----> r
//
//                     s                                    s
//                     ^                                    ^
//          |          |                          |         |          |
//          |          |                          |         |          |
//          |          |                          |         |          |
//        4 +          +-----> r    or          6 +         +-----> r  + 4
//          |          3                          |         5          |
//          |                                     |                    |
//          |                                     |                    |
//        1 + ---------+----------              1 +---------+----------+
//                     2                                    2          3
//
//            INFINITE 4-noded                      INFINITE 6-noded
//
//
//  The 4-noded and 6-noded Lagrangian isoparametric elements are two
//  dimensional mapped infinite elements.
//
//  The 6-noded Lagrangian isoparametric element is a single infinite element.
//  This term is used to denote elements which extend to infinity in one
//  direction only. This element may be regarded as a 9-noded Lagrangian
//  element in which the three nodes associated with s = +1 are positioned at
//  infinity.
//
//  The 4-noded is a doubly infinite element, which means that this element
//  extends to infinite in two directions.
// -------------------------------------------------------------------------

#ifndef _SHPINF_H
#define _SHPINF_H

#include "shpplane.h"

// -------------------------------------------------------------------------
// Infinite (Finite) Element Shape class:
//
class cShapePlaneInf : public cShapePlane
{
 public:
                 cShapePlaneInf(void);
   virtual      ~cShapePlaneInf(void);
   virtual bool  Isoparametric (void) { return(false); }
   virtual void  MapFunc       (sNatCoord, double *) = 0;
   virtual void  ShpFunc       (sNatCoord, double *) = 0;
   virtual void  DrvMapRST     (sNatCoord, sNatCoord *) = 0;
   virtual void  DrvShpRST     (sNatCoord, sNatCoord *) = 0;
   virtual int   GetEdge       (int *, eShpType *, int *, cNode **) = 0;

};


// ------------------------------------------------------------------------
// Definition of the L6Inf shape class:
//
class cShapeL6Inf : public virtual cShapePlaneInf
{
 public:
                cShapeL6Inf (void);
  virtual      ~cShapeL6Inf (void);
          void  NodeNatCoord(sNatCoord *);
          void  MapFunc     (sNatCoord, double *);
          void  ShpFunc     (sNatCoord, double *);
          void  DrvMapRST   (sNatCoord, sNatCoord *);
          void  DrvShpRST   (sNatCoord, sNatCoord *);
          int   GetEdge     (int *, eShpType *, int *, cNode **);
};

/*
// ------------------------------------------------------------------------
// Definition of the L4Inf shape class:
//
class cShapeL4Inf : public virtual cShapePlaneInf
{
 public:
                cShapeL4Inf (void);
  virtual      ~cShapeL4Inf (void);
          void  NodeNatCoord(sNatCoord *);
          void  MapFunc     (sNatCoord, double *);
          void  ShpFunc     (sNatCoord, double *);
          void  DrvMapRST   (sNatCoord, sNatCoord *);
          void  DrvShpRST   (sNatCoord, sNatCoord *);
          int   GetEdge     (int *, eShpType *, int *, cNode **);
};
*/

#endif
