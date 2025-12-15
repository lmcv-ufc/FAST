// -------------------------------------------------------------------------
// shpline.h - file containing the definition of the Line Shape class.
// -------------------------------------------------------------------------
// Remarks:
// -------------------------------------------------------------------------
//
// There are two line shape sub-classes:
//
//   cShapeLine2D
//     Planar version for line shapes lying in the XY plane.
//     The Jacobian matrix is computed as follows:
//      *The first row is formed by the derivatives of the position vector
//       with respect to the line parametric coordinate, i.e., the line
//       tangent vector.
//      *The second row is a direction in the plane which is normal to the
//       line tangent direction, i.e., the line normal vector.
//     The LocalSys matrix is [R] = [{t} {n} {ez}], where {t} is the unit
//     tangent vector, {n} is the unit normal vector and {ez} is the unit
//     z-vector.
//
//    cShapeLine3D
//      General, curved line shapes in space.
//      The Jacobian matrix is computed as follows:
//      *The first row is formed by the derivatives of the position vector
//       with respect to the line parametric coordinate, i.e., the line
//       tangent vector
//      *The second row is any vector which is found perdendicular to the
//       tangent vector. The first choice is a direction which also
//       perpendicular to the z-axis, the second is a direction
//       perpendicular to the y-axis.
//      *The third row is obtained by the cross product of the first
//       two rows.
//     The LocalSys matrix is [R] = [{u} {v} {w}], where {u} is the unit
//     tangent vector, {v} = {ez} x {u} or {ey} x {u} and {w} = {u} x {v}.
//
//    In both cases the returned determinant is equal to the size of the
//    vector (dx/dr,dy/dr,dz/dr), that is, equal to the length of
//    infinitesimal tangent element.
//
// -------------------------------------------------------------------------
// Each of these two classes can use the following line types:
//   Line2 => 2-noded linear line shape.
//   Line3 => 3-noded quadratic line shape.
// -------------------------------------------------------------------------
//
// The numbering scheme of the 2-noded line finite element shape is:
//
//                    s
//                    ^
//                    |    (natural coord. system)
//                    |
//                    -----> r
//
//            1 +-----------+ 2
//                  LINE2
// -------------------------------------------------------------------------
//
// The numbering scheme of the 3-noded line finite element shape is:
//
//                    s
//                    ^
//                    |    (natural coord. system)
//                    |
//                    -----> r
//
//            1 +-----+-----+ 3
//                    2            LINE3
// -------------------------------------------------------------------------

#ifndef _SHPLINE_H
#define _SHPLINE_H

#include "shape.h"

// ------------------------------------------------------------------------
// Definition of the Bar Shape class:
//
class cShapeBar : public cShape
{
 public:
                cShapeBar(void);
  virtual      ~cShapeBar(void);
          int   NumDim(void) { return 1; }
  virtual void  NodeNatCoord(sNatCoord *) { }
  virtual void  ShpFunc(sNatCoord, double *) { }
  virtual void  DrvShpRST(sNatCoord, sNatCoord *) { }
          void  DrvShpXYZ(sNatCoord, sNodeCoord *, double *,
                           sNodeCoord *) { }
  virtual void  LocalSys(sNatCoord, sNodeCoord *, cMatrix &) { }
  virtual int   GetEdge(int *, eShpType *, int *, cNode **) { return 0; }
};


// ------------------------------------------------------------------------
// Definition of the Line 2D Shape class:
//
class cShapeLine2D : public cShape
{
 public:
                cShapeLine2D(void);
  virtual      ~cShapeLine2D(void);
          int   NumDim(void) { return 1; }
  virtual void  NodeNatCoord(sNatCoord *) = 0;
  virtual void  ShpFunc(sNatCoord, double *) = 0;
  virtual void  DrvShpRST(sNatCoord, sNatCoord *) = 0;
  virtual void  DrvShpXYZ(sNatCoord, sNodeCoord *, double *,
                           sNodeCoord *);
  virtual void  LocalSys(sNatCoord, sNodeCoord *, cMatrix &);
  virtual int   GetEdge(int *, eShpType *, int *, cNode **);
};

// ------------------------------------------------------------------------
// Two-dimensional Line2 shape class:
//
class cShapeLine2 : public virtual cShapeLine2D
{
 public:
                cShapeLine2(void);
  virtual      ~cShapeLine2(void);
  virtual void  NodeNatCoord(sNatCoord *);
  virtual void  ShpFunc(sNatCoord, double *);
  virtual void  DrvShpRST(sNatCoord, sNatCoord *);
};

// ------------------------------------------------------------------------
// Two-dimensional Line3 shape class:
//
class cShapeLine3 : public virtual cShapeLine2D
{
 public:
                cShapeLine3(void);
  virtual      ~cShapeLine3(void);
  virtual void  NodeNatCoord(sNatCoord *);
  virtual void  ShpFunc(sNatCoord, double *);
  virtual void  DrvShpRST(sNatCoord, sNatCoord *);
};


// -------------------------------------------------------------------------
// Definition of the Line 3D Shape class:
//
class cShapeLine3D : public virtual cShapeLine2D
{
 public:
                cShapeLine3D(void);
  virtual      ~cShapeLine3D(void);
  virtual void  DrvShpXYZ(sNatCoord, sNodeCoord *, double *,
                          sNodeCoord *);
  virtual void  LocalSys (sNatCoord, sNodeCoord *, cMatrix &);
};

// -------------------------------------------------------------------------
// Three-dimensional Line2 Shape class:
//
class cShape3DLine2 : public virtual cShapeLine3D,
                      public virtual cShapeLine2
{
 public:
           cShape3DLine2(void);
  virtual ~cShape3DLine2(void);
};

// -------------------------------------------------------------------------
// Three-dimensional Line3 Shape class:
//
class cShape3DLine3 : public virtual cShapeLine3D,
                      public virtual cShapeLine3
{
 public:
           cShape3DLine3(void);
  virtual ~cShape3DLine3(void);
};

#endif
