// -------------------------------------------------------------------------
// shpsurf.h - Definitions of curve shape class.
// -------------------------------------------------------------------------
// Description:
//
// This class is intended to deal with curved area shapes in space. Thus,
// each derived shape has two parametric coordinates (r,s) which are mapped
// on three cartesian coordinates (x,y,z). This class is used in the compu-
// tation of equivalent forces in the loaded faces of solid elements.
//
// -------------------------------------------------------------------------
// Public methods:
// -------------------------------------------------------------------------
//
// void DerivXYZ( double **invjac,
//                sDerivNat *derivshp_rst,
//                sDerivCart *deriv_xyz )
//
//      invjac        -  inverse of the jacobian matrix          (  in )
//      derivshp_rst  -  local shape functions derivatives       (  in )
//      deriv_xyz     -  global shape functions derivatives      ( out )
//
//                 -     -                     -     -
//                | dN/dx |                   | dN/dr |
//                |       |                   |       |
//                | dN/dy | =  [invjac] *     | dN/ds |
//                |       |                   |       |
//                | dN/dz |                   |   0   | , where dN/dt = 0.
//                 -     -                     -     -
//
// This method evaluates the global shape functions derivatives at a given
// integration point.
//
// -------------------------------------------------------------------------
//
// void Jacobian( sDerivNat  *derivmap_rst,
//                double     *mapfunc,
//                sNodeCoord *coord,
//                double     *detjac,
//                double    **jac )
//                double    **invjac )
//
//      derivmap_rst  -  local map functions derivatives         (  in )
//      mapfunc       -  local map functions                     (  in )
//      coord         -  nodal coordinates                       (  in )
//      detjac        -  determinant of the jacobian matrix      ( out )
//      jac           -  jacobian matrix                         ( out )
//      invjac        -  inverse of the jacobian matrix          ( out )
//
// This method evaluates the jacobian matrix, its inverse and its
// determinant. The Jacobian matrix is computed as follows:
//    The first two rows are the derivatives of the position vector w.r.t.
// the parametric coordinates. It should be noted that vectors,
// (dx/dr,dy/dr,dz/dr) and (dx/ds,dy/ds,dz/ds), are tangent to the
// element surface.
//    The third row is obtained by the cross product of the first two rows
// and represents the normal vector.
//    The returned determinant is computed as the length of the normal
// vector and represents the infinitesimal area in parametric
// coordinates.
//
//                    -                     -
//                   | dx/dr   dy/dr   dz/dr |
//                   |                       |
//          [jac] =  | dx/ds   dy/ds   dz/ds |
//                   |                       |
//                   |   n1      n2      n3  |
//                    -                     -
//          where,
//
//               n1 = (dy/dr * dz/ds) - (dz/dr * dy/ds)
//               n2 = (dz/dr * dx/ds) - (dx/dr * dz/ds)
//               n3 = (dx/dr * dy/ds) - (dy/dr * dx/ds)
//
//               det_jac = sqrt( n1*n1 + n2*n2 + n3*n3 )
//
// -------------------------------------------------------------------------
//
// void LocalSys( sDerivNat  *derivmap_rst,
//                sNodeCoord *coord,
//                double    **rot )
//
//      derivmap_rst  -  local map functions derivatives         (  in )
//      coord         -  element coordinates                     (  in )
//      rot           -  rotation matrix                         ( out )
//
// This method computes the rotation matrix of a local system of coordina-
// tes evaluated at a point on the area. The rotation matrix is  formed as
// follows:
//    Each column of the rotation matrix is formed by the coordinates of
// each local system unity vector in the global system.
//    The non-normalized coordinates of the local system directions
// are considered as:
//
//           v1 = ( dx/dr, dy/dr, dz/dr )
//           v2 = ( dx/ds, dy/ds, dz/ds )
//           v3 = (   n1 ,   n2 ,   n3  )
//
//           where,
//               n1 = (dy/dr * dz/ds) - (dz/dr * dy/ds)
//               n2 = (dz/dr * dx/ds) - (dx/dr * dz/ds)
//               n3 = (dx/dr * dy/ds) - (dy/dr * dx/ds)
//
// -------------------------------------------------------------------------

#ifndef	_SHPSURF_H
#define	_SHPSURF_H

#include "shpplane.h"

// ------------------------------------------------------------------------
// Surface Shape class
//
class cShapeSurf : public virtual cShapePlane
{
  public:
               cShapeSurf(void);
  virtual     ~cShapeSurf(void);
  virtual void DrvShpXYZ (sNatCoord, sNodeCoord *, double *, sNodeCoord *);
  virtual void LocalSys  (sNatCoord, sNodeCoord *, cMatrix &);
  virtual void tVector   (sNatCoord, sNodeCoord *, cVector &);
};

// -------------------------------------------------------------------------
// Curve T3 Shape class:
//
class cShapeT3Surf : public virtual cShapeSurf,
                     public virtual cShapeT3
{
 public:
              cShapeT3Surf(void);
  virtual    ~cShapeT3Surf(void);
  virtual int GetEdge(int *, eShpType *, int *, cNode **);
};

// -------------------------------------------------------------------------
// Curve T6 Shape class:
//
class cShapeT6Surf : public virtual cShapeSurf,
                     public virtual cShapeT6
{
 public:
              cShapeT6Surf(void);
  virtual    ~cShapeT6Surf(void);
  virtual int GetEdge(int *, eShpType *, int *, cNode **);
};

// -------------------------------------------------------------------------
// Curve Q4 Shape class:
//
class cShapeQ4Surf : public virtual cShapeSurf,
                     public virtual cShapeQ4
{
 public:
              cShapeQ4Surf(void);
  virtual    ~cShapeQ4Surf(void);
  virtual int GetEdge(int *, eShpType *, int *, cNode **);
};

// -------------------------------------------------------------------------
// Curve Q8 Shape class:
//
class cShapeQ8Surf : public virtual cShapeSurf,
                     public virtual cShapeQ8
{
 public:
              cShapeQ8Surf(void);
  virtual    ~cShapeQ8Surf(void);
  virtual int GetEdge(int *, eShpType *, int *, cNode **);
};

// -------------------------------------------------------------------------
// Curve Q9 Shape class:
//
class cShapeQ9Surf : public virtual cShapeSurf,
                     public virtual cShapeQ9
{
 public:
              cShapeQ9Surf(void);
  virtual    ~cShapeQ9Surf(void);
  virtual int GetEdge(int *, eShpType *, int *, cNode **);
};

#endif

