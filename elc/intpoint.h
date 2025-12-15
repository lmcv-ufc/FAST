// ------------------------------------------------------------------------
// intpoint.h - file containing the definition of the Integration Point
//              class.
// ------------------------------------------------------------------------
//
// The Integration Point class basically stores the coordinates and weights
// of the integration points used in the numerical integration.
//
// -------------------------------------------------------------------------
// Static public methods:
// -------------------------------------------------------------------------
//
// static cIntPoint *CreateIntPoints(eShpType type, int order, int *n)
//
//  type   -  shape type                                              (in)
//  order  -  integration order                                       (in)
//  n      -  number of integration points                           (out)
//
// This method creates an array of integration points corresponding to
// given shape and order.
//
// -------------------------------------------------------------------------
// Shared methods (implemented only in the base class):
// -------------------------------------------------------------------------
//
// double GetWeight(void)
//
// This method returns the integration weight.
// -------------------------------------------------------------------------
//
// sNatCoord GetCoord(void)
//
// This method returns the parametric coordinates of the integration point.
//
// -------------------------------------------------------------------------
//
// ref: Dunavant, D. A. High degree efficient symmetrical Gaussian quadrature
//	rules for the triangle. International Journal for Numerical Methods in 
//      Engineering, vol. 21, issue 6, 1985, pp. 1129-1148.
// 
// ref: James Lyness, Dennis Jespersen,
//      Moderate Degree Symmetric Quadrature Rules for the Triangle,
//      Journal of the Institute of Mathematics and its Applications,
//      Volume 15, Number 1, February 1975, pages 19-32. 
//
// ref: Gellert, M., & Harbord, R. (1991). Moderate degree cubature formulas 
//      for 3-D tetrahedral finite-element approximations. Communications in 
//      Applied Numerical Methods, 7(6), 487–495. doi:10.1002/cnm.1630070609
// -------------------------------------------------------------------------

#ifndef _INTPOINT_H
#define _INTPOINT_H

#include "shape.h"

using namespace std;

// -------------------------------------------------------------------------
// Forward declarations:
//
#include <iosfwd>

// -------------------------------------------------------------------------
// Quadrature type:
//
typedef enum
{
  GAUSS,         // n-dimensional Gauss integration (tensor product)
  LOBATTO,       // n-dimensional Lobatto integration (tensor product)
  NEWTON_COTES,  // n-dimensional Newton-Cotes integration (tensor product)
  DUNAVANT,      // triangle integration - Dunavant (1985)
  LYNESS,        // triangle integration - Lyness and Jespersen (1975)
  TETRAHEDRA     // tetrahedra integration
} eQuadType;

// -------------------------------------------------------------------------
// Auxiliary functions:
//

istream& operator>> (istream &, eQuadType &);

// ------------------------------------------------------------------------
// Definition of the Integration Point class:
//
class cIntPoint
{
 protected:
  double    wgt;    // Integration weight
  sNatCoord pnt;    // Point coordinate (r,s,t)

 public:
  static cIntPoint *CreateIntPoints(eTopType , int, int[], int *);
  static cIntPoint *CreateLinePoints(int, eQuadType, int *);
  static cIntPoint *CreateTriaPoints(int, eQuadType, int *);
  static cIntPoint *CreateQuadPoints(int[], eQuadType, int *);
  static cIntPoint *CreateBrickPoints(int[], eQuadType, int *);
  static cIntPoint *CreateWedgePoints(int[], eQuadType, int *);
  static cIntPoint *CreateTetPoints(int, eQuadType, int *);
  static void       SetTrianQuad(eQuadType);
//  static cIntPoint *CreateLobattoLinePoints(int, int *);
//  static cIntPoint *CreateNCLinePoints(int, int *);
//  static cIntPoint *CreateNCQuadPoints(int, int *);

                    cIntPoint(void);
  virtual          ~cIntPoint(void);
         double     GetWeight(void){ return wgt; }
         sNatCoord  GetCoord(void) { return pnt; }
         void       SetWeight(double w) { wgt = w; }
         void       SetRCoord(double r) { pnt.r = r; }
         void       SetSCoord(double s) { pnt.s = s; }
         void       SetTCoord(double t) { pnt.t = t; }
};

#endif
