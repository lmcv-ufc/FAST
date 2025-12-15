// -------------------------------------------------------------------------
// shape.h - file containing the definition of the Shape class.
// -------------------------------------------------------------------------
//
// The Shape class handle the element features related to the field
// and geometric (shape) interpolation, as the computation of shape
// functions and their derivatives.
//
// cShape
// |-- cShapeLine2D => LINE2, LINE3
// |   |-- cShapeLine3D => LINE2_3D, LINE3_3D
// |
// |-- cShapePlane => T3, T6, Q4, Q8, Q9, IGA2D.
// |   |-- cShapeSurf => T3_SURF, T6_SURF, Q4_SURF, Q8_SURF, Q9_SURF
// |
// |-- cShapeSolid => BRICK8, BRICK20, TET4, TET10
//
// -------------------------------------------------------------------------
// Static public methods:
// -------------------------------------------------------------------------
//
// static cShape *CreateShape(eShpType type)
//
//  type   -  given shape type                                        (in)
//
// This method creates a new shape object from the given shape type and
// returns a pointer to the created shape.
// -------------------------------------------------------------------------
//
// static int MaxEdgeNodes(void)
//
// This method returns the maximum number of node in a element edge.
//
// -------------------------------------------------------------------------
// Shared methods (implemented only in the base class):
// -------------------------------------------------------------------------
//
// eShpType GetType(void)
//
// This method returns the shape type.
// -------------------------------------------------------------------------
//
// eTopType GetTopologyType(void)
//
// This method returns the shape topology type.
// -------------------------------------------------------------------------
//
// int GetNumNode(void)
//
// This method returns the number of shape nodes.
// -------------------------------------------------------------------------
//
// cNode *GetNode(int i)
//
//  i      -  node index                                              (in)
//
// This method returns the pointer to the given shape node.
// -------------------------------------------------------------------------
//
// void Read(void)
//
// This method returns the shape nodes.
// -------------------------------------------------------------------------
//
// void SetNodes(cNode **node)
//
//  node   -  array of nodes (pointers)                               (in)
//
// This method stores the given nodes in the shape connectivity.
//
// -------------------------------------------------------------------------
//
// void NodalCoord(sNodeCoord *coord)
//
//  coord  -  nodal coordinates                                      (out)
//
// This method returns the cartesian coordinates of the shape nodes.
//
// -------------------------------------------------------------------------
//
// void UpdatedNodalCoord(sNodeCoord *coord)
//
//  coord  -  updated nodal coordinates                              (out)
//
// This method returns the updated cartesian coordinates of the shape nodes.
//
// -------------------------------------------------------------------------
// Pure virtual methods:
// -------------------------------------------------------------------------
//
// virtual int NumDim(void)
//
// This method returns the number of shape dimensions (1, 2 or 3) in the
// parametric space. This value correspond to the number of parametric
// coordinates required to describe the shape.
// -------------------------------------------------------------------------
//
// virtual void NodeNatCoord(sNatCoord *c)
//
//  c      -  parametric coordinates of the shape nodes              (out)
//
// This method returns the parametric coordinates of the shape nodes.
// -------------------------------------------------------------------------
//
// virtual bool Isoparametric(void)
//
// This method returns true is the shape is isoparametric (default case) and
// false otherwise.
// -------------------------------------------------------------------------
//
// virtual void ShpFunc(sNatCoord p, double *N)
//
//  p      -  parametric coordinates of the given point               (in)
//  N      -  shape functions                                        (out)
//
// This method computes the shape functions at the given point.
// -------------------------------------------------------------------------
//
// virtual void MapFunc(sNatCoord p, double *M)
//
//  p      -  parametric coordinates of the given point               (in)
//  M      -  mapping functions                                      (out)
//
// This method computes the mapping functions at the given point. For
// isoparametric shapes (default case) this function calls ShpFunc.
// -------------------------------------------------------------------------
//
// virtual void DrvShpRST(sNatCoord p, double *dNrst)
//
//  p      -  parametric coordinates of the given point               (in)
//  dNrst  -  derivatives of the shape functions                     (out)
//
// This method computes the derivatives of the shape functions w.r.t. the
// parametric coordinates at a given point.
// -------------------------------------------------------------------------
//
// virtual void DrvMapRST(sNatCoord p, double *dMrst)
//
//  p      -  parametric coordinates of the given point               (in)
//  dMrst  -  derivatives of the mapping functions                   (out)
//
// This method computes the derivatives of the mapping functions w.r.t. the
// parametric coordinates at a given point. For isoparametric shapes
// (default case) this function calls DrvShpRST.
// -------------------------------------------------------------------------
//
// virtual void DrvShpXYZ(sNatCoord p, sNodeCoord *coord, double *detJ,
//                        sNodeCoord *dNxyz)
//
//  p      -  parametric coordinates of the given point               (in)
//  coord  -  nodal coordinates                                       (in)
//  detJ   -  determinant of the Jacobian matrix                     (out)
//  dNxyz  -  derivatives of the shape functions                     (out)
//
// This method computes the derivatives of the shape functions w.r.t. the
// cartesian coordinates and the determinant of the Jacobian matrix at a
// given point.
// -------------------------------------------------------------------------
//
// virtual void LocalSys(sNatCoord p, sNodeCoord *coord, cMatrix &R)
//
//  p      -  parametric coordinates of the given point               (in)
//  coord  -  nodal coordinates                                       (in)
//  R      -  rotation matrix                                        (out)
//
// This method computes the rotation matrix from the local system to the
// global one at a given point. It can be used when loads are defined in
// local system (e.g. pressure loads are always normal to a surface). It
// returns a 3x3 rotation matrix by default, but can also return a 2x2
// rotation matrix for plane problems.
// -------------------------------------------------------------------------
//
// virtual void RMatrix(sNatCoord p, sNodeCoord *coord, double *detJ,
//                      cMatrix &R)
//
//  p      -  parametric coordinates of the given point               (in)
//  coord  -  nodal coordinates                                       (in)
//  detJ   -  Jacobian determinant matrix                            (out)
//  R      -  rotation matrix                                        (out)
//
// This method computes the Jacobian and the rotation matrix from the local
// system to the global one at a given point. It returns a 3x3 rotation 
// matrix by default, but can also return a 2x2 rotation matrix for plane
// problems. This method should be used only for interface elements.
// -------------------------------------------------------------------------
//
// virtual void tVector(sNatCoord p, sNodeCoord *coord, cVector &v)
//
//  p      -  parametric coordinates of the given point               (in)
//  coord  -  nodal coordinates                                       (in)
//  v      -  vector in parametric "t" direction                     (out)
//
// This method computes the unit vector in the parametric "t" direction,
// which corresponds to the out-of-plane (normal) direction.
// -------------------------------------------------------------------------
//
// virtual int GetEdge(int *corner, eShpType *type, int *nnode,cNode **conn)
//
//  corner -  label of the corner nodes                               (in)
//  type   -  shape type                                             (out)
//  nnode  -  number of shape nodes                                  (out)
//  conn   -  shape connectivity (nodes)                             (out)
//
// This method returns 1/0 if the given corner nodes belongs or not to the
// shape. In the former case it also returns the data require to create
// a valid shape corresponding to the given edge.
// -------------------------------------------------------------------------
//
// virtual int GetFace(int *corner, eShpType *type, int *nnode,cNode **conn)
//
//  corner -  label of the corner nodes                               (in)
//  type   -  shape type                                             (out)
//  nnode  -  number of shape nodes                                  (out)
//  conn   -  shape connectivity (nodes)                             (out)
//
// This method returns 1/0 if the given corner nodes belongs or not to the
// shape. In the former case it also returns the data require to create
// a valid shape corresponding to the given face. This method is valid only
// for solid shapes.
// -------------------------------------------------------------------------

#ifndef _SHAPE_H
#define _SHAPE_H

#include "node.h"

// ------------------------------------------------------------------------
// Forward declarations:
//
class cVector;
class cMatrix;

// -------------------------------------------------------------------------
// Auxiliary types:
//
typedef struct
{
  double r,s,t;
} sNatCoord;

typedef struct
{
  double rr,ss,rs;
} sNatDrv;

// -------------------------------------------------------------------------
// Shape topology types:
//
typedef enum
{
  LINE_TOPOLOGY,
  TRIANGULAR_TOPOLOGY,
  QUADRILATERAL_TOPOLOGY,
  TETRAHEDRAL_TOPOLOGY,
  WEDGE_TOPOLOGY,
  HEXAHEDRAL_TOPOLOGY
} eTopType;

// -------------------------------------------------------------------------
// Shape types:
//
typedef enum
{
  SHAPE_BAR,           // Bar (truss and beam elements)
  SHAPE_L2_2D,         // Linear 1D element in XY plane
  SHAPE_L3_2D,         // Quadratic 1D element in XY plane
  SHAPE_L2_3D,         // Linear 1D element in 3D space
  SHAPE_L3_3D,         // Quadratic 1D element in 3D space
  SHAPE_BSP_PLCURV,    // Isogeometric BSpline curve element in XY plane
  SHAPE_BSP_CURV,      // Isogeometric BSpline curve element in 3D space
  SHAPE_T3,            // Linear triangle
  SHAPE_T3_SURF,       // Linear triangle in 3D space
  SHAPE_T3_SHELL,      // Linear triangle in 3D space for shell elem.
  SHAPE_T6,            // Quadratic triangle
  SHAPE_T6_SURF,       // Quadratic triangle in 3D space
  SHAPE_T6_SHELL,      // Quadratic triangle in 3D space for shell elem.
  SHAPE_Q4,            // Linear quadrilateral
  SHAPE_Q4_SURF,       // Linear quadrilateral in 3D space
  SHAPE_Q8,            // Quadratic (serendipity) quadrilateral
  SHAPE_Q8_SURF,       // Quadratic quadrilateral in 3D space
  SHAPE_Q8_SHELL,      // Quadratic quadrilateral in 3D space for shell elem.
  SHAPE_Q9,            // Quadratic (lagrangian) quadrilateral
  SHAPE_Q9_SURF,       // Quadratic (lagrangian) quadrilateral in 3D space
  SHAPE_Q9_SHELL,      // Quadratic (lagrangian) quadrilateral in 3D space for shell elem.
  SHAPE_BEZ_PLCURV,    // Isogeometric Bézier curve element in XY plane
  SHAPE_BEZ_CURV,      // Isogeometric Bézier curve element in 3D space
  SHAPE_BEZ_PLSURF,    // Isogeometric Bézier surface element in XY plane
  SHAPE_BEZ_SURF,      // Isogeometric Bézier surface element in 3D space
  SHAPE_BEZ_SHELL,     // Isogeometric Bézier surface element in 3D space for shell elem.
  SHAPE_BEZTRI_PLSURF, // Isogeometric Bézier triangle surface element in XY plane
  SHAPE_BEZTRI_SURF,   // Isogeometric Bézier triangle surface element in 3D space
  SHAPE_BEZTRI_SHELL,  // Isogeometric Bézier triangle surface element in 3D space for shell elem.
  SHAPE_BSP_PLSURF,    // Isogeometric BSpline surface element in XY plane
  SHAPE_BSP_SURF,      // Isogeometric BSpline surface element in 3D space
  SHAPE_BSP_SHELL,     // Isogeometric BSpline surface element in 3D space for shell elem.
  SHAPE_BRICK8,        // Linear brick (hexahedral) element
  SHAPE_BRICK20,       // Quadratic serendipity brick (hexahedra) element
  SHAPE_TET4,          // Linear tetrahedral element
  SHAPE_TET10,         // Quadratic tetrahedral element
  SHAPE_INTERF_L2,     // Linear line interface element in XY plane
  SHAPE_INTERF_L3,     // Quadratic line interface element in XY plane
  SHAPE_INTERF_Q4,     // Linear quadrilateral interface element
  SHAPE_INTERF_Q8,     // Quadratic quadrilateral interface element
  SHAPE_BSP_INTERF_LINE,  // Isogeometric BSpline line interface element
  SHAPE_BSP_SOLID,     // Isogeometric BSpline solid element
  SHAPE_L6_INF         // Mapped infinite element (6-node Lagrangian)
} eShpType;

// -------------------------------------------------------------------------
// Edge structure:
//
const int MAX_EDGE_NODE = 3;

typedef struct _edge
{
  eShpType type;               // Shape type
  int nnode;                   // Number of nodes
  int node[MAX_EDGE_NODE];     // Conectivity
} tEdge;


// -------------------------------------------------------------------------
// Face structure:
//
const int MAX_FACE_CORNER = 4;
const int MAX_FACE_NODE = 9;

typedef struct _face
{
  eShpType type;               // Shape type
  int ncorner;                 // Number of face corners
  int corner[MAX_FACE_CORNER]; // Corner nodes
  int nnode;                   // Number of face nodes
  int node[MAX_FACE_NODE];     // Conectivity
} tFace;


// ------------------------------------------------------------------------
// Definition of the Shape class:
//
class cShape
{
 protected:
  eShpType  Type;       // Shape type
  eTopType  TopType;    // Shape topology type
  int       NumNode;    // Number of element nodes
  cNode   **ElmNode;    // Element nodes

 protected:
  void   ShpMatZero(int n, double [3][3]);
  int    ShpMatInv(int n, double [3][3], double [3][3], double *det);
  double ShpVecLen3D(double *);
  void   ShpNormVec3D(double *);
  void   ShpCrossProd(double *, double *, double *);
  int    GetShapeEdge(int, tEdge *, int *, eShpType *, int *, cNode **);
  int    GetShapeFace(int, tFace *, int *, eShpType *, int *, cNode **);

 public:
  static cShape    *CreateShape(eShpType);
  static int        MaxEdgeNodes(void)    { return 3;          }

                    cShape(void);
  virtual          ~cShape(void);
          eShpType  GetType(void)         { return Type;         }
          eTopType  GetTopologyType(void) { return TopType;      }
          int       GetNumNode(void)      { return NumNode;      }
          cNode    *GetNode(int i)        { return ElmNode[i];   }
  virtual int       GetNumShp(void)       { return NumNode;      }
  virtual int       GetNumMap(void)       { return GetNumShp( ); }
  virtual void      Read(void);
  virtual void      SetNodes(cNode **, int size = 0);
  virtual void      SetNodalAxes(sNodeAxes*) { }
  virtual void      UpdNodalAxes(sNodeAxes*) { }
  virtual void      NodalCoord(sNodeCoord *);
          void      UpdatedNodalCoord(sNodeCoord *);
          void      Evaluate(sNatCoord,double*,sNodeCoord*,sNodeCoord&);
  virtual int       NumDim(void) = 0;
  virtual void      NodeNatCoord(sNatCoord *) = 0;
  virtual bool      Isoparametric(void) { return(true); }
  virtual void      ShpFunc(sNatCoord, double *) = 0;
  virtual void      MapFunc(sNatCoord p, double *M) { ShpFunc(p, M); }
  virtual void      DrvShpRST(sNatCoord, sNatCoord *) = 0;
  virtual void      DrvShpRST(sNatCoord, sNatCoord *, sNatDrv *) { }
  virtual void      DrvMapRST(sNatCoord p, sNatCoord *dM) { DrvShpRST(p, dM); };
  virtual void      DrvShpXYZ(sNatCoord, sNodeCoord *, double *,
                               sNodeCoord *) = 0;
  virtual void      DrvShpXYZ(sNatCoord, sNodeCoord *, double *,
                               sNodeCoord *, sNodeDrv *) { }
  virtual void      NMatrix(int,double*,cMatrix&);
  virtual void      LocalSys(sNatCoord, sNodeCoord *, cMatrix &);
  virtual void      RMatrix(sNatCoord, sNodeCoord *, double *, cMatrix &);
  virtual void      tVector(sNatCoord, sNodeCoord *, cVector &);
  virtual int       GetEdge(int *, eShpType *, int *, cNode **) = 0;
  virtual int       GetFace(int *, eShpType *, int *, cNode **) { return 0; }
  virtual cShape   *GetBaseShp(void) { return this; }
};

#endif
