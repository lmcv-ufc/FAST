// -------------------------------------------------------------------------
// shapeiga.h - file containing the definition of the Shape class.
// -------------------------------------------------------------------------
//
// The Shape IGA class stores data and implements methods required by
// isogeometric elements.
//
// cShapeIGA
// |-- cShapeBsp
// |   |-- cShapePlBspCurve
// |   |-- cShapeBspCurve
// |   |-- cShapePlBspSurf
// |   |-- cShapeBspSurf
// |   |-- cShapeBspSolid
// |
// |-- cShapeBezCurve
// |-- cShapePlBezTrianSurf
// |-- cShapePlBezSurf
// |
// |-- cShapeTSpline (not implemented)
// 
// -------------------------------------------------------------------------
// Protected methods:
// -------------------------------------------------------------------------
//
// int InitEqvMshQuad(eTopType t)
//
//  t      -  given topology type                                  (in)
//
// This method initialize respective equivalent mesh quadrature object
// associated with given topology type. 
//
// -------------------------------------------------------------------------
// Public methods:
// -------------------------------------------------------------------------
//
// void ReadPatPnt(void)
//
// This method read patch point list.
// -------------------------------------------------------------------------
//
// void ReadEqvMshDiv(void)
//
// This method read equivalent mesh division in each parametric dimension.
// -------------------------------------------------------------------------
//
// void Destroy(void)
//
// this method destroy equivalent mesh data.
// -------------------------------------------------------------------------
//
// void LoadEqvMsh(void)
//
// This method load equivalent mesh data (nodes and elements).
// -------------------------------------------------------------------------
//
// sEqvElmQuad& GetEqvElmPnt(eTopType t)
//
//  t      -  given topology type                                  (in)
//
// This method return the reference of respective equivalent mesh quadrature
// object associated with given topology type. The returned object is persistent
// (singleton pattern).
// -------------------------------------------------------------------------
//
// void PrintEqvMshDispl(std::ostream &out)
//
//   out - given output stream object                              (out)
//
// This method write equivalent mesh displacements in given output stream
// object.
// -------------------------------------------------------------------------
//
// void PrintEqvMshData(std::ostream &out)
//
//   out - given output stream object                              (out)
//
// This method write equivalent mesh data in given output stream object.
// -------------------------------------------------------------------------
//
// void ReadIGA(void)
//
// This method read isogeomtric shape data.
// -------------------------------------------------------------------------
//
// cPatch* GetPatch(void)
//
// This method return shape patch pointer.
// -------------------------------------------------------------------------
//
// cPatch* GetPatEdge(int f, int e, eShpType &t, int s[])
//
//   f - given face location                              (in)
//   e - given edge location                              (in)
//   t - shape type                                       (out)
//   s - span in each parametric dimension                (out)
//
// This method return curve patch corresponding to give face and edge locations.
// The face locations are Front, Back, Left, Right, Above, Below. The edge
// locations are Top, Bottom, Left and Right.
// -------------------------------------------------------------------------
//
// cPatch* GetPatFace(int f, eShpType &t, int s[])
//
//   f - given face location                              (in)
//   t - shape type                                       (out)
//   s - span in each parametric dimension                (out)
//
// This method return surface patch corresponding to give face.
// -------------------------------------------------------------------------
//
// void Init(void)
//
// This method initialize nodal (control points) incidence from patch control
// point vector.
//
// -------------------------------------------------------------------------
// Wrapper functions:
// -------------------------------------------------------------------------
//
// void ReadPatch(void) 
//
// This method call ReadPatch from cPatch class.
// -------------------------------------------------------------------------
//
// void InitIGA(void)
//
// This method set custom control points reading function.
//
// -------------------------------------------------------------------------

#ifndef _SHPIGA_H
#define _SHPIGA_H

#include "ctrl.h"
#include "patch.h"
#include "shape.h"
#include "shpline.h"
#include "shpplane.h"
#include "shpsurf.h"
#include "shpsolid.h"
#include "kdtree.h"

using namespace std;

#include <vector>
#include <map>

#include <iosfwd>
class cElement;

// -------------------------------------------------------------------------
// Edge types:
//
typedef enum
{
  EDGE_TOP,
  EDGE_BOTTOM,
  EDGE_LEFT,
  EDGE_RIGHT
} eEdgeType;

ifstream& operator>> (ifstream&,eEdgeType&);

// -------------------------------------------------------------------------
// Face types:
//
typedef enum
{
  FACE_FRONT,
  FACE_BACK,
  FACE_LEFT,
  FACE_RIGHT,
  FACE_ABOVE,
  FACE_BELOW
} eFaceType;

istream& operator>> (istream&,eFaceType&);

// -------------------------------------------------------------------------
// Equivalent node structure:
//
struct sEqvNode
{
  int label;
  int elm;
  sNatCoord  p;
  sNodeCoord coord;

  sEqvNode(void) { }
  sEqvNode(const sEqvNode&);
};

// -------------------------------------------------------------------------
// Equivalent element structure:
//
struct sEqvElm
{
  static int        buffsize;  // Node vector buffer size.
         int        nn;        // Number of nodes.
         int        label;     // Element label.
	 cElement*  parentelm; // Parent element.
         sEqvNode **nodes;     // Element nodes.

  sEqvElm(void);
  sEqvElm(const sEqvElm&);
 ~sEqvElm(void);

 void operator=(const sEqvElm&);
};

// -------------------------------------------------------------------------
// Equivalent element quadrature structure:
//
struct sEqvElmQuad
{
  int NumPnt;
  int NumElm;
  int NumElmPnt;

  std::vector<sNatCoord> coords;
};

// -------------------------------------------------------------------------
// Struct patch point:
//
struct sPatPnt
{
  int    patid;
  double parcoord[3];
};

// -------------------------------------------------------------------------
// Typedefs simplifications:
//
typedef std::vector<sEqvElm>                      sEqvElmStdVec;
typedef std::vector<sEqvNode>                     sEqvNodeStdVec;
typedef std::map<cPatch*,cElement**>              PatElmStdMap;
typedef std::map<int,cPatch*>                     IsoPatStdMap;
typedef std::vector<sPatPnt>                      PatPntStdVec;
typedef std::vector<cElement*>                    ElmStdVec;
typedef std::map<eTopType, std::vector<sEqvElm> > EqvElmStdMap;
typedef cKDTree<std::vector<sEqvNode*> >          cEqvNodeKDT; 

// -------------------------------------------------------------------------
// Isogeometric Shape class:
//
class cShapeIGA 
{
 protected:
  static PatElmStdMap     PatElmMap;
  static IsoPatStdMap     IsoPatMap;
  static PatPntStdVec     PatPntVec;
  static cEqvNodeKDT      EqvNodeTree;
  static ElmStdVec        ParElmVec;
  static EqvElmStdMap     EqvElmMap;
  static double           EqvMshTol;
  static int              NumEqvNode;
  static int              NumEqvElm;
  static int              EqvMshSL;
  static int              EqvMshDivRST[3];
                          
  static void             InitEqvMshQuad(eTopType);

 public:
                          cShapeIGA(void) { }
  virtual                ~cShapeIGA(void) { }

  static  void            ReadPatPnt(void);
  static  void            ReadEqvMshDiv(void);

  static  int             GetNumEqvNode(void) { return NumEqvNode; }
  static  int             GetNumEqvElm(void)  { return NumEqvElm;  }
  static  int             GetEqvMshDisc(int d)  { return EqvMshDivRST[d];  }
  static  sEqvElmQuad&    GetEqvElmPnt(eTopType);
  static  void            LoadEqvMsh(void);
  static  void            Destroy(void);
  static  void            PrintEqvMshDispl(std::ostream&);
  static  void            PrintEqvMshTemp(std::ostream&);
  static  void            PrintEqvMshData(std::ostream&);
  static  void            PrintMode(int,double,cVector&,ostream&,string,string);

  // Virtual functions.
  virtual void            ReadIGA(void) { }
  virtual cPatch*         GetPatch(void) { return 0; }
  virtual cPatch*         GetPatEdge(int,int,eShpType&,int[]) {return 0;}
  virtual cPatch*         GetPatFace(int,eShpType&,int[]) {return 0;}
  virtual void            Init(void) = 0;

  // Define const object access functions. 
  static const PatElmStdMap& GetPatElmMap(void)    { return PatElmMap; }
  static const double&       GetEqvMshTol(void)    { return EqvMshTol; }
  static const PatPntStdVec& GetPatPntVec(void)    { return PatPntVec; }
  static const ElmStdVec&    GetParentElmVec(void) { return ParElmVec; }
  static const EqvElmStdMap& GetEqvElmStdMap(void) { return EqvElmMap; }
};

// Global functions.
void          ReadPatch(void);
void          InitIGA(void);
std::istream& ReadCP(std::istream&,cPatch&);
std::string   GetEqvElmTag(const eTopType&);

#endif
