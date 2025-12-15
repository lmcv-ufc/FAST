// -------------------------------------------------------------------------
// load.h - file containing the definition of the cLoad class.
// -------------------------------------------------------------------------
//
// The cLoad class defines the general behavior of a external load in a
// finite element system and implements the methods which are common for the
// different load types. The specific features of each load are defined
// by a set of pure virtual methods implemented in the derived classes.
// The cLoad class is also responsible for the storage of all loads
// of the current mesh.
//
// Each load can vary in time (t) and space (x). In order to handle these
// features each load f is computed as:
//
//   f(x, t) = q(x)*h(t)
//
// where q(x) represents the variation in space (nodal/element) and h(t) is
// the associated time function. Since the neutral file description does
// not include time functions, it was adopted the concept of a current time
// function which, after read from the input file, is associated with all
// the next load definitions. In order to handle static cases a constant
// unit function is defined as the default time function.
//
// cLoad
// |-- cNodalLoad
// |-- cElemLoad => cParConcLoad, cPatConcLoad, cBodyLoad, cEdgeUnifLoad,
// |                cEdgeInifLoadIGA, cFaceUnifLoad, cFaceUnifLoadIGA,
// |                cPlFrameUnifLoad, cFrame3DUnifLoad, cGridUnifLoad
// |                cDonnellUnifLoad
//
// Each load object has a pointer to its parent entity (node or element).
// Moreover, for the element nodes also has an associated shape object which
// is automatically created by the program based on the parent element and
// on the load type: edge, face or body forces. Using these variables are
// used in the computation of the equivalent nodal forces and its summation
// in the global vector.
//
// -------------------------------------------------------------------------
// Static public methods:
// -------------------------------------------------------------------------
//
// static void ReadNodalLoads(void)
//
// This method creates and reads the nodal loads.
// -------------------------------------------------------------------------
//
// static void ReadBodyForces(void)
//
// This method creates and reads the body forces.
// -------------------------------------------------------------------------
//
// static void ReadLineUnifLoads(void)
//
// This method creates and reads the uniforme edge (line) loads.
// -------------------------------------------------------------------------
//
// static void ReadFaceUnifLoads(void)
//
// This method creates and reads the uniform face (area) loads.
// -------------------------------------------------------------------------
//
// static void ReadPlFrameUnifLoads(void)
//
// This method creates and reads the plane frame uniform loads.
// -------------------------------------------------------------------------
//
// static void ReadFrame3DUnifLoads(void)
//
// This method creates and reads the three-dimensional frame uniform loads.
// -------------------------------------------------------------------------
//
// static void Destroy(void)
//
// This method destroys all stored loads.
// -------------------------------------------------------------------------
//
// static int GetMaxDof(void)
//
// This method returns the maximum number of load dofs.
// -------------------------------------------------------------------------
//
// static cLoad *GetHead(void)
//
// This method returns the head of the load list.
//
// -------------------------------------------------------------------------
// Shared methods (implemented only in the base class):
// -------------------------------------------------------------------------
//
// cLoad *GetNext(void)
//
// This method returns next load in the list.
//
// -------------------------------------------------------------------------
// Pure virtual methods:
// -------------------------------------------------------------------------
//
// virtual void ExtForce(double t, cVector &f)
//
//  t      -  time                                                    (in)
//  f      -  equivalent force vector                                (out)
//
// This method computes the nodal force vector equivalent to the applied
// load at the given time.
// -------------------------------------------------------------------------
//
// virtual void AddGlobVec(cVector &felm, cVector &fglob)
//
//  felm   -  element load vector                                     (in)
//  fglob  -  global load vector                                  (in/out)
//
// This method adds the given element load vector to the global load vector.
// -------------------------------------------------------------------------
//
// bool HasSupp(void)
//
// This method returns true if the load element has a supported node and
// false otherwise.
// -------------------------------------------------------------------------
//
// void AddReactions(cVector &felm, cMatrix &R)
//
//   felm - external force vector                                      (in)
//   R    - nodal reactions                                        (in/out)
//
// This method adds the external forces (with negative sign) to the nodal
// support reactions: R = g - f, where R is a matrix (nnode x 6).
// -------------------------------------------------------------------------

#ifndef _LOAD_H
#define _LOAD_H

#include "vec.h"
#include "node.h"
#include "shape.h"
#include "gbldef.h"

// ------------------------------------------------------------------------
// Forward declarations:
//
class cElement;
class cMatrix;
class cTimeFunc;

// -------------------------------------------------------------------------
// Load types:
//
typedef enum
{
  NODAL_LOAD,       // Nodal loads
  PARCONC_LOAD,     // Concentrate load defined by element parametric coordinates
  BODY_LOAD,        // Body (domain) loads in parametric elements
  LINE_LOAD,        // Line (edge) loads in parametric elements
  LINE_PRESS_LOAD,  // Line (edge) pressure loads in parametric elements
  FACE_LOAD,        // Face (area) loads in parametric elements
  FACE_PRESS_LOAD,  // Face (area) pressure loads in parametric elements
  PLFRAME_LOAD,     // Plane frame element loads
  GRID_LOAD,        // Grid element load 
  FRAME3D_LOAD,     // Three-dimensional frame element loads
  DONNELL_LOAD      // Donnell shell element loads
} eLoadType;

// -------------------------------------------------------------------------
// Definition of the load class:
//
class cLoad
{
 private:
  static  cLoad  *Head;
  static  cLoad  *Tail;

 protected:
          cLoad     *Next;
          cTimeFunc *Func;
          eLoadType  Type;
	  bool       IsFollower;

 public:
  static  void   ReadNodalLoads(void);
  static  void   ReadParConcForces(void);
  static  void   ReadPatConcForces(void);
  static  void   ReadBodyForces(void);
  static  void   ReadShellUnifLoads(void);
  static  void   ReadLineUnifLoads(void);
  static  void   ReadLineUnifMoments(void);
  static  void   ReadLineGenLoads(void);
  static  void   ReadFaceUnifLoads(void);
  static  void   ReadLineUnifIgaLoads(void);
  static  void   ReadLineUnifIgaMoments(void);
  static  void   ReadFaceUnifIgaLoads(void);
  static  void   ReadLineGenIgaLoads(void);
  static  void   ReadPlFrameUnifLoads(void);
  static  void   ReadFrame3DUnifLoads(void);
  static  void   ReadGridUnifLoads(void);
  static  void   ReadLineFieldLoads(void);
  static  void   ReadDonnellUnifLoads(void);
  static  void   EqvForces(cElement *, double, cVector &);
  static  void   Destroy(void);
  static  int    GetMaxDof(void) { return 500; } //TODO discutir com evandro
  static  cLoad *GetHead(void) { return Head; }
                 cLoad(void);
  virtual       ~cLoad(void);
          cLoad *GetNext(void) { return Next; }
  virtual void   ExtForce(double, cVector &) = 0;
  virtual void   AddGlobVec(cVector &, cVector &) = 0;
  virtual bool   HasSupp(void) = 0;
  virtual void   AddReactions(cVector &, cMatrix &) = 0;
};


// -------------------------------------------------------------------------
// Definition of the nodal load class:
//
class cNodalLoad : public cLoad
{
 protected:
  cNode  *Node;
  double  LoadVal[6];

 public:
                cNodalLoad(cNode *, double, double, double,
                           double, double, double);
  virtual       ~cNodalLoad(void);
          void   ExtForce(double, cVector &);
          void   AddGlobVec(cVector &, cVector &);
          bool   HasSupp(void) { return Node->HasSupp( ); }
          void   AddReactions(cVector &, cMatrix &);
};


// -------------------------------------------------------------------------
// Definition of the element load class:
//
class cElemLoad : public cLoad
{
 public:                // TEMPORARIO
  cElement  *Elem;      // Parent element
  cShape    *Shape;     // Load shape

 public:
                 cElemLoad(void);
  virtual       ~cElemLoad(void);
  virtual void   ExtForce(double, cVector &);
  virtual void   AddGlobVec(cVector &, cVector &);
  virtual bool   HasSupp(void);
  virtual void   AddReactions(cVector &, cMatrix &);
  virtual void   Read(void) = 0;
  virtual void   GetVal(sNatCoord, sNodeCoord *, double *) = 0;

 protected:
          int    GetNumDofs(void);
          void   GetDofs(int *);
};


// -------------------------------------------------------------------------
// Definition of the parameter concentrate force class:
//
class cParConcLoad : public cElemLoad
{
 protected:
   cVector LoadVal;
   double ParVal[3];

 public:
                 cParConcLoad(void);
  virtual       ~cParConcLoad(void);
  virtual void   Read(void);
          void   GetVal(sNatCoord, sNodeCoord *, double *) { }
  virtual void   ExtForce(double, cVector &);
};

// -------------------------------------------------------------------------
// Definition of the body force class:
//
class cBodyLoad : public cElemLoad
{
 protected:
   cVector LoadVal;

 public:
                 cBodyLoad(void);
  virtual       ~cBodyLoad(void);
          void   Read(void);
          void   GetVal(sNatCoord, sNodeCoord *, double *);
};

// -------------------------------------------------------------------------
// Definition of the body force class:
//
class cShellUnifLoad : public cElemLoad
{
 protected:
   int     Local;            // Flag for loads in the local system (0/1)
   cVector LoadVal;

 public:
                 cShellUnifLoad(void);
  virtual       ~cShellUnifLoad(void);
          void   Read(void);
          void   GetVal(sNatCoord, sNodeCoord *, double *);
};

// -------------------------------------------------------------------------
// Definition of the edge uniform force class:
//
class cEdgeUnifLoad : public cElemLoad
{
 protected:
  int     Local;            // Flag for loads in the local system (0/1)
  cVector LoadVal;

 public:
                 cEdgeUnifLoad(void);
  virtual       ~cEdgeUnifLoad(void);
  virtual void   Read(void);
  virtual void   GetVal(sNatCoord, sNodeCoord *, double *);
};


// -------------------------------------------------------------------------
// Definition of the edge uniform moment class:
//
class cEdgeUnifMoment : public cEdgeUnifLoad
{
 public:
                 cEdgeUnifMoment(void);
  virtual       ~cEdgeUnifMoment(void);
  virtual void   GetVal(sNatCoord, sNodeCoord *, double *);
};

// -------------------------------------------------------------------------
// Definition of the edge field force class:
//
class cEdgeFieldLoad : public cElemLoad
{
 public:
                 cEdgeFieldLoad(void);
  virtual       ~cEdgeFieldLoad(void);
  virtual void   Read(void);
  virtual void   GetVal(sNatCoord, sNodeCoord *, double *);
};

// -------------------------------------------------------------------------
// Definition of the edge general force class:
//
class cEdgeGenLoad : public cElemLoad
{
 protected:
  sExpEval expr;      // Expression evaluator.

  double x, y, z;     // Input/Output variables used
  double fx, fy, fz;  // in the expression evaluator.

 public:
                 cEdgeGenLoad(void);
  virtual       ~cEdgeGenLoad(void);
  virtual void   Read(void);
  virtual void   GetVal(sNatCoord, sNodeCoord *, double *);
};


// -------------------------------------------------------------------------
// Definition of the edge uniform pressure force class:
//
class cEdgeUnifPressLoad : public cElemLoad
{
 protected:
  double PressVal;          // Pressure value

 public:
                 cEdgeUnifPressLoad(void);
  virtual       ~cEdgeUnifPressLoad(void);
  virtual void   Read(void);
  virtual void   GetVal(sNatCoord, sNodeCoord *, double *);
};

// -------------------------------------------------------------------------
// Definition of the face uniform force class:
//
class cFaceUnifLoad : public cElemLoad
{
 protected:
  int     Local;            // Flag for loads in the local system (0/1)
  cVector LoadVal;

 public:
                 cFaceUnifLoad(void);
  virtual       ~cFaceUnifLoad(void);
  virtual void   Read(void);
  virtual void   GetVal(sNatCoord, sNodeCoord *, double *);
};

// -------------------------------------------------------------------------
// Definition of the face uniform pressure force class:
//
class cFaceUnifPressLoad : public cElemLoad
{
 protected:
  double PressVal;          // Pressure value

 public:
                 cFaceUnifPressLoad(void);
  virtual       ~cFaceUnifPressLoad(void);
  virtual void   Read(void);
  virtual void   GetVal(sNatCoord, sNodeCoord *, double *);
};

// -------------------------------------------------------------------------
// Definition of the plane frame uniform element load class:
//
class cPlFrameUnifLoad : public cElemLoad
{
 protected:
  int     Local;            // Flag for loads in the local system (0/1)
  cVector LoadVal;

 public:
                 cPlFrameUnifLoad(void);
  virtual       ~cPlFrameUnifLoad(void);
          void   ExtForce(double, cVector &);
          void   Read(void);
          void   GetVal(sNatCoord, sNodeCoord *, double *);
};

// -------------------------------------------------------------------------
// Definition of grid uniform element load class:
//
class cGridUnifLoad : public cElemLoad
{
 protected:
  int       Local;             // Flag for loads in the local system (0/1)  
  cVector   LoadVal;
    
 public:
                 cGridUnifLoad(void);
  virtual       ~cGridUnifLoad(void);
          void   Read(void);
          void   ExtForce(double, cVector &);
          void   GetVal(sNatCoord, sNodeCoord *, double *);
};

// -------------------------------------------------------------------------
// Definition of the three-dimensional frame uniform element load class:
//
class cFrame3DUnifLoad : public cElemLoad
{
 protected:
  int     Local;            // Flag for loads in the local system (0/1)
  cVector LoadVal;
  
 public:
                 cFrame3DUnifLoad(void);
  virtual       ~cFrame3DUnifLoad(void);
          void   ExtForce(double, cVector &);
          void   Read(void);
          void   GetVal(sNatCoord, sNodeCoord *, double *);
};

// -------------------------------------------------------------------------
// Definition of the cylindrical shell uniform load class:
//
class cDonnellUnifLoad : public cElemLoad
{
 protected:
  int     Cartesian;        // Flag for loads in the cartesian system (0/1)
  cVector LoadVal;

 public:
                 cDonnellUnifLoad(void);
  virtual       ~cDonnellUnifLoad(void);
          void   Read(void);
          void   GetVal(sNatCoord, sNodeCoord *, double *);
};

#endif
