// -------------------------------------------------------------------------
// heat.h - file containing the definition of the cHeat class.
// -------------------------------------------------------------------------
//
// The cLoad class defines the general behavior of a heat load in a
// finite element system and implements the methods which are common for the
// different heat generation types. The specific features of each heat are defined
// by a set of pure virtual methods implemented in the derived classes.
// The cHeat class is also responsible for the storage of all external heat
// generation of the current mesh.
//
// Each heat generation can vary in time (t) and space (x). In order to handle these
// features each heat vector f is computed as:
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
// |-- cNodalSource
// |-- cElemHeat
// |   |--cBodySource
// |   |--cEdgeUnifFlux
// |   |--cFaceUnifFlux
//
// Each heat object has a pointer to its parent entity (node or element).
// Moreover, for the element nodes also has an associated shape object which
// is automatically created by the program based on the parent element and
// on the heat generation type: edge, face or body source. Using these variables are
// used in the computation of the equivalent nodal heat and its summation
// in the global vector.
//
// -------------------------------------------------------------------------
// Static public methods:
// -------------------------------------------------------------------------
//
// static void ReadNodalSources(void)
//
// This method creates and reads the nodal sources.
// -------------------------------------------------------------------------
//
// static void ReadBodySources(void)
//
// This method creates and reads the body sources.
// -------------------------------------------------------------------------
//
// static void ReadLineUnifFluxes(void)
//
// This method creates and reads the uniforme edge (line) fluxes.
// -------------------------------------------------------------------------
//
// static void ReadFaceUnifFluxes(void)
//
// This method creates and reads the uniform face (area) fluxes.
// -------------------------------------------------------------------------
//
// static void Destroy(void)
//
// This method destroys all stored external heat generators.
// -------------------------------------------------------------------------
//
// static int GetMaxDof(void)
//
// This method returns the maximum number of flux dofs.
// -------------------------------------------------------------------------
//
// static cLoad *GetHead(void)
//
// This method returns the head of the flux list.
//
// -------------------------------------------------------------------------
// Shared methods (implemented only in the base class):
// -------------------------------------------------------------------------
//
// cLoad *GetNext(void)
//
// This method returns next heat generator in the list.
//
// -------------------------------------------------------------------------
// Pure virtual methods:
// -------------------------------------------------------------------------
//
// virtual void ExtHeat(double t, cVector &f)
//
//  t      -  time                                                    (in)
//  f      -  equivalent heat vector                                  (out)
//
// This method computes the nodal heat vector equivalent to the applied
// heat at the given time.
// -------------------------------------------------------------------------
//
// virtual void AddGlobVec(cVector &felm, cVector &fglob)
//
//  felm   -  element heat vector                                     (in)
//  fglob  -  global heat vector                                  (in/out)
//
// This method adds the given element heat vector to the global load vector.
// -------------------------------------------------------------------------

#ifndef _HEATFLUX_H
#define _HEATFLUX_H

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
  NODAL_SOURCE,     // Nodal source.
  BODY_SOURCE,      // Body (domain) source in parametric elements
  LINE_FLUX,        // Line (edge) flux in parametric elements
  FACE_FLUX         // Face (area) flux in parametric elements
} eHeatType;

// -------------------------------------------------------------------------
// Definition of the heat load class:
//
class cHeat
{
 private:
  static  cHeat  *Head;
  static  cHeat  *Tail;

 protected:
          cHeat     *Next;
          cTimeFunc *Func;
          eHeatType  Type;

 public:
  static  void   ReadNodalSources(void);
  static  void   ReadBodySources(void);
  static  void   ReadLineUnifFluxes(void);
  static  void   ReadLineGenFluxes(void);
  static  void   ReadBodyGenSources(void);
  //static  void       ReadFaceUnifFluxes(void);
  static  void   ReadFieldBodySources(void);
  static  void   ReadLineFieldFluxes(void);
  static  void   EqvHeat(cElement*,double,cVector&);
  static  void   Destroy(void);
  static  int    GetMaxDof(void) { return 500; } //TODO discutir com evandro
  static  cHeat *GetHead(void) { return Head; }
                 cHeat(void);
  virtual       ~cHeat(void);
          cHeat *GetNext(void) { return Next; }
  virtual void   ExtHeat(double,cVector&) = 0;
  virtual void   AddGlobVec(cVector&,cVector&) = 0;
};


// -------------------------------------------------------------------------
// Definition of the nodal load class:
//
class cNodalSource : public cHeat
{
 protected:
  cNode  *Node;
  double  Value;

 public:
                 cNodalSource(cNode *,double);
  virtual       ~cNodalSource(void);
          void   ExtHeat(double, cVector &);
          void   AddGlobVec(cVector &, cVector &);
};


// -------------------------------------------------------------------------
// Definition of the element flux class:
//
class cElemHeat : public cHeat
{
 public:                // TEMPORARIO
  cElement  *Elem;      // Parent element
  cShape    *Shape;     // Load shape

 public:
                 cElemHeat(void);
  virtual       ~cElemHeat(void);
  virtual void   ExtHeat(double, cVector &);
  virtual void   AddGlobVec(cVector &, cVector &);
  virtual void   Read(void) = 0;
  virtual double GetVal(sNatCoord, sNodeCoord *) = 0;

 protected:
          int    GetNumDofs(void);
          void   GetDofs(int *);
};

// -------------------------------------------------------------------------
// Definition of the body force class:
//
class cBodySource : public cElemHeat
{
 protected:
   double Value;

 public:
                  cBodySource(void);
  virtual        ~cBodySource(void);
          void    Read(void);
          double  GetVal(sNatCoord, sNodeCoord*);
};

// -------------------------------------------------------------------------
// Definition of the edge uniform force class:
//
class cEdgeUnifFlux : public cElemHeat
{
 protected:
  double  FluxVal;             // Constant heat flux.

 public:
                 cEdgeUnifFlux(void);
  virtual       ~cEdgeUnifFlux(void);
  virtual void   Read(void);
  virtual double GetVal(sNatCoord, sNodeCoord *);
};

// -------------------------------------------------------------------------
// Definition of the edge general flux class:
//
class cEdgeGenFlux : public cElemHeat
{
 protected:
  sExpEval expr;      // Expression evaluator.

  double x, y, z;     // Input/Output variables used
  double f;           // in the expression evaluator.

 public:
                 cEdgeGenFlux(void);
  virtual       ~cEdgeGenFlux(void);
  void   Read(void);
  double GetVal(sNatCoord, sNodeCoord *);
};

// -------------------------------------------------------------------------
// Definition of the general body heat source class:
//
class cBodyGenSource : public cElemHeat
{
 protected:
  sExpEval expr;      // Expression evaluator.

  double x, y, z;     // Input/Output variables used
  double f;           // in the expression evaluator.

 public:
                 cBodyGenSource(void);
  virtual       ~cBodyGenSource(void);
  void   Read(void);
  double GetVal(sNatCoord, sNodeCoord *);
};


// -------------------------------------------------------------------------
// Definition of the body source field class:
//
class cBodySrcField: public cElemHeat
{
 public:
                  cBodySrcField(void);
  virtual        ~cBodySrcField(void);
          void    Read(void);
          double  GetVal(sNatCoord, sNodeCoord*);
};

// -------------------------------------------------------------------------
// Definition of the edge field force class:
//
class cEdgeFieldFlux : public cElemHeat
{
 public:
                 cEdgeFieldFlux(void);
  virtual       ~cEdgeFieldFlux(void);
  virtual void   Read(void);
  virtual double GetVal(sNatCoord, sNodeCoord *);
};

// -------------------------------------------------------------------------
// Definition of the face uniform force class:
//
/*
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
*/

// -------------------------------------------------------------------------
// Definition of the face uniform pressure force class:
//
/*
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
*/

#endif
