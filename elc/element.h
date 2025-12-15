// -------------------------------------------------------------------------
// element.h - file containing the definition of the Element class.
// -------------------------------------------------------------------------
//
// The Element class defines the general behavior of a finite element in
// the analysis program and implements the methods which are common for
// the different element types. The specific features of each element are
// defined by a set of pure virtual methods implemented in the derived
// classes. The Element class is also responsible for the storage of all
// elements of the mesh.
//
// Element
// |-- PlTruss
// |   |-- PlTrussTL
// |   |-- PlTrussCR
// |-- Truss3D
// |   |-- Truss3DTL
// |   |-- Truss3DCR
// |-- PlFrame (see elmplfrm.h for complete tree)
// |   |-- PlFrameTL
// |   |-- PlFrameCR
// |   |-- Grid
// |-- Frame3D
// |   |-- Frame3DTL
// |   |-- Frame3DCR
// |   |-- Frame3DCRTL
// |-- Parametric
// |   |-- ParametricTL
// |   |-- Interface
//
// -------------------------------------------------------------------------
// Static public methods:
// -------------------------------------------------------------------------
//
// static void ReadGeoStat(void)
//
// This method reads the parameters of the geostatic stress state.
// -------------------------------------------------------------------------
//
// static void ReadIntOrd(void)
//
// This method reads the different integration orders.
// -------------------------------------------------------------------------
//
// static void ReadNumElm(void)
//
// This method reads the total number of finite elements.
// -------------------------------------------------------------------------
//
// static void ReadBeamOrt(void)
//
// This method reads the beam orientation vectors, corresponding to the
// direction cosines of the cross-section y-axis (largest inertia).
// -------------------------------------------------------------------------
//
// static void ReadPlTruss(void)
// static void ReadPlTrussTL(void)
// static void ReadPlTrussCR(void)
// static void ReadTruss3D(void)
// static void ReadTruss3DTL(void)
// static void ReadTruss3DCR(void)
// static void ReadFrame3D(void)
// static void ReadFrame3DTL(void)
// static void ReadFrame3DCR(void)
// static void ReadFrame3DCRTL(void)
// static void ReadPlTruss(void)
// static void ReadPlTrussTL(void)
// static void ReadPlTrussCR(void)
// static void ReadPlFrame(void)
// static void ReadCoPlFrame(void)
// static void ReadEnCoPlFrame(void)
// static void ReadPlframeCR1(void)
// static void ReadPlframeCR2(void)
// static void ReadPlframeTL1(void)
// static void ReadPlframeTL2(void)
//
// Each of these methods reads and creates bar elements.
// -------------------------------------------------------------------------
//
// static void ReadPlStressQ4(void)
// static void ReadPlStressQ4TL(void)
// static void ReadPlStressQ8(void)
// static void ReadPlStressQ8TL(void)
// static void ReadPlStressQ9(void)
// static void ReadPlStressQ9TL(void)
// static void ReadPlStressT3(void)
// static void ReadPlStressT3TL(void)
// static void ReadPlStressT6(void)
// static void ReadPlStressT6TL(void)
// static void ReadPlStrainQ4(void)
// static void ReadPlStrainQ4TL(void)
// static void ReadPlStrainQ8(void)
// static void ReadPlStrainQ8TL(void)
// static void ReadPlStrainQ9(void)
// static void ReadPlStrainQ9TL(void)
// static void ReadPlStrainT3(void)
// static void ReadPlStrainT3TL(void)
// static void ReadPlStrainT6(void)
// static void ReadPlStrainT6TL(void)
// static void ReadAxiSymL6Inf(void)
// static void ReadAxiSymQ4(void)
// static void ReadAxiSymQ4TL(void)
// static void ReadAxiSymQ8(void)
// static void ReadAxiSymQ8TL(void)
// static void ReadAxiSymQ9(void)
// static void ReadAxiSymQ9TL(void)
// static void ReadAxiSymT3(void)
// static void ReadAxiSymT3TL(void)
// static void ReadAxiSymT6(void)
// static void ReadAxiSymT6TL(void)
// static void ReadBrick8(void)
// static void ReadBrick8TL(void)
// static void ReadBrick20(void)
// static void ReadBrick20TL(void)
// static void ReadTet4(void)
// static void ReadTet4TL(void)
// static void ReadTet10(void)
// static void ReadTet10TL(void)
// static void ReadThickPltQ4(void)
// static void ReadHSDTPlateBsp(void);
// static void ReadHSDTPlateBspTL(void);
// static void ReadThickPltQ8(void)
// static void ReadShallowShellQ8(void)
// static void ReadShallowShellQ8TL(void)
// static void ReadShellQ8(void)
// static void ReadDonnellShellQ8(void)
// static void ReadDonnellShellQ8TL(void)
// static void ReadDonnellShellBsp(void);
// static void ReadDonnellShellBspTL(void);
//
// Each of these methods reads and creates parametric elements described
// by the combination of an analysis model (e.g. plane stress) and a shape
// (e.g. Q8).
// -------------------------------------------------------------------------
//
// static void ReadInterfL2(void)
// static void ReadInterfL3(void)
// static void ReadInterfQ4(void)
// static void ReadInterfQ8(void)
//
// Each of these methods reads and creates interface parametric elements
// described by the combination of an analysis model (e.g. interface 3D) and
// a shape (e.g. Q8).
// -------------------------------------------------------------------------
//
// static void Destroy(void)
//
// This method destroys all the elements of the current model. It should be
// called only at the end of the program.
// -------------------------------------------------------------------------
//
// static int GetLargeRot3D(void)
//
// This method retuns (1) for elements able to treat large rotation in the
// space and (0) otherwise. It can be used in the control algorithms.
// -------------------------------------------------------------------------
//
// static int GetMaxNode(void)
//
// This method returns the maximum number of element nodes for the current
// mesh. It can be used in the allocation of temporary vectors storing
// element nodal coordinates or similar quantities.
// -------------------------------------------------------------------------
//
// static int GetMaxDof(void)
//
// This method returns the maximum number of element dofs for the current
// mesh. It can be used in the allocation of temporary arrays.
// -------------------------------------------------------------------------
//
// static int GetMaxIntPnt(void)
//
// This method returns the maximum number of integration points of an
// element in the current mesh. It can be used in the allocation of
// temporary arrays.
// -------------------------------------------------------------------------
//
// static int GetNumSclLab(void)
//
// This method returns the total number of scalar response labels (e.g.
// FORCE_X, STRESS_XY, ...) which should be printed in the output file.
// -------------------------------------------------------------------------
//
// static int GetNumElm(void)
//
// This method returns the total number of elements of the current model.
// -------------------------------------------------------------------------
//
// static cElement *GetElm(int i)
//
//   i - element index                                                (in)
//
// This method returns a pointer to the element of index i.
// -------------------------------------------------------------------------
//
// static cElement *FindElm(int label)
//
//   label - element id                                               (in)
//
// This method returns a pointer to the element associated with the given
// label.
// -------------------------------------------------------------------------
//
// static void GetVecSclLab(int *label)
//
//   label - scalar response labels                                  (out)
//
// This method stores (in the given vector) the indices of the scalar
// response labels of the current model.
//
// -------------------------------------------------------------------------
// Shared methods (implemented only in the base class):
// -------------------------------------------------------------------------
//
// eElmType GetType(void)
//
// This method returns the element type.
// -------------------------------------------------------------------------
//
// int GetLabel(void)
//
// This method returns the element label (id).
// -------------------------------------------------------------------------
//
// cShape *GetShape(void)
//
// This method returns a pointer to the element shape.
// -------------------------------------------------------------------------
//
// int GetNumNode(void)
//
// This method returns the number of element nodes.
// -------------------------------------------------------------------------
//
// int GetNumElmDof(void)
//
// This method returns the number of element dofs.
// -------------------------------------------------------------------------
//
// void GetElmDofs(int *dof)
//
//   dof - vector of element dofs                                      (in)
//
// This method returns the element dofs, i.e. the equations associated
// with the element nodes. It is important to note that dof[i] == 0
// indicates that this dof is inactive.
// -------------------------------------------------------------------------
//
// void ActivateDof(void)
//
// This method activates the dofs of the element nodes. It is important to
// note that only the nodal dofs associated with the element formulation
// are activated. Thus, for a node linked only to plane truss elements,
// only dof[0] (u) and dof[1] (v) are activated.
// -------------------------------------------------------------------------
//
// void AddGlobVec(cVector &elmvec, cVector &glbvec)
//
//   elmvec - element vector                                           (in)
//   glbvec - global vector                                        (in/out)
//
// This method adds the given element vector at the appropriate places of
// the given global vector. It is used in the assembly of the global
// internal force vector.
// -------------------------------------------------------------------------
//
// void AddGlobMat(cMatrix &elmmat, cSysMatrix &glbmat)
//
//   elmmat - element matrix                                           (in)
//   glbmat - global matrix                                        (in/out)
//
// This method adds the given element matrix at the appropriate places of
// the given global matrix. It is used in the assembly of the global
// stiffness matrix.
// -------------------------------------------------------------------------
//
// void GlobToElm(cVector &glbvec, cVector &elmvec)
//
//   glbvec - global vector                                            (in)
//   elmvec - element vector                                          (out)
//
// This method extracts the displacements of the element nodes from the
// given global nodal displacement vector. It assigns 0.0 to the fixed
// displacements.
// -------------------------------------------------------------------------
//
// void NodalDispl(cVector &u)
//
//   u - nodal displacements                                          (out)
//
// This method returns the element nodal displacements.
// -------------------------------------------------------------------------
//
// bool HasSupp(void)
//
// This method returns true if the element has a supported node and false
// otherwise.
// -------------------------------------------------------------------------
//
// void AddReactions(cVector &gelm, cMatrix &R)
//
//   gelm - element internal force vector                              (in)
//   R    - nodal reactions                                        (in/out)
//
// This method adds the element internal forces to the nodal support
// reactions: R = g - f, where R is a matrix (nnode x 6).
//
// -------------------------------------------------------------------------
// Virtual methods:
// -------------------------------------------------------------------------
//
// virtual void CalcRotMat(void)
//
// This method computes the rotation matrix (3x3), which elements are the
// direct cosines of the undeformed configuration.
// -------------------------------------------------------------------------
//
// virtual void UpdateItera(cVector &du)
//
// This method updates the element internal variables. It should be called
// after the correction of each nonlinear analysis iteration.
// -------------------------------------------------------------------------
//
// virtual void UpdateState(void)
//
// This method updates the element internal variables. It should the called
// after the convergence of each nonlinear analysis step.
// -------------------------------------------------------------------------
//
// virtual cAnModel *GetAnModel(void)
//
// This method returns a pointer to the element analysis model or 0 if the
// element does not have one.
// -------------------------------------------------------------------------
//
// virtual int GetIntOrd(void)
//
// This method returns the index of the element integration order.
// -------------------------------------------------------------------------
//
// virtual int GetNumIntPnt(void)
//
// This method returns the total number of integration points used in the
// integration of the stiffness matrix.
// -------------------------------------------------------------------------
//
// virtual void IntPntTemp(cVector &Temp)
//
//  Tempe -  Temperature                                             (out)
//
// This method computes the temperatures at the integration points
// for the current temperature field. For a generic element.
// -------------------------------------------------------------------------
//
// virtual void IntPntStress(cMatrix &Strain, cMatrix &Stress)
//
//  Strain -  strains                                                (out)
//  Stress -  stresses                                               (out)
//
// This method computes the strains and stresses at the integration points
// for the current nodal displacements. For a generic element S[i][j] the
// index "i" is related to the integraton point and the index "j" is related
// to the stress/strain component.
// -------------------------------------------------------------------------
//
// virtual void PntStress(sNatCoord *Pnts, int NumPnt, cMatrix &Strain,
//                        cMatrix &Stress, bool opt = false)
//
//  Pnts     -  input points in matrix format (nx3)                     (in)
//  Strain   -  strains                                                (out)
//  Stress   -  stresses                                               (out)
//  opt      -  optimization flag                                       (in)
//
// This method computes the strains and stresses at given parametric location
// for current nodal displacements. For a generic element S[i][j] the
// index "i" is related to the input point and the index "j" is related
// to the stress/strain component. In stress computation of each point, the
// section analysis object of the nearest integration point is considered.
// When the optimization flag is true, the program reuses the same
// section analysis indexes computed in the last call, if the element shape and
// integration quadrature do not change.
// -------------------------------------------------------------------------
//
// virtual void InitialStresses(sNatCoord p, cVector &str)
//
//  p      -  parametric coords of the integration point              (in)
//  str    -  stresses                                               (out)
//
// This method computes the initial stresses at a given integration point.
// For geostatic problems the initial stress field depend on the depth of
// point (amount of soil cover).
// -------------------------------------------------------------------------
//
// virtual double VolCoeff(double t,int nn,double *mapfnc,sNodeCoord *coord)
//
//  t      -  thickness                                               (in)
//  nn     -  number of element nodes                                 (in)
//  mapfnc -  mapping functions                                       (in)
//  coord  -  nodal coordinates                                       (in)
//
// This method returns the coefficient that multiply the Jacobian
// determinant (|J|) in order to obtain the differential volume (dV) used
// in the integration of the element matrices. Some examples:
// plane stress  dV = t*|J| drds,  it returns the thickness t.
// axisymemtric  dV = 2*PI*R*|J| drds, it returns 2*PI*R.
// solids        dV = |J| drdsdt, it returns 1.
// -------------------------------------------------------------------------
//
// virtual int GetNumElmDof(void)
//
// This method returns the number of element degrees of freedom.
// -------------------------------------------------------------------------
//
// virtual void GetElmDofs(int *dof)
//
//  dof    -  degrees of freedom                                     (out)
//
// This method returns the element degrees of freedom.
//
// -------------------------------------------------------------------------
// Pure virtual methods:
// -------------------------------------------------------------------------
//
// virtual void Read(void)
//
// This method should read the element attributes (material, section, ...)
// and the element nodes.
// -------------------------------------------------------------------------
//
// virtual int GetNumDofNode(void)
//
// This method returns the number of degrees of freedom of each element
// node.
// -------------------------------------------------------------------------
//
// virtual int GetActDir(int *dir)
//
//  dir    -  array of the active directions (dofs)                  (out)
//
// This method returns the active directions (degrees of freedom) of the
// element nodes. Thus, for a space truss element the active directions
// are dir[0] = dir[1] = dir[2] = 1, and dir[3] = dir[4] = dir[5] = 0,
// since only the displacements u,w,w are active.
// -------------------------------------------------------------------------
//
// virtual int GetNumStrCmp(void)
//
// This method returns the number of element stress components.
// -------------------------------------------------------------------------
//
// virtual void GetStrLabels(int *label)
//
//  label  -  array of response labels                               (out)
//
// This method returns the indices of the element stress component labels
// in the global vector VEC_SCL_LAB.
// -------------------------------------------------------------------------
//
// virtual int IntForce(cVector &g)
//
//  g      -  element internal force vector                          (out)
//
// This method computes the element internal force vector for the current
// nodal displacements. It returns (1) for succesfull computation and (0)
// otherwise.
// -------------------------------------------------------------------------
//
// virtual void StiffMat(cMatrix &K)
//
//  K      -  element tangent stiffness matrix                       (out)
//
// This method evaluates the tangent stiffness matrix for the current nodal
// displacements.
// -------------------------------------------------------------------------
//
// virtual void MassMat(cMatrix &M)
//
//  M      -  element mass matrix                                    (out)
//
// This method evaluates the element mass matrix.
// -------------------------------------------------------------------------
//
// virtual void GeomStiff(cMatrix &G)
//
//  G      -  element geometric stiffness matrix                     (out)
//
// This method evaluates the element geometric stiffness matrix used in
// linearized (eigenvalue) stability analysis.
// -------------------------------------------------------------------------
//
// virtual void NodalStress(cMatrix &S)
//
//  S      -  element nodal stresses                                 (out)
//
// This method computes the stresses at each element node for the current
// nodal displacements. For a generic element S[i][j] the index "i"
// is related to the element node and the index "j" is related to the
// stress component.
// -------------------------------------------------------------------------

#ifndef _ELEMENT_H
#define _ELEMENT_H

#include "gbldef.h"
#include "node.h"
#include "shape.h"

// ------------------------------------------------------------------------
// Forward declarations:
//
class cVector;
class cMatrix;
class cSysMatrix;
class cAnModel;
class cSection;
class cIntPoint;

// -------------------------------------------------------------------------
// Auxiliary types:
//
typedef struct
{
  int K[3];   // Order for integration of the stiffness matrix
  int M[3];   // Order for integration of the mass matrix
  int type;   // Integration quadrature type
} sIntOrd;

typedef struct
{
  double top;     // Vertical coordinate of the top of the layer
  double gamma;   // Specific weight
  double k0;      // At-rest lateral earth pressure coefficient
} sGeoStat;

typedef struct
{
  int idort;      // Orientation ID
  double e0y[3];  // Direction cosine of the y-axis {e0} = {ly, my, ny}
}sBeamOrt;

typedef enum
{
  DIRECT_STIFFNESS,
  STATICS
} eOutConv;


// -------------------------------------------------------------------------
// Type definitions:
//
typedef struct _tempdata  // Temperature data
{
 int label;         // Element or Node label
 double tinf;       // Temperature variation in the bottom face
 double tsup;       // Temperature variation in the top face

} sTempData;

inline int sTempDataCmp(const void *e1, const void *e2)
{
  return(((sTempData*)e1)->label - ((sTempData*)e2)->label);
}

// -------------------------------------------------------------------------
// Element types:
//
typedef enum
{
  PLTRUSS,          // Plane truss
  PLTRUSSTL,        // Total Lagrangian plane truss
  PLTRUSSCR,        // Corotational plane truss
  TRUSS3D,          // 3D truss
  TRUSS3DTL,        // Total Lagrangian 3D truss
  TRUSS3DCR,        // Corotational 3D truss
  PLFRAME,          // Plane frame
  PLFRAMENI,        // Plane frame for materially nonlinear analysis (MNA)
  PLFRAMECR1,       // Corotational Plane frame
  PLFRAMECR1NI,     // Corotational Plane frame for MNA
  PLFRAMECR2,       // Enhanced Corotational Plane frame
  PLFRAMECR2NI,     // Enhanced Corotational Plane frame for MNA
  PLFRAMETL1,       // Total Lagrangian Plane frame
  PLFRAMETL2,       // Total Lagrangian Plane frame (shallow arc)
  FRAME3D,          // Linear 3D frame element (small displacements)
  FRAME3DTL,        // Total Lagrangian 3D frame element (moderate rot.)
  FRAME3DCR,        // Co-rotational 3D frame
  FRAME3DCRTL,      // Co-rotational 3D frame with local TL strains
  GRID,             // Grid Element.
  PARAMETRIC,       // Isoparametric element
  PARAMETRICTL,     // Total Lagrangian Isoparametric element
  DEGSHELL2D,       // Degenerated shell 2D.
  INTERFACE         // Interface element
} eElmType;

// -------------------------------------------------------------------------
// Definition of the element base class:
//
class cElement
{
 private:
  static  int        NumElm;       // Number of elements
  static  cElement **VecElm;       // Element array

 protected:
  static  int        NumIntOrd;    // Number of integration orders
  static  sIntOrd   *VecIntOrd;    // Array of integration orders
  static  int        MaxNode;      // Max. number of elm. nodes
  static  int        MaxDof;       // Max. number of elm. dofs
  static  int        MaxIntPnt;    // Max. number of elm. integ. points
  static  int        NumSclLab;    // Number of scalar labels
  static  int       *VecSclLab;    // Vector of scalar labels
  static  int        NumNodSclLab; // Number of node scalar labels
  static  int       *VecNodSclLab; // Vector of node scalar labels
  static  int        NumElemTemp;  // Number of elms with temp. var.
  static  int        IniStrFlag;   // Flag for initial stresses
  static  eAxis      VertAxis;     // Axis for vertical geostatic stresses
  static  int        NumSoil;      // Number of element layers
  static  sGeoStat  *VecSoil;      // Array of soil data for geostatic stresses
  static  eOutConv   OutputConv;   // Output convention for frame forces
  static  int        NumBeamOrt;   // Number of beam orientations
  static  sBeamOrt  *VecBeamOrt;   // Array of the beam orientations
  static  int        LargeRot3D;   // Flag for large rotations in 3D
          eElmType   Type;         // Element type
          int        Label;        // Element label
          int        BeamOrtIdx;   // Beam orientation index (Precisa aqui?)
          cSection  *Section;      // Element section
          cShape    *Shape;        // Element shape (nodes, shp funcs, ...)
 static   sTempData *ElemTemp;    // Element temperatures
          void       UpdGlbVars(void);


 public:
  static  void       ReadGeoVertAxis(void);
  static  void       ReadGeoStat(void);
  static  void       ReadInitialStress(void);
  static  void       ReadOutputConv(void);
  static  void       ReadBeamOrt(void);
  static  void       ReadIntOrd(void);
  static  void       ReadIntType(void);
  static  void       ReadNumElm(void);
  static  void       ReadPlTruss(void);
  static  void       ReadPlTrussTL(void);
  static  void       ReadPlTrussCR(void);
  static  void       ReadTruss3D(void);
  static  void       ReadTruss3DTL(void);
  static  void       ReadTruss3DCR(void);
  static  void       ReadPlFrame(void);
  static  void       ReadPlFrameCR1(void);
  static  void       ReadPlFrameCR2(void);
  static  void       ReadPlFrameTL1(void);
  static  void       ReadPlFrameTL2(void);
  static  void       ReadGrid(void);
  static  void       ReadFrame3D(void);
  static  void       ReadFrame3DTL(void);
  static  void       ReadFrame3DCR(void);
  static  void       ReadFrame3DCRTL(void);
  static  void       ReadPlStressQ4(void);
  static  void       ReadPlStressQ4TL(void);
  static  void       ReadPlStressQ8(void);
  static  void       ReadPlStressQ8TL(void);
  static  void       ReadPlStressQ9(void);
  static  void       ReadPlStressQ9TL(void);
  static  void       ReadPlStressT3(void);
  static  void       ReadPlStressT3TL(void);
  static  void       ReadPlStressT6(void);
  static  void       ReadPlStressT6TL(void);
  static  void       ReadPlStressBezTri(void);
  static  void       ReadPlStressBezTriTL(void);
  static  void       ReadPlStressBezSurf(void);
  static  void       ReadPlStressBezSurfTL(void);
  static  void       ReadPlStressBsp(void);
  static  void       ReadPlStressBspTL(void);
  static  void       ReadPlStrainQ4(void);
  static  void       ReadPlStrainQ4TL(void);
  static  void       ReadPlStrainQ8(void);
  static  void       ReadPlStrainQ8TL(void);
  static  void       ReadPlStrainQ9(void);
  static  void       ReadPlStrainQ9TL(void);
  static  void       ReadPlStrainT3(void);
  static  void       ReadPlStrainT3TL(void);
  static  void       ReadPlStrainT6(void);
  static  void       ReadPlStrainT6TL(void);
  static  void       ReadPlStrainBezTri(void);
  static  void       ReadPlStrainBezTriTL(void);
  static  void       ReadPlStrainBezSurf(void);
  static  void       ReadPlStrainBezSurfTL(void);
  static  void       ReadPlStrainBsp(void);
  static  void       ReadPlStrainBspTL(void);
  static  void       ReadPlThermQ4(void);
  static  void       ReadPlThermQ8(void);
  static  void       ReadPlThermT3(void);
  static  void       ReadPlThermT6(void);
  static  void       ReadPlThermBezTri(void);
  static  void       ReadPlThermBezSurf(void);
  static  void       ReadAxiSymL6Inf(void);
  static  void       ReadAxiSymQ4(void);
  static  void       ReadAxiSymQ4TL(void);
  static  void       ReadAxiSymQ8(void);
  static  void       ReadAxiSymQ8TL(void);
  static  void       ReadAxiSymQ9(void);
  static  void       ReadAxiSymQ9TL(void);
  static  void       ReadAxiSymT3(void);
  static  void       ReadAxiSymT3TL(void);
  static  void       ReadAxiSymT6(void);
  static  void       ReadAxiSymT6TL(void);
  static  void       ReadBrick8(void);
  static  void       ReadBrick8TL(void);
  static  void       ReadBrick20(void);
  static  void       ReadBrick20TL(void);
  static  void       ReadTet4(void);
  static  void       ReadTet4TL(void);
  static  void       ReadTet10(void);
  static  void       ReadTet10TL(void);
  static  void       ReadSolidBsp(void);
  static  void       ReadSolidBspTL(void);
  static  void       ReadThickPltT3(void);
  static  void       ReadThickPltT6(void);
  static  void       ReadThickPltQ4(void);
  static  void       ReadHSDTPlateBsp(void);
  static  void       ReadHSDTPlateBspTL(void);
  static  void       ReadThickPltQ8(void);
  static  void       ReadShallowShellQ8(void);
  static  void       ReadShallowShellQ8TL(void);
  static  void       ReadShallowShellT6(void);
  static  void       ReadShallowShellT6TL(void);
  static  void       ReadShallowShellBezTrian(void);
  static  void       ReadShallowShellBezTrianTL(void);
  static  void       ReadShallowShellBezSurf(void);
  static  void       ReadShallowShellBezSurfTL(void);
  static  void       ReadShallowShellBsp(void);
  static  void       ReadShallowShellBspTL(void);
  static  void       ReadDonnellShellQ8(void);
  static  void       ReadDonnellShellQ8TL(void);
  static  void       ReadDonnellShellBsp(void);
  static  void       ReadDonnellShellBspTL(void);
  static  void       ReadShellT3(void);
  static  void       ReadShellT6(void);
  static  void       ReadShellQ8(void);
  static  void       ReadShellQ8TL(void);
  static  void       ReadShellQ9(void);
  static  void       ReadShellBezTrian(void);
  static  void       ReadShellBezTrianTL(void);
  static  void       ReadShellBez(void);
  static  void       ReadShellBezTL(void);
  static  void       ReadShellBsp(void);
  static  void       ReadShellBspTL(void);
  static  void       ReadInterfL2(void);
  static  void       ReadInterfL3(void);
  static  void       ReadInterfQ4(void);
  static  void       ReadInterfQ8(void);
  static  void       ReadInterfLineBsp(void);
  static  void       ReadElemTemp(void);
  static  void       Destroy(void);
  static  int        GetLargeRot3D(void) { return LargeRot3D; }
  static  int        GetMaxNode(void){ return MaxNode; }
  static  int        GetMaxDof(void) { return MaxDof; }
  static  int        GetMaxIntPnt(void) { return MaxIntPnt; }
  static  int        GetNumSclLab(void) { return NumSclLab; }
  static  int        GetNumNodSclLab(void) { return NumNodSclLab; }
  static  int        GetNumIntOrd(void) { return NumIntOrd; }
  static  int        GetNumElm(void) { return NumElm; }
  static  cElement  *GetElm(int i)   { return VecElm[i]; }
  static  cElement  *FindElm(int);
  static  void       GetVecSclLab(int *);
  static  void       GetVecNodSclLab(int *);
  static  int        ExistThermalLoad  (void);
          sTempData *GetTemp(void);

                     cElement(int);
  virtual           ~cElement(void);
          eElmType   GetType(void) { return Type; }
          int        GetLabel(void) { return Label; }
          cShape    *GetShape(void) { return Shape; }
          int        GetNumNode(void) { return Shape->GetNumNode( ); }
          void       ActivateDof(void);
          void       AddGlobVec(cVector &, cVector &);
          void       AddGlobMat(cMatrix &, cSysMatrix *);
          void       GlobToElm(cVector &, cVector &);
          void       NodalDispl(cVector &);
          void       NodalTemp(cVector &);
          void       NodalDu(cVector &, cVector &);
          bool       HasSupp(void);
          void       AddReactions(cVector &, cMatrix &);
  virtual void       CalcRotMat(cMatrix &) { }
  virtual void       UpdateItera(cVector &) { }
  virtual void       UpdateState(void) { }
  virtual cAnModel  *GetAnModel(void) { return 0; }
  virtual cSection  *GetSection(void) { return Section; }
  virtual int        GetIntOrd(void) { return 0; }
  virtual int        GetIntOrdIdx(void) { return -1; }
  virtual int        GetIntType(void) { return 0; }
  virtual int        GetIntOrdZ(void) { return 0; }  
  virtual int        GetNumIntPnt(void) { return 0; }
  virtual void       GetIntPnt(int,sNatCoord&,double&) { }
  virtual void       IntPntTemp(cVector&) { }
  virtual void       IntPntStress(cMatrix &, cMatrix &) { }
  virtual void       IntPntFlux(cMatrix &) { }
  virtual void       CalcStrInt(double &, cVector &) { } 
  virtual void       PntStress(sNatCoord*,int,cMatrix&,cMatrix&,bool opt = 0) { }
  virtual double     IntFunc(int,cIntPoint*,cVector&,double&){ return 0.0; }
  virtual double     IntFunc(cVector&,double&){ return 0.0; }
  virtual void       ReadIniStr(void) { }
  virtual void       InitialStresses(sNatCoord, cVector &) { }
  virtual double     VolCoeff(double, int, double *, sNodeCoord *) { return 1.0; }
  virtual int        GetNumElmDof(void);
  virtual void       GetElmDofs(int *);
  virtual void       Read(void) = 0;
  virtual int        GetNumDofNode(void) = 0;
  virtual void       GetActDir(int *) = 0;
  virtual int        GetNumStrCmp(void) = 0;
  virtual int        GetNumNodStrCmp(void)     { return GetNumStrCmp( ); }
  virtual void       GetStrLabels(int *) = 0;
  virtual void       GetNodStrLabels(int *l) { GetStrLabels(l); }
  virtual int        IntForce(cVector &) = 0;
  virtual void       StiffMat(cMatrix &) = 0;
  virtual void       ConducMat(cMatrix &){ }
  virtual void       MassMat(cMatrix &) = 0;
  virtual void       GeomStiff(cMatrix &) = 0;
  virtual void       NodalStress(cMatrix &) = 0;
};


#endif
