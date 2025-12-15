// -------------------------------------------------------------------------
// node.h - file containing the definition of node class.
// -------------------------------------------------------------------------
//
// The cNode class defines the general behavior of a finite element node and
// handles a series of attributes (boundary conditions, springs, etc.)
// associated with the nodes of the FE model.
//
// -------------------------------------------------------------------------
// Static public methods:
// -------------------------------------------------------------------------
//
// static void ReadNumNode(void)
//
// This method reads the number of nodes and allocates the memory to store
// the nodes of a FE model.
// -------------------------------------------------------------------------
//
// static void ReadCoord(void)
//
// This method reads the coordinates and creates each node of the FE model.
// -------------------------------------------------------------------------
//
// static void ReadSolverOrder(void)
//
// This method reads the optimal nodal ordering to be used in the assembly
// of the stiffness matrix.
// -------------------------------------------------------------------------
//
// static void ReadSupp(void)
//
// This method reads the nodal supports (boundary conditions) of the FE
// model.
// -------------------------------------------------------------------------
//
// static void ReadNodeMass(void)
//
// This method reads the nodal (concentrated) masses.
// -------------------------------------------------------------------------
//
// static void ReadNodeConstraint(void)
//
// This method reads the nodal constraints (master-slave).
// -------------------------------------------------------------------------
//
// static void Destroy(void)
//
// This method destroys each mesh node and release the allocated memory. It
// shound be called only after the end of the analysis.
// -------------------------------------------------------------------------
//
// static int GetNumSupp(void)
//
// This method returns the number of supports of the FE model.
// -------------------------------------------------------------------------
//
// static sNodeSupp GetSupp(int i)
//
//   i - index of the required support                                 (in)
//
// This method returns the required support from the given index.
// -------------------------------------------------------------------------
//
// static int GetNumNodeMass(void)
//
// This method returns the number of nodal masses of the FE model.
// -------------------------------------------------------------------------
//
// static sNodeMass GetMass(int i)
//
//   i - index of the required nodal mass                              (in)
//
// This method returns the required nodal mass from the given index.
// -------------------------------------------------------------------------
//
// static int GetNumNode(void)
//
// This method returns the number of nodes of the FE model.
// -------------------------------------------------------------------------
//
// static cNode *GetNode(int i)
//
//   i - node index                                                    (in)
//
// This method returns a pointer to the required node from the given index.
// It is generally used to loop through all the model nodes.
// -------------------------------------------------------------------------
//
// static cNode *GetNodeInSolverOrder(int i)
//
//   i - node index                                                    (in)
//
// This method returns a pointer to the required node according to the
// optimal ordering that minimize the memory used to store the global
// stiffness matrix.
// It is generally used to loop through all the model nodes.
// -------------------------------------------------------------------------
//
// static int GetNodeIndexInSolverOrder(int label)
//
//   label - node label (id)                                           (in)
//
// This method returns the node index (position) in the solver order
// (standard or optimized).
// -------------------------------------------------------------------------
// static cNode *FindNode(int label)
//
//   label - node label (id)                                           (in)
//
// This method returns a pointer to the node associated from the given
// label. It is  used when it is necessary to directly access a node whose
// label is known.
// -------------------------------------------------------------------------
//
// static sNodeConstr* FindConstraint(int label)
//
//   label - node label (id)                                           (in)
//
// This method check if the node is linked to a master node and returns the
// constraint pointer. It returns 0 if the node is not linked.
//
// -------------------------------------------------------------------------
// Public methods:
// -------------------------------------------------------------------------
//
// int GetLabel(void)
//
// This method returns the nodal label.
// -------------------------------------------------------------------------
//
// sNodeCoord GetCoord(void)
//
// This method returns the nodal coordinates.
// -------------------------------------------------------------------------
//
// int GetDof(int i)
//
//   i - index of the nodal dof                                        (in)
//
// This method returns the nodal dof corresponding to the given index.
// -------------------------------------------------------------------------
//
// double GetDispl(int i)
//
//   i - index of the nodal displacement                               (in)
//
// This method returns the nodal displacement corresponding to the given
// index.
// -------------------------------------------------------------------------
//
// double GetTemp(void)
//
// This method returns the nodal temperature.
// -------------------------------------------------------------------------
//
// double GetDofVal(int i)
//
//   i - index of the nodal dof                                        (in)
//
// This method returns the nodal dof value corresponding to the given
// index. The following dofs order are assumed:
//
//   0 - displacement u. 
//   1 - displacement v. 
//   2 - displacement w. 
//   3 - rotation in x axis. 
//   4 - rotation in y axis. 
//   5 - rotation in z axis. 
//   6 - warping. 
//   7 - temperature.
//
// -------------------------------------------------------------------------
//
// void SetDof(int i, int d)
//
//   i - index of the nodal dof                                        (in)
//   d - given value                                                   (in)
//
// This method assigns the given value to the dof corresponding to the
// given index. The value 0 indicates a fixed dof, while a positive value
// correspond to the label of the associated equation.
// -------------------------------------------------------------------------
//
// void SetDispl(int i, double d)
//
//   i - index of the nodal dof                                        (in)
//   d - given displacement                                            (in)
//
// This method assigns the given displacement value to the dof
// corresponding to the  given index.
// -------------------------------------------------------------------------
//
// void SetTemp(double d)
//
//   d - given temperature                                             (in)
//
// This method assigns the given nodal temperature.
// -------------------------------------------------------------------------
//
// void SetSupp(void)
//
// This method change the value of the Supp flag to true, indicating that
// the node has a constrained dof (support, spring or prescribed).
// -------------------------------------------------------------------------
//
// bool HasSupp(void)
//
// This method returns the Supp flag, indicating if the node has a
// constrained dof (true) or not (false).
// -------------------------------------------------------------------------

#ifndef _NODE_H
#define _NODE_H

#include <vector>
#include <map>
#include <vec.h>

// -------------------------------------------------------------------------
// Forward declarations:
//
class cNode;

// -------------------------------------------------------------------------
// Auxiliary types:
//
struct sNodeAxes
{
  cVector V1, V2, V3;

  sNodeAxes( ) : V1(3), V2(3), V3(3) { }
};

typedef struct
{
  double     x,y,z;  // Coordinates
  double    *thk;    // Shell thickness
  sNodeAxes *axes;   // Shell updated axes
  sNodeAxes *iaxes;  // Shell initial axes
  cVector   *pa;     // Shell d{p}/dalpha
  cVector   *pb;     // Shell d{p}/dbeta
  cVector   *paa;    // Shell d2{p}/dalpha2
  cVector   *pab;    // Shell d2{p}/dalphadbeta	   
  cVector   *pbb;    // Shell d2{p}/dbeta2
} sNodeCoord;

typedef struct
{
  double xx,yy,xy;
} sNodeDrv;

typedef struct
{
  cNode *node;
  int optim;
} sOptOrder;

typedef struct
{
  cNode *node;
  int dx,dy,dz,rx,ry,rz;
} sNodeSupp;

typedef struct
{
  cNode *node;
  double mx,my,mz,mrx,mry,mrz;
} sNodeMass;

typedef struct
{
  cNode *master;
  bool  doflink[8];
} sNodeConstr;

typedef struct
{
  int    node;
  int    dir;
  double coef;
} sMPCTerm;

typedef struct
{
  std::vector<sMPCTerm> EqTerm;
  double b;
} sNodeMPC;

typedef std::map<int, sNodeConstr*> mNodeConstr;
typedef std::vector<sNodeMPC>       vNodeMPC;

// -------------------------------------------------------------------------
// Definition of node class:
//
class cNode
{
 private:
  static int          NumNode;      // Number of mesh nodes
  static cNode**      VecNode;      // Vector of mesh nodes
  static int          NumW;         // Number of control points weights
  static double*      VecW;         // Vector of control points weights
  static int          NumThk;       // Number of node thickness
  static double*      VecThk;       // Vector of node thickness
  static double       DefThk;       // Default node thickness
  static int          NumNormal;    // Number of node normal vector
  static cVector     *VecNormal;    // Vector of node normal vector
  static sOptOrder*   VecSolOrd;    // Vector of nodes in optimized order
  static int          NumSupp;      // Number of supported nodes
  static sNodeSupp*   VecSupp;      // Vector of supported nodes
  static int          NumMass;      // Number of nodal masses
  static sNodeMass*   VecMass;      // Vector of nodal masses
  static mNodeConstr  MapNodeConstr;// Map with node constraints
  static vNodeMPC     VecNodeMPC;   // Vector of multi-point constraints

 protected:
         int          Label;        // Node label
         bool         Supp;         // Flag for support or presc. displ.
         sNodeCoord   Coord;        // Node coordinates
         int          VecDof[8];    // Nodal dofs
         double       VecDispl[6];  // Nodal displacements
         double       Temp;         // Nodal temperature

 public:
  static void         ReadNumNode(void);
  static void         ReadCoord(void);
  static void         ReadCtrlPntW(void);
  static void         ReadDefThk(void);
  static void         ReadThk(void);
  static void         ReadNodeNormal(void);
  static void         ReadSolverOrder(void);
  static void         ReadSupp(void);
  static void         ReadNodeMass(void);
  static void         ReadNodeConstraint(void);
  static void         ReadNodeMPC(void);  
  static void         ReadSmoothNormal(void);
  static void         Destroy(void);
  static int          GetNumSupp(void)   { return NumSupp;    } 
  static sNodeSupp    GetSupp(int i)     { return VecSupp[i]; }
  static int          GetNumMass(void)   { return NumMass;    }
  static sNodeMass    GetMass(int i)     { return VecMass[i]; }
  static int          GetNumNode(void)   { return NumNode;    }
  static cNode*       GetNode(int i)     { return VecNode[i]; }
  static int          GetNumNormal(void) { return NumNormal;  }
  static cNode*       GetNodeInSolverOrder(int);
  static int          GetNodeIndexInSolverOrder(int);
  static cNode*       FindNode(int);
  static sNodeConstr* FindConstraint(int);
  static int          GetNumMPC(void);
  static void         GetMPCData(cMatrix&, cVector&);    

                      cNode(int, double, double, double);
                     ~cNode(void);
         int          GetLabel(void)  { return Label; }
         sNodeCoord   GetCoord(void)  { return Coord; }
         double       GetW(void);
         double       GetThk(void);
         cVector&     GetNormal(void) { return VecNormal[Label-1]; }
         int          GetDof(int i)   { return VecDof[i]; }
         double       GetDofVal(int i);
         double       GetDispl(int i) { return VecDispl[i]; }
         double       GetTemp(void)   { return Temp; } 
         bool         HasSupp(void)   { return Supp; }
         void         SetSupp(void)   { Supp = true; }
         void         SetDof(int i, int d)      { VecDof[i] = d;   }
         void         SetDispl(int i, double d) { VecDispl[i] = d; }
         void         SetTemp(double t)         { Temp = t;        }
};

#endif
