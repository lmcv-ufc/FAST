// -------------------------------------------------------------------------
// node.cpp - implementation of node class.
// -------------------------------------------------------------------------
// Created:      10-Nov-2000     Evandro Parente Junior
//
// Modified:     02-Dec-2000     Evandro Parente Junior
//               Implementation of nodal springs.
//
// Modified:     08-Aug-2011     Evandro Parente Junior
//               Implementation of nodal constraints (master-slave).
//
// Modified:     08-Apr-2015     Elias Saraiva Barroso
//               Input/output operations using c++ streams.
//
// Modified:     04-Aug-2016     Elias Saraiva Barroso
//               Creation of control point weights data structure and methods
//               in cNode class.
//
// Modified:     23-Jul-2024     Elias Saraiva Barroso
//               Implementation of MPCs.
//
// Modified:     31-Oct-2025     Evandro Parente Junior
//               Bug correction in GetNodeIndexInSolverOrder affecting nodal
//               constraints.
// -------------------------------------------------------------------------

#include <iostream>
#include <stdlib.h>

using namespace std;

#include "node.h"
#include "gblvar.h"

// -------------------------------------------------------------------------
// Static variables:
//
int         cNode :: NumNode   = 0;
cNode**     cNode :: VecNode   = 0;
int         cNode :: NumW      = 0;         
double*     cNode :: VecW      = 0;         
int         cNode :: NumThk    = 0;      
double*     cNode :: VecThk    = 0;      
double      cNode :: DefThk    = 1;      
int         cNode :: NumNormal = 0;      
cVector*    cNode :: VecNormal = 0;      
sOptOrder*  cNode :: VecSolOrd = 0;
int         cNode :: NumSupp   = 0;
sNodeSupp*  cNode :: VecSupp   = 0;
int         cNode :: NumMass   = 0;
sNodeMass*  cNode :: VecMass   = 0;
mNodeConstr cNode :: MapNodeConstr;
vNodeMPC    cNode :: VecNodeMPC;

// -------------------------------------------------------------------------
// Public methods:
//

// ============================== ReadNumNode ==============================

void cNode :: ReadNumNode(void)
{
  // Read the number of mesh nodes.

  if (!(in >> NumNode) || (NumNode <= 0))
  {
    cout << "Error in the input of the number of nodes!\n";
    exit(0);
  }

  // Alloc the array of nodes.

  VecNode = new cNode*[NumNode];
}

// =============================== ReadCoord ===============================

void cNode :: ReadCoord(void)
{
  // Read the number of mesh nodes.

  int n;
  if (!(in >> n) || (NumNode != n))
  {
    cout << "Error in the input of the number of nodal coordinates!\n";
    exit(0);
  }

  // Create and read each node.

  int id;
  double x,y,z;
  for (int i = 0; i < NumNode; i++)
  {
    if (!(in >> id) || !(in >> x) || !(in >> y) || !(in >> z))
    {
      cout << "Error in the input of the coordinates of node " << i+1 << "!\n";
      exit(0);
    }
    VecNode[i] = new cNode(i+1, x, y, z);
  }
}

// =============================== ReadCtrlPntW ============================

void cNode :: ReadCtrlPntW(void)
{
  // Check if weights are instantiated.

  if (NumW != 0 || NumNode == 0)
  {
    cout << "Error in the input of the control points weights!\n";
    exit(0);
  }

  // Create the vector of control points weights and set its components to 1.
  
  NumW = NumNode;
  VecW = new double [NumW];
  for(int i = 0; i < NumW; ++i)
    VecW[i] = 1.0;

  // Read the number of input weights

  int n;
  if (!(in >> n) || (n < 1) || (n > NumNode))
  {
    cout << "Error in the input of the number of control points weights!\n";
    exit(0);
  }

  // Read each weight.

  int    id; // Ctrl point id;
  for(int i = 0; i < n; ++i)
  {
    if (!(in >> id) || (id < 1) || (id > NumNode))
    {
      cout << "Error in the input of control point index!\n";
      exit(0);
    }

    if (!(in >> VecW[id-1]))
    {
      cout << "Error in the input of control point weight!\n";
      exit(0);
    }
  }
}

// =============================== ReadDefThk ==============================

void cNode :: ReadDefThk(void)
{
  // Read default node thickness.

  if (!(in >> DefThk))
  {
    cout << "Error in the input of the default node thickness!\n";
    exit(0);
  }
}

// =============================== ReadThk =================================

void cNode :: ReadThk(void)
{
  // Read the number of input node thickness.

  if (!(in >> NumThk) || (NumThk < 1) || ( NumThk > NumNode))
  {
    cout << "Error in the input of the number of node thickness!\n";
    exit(0);
  }
  VecThk = new double [NumThk];

  // Read each thickness.

  int    id;
  for(int i = 0; i < NumThk; ++i)
  {
    if (!(in >> id) || (id < 1) || (id > NumNode))
    {
      cout << "Error in the input of node thickness index!\n";
      exit(0);
    }

    if (!(in >> VecThk[id-1]))
    {
      cout << "Error in the input of node thickness!\n";
      exit(0);
    }
  }
}

// =============================== ReadNodeNormal ==========================

void cNode :: ReadNodeNormal(void)
{
  // Check if node normal vector exist.

  if (!NumNormal)
  {
    NumNormal = NumNode;
    VecNormal = new cVector [NumNode];
  } 

  // Read the number of input node normal.

  int n;
  if (!(in >> n) || (n < 1) || (n > NumNode))
  {
    cout << "Error in the input of the number of node axis V3!\n";
    exit(0);
  }

  // Read each vector.

  int    id;
  double v1,v2,v3;
  for(int i = 0; i < n; ++i)
  {
    if (!(in >> id) || (id < 1) || (id > NumNode))
    {
      cout << "Error in the input of node axis V3!\n";
      exit(0);
    }

    if (!(in >> v1) || !(in >> v2) || !(in >> v3))
    {
      cout << "Error in the input of node axis V3!\n";
      exit(0);
    }
    VecNormal[id-1].Resize(3);
    VecNormal[id-1][0] = v1;
    VecNormal[id-1][1] = v2;
    VecNormal[id-1][2] = v3;
  }
}

// =============================== ReadSolverOrder ================================

static int CompOpt(const void *a, const void *b)
{
  return(((sOptOrder *)a)->optim - ((sOptOrder *)b)->optim);
}

void cNode :: ReadSolverOrder(void)
{
  // Alloc the number of nodes.

  int n;
  if (!(in >> n) || (n != NumNode))
  {
    cout << "Error in the input of the number of nodes in solver order!\n";
    exit(0);
  }
  VecSolOrd = new sOptOrder[NumNode];

  // Read each node.

  int opt;
  for (int i = 0; i < NumNode; i++)
  {
    if (!(in >> opt))
    {
      cout << "Error in the input of nodes in solver order!\n";
      exit(0);
    }

    if (opt < 1 || opt > NumNode)
    {
      cout << "Invalid node number " << opt << " in solver order)!\n";
      exit(0);
    }
    VecSolOrd[i].optim = opt;
    VecSolOrd[i].node  = VecNode[i];
  }

  // Sort the vector in ascending order.

  qsort(VecSolOrd, NumNode, sizeof(sOptOrder), CompOpt);
//  for (int i = 0; i < NumNode; i++)
//    cout << "Opt(%d) = %d\n", VecSolOrd[i].optim, VecSolOrd[i].node->Label;
}

// =============================== ReadSupp ================================

void cNode :: ReadSupp(void)
{
  // Read the number of mesh supports.

  if (!(in >> NumSupp) || (NumSupp < 1))
  {
    cout << "Error in the input of the number of supports!\n";
    exit(0);
  }

  // Alloc the array of supports.

  VecSupp = new sNodeSupp[NumSupp];

  // Read each support.

  int id,dx,dy,dz,rx,ry,rz;
  for (int i = 0; i < NumSupp; i++)
  {
    if (!(in >> id) || !(in >> dx) || !(in >> dy) || !(in >> dz) ||
        !(in >> rx) || !(in >> ry) || !(in >> rz))
    {
      cout << "Error in the input of the support " << i+1 << "!\n";
      exit(0);
    }

    VecSupp[i].node = FindNode(id);
    if (!VecSupp[i].node)
    {
      cout << "Error in the input of the support " << i+1 << " (invalid node)!\n";
      exit(0);
    }

    VecSupp[i].node->SetSupp( );
    VecSupp[i].dx = dx;
    VecSupp[i].dy = dy;
    VecSupp[i].dz = dz;
    VecSupp[i].rx = rx;
    VecSupp[i].ry = ry;
    VecSupp[i].rz = rz;
  }
}

// ============================= ReadNodeMass ==============================

void cNode :: ReadNodeMass(void)
{
  // Read the number of nodal masses.

  if (!(in >> NumMass) || (NumMass < 1))
  {
    cout << "Error in the input of the number of nodal masses!\n";
    exit(0);
  }

  // Alloc the array of nodal masses.

  VecMass = new sNodeMass[NumMass];

  // Read each mass.

  int id;
  double mx,my,mz,mrx,mry,mrz;
  for (int i = 0; i < NumMass; i++)
  {
    if (!(in >> id) || !(in >> mx) || !(in >> my) || !(in >> mz) ||
        !(in >> mrx) || !(in >> mry) || !(in >> mrz))
    {
      cout << "Error in the input of the nodal mass " << i+1 << "!\n";
      exit(0);
    }

    VecMass[i].node = FindNode(id);
    if (!VecMass[i].node)
    {
      cout << "Error in the input of the nodal mass " << i+1 << "(invalid node)!\n";
      exit(0);
    }

    VecMass[i].mx  = mx;
    VecMass[i].my  = my;
    VecMass[i].mz  = mz;
    VecMass[i].mrx = mrx;
    VecMass[i].mry = mry;
    VecMass[i].mrz = mrz;
  }
}

// ========================== ReadNodeConstraint ===========================

void cNode :: ReadNodeConstraint(void)
{
  // Read the number of node constraints.

  int n;
  if (!(in >> n) || (n < 1))
  {
    cout << "Error in the input of the number of node constraints!\n";
    exit(0);
  }

  // Read each constraint.

  int no1,no2;
  for (int i = 0; i < n; i++)
  {
    if (!(in >> no1) || !(in >> no2))
    {
      cout << "Error in the input of the node constraint "<< i+1 << "!\n";
      exit(0);
    }

    if ((no1 < 1) || (no2 < 1) || (no1 == no2))
    {
      cout << "Error in the input of the node constrant " << i+1 << " (invalid nodes)!\n";
      exit(0);
    }

    // Get node pointers.

    cNode *node1 = FindNode(no1);
    cNode *node2 = FindNode(no2);
    if (!node1 || !node2)
    {
      cout << "Error in the input of the node constrant " << i+1 << " (invalid node)!\n";
      exit(0);
    }

    // Read linked dofs.

    int dof[6];
    if (!(in >> dof[0]) || !(in >> dof[1]) || !(in >> dof[2]) ||
        !(in >> dof[3]) || !(in >> dof[4]) || !(in >> dof[5]))
    {
      cout << "Error in the input of the node constrant " << i+1 << "(invalid dofs)!\n";
      exit(0);
    }

    // Choose the master/slave pair. The master is always the first node in
    // the solver order (standard or optimized).

    int slave = no1;
    cNode *master = node2;
    if (GetNodeIndexInSolverOrder(no2) > GetNodeIndexInSolverOrder(no1))
    {
      slave  = no2;
      master = node1;
    }
    cout << "Node Constraint: master " << master->GetLabel( ) << " slave " << slave << "\n";

    // Assembly the constraint data.

    sNodeConstr *constraint = new sNodeConstr;
    constraint->master = master;
    for (int j = 0; j < 6; j++)
    {
      if (dof[j])
        constraint->doflink[j] = true;
      else
        constraint->doflink[j] = false;
    }

    // Store the constraint at the map structure.

    MapNodeConstr[slave] = constraint;
  }
}

// ============================== ReadNodeMPC ==============================

void cNode :: ReadNodeMPC(void)
{
  // Read the number of constraint equations.

  int n;
  if (!(in >> n) || (n < 1))
  {
    cout << "Error in the input of the number of linear multi-point constraint equations!\n";
    exit(0);
  }
  VecNodeMPC.resize(n);

  // Read each equation data.

  for (int i = 0; i < n; i++)
  {
    // Read the number of terms.

    int nt;
    if (!(in >> nt) || (nt < 1))
    {
      cout << "Error in the input of the MPC equation number of terms (eq "<< i+1 << ")!\n";
      exit(0);
    }
    VecNodeMPC[i].EqTerm.resize(nt);

    // Read each equation term.

    int nid, ndir;
    double ncoef;
    for (int j = 0; j < nt; j++)
    {
      if (!(in >> nid) || !(in >> ndir) || !(in >> ncoef))
      {
        cout << "Error in the input of the MPC equation term (eq "<< i+1 << ")!\n";
        exit(0);
      }
      VecNodeMPC[i].EqTerm[j].node = nid;
      VecNodeMPC[i].EqTerm[j].dir  = ndir;
      VecNodeMPC[i].EqTerm[j].coef = ncoef;
    }

    // Read the constant term.

    if (!(in >> VecNodeMPC[i].b))
    {
      cout << "Error in the input of the MPC equation constant term (eq "<< i+1 << ")!\n";
      exit(0);
    }
  }
}

// ================================ Destroy ================================

void cNode :: Destroy(void)
{
  // Destroy each node.

  for (int i = 0; i < NumNode; i++) delete VecNode[i];

  // Release the array of nodes.

  delete []VecNode;

  // Release the other arrays: order, supports, springs and masses.

  delete []VecSolOrd;
  delete []VecSupp;
  delete []VecMass;

  // INCLUIR AS CONSTRAINTS!!!
}

// =============================== FindNode ================================

cNode* cNode :: FindNode(int label)
{
  if (label < 1 || label > NumNode) return(0);

  return(VecNode[label-1]);
}

// ========================= GetNodeInSolverOrder ==========================

cNode* cNode :: GetNodeInSolverOrder(int i)
{
  // Return the node according to the optimum ordering.

  if (VecSolOrd) return(VecSolOrd[i].node);

  // Return the node in the standard ordering.

  return(VecNode[i]);
}

// ======================= GetNodeIndexInSolverOrder =======================

int cNode :: GetNodeIndexInSolverOrder(int label)
{
  // Return the node according to the optimum ordering.

  if (VecSolOrd)
    for (int i = 0; i < NumNode; i++)
      if (VecSolOrd[i].node->Label == label) return(i);

  // Return the node in the standard ordering.

  return(label-1);
}

// ============================ FindConstraint =============================

sNodeConstr* cNode :: FindConstraint(int slave)
{
  // Search for the given slave node.

  mNodeConstr::iterator slvit = MapNodeConstr.find(slave);

  // Return the associated constraint (if found).

 // cout << "slave = %d\n", slave;
  if (slvit != MapNodeConstr.end( )) return((*slvit).second);
 // cout << "slave = %d not found \n", slave;

  return(0);
}

// ============================ GetNumMPC ==================================

int cNode :: GetNumMPC( )
{
  return VecNodeMPC.size( );
}

// ============================ GetMPCData =================================

void cNode :: GetMPCData(cMatrix &Q, cVector &b)
{
  Q.Zero( );	
  b.Zero( );	
  for (int i = 0; i < VecNodeMPC.size( ); i++)
  {
    for (int j = 0; j < VecNodeMPC[i].EqTerm.size( ); j++)
    {
      int node = VecNodeMPC[i].EqTerm[j].node - 1;
      int dir  = VecNodeMPC[i].EqTerm[j].dir - 1;
      int dof  = GetNode(node)->GetDof(dir) - 1;

      if (dof < 0)
      {
        printf("Invalid constrained dof: node = %d dir = %d\n", node, dir+1);
        exit(0);
      }

      Q[i][dof] = VecNodeMPC[i].EqTerm[j].coef;
    }
    b[i] = VecNodeMPC[i].b;
  }
  //cout << "[Q] = \n" << scientific;
  //Q.Print( );
  //cout << "{b} = \n";
  //b.Print( );
}

// ================================= cNode =================================

cNode :: cNode(int id, double x, double y, double z)
{
  // Store the node label.

  Label = id;

  // Support flag.

  Supp = false;

  // Store the node coordinates.

  Coord.x = x;
  Coord.y = y;
  Coord.z = z;

  // Initialize the node data.

  for (int i = 0; i < 6; i++) VecDispl[i] = 0.0;
  for (int i = 0; i < 8; i++) VecDof[i]   = 0;
  Temp = 0.0;
}

// ================================ ~cNode =================================

cNode :: ~cNode(void)
{
}

// ================================ GetW ===================================

double cNode :: GetW(void)
{
  if (NumW == 0)
    return 1.0;
  else
    return VecW[Label-1];
}

// ================================ GetThk =================================

double cNode :: GetThk(void)
{
  if (NumThk == 0)
    return DefThk; 
  else
    return VecThk[Label-1];
}

// ================================ GetDofVal ==============================

double cNode :: GetDofVal(int i)
{
  if (i < 6)
    return VecDispl[i];
  else if (i == 7)
    return Temp;
}

// ======================================================= End of file =====
