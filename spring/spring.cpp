// -------------------------------------------------------------------------
// spring.cpp - implementation of the spring class.
// -------------------------------------------------------------------------
// Created:      10-Ago-2011     Evandro Parente Junior
//
// Modified:
// -------------------------------------------------------------------------

#include <iostream>
#include <math.h>
#include <stdlib.h>

using namespace std;

#include "spring.h"
#include "sprprop.h"
#include "node.h"
#include "vec.h"
#include "mat.h"
#include "utl.h"
#include "gblvar.h"
#include "gbldef.h"

// -------------------------------------------------------------------------
// Static variables:
//
int           cSpringSupp :: NumSpring = 0;
cSpringSupp** cSpringSupp :: VecSpring = 0;
int           cSpringConn :: NumSpring = 0;
cSpringConn** cSpringConn :: VecSpring = 0;

// -------------------------------------------------------------------------
// Public methods:
//

// ================================ cSpring ================================

cSpring :: cSpring(void)
{
}

// =============================== ~cSpring ================================

cSpring :: ~cSpring(void)
{
}

// ============================== AddGlobVec ===============================

void cSpring :: AddGlobVec(cVector &sprvec, cVector &glbvec)
{
  // Get spring dofs.

  int ndof = GetNumDofs( );
  int *dof = new int[ndof];
  GetDofs(dof);

  // Add spring vector to the global vector.

  for (int i = 0; i < ndof; i++) if (dof[i])
  {
    glbvec[dof[i]-1] += sprvec[i];
  }

  // Release the dof array.

  delete []dof;
}

// ============================== AddGlobMat ===============================

void cSpring :: AddGlobMat(cMatrix &sprmat, cSysMatrix *glbmat)
{
  // Get spring dofs.

  int ndof = GetNumDofs( );
  int *dof = new int[ndof];
  GetDofs(dof);

  // Add spring matrix to the global matrix.

  for (int i = 0; i < ndof; i++) if (dof[i])
  {
    for (int j = 0; j < ndof; j++) if (dof[j])
    {
      // Add if element is below/at the main diagonal or the matrix is
      // not symmetric.

      if ((dof[i] >= dof[j]) || (glbmat->Symmetric( ) == 0))
      {
        glbmat->Add(dof[i]-1, dof[j]-1, sprmat[i][j]);
      }
    }
  }

  // Release the dof array.

  delete []dof;
}

// -------------------------------------------------------------------------
// Class cSpringSupp:
// -------------------------------------------------------------------------

// ============================== ReadSpringSupp ===========================

void cSpringSupp :: ReadSpringSupp(void)
{
  // Read the number of spring supports.

  if (!(in >> NumSpring) || (NumSpring < 1))
  {
    cout << "Error in the input of the number of spring supports!\n";
    exit(0);
  }

  // Alloc the array of spring supports.

  VecSpring = new cSpringSupp*[NumSpring];

  // Read each spring.

  int id,nodeid,dir,propid;
  for (int i = 0; i < NumSpring; i++)
  {
    if (!(in >> id) || !(in >> nodeid) || !(in >> dir) || !(in >> propid))
    {
      cout << "Error in the input of the spring support " << i+1 << "!\n";
      exit(0);
    }

    cNode *node = cNode::FindNode(nodeid);
    if (!node)
    {
      cout << "Error in the input of the spring support " << i+1 << "(invalid node)!\n";
      exit(0);
    }

    cSpringProp *prop = cSpringProp::GetSpringProp(propid);
    if (!prop)
    {
      cout << "Error in the input of the spring support " << i+1 << " (invalid property)!\n";
      exit(0);
    }

    VecSpring[i] = new cSpringSupp(id, node, dir, prop);
  }
}

// ================================ Destroy ================================

void cSpringSupp :: Destroy(void)
{
  // Destroy each node.

  for (int i = 0; i < NumSpring; i++) delete VecSpring[i];

  // Release the array of nodes.

  delete []VecSpring;
}

// ============================== cSpringSupp ==============================

cSpringSupp :: cSpringSupp(int label, cNode *node, int dir,cSpringProp *prop) :
               cSpring( )
{
  Label = label;
  Node  = node;
  Dir   = dir;
  Prop  = prop;
  Node->SetSupp( );  // Turn on the support flag
}

// ============================= ~cSpringSupp ==============================

cSpringSupp :: ~cSpringSupp(void)
{
}

// ================================ GetDof =================================

void cSpringSupp :: GetDofs(int *dof)
{
  dof[0] = Node->GetDof(Dir-1);
}

// =============================== IntForce ================================

void cSpringSupp :: IntForce(cVector &g)
{
  double u = Node->GetDispl(Dir-1);

  g[0] = Prop->GetForce(u);
}

// =============================== StiffMat ================================

void cSpringSupp :: StiffMat(cMatrix &K)
{
  double u = Node->GetDispl(Dir-1);

  K[0][0] = Prop->GetStiff(u);
}

// -------------------------------------------------------------------------
// Class cSpringConn:
// -------------------------------------------------------------------------

// ============================== ReadSpringConn ===========================

void cSpringConn :: ReadSpringConn(void)
{
  // Read the number of spring connections.
  if (!(in >> NumSpring) || (NumSpring < 1))
  {
    cout << "Error in the input of the number of spring supports!\n";
    exit(0);
  }

  // Alloc the array of spring connections.

  VecSpring = new cSpringConn*[NumSpring];

  // Read each spring.

  int id,nodeid1,nodeid2,dir,propid;
  for (int i = 0; i < NumSpring; i++)
  {
    if (!(in >> id) || !(in >> nodeid1) || !(in >> nodeid2) || !(in >> dir) ||
        !(in >> propid))
    {
      cout << "Error in the input of the spring support " << i+1 << "!\n";
      exit(0);
    }

    cNode *node1 = cNode::FindNode(nodeid1);
    cNode *node2 = cNode::FindNode(nodeid2);

    if (!node1 || !node2)
    {
     cout << "Error in the input of the spring connection " << i+1 << " (invalid node)!\n";
      exit(0);
    }
    cSpringProp *prop = cSpringProp::GetSpringProp(propid);
    if (!prop)
    {
      cout << "Error in the input of the spring connection " << i+1 << " (invalid property)!\n";
      exit(0);
    }
    // Create and store the spring connection object.

    VecSpring[i] = new cSpringConn(id, node1, node2, dir, prop);
  }
}

// ================================ Destroy ================================

void cSpringConn :: Destroy(void)
{
  // Destroy each node.

  for (int i = 0; i < NumSpring; i++) delete VecSpring[i];

  // Release the array of nodes.

  delete []VecSpring;
}

// ============================== cSpringConn ==============================

cSpringConn :: cSpringConn(int label, cNode *node1, cNode *node2, int dir,cSpringProp *prop) :
               cSpring( )
{
  Label = label;
  Node1 = node1;
  Node2 = node2;
  Dir   = dir;
  Prop  = prop;
}

// ============================= ~cSpringConn ==============================

cSpringConn :: ~cSpringConn(void)
{
}

// ================================ GetDof =================================

void cSpringConn :: GetDofs(int *dof)
{
  dof[0] = Node1->GetDof(Dir-1);
  dof[1] = Node2->GetDof(Dir-1);
}

// =============================== IntForce ================================

void cSpringConn :: IntForce(cVector &g)
{
  double u1 = Node1->GetDispl(Dir-1);
  double u2 = Node2->GetDispl(Dir-1);
  double d  = u2 - u1;                  // Relative displacement
  double f  = Prop->GetForce(d);

  g[0] = -f;
  g[1] =  f;
}

// =============================== StiffMat ================================

void cSpringConn :: StiffMat(cMatrix &K)
{
  double u1 = Node1->GetDispl(Dir-1);
  double u2 = Node2->GetDispl(Dir-1);
  double d  = u2 - u1;                  // Relative displacement
  double k  = Prop->GetStiff(d);

  K[0][0] =  k;
  K[1][1] =  k;
  K[0][1] = -k;
  K[1][0] = -k;
}

// ======================================================= End of file =====
