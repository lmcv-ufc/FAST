// -------------------------------------------------------------------------
// kdtree.tpp - implementation of kdtree data structure.
// -------------------------------------------------------------------------
// Created:      30-Jan-2018     Elias Saraiva Barroso
//
// Modified:
//
// -------------------------------------------------------------------------

#include <iostream>
#include <cmath>
#include <algorithm>

// -------------------------------------------------------------------------
// cKDTree class:
// -------------------------------------------------------------------------

// ============================= cKDTree ===================================

template <class N>
cKDTree<N> :: cKDTree( )
{
  datavec = 0;
  nodes   = 0;

  Tol = 1e-6;
}

// ============================= ~cKDTree ==================================

template <class N>
cKDTree<N> :: ~cKDTree( )
{
  if (nodes)
    delete [] nodes;

  if (datavec)
    delete [] datavec;
}

// ============================= IsReady ===================================

template <class N>
bool cKDTree<N> :: IsReady(void)
{
  // Check tree data.
  if (!nodes || size < 1)
    return false;

  // Check leaves data.
  if (!datavec || numleaves < 1)

  // Check if each node is initialized.
  for(int i = 0; i < size; ++i)
    if (nodes[i].dim < 0)
      return false;

  return true;
}

// ============================= Load ======================================

template <class N>
void cKDTree<N> :: Load(const int &d)
{
  // Release old data.
  Release( );

  // Setup depth, tree size and number of leaves.
  depth     = d;
  size      = round(pow(2,d+1)) - 1;
  numleaves = round(pow(2,d+1));

  // Alloc memory for tree and leafdata.
  nodes    = new sKDTreeNode[size];
  datavec  = new N [numleaves];
}

// ============================= Release ===================================

template <class N>
void cKDTree<N> :: Release( )
{
  if (nodes)
  {
    delete [] nodes;
    nodes = 0;
  }

  if (datavec)
  {
    delete [] datavec;
    datavec = 0;
  }

  size = numleaves = 0;
}

// ============================= GetNodeHeight =============================

template <class N>
int cKDTree<N> :: GetNodeHeight(int i) const
{
  int d = 0;

  while (i > 0)
  {
    i = GetParent(i);
    ++d;
  }

  return d;
}

// ============================= SetNode ===================================

template <class N>
void cKDTree<N> :: SetNode(const int &i, const int &d, const double &v)
{
  nodes[i].dim = d;
  nodes[i].val = v;
}

// ============================= SetTolerance ==============================

template <class N>
void cKDTree<N> :: SetTolerance(const double &tol)
{
  Tol = tol;
}

// ============================= GetNode ===================================

template <class N>
sKDTreeNode* cKDTree<N> :: GetNode(const int &i) const
{
  return &(nodes[i]);
}

template <class N>
int cKDTree<N> :: GetNode(double pos[]) const
{
  int  i = 0;
  bool maxdepth = false;

  while (1)
  {
    if (GetLeft(i) >= size)
      maxdepth = true;

    if ((pos[nodes[i].dim] + Tol) < nodes[i].val)
      i = GetLeft(i);
    else
      i = GetRight(i);

    if (maxdepth) return i;
  }
}

// ============================= GetNodeData ===============================

template <class N>
N& cKDTree<N> :: GetNodeData(double pos[]) const
{
  int i = GetNode(pos);
  return datavec[(i-size)];
}

template <class N>
N& cKDTree<N> :: GetNodeData(const int &il) const
{
  return datavec[il];
}

// ======================================================= End of file =====
