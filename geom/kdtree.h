// -------------------------------------------------------------------------
// Copyright (c) 2015 LMCV/UFC
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without 
// modification, are permitted provided that the following conditions
// are met:
//
// 1. Redistributions of source code must retain the above copyright 
//    notice, this list of conditions and the following disclaimer.
// 2. Redistributions in binary form must reproduce the above copyright 
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS 
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT 
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
// PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT 
// HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, 
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED 
// TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR 
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF 
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING 
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS 
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// -------------------------------------------------------------------------
//
// The template class cKDTree implements a static binary decom

#ifndef _KDTREE_H
#define _KDTREE_H

#include <vector>

// -------------------------------------------------------------------------
// Struct sKDTreeNode:
//
struct sKDTreeNode
{
  int    dim;
  double val;

  sKDTreeNode(void) : dim(-1), val(0.0) { }
};

// -------------------------------------------------------------------------
// Template class cKDTree:
//
template <class NodeData>
class cKDTree
{
 protected:
  int depth;              // Binary tree max depth.
  int size;               // Number of nodes.
  int numleaves;          // Number of leaves.
  double Tol;             // Tolerance used in floops.
  NodeData    *datavec;   // Vector of node data (only in leaves).
  sKDTreeNode *nodes;     // Tree nodes stored in a vector.

 public:

               cKDTree(void);
              ~cKDTree(void);

         void  Load(const int&);
         void  Release(void);
         bool  IsReady(void);
         void  SetTolerance(const double&);
         void  SetNode(const int&,const int&,const double&);
  inline int   GetLeft(const int &i) const { return (i * 2 + 1); }
  inline int   GetRight(const int &i) const { return (i * 2 + 2); }
  inline int   GetParent(const int &i) const { return ((i - 1) / 2); }

  sKDTreeNode* GetNode(const int &) const;
  int          GetNode(double pos[]) const;
  NodeData&    GetNodeData(double pos[]) const;
  NodeData&    GetNodeData(const int&) const;
  const int&   GetTreeSize(void) const { return size; }
  const int&   GetLeafSize(void) const { return numleaves; }
        int    GetNodeHeight(int i) const;
};

// -------------------------------------------------------------------------
// Include template implementation file:
//
#include "kdtree.tpp"

// ======================================================= End of file =====

#endif
