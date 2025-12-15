// -------------------------------------------------------------------------
// subpat.cpp - implementation of sub-patch routines.
// -------------------------------------------------------------------------
// Copyright (c) 2018 LMCV/UFC
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
// Created:      04-Jun-2018    Elias Saraiva Barroso
//
// -------------------------------------------------------------------------

#include "subpat.h"
#include "patch.h"
#include "bspbasis.h"

using namespace std;
using namespace nGeomAlg;

// -------------------------------------------------------------------------
// nCurvAlg namespace:
//

// ============================= CreateSubPatch ============================

cBspCurve* nCurvAlg :: CreateSubPatch(cBspPat *pat, const int &f, const int &e)
{
  // This tensor of indexes maps the r,s,t parameters according to given face 
  // and edge of the patch. The parameter values 0 and 1 means the
  // beginning and ending of the patch in this direction. The value 2 means
  // this parametric value is used in the definition of the isopatch curve. For
  // example, bottom edge of front face is {2,0,0}.

                     // Edge: Top    Bot      Lef     Rig       // Faces:
  static int inds[][4][3] = {{{2,1,0},{2,0,0},{0,2,0},{1,2,0}},  //  Front
                             {{2,1,1},{2,0,1},{1,2,1},{0,2,1}},  //  Back
                             {{0,1,2},{0,0,2},{0,2,1},{0,2,0}},  //  Left
                             {{1,1,2},{1,0,2},{1,2,0},{1,2,1}},  //  Right
                             {{2,1,1},{2,1,0},{0,1,2},{1,1,2}},  //  Above
                             {{2,0,0},{2,0,1},{0,0,2},{1,0,2}}}; //  Below


  // Auxiliary data. 
  cBSplineBasis *bas = 0;  // B-Spline basis of the curve.
  int            pd = 0;   // Parametric variable index used in the new curve.
  int            ncp = 1;  // Number of control points in the the new curve.
  int cnr[3][2] = {{0,0},  // Incidence of this patch control points in the
                   {0,0},  // new curve.        
      	           {0,0}}; 
                            
  // Evaluate the new curve incidence in this patch.
  int nbas;
  for(int pvar = 0; pvar < pat->getNumParVar( ); ++pvar)
  {
    nbas = pat->getBasis(pvar).getNumBas( );
  
    if (inds[f][e][pvar] == 2) // Variables used in new patch.
    {
      pd = pvar;
      cnr[pvar][1] = nbas;
      bas = new cBSplineBasis(pat->getBasis(pvar));
      ncp *= nbas;
    }
    else  // Fixed variables.
    {
      if (inds[f][e][pvar] == 0)
        cnr[pvar][0] = cnr[pvar][1] = 0;
      else
        cnr[pvar][0] = cnr[pvar][1] = nbas - 1;
    }
  }

  // Create B-Spline curve object.
  cBspCurve *curv = new cBspCurve(bas);
  curv->setRational(pat->isRational( ));
  curv->setRelCtrlPntFlag(false);
  
  // Loop over this patch control net and store curve control polygon.
  sCtrlPnt** CtrlPol = new sCtrlPnt*[curv->getNumBas( )];

  int counter = 0;
  int id;
  int pid[3] = {cnr[0][0],cnr[1][0],cnr[2][0]};

  for(pid[pd] = cnr[pd][0]; pid[pd] < cnr[pd][1]; ++pid[pd])
  {
    id = pat->getCPInd(pid[0],pid[1],pid[2]);
    CtrlPol[counter++] = (pat->getCtrlPnt(id));
  }

  curv->setCtrlPnt(CtrlPol,curv->getNumBas( ));
  
  // Release memory
  delete [] CtrlPol;

  return curv;
}

// -------------------------------------------------------------------------
// nSurfAlg namespace:
//

// ============================= CreateSubPatch ============================

cBspSurface* nSurfAlg :: CreateSubPatch(cBspPat *pat, const int &f)
{
  // This tensor of indexes maps the r,s,t parameters according to given face 
  // and edge of the patch. The parameter values 0 and 1 means the
  // beginning and ending of the patch in this direction. The value 2 means
  // this parametric value is used in the definition of the isopatch surface.
  // For example, front face is {2,2,0}.

                                    // Faces:
  static int inds[][3] = {{2,2,0},  // Front
                          {2,2,1},  // Back
                          {0,2,2},  // Left
                          {1,2,2},  // Right
                          {2,1,2},  // Above
                          {2,0,2}}; // Below

  // Auxiliary data. 
  cBSplineBasis *bas[2];   // B-Spline basis of the surface.
  int             pd[2];   // Parametric index used in the new surface.
  int            ncp = 1;  // Number of control points in the the new surface.
  int cnr[3][2] = {{0,0},  // Incidence of this patch control points in the
                   {0,0},  // new surface.        
		   {0,0}}; 
  
  // Evaluate the new surface incidence in this patch.
  int nbas[3];
  int basid = 0;
  for(int pvar = 0; pvar < pat->getNumParVar( ); ++pvar)
  {
    nbas[pvar] = pat->getBasis(pvar).getNumBas( );
  
    if (inds[f][pvar] == 2) // Variables used in new patch.
    {
      pd[basid] = pvar;
      cnr[pvar][1] = nbas[pvar];
      bas[basid++] = new cBSplineBasis(pat->getBasis(pvar));
      ncp *= nbas[pvar];
    }
    else  // Fixed variables.
    {
      if (inds[f][pvar] == 0)
        cnr[pvar][0] = cnr[pvar][1] = 0;
      else
        cnr[pvar][0] = cnr[pvar][1] = nbas[pvar] - 1;
    }
  }

  // Create B-Spline surface object.
  cBspSurface *surf = new cBspSurface(bas[0],bas[1]);
  surf->setRational(pat->isRational( ));
  surf->setRelCtrlPntFlag(false);
  
  // Loop over this patch control net and store surface control net.
  sCtrlPnt** CtrlPol = new sCtrlPnt* [surf->getNumBas( )];

  int counter = 0;
  int id;
  int sr = pd[0], ss = pd[1];
  int pid[3] = {cnr[0][0],cnr[1][0],cnr[2][0]};

  for(pid[ss] = cnr[ss][0]; pid[ss] < cnr[ss][1]; ++pid[ss])
    for(pid[sr] = cnr[sr][0]; pid[sr] < cnr[sr][1]; ++pid[sr])
    {
      id = pat->getCPInd(pid[0],pid[1],pid[2]);
      CtrlPol[counter++] = pat->getCtrlPnt(id);
    }

  surf->setCtrlPnt(CtrlPol,surf->getNumBas( ));

  // Release memory
  delete [] CtrlPol;

  return surf;
}

// ======================================================= End of file =====
