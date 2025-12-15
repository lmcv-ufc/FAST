// -------------------------------------------------------------------------
// loadiga.h - file containing the definition of isogeometric loads.
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

#ifndef _LOADIGA_H
#define _LOADIGA_H

#include "load.h"

// -------------------------------------------------------------------------
// Definition of the patch concentrate force class:
//
class cPatConcLoad : public cParConcLoad 
{
 public:
                 cPatConcLoad(void);
  virtual       ~cPatConcLoad(void);
          void   Read(void);
};

// -------------------------------------------------------------------------
// Definition of the edge uniform force class for IGA problems:
//
class cEdgeUnifLoadIGA : public cEdgeUnifLoad
{
 protected:
  int facetype;
  int edgetype;

 public:
                 cEdgeUnifLoadIGA(void);
  virtual       ~cEdgeUnifLoadIGA(void);
  void           Read(void);
};

// -------------------------------------------------------------------------
// Definition of the edge uniform moment class for IGA problems:
//
class cEdgeUnifMomentIGA : public cEdgeUnifMoment
{
 protected:
  int facetype;
  int edgetype;

 public:
                 cEdgeUnifMomentIGA(void);
  virtual       ~cEdgeUnifMomentIGA(void);
  void           Read(void);


};

// -------------------------------------------------------------------------
// Definition of the face uniform force class:
//
class cFaceUnifLoadIGA : public cFaceUnifLoad
{
 protected:
  int facetype;

 public:
                 cFaceUnifLoadIGA(void);
  virtual       ~cFaceUnifLoadIGA(void);
  virtual void   Read(void);
};

// -------------------------------------------------------------------------
// Definition of the edge general force class:
//
class cEdgeGenLoadIGA : public cEdgeGenLoad
{
 protected:
  int facetype;
  int edgetype;

 public:
                 cEdgeGenLoadIGA(void);
  virtual       ~cEdgeGenLoadIGA(void);
  virtual void   Read(void);
};

#endif
