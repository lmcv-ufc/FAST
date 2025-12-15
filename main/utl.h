// -------------------------------------------------------------------------
// utl.h - This module contains utilitary functions.
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

#ifndef _UTL_H
#define _UTL_H

//#include <stdio.h>
//#include <stdlib.h>
#include <vector>
//#include <iostream>
#include <string>

#include <iosfwd>

class Utl
{
 protected:
  static std::vector<std::string> ExpVec;

 public:
  static int         NextLabel   (std::istream &fp, char *s);
  static int         ReadString  (std::istream &fp, char *s);
  static int         CountString (const char*);
  static int         RandInt     (int floor, int ceiling);
  static int         GetNumExp   (void)   { return ExpVec.size( ); } 
  static double      RandDec     (void);
  static double      RandDouble  (double floor, double ceiling);
  static std::string  GetFilePath(const std::string&);
  static double      Min         (double val1, double val2);
  static double      Max         (double val1, double val2);
  static void        Exit        (std::string messeger);
  static void        ReadExp     (void);
  static std::string GetExp      (int id) { return ExpVec[id-1]; } 
};


#endif


