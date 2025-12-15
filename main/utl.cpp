// -------------------------------------------------------------------------
// utl.cpp - utilitary functions.
// -------------------------------------------------------------------------
// Copyright (c) 2013 LMCV/UFC
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


#include <cstdlib>
#include <math.h>
#include <iostream>
#include <cctype>
#include <sstream>

using namespace std;

#include "utl.h"
#include "gblvar.h"

vector<string> Utl :: ExpVec;

// -------------------------------------------------------------------------
// Public functions:
//

// ============================ ReadExp ====================================

void Utl :: ReadExp(void)
{
  // Read the number expressions.
  int nexp;

  if (!(in >> nexp) || (nexp <= 0))
  {
    cout << "Error in the input of the number of expressions!\n";
    exit(0);
  }
  ExpVec.resize(nexp);

  // Read expressions.

  int id;
  for (int i = 0; i < nexp; i++)
  {
    // Read expression label.
    if (!(in >> id) || (id <=0) || (id > nexp))
    {
      cout << "Error in the input of expression label!\n";
      exit(0);
    }

    // Read input code.
    char aux[500];
    if (!ReadString(in,aux))
    {
      cout << "Error in the input of expression label!\n";
      exit(0);
    }
    
    if (aux[0] == 'd')        // Direct input.
    {
      if (!ReadString(in,aux))
      {
        cout << "Error in the input of expression " << id << " !\n";
        exit(0);
      }

      ExpVec[id-1] = string(aux);
    }
    else if (aux[0] == 'f')   // Input file.
    {
      // Read file name.
      if (!ReadString(in,aux))
      {
        cout << "Error in the input of expression " << id << " file name !\n";
        exit(0);
      }
      string fpath = GetFilePath(fname) + string(aux);
      ifstream file(fpath.c_str( ));
 
      if (!file.is_open( ))
      {
        cout << "Error in the input of expression " << id << "!\n";
	cout << "Unable to open file: " << fpath << endl;
        exit(0);
      }

      stringstream ss;
      ss << file.rdbuf( );
      ExpVec[id-1] = ss.str( );
      //cout << ss.str( ) << endl;
    }
    else
    {
      cout << "Invalid expression input code.\n";
      exit(0);
    }
  }
}

// ============================ UtlNextLabel ============================

int Utl :: NextLabel(istream &fp, char *s)
{
  while (fp.get() != '%')
    if (fp.eof())
      return(0);

  fp >> s;
  return(1);
}

// ==============================  ReadString  ==========================

int Utl :: ReadString(istream &fp, char *s)
{
  int i,c;
 
  while ((c = fp.get()) != '\'')
    if (fp.eof())
      return(0);

  for (i = 0; (c = fp.get()) != '\''; i++)
  {
    if (fp.eof()) return(0);
    s[i] = c;
  }

  s[i] = '\0';
  return(1);
}

// ==============================  CountString  =======================

int Utl :: CountString(const char* str)
{
  if (str == 0)
    return 0; 

  bool inSpaces = true;
  int numWords = 0;

  string line(str);

  while (*str != 0)
  {
    if (std::isspace(*str))
      inSpaces = true;
    else if (inSpaces)
    {		
      numWords++;
      inSpaces = false;
    }
    
    ++str;
  }

  return numWords;
}

// ============================== RandInt =============================

int Utl :: RandInt(int floor, int ceiling)
{
  return floor + rand() % (ceiling-floor+1);
}

// ============================== RandDec =============================

double Utl :: RandDec(void)
{
  return rand() / (RAND_MAX+1.0);
}

// ============================ RandDouble ============================

double Utl :: RandDouble(double floor, double ceiling)
{
  double f = (double)rand() / RAND_MAX;
  return floor + f * (ceiling - floor);
}

// ================================ Min ===============================

double Utl :: Min(double val1, double val2)
{
  if (val1 < val2)
    return val1;
  else
    return val2;
}

// ================================ Max ===============================

double Utl :: Max(double val1, double val2)
{
  if (val1 > val2)
    return val1;
  else
    return val2;
}

// ================================ FilePath ===============================

string Utl :: GetFilePath(const string& str)
{
  size_t found = str.find_last_of("/\\");
  return (found != string::npos) ? str.substr(0,found+1) : string("");
}

// ================================ Exit ==============================

void Utl :: Exit(string messeger)
{
  cout << "\n" << messeger << "\n";
  exit(0);
}

// ================================================= End of file ========
