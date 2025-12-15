// ------------------------------------------------------------------------
// gblvar.h - global variables.
// ------------------------------------------------------------------------

#ifndef _GBLVAR_H
#define _GBLVAR_H

#include <cstdio>
#include <fstream>
#include <string>

extern std::string    fname;     // Input file name
extern std::ifstream  in;        // Input file
extern std::ofstream  out;       // Output file
extern std::ofstream  outem;     // Output file for IGA equivalent mesh
extern bool           Feedback;  // Flag for printed feedback
extern bool           Printing;  // Printing status

#endif
