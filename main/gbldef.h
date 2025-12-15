// ------------------------------------------------------------------------
// gbldef.h - global definitions.
// ------------------------------------------------------------------------

#ifndef _GBLDEF_H
#define _GBLDEF_H

// ------------------------------------------------------------------------
// Cartesian axes:
//
typedef enum
{
  X_AXIS = 1,
  Y_AXIS,
  Z_AXIS
} eAxis;

// ------------------------------------------------------------------------
// Response labels:
//
enum
{
  FORCE_X = 0,
  FORCE_Y,
  FORCE_Z,
  FORCE_XY,
  FORCE_XZ,
  FORCE_YZ,
  MOMENT_X,
  MOMENT_Y,
  MOMENT_Z,
  MOMENT_XY,
  STRESS_XX,
  STRESS_YY,
  STRESS_ZZ,
  STRESS_XY,
  STRESS_XZ,
  STRESS_YZ,
  STRAIN_XX,
  STRAIN_YY,
  STRAIN_ZZ,
  STRAIN_XY,
  STRAIN_XZ,
  STRAIN_YZ,
  P_XX,
  P_YY,
  P_XY,
  R_XZ,
  R_YZ,
  HEAT_FLUX_X,
  HEAT_FLUX_Y,
  HEAT_FLUX_Z,
  MAX_SCL_LAB     // Total number of scalar labels
};

typedef char tLabel[14];

const tLabel VEC_SCL_LAB[] =
{
  "FORCE_X",
  "FORCE_Y",
  "FORCE_Z",
  "FORCE_XY",
  "FORCE_XZ",
  "FORCE_YZ",
  "MOMENT_X",
  "MOMENT_Y",
  "MOMENT_Z",
  "MOMENT_XY",
  "STRESS_XX",
  "STRESS_YY",
  "STRESS_ZZ",
  "STRESS_XY",
  "STRESS_XZ",
  "STRESS_YZ",
  "STRAIN_XX",
  "STRAIN_YY",
  "STRAIN_ZZ",
  "STRAIN_XY",
  "STRAIN_XZ",
  "STRAIN_YZ",
  "P_XX",
  "P_YY",
  "P_XY",
  "R_XZ",
  "R_YZ",
  "HEAT_FLUX_X",
  "HEAT_FLUX_Y",
  "HEAT_FLUX_Z"
};

enum
{
  STRESSES = 0,
  STRAINS,
  PRINCSTRESSES
};


// ------------------------------------------------------------------------
// Mathematical Expression Toolkit Library 
//
#include <string>

#ifdef _EXPR_
#include "exprtk.hpp"
#endif

typedef struct
{
 #ifdef _EXPR_
  exprtk::symbol_table<double> table;
  exprtk::expression<double>   exp;
  exprtk::parser<double>       parser;

  void register_symbol_table(void)
  {
    exp.register_symbol_table(table);
  }

  int compile(std::string str)
  {
    return parser.compile(str,exp);
  }

  double evaluate(void)
  {
    return exp.value( );
  }

  void AddVar(const std::string tvar, double &var)
  {
    table.add_variable(tvar,var);
  }
 #else
  void   register_symbol_table(void)       {             }
  int    compile(std::string str)          { return 0;   }
  double evaluate(void)                    { return 0.0; }
  void   AddVar(const std::string,double&) {             }
 #endif

} sExpEval;

// ------------------------------------------------------------------------
// Useful constants:
//
const double PI = 3.141592653589793;
const double PENALTY_STIFFNESS = 1.0e30;

// ------------------------------------------------------------------------
// Useful functions:
//
template <typename T>
inline T MAX(T x, T y)
{
  if (x < y)
    return y;
  else
    return x;
}

template <typename T>
inline T MIN(T x, T y)
{
  if (x > y)
    return y;
  else
    return x;
}

template <typename T>
inline T Sign(T x)
{
  if (x < 0)
    return -1;
  else
    return 1;
}

#endif
