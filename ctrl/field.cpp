#include<iostream>
#include <cmath>

using namespace std;

#include "field.h"
#include "gbldef.h"
#include "gblvar.h"
#include "utl.h"

cField* cField :: f =  new cField;
//cField* cField :: f =  new cPlateHoleHeat;


cField :: cField( )
{
  displ_flag  = false; 
  strs_flag   = false; 
  strn_flag   = false; 
  temp_flag   = false; 
}

// ========================== CompileDisplField ============================

int cField :: CompileDisplField(int exprid)
{
  // Setup variables table.
  expr_displ.AddVar("x",z);
  expr_displ.AddVar("y",y);
  expr_displ.AddVar("z",z);
  expr_displ.AddVar("u",displ[0]);
  expr_displ.AddVar("v",displ[1]);
  expr_displ.AddVar("w",displ[2]);
  expr_displ.AddVar("rx",displ[3]);
  expr_displ.AddVar("ry",displ[4]);
  expr_displ.AddVar("rz",displ[5]);

  // Compile input expression.
  expr_displ.register_symbol_table( );
  return expr_displ.compile(Utl :: GetExp(exprid));
}

// ========================== CompileStrainField ===========================

int cField :: CompileStrainField(int exprid)
{
  // Setup variables table.
  expr_strn.AddVar("x",x);
  expr_strn.AddVar("y",y);
  expr_strn.AddVar("z",z);
  expr_strn.AddVar("sxx",strn[0]);
  expr_strn.AddVar("syy",strn[1]);
  expr_strn.AddVar("szz",strn[2]);
  expr_strn.AddVar("sxy",strn[3]);
  expr_strn.AddVar("sxz",strn[4]);
  expr_strn.AddVar("syz",strn[5]);

  // Compile input expression.
  expr_strn.register_symbol_table( );
  return expr_strn.compile(Utl :: GetExp(exprid));
}

// ========================== CompileStressField ===========================

int cField :: CompileStressField(int exprid)
{
  // Setup variables table.
  expr_strs.AddVar("x",x);
  expr_strs.AddVar("y",y);
  expr_strs.AddVar("z",z);
  expr_strs.AddVar("sxx",strs[0]);
  expr_strs.AddVar("syy",strs[1]);
  expr_strs.AddVar("szz",strs[2]);
  expr_strs.AddVar("sxy",strs[3]);
  expr_strs.AddVar("sxz",strs[4]);
  expr_strs.AddVar("syz",strs[5]);

  // Compile input expression.
  expr_strs.register_symbol_table( );
  return expr_strs.compile(Utl :: GetExp(exprid));
}

// ========================== CompileTempField =============================

int cField :: CompileTempField(int exprid)
{
  // Setup variables table.
  expr_temp.AddVar("x",x);
  expr_temp.AddVar("y",y);
  expr_temp.AddVar("z",z);
  expr_temp.AddVar("temp",temp);

  // Compile input expression.
  expr_temp.register_symbol_table( );
  return expr_temp.compile(Utl :: GetExp(exprid));
}

// ============================ ReadDisplField =============================

void cField:: ReadDisplField( )
{
  // Read input expression label.
  int exprid;
  if (!(in >> exprid) || (exprid <= 0) || (exprid > Utl :: GetNumExp( )))
  {
    cout << "Invalid expression label: " << exprid << "!\n";
    exit(0);
  }

  // Compile displacement field.
  if (!GetInstance( )->CompileDisplField(exprid))
  {
    cout << "Error in compilation of the expression " << exprid;
    cout << " for displacement field!" << endl;
    exit(0);
  }

  GetInstance( )->displ_flag = true;
}

// ============================ ReadStrainField ============================

void cField:: ReadStrainField( )
{
  // Read input expression label.
  int exprid;
  if (!(in >> exprid) || (exprid <= 0) || (exprid > Utl :: GetNumExp( )))
  {
    cout << "Invalid expression label: " << exprid << "!\n";
    exit(0);
  }

  // Compile displacement field.
  if (!GetInstance( )->CompileStrainField(exprid))
  {
    cout << "Error in compilation of the expression " << exprid;
    cout << " for strain field!" << endl;
    exit(0);
  }

 GetInstance( )->strn_flag = true;
}

// ============================ ReadStressField ============================

void cField:: ReadStressField( )
{
  // Read input expression label.
  int exprid;
  if (!(in >> exprid) || (exprid <= 0) || (exprid > Utl :: GetNumExp( )))
  {
    cout << "Invalid expression label: " << exprid << "!\n";
    exit(0);
  }

  // Compile displacement field.
  if (!GetInstance( )->CompileStressField(exprid))
  {
    cout << "Error in compilation of the expression " << exprid;
    cout << " for stress field!" << endl;
    cout << Utl :: GetExp(exprid) << endl;
    exit(0);
  }

  GetInstance( )->strs_flag = true;
}

// ============================ ReadTempField ==============================

void cField:: ReadTempField( )
{
  // Read input expression label.
  int exprid;
  if (!(in >> exprid) || (exprid <= 0) || (exprid > Utl :: GetNumExp( )))
  {
    cout << "Invalid expression label: " << exprid << "!\n";
    exit(0);
  }

  // Compile displacement field.
  if (!GetInstance( )->CompileTempField(exprid))
  {
    cout << "Error in compilation of the expression " << exprid;
    cout << " for temperature field!" << endl;
    exit(0);
  }

  GetInstance( )->temp_flag = true;
}

void cField :: ReadField(void)
{
  // Create all possible solution fields.
  static cInfSheetHole  f1;
  static cPlateHoleHeat f2;

  // Read the field label.
  char label[100];
  if (!Utl::ReadString(in, label))
  {
    cout << "Error in the input of the field label!\n";
    exit(0);
  }

  // Create the appropriate algorithm.
  if (string(label) == "InfSheetHole")
    f = &f1;
  else if (string(label) == "PlateHoleHeat")
    f = &f2;
  else
  {
    cout << "Unknown field: " << label << "\n";
    exit(0);
  }
}


void cField :: Displ(double x, double y, double z, double dvec[])
{
  // Store variables.
  this->x = x;  
  this->y = y;  
  this->z = z;  

  // Evaluate expression.
  expr_displ.evaluate( );

  // Store evaluated displacements.
  for(int i = 0; i < 6; ++i) 
    dvec[i] = displ[i];
}

void cField :: Stress(double x, double y, double z, double svec[])
{
  // Store variables.
  this->x = x;  
  this->y = y;  
  this->z = z;

  // Evaluate expression.
  expr_strs.evaluate( );

  // Store evaluated displacements.
  for(int i = 0; i < 6; ++i) 
    svec[i] = strs[i];
}

void cField :: Strain(double x, double y, double z, double svec[])
{
  // Store variables.
  this->x = x;  
  this->y = y;  
  this->z = z;  

  // Evaluate expression.
  expr_strn.evaluate( );

  // Store evaluated displacements.
  for(int i = 0; i < 6; ++i) 
    svec[i] = strn[i];
}

double cField :: Temp(double x, double y, double z)
{
  // Store variables.
  this->x = x;  
  this->y = y;  
  this->z = z;  

  // Evaluate expression.
  expr_temp.evaluate( );

  return temp;
}

cInfSheetHole :: cInfSheetHole( )
{
  R = 1.0;
  T = 10;
}

void cInfSheetHole :: Sig(double x, double y, double z, double sv[])
{
  // Evaluate polar coordinates.
  double PI    = 3.14159265359;
  //double theta = (x < -1e-6) ? atan(y/x) : 1.5 * PI; 
  double theta = atan2(y,x); // = (x < -1e-6) ? atan(y/x) : 1.5 * PI; 
  double r     = sqrt(x*x+y*y);

  if (fabs(theta) > 5)
    cout << "theta  " << theta << endl;

  // Evaluate stresses in polar coordinates.
  double sp[3];
  double  c2t = cos(theta*2);
  double  s2t = sin(theta*2);
  double  a   = (R*R)/(r*r);
       sp[0]  = T/2.0 * (1.0 - a) + T/2.0 * (1.0 - 4*a + 3*a*a) * c2t;
       sp[1]  = T/2.0 * (1.0 + a) - T/2.0 * (1.0 + 3*a*a) * c2t;
       sp[2]  = -T/2.0 * (1.0 + 2*a - 3*a*a) * s2t;

  // Evaluate stresses in cartesian system.
  double c  = cos(theta);
  double c2 = c*c;
  double s  = sin(theta);
  double s2 = s*s;

  for(int i = 0; i < 9; ++i) sv[i] = 0.0;

  sv[0] = sp[0] * c2 + sp[1] * s2 - 2 * sp[2] * s * c;
  sv[4] = sp[0] * s2 + sp[1] * c2 + 2 * sp[2] * s * c;
  sv[1] = (sp[0] - sp[1]) * s * c + sp[2] * (c2 - s2);
  sv[3] = sv[1];

/*
  cout << "x " << x  << ", y " << y << endl;
  cout << "r " << r  << ", theta " << theta << endl;
  for(int i = 0; i < 3; ++i)
    for(int j = 0; j < 3; ++j)
      cout << sv[i*3+j] << (j==2) ? "\n" : " ";
  cout << endl;
  */
}

cPlateHoleHeat :: cPlateHoleHeat( )
{
  a = 0.3;
}

double cPlateHoleHeat :: Temp(double x, double y, double z)
{
  /* Problem data              */
  double ro = 2;     // Outer radius.
  double ri = 1;     // Inner radius.
  double uo = 10;    // Outer temperature.
  double ui = 20;    // Inner temperature.

  /* Evaluate polar coordinates             */
  double r     = sqrt(x*x+y*y);

  /* Evaluate temperature */
  return (uo*log(r/ri)-ui*log(r/ro)) / log(ro/ri); 

  //return x*x + y*y - 2.0*a*sqrt(x*x+y*y) + a*a;
}

double cPlateHoleHeat :: Source(double x, double y, double z)
{
  return 2.0 * a / sqrt(x*x+y*y) - 4.0;
}

void cPlateHoleHeat :: Flux(double x, double y, double z, double f[])
{
  double sq = sqrt(x*x+y*y);

  f[0] = 2.0 * x - 2.0 * a * x / sq;
  f[1] = 2.0 * y - 2.0 * a * y / sq;
  f[2] = 0.0;
}


