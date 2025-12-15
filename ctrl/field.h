
#ifndef _FIELD_H
#define _FIELD_H

#include "field.h"
#include <cmath>
#include "gbldef.h"


class cField
{
 private:

 protected:
  cField(void);
  static cField *f;

  // Expression evaluator variables.
  sExpEval expr_displ;
  sExpEval expr_strs;
  sExpEval expr_strn;
  sExpEval expr_temp;
  //symbol_table  tab_d, tab_strs, tab_strn, tab_temp;
  //expression    exp_d, exp_strs, exp_strn, exp_temp;
  //parser        pars;

  // Input variables used in the expression evaluator.
  double x, y, z;     // Input cartesian coordinates.

  // Output variables used in the expression evaluator.
  double displ[6];    // Displacement vector.
  double strs[6];     // Stress vector.
  double strn[6];     // Strain vector.  
  double temp;        // Temperature.
  
  bool displ_flag;
  bool strs_flag;
  bool strn_flag;
  bool temp_flag;

 public: 
  static void    ReadDisplField(void);
  static void    ReadStrainField(void);
  static void    ReadStressField(void);
  static void    ReadTempField(void);
  static void    ReadField(void);

         bool    hasDisplField(void)  { return displ_flag;  }
         bool    hasStressField(void) { return strs_flag;   }
         bool    hasStrainField(void) { return strn_flag;   }
         bool    hasTempField(void)   { return temp_flag;   }

  static cField* GetInstance(void)    { return f;           }


  int CompileDisplField(int);
  int CompileStrainField(int);
  int CompileStressField(int);
  int CompileTempField(int);

  virtual void   Displ(double,double,double,double[]);
  virtual void   Stress(double,double,double,double[]);
  virtual void   Strain(double,double,double,double[]);
  virtual double Temp(double,double,double);


  virtual void   Sig(double,double,double,double[])   { } 
  virtual void   Flux(double,double,double,double[])  { } 
  virtual double Source(double, double, double)       { return 0; }
};

class cInfSheetHole : public cField
{
 protected:
  double R;
  double T;
 
 public:
  cInfSheetHole(void);
  void Sig(double,double,double,double[]);
  
};

class cPlateHoleHeat : public cField
{
 protected:
  double a;
 
 public:
  cPlateHoleHeat(void);
  double Temp(double,double,double);
  void   Flux(double,double,double,double[]);
  double Source(double,double,double);
};


#endif
