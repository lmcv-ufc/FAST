// -------------------------------------------------------------------------
// timefunc.h - file containing the definition of the TimeFunction class.
// -------------------------------------------------------------------------
//
// The TimeFunction class handle the time dependent characteristic of the
// applied loads.
//
// cTimeFunc
// |-- cFuncConst => constant time function.
// |-- cFuncHarmonic => harmonic (sine) time function.
// |-- cFuncHalfSine => repeated half-sine time function.
// |-- cFuncTable => piecewise linear time function.
//
// -------------------------------------------------------------------------
// Static public methods:
// -------------------------------------------------------------------------
//
// static void ReadTimeFuncConst(void)
//
// This method creates and reads a constant time function.
// -------------------------------------------------------------------------
//
// static void ReadTimeFuncHarmonic(void)
//
// This method creates and reads a harmonic time function.
// -------------------------------------------------------------------------
//
// static void ReadTimeFuncHalfSine(void)
//
// This method creates and reads a time function with the shape of a
// repeated half-sine function.
// -------------------------------------------------------------------------
//
// static void ReadTimeFuncTable(void)
//
// This method creates and reads a piecewise linear (tabular) time function.
// -------------------------------------------------------------------------
//
// static void Destroy(void)
//
// This method destroys all stored time functions.
// -------------------------------------------------------------------------
//
// static cTimeFunc *GetCurr(void)
//
// This method returns the current time function, which is generally the
// last one to be created in the time when this method is called.
//
// -------------------------------------------------------------------------
// Pure virtual methods:
// -------------------------------------------------------------------------
//
// virtual void Read(void)
//
// This method reads the function parameters.
// -------------------------------------------------------------------------
//
// virtual double GetVal(double t)
//
//  t      -  time                                                    (in)
//
// This method returns the function value at the given time.
//
// -------------------------------------------------------------------------

#ifndef _TIMEFUNC_H
#define _TIMEFUNC_H

// -------------------------------------------------------------------------
// Definition of the Time Function class:
//
class cTimeFunc
{
 private:
  static  cTimeFunc  *Head;
  static  cTimeFunc  *Tail;
  static  cTimeFunc  *Curr;

 protected:
          cTimeFunc  *Next;

 public:
  static  void       ReadTimeFuncConst(void);
  static  void       ReadTimeFuncHarmonic(void);
  static  void       ReadTimeFuncTable(void);
  static  void       ReadTimeFuncHalfSine(void);
  static  void       Destroy(void);
  static  cTimeFunc *GetCurr(void) { return Curr; }
                     cTimeFunc(void);
  virtual           ~cTimeFunc(void);
  virtual void       Read(void) = 0;
  virtual double     GetVal(double) = 0;
};


// -------------------------------------------------------------------------
// Definition of the Constant Time Function class:
//
class cFuncConst : public cTimeFunc
{
 protected:
  double Val;

 public:
                  cFuncConst(void);
                  cFuncConst(double);
  virtual        ~cFuncConst(void);
          void    Read(void);
          double  GetVal(double) { return Val; }
};

// -------------------------------------------------------------------------
// Definition of the Harmonic (sine) Time Function class:
//
class cFuncHarmonic : public cTimeFunc
{
 protected:
  double Amplit;  // Amplitude
  double Period;  // Period of the sine funcion
  double Phase;   // Phase angle (radians)
  double T0,T1;   // Time interval where the function is non-zero

 public:
                  cFuncHarmonic(void);
                  cFuncHarmonic(double, double, double, double, double);
  virtual        ~cFuncHarmonic(void);
          void    Read(void);
          double  GetVal(double);
};

// -------------------------------------------------------------------------
// Definition of the Harmonic (half-sine) Time Function class:
//
class cFuncHalfSine : public cTimeFunc
{
 protected:
  double Amplit;  // Amplitude
  double Period;  // Period of the sine funcion
  double Rest;    // Rest time

 public:
                  cFuncHalfSine(void);
                  cFuncHalfSine(double, double, double);
  virtual        ~cFuncHalfSine(void);
          void    Read(void);
          double  GetVal(double);
};

// -------------------------------------------------------------------------
// Definition of the Table Time Function class:
//
typedef struct
{
  double t;
  double f;
} sTablePnt;

class cFuncTable : public cTimeFunc
{
 protected:
  int        NumPnt;
  sTablePnt *Pnt;

 public:
                  cFuncTable(void);
                  cFuncTable(int, double *, double *);
  virtual        ~cFuncTable(void);
          void    Read(void);
          double  GetVal(double);
};

#endif
