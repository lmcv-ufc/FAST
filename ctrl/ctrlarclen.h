// -------------------------------------------------------------------------
// ctrlarclen.h - file containing the definition of the cCtrlArcLen class.
// -------------------------------------------------------------------------
//
// The cCtrlArcLen class implements the analysis of nonlinear static
// structures by arc-length method with Newton-Raphson iterations.
//
// -------------------------------------------------------------------------

#ifndef _CTRLARCLEN_H
#define _CTRLARCLEN_H

#include "ctrl.h"
#include "ctrlpath.h"

// -------------------------------------------------------------------------
// Definition of cCtrlArcLength base class:
//
class cCtrlArcLength : public cCtrlPath
{
 protected:
  static  double   LoadCoeff;  // Coefficient of load terms (psi)
          double   ArcLen;     // Current arc length
          double   Dlf;        // Step load increment
  static  bool     Restart;    // Flag for restarts

  virtual void     SetArcLength(cVector &);
  virtual double   FirstIncLoadFac(cVector &, cVector &);

 public:
  static  void     ReadData(void);
                   cCtrlArcLength(void);
                  ~cCtrlArcLength(void);
};

// -------------------------------------------------------------------------
// Definition of cIniOrtArcLength class:
//
class cIniOrtArcLength: public cCtrlArcLength
{
 public:
                   cIniOrtArcLength(void);
                  ~cIniOrtArcLength(void);
          double   IncLoadFac(cVector &, cVector &, cVector &);
};

// -------------------------------------------------------------------------
// Definition of cUpOrtArcLength class:
//
class cUpOrtArcLength: public cCtrlArcLength
{
 public:
                   cUpOrtArcLength(void);
                  ~cUpOrtArcLength(void);
          double   IncLoadFac(cVector &, cVector &, cVector &);
};

// -------------------------------------------------------------------------
// Definition of cConLinArcLength class:
//
class cConLinArcLength: public cCtrlArcLength
{
 public:
                   cConLinArcLength(void);
                  ~cConLinArcLength(void);
          double   IncLoadFac(cVector &, cVector &, cVector &);
};

// -------------------------------------------------------------------------
// Definition of cCilArcLength class:
//
class cCilArcLength: public cCtrlArcLength
{
 public:
                   cCilArcLength(void);
                  ~cCilArcLength(void);
          double   IncLoadFac(cVector &, cVector &, cVector &);
};

#endif
