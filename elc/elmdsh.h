// ------------------------------------------------------------------------
// elmdsh.h - class to handle degenerated shell elements.
// ------------------------------------------------------------------------
//
// ------------------------------------------------------------------------

#ifndef _ELMDSH_H
#define _ELMDSH_H

#include "elmparam.h"

// ------------------------------------------------------------------------
// Definition of the Degenerated Shell Element class:
//
class cElmDegShell : public cElmParam
{
 public:
  static void       CompDegShlAxes(void);
  static void       PrintNormalVec(int,cVector*,cVector*,cVector*);

                    cElmDegShell(int, eShpType, eAnmType);
  virtual          ~cElmDegShell(void);
  
  virtual void      NodalStress(cMatrix &);
};

// ------------------------------------------------------------------------
// Definition of the Total Lagrangian Degenerated Shell Element class:
//
class cElmDegShellTL : public cElmParamTL
{
 protected:
          void      GreenStrain(int, double *, cVector &, sNodeCoord *, sNodeCoord *, cVector &);
          void      KgnlMatrix(int, double *, sNodeCoord *, sNodeCoord *, cVector &, cVector &, cMatrix &);
          
 public:
  static  void      ReadUpdateOption(void);	  

                    cElmDegShellTL(int, eShpType, eAnmType);
  virtual          ~cElmDegShellTL(void);
  
  virtual void      NodalStress(cMatrix &);
  
          int       IntForce(cVector &);
          void      StiffMat(cMatrix &);
          void      IntPntStress(cMatrix &, cMatrix &);
          void      UpdateItera(cVector &);
};

#endif
