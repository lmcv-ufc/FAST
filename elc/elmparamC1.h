// ------------------------------------------------------------------------
// elmparamC1.h - class to handle parametric elements with C1 continuity.
// ------------------------------------------------------------------------

#ifndef _ELMPARAMC1_H
#define _ELMPARAMC1_H

#include <string>

#include "elmparam.h"

// ------------------------------------------------------------------------
// Definition of the Parametric Element C1 class:
//
class cElmParamC1 : public cElmParam
{
 public:
  static  void      ReadElmBSP(std::string, eShpType, eAnmType, bool tl = 0);

                    cElmParamC1(int, eShpType, eAnmType);
  virtual          ~cElmParamC1(void) {}
  
  virtual int       IntForce(cVector &);
  virtual void      StiffMat(cMatrix &);
  virtual void      MassMat(cMatrix &);
  virtual void      GeomStiff(cMatrix &);
  virtual void      IntPntStress(cMatrix &, cMatrix &);
  virtual void      PntStress(sNatCoord*, int, cMatrix&, cMatrix&, bool opt = 0);
};

// ------------------------------------------------------------------------
// Definition of the Parametric Total Lagrangian Element C1 class:
//
class cElmParamC1TL : public cElmParamC1
{
 protected:
          void      GreenStrain(cMatrix &, cMatrix &, cVector &, cVector &);

 public:
   static  void     ReadElmBSP(std::string,eShpType,eAnmType, bool tl = 0);

                    cElmParamC1TL(int, eShpType, eAnmType);
  virtual          ~cElmParamC1TL(void) {}

  virtual int       IntForce(cVector &);
  virtual void      StiffMat(cMatrix &);
  virtual void      IntPntStress(cMatrix &, cMatrix &);
  virtual void      PntStress(sNatCoord*, int, cMatrix&, cMatrix&, bool opt = 0);
};

#endif
