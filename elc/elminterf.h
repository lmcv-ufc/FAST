// ------------------------------------------------------------------------
// elminterf.h - class to handle parametric interface elements.
// ------------------------------------------------------------------------

#ifndef _ELMINTERF_H
#define _ELMINTERF_H

#include "elmparam.h"

// ------------------------------------------------------------------------
// Definition of the Parametric Interface Element class:
//
class cElmInterf : public cElmParam
{
 public:
                cElmInterf(int, eShpType, eAnmType);
  virtual      ~cElmInterf(void);
  virtual int   IntForce(cVector &);
  virtual void  StiffMat(cMatrix &);
  virtual void  MassMat(cMatrix &);
  virtual void  GeomStiff(cMatrix &);
  virtual void  NodalStress(cMatrix &);
  virtual void  IntPntStress(cMatrix &, cMatrix &);
};

#endif

