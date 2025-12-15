// ------------------------------------------------------------------------
// elmparam.h - class to handle parametric elements.
// ------------------------------------------------------------------------
//
// The ElmParam class defines the behavior parametric finite elements.
// In order to allow a generic implementation, each parametric element has
// pointers to an Analysis Model and to a Shape object. The former handles
// the features related to the differential equation/variational statement
// governing the problem (e.g. number dofs/node and form of constitutive
// and strain-displacement matrices) while the latter handles the field and
// geometry interpolation aspects (e.g. element topology, connectivity and
// shape functions).
//
// Each parametric element also has an array Integration Points and an
// array of Constitutive Models, where each Constitutive Model corresponds
// to one Integration Point. The Constitutive Model is resposible by the
// computation of the stress vector and the tangent constitutive matrix at
// an Integration Point, while the Integration Point stores the coordinates
// and weights used in the numerical integration of the stiffness matrix
// and the internal force vector. The numerical integration of the mass
// matrix is performed using a different (generally higher) integration
// order and its Integration Points are not stored, but created/destroyed
// inside the MassMat method.
//
// The nodal stresses of parametric elements are obtained by the
// extrapolation of stresses computed at the integration points used in
// numerical integration of the stiffness matrix. Currenctly, this
// extrapolation (smoothing) is carried out building a linear approximation
// of the stress field (in parametric coordinates) by using the least
// square method.
// ------------------------------------------------------------------------

#ifndef _ELMPARAM_H
#define _ELMPARAM_H

#include <string>
#include <vector>

#include "element.h"
#include "anmodel.h"
#include "cmodel.h"
#include "intpoint.h"
#include "mat.h"
#include "section.h"
#include "secanalysis.h"

// ------------------------------------------------------------------------
// Definition of the Parametric Element class:
//
class cElmParam : public cElement
{
 protected:
  static cMatrix       TR;          // Stress extrapolation matrix
  static int           PrvOrdIdx;   // Int. order of previous element
  static eShpType      PrvShpType;  // Shape type of the previous element
         int           OrdIdx;      // Integration order index
         int           NumIntPnt;   // Number of integration points
         cAnModel     *Model;       // Pointer to element analysis model
         cIntPoint    *IntPnt;      // Array of integration points
         cSecAnalysis **SecAn;      // Array of section analyses
         cVector       IniStr;      // Initial element stresses
         cMatrix       StrIpt;      // Integration point stress

 protected:
          void      CompTRMatrix(int n = 0, cIntPoint *ipt = 0);
          void      NearestIntPnt(int, sNatCoord *, int, cIntPoint *, int, 
                                  vector<int> &);

 public:

  template <typename T>
  static  void ReadElm(string,eShpType,eAnmType);

  static  void      ReadElmFEM(string,eShpType,eAnmType,bool tl = 0);
  static  void      ReadElmBSP(std::string,eShpType,eAnmType, bool tl = 0);


                    cElmParam(int, eShpType, eAnmType);
  virtual          ~cElmParam(void);
          void      UpdateState(void);
          void      Read(void);
          void      ReadIniStr(void);
          int       GetNumIntPnt(void) { return NumIntPnt; }
          int       GetIntOrd(void);
          int       GetIntOrdIdx(void) { return OrdIdx; }
          int       GetIntType(void){ return VecIntOrd[OrdIdx].type;}
          int       GetIntOrdZ(void){ return VecIntOrd[OrdIdx].K[2];}
          cAnModel *GetAnModel(void) { return Model; }
          cSection *GetSection(void) { return Section; }
          int       GetNumDofNode(void) { return Model->GetNumDofNode( ); }
          int       GetNumStrCmp(void) { return Model->GetDimBMatrix( ); } // Evandro: verificar.
          int       GetNumNodStrCmp(void) { return Model->GetNumNodStrCmp( ); }
	  void      GetIntPnt(int,sNatCoord&,double&);
          void      GetStrLabels(int *lab) { Model->GetStrLabels(lab); }
          void      GetNodStrLabels(int *lab) { Model->GetNodStrLabels(lab); }
          void      GetActDir(int *dir) { Model->GetActDir(dir); }
          double    VolCoeff(double t, int n, double *mapfnc, sNodeCoord *coord)
                    { return Model->VolCoeff(t, n, mapfnc, coord); }
          double    IntFunc(int,cIntPoint*,cVector&,double&);
          double    IntFunc(cVector&,double&);
  virtual int       IntForce(cVector &);
  virtual void      StiffMat(cMatrix &);
  virtual void      ConducMat(cMatrix &);
  virtual void      MassMat(cMatrix &);
  virtual void      GeomStiff(cMatrix &);
  virtual void      NodalStress(cMatrix &);
  virtual void      IntPntTemp(cVector&);
  virtual void      IntPntStress(cMatrix &, cMatrix &);
  virtual void      IntPntFlux(cMatrix &);
  virtual void      PntStress(sNatCoord*,int,cMatrix&,cMatrix&,bool opt = 0);
  virtual void      CalcStrInt(double &, cVector &);
  virtual void      InitialStresses(sNatCoord, cVector &);
  virtual  int      ThermalEffects( void ) { return 1; }
  virtual  int      GetStrnTemp(sNatCoord,double *,cVector & );
};

// ------------------------------------------------------------------------
// Definition of the Parametric Total Lagrangian Element class:
//
class cElmParamTL : public cElmParam
{
 protected:
          void      GreenStrain(cMatrix &, cMatrix &, cVector &, cVector &);

 public:

                    cElmParamTL(int, eShpType, eAnmType);
  virtual          ~cElmParamTL(void) {}

          int       IntForce(cVector &);
          void      StiffMat(cMatrix &);
          void      IntPntStress(cMatrix &, cMatrix &);
};

#endif
