// ------------------------------------------------------------------------
// secbar.h - definition of cross-section class for bar elements.
// ------------------------------------------------------------------------
//
// The SecBar class defines the general behavior of a cross-section in the 
// structural analysis using bar (truss and beam/frame) elements and 
// implements methods shared by different section types. The specific 
// features of the each section type are defined by pure virtual methods 
// implemented in the derived classes. The BarSection class is also 
// responsible for the storage of the sections read from the input file.
// Class hierarchy:
//
// SecBar
// |-- SecBarHom
// |   |-- SecBarGen
// |   |-- SecBarRect
// |   |-- SecBarCirc
// |   |-- SecBarTube
// |-- SecBarB3DCoupled
// |-- SecBarRC
// |   |-- SecRCRect
//
// The geometric properties, stress resultants and section stiffness are
// always referenced to the Centroid (C) of the cross-section. Depending
// on the materials/constitutive models of the section the integration
// of stresses and stiffnesses can be carried out before the analysis (pre-
// integration) or during the analysis. The former is used for sections 
// with linear elastic materials and latter is used in the other cases.
//
// SecBarGen corresponds to a cross-section of unknown shape whose geometric
// properties are given in the input file. This type of section cannot 
// be used in the materially nonlinear analysis, since the integration is
// always carried out before the analysis (pre-integration).
//
// SecBarB3DCoupled should be used in the materially linear analysis of 
// laminated composite beams with fully coupled (4x4) constitutive matrix.
// The geometric properties and section constitutive matrix are read from
// the input file, thus the associated (elastic) material is used only for
// computation of mass and weight.
//
// SecBarHom corresponds to a cross-section of given shape and only one 
// (isotropic) material. It can be used in both linear and nonlinear 
// materials. Section integration is carried out before the analysis (pre-
// integration) for linear materials and during the analysis for the other
// cases.
//
// SecBarRC corresponds to a reinforced concrete (RC) cross-section of given
// shape. This section is composed of an homogeneous material representing
// the concrete and a set of rebars. Different materials/constitutive 
// models can be used to model the concrete and the rebars. Integration is
// always carried out during the analysis. The position of rebars are defi-
// ned in a local (x, y) coordinate system whose origin depend on the
// section shape. For double symmetric sections (e.g. rectangles and
// circles) the origin of this system is on the Centroid (C) of concrete 
// section.
// ------------------------------------------------------------------------

#ifndef _SECBAR_H
#define _SECBAR_H

#include "section.h"
#include "secanalysis.h"

// ------------------------------------------------------------------------
// Forward declarations:
//
class cVector;
class cMatrix;
class cMaterial;
class cConstModel;

// ------------------------------------------------------------------------
// Definition of the beam section class:
//
class cSecBar : public cSection
{
 protected:
  double A;          // Area
  double Ay,Az;      // Shear areas (Timoshenko beam)
  double Iy,Iz;      // Moments of inertia
  double Jt;         // Torsion contant
  double Cw;         // Warping constant

 public:
                  cSecBar(int);
  virtual        ~cSecBar(void);
  virtual double  GetA(void) { return A; }
  virtual double  GetAy(void) { return Ay; }
  virtual double  GetAz(void) { return Az; }
  virtual double  GetIy(void) { return Iy; }
  virtual double  GetIz(void) { return Iz; }
  virtual double  GetJt(void) { return Jt; }
  virtual double  GetCw(void) { return Cw; }
  virtual void    Read(void) = 0;
  virtual void    SecStiffness2D(cMatrix &){ }
  virtual void    SecStiffness3D(cMatrix &){ }
};


// ------------------------------------------------------------------------
// Definition of the Homogeneous section class:
//
class cSecBarHom : public cSecBar
{
 protected:
  cMaterial *Mat;    // Section material

 public:
                     cSecBarHom(int);
  virtual           ~cSecBarHom(void);
          cMaterial *GetMaterial(void) { return Mat; }
          bool       PreIntegrated(void);
          void       Read(void);
  virtual void       ReadGeom(void) = 0;
  virtual void       SecStiffness2D(cMatrix &);
  virtual void       SecStiffness3D(cMatrix &);
};


// ------------------------------------------------------------------------
// Definition of the General section class:
//
class cSecBarGen : public cSecBarHom
{
 public:
                cSecBarGen(int);
  virtual      ~cSecBarGen(void);
          bool  PreIntegrated(void) { return true; }
          void  ReadGeom(void);
};


// ------------------------------------------------------------------------
// Definition of the Rectangular section class:
//
class cSecBarRect : public cSecBarHom
{
 friend class cSecAnPlFrameRectFiber;
 friend class cSecAnPlFrameRectGauss;

 protected:
  double b;      // Width
  double h;      // Height

 public:
                cSecBarRect(int);
  virtual      ~cSecBarRect(void);
          void  ReadGeom(void);
};


// ------------------------------------------------------------------------
// Definition of the Circular section class:
//
class cSecBarCirc : public cSecBarHom
{
 protected:
  double R;      // Radius

 public:
                cSecBarCirc(int);
  virtual      ~cSecBarCirc(void);
          void  ReadGeom(void);
};


// ------------------------------------------------------------------------
// Definition of the Tubular section class:
//
class cSecBarTube : public cSecBarHom
{
 protected:
  double R;      // Outer radius
  double t;      // Thickness

 public:
                cSecBarTube(int);
  virtual      ~cSecBarTube(void);
          void  ReadGeom(void);
};


// ------------------------------------------------------------------------
// Definition of the Fully Coupled 3D Beam section class:
//
class cSecBarB3DCoupled : public cSecBar
{
 protected:
  cMaterial *Mat;   // Section material
  cMatrix    C;     // Constitutive matrix (4x4)

 public:
                cSecBarB3DCoupled(int);
  virtual      ~cSecBarB3DCoupled(void);
          bool  PreIntegrated(void) { return true; }
          void  Read(void);
          void  SecStiffness2D(cMatrix &);
          void  SecStiffness3D(cMatrix &);
};


// ------------------------------------------------------------------------
// Definition of the Reinforced Concrete section class:
//
class cSecBarRC : public cSecBar
{
 protected:
  int         NumBar;   // Number of rebars
  double     *XBar;     // X-coordinate of each rebar
  double     *YBar;     // Y-coordinate of each rebar
  double     *ABar;     // Area of each rebar
  cMaterial  *Conc;     // Concrete material
  cMaterial **MatBar;   // Material of each rebar

 public:
                      cSecBarRC(int);
  virtual            ~cSecBarRC(void);
          cMaterial  *GetMaterial(void) { return Conc; }
          bool        PreIntegrated(void) { return false; }
          void        Read(void);
  virtual void        ReadGeom(void) = 0;
  virtual void        GetCentroid(double *, double *) = 0;
};


// ------------------------------------------------------------------------
// Definition of the Rectangular Reinforced Concrete section class:
//
class cSecBarRCRect : public cSecBarRC
{
 friend class cSecAnPlFrameRCRect;

 protected:
  double b;      // Width
  double h;      // Height

 public:
                cSecBarRCRect(int);
  virtual      ~cSecBarRCRect(void);
          void  ReadGeom(void);
          void  GetCentroid(double *, double *);
};

#endif
