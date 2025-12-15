// ------------------------------------------------------------------------
// material.h - file containing the definition of the Material class.
// ------------------------------------------------------------------------
//
// The Material class basically reads the material parameters from the
// input file. These parameters are stored and can be queried by others
// when these parameters are necessary.
//
// Material
// |-- Elastic (Isotropic)
// |   |-- Orthotropic
// |   |-- ElastNonlin
// |   |-- FGM
// |   |-- ConcreteEuroCEB
// |   |-- ConcreteNBRCEB
// |-- ViscoElastic
// |-- PlastIsoLinHard
// |-- PlastLinHard
// |-- PlastIsoNonlinHard
// |-- PlastIsoExpHard
// |-- ProgFailHashin
// |-- ProgFailTsai
// |-- ProgFailEngelstad
// |-- ProgFailEngelstad
// |-- CZMModeIBilinear
// |-- CZMModeIExponential
// |-- Thermal Isotropic
// |-- Thermal Orthotropic
//
// -------------------------------------------------------------------------
// Static public methods:
// -------------------------------------------------------------------------
//
// static void ReadNumMat(void)
//
// This method reads the number of materials and creates an array to store
// these materials.
// -------------------------------------------------------------------------
//
// static void ReadElastic(void)
// static void ReadElastOrtho(void)
// static void ReadElastNonlin(void)
// static void ReadFGM(void)
// static void ReadConcEuroCEB(void)
// static void ReadConcNBRCEB(void)
// static void ReadViscoElastic(void)
// static void ReadPlastIsoLinHard(void);
// static void ReadPlastLinHard(void)
// static void ReadPlastIsoNonlinHard(void)
// static void ReadPlastIsoExpHard (void)
// static void ReadDamageMazars(void)
// static void ReadDamageMuModel(void)
// static void ReadLeeFenves(void)
// static void ReadProgFailHashin(void)
// static void ReadProgFailTsai(void)
// static void ReadProgFailEngelstad(void)
// static void ReadCZMModeIBilinear(void)
// static void ReadCZMModeIExponential(void)
//
// Create and read each type of material.
// -------------------------------------------------------------------------
//
// static void ReadDensity(void)
//
// This method reads the mass density of the given materials.
// -------------------------------------------------------------------------
//
// static void Destroy(void)
//
// This method destroys all stored materials.
// -------------------------------------------------------------------------
//
// static cMaterial *GetMaterial(int id)
//
//  id     -  material label                                          (in)
//
// This method returns a pointer to the material corresponding to the given
// label.
//
// -------------------------------------------------------------------------
// Shared methods (implemented only in the base class):
// -------------------------------------------------------------------------
//
// eMatType GetType(void)
//
// This method returns the material type.
//
// -------------------------------------------------------------------------
//
// int GetLabel(void)
//
// This method returns the material label.
//
// -------------------------------------------------------------------------
//
// double GetDensity(void)
//
// This method returns the material mass density.
//
// -------------------------------------------------------------------------
// Pure virtual methods:
// -------------------------------------------------------------------------
//
// virtual int NumParam(void)
//
// This method returns the number of material parameters.
// -------------------------------------------------------------------------
//
// virtual void GetParam(double *param)
//
//  param  -  array of parameters                                    (out)
//
// This method returns all the material parameters in the given array.
// -------------------------------------------------------------------------
//
// virtual void Read(void)
//
// This method reads the material parameters.
// -------------------------------------------------------------------------

#include "cmodel1d.h"

#ifndef _MATERIAL_H
#define _MATERIAL_H

// -------------------------------------------------------------------------
// Foward declaration:
//
class cMatProp;
class cMatMec;
class cMatTherm;

// -------------------------------------------------------------------------
// Mechanical Material types:
//
typedef enum
{
  MAT_ELASTIC_ISOTROPIC,     // Elastic isotropic material
  MAT_ELASTIC_ORTHOTROPIC,   // Elastic orthotropic material
  MAT_ELASTIC_NONLINEAR,     // Nonlinear elastic isotropic material
  MAT_FGM,                   // Functionally Graded Material
  MAT_FGM_TTO,               // Functionally Graded Material with TTO model
  MAT_CONCRETE_EUROCEB,      // Nonlinear elastic EuroCEB concrete material
  MAT_CONCRETE_NBRCEB,       // Nonlinear elastic NBRCEB concrete material
  MAT_VISCOELASTIC,          // Viscoelastic (isotropic) material
  MAT_PLASTIC_ISOLINHARD,    // Elastoplastic - Von Mises Multidimensional Lin. Hardening
  MAT_PLASTIC_LINHARD,       // Elastoplastic - iso/kin linear hardening
  MAT_PLASTIC_ISONONLINHARD, // Elastoplastic - iso. nonlinear hardening
  MAT_PLASTIC_ISOEXPHARD,    // Elastoplastic - iso. exponential hardening
  MAT_DAMAGE_MAZARS,         // Damage material - Mazars (1985)
  MAT_DAMAGE_MU_MODEL,       // Damage material - Mazars (2015)
  MAT_LEE_FENVES,            // Plastic Damage mat - Lee & Fenves (1998)
  MAT_PROGFAIL_HASHIN,       // FRC with prog. failure - Hashin
  MAT_PROGFAIL_TSAI,         // FRC with prog. failure - Tsai
  MAT_PROGFAIL_ENGELSTAD,    // FRC with prog. failure - Engelstad
  MAT_CZM_MODEI_BILINEAR,    // Cohesive Zone Bilinear Law (Mode I)
  MAT_CZM_MODEI_EXPONENTIAL  // Cohesive Zone Exponential Law (Mode I)
} eMatMecType;

// -------------------------------------------------------------------------
// Thermal Material types:
//
typedef enum
{
  MAT_THERMAL_ISOTROPIC,       // Thermal isotropic material
  MAT_THERMAL_ORTHOTROPIC      // Thermal orthotropic material
} eMatThermType;

// ------------------------------------------------------------------------
// Definition of the Material class:
//
class cMaterial
{
 private:
  static int         NumMat;
  static cMaterial  *VecMat;

 protected:
  cMatMec   *MecProp;   // Mechanical properties
  cMatTherm *ThermProp; // Thermal properties
  int        Label;     // Material label
  cVector    Density;   // Mass density

         void       SetMecProp(cMatMec*);
         void       SetThermProp(cMatTherm*);

 public:
  static void       ReadNumMat(void);
  static void       ReadElastIso(void);
  static void       ReadElastOrtho(void);
  static void       ReadElastNonlin(void);
  static void       ReadFGM(void);
  static void       ReadFGMTTO(void);  
  static void       ReadViscoElastic(void);
  static void       ReadConcEuroCEB(void);
  static void       ReadConcNBRCEB(void);
  static void       ReadPlastIsoLinHard(void);
  static void       ReadPlastLinHard(void);
  static void       ReadPlastIsoNonlinHard(void);
  static void       ReadPlastIsoExpHard(void);
  static void       ReadDamageMazars(void);
  static void       ReadDamageMuModel(void);
  static void       ReadLeeFenves(void);
  static void       ReadProgFailHashin(void);
  static void       ReadProgFailTsai(void);
  static void       ReadProgFailEngelstad(void);
  static void       ReadCZMModeIBilinear(void);
  static void       ReadCZMModeIExponential(void);
  static void       ReadThermalIso(void);
  static void       ReadThermalOrtho(void);
  static void       ReadDensity(void);
  static void       ReadThermExpan(void);
  static void       Destroy(void);
  static cMaterial *GetMaterial(int);

                     cMaterial(void);
                    ~cMaterial(void);
          int        GetLabel(void)           { return Label;      }
          double     GetDensity(void)         { return Density[0]; }
          void       GetDensity(cVector &rho) { rho = Density;     }
          cMatMec*   GetMecProp(void)         { return MecProp;    }
          cMatTherm* GetThermProp(void)       { return ThermProp;  }

  friend class cMatProp;
};

// ------------------------------------------------------------------------
// Definition of the Material properties interface:
//
class cMatProp
{
 protected:
  cMaterial *mat;

  cVector& GetDensity(void) { return mat->Density; }

 public:
                       cMatProp(void) { }
  virtual             ~cMatProp(void) { }
  virtual int          NumParam(void) = 0;
  virtual void         GetParam(double *) = 0;
  virtual void         Read(void) = 0;
  virtual void         ReadAlpha(void) { }

  friend class cMaterial;
};

// ------------------------------------------------------------------------
// Definition of the Material Mechanical Properties class:
//
class cMatMec : public cMatProp
{
 protected:
  eMatMecType Type;    // Material type

 public:

                       cMatMec(void){ }
  virtual             ~cMatMec(void){ }
          eMatMecType  GetType(void) { return Type; }
};


// ------------------------------------------------------------------------
// Definition of the Material Thermal Properties class:
//
class cMatTherm : public cMatProp
{
 protected:
  eMatThermType Type;    // Material type

 public:

                        cMatTherm(void){ }
  virtual              ~cMatTherm(void){ }
          eMatThermType GetType(void) { return Type; }
};


// ------------------------------------------------------------------------
// Definition of the Elastic Isotropic Material class:
//
class cMatElastIso : public cMatMec
{
 protected:
  double E;     // Elastic modulus
  double Nu;    // Poisson's ratio
  double Alpha; // Thermal expansion coefficient

 public:
                  cMatElastIso(void);
  virtual        ~cMatElastIso(void);
           void   Read(void);
           void   ReadAlpha(void);
           void   GetParam(double *);
           int    NumParam(void) { return 3; }
};


// ------------------------------------------------------------------------
// Definition of the Elastic Orthotropic Material class:
//
class cMatElastOrtho : public cMatMec
{
 protected:
  double E1,E2,E3;       // Elastic moduli in the material directions
  double Nu12,Nu13,Nu23; // Poisson's ratios
  double G12,G13,G23;    // Shear moduli
  double Alpha1,Alpha2;  // Thermal expansion coefficients

 public:
                  cMatElastOrtho(void);
  virtual        ~cMatElastOrtho(void);
           void   Read(void);
           void   ReadAlpha(void);
           void   GetParam(double *);
            int   NumParam(void) { return 11; }
};


// ------------------------------------------------------------------------
// Definition of the Nonlinear Elastic Isotropic Material class:
//
class cMatElastNonlin : public cMatMec
{
 friend class cMod1DElastNonlin;

 protected:
  int    NumPnt; // Number of points
  double *Eps;   // Strains
  double *Sig;   // Stresses
  double Nu;     // Poisson's ratio

 public:
                cMatElastNonlin(void);
  virtual      ~cMatElastNonlin(void);
          bool  IsNonlinear(void) { return true; }
          void  Read(void);
          void  GetParam(double *);
          int   NumParam(void) { return(2*NumPnt + 1); }
};


// ------------------------------------------------------------------------
// Definition of the Functionally Graded Material class:
//
class cMatFGM : public cMatMec
{
 protected:
  int    HomMethod;     // Homogenization method (1 = RoM, 2 = Mori-Tanaka)
  double E1,E2;         // Elastic moduli
  double Nu1,Nu2;       // Poisson's ratios
  double Alpha1,Alpha2; // Thermal expansion coefficients

 public:
                  cMatFGM(void);
  virtual        ~cMatFGM(void);
           void   Read(void);
           void   ReadAlpha(void);
           void   GetParam(double *);
           int    NumParam(void) { return 7; }
};

// ------------------------------------------------------------------------
// Definition of the FGM TTO Material Class:
//
class cMatFGMTTO : public cMatMec
{
 protected:
  int    HomMethod;     // Homogenization method (1 = RoM, 2 = Mori-Tanaka)
  double E1,E2;         // Elastic moduli
  double Nu1,Nu2;       // Poisson's ratios
  double SigYm;         // Matrix material yield stress
  double Hm;            // Matrix material tangent modulus
  double qtto;          // Stress transfer coefficient from TTO method

 public:
                  cMatFGMTTO(void);
  virtual        ~cMatFGMTTO(void);
           void   Read(void);
           void   GetParam(double *);
           int    NumParam(void) { return 8; }
};

// ------------------------------------------------------------------------
// Definition of the Prony series.
//
typedef struct
{
  double Eoo;
  int n;
  double *E;
  double *rho;
} sProny;

// ------------------------------------------------------------------------
// Definition of the ViscoElastic Material class:
//
class cMatViscoElastic : public cMatMec
{
 protected:
  sProny  RelMod;  // Relaxation modulus: E(t) = Eoo + E(i)*exp(-t/rho(i))
  double  Nu;      // Poisson's ratio

 public:
                cMatViscoElastic(void);
  virtual      ~cMatViscoElastic(void);
          void  Read(void);
          int   NumParam(void) { return(2*RelMod.n + 2); }
          void  GetParam(double *);
};


// ------------------------------------------------------------------------
// Definition of the EuroCEB Concrete Material class:
//
class cMatConcEuroCEB : public cMatMec
{
 friend class cMod1DConcEuroCEB;

 protected:
  double fcm;    // Compressive strength of concrete
  double ec1;    // Strain at peak stress for compression of concrete
  double Ecm;    // Secant modulus from the origin to 0.4fcm (compresion)
  double fct;    // Tensile strength of concrete
  double Eci;    // Initial tangent modulus (tension)
  double ey;     // Yield strain of steel
  double Es;     // Young's modulus of steel
  double Rho;    // Reinforcement ratio

 public:
                cMatConcEuroCEB(void);
  virtual      ~cMatConcEuroCEB(void);
          bool  IsNonlinear(void) { return true; }
          void  Read(void);
          void  GetParam(double *);
          int   NumParam(void) { return 8; }
};


// ------------------------------------------------------------------------
// Definition of the NBRCEB Concrete Material class:
//
class cMatConcNBRCEB : public cMatMec
{
 friend class cMod1DConcNBRCEB;

 protected:
  double fc;     // Compressive strength
  double fct;    // Tensile strength
  double Eci;    // Initial tangent modulus for tensile strains
  double ey;     // Yield strain of steel
  double Es;     // Young's modulus of steel
  double Rho;    // Reinforcement ratio

 public:
                cMatConcNBRCEB(void);
  virtual      ~cMatConcNBRCEB(void);
          bool  IsNonlinear(void) { return true; }
          void  Read(void);
          void  GetParam(double *);
          int   NumParam(void) { return 6; }
};

// ------------------------------------------------------------------------
// Definition of the ElastoPlastic Material (von Mises) with Linear
// Isotropic Hardening class:
//
class cMatPlastIsoLinHard : public cMatMec
{
 protected:
  double E;     // Elastic modulus
  double Nu;    // Poisson's ratio
  double SigY0; // Initial yield Stress
  double H;     // Isotropic hardening modulus

 public:
                  cMatPlastIsoLinHard(void);
  virtual        ~cMatPlastIsoLinHard(void);
           void   Read(void);
           void   GetParam(double *);
           int    NumParam(void) { return 4; }
};

// ------------------------------------------------------------------------
// Definition of the ElastoPlastic Material with Linear Isotropic and
// Kinematic Hardening class:
//
class cMatPlastLinHard : public cMatMec
{
 friend class cMod1DPlastLinHard;

 protected:
  double E;     // Elastic modulus
  double Nu;    // Poisson's ratio
  double SigY;  // Yield Stress
  double K;     // Isotropic hardening modulus
  double H;     // Kinematic hardening modulus

 public:
                cMatPlastLinHard(void);
  virtual      ~cMatPlastLinHard(void);
          void  Read(void);
          void  GetParam(double *);
          int   NumParam(void) { return 5; }
};


// ------------------------------------------------------------------------
// Definition of the ElastoPlastic Material with Nonlinear Isotropic
// Hardening class:
//
class cMatPlastIsoNonlinHard : public cMatMec
{
 friend class cMod1DPlastIsoNonlinHard;

 protected:
  double E;     // Elastic modulus
  double Nu;    // Poisson's ratio
  double SigY;  // Yield Stress
  double K;     // Plastic Modulus
  double ExpH;  // Plastic Exponent

 public:
                cMatPlastIsoNonlinHard(void);
  virtual      ~cMatPlastIsoNonlinHard(void);
          void  Read(void);
          void  GetParam(double *);
          int   NumParam(void) { return 5; }
};


// ------------------------------------------------------------------------
// Definition of the ElastoPlastic Material with Exponential Isotropic
// Hardening class:
//
class cMatPlastIsoExpHard : public cMatMec
{
 friend class cMod1DPlastIsoExpHard;

 protected:
  double E;     // Elastic modulus
  double Nu;    // Poisson's ratio
  double SigY;  // Yield Stress
  double Mu;    // Hardening parameter

 public:
                cMatPlastIsoExpHard(void);
  virtual      ~cMatPlastIsoExpHard(void);
          void  Read(void);
          void  GetParam(double *);
          int   NumParam(void) { return 4; }
};


// ------------------------------------------------------------------------
// Definition of the Mazars Damage Model for concrete class:
//
class cMatDamageMazars : public cMatMec
{
 friend class cMod1DDamageMazars;

 protected:
  double E;     // Elastic Modulus
  double Ac;    // Ac parameter
  double Bc;    // Bc parameter
  double At;    // At parameter
  double Bt;    // Bt parameter
  double epsd0; // epsd0 parameter
  double Nu;    // Poisson's ratio

 public:
                cMatDamageMazars(void);
  virtual      ~cMatDamageMazars(void);
          void  Read(void);
          void  GetParam(double *);
          int   NumParam(void) { return 7; }
};


// ------------------------------------------------------------------------
// Definition of the Damage Mu-Model for concrete class:
//
class cMatDamageMuModel : public cMatMec
{
 friend class cMod1DDamageMuModel;

 protected:
  double E;     // Elastic Modulus
  double Ac;    // Ac parameter
  double Bc;    // Bc parameter
  double At;    // At parameter
  double Bt;    // Bt parameter
  double epsdt0; // epsdt0 parameter
  double epsdc0; // epsdc0 parameter
  double Nu;    // Poisson's ratio

 public:
                cMatDamageMuModel(void);
  virtual      ~cMatDamageMuModel(void);
          void  Read(void);
          void  GetParam(double *);
          int   NumParam(void) { return 8; }
};


// ------------------------------------------------------------------------
// Definition of the Plastic Damage Model for concrete by Lee & Fenves class:
//
class cMatLeeFenves : public cMatMec
{
 friend class cMod1DLeeFenves;

 protected:
  double E;     // Elastic Modulus
  double fb0;   // biaxial compressive strength
  double fc0;   // uniaxial compressive strength
  double ft0;   // uniaxial tensile strength
  double a_c;   // ac parameter
  double b_c;   // bc parameter
  double d_c;   // dc parameter
  double a_t;   // at parameter
  double b_t;   // bt parameter
  double d_t;   // dt parameter
  double Gt;    // Fracture energy
  double Gc;    // Crushing energy
  double lt;    // Characteristic length for tensile
  double lc;    // Characteristic length for compression
  double s0;    // alpha_p parameter

 public:
                cMatLeeFenves(void);
  virtual      ~cMatLeeFenves(void);
          void  Read(void);
          void  GetParam(double *);
          int   NumParam(void) { return 15; }
};


// ------------------------------------------------------------------------
// Definition of the ProgFailHashin Material class:
//
class cMatProgFailHashin : public cMatMec
{
 protected:
  double E1,E2,E3;                      // Elastic moduli
  double Nu12,Nu13,Nu23;                // Poisson's ratios
  double G12,G13,G23;                   // Shear moduli
  double Xt,Xc,Yt,Yc,Zt,Zc,S12,S13,S23; // Strength parameters
  double Alpha;                         // Degradation parameter

 public:
                cMatProgFailHashin(void);
  virtual      ~cMatProgFailHashin(void);
          void  Read(void);
          int   NumParam(void) { return 19; }
          void  GetParam(double *);
};


// ------------------------------------------------------------------------
// Definition of the ProgFailTsai Material class:
//
class cMatProgFailTsai : public cMatMec
{
 protected:
  double E1,E2,E3;                      // Elastic moduli
  double Nu12,Nu13,Nu23;                // Poisson's ratios
  double G12,G13,G23;                   // Shear moduli
  double Xt,Xc,Yt,Yc,Zt,Zc,S12,S13,S23; // Strength parameters

 public:
                cMatProgFailTsai(void);
  virtual      ~cMatProgFailTsai(void);
          void  Read(void);
          int   NumParam(void) { return 18; }
          void  GetParam(double *);
};


// ------------------------------------------------------------------------
// Definition of the ProgFailEngelstad Material class:
//
class cMatProgFailEngelstad : public cMatMec
{
 protected:
  double E1,E2,E3;                      // Elastic moduli
  double Nu12,Nu13,Nu23;                // Poisson's ratios
  double G12,G13,G23;                   // Shear moduli
  double Xt,Xc,Yt,Yc,Zt,Zc,S12,S13,S23; // Strength parameters
  double Alpha;                         // Degradation parameter

 public:
                cMatProgFailEngelstad(void);
  virtual      ~cMatProgFailEngelstad(void);
          void  Read(void);
          int   NumParam(void) { return 19; }
          void  GetParam(double *);
};


// ------------------------------------------------------------------------
// Definition of the CZMModeIBilinear Material class:
//
class cMatCZMModeIBilinear : public cMatMec
{
 protected:
  double Ft;    // Tensile strength
  double GIc;   // Fracture energy
  double K0;    // Penalty stiffness

 public:
                cMatCZMModeIBilinear(void);
  virtual      ~cMatCZMModeIBilinear(void);
          void  Read(void);
          int   NumParam(void) { return 3; }
          void  GetParam(double *);
};


// ------------------------------------------------------------------------
// Definition of the CZMModeIExponential Material class:
//
class cMatCZMModeIExponential : public cMatMec
{
 protected:
  double Ft;    // Tensile strength
  double GIc;   // Fracture energy

 public:
                cMatCZMModeIExponential(void);
  virtual      ~cMatCZMModeIExponential(void);
          void  Read(void);
          int   NumParam(void) { return 2; }
          void  GetParam(double *);
};


// ------------------------------------------------------------------------
// Definition of the Thermal Isotropic Material class:
//
class cMatThermalIso : public cMatTherm
{
 protected:
  double k;   // Thermal conductivity

 public:
                  cMatThermalIso(void);
  virtual        ~cMatThermalIso(void);
           void   Read(void);
           void   GetParam(double *);
           int    NumParam(void) { return 1; }
};


// ------------------------------------------------------------------------
// Definition of the Thermal Orthotropic Material class:
//
class cMatThermalOrtho : public cMatTherm
{
 protected:
  double k1,k2,k3;       // Thermal conductivities

 public:
                  cMatThermalOrtho(void);
  virtual        ~cMatThermalOrtho(void);
           void   Read(void);
           void   GetParam(double *);
            int   NumParam(void) { return 3; }
};


#endif
