// -------------------------------------------------------------------------
// section.h - file containing the definition of the Section class.
// -------------------------------------------------------------------------
//
// The Section class deals with the different element section types, for
// both parametric and bar elements. Each section reads and stores specific
// data (geometry and materials) used by the cSecAnalysis class methods,
// allowing the analysis of homogeneous and composite (e.g. FRC laminates
// and reinforced concrete) structures.
//
// Section
// |-- SecHomIso2D => Plane Stress, Strain, Axisymmetric, Plate, and Shell
// |-- SecHomIso3D => Solid and DegShell
// |-- SecHomOrtho2D => Plane Stress
// |-- SecHomOrtho3D => Solid
// |-- SecFGM2D => Plane Stress, Strain, and Axisymmetric 
// |-- SecFGM3D => Solid
// |-- SecHSDTPlate
// |-- SecHSDTFGMPlate
// |-- SecFGMShell
// |-- SecGeneralShell
// |-- SecLaminated2D => Plane Stress and Plate (symmetric layup)
// |-- SecLaminatedShell
// |-- SecLaminatedSolid
// |-- SecInterface2D
// |-- SecInterface3D
// |-- SecBar  (see secbar.h)
//
// Important: 
// 1) Shallow and Donnell shell elements use 2D sections.
// 2) Degenerated solid shell elements use 3D sections.
// -------------------------------------------------------------------------
// Static public methods:
// -------------------------------------------------------------------------
//
// static void ReadNumSec(void)
//
// This method reads the number of sections and creates an array to store
// them.
// -------------------------------------------------------------------------
//
// static void ReadHomIso2D(void)
// static void ReadHomIso3D(void)
// static void ReadHomOrtho2D(void)
// static void ReadHomOrtho3D(void)
// static void ReadFGM2D(void)
// static void ReadFGM3D(void)
// static void ReadHSDTPlate(void)
// static void ReadHSDTFGMPlate(void)
// static void ReadFGMShell(void)
// static void ReadGeneralShell(void)
// static void ReadLaminated2D(void)
// static void ReadLaminatedShell(void)
// static void ReadLaminatedSolid(void)
// static void ReadInterface2D(void)
// static void ReadInterface3D(void)
// static void ReadSecGen(void)
// static void ReadSecB3DCoupled(void)
// static void ReadSecRect(void)
// static void ReadSecCirc(void)
// static void ReadSecTube(void)
// static void ReadSecRCRect(void)
//
// Each of these methods reads and creates sections.
// -------------------------------------------------------------------------
//
// static void Destroy(void)
//
// This method destroys all stored sections.
// -------------------------------------------------------------------------
//
// static cSection *GetSection(int id)
//
//  id     -  section label                                          (in)
//
// This method returns a pointer to the section corresponding to the given
// label.
//
// -------------------------------------------------------------------------
// Shared methods (implemented only in the base class):
// -------------------------------------------------------------------------
//
// eSecType GetType(void)
//
// This method returns the section type.
// -------------------------------------------------------------------------
//
// int GetLabel(void)
//
// This method returns the section label.
//
// -------------------------------------------------------------------------
// Virtual methods:
// -------------------------------------------------------------------------
//
// bool PreIntegrated(void)
//
// This function returns true if the section properties are computed before
// the analysis (pre-integration) and false if these properties are computed
// during the analysis.
// -------------------------------------------------------------------------
//
// int GetNumLam(void)
//
// This function returns the number of section laminas (plies) for laminated
// sections and 0 for homogeneous sections.
// -------------------------------------------------------------------------
//
// double GetThickness(void)
//
// This function returns the total thickness for 2D and shell sections and
// 1.0 for 3D solid sections.
// -------------------------------------------------------------------------
//
// double GetDensity(void)
//
// This function returns the material density (specific mass) for homogenous
// sections and the average material density for laminated sections.
// -------------------------------------------------------------------------
//
// cMaterial *GetMaterial(void)
//
// This function returns the section material for homogenous sections and
// 0 (null pointer) for laminated sections.
// -------------------------------------------------------------------------
//
// sLamina *GetLayup(void)
//
// This function returns the layup (laminate data) for laminated sections
// and 0 (null pointer) for homogeneous sections.
// -------------------------------------------------------------------------
//
// double GetOrientationAngle(void)
//
// This method returns the section orientation angle (i.e. angle between
// the material x-axis) for homogenous orthotropic plane (2D/plate) sections
// and 0.0 for other sections.
// -------------------------------------------------------------------------
//
// void GetOrientationVector(cVector &v)
//
//  v   -  orientation vector                                        (out)
//
// This method returns the section orientation vector corresponding to
// the material x-axis for homogenous orthotropic sections and the laminate
// x-axis for laminated shell sections. It returns the vector {1, 0, 0}
// for other sections.
//
// -------------------------------------------------------------------------
// Pure virtual methods:
// -------------------------------------------------------------------------
//
// void Read(void)
//
// This method reads the section parameters.
// -------------------------------------------------------------------------

#ifndef _SECTION_H
#define _SECTION_H

#include <map>
#include "material.h"
#include "vec.h"

// -------------------------------------------------------------------------
// Forward declarations:
//
class cElement;
class cMatrix;

// -------------------------------------------------------------------------
// Auxiliary types:
//
typedef struct
{
  double     Thk;    // Thickness
  double     Ang;    // Orientation (in degrees)
  cMaterial *Mat;    // Material
} sLamina;

typedef enum
{
  SEC_DEFAULT,      // Integration chosen by the program
  SEC_FIBER,        // Fiber integration
  SEC_GAUSS         // Gauss integration
} eSecIntMet;

typedef enum
{
 FSDT,  // First-Order Shear Deformation Theory 
 TSDT,  // Third-Order Shear Deformation Theory
 ESDT,  // Exponential Shear Deformation Theory
 SSDT,  // Sinusoidal Shear Deformation Theory
 SHSDT  // Soldatos Hyperbolic Shear Deformation Theory
} eHSDT;

// -------------------------------------------------------------------------
// Section types:
//
typedef enum
{
  SEC_HOM_ISO_2D,       // Isotropic section
  SEC_HOM_ISO_3D,       // Isotropic section for 3D parametric elms
  SEC_HOM_ORTHO_2D,     // Orthotropic section for 2D parametric elms
  SEC_HOM_ORTHO_3D,     // Orthotropic section for 3D parametric elms
  SEC_FGM_2D,           // FGM section for 2D parametric elms
  SEC_FGM_3D,           // FGM section for 3D parametric elms
  SEC_FGM_SHELL,        // FGM section for shells
  SEC_HSDT_PLATE,       // HSDT section for homogeneous materials
  SEC_HSDT_FGM_PLATE,   // HSDT section for FGM materials
  SEC_GENERAL_SHELL,    // General section for shells
  SEC_LAMINATED_2D,     // Laminated section for 2D parametric elms
  SEC_LAMINATED_SHELL,  // Laminated section for shells
  SEC_LAMINATED_SOLID,  // Laminated section for solids
  SEC_INTERFACE_2D,     // Section for 2D interface elements
  SEC_INTERFACE_3D,     // Section for 3D interface elements
  SEC_BAR_GEN,          // General section with given geometric properties
  SEC_BAR_B3DCOUPLED,   // 3D beam section with 4x4 coupled stiffness
  SEC_BAR_HOM_RECT,     // Rectangular cross section (homogeneous)
  SEC_BAR_HOM_CIRCLE,   // Circular cross section (homogeneous)
  SEC_BAR_HOM_TUBE,     // Tubular cross section (homogeneous)
  SEC_BAR_RC_RECT       // Rectangular reinforced concrete cross-section
} eSecType;


// ------------------------------------------------------------------------
// Definition of the Section class:
//
class cSection
{
 private:
  static int        NumSec;    // Number of sections
  static cSection **VecSec;    // Array of sections

 protected:
  static eSecIntMet IntMethod; // Integration method
  static int        NumIntX;   // Number of integ. points (X-axis)
  static int        NumIntY;   // Number of integ. points (Y-axis)
         eSecType   Type;      // Section type
         int        Label;     // Section label

 public:
  static  void       ReadNumSec(void);
  static  void       ReadHomIso2D(void);
  static  void       ReadHomIso3D(void);
  static  void       ReadHomOrtho2D(void);
  static  void       ReadHomOrtho3D(void);
  static  void       ReadFGM2D(void);
  static  void       ReadFGM3D(void);
  static  void       ReadHSDTPlate(void);
  static  void       ReadHSDTFGMPlate(void);
  static  void       ReadFGMShell(void);
  static  void       ReadGeneralShell(void);
  static  void       ReadLaminated2D(void);
  static  void       ReadLaminatedShell(void);
  static  void       ReadLaminatedSolid(void);
  static  void       ReadInterface2D(void);
  static  void       ReadInterface3D(void);
  static  void       ReadSecGen(void);
  static  void       ReadSecB3DCoupled(void);
  static  void       ReadSecRect(void);
  static  void       ReadSecCirc(void);
  static  void       ReadSecTube(void);
  static  void       ReadSecRCRect(void);
  static  void       ReadSecInteg(void);
  static  void       Destroy(void);
  static  eSecIntMet IntegrationMethod(void) { return IntMethod; }
  static  int        NumIntPntX(void) { return NumIntX; }
  static  int        NumIntPntY(void) { return NumIntY; }
  static  cSection  *GetSection(int);

                     cSection(int);
  virtual           ~cSection(void);
          eSecType   GetType(void) { return Type; }
          int        GetLabel(void) { return Label; }
  virtual bool       PreIntegrated(void) { return false; }
  virtual int        GetNumLam(void) { return 0; }
  virtual int        GetNumSecPnt(void) { return 0; }
  virtual double     GetThickness(void) { return 1.0; }
  virtual double     GetDensity(void) { return 0.0; }
  virtual double     GetOrientationAngle(void) { return 0.0; }
  virtual cMaterial *GetMaterial(void) { return 0; }
  virtual sLamina   *GetLayup(void) { return 0; }
  virtual void       GetOrientationVector(cVector &);
  virtual void       Read(void) = 0;
};


// ------------------------------------------------------------------------
// Definition of the 2D Homogeneous Isotropic Section class:
//
class cSecHomIso2D : public cSection
{
 protected:
  double     Thk;   // Section thickness
  cMaterial *Mat;   // Section material
  
 public:
                     cSecHomIso2D(int);
  virtual           ~cSecHomIso2D(void);
          double     GetThickness(void) { return Thk; }
          double     GetDensity(void) { return Mat->GetDensity( ); }
          cMaterial *GetMaterial(void) { return Mat; }
          void       Read(void);
};


// ------------------------------------------------------------------------
// Definition of the 3D Homogeneous Isotropic Section class:
//
class cSecHomIso3D : public cSection
{
 protected:
  cMaterial *Mat;   // Section material
  
 public:
                     cSecHomIso3D(int);
  virtual           ~cSecHomIso3D(void);
          double     GetDensity(void) { return Mat->GetDensity( ); }
          cMaterial *GetMaterial(void) { return Mat; }
          void       Read(void);
};


// ------------------------------------------------------------------------
// Definition of the 2D Elastic Orthotropic Section class:
//
class cSecHomOrtho2D : public cSection
{
 protected:
  double     Thk;   // Section thickness
  double     Ang;   // Orientation angle (radians)
  cMaterial *Mat;   // Section material
  
 public:
                     cSecHomOrtho2D(int);
  virtual           ~cSecHomOrtho2D(void);
          double     GetThickness(void) { return Thk; }
          double     GetOrientationAngle(void) { return Ang; }
          double     GetDensity(void) { return Mat->GetDensity( ); }
          cMaterial *GetMaterial(void) { return Mat; }
          void       GetOrientationVector(cVector &);
          void       Read(void);
};


// ------------------------------------------------------------------------
// Definition of the 3D Elastic Orthotropic Section class:
//
class cSecHomOrtho3D : public cSection
{
 protected:
  cMaterial *Mat;     // Section material
  cVector    OriVec;  // Orientation vector (material x-axis)
  
 public:
                     cSecHomOrtho3D(int);
  virtual           ~cSecHomOrtho3D(void);
          double     GetDensity(void) { return Mat->GetDensity( ); }
          cMaterial *GetMaterial(void) { return Mat; }
          void       GetOrientationVector(cVector &v) { v = OriVec; }
          void       Read(void);
};


// ------------------------------------------------------------------------
// Definition of the FGM 2D Section class:
//
typedef std::map <int, double> VFNodeMap;
typedef std::pair<int, double> VFNodePair;
class cSecFGM2D : public cSection
{
 friend class cSecAnFGM2D;

 protected:
  double     Thk;      // Section thickness
  cMaterial *Mat;      // Section material
  int        NumNode;  // Number of nodes in the FGM section
  VFNodeMap  SecVF;    // Map (node, v2) containing the section nodes

 public:
                     cSecFGM2D(int);
  virtual           ~cSecFGM2D(void);
          double     GetThickness(void) { return Thk; }
          double     GetDensity(void) { return Mat->GetDensity( ); }
          cMaterial *GetMaterial(void) { return Mat; }
          void       Read(void);
};


// ------------------------------------------------------------------------
// Definition of the FGM 3D Section class:
//
class cSecFGM3D : public cSection
{
 friend class cSecAnFGM3D;

 protected:
  cMaterial *Mat;      // Section material
  int        NumNode;  // Number of nodes in the FGM section
  VFNodeMap  SecVF;    // Map (node, v2) containing the section nodes

 public:
                     cSecFGM3D(int);
  virtual           ~cSecFGM3D(void);
          double     GetDensity(void) { return Mat->GetDensity( ); }
          cMaterial *GetMaterial(void) { return Mat; }
          void       Read(void);
};


// ------------------------------------------------------------------------
// Definition of the HSDT Plate Section class:
//
class cSecHSDTPlate : public cSection
{
 friend class cSecAnHSDTPlate;

 protected:
  eHSDT      Theory;   // Plate theory
  double     Thk;      // Section thickness
  cMaterial *Mat;      // Section material

 public:
                     cSecHSDTPlate(int);
  virtual           ~cSecHSDTPlate(void);
          double     GetThickness(void) { return Thk; }
          double     GetDensity(void) { return Mat->GetDensity( ); }
          cMaterial *GetMaterial(void) { return Mat; }
          void       Read(void);
};


// ------------------------------------------------------------------------
// Definition of the HSDT FGM Plate Section class:
//
class cSecHSDTFGMPlate : public cSection
{
 friend class cSecAnHSDTFGMPlate;

 protected:
  eHSDT      Theory;   // Plate theory
  int        NumIpt;   // Integration points
  double     ExpN;     // Volume fraction exponent
  double     Thk;      // Section thickness
  cMaterial *Mat;      // Section material

 public:
                     cSecHSDTFGMPlate(int);
  virtual           ~cSecHSDTFGMPlate(void);
          int        GetNumSecPnt(void) { return NumIpt; }
          double     GetThickness(void) { return Thk; }
          double     GetDensity(void) { return Mat->GetDensity( ); }
          cMaterial *GetMaterial(void) { return Mat; }
          void       Read(void);
};


// ------------------------------------------------------------------------
// Definition of the FGM Shell Section class:
//
class cSecFGMShell : public cSection
{
 friend class cSecAnFGMShell;

 protected:
  int        NumIpt;   // Integration points
  double     ExpN;     // Volume fraction exponent
  double     Thk;      // Section thickness
  cMaterial *Mat;      // Section material

 public:
                     cSecFGMShell(int);
  virtual           ~cSecFGMShell(void);
          int        GetNumSecPnt(void) { return NumIpt; }
          double     GetThickness(void) { return Thk; }
          double     GetDensity(void) { return Mat->GetDensity( ); }
          cMaterial *GetMaterial(void) { return Mat; }
          void       Read(void);
};


// ------------------------------------------------------------------------
// Definition of the General Shell Section class:
//
class cSecGeneralShell : public cSection
{
 friend class cSecAnGenShell;

 protected:
  cMatrix C;  // ABDG matrix
  cMatrix Mb; // Section mass matrix

 public:
                     cSecGeneralShell(int);
  virtual           ~cSecGeneralShell(void);
          void       Read(void);
};


// ------------------------------------------------------------------------
// Definition of the Laminated 2D Section class:
//
class cSecLaminated2D : public cSection
{
 protected:
  int      NumLam;  // Number of plies
  double   LamThk;  // Total laminate thickness
  sLamina *Layup;   // Laminate data
  
 public:
                   cSecLaminated2D(int);
  virtual         ~cSecLaminated2D(void);
          int      GetNumLam(void) { return NumLam; }
          double   GetThickness(void) { return LamThk; }
          double   GetDensity(void);
          sLamina *GetLayup(void) { return Layup; }
          void     Read(void);
};


// ------------------------------------------------------------------------
// Definition of the Laminated Shell Section class:
//
class cSecLaminatedShell : public cSecLaminated2D
{
 protected:
  cVector  OriVec;  // Orientation vector (laminate x-axis)
  
 public:
                   cSecLaminatedShell(int);
  virtual         ~cSecLaminatedShell(void);
          void     GetOrientationVector(cVector &v) { v = OriVec; }
          void     Read(void);
};


// ------------------------------------------------------------------------
// Definition of the Laminated Solid Section class:
//
class cSecLaminatedSolid : public cSection
{
 protected:
  int      NumLam;  // Number of plies
  sLamina *Layup;   // Laminate data
  cVector  OriVec;  // Orientation vector (laminate x-axis)

 public:
                   cSecLaminatedSolid(int);
  virtual         ~cSecLaminatedSolid(void);
          int      GetNumLam(void) { return NumLam; }
          double   GetDensity(void);
          sLamina *GetLayup(void) { return Layup; }
          void     GetOrientationVector(cVector &v) { v = OriVec; }
          void     Read(void);
};


// ------------------------------------------------------------------------
// Definition of the Interface 2D Section class:
//
class cSecInterface2D : public cSection
{
 protected:
  double     Thk;   // Section thickness
  cMaterial *Mat;   // Section material

 public:
                     cSecInterface2D(int);
  virtual           ~cSecInterface2D(void);
          double     GetThickness(void) { return Thk; }
          double     GetDensity(void) { return 0.0; }
          cMaterial *GetMaterial(void) { return Mat; }
          void       Read(void);
};


// ------------------------------------------------------------------------
// Definition of the Interface 3D Section class:
//
class cSecInterface3D : public cSection
{
 protected:
  cMaterial *Mat;   // Section material

 public:
                     cSecInterface3D(int);
  virtual           ~cSecInterface3D(void);
          double     GetDensity(void) { return 0.0; }
          cMaterial *GetMaterial(void) { return Mat; }
          void       Read(void);
};

#endif
