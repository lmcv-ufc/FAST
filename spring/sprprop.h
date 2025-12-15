// -------------------------------------------------------------------------
// sprprop.h - file containing the definition of Spring Property class.
// -------------------------------------------------------------------------
//
// The SpringProp class reads and stores the parameters defining the
// force x deformation (displacement) curve of each spring. In addtion it
// provides methods to compute the force and stiffness of the spring from a
// given deformation.
//
// SpringProp
// |-- Linear
// |-- Nonlinear (piecewise linear)
//
// -------------------------------------------------------------------------
// Static public methods:
// -------------------------------------------------------------------------
//
// static void ReadSpringProp(void)
//
// This method reads the spring properties.
// -------------------------------------------------------------------------
//
// static void Destroy(void)
//
// This method destroys each spring and release the allocated memory. It
// shound be called only after the end of the analysis.
// -------------------------------------------------------------------------
//
// static cSpringProp *GetSpringProp(int label)
//
//   label - label of the spring property                              (in)
//
// This method returns the required spring property from the given label.
// -------------------------------------------------------------------------
// Public methods:
// -------------------------------------------------------------------------
//
// int GetLabel(void)
//
// This method returns the spring property label.
// -------------------------------------------------------------------------
// Pure virtual methods:
// -------------------------------------------------------------------------
//
// virtual void Read(void)
//
// This method reads the spring data (e.g stiffness or load-displ curve).
// -------------------------------------------------------------------------
//
// virtual double GetForce(double u)
//
//   u - spring displacement/deformation                               (in)
//
// This method evaluates the spring force due to the given displacement (or
// deformation).
// -------------------------------------------------------------------------
//
// virtual double GetStiff(double u)
//
//   u - spring displacement/deformation                               (in)
//
// This method evaluates the spring tangent stiffness associated with the
// given displacement (or deformation).
// -------------------------------------------------------------------------

#ifndef _SPRPROP_H
#define _SPRPROP_H

// -------------------------------------------------------------------------
// Definition of SpringProp class:
//
class cSpringProp
{
 private:
  static  int           NumSprProp;    // Number of spring properties
  static  cSpringProp** VecSprProp;    // Vector of spring properties

 protected:
          int           Label;         // Spring property label

 public:
  static  void          ReadSpringProp(void);
  static  void          Destroy(void);
  static  cSpringProp*  GetSpringProp(int);
                        cSpringProp(void);
  virtual              ~cSpringProp(void);
          int           GetLabel(void) { return Label; }
  virtual void          Read(void) = 0;
  virtual double        GetForce(double) = 0;
  virtual double        GetStiff(double) = 0;
};


// -------------------------------------------------------------------------
// Definition of SpringPropLinear class:
//
class cSpringPropLinear : public cSpringProp
{
 protected:
  double K;         // Spring stiffness

 public:
          cSpringPropLinear(int);
         ~cSpringPropLinear(void);
  void    Read(void);
  double  GetForce(double);
  double  GetStiff(double) { return K; }
};


// -------------------------------------------------------------------------
// Definition of cSpringPropNonlin class:
//
class cSpringPropNonlin : public cSpringProp
{
 protected:
  int     NumPnt;    // Number of points (displ x force)
  double *Displ;     // Displacement values
  double *Force;     // Force values

 public:
          cSpringPropNonlin(int);
         ~cSpringPropNonlin(void);
  void    Read(void);
  double  GetForce(double);
  double  GetStiff(double);
};


#endif
