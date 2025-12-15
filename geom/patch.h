// -------------------------------------------------------------------------
// patch.h - file containing the definition of the cPatch2D class.
// -------------------------------------------------------------------------
// Copyright (c) 2015 LMCV/UFC
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without 
// modification, are permitted provided that the following conditions
// are met:
//
// 1. Redistributions of source code must retain the above copyright 
//    notice, this list of conditions and the following disclaimer.
// 2. Redistributions in binary form must reproduce the above copyright 
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS 
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT 
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
// PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT 
// HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, 
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED 
// TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR 
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF 
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING 
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS 
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// -------------------------------------------------------------------------
//
// The cPatch class abstract parametric geometry entity (curves, surfaces and
// solids). The rational propertie is hold in base class and influences
// implementation of parametric location evaluation. Following entities are
// avaliable:
//
// Patch
// |-- BezCurve
// |-- BezTriangle
// |-- BezSurface
// |-- BspPat
// |   |-- BspCurve
// |   |-- BspSurface
// |   |-- BspSolid
//
// -------------------------------------------------------------------------
// Private methods:
// -------------------------------------------------------------------------
//
// void PopBack(void)
//
//   type - given eKnotType type                                        (in)
//   type - given cBSplineBasis *bsp                                    (in)
//
// Remove the last element of the global patch linked-list.
//
// -------------------------------------------------------------------------
// Public methods:
// -------------------------------------------------------------------------
//
// cPatch* CreatePatch(ePatType type)
//
//   type - given ePatType type                                        (in)
//
// This method creates a cPatch object according with given ePatType and return
// its pointer.
// -------------------------------------------------------------------------
//
// void DeletePatch(int label)
//
//   label - given patch label                                         (in)
//
// This method destroys the patch with given label in the global patch list. If
// target patch is not at the end of the list, its position will be swapped with
// the last element and a PopBack( ) operation is performed. Hence, this
// operation may change the position and label of the last element of the list.
// -------------------------------------------------------------------------
//
// void InsertPatch(cPatch *pat)
//
//   pat - given patch                                             (in)
//
// This method insert given patch at the end of global patch list.
// -------------------------------------------------------------------------
//
// void Destroy(void)
//
// This destroys all patches of global list. 
// -------------------------------------------------------------------------
//
// void FindPatch(int i)
//
//   i - given patch label                                         (in)
//
// This method return the patch with given label. A null pointer is returned if
// no patch is found.
// -------------------------------------------------------------------------
//
// void ReadPatch(std::istream &in)
//
//   in - given input stream object                                (in)
//
// This method read a set of patches from given input stream and insert them in
// the global patch list. The list is destroyed before process reading
// operations.
// -------------------------------------------------------------------------
//
// void WritePatch(std::ostream &out)
//
//   out - given output stream object                              (in)
//
// This method writes global patch list in given output stream.
// -------------------------------------------------------------------------
//
// int getNumPatch(void)
//
// This method return the number of patches in global patch list.
// -------------------------------------------------------------------------
//
// void setCtrlPntReadFunc(PatWriteFunc func)
//
// This method register control point read function processed by patches at
// reading phase.
// -------------------------------------------------------------------------
//
// void setCtrlPntWriteFunc(PatWriteFunc func)
//
// This method register control point write function processed by patches at
// writting phase.
// -------------------------------------------------------------------------
//
// int GetNumCP(void) const
//
// This method return the number of control points.
// -------------------------------------------------------------------------
//
// bool isRational(void) const
//
// This method return the Rational flag.
// -------------------------------------------------------------------------
//
// void setRational(bool rat)
//
//   rat - given output stream object                              (in)
//
// This method sets the state of Rational flag.
// -------------------------------------------------------------------------
//
// void setLabel(int id)
//
//   id - given label                                              (in)
//
// This method sets the patch label.
// -------------------------------------------------------------------------
//
// void setRelCtrlPntFlag(bool d)
//
//   d - given flag value.                                             (in)
//
// This method sets the release control point memory flag. (default is true).
// -------------------------------------------------------------------------
//
// bool getLabel(void) const
//
// This method return the patch label.
// -------------------------------------------------------------------------
//
// ePatType getType(void) const
//
// This method return the patch type.
// -------------------------------------------------------------------------
//
// void getNumBas(void) const
//
// This method return total number of basis functions.
// -------------------------------------------------------------------------
//
// sCtrlPnt* getCtrlPnt(const int &id) const
//
//   id - given control point label                                (in)
//
// This method return the idth control point pointer.
// -------------------------------------------------------------------------
//
// void getCtrlPnt(sCtrlPnt* &cpvec, bool alloc) const
//
//   cpvec - Dynamic array of control points                       (out)
//   alloc - memory allocation flag                                (in)
//
// This method return the control points. The alloc flag is used alloc memory
// for dynamic array. Default value is false.
// -------------------------------------------------------------------------
//
// void getCtrlPnt(sCtrlPnt*i* &cpvec, bool alloc) const
//
//   cpvec - Dynamic array of control points pointers              (out)
//   alloc - memory allocation flag                                (in)
//
// This method return the control points references. The alloc flag is used to
// alloc memory for dynamic array. Default value is false.
// -------------------------------------------------------------------------
//
// void setCtrlPnt(const int &id, sCtrlPnt *cp, bool destroy) 
//
//   id      - control point label                                 (in)
//   cp      - given control point pointer                         (in)
//   destroy - flag to destroy old data                            (in)
//
// This method stores given control point pointer in idth position of control
// point vector. The destroy flag is used to delete the old control point
// (default value is true).
// -------------------------------------------------------------------------
//
// void setCtrlPnt(sctrlpnt** cpvec, const int &s, bool destroy, bool inv) 
//
//   cpvec   - dynamic array of control points pointers            (in)
//   s       - size of given dynamic vector                        (in)
//   destroy - flag to destroy old data                            (in)
//   inv     - flag to reverse assign                              (in)
//
// This method stores given array of control points pointers. The destroy flag is used to
// delete the old data (default value is true). The inv flag is used to assign
// control points in reverse order (default value is false).
// -------------------------------------------------------------------------
//
// void setCtrlPnt(sctrlpnt* cpvec, const int &s, bool destroy, bool inv) 
//
//   cpvec   - dynamic array of control points pointers            (in)
//   s       - size of given dynamic vector                        (in)
//   destroy - flag to destroy old data                            (in)
//   inv     - flag to reverse assign                              (in)
//
// This method stores a copy of given array of control points. The destroy flag is used to
// delete the old data (default value is true). The inv flag is used to assign
// control points in reverse order (default value is false).
// -------------------------------------------------------------------------
//
// void getParVarLimits(const int &dim, double &inf, double &sup) const 
//
//   dim - Parametric dimension (e.g. r, s and t)                  (in)
//   inf - inferior limit                                          (in)
//   sup - superior limit                                          (in)
//
// This method return inferior and superior limits of given parametric
// dimension.
// -------------------------------------------------------------------------
//
// int getSpanID(int spanR, int spanS, int spanT)
//
//   spanR  - span index in r direction                      (in)
//   spanS  - span index in s direction                      (in)
//   spanT  - span index in t direction                      (in)
//
// This method return the span index based in span id of each parametric
// direction. spanS and spanT have default value zero.
// -------------------------------------------------------------------------
//
// int getSpanID(double r, double s, double t)
//
//   r,s,t - parametric location                                   (in)
//
// This method return the span index of given parametric location. In case of
// location at boundary of knot spans, the span with smaller index is returned.
// -------------------------------------------------------------------------
//
// int getSpanID(double r, double s, double t, int &size, int inds[])
//
//   r,s,t - parametric location                                   (in)
//   size  - number of span indexes found                          (out)
//   inds  - span indexes found                                    (out)
//
// This method return all span indexes of given parametric location. Note that
// overload versions with s = 0 and t = 0 are avaliable.
// -------------------------------------------------------------------------
//
// void getSpanLimits(int sid, doubl ParVar[][2]) const 
//
//   sid    - given span index                                     (in)
//   ParVar - span limits in each parametric dimension             (in)
//
// This method return inferior and superior limits of given span index for each
// parametric dimension.
// -------------------------------------------------------------------------
//
// sCtrlPnt Evaluate(double r) const 
//
//   r - parametric location                                       (in)
//
// This method evaluate given parametric location in curves patches. 
// -------------------------------------------------------------------------
//
// sCtrlPnt Evaluate(double r, double s) const 
//
//   r,s - parametric location                                     (in)
//
// This method evaluate given parametric location in surface patches. 
// -------------------------------------------------------------------------
//
// sCtrlPnt Evaluate(double r, double s, double t) const 
//
//   r,s,t - parametric location                                   (in)
//
// This method evaluate given parametric location in solid patches. 
// -------------------------------------------------------------------------
//
// void Read(std::istream &in)
//
//   in - given input stream object                                (in)
//
// This method read patch data from given input stream object.
// -------------------------------------------------------------------------
//
// void Write(std::ostream &out)
//
//   out - given output stream object                              (out)
//
// This method write patch data in given output stream object.
//
// -------------------------------------------------------------------------
// Friend Methods:
// -------------------------------------------------------------------------
//
// istream& operator>>(istream &in, cPatch &pat)
//
//   in  - given istream object                                    (in/out)
//   pat - given patch object                                      (in)
//
// This method read patch data from given input stream object.
// -------------------------------------------------------------------------
//
// ostream& operator<<(ostream &out, cPatch &pat)
//
//   out  - given ostream object                                   (in/out)
//   pat - given patch object                                      (out)
//
// This method write patch data in given output stream object.
// -------------------------------------------------------------------------


#ifndef _PATCH_H
#define _PATCH_H

// -------------------------------------------------------------------------
// Forward declarations:
//
#include <iosfwd>

class  cPatch;
class  cBSplineBasis;
struct sCtrlPnt;

// -------------------------------------------------------------------------
// Patch types:
//
typedef enum
{
  BEZ_CURVE,
  BEZ_SURFACE,
  BEZ_TRIAN,
  BSPLINE_CURVE,
  BSPLINE_SURFACE,
  BSPLINE_SOLID
} ePatType;

std::istream& operator>> (std::istream&,ePatType&);
std::ostream& operator<< (std::ostream&,ePatType&);

// -------------------------------------------------------------------------
// Control point I/O function pointer definition:
//
typedef std::ostream& (*PatWriteFunc) (std::ostream&,const cPatch&);
typedef std::istream& (*PatReadFunc)  (std::istream&,cPatch&);

// -------------------------------------------------------------------------
// Basic foward linked-list for global patch vector.
//
struct sPatCell
{
  sPatCell *nxt;
  cPatch   *pat;

  sPatCell( ); 
  sPatCell(cPatch*);
  void Turn(sPatCell*);
};

// -------------------------------------------------------------------------
// Definition of cPatch class:
//
class cPatch
{
 private:
  static int        NewPatLabel;            // Correct new patch label;
  static int        NumPatch;               // Number of patches.  
  static sPatCell*  ListPatHead;            // Global list of patches (head).
  static sPatCell*  ListPatTail;            // Global list of patches (tail).

  static void       PopBack( ); 

 protected:
  static PatWriteFunc WriteCtrlPntFunc;
  static PatReadFunc  ReadCtrlPntFunc;

         int          label;                // Patch label.
	 int          NumCP;                // Number of control points.
         int          NumParVar;            // Number of parametric variables.
         ePatType     type;                 // Patch type.
         sCtrlPnt**   CP;                   // Patch control points.
         bool         Rational;             // Rational flag
         bool         RelCtrlPntMem;        // Flag to destroy control 
                                            // points when patch is destroyed.


 public:
  static   cPatch*   CreatePatch(ePatType);
  static   void      DeletePatch(int);
  static   void      InsertPatch(cPatch *pat);
  static   void      Destroy(void);
  static   cPatch*   FindPatch(int i);
  static   void      ReadPatch(std::istream&);
  static   void      WritePatch(std::ostream&);
  static   int       getNumPatch(void) { return NumPatch; }
           
  static   void      setCtrlPntWriteFunc(PatWriteFunc);
  static   void      setCtrlPntReadFunc(PatReadFunc);


		     cPatch(void);
  virtual            ~cPatch(void);

          int        getNumCP(void) const        { return NumCP;     }
          bool       isRational(void) const      { return Rational;  }
          int        getLabel(void) const        { return label;     }
	  int        getNumParVar(void) const    { return NumParVar; }
          ePatType   getType(void) const         { return type;      }
          sCtrlPnt*  getCtrlPnt(const int &id) const;
          void       getCtrlPnt(sCtrlPnt*&,  bool alloc = 0) const;  
          void       getCtrlPnt(sCtrlPnt**&, bool alloc = 0) const; 
          void       getSpanID(double,double,int&,int[]) const;
          void       getSpanID(double,int&,int[]) const;
          void       setCtrlPnt(const int&,sCtrlPnt*, bool d = 1);
          void       setCtrlPnt(sCtrlPnt*,const int&, bool d = 1, bool inv = 0);
          void       setCtrlPnt(sCtrlPnt**, const int&, bool d = 1, bool inv = 0);
          void       setRational(bool rat)       { Rational = rat;   }
          void       setLabel(int id)            { label = id;       }
	  void       setRelCtrlPntFlag(bool d)   { RelCtrlPntMem = d; }

     // Virtual methods.
     virtual int      getSpanID(int spanR, int spanS = 0, int spanT = 0) const;
     virtual int      getSpanID(double r, double s = 0, double t = 0) const;
     virtual void     getSpanID(double,double,double,int&,int[]) const;
     virtual void     getSpanLimits(int,double[][2]) const;
     virtual int      getNumSpan(void) { return 1; }
     virtual void     Read(std::istream&)  { }
     virtual void     Write(std::ostream&) const { }
     virtual void     getParVarLimits(const int&,double&,double&) const;
     virtual sCtrlPnt Evaluate(double) const;
     virtual sCtrlPnt Evaluate(double,double) const;
     virtual sCtrlPnt Evaluate(double,double,double) const;
     virtual int      getNumBas(void) const = 0;
};

// -------------------------------------------------------------------------
// Define I/O c++ operator:
//
std::istream& operator>>(std::istream&,cPatch&);
std::ostream& operator<<(std::ostream&,const cPatch&);

// -------------------------------------------------------------------------
// Definition of cBezCurvPat class:
//
class cBezCurve : public cPatch
{
 protected:
  int Deg;
 
 public:
                          cBezCurve(void);
                          cBezCurve(const int&, sCtrlPnt**);
  virtual                ~cBezCurve(void);

            void setDegree(const int &deg) {Deg = deg;}
            int getDegree(void) { return Deg;}
            int getNumBas(void) { return Deg+1;}
            sCtrlPnt        Evaluate(double) const;
/*    TODO: Passar para lib geom alg
  void          DegreeElevation(const int&,std::vector<sCtrlPnt>&);
  void          DegreeElevation(const int&,int &, sCtrlPnt *&);
*/



  void      Read(std::istream&);
  void      Write(std::ostream&) const;

};

// -------------------------------------------------------------------------
// Definition of Bezier surface class:
//
class cBezSurface : public cPatch
{
 protected:
  int Deg[2];

 public:
                          cBezSurface(void);
                          cBezSurface(int,int,sCtrlPnt**);
  virtual                ~cBezSurface(void);

            void            setDegree(const int &id, const int &deg) {Deg[id] = deg;}
            int             getDegree(const int &id) {return Deg[id];}
            int             getNumBas(void) const;
            sCtrlPnt        Evaluate(double,double) const;

 // void          DegreeElevation(const int&,std::vector<sCtrlPnt>&);
 // void          DegreeElevation(const int&,int &, sCtrlPnt *&);



  void      Read(std::istream&);
  void      Write(std::ostream&) const;
};


// -------------------------------------------------------------------------
// Definition of cBezTrianPat class:
//

class cBezTriangle : public cPatch
{
 protected:
  int Deg;
 
 public:

                          cBezTriangle(void);
  virtual                ~cBezTriangle(void);
           sCtrlPnt      Evaluate(double,double) const;


           int getDegree(void) const {return Deg;}
           int getNumBas(void) const;

/*    TODO: Passar para lib geom alg
  void      DegreeElevation(const int&,int&,sCtrlPnt*&);
  void      DegreeElevation(const int&,std::vector<sCtrlPnt>&);
*/  
  void      setDegree(const int &deg) { Deg = deg; }
  void      getSideCtrlPntID(int,int*) const;
  void      getInnerCtrlPntID(int&,int*) const;

//  void      Read(std::istream&);
//  void      Write(std::ostream&) const;
};


// -------------------------------------------------------------------------
// Definition Bézier extraction matrix structure:
//
struct sExtMatData
{
  int deg;           // Degree.
  int nspan;         // Number of knotspan.
  double ***C;       // Extraction matrices.
};

// -------------------------------------------------------------------------
// Definition of B-Spline patch class:
//
class cBspPat : public cPatch
{
 protected:
  cBSplineBasis   **bsp;      // B-Spline basis. 
  sExtMatData      *ExtMat;   // Bézier Extraction Matrices.

 public:
                   cBspPat(void);
 virtual          ~cBspPat(void);
          int      getCPInd(int,int j = 0, int k = 0) const;
	  void     getBasInd(const int&, const int&,int*) const;

          void     setDegree(const int&,const int&);
          void     setKnots(const int&,const int&,double*);

          void     LoadExtMat(void);
	  void     ReleaseExtMat(void);
	  double** getExtMat(const int&,const int&) const;
          void     getParVarLimits(const int&,double&,double&) const;
          int      getNumSpan(void) const;
          void     getSpanLimits(int,double[][2]) const;
          int      getSpanID(int spanR, int spanS = 0, int spanT = 0) const;
          int      getSpanID(double r, double s = 0, double t = 0) const;
          void     getSpanID(double,double,double,int&,int[]) const;
          int      getNumBas(void) const;
          void     Read(std::istream &);

  const cBSplineBasis& getBasis(const int&) const;
};

// -------------------------------------------------------------------------
// Definition of B-Spline curve class:
//

class cBspCurve : public cBspPat
{
 public:
                          cBspCurve(void);
                          cBspCurve(cBSplineBasis*);
  virtual                ~cBspCurve(void);

  sCtrlPnt  Evaluate(double) const;

/*    TODO: Passar para lib geom alg
  sCtrlPnt  CompDervFD(double,double dr = 0.0001);
  double  Length(double t1, double t2);
  */

 // void    KnotIns(const double &u, const int &t)
  const cBSplineBasis& getBasis(void) const;
  void    setDegree(const int &d)            {cBspPat::setDegree(0,d);}
  void    setKnots(const int &ks,double *kv) {cBspPat::setKnots(0,ks,kv);}
  int     getNumBas(void) const;


/*    TODO: Passar para lib geom alg
  bool    KnotInsertion(double,int,sCtrlPnt*&,int&,double*&,int&);
  bool    KnotRemoval(double,int,sCtrlPnt*&,int&,double*&,int&);
  bool    DegreeElevation(int,int&,double*&,int&,sCtrlPnt **&);

  void CompBezCtrlPnt(int &nc, int &cs, sCtrlPnt*** &cpvec);
*/  

//  void  Read(std::istream &);
  void  Write(std::ostream &) const;
};


// -------------------------------------------------------------------------
// Definition of B-Spline surface class:
//

class cBspSurface : public cBspPat
{
 public:
                          cBspSurface(void);
                          cBspSurface(cBSplineBasis*,cBSplineBasis*);
  virtual                ~cBspSurface(void);


//  virtual void            Read(ifstream&);

   //        sCtrlPnt  Evaluate(double,double);
  //virtual sNodeCoord      Evaluate(double,double,bool deform = false);

};

// -------------------------------------------------------------------------
// Definition of B-Spline solid class:
//

class cBspSolid : public cBspPat
{
 public:
                          cBspSolid(void);
                          cBspSolid(cBSplineBasis*,cBSplineBasis*,cBSplineBasis*);
  virtual                ~cBspSolid(void);
};


#endif
