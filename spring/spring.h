// -------------------------------------------------------------------------
// spring.h - file containing the definition of Spring class.
// -------------------------------------------------------------------------
//
// The Spring class define the behavior of elastic springs in a FE program.
// This class is very close conceptually to the Element class. The spring
// mechanical behavior is described by the associated SpringProp class,
// which is strongly related to the Material and ConstitutiveModel classes.
// The Spring class reads and stores the data (nodes/directions/properties)
// of each spring and provides methods to compute the internal force vector
// and stiffness matrix of each spring.
//
// Spring
// |-- SpringSupport (1-DOF)
// |-- SpringConnection (2-DOF)
//
// -------------------------------------------------------------------------
// Static public methods:
// -------------------------------------------------------------------------
//
// static void ReadSpringSupp(void)
//
// This method reads the elastic supports (1-dof springs).
// -------------------------------------------------------------------------
//
// static void ReadSpringConn(void)
//
// This method reads the elastic (semi-rigid) connections (2-dof springs).
// -------------------------------------------------------------------------
//
// static void Destroy(void)
//
// This method destroys each spring and release the allocated memory. It
// shound be called only after the end of the analysis.
// -------------------------------------------------------------------------
//
// static int GetNumSpring(void)
//
// This method returns the number of springs of the FE model.
// -------------------------------------------------------------------------
//
// static cSpringSupp *GetSpringSupp(int i)
//
//   i - index of the spring support                                   (in)
//
// This method returns the required elastic support from the given index.
// -------------------------------------------------------------------------
//
// static cSpringConn *GetSpringConn(int i)
//
//   i - index of the spring connection                                (in)
//
// This method returns the required elastic connection from the given index.
//
// -------------------------------------------------------------------------
// Public methods:
// -------------------------------------------------------------------------
//
// int GetLabel(void)
//
// This method returns the spring label.
// -------------------------------------------------------------------------
//
// void AddGlobVec(cVector &sprvec, cVector &glbvec)
//
//   sprvec - spring vector                                            (in)
//   glbvec - global vector                                        (in/out)
//
// This method adds the given spring vector at the appropriate places of
// the given global vector. It is used in the assembly of the global
// internal force vector.
// -------------------------------------------------------------------------
//
// void AddGlobMat(cMatrix &sprmat, cSysMatrix &glbmat)
//
//   sprmat - spring matrix                                            (in)
//   glbmat - global matrix                                        (in/out)
//
// This method adds the given spring matrix at the appropriate places of
// the given global matrix. It is used in the assembly of the global
// stiffness matrix.
//
// -------------------------------------------------------------------------
// Pure virtual methods:
// -------------------------------------------------------------------------
//
// int GetNumDofs(void)
//
// This method returns the number of spring dofs.
// -------------------------------------------------------------------------
//
// void GetDofs(int *dof)
//
//   dof - vector of spring dofs                                      (out)
//
// This method returns the vector of nodal dofs.
// -------------------------------------------------------------------------
//
// virtual void IntForce(cVector &g)
//
//   g - spring internal force vector                                 (out)
//
// This method computes the internal force vector for the current nodal
// displacements.
// -------------------------------------------------------------------------
//
// virtual void StiffMat(cMatrix &K)
//
//   K - spring tangent stiffness matrix                              (out)
//
// This method evaluates the tangent stiffness matrix for the current nodal
// displacements.
// -------------------------------------------------------------------------
//

#ifndef _SPRING_H
#define _SPRING_H

// -------------------------------------------------------------------------
// Forward declarations:
//
class cSpringProp;
class cSpringSupp;
class cSpringConn;
class cNode;
class cVector;
class cMatrix;
class cSysMatrix;

// -------------------------------------------------------------------------
// Definition of Spring class:
//
class cSpring
{
 protected:
          int          Label;        // Spring label
          int          Dir;          // Spring direction (1 to 6)
          cSpringProp* Prop;         // Spring property

 public:
  static  int          GetMaxDof(void) { return 2; }
                       cSpring(void);
  virtual             ~cSpring(void);
          int          GetLabel(void) { return Label; }
          void         AddGlobVec(cVector &, cVector &);
          void         AddGlobMat(cMatrix &, cSysMatrix *);
  virtual int          GetNumDofs(void) = 0;
  virtual void         GetDofs(int *) = 0;
  virtual void         IntForce(cVector &) = 0;
  virtual void         StiffMat(cMatrix &) = 0;
};


// -------------------------------------------------------------------------
// Definition of SpringSupp class:
//
class cSpringSupp : public cSpring
{
 private:
  static  int           NumSpring;    // Number of springs
  static  cSpringSupp** VecSpring;    // Vector of springs

 protected:
          cNode*        Node;         // Spring node

 public:
  static  void          ReadSpringSupp(void);
  static  void          Destroy(void);
  static  int           GetNumSpring(void) { return NumSpring; }
  static  cSpringSupp*  GetSpring(int i)   { return VecSpring[i]; }
                        cSpringSupp(int, cNode *, int, cSpringProp *);
                       ~cSpringSupp(void);
          int           GetNumDofs(void)   { return 1; }
          void          GetDofs(int *);
          void          IntForce(cVector &);
          void          StiffMat(cMatrix &);
};


// -------------------------------------------------------------------------
// Definition of SpringConn class:
//
class cSpringConn : public cSpring
{
 private:
  static  int           NumSpring;    // Number of springs
  static  cSpringConn** VecSpring;    // Vector of springs

 protected:
          cNode*        Node1;         // Spring node (vetor)
          cNode*        Node2;
 public:
  static  void          ReadSpringConn(void);
  static  void          Destroy(void);
  static  int           GetNumSpring(void) { return NumSpring; }
  static  cSpringConn*  GetSpring(int i)   { return VecSpring[i]; }
                        cSpringConn(int, cNode *, cNode *, int, cSpringProp *);
                       ~cSpringConn(void);
          int           GetNumDofs(void)   { return 2; }
          void          GetDofs(int *);
          void          IntForce(cVector &);
          void          StiffMat(cMatrix &);
};

#endif
