// -------------------------------------------------------------------------
// prescdof.h - file containing the definition of the cPrescDof class.
// -------------------------------------------------------------------------
//
// The cPrescDof class reads and stores the parameters defining the
// prescribed degrees of freedom (dofs) of the finite element model. These 
// parameters can be queried by other classes when they are needed.
//
// Each prescribed dof value can vary in time (t) using a different 
// time function h(t):
//
//   dofval(t) = val*h(t)
//
// where val corresponds to the value read from the input file.
//
// Since the neutral file description does not include time functions, it 
// was adopted the concept of a current time function which, after read 
// from the input file, is associated with all the next loads and prescribed
// displacements. In order to handle static cases a constant unit function 
// is defined as the default time function.
//
// -------------------------------------------------------------------------
// Static public methods:
// -------------------------------------------------------------------------
//
// static void ReadPrescDofs(void)
//
// This method creates and reads the prescribed dofs.
// -------------------------------------------------------------------------
//
// static void Destroy(void)
//
// This method destroys all stored dofs.
// -------------------------------------------------------------------------
//
// static int GetNumPrescDofs(void)
//
// This method returns the number of prescribed dofs.
// -------------------------------------------------------------------------
//
// static void GetPrescDofs(int *prdof)
//
//  prdof  -  prescribed dofs                                        (out)
//
// This method returns an array containing the prescribed dofs.
// -------------------------------------------------------------------------
//
// static void GetPrescVals(double t, cVector &prval)
//
//  t      -  time                                                    (in)
//  prval  -  value of the prescribed dof                            (out)
//
// This method returns the values of prescribed dof at the given time.
//
// -------------------------------------------------------------------------
//
// static void CurrNodalVals(cVector &crval)
//
//  crval  -  current values stored at nodes                         (out)
//
// This method returns the current displacements stored at nodes with
// prescribed dofs. These values may be different from the prescribed values
// during the iterative processo, but equal to them after convergence.
//
// -------------------------------------------------------------------------
// Public methods:
// -------------------------------------------------------------------------
//
// double Evaluate(double t)
//
//  t      -  time                                                    (in)
//
// This method returns the prescribed dof value at the given time.
// -------------------------------------------------------------------------
//

#ifndef _PRESCDOF_H
#define _PRESCDOF_H

// ------------------------------------------------------------------------
// Forward declarations:
//
class cNode;
class cVector;
class cTimeFunc;

// -------------------------------------------------------------------------
// Definition of the Prescribed Dof class:
//
class cPrescDof
{
 private:
  static cPrescDof *Head;         // List head
  static cPrescDof *Tail;         // List tail
  static int         NumPrescDof; // List size

 protected:
         int        Dir;          // Dof. direction (1 = u, 2 = v, ...)
         double     DofVal;       // Base value 
         cNode     *Node;         // Node of prescribed dof
         cTimeFunc *Func;         // Time function
         cPrescDof *Next;         // Next dof

 public:
  static void    ReadPrescDofs(void);
  static void    ReadPrescTemps(void);
  static int     GetNumPrescDofs(void) { return NumPrescDof; }
  static void    GetPrescDofs(int *);
  static void    GetPrescVals(double, cVector &);
  static void    CurrNodalVals(cVector &);
  static void    Destroy(void);
                 cPrescDof(cNode *, int, double);
                ~cPrescDof(void);
         double  Evaluate(double);
};

#endif
