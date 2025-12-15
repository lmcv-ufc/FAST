
#ifndef _EIGEN_H
#define _EIGEN_H

// -------------------------------------------------------------------------
// Foward declaration:
//
class cVector;
class cSysMatrix;

// -------------------------------------------------------------------------
// Eigen algorithm  types:
//
enum eEigAlgType
{
  INVERSE_POWER,
  SHIFTED_POWER,
  GENERALIZED_JACOBI,
  SUBSPACE_ITERATION
};

istream& operator>> (istream&,eEigAlgType&);
ostream& operator<< (ostream&,const eEigAlgType&);

// -------------------------------------------------------------------------
// Eigen value and vector pair struct:
//
struct sEigPair
{
  double   lbd;
  cVector *EigVec;

  sEigPair(int);
 ~sEigPair(void);

  static bool ComparePair(sEigPair*,sEigPair*);
  void operator= (const sEigPair&);

};

// ------------------------------------------------------------------------
// Definition of abstract generalized eigen value problem solver class:
//
class cGenEigenProb
{
 protected:
          int             NumEigPair;
	  bool            CompVec;     // Option to compute eigen vectors.
	  sEigPair      **EigPairVec;  // Vector of eigen pairs.

	  void   AllocMemory(const int&);
	  void   ReleaseMemory( );

 public:
  static  cGenEigenProb*  CreateEigenProbAlg(eEigAlgType);


                 void setNumPair(const int&);
                 void setCompVec(const bool&);
                 int  getNumPair(void);
                 sEigPair*  getPair(const int&);
                 cGenEigenProb(void);
  virtual       ~cGenEigenProb(void);
  virtual int    Solver(cSysMatrix*,cSysMatrix*,int,double) = 0;
};

// ------------------------------------------------------------------------
// Definition of generalized eigen inverse power iteration solver class:
//
class cGenEigInvPower : public cGenEigenProb
{
 public:
  cGenEigInvPower(void);
  ~cGenEigInvPower(void);
  int  Solver(cSysMatrix*,cSysMatrix*,int,double);
};

// ------------------------------------------------------------------------
// Definition of generalized eigen shifted power iteration solver class:
//
class cGenEigShfPower : public cGenEigInvPower
{
 protected:
  static double shf;       // Shift factor.

 public:
  static void ReadShiftFactor(void);

  cGenEigShfPower(void);
  ~cGenEigShfPower(void);
  int  Solver(cSysMatrix*,cSysMatrix*,int,double);
};

// ------------------------------------------------------------------------
// Definition of Subspaces Iteration Method
//
class cGenEigSubSpaces : public cGenEigenProb
{
 public:
  cGenEigSubSpaces(void);
  ~cGenEigSubSpaces(void);
  int  Solver(cSysMatrix*,cSysMatrix*,int,double);
};

// ------------------------------------------------------------------------
// Definition of generalized eigen Jacobi solver class:
//
class cGenEigJacobi : public cGenEigenProb
{
 public:
  cGenEigJacobi(void);
  ~cGenEigJacobi(void);
  int  Solver(cSysMatrix*,cSysMatrix*,int,double);
};


#endif
