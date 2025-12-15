// -------------------------------------------------------------------------
// eig.cpp - implementation of the generalized eigen value problem solver.
// -------------------------------------------------------------------------
//
// Created:      17-Aug-2017     Elias Saraiva Barroso
//
// Modified:     7-May-2021      Elias Saraiva Barroso
//               Implementation of Shift Iteration method.
//
// Modified:     11-May-2021     Elias Saraiva Barroso
//               Implementation of the Subspaces Iteration method.
// -------------------------------------------------------------------------

#include <iostream>
#include <iomanip>
#include <string>
#include <cmath>
#include <cstdlib>
#include <algorithm>

using namespace std;

#include "eig.h"
#include "sysmat.h"
#include "vec.h"
#include "utl.h"
#include "gblvar.h"

// -------------------------------------------------------------------------
// Auxiliary functions:
//

// ============================= operator>> (eEigAlgType) ================================

istream& operator>> (istream &in, eEigAlgType &type)
{
  // Read the algorithm label.
  char label[100];
  if (!Utl::ReadString(in, label))
  {
    cout << "Algorithm label must be enclosed by double quotes (\"Identifier\").\n";
    in.setstate(ios::failbit);
    return in;
  }

  if (string(label) == "InversePower")
    type = INVERSE_POWER;
  else if (string(label) == "ShiftedPower")
    type = SHIFTED_POWER;
  else if (string(label) == "GeneralizedJacobi")
    type = GENERALIZED_JACOBI;
  else if (string(label) == "SubspaceIteration")
    type = SUBSPACE_ITERATION;
  else
  {
    cout << "Unknown eigen algorithm: " << string(label) << "\n";
    in.setstate(ios::failbit);
    return in;
  }

  return in;
}

// ============================= operator<< (eEigAlgType) ================================

ostream& operator<< (ostream &out, const eEigAlgType &type)
{
  if (type == INVERSE_POWER)
    out << "InversePower";
  else if (type == SHIFTED_POWER)
    out << "ShiftedPower";
  else if (type == GENERALIZED_JACOBI)
    out << "GeneralizedJacobi";
  else if (type == SUBSPACE_ITERATION)
    out << "SubspaceIteration";

  return out;
}

// -------------------------------------------------------------------------
// Struct sEigPairs:
//

// ============================= sEigPair =================================

sEigPair :: sEigPair(int size) : lbd(0.0)
{
  EigVec = new cVector(size);
  EigVec->Zero( );
}

// ============================= ~sEigPair ================================

sEigPair :: ~sEigPair( )
{
  delete EigVec;
}

// ============================= operator= =================================

void sEigPair :: operator= (const sEigPair &p)
{
  lbd    = p.lbd;
  EigVec = p.EigVec;
}

// ============================= ComparePair ===============================

bool sEigPair :: ComparePair(sEigPair *p1, sEigPair *p2)
{
  return (p1->lbd < p2->lbd);
}

// -------------------------------------------------------------------------
// Class cGenEigenProb:
//

// -------------------------------------------------------------------------
// Private methods:
//

// ============================= AllocMemory ===============================

void cGenEigenProb :: AllocMemory(const int &size)
{
  if (NumEigPair > 0 && size > 0) 
  {
    // Release old memory.
    ReleaseMemory( );

    // Alloc memory for the vector of eigen value/vector pairs.
    EigPairVec = new sEigPair* [NumEigPair];
    for(int i = 0; i < NumEigPair; ++i)
      EigPairVec[i] = new sEigPair(size);
  }
}

// ============================= ReleaseMemory =============================

void cGenEigenProb :: ReleaseMemory(void)
{
  if (EigPairVec)
  {
    for(int i = 0; i < NumEigPair; ++i) delete EigPairVec[i];
    delete [] EigPairVec;
  }
}

// -------------------------------------------------------------------------
// Public methods:
//

// ============================= CreateEigenProbAlg ========================

cGenEigenProb* cGenEigenProb :: CreateEigenProbAlg(eEigAlgType type)
{
  // Create the appropriate algorithm.
  if (type == INVERSE_POWER)
    return new cGenEigInvPower( );
  else if (type == SHIFTED_POWER)
    return new cGenEigShfPower( );
  else if (type == GENERALIZED_JACOBI)
    return new cGenEigJacobi( );
  else if (type == SUBSPACE_ITERATION)
    return new cGenEigSubSpaces( );

  return 0;
}

// ============================= cGenEigenProb =============================

cGenEigenProb :: cGenEigenProb( )
{
  EigPairVec = 0;
  NumEigPair = 1;
  CompVec    = true;
}

// ============================= ~cGenEigenProb ============================

cGenEigenProb :: ~cGenEigenProb( )
{
  ReleaseMemory( );
}

// ============================= setCompVec ================================

void cGenEigenProb :: setCompVec(const bool &flag)
{
  CompVec = flag;
}

// ============================= setNumPair ================================

void cGenEigenProb :: setNumPair(const int &np)
{
  if (np > 0)
    NumEigPair = np;
  else
    cout << "Invalid input in cGenEigenProb :: setNumPair, " << np << ".\n";
}

// ============================= getNumPair ================================

int cGenEigenProb :: getNumPair( )
{
  return NumEigPair;
}

// ============================= getPair ===================================

sEigPair* cGenEigenProb :: getPair(const int &id)
{
  if (id >= 0 && id < NumEigPair)
    return EigPairVec[id];
  else
  {
    cout << "Invalid identifier in cGenEigenProb :: GetPair!\n";
    return 0;
  }
}

// -------------------------------------------------------------------------
// Class cGenEigInvPower:
//

// ============================= cGenEigInvPower ===========================

cGenEigInvPower :: cGenEigInvPower( )
{
}

// ============================= ~cGenEigInvPower ==========================

cGenEigInvPower :: ~cGenEigInvPower( )
{
}

// ============================= Solver ====================================

int cGenEigInvPower :: Solver(cSysMatrix *A, cSysMatrix *B, int maxit, double tol)
{
  assert(A->Dim( ) == B->Dim( ) && A->Dim( ) > 0);

  // Auxiliary data.
  int    curriter;
  double dot1,dot2;
  double lbd1,lbd2 = 0.0;
  double err = 0.0;

  // Initialization

  int numeq = A->Dim( );
  cVector v(numeq);
  cVector w(numeq);

  for (int i = 0; i < numeq; i++) v[i] = pow(-1.0,(i % 2))*rand( );
  v /= v.Length( );

  w = B*v;
  lbd1 = 1.0;

  // Set the number of Eigen pairs.
  
  NumEigPair = 1;
  AllocMemory(numeq);

  // Inverse iteration loop

  for (curriter = 1; curriter <= maxit; ++curriter)
  {
    A->Solve(w, v);
    dot1 = v*w;

    w = B*v;

    dot2 = v*w;
    lbd2 = dot1/dot2;
    err  = fabs(lbd2 - lbd1)/(1.0 + fabs(lbd2));

    dot2 = sqrt(fabs(dot2));
    w   /= dot2; // w = dot2/w
    lbd1 = lbd2;

    if (Feedback)
    {
      cout << "Iter: "  << setw(2) << curriter;
      cout << "  Err: " << scientific  << setprecision(4) << err;
      cout << "  LF: "  << fixed << showpos << lbd2 << endl << noshowpos;
    }
    if (err < tol) break;
  }

  // Return computed values

  if (err < tol)
  {
    v /= dot2;

    EigPairVec[0]->lbd       =  lbd2;
    *(EigPairVec[0]->EigVec) =  v;

    return(1);
  }

  return(0);
}

// -------------------------------------------------------------------------
// Class cGenEigShfPower:
//

double cGenEigShfPower :: shf = 0.0;

#include "gblvar.h"

// ============================= ReadShiftFactor ===========================

void cGenEigShfPower :: ReadShiftFactor( )
{
  if (!(in >> shf))
  {
    cout << "Error in the input of the number of steps!\n";
    exit(0);
  }
}

// ============================= cGenEigShfPower ===========================

cGenEigShfPower :: cGenEigShfPower( )
{
}

// ============================= ~cGenEigShfPower ==========================

cGenEigShfPower :: ~cGenEigShfPower( )
{
}

// ============================= Solver ====================================

int cGenEigShfPower :: Solver(cSysMatrix *A, cSysMatrix *B, int maxit, double tol)
{
  assert(A->Dim( ) == B->Dim( ) && A->Dim( ) > 0);

  // Add shift factor to A matrix.
  A->AddMat(-shf,B);
  
  // Process inverse power method.
  int code = this->cGenEigInvPower::Solver(A,B,maxit,tol);

  // Add shift factor to eigen value.
  if (code)
    EigPairVec[0]->lbd += shf;

  return code;
}

// -------------------------------------------------------------------------
// Class cGenEigSubSpaces:
//

// ============================= cGenEigSubSpaces ===========================

cGenEigSubSpaces :: cGenEigSubSpaces( )
{
}

// ============================= ~cGenEigSubSpaces ==========================

cGenEigSubSpaces :: ~cGenEigSubSpaces( )
{
}

// ============================= Solver ====================================

int cGenEigSubSpaces :: Solver(cSysMatrix *A, cSysMatrix *B, int maxit, double tol)
{
  assert(A->Dim( ) == B->Dim( ) && A->Dim( ) > 0);

  bool conv = false;
  cGenEigJacobi jacobi;
  int npair = max(2*NumEigPair,NumEigPair+8);
  int profile[npair];
  for(int i = 0; i < npair; ++i)
  profile[i] = 0;

  cSymSkylMatrix K(npair,profile), M(npair,profile);

  int numeq = B->Dim( );

  AllocMemory(numeq);

  cVector *X = new cVector [npair]; 
  for(int i = 0; i < npair; ++i) X[i].Resize(numeq);    

  // Set initial vectors.
  
  for (int j = 0; j < numeq; j++)  X[0][j] = B->Get(j,j); // X[0] = Dig(B).

  // Set ei vectors. 
  pair<double,int> *inds =  new pair<double,int> [numeq];
 
  for (int i = 0; i < numeq; i++)
  {
    inds[i].first  = A->Get(i,i) / B->Get(i,i);
    inds[i].second = i;
  }
  sort(inds,inds+numeq);

  for(int i = 0; i < npair; ++i)
  {
    X[i].Zero( );
    X[i][inds[i].second] = 1.0;
  }

  // Random vector (last vector).
  for (int j = 0; j < numeq; j++) X[npair-1][j] = pow(-1.0,(j % 2))*rand( );
  X[npair-1] /= X[npair-1].Length( );

  // Setup initial vectors.
  cVector aux(numeq), aux2(numeq), aux3(numeq);
  cMatrix auxM(numeq,npair); 
  
  for (int iter = 0; iter < maxit; ++iter)
  {
    // Solve K Y = M X.
    for(int i = 0; i < npair; ++i)
    {
      if (iter==0) 
        aux = X[i];   // Using starting vectors.
      else 
        aux = B * X[i];

      A->Solve(aux,X[i]);
    }

    // Find projected K.
    K.Zero( );
    for(int i = 0; i < npair; ++i)
    {
      aux = A * X[i];
      for(int j = 0; j < npair; ++j)
        K.Add(j,i,X[j]*aux);
    }

    // Find projected M.
    M.Zero( );
    for(int i = 0; i < npair; ++i)
    {
      aux = B * X[i];
      for(int j = 0; j < npair; ++j)
        M.Add(j,i,X[j]*aux);
    }
   
    // Evaluate eigen value solution in {npair,npair} space.
    bool fb = Feedback;
    Feedback = false;
    int code = jacobi.Solver(&K,&M,maxit,tol);
    if (!code) return 0;
    Feedback = fb;

    // Evaluate updated solution.
    auxM.Zero( );
    for(int i = 0; i < numeq; ++i)
      for(int j = 0; j < npair; ++j)
        for(int k = 0; k < npair; ++k)
          auxM(i,j) += X[k][i] * (*jacobi.getPair(j)->EigVec)[k];

    for(int i = 0; i < numeq; ++i)
      for(int j = 0; j < npair; ++j)
        X[j][i] = auxM[i][j];

    // Evaluate convergence.
    conv = true;
    double err, maxerr = -10e10;

    for(int i = 0; i < NumEigPair; ++i)
    {
      double nom   = pow(jacobi.getPair(i)->lbd,2.0);
      double denom = pow(jacobi.getPair(i)->EigVec->Length( ),2.0);
      err = sqrt(fabs(1.0-nom/denom));

      maxerr = max(err,maxerr);
      if (err - tol > 0) conv = false;
    }

    if (Feedback)
    {
      cout << "Iter: "  << setw(2) << iter+1;
      cout << "  Err: " << scientific  << setprecision(4) << maxerr << endl;
    }

    if (conv)
    {
      for(int i = 0; i < NumEigPair; ++i)
      {
        EigPairVec[i]->lbd        = jacobi.getPair(i)->lbd;
        *(EigPairVec[i]->EigVec)  = X[i];
      }

      break;
    } 
  }

  return conv;
}


// -------------------------------------------------------------------------
// Class cGenEigJacobi:
//

// ============================= cGenEigJacobi =============================

cGenEigJacobi :: cGenEigJacobi( )
{
}

// ============================= ~cGenEigJacobi ============================

cGenEigJacobi :: ~cGenEigJacobi( )
{
}

// ============================= Solver ====================================

int cGenEigJacobi :: Solver(cSysMatrix *A, cSysMatrix *B, int maxit, double tol)
{
  assert(A->Dim( ) == B->Dim( ) && A->Dim( ) > 0);

  // Auxiliary data.
  bool conv = false;
  int  n = A->Dim( );
  int i,j,k;
  double a_ii, a_jj, a_ij, b_ii, b_jj, b_ij, c_ii, c_jj, c, x;
  double gamma, alpha;
  double aux;
  double err;

  // Set the number of Eigen pairs.
  NumEigPair = n;
  AllocMemory(n);

  for(i = 0; i < n; ++i)
    (*(EigPairVec[i]->EigVec))(i) = 1.0;

  cVector currlvec(n), lvec(n);
  for(i = 0; i < n; ++i)
    currlvec(i) = A->Get(i,i);

  for(int curriter = 1; curriter < maxit; ++curriter)
  { 
    // Update vector of eigen values.
    lvec = currlvec;

    // Process sweep.
    for(i = 0; i < n-1; ++i)
      for(j = i+1; j < n; ++j)
      {
        //if (fabs(A->Get(i,j)) < tol && fabs(B->Get(i,j)) < tol)
	//  continue;

        // Get matrices components.
        a_ii = A->Get(i,i);   b_ii = B->Get(i,i);
        a_jj = A->Get(j,j);   b_jj = B->Get(j,j);
        a_ij = A->Get(i,j);   b_ij = B->Get(i,j);


	// Go to the next element
        if (sqrt((a_ij*a_ij)/(a_ii * a_jj)) <= tol)
          if (sqrt((b_ij*b_ij) /(b_ii * b_jj)) <= tol)
	    continue;

	c_ii = a_ii * b_ij - b_ii * a_ij;
	c_jj = a_jj * b_ij - b_jj * a_ij;
	c    = a_ii * b_jj - b_ii * a_jj;

	//if (((fabs(a_ii/b_ii - a_jj/b_jj) < tol) && (fabs(a_ii/b_ii - a_ij/b_ij) < tol)) || (fabs(c)-tol < 0))
	//if (fabs(c) - tol < 0)
	if ((fabs(a_ii/b_ii - a_jj/b_jj) < tol) && (fabs(a_ii/b_ii - a_ij/b_ij) < tol))
	{
	  alpha = 0.0;
	  gamma = -a_ij/a_jj;
	}
	else
	{
	  if (c > 0.0)
  	    x = c / 2.0 + sqrt(pow(c/2.0,2.0) + c_ii*c_jj);
	  else
	    x = c / 2.0 - sqrt(pow(c/2.0,2.0) + c_ii*c_jj);
	  
	  // Evaluate components of equation (11.89).
	  gamma = - c_ii / x;
	  alpha =   c_jj / x;
        }
 
        // =================================== //
        //          Update matrices            // 
        // =================================== //

        // Evaluate elements (i,i), (j,j) and (i,j), in matrices A and B.
	A->Add(i,i,gamma * (2.0 * a_ij + gamma * a_jj));
	A->Add(j,j,alpha * (2.0 * a_ij + alpha * a_ii));
	A->Add(j,i,alpha * (a_ii + gamma * a_ij) + gamma * a_jj);

	B->Add(i,i,gamma * (2.0 * b_ij + gamma * b_jj));
	B->Add(j,j,alpha * (2.0 * b_ij + alpha * b_ii));
	B->Add(j,i,alpha * (b_ii + gamma * b_ij) + gamma * b_jj);
	
        // Combine rows A(i,:) and A(j,:).
	// Index k from 0 to i-1:
        for(k = 0; k < i; ++k) // 
        {
          aux    = A->Get(i,k) * alpha; 
          A->Add(i,k,A->Get(j,k) * gamma); 
          A->Add(j,k,aux);

          aux    = B->Get(i,k) * alpha; 
          B->Add(i,k,B->Get(j,k) * gamma); 
          B->Add(j,k,aux);
        }

	// Index k from i+1 to j-1:
        for(k = i+1; k < j; ++k) 
        {
          aux    = A->Get(k,i) * alpha; 
          A->Add(k,i,A->Get(j,k) * gamma); 
          A->Add(j,k,aux);

          aux    = B->Get(k,i) * alpha; 
          B->Add(k,i,B->Get(j,k) * gamma); 
          B->Add(j,k,aux);
        }

	// Index k from j+1 to n-1:
        for(k = j+1; k < n; ++k) 
        {
          aux    = A->Get(k,i) * alpha; 
          A->Add(k,i,A->Get(k,j) * gamma); 
          A->Add(k,j,aux);

          aux    = B->Get(k,i) * alpha; 
          B->Add(k,i,B->Get(k,j) * gamma); 
          B->Add(k,j,aux);
        }

	// Update Eigen Vectors.
	currlvec = *(EigPairVec[i]->EigVec);
        *(EigPairVec[i]->EigVec) += gamma * *(EigPairVec[j]->EigVec);
        *(EigPairVec[j]->EigVec) += alpha * currlvec;
      }

    // Evaluate convergence criteria.
    conv = true;

    cVector AuxVec(n);

    for(i = 0; i < n; ++i)
      currlvec(i) = A->Get(i,i)/B->Get(i,i);

    for(i = 0; i < n; ++i)
    {
      AuxVec(i) = (lvec(i) - currlvec(i))/currlvec(i);

      if (fabs((lvec(i) - currlvec(i))/currlvec(i)) - tol > 0.0)         
	conv = false;
    }

    err = AuxVec.Length( );

    // Test off-diagonal elements.
    if (conv)
      for(i = 1; i < n; ++i)
	if (conv)
          for(j = 1; j < i; ++j)
          {
            // Get matrices components.
            a_ii = A->Get(i,i);   b_ii = B->Get(i,i);
            a_jj = A->Get(j,j);   b_jj = B->Get(j,j);
            a_ij = A->Get(i,j);   b_ij = B->Get(i,j);
          
            //if (pow(a_ij/(a_ii * a_jj),2.0) - tol > 0.0)
            if (sqrt((a_ij*a_ij)/(a_ii * a_jj)) > tol)
            {
              conv = false;
    //ifff(!conv) cout << "conv2: " << pow(a_ij/(a_ii * a_jj),2.0) << endl;
              break;
            }
            //else if (pow(b_ij/(b_ii * b_jj),2.0) - tol > 0.0)
            else if (sqrt((b_ij*b_ij) /(b_ii * b_jj)) > tol)
            {
              conv = false;
    //ifff(!conv) cout << "conv3: " << pow(b_ij/(b_ii * b_jj),2.0) << endl;
              break;
            }
          }

    if (Feedback)
    {
      cout << "Iter: "  << setw(2) << curriter;
      cout << "  Err: " << scientific  << setprecision(4) << err << endl;
      cout << fixed;
    }
    if (conv && curriter > 1) break;
  } 

  if (conv) 
  {
    // Finish the computation of eigen values and vectors (Eqs. 11.75, 11.76).
    for(i = 0; i < n; ++i)
    {
      EigPairVec[i]->lbd        = currlvec(i);// / B->Get(i,i);
      *(EigPairVec[i]->EigVec) *= 1.0/sqrt(B->Get(i,i));
    }

    // Sort eigen pairs.
    sort(EigPairVec,EigPairVec+n,sEigPair::ComparePair);

    return (1);
  }
  else return (0);
}

