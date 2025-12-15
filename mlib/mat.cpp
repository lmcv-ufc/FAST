// -------------------------------------------------------------------------
// mat.cpp - implementation of the matrix class.
// -------------------------------------------------------------------------
// Created:  14-Jan-1998    Evandro Parente Junior
//
// Modified: 02-Aug-2000    Evandro Parente Junior
//           Use of new/delete instead of calloc/free.
//
// Modified: 21-Sep-2001    Evandro Parente Junior
//           Use of compositors to increase the effic. of some operators.
//
// Modified: 08-Apr-2015    Elias Saraiva Barroso
//           Input/output operations using c++ streams.
//
// Modified: 17-May-2025    Evandro Parente Junior
//           Creation of method Set for submatrices.
// -------------------------------------------------------------------------

#include <math.h>
#include <iostream>
#include <iomanip>
#include <assert.h>

using namespace std;

#include "mat.h"
#include "vec.h"

// -------------------------------------------------------------------------
// Private methods:
//

// ================================ AllocMem ===============================

void cMatrix :: AllocMem(void)
{
  val = new double*[nrow];

  for (int i = 0; i < nrow; i++) val[i] = new double[ncol];
}

// =============================== ReleaseMem ==============================

void cMatrix :: ReleaseMem(void)
{
  for (int i = 0; i < nrow; i++) delete []val[i];

  delete []val;
}

// -------------------------------------------------------------------------
// Public methods:
//

// ============================== cMatrix ==================================

cMatrix :: cMatrix(void)
{
  nrow = 0;
  ncol = 0;
  val  = 0;
}

// ============================== cMatrix ==================================

cMatrix :: cMatrix(int n, int m)
{
  assert(n > 0 && m > 0);
  nrow = n;
  ncol = m;
  AllocMem( );
}

// ============================== cMatrix ==================================

cMatrix :: cMatrix(int n, int m, double **mat)
{
  assert(n > 0 && m > 0);
  nrow = n;
  ncol = m;
  AllocMem( );
  MatAssign(nrow, ncol, mat, val);
}

// ============================== cMatrix ==================================

cMatrix :: cMatrix(const cMatrix &mat)
{
  nrow = mat.nrow;
  ncol = mat.ncol;
  AllocMem( );
  MatAssign(nrow, ncol, mat.val, val);
}

// ============================== cMatrix ==================================

cMatrix :: cMatrix(const sSclMat &op)
{
  nrow = op.mat.nrow;
  ncol = op.mat.ncol;
  AllocMem( );
  MatMult(nrow, ncol, op.scl, op.mat.val, val);
}

// ============================== cMatrix ==================================

cMatrix :: cMatrix(const sAddMat &op)
{
  assert(op.lmat.nrow == op.rmat.nrow && op.lmat.ncol == op.rmat.ncol);
  nrow = op.lmat.nrow;
  ncol = op.lmat.ncol;
  AllocMem( );
  MatAdd(nrow, ncol, op.lmat.val, op.rmat.val, val);
}

// ============================== cMatrix ==================================

cMatrix :: cMatrix(const sSubMat &op)
{
  assert(op.lmat.nrow == op.rmat.nrow && op.lmat.ncol == op.rmat.ncol);
  nrow = op.lmat.nrow;
  ncol = op.lmat.ncol;
  AllocMem( );
  MatSub(nrow, ncol, op.lmat.val, op.rmat.val, val);
}

// ============================== cMatrix ==================================

cMatrix :: cMatrix(const sMulMat &op)
{
  assert(op.lmat.ncol == op.rmat.nrow);
  nrow = op.lmat.nrow;
  ncol = op.rmat.ncol;
  AllocMem( );
  MatMult(nrow, ncol, op.lmat.ncol, op.lmat.val, op.rmat.val, val);
}

// ============================== ~cMatrix =================================

cMatrix :: ~cMatrix(void)
{
  ReleaseMem( );
}

// =============================== Resize ==================================

void cMatrix :: Resize(int n, int m)
{
  ReleaseMem( );
  nrow = n;
  ncol = m;
  AllocMem( );
}

// =============================== Print ===================================

void cMatrix :: Print(void)
{
  for (int i = 0; i < nrow; i++)
  {
    for (int j = 0; j < ncol; j++) cout << val[i][j] << "  ";
    cout << "\n";
  }
}

// ================================ Set ====================================

void cMatrix :: Set(int r, int c, cMatrix &A)
{
  for (int i = 0; i < A.nrow; i++)
  {
    for (int j = 0; j < A.ncol; j++) val[r+i][c+j] = A[i][j];
  }
}

// ================================= DecompLU ==============================

int cMatrix :: DecompLU(void)
{
  return MatDecompLU(nrow, val);
}

// ================================= SolveLU ===============================

int cMatrix :: SolveLU(cVector &b, cVector &x)
{
  x = b;
  MatSolveLU(nrow, val, x.Val( ));
  return(1);
}

// ================================= Solve =================================

int cMatrix :: Solve(cVector &b, cVector &x)
{
  int dec = MatDecompLU(nrow, val);
  if (!dec) return(0);

  x = b;
  MatSolveLU(nrow, val, x.Val( ));
  return(1);
}

// ============================== CompInverse ==============================

int cMatrix :: CompInverse(cMatrix &B)
{
  assert(nrow == ncol && B.nrow == B.ncol && nrow == B.nrow);
  cMatrix A(*this);  // Create a temporary copy

  int dec = MatDecompLU(nrow, A.val);
  if (!dec) return(0);

  // Compute the inverse matrix solving [A][B] = [I].

  double *e = new double[nrow];
  for (int j = 0; j < ncol; j++)
  {
    VecZero(nrow, e);
    e[j] = 1.0;
    MatSolveLU(nrow, A.val, e);
    if (!dec) break;
    for (int i = 0; i < ncol; i++) B.val[i][j] = e[i];
  }
  delete []e;
  return(1);
}

// ============================= operator= =================================

#if 0
cMatrix& cMatrix :: operator=(const cMatrix &mat)
{
  // Check asssignment to self.

  if (this == &mat) return(*this);

  // Handle different dimensions.

  if (nrow != mat.nrow || ncol != mat.ncol)
  {
    ReleaseMem( );
    nrow = mat.nrow;
    ncol = mat.ncol;
    AllocMem( );
  }

  // Copy the given matrix.

  MatAssign(nrow, ncol, mat.val, val);

  // Allow multiple assignment.

  return(*this);
}
#endif


// ================================================ End of file ============
