// -------------------------------------------------------------------------
// elmplfrm3dcr.cpp - implementation of co-rotational 3D frame element
//                    classes.
// -------------------------------------------------------------------------
// Created:      14-Nov-2012     Luiz Antonio Taumaturgo Mororo
//
// Modified:     28-Oct-2013     Evandro Parente Junior
//               Compatibilization with SecBar class.
// -------------------------------------------------------------------------

#include <math.h>

#include "elmfrm3dcr.h"
#include "elmfrm3d.h"
#include "node.h"
#include "material.h"
#include "secbar.h"
#include "shpline.h"
#include "intpoint.h"
#include "ctrl.h"
#include "load.h"
#include "vec.h"
#include "mat.h"
#include "gblvar.h"
#include "gbldef.h"

// -------------------------------------------------------------------------
// cFrame3DEICR class:
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Public methods:
//

// ============================== cFrame3DEICR =============================

cFrame3DEICR :: cFrame3DEICR(int id, cSection *sec) : cFrame3D(id, sec)
{
  LargeRot3D = 1;
  
  delta_d.Resize(12);
  delta_d.Zero( );
  
  R_1.Resize(3, 3);
  R_1.Zero( );
  R_2.Resize(3, 3);
  R_2.Zero( );
  
  for (int i = 0; i < 3; i++)
  {
    R_1[i][i] = 1;
    R_2[i][i] = 1;
  }
}

// ============================= ~cFrame3DEICR =============================

cFrame3DEICR :: ~cFrame3DEICR(void)
{
}

// =============================== IntForce ================================

int cFrame3DEICR :: IntForce(cVector &g)
{
  // Get the node coordinates.
  
  sNodeCoord coord[2];
  Shape->NodalCoord(coord);
  
  // Get the nodal displacement.
  
  cVector u(12);
  NodalDispl(u);
  
  // Compute the element length.
  
  double L0, L;
  CalcL0(coord, &L0);
  CalcL(coord, u, &L);
  
  // Compute the [T0] matrix.
  
  cMatrix T0(3, 3);
  T0Mat(coord, T0);
  
  // Compute the [Tr] matrix.
  
  cMatrix Tr(3, 3);
  TrMat(coord, u, Tr);
  
  // Compute the corotated nodal displacement
  // => {d}_c = {{0}, {theta}_1_c, (L-L0), 0, 0, {theta}_2_c}.
  
  cVector d_c(12);
  CorNodalDisp(L0, L, T0, Tr, d_c);

  cVector theta_1_c(3);
  cVector theta_2_c(3);
  for (int i = 0; i < 3; i++)
  {
    theta_1_c[i] = d_c[i+3];
    theta_2_c[i] = d_c[i+9];
  }
      
  // Compute the local internal force vector => {f}_c_cls.
  
  cVector f_c_cls(12);
  LocIntForce(f_c_cls);
  
  // Assembly [H]_c matrix.
  
  cMatrix H_c(12, 12);
  HMat(theta_1_c, theta_2_c, H_c);
  
  // Compute the {f}_c vector => {f}_c = [H]t * {f}_c_cls.
  
  cVector f_c(12);
  f_c = t(H_c)*f_c_cls;
  
  // Assembly [P] matrix.
  
  cMatrix PSI(12, 3);
  cMatrix LAMBDA(3, 12);
  cMatrix P(12, 12);
  double eta;
  
  PMat(T0, Tr, L, &eta, PSI, LAMBDA, P);
  
  // Compute the {f}_e vector => {f}_e = [P]t * {f}_c.
  
  cVector f_e(12);
  f_e = t(P)*f_c;
  
  // Assembly [G] matrix of coordinates transformation.
  
  cMatrix G(12, 12);
  GetTrnMat(Tr, G);
  
  // Compute the internal force vector {g}.
  
  g = G*f_e;

  return (1);
}

// =============================== StiffMat ================================

void cFrame3DEICR :: StiffMat(cMatrix &Kt)
{
  // Get the node coordinates.
  
  sNodeCoord coord[2];
  Shape->NodalCoord(coord);
  
  // Get the nodal displacement.
  
  cVector u(12);
  NodalDispl(u);
    
  // Compute the element length.
  
  double L0, L;
  CalcL0(coord, &L0);
  CalcL(coord, u, &L);
  
  // Compute the [T0] matrix.
  
  cMatrix T0(3, 3);
  T0Mat(coord, T0);
    
  // Compute the [Tr] matrix.
  
  cMatrix Tr(3, 3);
  TrMat(coord, u, Tr);
  
  // Compute the corotated nodal displacement
  // => {d}_c = {{0}, {theta}_1_c, (L-L0), 0, 0, {theta}_2_c}.
  
  cVector d_c(12);
  CorNodalDisp(L0, L, T0, Tr, d_c);
  
  cVector theta_1_c(3);
  cVector theta_2_c(3);
  for (int i = 0; i < 3; i++)
  {
    theta_1_c[i] = d_c[i+3];
    theta_2_c[i] = d_c[i+9];
  }
  
  // Get the local stiff. matrix [k]_c_cls.
  
  cMatrix k_c_cls(12, 12);
  LocStiffMat(k_c_cls);
  
  // Compute the {f}_c_cls vector => {f}_c_cls = [k]_c_cls*{d}_c.
  
  cVector f_c_cls(12);
  LocIntForce(f_c_cls);
  
  // Assembly => {n}_1_cls, {m}_1_cls, {n}_2_cls, {m}_2_cls.
  
  cVector n_1_cls(3);
  cVector m_1_cls(3);
  cVector n_2_cls(3);
  cVector m_2_cls(3);
  
  for (int i = 0; i < 3; i++)
  {
    n_1_cls[i] = f_c_cls[i];
    m_1_cls[i] = f_c_cls[i+3];
    n_2_cls[i] = f_c_cls[i+6];
    m_2_cls[i] = f_c_cls[i+9];
  }
  
  // Assembly [H]_c matrix.
  
  cMatrix H_c(12, 12);
  HMat(theta_1_c, theta_2_c, H_c);
  
  cMatrix H_ct(12, 12);    // [H]_c transpose
  for (int i = 0; i < 12; i++)
    for (int j = 0; j < 12; j++)
      H_ct[i][j] = H_c[j][i];
  
  // Assembly [O]_c matrix.
  
  cMatrix O_c(12, 12);
  OMat(theta_1_c, m_1_cls, theta_2_c, m_2_cls, O_c);
  
  // Compute the {f}_c vector => {f]_c = [H]_c t * {f}_c_cls.
  
  cVector f_c(12);
  f_c = t(H_c)*f_c_cls;
  
  // Assembly {m}_1_c and {m}_2_c vectors.
  
  cVector m_1_c(3);
  cVector m_2_c(3);
  for (int i = 0; i < 3; i++)
  {
    m_1_c[i] = f_c[i+3];
    m_2_c[i] = f_c[i+9];
  }
  
  // Compute the [k]_c matrix => [k]_c = [k]_c_1 + [k]_c_2.
  
  cMatrix k_c_1(12, 12);
  cMatrix k_c_1Aux(12, 12);
  k_c_1Aux = H_ct*k_c_cls;
  k_c_1 = k_c_1Aux*H_c;
  
  cMatrix k_c_2(12, 12);
  k_c_2 = O_c*H_c;
  
  cMatrix k_c(12, 12);
  k_c = k_c_1 + k_c_2;
  
  // Assembly eta, [PSI], [LAMBDA] and [P] matrices.
  
  cMatrix PSI(12, 3);
  cMatrix LAMBDA(3, 12);
  cMatrix P(12, 12);
  double eta;
  
  PMat(T0, Tr, L, &eta, PSI, LAMBDA, P);
  
  cMatrix LAMBDAt(12, 3);           // [LAMBDA]t transpose
  for (int i = 0; i < 12; i++)
    for (int j = 0; j < 3; j++)
      LAMBDAt[i][j] = LAMBDA[j][i];
  
  cMatrix Pt(12,12);               // [P]t transpose
  for (int i = 0; i < 12; i++)
    for (int j = 0; j < 12; j++)
      Pt[i][j] = P[j][i];
  
  // Compute the {f}_e vector => {f}_e = [P]t * {f}_c.
  
  cVector f_e(12);
  f_e = Pt*f_c;
  
  // Assembly => {n}_1_e, {m}_1_e, {n}_2_e, {m}_2_e.
  
  cVector n_1_e(3);
  cVector m_1_e(3);
  cVector n_2_e(3);
  cVector m_2_e(3);
  
  for (int i = 0; i < 3; i++)
  {
    n_1_e[i] = f_e[i];
    m_1_e[i] = f_e[i+3];
    n_2_e[i] = f_e[i+6];
    m_2_e[i] = f_e[i+9];
  }
  
  // Assembly [H1] and [H2] matrices.
  
  cMatrix H1(12, 3);
  cMatrix H2(12, 3);
  
  H1Mat(n_1_e, n_2_e, H1);
  H2Mat(n_1_e, m_1_e, n_2_e, m_2_e, H2);
  
  cMatrix H1t(3, 12);             // [H1]t transpose
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 12; j++)
      H1t[i][j] = H1[j][i];
  
  // Assembly [Y] matrix.
  
  cMatrix Y(12, 12);
  YMat(m_1_c, m_2_c, L, eta, Y);
  
  // Compute the [k]_t1_e matrix => [k]_t1_e = [P]t * [k]_c * [P].
  
  cMatrix k_t1_e(12, 12);
  cMatrix k_t1_eAux(12, 12);
  
  k_t1_eAux = Pt*k_c;
  k_t1_e = k_t1_eAux*P;
  
  // Compute the [k]_t2_e matrix => [k]_t2_e = -[LAMBDA]t*[H1]t*[P] - [Y].
  
  cMatrix k_t2_e(12, 12);
  cMatrix k_t2_eAux(12, 12);
  
  k_t2_eAux = -1*LAMBDAt*H1t;
  k_t2_e  = k_t2_eAux*P;
  k_t2_e -= Y;
  
  // Compute the [k]_t3_e matrix => [k]_t3_e = -[H2]*[LAMBDA].
  
  cMatrix k_t3_e(12, 12);
  k_t3_e = -1*H2*LAMBDA;
  
  // Compute the [k]_t_e atrix => [k]_t_e = [k]_t1_e + [k]_t2_e + [k]_t3_e.
  
  cMatrix k_t_e(12, 12);
  k_t_e  = k_t1_e + k_t2_e;
  k_t_e += k_t3_e;
  
  // Get [G] matrix of coordinates transformation.
  
  cMatrix G(12, 12);
  GetTrnMat(Tr, G);
  
  cMatrix Gt(12,12);
  for (int i = 0; i < 12; i++)
    for (int j = 0; j < 12; j++)
      Gt[i][j] = G[j][i];
  
  // Compute [K] matrix => [K] = [G] * [k]_t_e * [G]t.
  
  cMatrix KAux(12, 12);
  KAux.Zero( );
  KAux = G*k_t_e;
  
  Kt.Zero( );
  Kt = KAux*Gt;
}

// ============================== NodalStress ==============================

void cFrame3DEICR :: NodalStress(cMatrix &S)
{
  // Compute the internal force vector.
  
  cVector gl(12);
  LocIntForce(gl);
  
  // Add the fixed-end forces (- equivforces).
  
  cVector fl(12);
  EqvForces(fl);
  gl -= fl;
  
  // Return nodal forces.

  if (OutputConv == DIRECT_STIFFNESS)
  {
    S[0][0] = gl[ 0];  // Initial node
    S[0][1] = gl[ 1];
    S[0][2] = gl[ 2];
    S[0][3] = gl[ 3];
    S[0][4] = gl[ 4];
    S[0][5] = gl[ 5];
    S[1][0] = gl[ 6];  // Final node
    S[1][1] = gl[ 7];
    S[1][2] = gl[ 8];
    S[1][3] = gl[ 9];
    S[1][4] = gl[10];
    S[1][5] = gl[11];
  }
  else 
  {
    S[0][0] = -gl[ 0];  // Initial node
    S[0][1] =  gl[ 1];
    S[0][2] =  gl[ 2];
    S[0][3] = -gl[ 3];
    S[0][4] = -gl[ 4];
    S[0][5] = -gl[ 5];
    S[1][0] =  gl[ 6];  // Final node
    S[1][1] = -gl[ 7];
    S[1][2] = -gl[ 8];
    S[1][3] =  gl[ 9];
    S[1][4] =  gl[10];
    S[1][5] =  gl[11];
  }
}

// ============================== UpdateItera ==============================

void cFrame3DEICR :: UpdateItera(cVector &du)
{
  // Extract local corrector vector from the global one.
  
  GlobToElm(du, delta_d);
  
  // Compute delta_rot_1 and delta_rot_2.
  
  cVector delta_rot_1(3);
  for (int i = 0; i < 3; i++)
    delta_rot_1[i] = delta_d[i+3];
  
  cVector delta_rot_2(3);
  for (int i = 0; i < 3; i++)
    delta_rot_2[i] = delta_d[i+9];
  
  // Update R_1 and R_2.
  
  cMatrix R_1_d(3, 3);
  RMat(delta_rot_1, R_1_d);
  
  cMatrix aux1(3, 3);
  aux1 = R_1_d*R_1;
  R_1  = aux1;
  
  cMatrix R_2_d(3, 3);
  RMat(delta_rot_2, R_2_d);
  
  cMatrix aux2(3, 3);
  aux2 = R_2_d*R_2;
  R_2  = aux2;
}

// -------------------------------------------------------------------------
// Protected methods:
//

// ================================== CalcL0 ===============================
//
// This method computes the length of a given bar in base configuration.
//
//   coord  - vector of nodal coordinates                              (in)
//   L0     - element length                                          (out)
//
void cFrame3DEICR :: CalcL0(sNodeCoord *coord, double *L0)
{
  double dx = coord[1].x - coord[0].x;
  double dy = coord[1].y - coord[0].y;
  double dz = coord[1].z - coord[0].z;
  
  *L0 = sqrt(dx*dx + dy*dy + dz*dz);
}

// ================================== CalcL ================================
//
// This method computes the length of a given bar in current configuration.
//
//   coord  - vector of nodal coordinates                              (in)
//   u      - vector of nodal displacement                             (in)
//   L      - element length                                          (out)
//
void cFrame3DEICR :: CalcL(sNodeCoord *coord, cVector &u, double *L)
{
  // Compute the current position vector of the node 1 => {x1} = {X1} + {u1}.
  
  double x1[3];
  x1[0] = coord[0].x + u[0];
  x1[1] = coord[0].y + u[1];
  x1[2] = coord[0].z + u[2];
  
  // Compute the current position vector of the node 2 => {x2} = {X2} + {u2}.
  
  double x2[3];
  x2[0] = coord[1].x + u[6];
  x2[1] = coord[1].y + u[7];
  x2[2] = coord[1].z + u[8];
  
  // Compute the vector => {x21} = {x2} - {x1}.
  
  double x21[3];
  VecSub(3, x2, x1, x21);     //  computes {x21}
  
  // Compute the length  |{x21}| = L module.
  
  *L = sqrt(x21[0]*x21[0] + x21[1]*x21[1] + x21[2]*x21[2]);
}

// ================================ Spurrier ===============================
//
// This method computes the pseudovector from a rotator matrix by mean
// Spurrier algorithm.
//
//   R      - rotation matrix (3x3)                                    (in)
//   vtheta - pseudovector                                            (out)
//
void cFrame3DEICR :: Spurrier(cMatrix &R, cVector &vtheta)
{
  double trR = R[0][0] + R[1][1] + R[2][2];       // trace of the R matrix
  
  double v[4] = {trR, R[0][0], R[1][1], R[2][2]}; 
  double m = v[0];                                
  double n;                                       
  double vq[3];                                   // quaternion
  double q = 0;                                   // quaternion
  double vw[3];                                   
  double w;                                       // length of vw vector
  
  // Verifying of the max value of the v vector => m = max({v}).
  
  for (int i = 1; i < 4; i++)
    if (m < v[i])
      m = v[i];
  
  // Compute the n, q, {q}.
  
  n = sqrt(1 + 2*m - trR);
  
  if (m == trR)
  {
    q     = n/2;
    vq[0] = (R[2][1] - R[1][2])/(2*n);
    vq[1] = (R[0][2] - R[2][0])/(2*n);
    vq[2] = (R[1][0] - R[0][1])/(2*n);
  }
  
  if (m == R[0][0])
  {
    q     = (R[2][1] - R[1][2])/(2*n);
    vq[0] = (n*n)/(2*n);
    vq[1] = (R[1][0] + R[0][1])/(2*n);
    vq[2] = (R[2][0] + R[0][2])/(2*n);
  }
  
  if (m == R[1][1])
  {
    q     = (R[0][2] - R[2][0])/(2*n);
    vq[0] = (R[0][1] + R[1][0])/(2*n);
    vq[1] = (n*n)/(2*n);
    vq[2] = (R[2][1] + R[1][2])/(2*n);
  }
  
  if (m == R[2][2])
  {
    q     = (R[1][0] - R[0][1])/(2*n);
    vq[0] = (R[0][2] + R[2][0])/(2*n);
    vq[1] = (R[1][2] + R[2][1])/(2*n);
    vq[2] = (n*n)/(2*n);
  }
  
  // Compute => {w} = 2 {q}/q.
  
  for (int i = 0; i < 3; i++)
    vw[i] = 2*vq[i]/q;
  
  // Compute w => w = ||{w}||.
  
  w = sqrt(vw[0]*vw[0] + vw[1]*vw[1] + vw[2]*vw[2]);
  
  // Compute {theta} vector.
  
  if (w == 0) // {theta} = {0}
  {
    for (int i = 0; i < 3; i++)
      vtheta[i] = 0.0;
  }
  else        // {theta} = 2 arctan(w/2) {w}/w
  {
    for (int i = 0; i < 3; i++)
      vtheta[i] = 2*vw[i]/w*atan(w/2);
  }
}

// ================================== T0Mat ================================
//
// This method computes the [T0] matrix. The [T0] elements are the
// direct cosines of the base configuration.
// [T0] = [{e1}_0 {e2}_0 {e3}_0]
//
//   coord  - vector of nodal coordinates                              (in)
//   [T0]   - [T0] matrix (3x3)                                       (out)
//
void cFrame3DEICR :: T0Mat(sNodeCoord *coord, cMatrix &T0)
{
  double e0x[3];  //  vector of the direct cosine of the bar
  double e0y[3];  //  vector of the direct cosine of the y-axis
  double e0z[3];  //  vector of the direct cosine of the z-axis
  
  // Compute the direct cosines.
  
  double dx = coord[1].x - coord[0].x;
  double dy = coord[1].y - coord[0].y;
  double dz = coord[1].z - coord[0].z;
  double L0 = sqrt(dx*dx + dy*dy + dz*dz);
  
  e0x[0] = dx/L0;  //  lx - direct cosine between bar and global x-axis
  e0x[1] = dy/L0;  //  mx - direct cosine between bar and global y-axis
  e0x[2] = dz/L0;  //  nx - direct cosine between bar and global z-axis
  
  VecCrossProd(e0x, VecBeamOrt[BeamOrtIdx].e0y, e0z);
  VecCrossProd(e0z, e0x, e0y);
  
  // Assembly the [T0] matrix.
  
  T0.Zero( );
  T0[0][0] = e0x[0];  T0[0][1] = e0y[0];  T0[0][2] = e0z[0];
  T0[1][0] = e0x[1];  T0[1][1] = e0y[1];  T0[1][2] = e0z[1];
  T0[2][0] = e0x[2];  T0[2][1] = e0y[2];  T0[2][2] = e0z[2];
}

// ================================== TrMat ================================
//
// This method computes the [Tr] matrix. The [Tr] elements are the
// direct cosines of the current configuration.
// [T0] = [{e1} {e2} {e3}]
//
//   coord  - vector of nodal coordinates                              (in)
//   u      - nodal displacements                                      (in)
//   [R1]   - nodal deformational matrix of the node 1                 (in)
//   [T0]   - [T0] matrix (3x3)                                       (out)
//
void cFrame3DEICR :: TrMat(sNodeCoord *coord, cVector &u, cMatrix &Tr)
{
  // Compute the current position vector of the node 1 => {x1} = {X1} + {u1}.
  
  double x1[3];
  x1[0] = coord[0].x + u[0];
  x1[1] = coord[0].y + u[1];
  x1[2] = coord[0].z + u[2];
  
  // Compute the current position vector of the node 2 => {x2} = {X2} + {u2}.
  
  double x2[3];
  x2[0] = coord[1].x + u[6];
  x2[1] = coord[1].y + u[7];
  x2[2] = coord[1].z + u[8];
  
  // Compute the vector => {x21} = {x2} - {x1}.
  
  double x21[3];
  VecSub(3, x2, x1, x21);
  
  //  Compute the |{x21}| = L module.
  
  double L = sqrt(x21[0]*x21[0] + x21[1]*x21[1] + x21[2]*x21[2]);
  
  // Compute the {e1} vector => {e1} = {x21}/L.
  
  double e1[3];
  for (int i = 0; i < 3; i++)
    e1[i] = x21[i]/L;
  
  // Compute the [T0] matrix.
  
  cMatrix T0(3,3);
  T0Mat(coord, T0);
  
  // Compute the {a2} vector => {a2} = [R1][T0]{0 1 0}t.
  
  double a2[3];
  double e02[3] = {0, 1, 0};        //  {e02}  = {0 1 0}t
  cMatrix AuxMat(3,3);              //  AuxMat = [R1][T0]
  AuxMat = R_1*T0;
  for (int i = 0; i < 3; i++)
  {
    a2[i] = 0.0;
    for (int j = 0; j < 3; j++)
    {
      a2[i] += AuxMat[i][j]*e02[j];
    }
  }
  
  // Compute the {e3} vector, {e3} = {e1}x{a2}/(|{e1}x{a2}|).
  
  double e3[3];
  double eaux3[3];
  VecCrossProd(e1, a2, eaux3);  // {eaux3} = {e1}x{a2}
  
  // m_eaux3 = |{eaux3}| = (|{e1}x{a2}|)
  
  double m_eaux3 = sqrt(eaux3[0]*eaux3[0] + eaux3[1]*eaux3[1] + eaux3[2]*eaux3[2]);
  
  for (int i = 0; i < 3; i++)   // {e3}
    e3[i] = eaux3[i]/m_eaux3;
  
  // Compute the {e2} vector, {e2} = {e3}x{e1}.
  
  double e2[3];
  VecCrossProd(e3, e1, e2);
  
  // Assembly [Tr] matrix.
  
  Tr.Zero( );
  Tr[0][0] = e1[0];  Tr[0][1] = e2[0];  Tr[0][2] = e3[0];
  Tr[1][0] = e1[1];  Tr[1][1] = e2[1];  Tr[1][2] = e3[1];
  Tr[2][0] = e1[2];  Tr[2][1] = e2[2];  Tr[2][2] = e3[2];
}

// ================================== SMat =================================
//
// This method computes the skew-symetric matrix S given a vector.
//
//   vec  - vector                                                     (in)
//   S    - Skew-symetric matrix (3x3)                                (out)
//
void cFrame3DEICR :: SMat(cVector &vec, cMatrix &S)
{
  S.Zero( );
  S[0][1] = -vec[2];
  S[0][2] =  vec[1];
  S[1][0] =  vec[2];
  S[1][2] = -vec[0];
  S[2][0] = -vec[1];
  S[2][1] =  vec[0];
}

// ================================== H1Mat ================================
//
// This method computes the auxiliary matrix given {n}_1_e and {n}_2_e
// vectors.
//
//   n_1_e  - vector                                                   (in)
//   n_2_e  - vector                                                   (in)
//   [H1]   - Skew-symetric matrix (12x3)                             (out)
//
//        |[S]_n_1_e|
//        |         |
//        |   [0]   |
// [H1] = |         |
//        |[S]_n_2_e|
//        |         |
//        |   [0]   |
//
void cFrame3DEICR :: H1Mat(cVector &n_1_e, cVector &n_2_e, cMatrix &H1)
{
  H1.Zero( );
  
  H1[0][1] = -n_1_e[2]; H1[0][2] =  n_1_e[1];
  H1[1][0] =  n_1_e[2]; H1[1][2] = -n_1_e[0];
  H1[2][0] = -n_1_e[1]; H1[2][1] =  n_1_e[0];
  
  H1[6][1] = -n_2_e[2]; H1[6][2] =  n_2_e[1];
  H1[7][0] =  n_2_e[2]; H1[7][2] = -n_2_e[0];
  H1[8][0] = -n_2_e[1]; H1[8][1] =  n_2_e[0];
  
}

// ================================== H2Mat ================================
//
// This method comutes the auxiliary matrix given {n}_1_e, {m}_1_e,
// {n}_2_e and {m}_2_e
// vectors.
//
//   n_1_e  - vector                                                   (in)
//   m_1_e  - vector                                                   (in)
//   n_2_e  - vector                                                   (in)
//   m_2_e  - vector                                                   (in)
//   [H2]   - Skew-symetric matrix (12x3)                             (out)
//
//        |[S]_n_1_e|
//        |         |
//        |[S]_m_1_e|
// [H2] = |         |
//        |[S]_n_2_e|
//        |         |
//        |[S]_m_2_e|
//
void cFrame3DEICR :: H2Mat(cVector &n_1_e, cVector &m_1_e, cVector &n_2_e,
                           cVector &m_2_e, cMatrix &H2)
{
  H2.Zero( );
  
  H2[0][1] = -n_1_e[2]; H2[0][2] =  n_1_e[1];
  H2[1][0] =  n_1_e[2]; H2[1][2] = -n_1_e[0];
  H2[2][0] = -n_1_e[1]; H2[2][1] =  n_1_e[0];
  
  H2[3][1] = -m_1_e[2]; H2[3][2] =  m_1_e[1];
  H2[4][0] =  m_1_e[2]; H2[4][2] = -m_1_e[0];
  H2[5][0] = -m_1_e[1]; H2[5][1] =  m_1_e[0];
  
  H2[6][1] = -n_2_e[2]; H2[6][2] =  n_2_e[1];
  H2[7][0] =  n_2_e[2]; H2[7][2] = -n_2_e[0];
  H2[8][0] = -n_2_e[1]; H2[8][1] =  n_2_e[0];
  
  H2[ 9][1] = -m_2_e[2]; H2[ 9][2] =  m_2_e[1];
  H2[10][0] =  m_2_e[2]; H2[10][2] = -m_2_e[0];
  H2[11][0] = -m_2_e[1]; H2[11][1] =  m_2_e[0];
}

// ================================== RMat =================================
//
// This method computes the [R] matrix (tensor) of rotation given a
// pseudovector.
//
//   vec  - pseudovector                                               (in)
//   [R]  - rotation matrix (3x3)                                     (out)
//
void cFrame3DEICR :: RMat(cVector &vec, cMatrix &R)
{
  double m_vec = sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
  
  cMatrix I(3,3);        // Identity matrix 3x3
  I.Zero( );
  I[0][0] = I[1][1] = I[2][2] = 1;
  
  cMatrix S_vec(3,3);    // [S]_vec matrix
  SMat(vec, S_vec);
  
  cMatrix S_vec_2(3,3);  // ([S]_vec)^2 matrix
  S_vec_2 = S_vec*S_vec;
  
  R.Zero( );
  if (m_vec != 0.000)
  {
    R += I;
    R += (sin(m_vec)/m_vec)*S_vec;
    R += (0.5*(sin(m_vec/2)/(m_vec/2)*sin(m_vec/2)/(m_vec/2)))*S_vec_2;
  }
  else
  {
    R += I;
  }
}

// ================================ AlfaMat ================================
//
// This method computes the auxiliary [Alfa] matrix given a vector.
//
//   vec - vector                                                      (in)
//   A   - [Alfa] matrix (3x3)                                        (out)
//
void cFrame3DEICR :: AlfaMat(cVector &vec, cMatrix &A)
{
  double m_vec = sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
  
  cMatrix I(3,3);        // Identity matrix 3x3
  I.Zero( );
  I[0][0] = I[1][1] = I[2][2] = 1;
  
  A.Zero( );
  if (m_vec != 0.000)
  {
    double ksi = (2*sin(m_vec) - m_vec*(1 + cos(m_vec)))/(2*m_vec*m_vec*sin(m_vec));
    
    cMatrix S_vec(3,3);    // [S]_vec matrix
    SMat(vec, S_vec);
    
    cMatrix S_vec_2(3,3);  // ([S]_vec)^2 matrix
    S_vec_2 = S_vec*S_vec;
    
    A += I;
    A -= 0.5*S_vec;
    A += ksi*S_vec_2;
  }
  else
  {
    A += I;
  }
}

// ================================ OmegaMat ===============================
//
// This method computes the auxiliary [Omega] matrix given two vectors.
//
//   vec1, vec2 - vectors                                              (in)
//   Omega      - [Omega] matrix (3x3)                                (out)
//
void cFrame3DEICR :: OmegaMat(cVector &vec1, cVector &vec2, cMatrix &Omega)
{
  double theta = sqrt(vec1[0]*vec1[0] + vec1[1]*vec1[1] + vec1[2]*vec1[2]);
  
  Omega.Zero( );
  if (theta != 0.000)
  {
    cMatrix I(3,3);        // Identity matrix 3x3
    I.Zero( );
    I[0][0] = I[1][1] = I[2][2] = 1;
    
    cMatrix S_vec1(3,3);      // [S]_vec1
    SMat(vec1, S_vec1);
    
    cMatrix S_vec1_2(3,3);    // ([S]_vec1)^2
    S_vec1_2 = S_vec1*S_vec1;
    
    cMatrix S_vec2(3,3);      // [S]_vec2
    SMat(vec2, S_vec2);
    
    double ksi = (2*sin(theta) - theta*(1 + cos(theta)))/(2*theta*theta*sin(theta));
    double mi  = (theta*(sin(theta) + theta) - 8*sin(theta/2)*sin(theta/2))/
    (4*theta*theta*theta*theta*sin(theta/2)*sin(theta/2));
    
    double vec1_tvec2 = 0.0;   // {vec1}t{vec2}
    for (int i = 0; i < 3; i++)
      vec1_tvec2 += vec1[i]*vec2[i];
    
    cMatrix vec1vec2_t(3,3);  // {vec1}{vec2}t
    vec1vec2_t.Zero( );
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
        vec1vec2_t[i][j] += vec1[i]*vec2[j];
    
    cMatrix vec2vec1_t(3,3);  // {vec2}{vec1}t
    vec2vec1_t.Zero( );
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
        vec2vec1_t[i][j] += vec2[i]*vec1[j];
    
    cMatrix Aux(3,3);         // [Aux] = ([S]_vec1)^2
    Aux = S_vec1_2*vec2vec1_t;
    
    // Assembly Omega matrix => [Omega].
    
    // [Omega] = -0.5[S]_vec2 + ksi*({vec1}t{vec2}[I] + {vec1}{vec2}t - 2{vec2}{vec1}t)
    //           + mi*([S]_vec1)^2*{vec2}{vec1}t
    
    Omega -= 0.5*S_vec2;
    Omega += (ksi*vec1_tvec2)*I;
    Omega += ksi*vec1vec2_t;
    Omega -= (2*ksi)*vec2vec1_t;
    Omega += mi*Aux;
  }
  else
  {
    // [Omega] = -0.5[S]_vec2
    
    cMatrix S_vec2(3,3);      // [S]_vec2
    SMat(vec2, S_vec2);
    Omega += -0.5*S_vec2;
  }
}

// ================================== HMat =================================
//
// This method computes the auxiliary [H] matrix given two vectors.
//
//   vec1, vec2 - vectors                                              (in)
//   H          - [H] matrix (12x12)                                  (out)
//
//       | [I]  [0 ]  [0]  [0 ] |
//       | [0]  [A1]  [0]  [0 ] |
// [H] = |                      |
//       | [0]  [0 ]  [I]  [0 ] |
//       | [0]  [0 ]  [0]  [A2] |
//
void cFrame3DEICR :: HMat(cVector &vec1, cVector &vec2, cMatrix &H)
{
  H.Zero( );
  
  // Compute the [H](1,1) = [H](3,3) = [I].
  
  for (int i = 0; i < 3; i++)
    H[i][i] = 1;
  
  for (int i = 6; i < 9; i++)
    H[i][i] = 1;
  
  // Compute the [H](2,2) = [A1], [H](4,4) = [A2].
  
  cMatrix A1(3,3);    // [A1]
  AlfaMat(vec1, A1);
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      H[i+3][j+3] = A1[i][j];
  
  cMatrix A2(3,3);    // [A2]
  AlfaMat(vec2, A2);
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      H[i+9][j+9] = A2[i][j];
}

// ================================== OMat =================================
//
// This method comutes the auxiliary matrix given four vectors.
//
//   vec1, vec2, vec3, vec4  - vector                                  (in)
//   [H2]   - Skew-symetric matrix (12x12)                            (out)
//
//       | [0]  [0  ]  [0]  [0  ] |
//       | [0]  [O_1]  [0]  [0  ] |
// [O] = |                        |
//       | [0]  [0  ]  [0]  [0  ] |
//       | [0]  [0  ]  [0]  [O_2] |
//
void cFrame3DEICR :: OMat(cVector &vec1, cVector &vec2, cVector &vec3, 
                          cVector &vec4, cMatrix &O)
{
  O.Zero( );
  
  cMatrix O_1(3,3);
  OmegaMat(vec1, vec2, O_1);
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      O[i+3][j+3] = O_1[i][j];
  
  cMatrix O_2(3,3);
  OmegaMat(vec3, vec4, O_2);
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      O[i+9][j+9] = O_2[i][j];
}

// ================================== YMat =================================

void cFrame3DEICR :: YMat(cVector &m_1_c, cVector &m_2_c, double L,
                          double eta, cMatrix &Y)
{
  double v_1 = m_1_c[0] + m_2_c[0];
  
  Y.Zero( );
  
  Y[2][5] = -v_1*(1 + eta*eta)/L;
  Y[4][5] =  v_1*(1 + eta*eta);
  Y[8][5] =  v_1*(1 + eta*eta)/L;
}

// ================================== PMat =================================

void cFrame3DEICR :: PMat(cMatrix &T0, cMatrix &Tr, double L, double *eta,
                          cMatrix &PSI,cMatrix &LAMBDA, cMatrix &P)
{
  cMatrix Trt(3,3);                 //  [Tr] transpose
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      Trt[i][j] = Tr[j][i];
  
  cMatrix Aux1(3,3);
  Aux1 = Trt*R_1;
  
  cMatrix Aux2(3,3);
  Aux2 = Aux1*T0;
  
  double e02[3] = {0, 1, 0};        //  {e02}  = {0 1 0}t
  
  cVector a_2_e(3);
  for (int i = 0; i < 3; i++)
  {
    a_2_e[i] = 0.0;
    for (int j = 0; j < 3; j++)
    {
      a_2_e[i] += Aux2[i][j]*e02[j];
    }
  }
  
  double q1 = a_2_e[0];
  double q2 = a_2_e[1];
  
  if (q2 != 0)
    *eta = q1/q2;
  else
    *eta = 0.0;
  
  // Assembly [PSI] matrix.
  
  PSI.Zero( );
  
  PSI[3][0] = PSI[ 4][1] = PSI[ 5][2] = 1.0;
  PSI[9][0] = PSI[10][1] = PSI[11][2] = 1.0;
  PSI[7][2] =  L;
  PSI[8][1] = -L;
  
  // Assembly [LAMBDA] matrix.
  
  LAMBDA.Zero( );
  
  LAMBDA[0][2] =  (*eta)/L;
  LAMBDA[0][3] =  1.0;
  LAMBDA[0][4] = -(*eta);
  LAMBDA[0][8] = -(*eta)/L;
  
  LAMBDA[1][2] =  1.0/L;
  LAMBDA[1][8] = -1.0/L;
  
  LAMBDA[2][1] = -1.0/L;
  LAMBDA[2][7] =  1.0/L;
  
  // Compute the [P] matrix => [P] = [I]_(12x12) - [PSI]*[LAMBDA].
  
  P.Zero( );
  P = PSI*LAMBDA;
  
  for (int i = 0; i < 12; i++)
    for (int j = 0; j < 12; j++)
    {
      if (i == j)
        P[i][j] = 1 - P[i][j];
      else
        P[i][j] = -P[i][j];
    }
}

// =============================== GetTrnMat ===============================
//
// This method computes the transformation matrix between the current and
// the global coordinate system.
//
//   [Tr] - matrix with direction cosines of current configuration     (in)
//   [G]  - transformation matrix (12x12)                             (out)
//
void cFrame3DEICR :: GetTrnMat(cMatrix &Tr, cMatrix &G)
{
  G.Zero( );
  
  G[0][0] = G[3][3] = G[6][6] = G[ 9][ 9] = Tr[0][0];
  G[0][1] = G[3][4] = G[6][7] = G[ 9][10] = Tr[0][1];
  G[0][2] = G[3][5] = G[6][8] = G[ 9][11] = Tr[0][2];
  
  G[1][0] = G[4][3] = G[7][6] = G[10][ 9] = Tr[1][0];
  G[1][1] = G[4][4] = G[7][7] = G[10][10] = Tr[1][1];
  G[1][2] = G[4][5] = G[7][8] = G[10][11] = Tr[1][2];
  
  G[2][0] = G[5][3] = G[8][6] = G[11][ 9] = Tr[2][0];
  G[2][1] = G[5][4] = G[8][7] = G[11][10] = Tr[2][1];
  G[2][2] = G[5][5] = G[8][8] = G[11][11] = Tr[2][2];
}

// ============================= CorNodalDisp ==============================
//
// This method computes the nodal vector of the corotated displacements in
// the local level.
//
//   L0   - element length in base configuration                       (in)
//   L    - element length in base current                             (in)
//   [Tr] - matrix with direction cosines of current configuration     (in)
//   d_c  - nodal vector of the corotated nodal displacements         (out)
//
void cFrame3DEICR :: CorNodalDisp(double L0, double L, cMatrix &T0,
                                  cMatrix &Tr, cVector &d_c)
{
  d_c.Zero( );
  
  // Assembly the transpose of [Tr] => [Tr]t.
  
  cMatrix Trt(3, 3);
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      Trt[i][j] = Tr[j][i];
  
  // Compute the axial deformation.
  
  d_c[6] = L - L0;
  
  // Compute the curvatures of the node 1.
  
  cMatrix R_1_c(3, 3);
  R_1_c = Trt*R_1*T0;
  
  cVector theta_1_c(3);
  Spurrier(R_1_c, theta_1_c);
  for (int i = 0; i < 3; i++)
    d_c[i+3] = theta_1_c[i];
  
  // Compute the curvatures of the node 2.
  
  cMatrix R_2_c(3, 3);
  R_2_c = Trt*R_2*T0;
  
  cVector theta_2_c(3);
  Spurrier(R_2_c, theta_2_c);
  for (int i = 0; i < 3; i++)
    d_c[i+9] = theta_2_c[i];
}


// -------------------------------------------------------------------------
// cFrame3DCR class:
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Public methods:
//

// ============================== cFrame3DCR ===============================

cFrame3DCR :: cFrame3DCR(int id, cSection *sec) : cFrame3DEICR(id, sec)
{
  Type = FRAME3DCR;
}

// ============================== ~cFrame3DCR ==============================

cFrame3DCR :: ~cFrame3DCR(void)
{
}

// ============================= LocIntForce ===============================

void cFrame3DCR :: LocIntForce(cVector &gl)
{
  // Get the nodal coordinates.
  
  sNodeCoord coord[2];
  Shape->NodalCoord(coord);
  
  // Get the nodal displacement.
  
  cVector u(12);
  NodalDispl(u);
  
  // Compute the element length.
  
  double L, L0;
  CalcL0(coord, &L0);
  CalcL(coord, u, &L);
  
  // Compute the [T0] and [Tr] matrix.
  
  cMatrix T0(3, 3);
  T0Mat(coord, T0);
  
  cMatrix Tr(3, 3);
  TrMat(coord, u, Tr);
    
  // Compute the corotated nodal displacement
  // => {d}_c = {{0}, {theta}_1_c, (L-L0), 0, 0, {theta}_2_c}.
  
  cVector d_c(12);
  CorNodalDisp(L0, L, T0, Tr, d_c);
  
  // Get local stiffness matrix [k]_c_cls.
  
  cMatrix k_c_cls(12, 12);
  LocStiffMat(k_c_cls);
  
  gl = k_c_cls*d_c;
}

// ============================= LocStiffMat ===============================

void cFrame3DCR :: LocStiffMat(cMatrix &Kl)
{
  // Get the nodal coordinates.
  
  sNodeCoord coord[2];
  Shape->NodalCoord(coord);
  
  // Compute the element length.
  
  double L0;
  CalcL0(coord, &L0);
  
  // Assembly the local stiffness matrix.
  
  if (Section->GetType( ) == SEC_BAR_B3DCOUPLED)
    CoupLocStiff(L0, Kl);
  else
    ConvLocStiff(L0, Kl);
}

// -------------------------------------------------------------------------
// Protected methods:
//

// ============================== ConvLocStiff =============================
//
// This method computes the stiffness matrix in the local system for conven-
// tional (i.e. diagonal/uncoupled) section matrix.
//
//   L  - element length                                               (in)
//   kl - local stiffness matrix (12x12)                              (out)
//
void cFrame3DCR :: ConvLocStiff(double L, cMatrix &kl)
{
  // Initialization.

  kl.Zero( );
  double L2 = L*L;
  double L3 = L2*L;

  // Material properties.

  cMatMec *mat = Section->GetMaterial( )->GetMecProp( );
  double matpar[2];
  mat->GetParam(matpar);
  double E  = matpar[0];
  double nu = matpar[1];
  double G  = 0.5*E/(1.0 + nu);
  
  // Section stiffness (uncoupled).

  double EA  = E*Section->GetA( );
  double EIz = E*Section->GetIz( );
  double EIy = E*Section->GetIy( );
  double GJ  = G*Section->GetJt( );
  
  // Elastic stiffness matrix.

  kl[ 0][ 0] = kl[ 6][ 6] =  EA/L;
  kl[ 0][ 6] = kl[ 6][ 0] = -EA/L;
  
  kl[ 1][ 1] =  12*EIz/L3;
  kl[ 1][ 5] = kl[ 5][ 1] =   6*EIz/L2;
  kl[ 1][ 7] = kl[ 7][ 1] = -12*EIz/L3;
  kl[ 1][11] = kl[11][ 1] =   6*EIz/L2;
  
  kl[ 2][ 2] =  12*EIy/L3;
  kl[ 2][ 4] = kl[ 4][ 2] =  -6*EIy/L2;
  kl[ 2][ 8] = kl[ 8][ 2] = -12*EIy/L3;
  kl[ 2][10] = kl[10][ 2] =  -6*EIy/L2;
  
  kl[ 3][ 3] = kl[ 9][ 9] =  GJ/L;
  kl[ 3][ 9] = kl[ 9][ 3] = -GJ/L;
  
  kl[ 4][ 4] =  4*EIy/L;
  kl[ 4][ 8] = kl[ 8][ 4] = 6*EIy/L2;
  kl[ 4][10] = kl[10][ 4] = 2*EIy/L;
  
  kl[ 5][ 5] =  4*EIz/L;
  kl[ 5][ 7] = kl[ 7][ 5] = -6*EIz/L2;
  kl[ 5][11] = kl[11][ 5] =  2*EIz/L;
  
  kl[ 7][ 7] =  12*EIz/L3;
  kl[ 7][11] = kl[11][7] = -6*EIz/L2;
  
  kl[ 8][ 8] =  12*EIy/L3;
  kl[ 8][10] = kl[10][8] = 6*EIy/L2;
  
  kl[10][10] =  4*EIy/L;
  
  kl[11][11] =  4*EIz/L;
}

// ============================== CoupLocStiff =============================
//
// This method computes the stiffness matrix in the local system for coupled
// (4x4) section matrix.
//
//   L  - element length                                               (in)
//   kl - local stiffness matrix (12x12)                              (out)
//
void cFrame3DCR :: CoupLocStiff(double L, cMatrix &kl)
{
  // Initialization.

  kl.Zero( );

  // Section stiffness - fully coupled 4x4 matrix.
  
  cMatrix Cb(4, 4);
  Section->SecStiffness3D(Cb);
  double C11 = Cb[0][0]; double C12 = Cb[0][1]; double C13 = Cb[0][2]; double C14 = Cb[0][3];
  double C21 = Cb[1][0]; double C22 = Cb[1][1]; double C23 = Cb[1][2]; double C24 = Cb[1][3];
  double C31 = Cb[2][0]; double C32 = Cb[2][1]; double C33 = Cb[2][2]; double C34 = Cb[2][3];
  double C41 = Cb[3][0]; double C42 = Cb[3][1]; double C43 = Cb[3][2]; double C44 = Cb[3][3];
    
  // Auxiliary variables (Maple code).  

  double t1 = 0.1e1 / L;
  double t2 = t1 * C11;
  double t3 = t1 * C14;
  double t4 = t1 * C13;
  double t5 = t1 * C12;
  double t6 = L * L;
  double t8 = 0.1e1 / t6 / L;
  double t10 = 0.12e2 * t8 * C22;
  double t12 = 0.12e2 * t8 * C23;
  double t13 = 0.1e1 / t6;
  double t15 = 0.6e1 * t13 * C23;
  double t17 = 0.6e1 * t13 * C22;
  double t19 = 0.12e2 * t8 * C32;
  double t21 = 0.12e2 * t8 * C33;
  double t23 = 0.6e1 * t13 * C33;
  double t25 = 0.6e1 * t13 * C32;
  double t26 = t1 * C41;
  double t27 = t1 * C44;
  double t28 = t1 * C43;
  double t29 = t1 * C42;
  double t30 = t1 * C31;
  double t31 = t1 * C34;
  double t32 = t1 * C33;
  double t33 = 0.4e1 * t32;
  double t34 = t1 * C32;
  double t35 = 0.4e1 * t34;
  double t36 = 0.2e1 * t32;
  double t37 = 0.2e1 * t34;
  double t38 = t1 * C21;
  double t39 = t1 * C24;
  double t40 = t1 * C23;
  double t41 = 0.4e1 * t40;
  double t42 = t1 * C22;
  double t43 = 0.4e1 * t42;
  double t44 = 0.2e1 * t40;
  double t45 = 0.2e1 * t42;
  
  // Assembly the [kl] matrix.
  
  kl[0][0] = t2;
  kl[0][1] = 0.0e0;
  kl[0][2] = 0.0e0;
  kl[0][3] = t3;
  kl[0][4] = -t4;
  kl[0][5] = t5;
  kl[0][6] = -t2;
  kl[0][7] = 0.0e0;
  kl[0][8] = 0.0e0;
  kl[0][9] = -t3;
  kl[0][10] = t4;
  kl[0][11] = -t5;
  kl[1][0] = 0.0e0;
  kl[1][1] = t10;
  kl[1][2] = t12;
  kl[1][3] = 0.0e0;
  kl[1][4] = -t15;
  kl[1][5] = t17;
  kl[1][6] = 0.0e0;
  kl[1][7] = -t10;
  kl[1][8] = -t12;
  kl[1][9] = 0.0e0;
  kl[1][10] = -t15;
  kl[1][11] = t17;
  kl[2][0] = 0.0e0;
  kl[2][1] = t19;
  kl[2][2] = t21;
  kl[2][3] = 0.0e0;
  kl[2][4] = -t23;
  kl[2][5] = t25;
  kl[2][6] = 0.0e0;
  kl[2][7] = -t19;
  kl[2][8] = -t21;
  kl[2][9] = 0.0e0;
  kl[2][10] = -t23;
  kl[2][11] = t25;
  kl[3][0] = t26;
  kl[3][1] = 0.0e0;
  kl[3][2] = 0.0e0;
  kl[3][3] = t27;
  kl[3][4] = -t28;
  kl[3][5] = t29;
  kl[3][6] = -t26;
  kl[3][7] = 0.0e0;
  kl[3][8] = 0.0e0;
  kl[3][9] = -t27;
  kl[3][10] = t28;
  kl[3][11] = -t29;
  kl[4][0] = -t30;
  kl[4][1] = -t25;
  kl[4][2] = -t23;
  kl[4][3] = -t31;
  kl[4][4] = t33;
  kl[4][5] = -t35;
  kl[4][6] = t30;
  kl[4][7] = t25;
  kl[4][8] = t23;
  kl[4][9] = t31;
  kl[4][10] = t36;
  kl[4][11] = -t37;
  kl[5][0] = t38;
  kl[5][1] = t17;
  kl[5][2] = t15;
  kl[5][3] = t39;
  kl[5][4] = -t41;
  kl[5][5] = t43;
  kl[5][6] = -t38;
  kl[5][7] = -t17;
  kl[5][8] = -t15;
  kl[5][9] = -t39;
  kl[5][10] = -t44;
  kl[5][11] = t45;
  kl[6][0] = -t2;
  kl[6][1] = 0.0e0;
  kl[6][2] = 0.0e0;
  kl[6][3] = -t3;
  kl[6][4] = t4;
  kl[6][5] = -t5;
  kl[6][6] = t2;
  kl[6][7] = 0.0e0;
  kl[6][8] = 0.0e0;
  kl[6][9] = t3;
  kl[6][10] = -t4;
  kl[6][11] = t5;
  kl[7][0] = 0.0e0;
  kl[7][1] = -t10;
  kl[7][2] = -t12;
  kl[7][3] = 0.0e0;
  kl[7][4] = t15;
  kl[7][5] = -t17;
  kl[7][6] = 0.0e0;
  kl[7][7] = t10;
  kl[7][8] = t12;
  kl[7][9] = 0.0e0;
  kl[7][10] = t15;
  kl[7][11] = -t17;
  kl[8][0] = 0.0e0;
  kl[8][1] = -t19;
  kl[8][2] = -t21;
  kl[8][3] = 0.0e0;
  kl[8][4] = t23;
  kl[8][5] = -t25;
  kl[8][6] = 0.0e0;
  kl[8][7] = t19;
  kl[8][8] = t21;
  kl[8][9] = 0.0e0;
  kl[8][10] = t23;
  kl[8][11] = -t25;
  kl[9][0] = -t26;
  kl[9][1] = 0.0e0;
  kl[9][2] = 0.0e0;
  kl[9][3] = -t27;
  kl[9][4] = t28;
  kl[9][5] = -t29;
  kl[9][6] = t26;
  kl[9][7] = 0.0e0;
  kl[9][8] = 0.0e0;
  kl[9][9] = t27;
  kl[9][10] = -t28;
  kl[9][11] = t29;
  kl[10][0] = t30;
  kl[10][1] = -t25;
  kl[10][2] = -t23;
  kl[10][3] = t31;
  kl[10][4] = t36;
  kl[10][5] = -t37;
  kl[10][6] = -t30;
  kl[10][7] = t25;
  kl[10][8] = t23;
  kl[10][9] = -t31;
  kl[10][10] = t33;
  kl[10][11] = -t35;
  kl[11][0] = -t38;
  kl[11][1] = t17;
  kl[11][2] = t15;
  kl[11][3] = -t39;
  kl[11][4] = -t44;
  kl[11][5] = t45;
  kl[11][6] = t38;
  kl[11][7] = -t17;
  kl[11][8] = -t15;
  kl[11][9] = t39;
  kl[11][10] = -t41;
  kl[11][11] = t43;
}


// -------------------------------------------------------------------------
// cFrame3DCRTL class:
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Public methods:
//

// ============================== cFrame3DCRTL =============================

cFrame3DCRTL :: cFrame3DCRTL(int id, cSection *sec) : cFrame3DEICR(id, sec)
{
  Type = FRAME3DCRTL;
}

// ============================== cFrame3DCRTL =============================

cFrame3DCRTL :: ~cFrame3DCRTL(void)
{
}

// ============================== LocIntForce ==============================

void cFrame3DCRTL :: LocIntForce(cVector &gl)
{
   // Get the nodal coordinates.
  
  sNodeCoord coord[2];
  Shape->NodalCoord(coord);

  // Get the nodal displacement.
  
  cVector u(12);
  NodalDispl(u);

  // Compute the element length.

  double L, L0;
  CalcL0(coord, &L0);
  CalcL(coord, u, &L);

  // Compute the [T0] and [Tr] matrix.

  cMatrix T0(3, 3);
  T0Mat(coord, T0);

  cMatrix Tr(3, 3);
  TrMat(coord, u, Tr);

  // Compute the corotated nodal displacement
  // => {d}_c = {{0}, {theta}_1_c, (L-L0), 0, 0, {theta}_1_c}.

  cVector d_c(12);
  CorNodalDisp(L0, L, T0, Tr, d_c);
  
  // Evaluate the internal force vector.

  if (Section->GetType( ) == SEC_BAR_B3DCOUPLED)
    CoupIntForce(L0, d_c, gl);
  else
    ConvIntForce(L0, d_c, gl);
}

// ============================== LocStiffMat ==============================

void cFrame3DCRTL :: LocStiffMat(cMatrix &Kl)
{
  // Get the nodal coordinates.
  
  sNodeCoord coord[2];
  Shape->NodalCoord(coord);

  // Get the nodal displacement.
  
  cVector u(12);
  NodalDispl(u);

  // Compute the element length.

  double L, L0;
  CalcL0(coord, &L0);
  CalcL(coord, u, &L);

  // Compute the [T0] and [Tr] matrix.

  cMatrix T0(3, 3);
  T0Mat(coord, T0);

  cMatrix Tr(3, 3);
  TrMat(coord, u, Tr);

  // Compute the corotated nodal displacement
  // => {d}_c = {{0}, {theta}_1_c, (L-L0), 0, 0, {theta}_1_c}.

  cVector d_c(12);
  CorNodalDisp(L0, L, T0, Tr, d_c);

  // Compute the [Kl] = int([Bbar]t[C][Bbar]).

  Kl.Zero( );
  if (Section->GetType( ) == SEC_BAR_B3DCOUPLED)
    CoupStiffMat(L0, d_c, Kl);
  else
    ConvStiffMat(L0, d_c, Kl);

  // Compute the normal force.

  double NL;
  if (Section->GetType( ) == SEC_BAR_B3DCOUPLED)
    NL = CoupNormal(L, d_c); // Inclui o L ???
  else
  {
    double N = ConvNormal(L0, d_c);
    NL = N*L;   
  }

  // Compute the "geometric" stiffness matrix => Kg = NL*[A]t.

  double A  = Section->GetA( );
  double Ip = Section->GetJt( ); // ??
  cMatrix At(12, 12);
  AMat(L0, Ip, A, At);

  cMatrix Kg(12, 12);
  Kg = NL*At;

  // Add the geometric stiffness matrix.
  
  Kl += Kg;
}

// -------------------------------------------------------------------------
// Protected methods:
//

// =================================== AMat ================================
//
// This method computes the auxiliary matrix [At] due to transverse
// displacements v and w, and rotation of the cross section effects.
//
//   L  - element length                                               (in)
//   Ip - polar moment of inertia of the cross section                 (in)
//   A  - area of the cross section                                    (in)
//   At - auxiliary matrix (12x12)                                    (out)
//
void cFrame3DCRTL :: AMat(double L, double Ip, double A, cMatrix & At)
{
  // Auxiliary variables.
  
  double t1 = L * L;
  double t2 = 0.1e1 / t1;
  double t3 = 0.6e1 / 0.5e1 * t2;
  double t5 = 0.1e1 / L / 0.10e2;
  double t8 = t2 * Ip / A;
  
  // Assembly the [A] matrix.
  
  At[0][0] = 0.0e0;
  At[0][1] = 0.0e0;
  At[0][2] = 0.0e0;
  At[0][3] = 0.0e0;
  At[0][4] = 0.0e0;
  At[0][5] = 0.0e0;
  At[0][6] = 0.0e0;
  At[0][7] = 0.0e0;
  At[0][8] = 0.0e0;
  At[0][9] = 0.0e0;
  At[0][10] = 0.0e0;
  At[0][11] = 0.0e0;
  At[1][0] = 0.0e0;
  At[1][1] = t3;
  At[1][2] = 0.0e0;
  At[1][3] = 0.0e0;
  At[1][4] = 0.0e0;
  At[1][5] = t5;
  At[1][6] = 0.0e0;
  At[1][7] = -t3;
  At[1][8] = 0.0e0;
  At[1][9] = 0.0e0;
  At[1][10] = 0.0e0;
  At[1][11] = t5;
  At[2][0] = 0.0e0;
  At[2][1] = 0.0e0;
  At[2][2] = t3;
  At[2][3] = 0.0e0;
  At[2][4] = -t5;
  At[2][5] = 0.0e0;
  At[2][6] = 0.0e0;
  At[2][7] = 0.0e0;
  At[2][8] = -t3;
  At[2][9] = 0.0e0;
  At[2][10] = -t5;
  At[2][11] = 0.0e0;
  At[3][0] = 0.0e0;
  At[3][1] = 0.0e0;
  At[3][2] = 0.0e0;
  At[3][3] = t8;
  At[3][4] = 0.0e0;
  At[3][5] = 0.0e0;
  At[3][6] = 0.0e0;
  At[3][7] = 0.0e0;
  At[3][8] = 0.0e0;
  At[3][9] = -t8;
  At[3][10] = 0.0e0;
  At[3][11] = 0.0e0;
  At[4][0] = 0.0e0;
  At[4][1] = 0.0e0;
  At[4][2] = -t5;
  At[4][3] = 0.0e0;
  At[4][4] = 0.2e1 / 0.15e2;
  At[4][5] = 0.0e0;
  At[4][6] = 0.0e0;
  At[4][7] = 0.0e0;
  At[4][8] = t5;
  At[4][9] = 0.0e0;
  At[4][10] = -0.1e1 / 0.30e2;
  At[4][11] = 0.0e0;
  At[5][0] = 0.0e0;
  At[5][1] = t5;
  At[5][2] = 0.0e0;
  At[5][3] = 0.0e0;
  At[5][4] = 0.0e0;
  At[5][5] = 0.2e1 / 0.15e2;
  At[5][6] = 0.0e0;
  At[5][7] = -t5;
  At[5][8] = 0.0e0;
  At[5][9] = 0.0e0;
  At[5][10] = 0.0e0;
  At[5][11] = -0.1e1 / 0.30e2;
  At[6][0] = 0.0e0;
  At[6][1] = 0.0e0;
  At[6][2] = 0.0e0;
  At[6][3] = 0.0e0;
  At[6][4] = 0.0e0;
  At[6][5] = 0.0e0;
  At[6][6] = 0.0e0;
  At[6][7] = 0.0e0;
  At[6][8] = 0.0e0;
  At[6][9] = 0.0e0;
  At[6][10] = 0.0e0;
  At[6][11] = 0.0e0;
  At[7][0] = 0.0e0;
  At[7][1] = -t3;
  At[7][2] = 0.0e0;
  At[7][3] = 0.0e0;
  At[7][4] = 0.0e0;
  At[7][5] = -t5;
  At[7][6] = 0.0e0;
  At[7][7] = t3;
  At[7][8] = 0.0e0;
  At[7][9] = 0.0e0;
  At[7][10] = 0.0e0;
  At[7][11] = -t5;
  At[8][0] = 0.0e0;
  At[8][1] = 0.0e0;
  At[8][2] = -t3;
  At[8][3] = 0.0e0;
  At[8][4] = t5;
  At[8][5] = 0.0e0;
  At[8][6] = 0.0e0;
  At[8][7] = 0.0e0;
  At[8][8] = t3;
  At[8][9] = 0.0e0;
  At[8][10] = t5;
  At[8][11] = 0.0e0;
  At[9][0] = 0.0e0;
  At[9][1] = 0.0e0;
  At[9][2] = 0.0e0;
  At[9][3] = -t8;
  At[9][4] = 0.0e0;
  At[9][5] = 0.0e0;
  At[9][6] = 0.0e0;
  At[9][7] = 0.0e0;
  At[9][8] = 0.0e0;
  At[9][9] = t8;
  At[9][10] = 0.0e0;
  At[9][11] = 0.0e0;
  At[10][0] = 0.0e0;
  At[10][1] = 0.0e0;
  At[10][2] = -t5;
  At[10][3] = 0.0e0;
  At[10][4] = -0.1e1 / 0.30e2;
  At[10][5] = 0.0e0;
  At[10][6] = 0.0e0;
  At[10][7] = 0.0e0;
  At[10][8] = t5;
  At[10][9] = 0.0e0;
  At[10][10] = 0.2e1 / 0.15e2;
  At[10][11] = 0.0e0;
  At[11][0] = 0.0e0;
  At[11][1] = t5;
  At[11][2] = 0.0e0;
  At[11][3] = 0.0e0;
  At[11][4] = 0.0e0;
  At[11][5] = -0.1e1 / 0.30e2;
  At[11][6] = 0.0e0;
  At[11][7] = -t5;
  At[11][8] = 0.0e0;
  At[11][9] = 0.0e0;
  At[11][10] = 0.0e0;
  At[11][11] = 0.2e1 / 0.15e2;
}

// ================================== AvMat ================================
//
// This method computes the auxiliary matrix [Av] due to transverse
// displacement v effects.
//
//   L  - element length                                               (in)
//   Av - auxiliary matrix (12x12)                                    (out)
//
void cFrame3DCRTL :: AvMat(double L, cMatrix &Av)
{
  Av.Zero( );
  double L2 = L*L;
  
  Av[1][1]   =  6/(5*L2);
  Av[1][5]   =  Av[5][1]  =  1/(10*L);
  Av[1][7]   =  Av[7][1]  = -6/(5*L2);
  Av[1][11]  =  Av[11][1] =  1/(10*L);
  
  Av[5][5]   =  0.13333333333;
  Av[5][7]   =  Av[7][5]  = -1/(10*L);
  Av[5][11]  =  Av[11][5] = -0.03333333333;
  
  Av[7][7]   =  6/(5*L2);
  Av[7][11]  =  Av[11][7] = -1/(10*L);
  
  Av[11][11] = 0.13333333333;
}

// ================================== AwMat ================================
//
// This method computes the auxiliary matrix [Aw] due to transverse
// displacement w effects.
//
//   L  - element length                                               (in)
//   Aw - auxiliary matrix (12x12)                                    (out)
//
void cFrame3DCRTL :: AwMat(double L, cMatrix &Aw)
{
  Aw.Zero( );
  double L2 = L*L;
  
  Aw[2][2]   =  6/(5*L2);
  Aw[2][4]   =  Aw[4][2]  = -1/(10*L);
  Aw[2][8]   =  Aw[8][2]  = -6/(5*L2);
  Aw[2][10]  =  Aw[10][2] = -1/(10*L);
  
  Aw[4][4]   =  0.13333333333;
  Aw[4][8]   =  Aw[8][4]  = 1/(10*L);
  Aw[4][10]  =  Aw[10][4] = -0.03333333333;
  
  Aw[8][8]   = 6/(5*L2);
  Aw[8][10]  = Aw[10][8]  = 1/(10*L);
  
  Aw[10][10] = 0.13333333333;
}

// ================================= AthetaMat =============================
//
// This method computes the auxiliary matrix [Atheta] due to twisting
// effects.
//
//   L      - element length                                           (in)
//   Ip     - polar moment of inertia of the cross section             (in)
//   A      - area of the cross section                                (in)
//   Atheta - auxiliary matrix (12x12)                                (out)
//
void cFrame3DCRTL :: AthetaMat(double L, double Ip, double A, cMatrix &Atheta)
{
  Atheta.Zero( );
  double L2 = L*L;
  
  Atheta[3][3] = Atheta[9][9] =  Ip/(L2*A);
  Atheta[3][9] = Atheta[9][3] = -Ip/(L2*A);
}

// ================================== B0Mat ================================
//
// This method computes the element linear strain-displacement matrix.
//
//   L  - element length                                               (in)
//   B0 - linear strain-displacement matrix (1x12)                    (out)
//
void cFrame3DCRTL :: B0Mat(double L, cMatrix &B0)
{
  B0.Zero( );
  B0[0][0] = -1/L;
  B0[0][6] =  1/L;
}

// ================================== BLvMat ===============================
//
// This method computes the [BLv] matrix (nonlinear) related to membrane
// strains. Obtained from transverse displacement v effects.
//
//   L   - element length                                              (in)
//   u   - nodal displacement vector                                   (in)
//   BLv - strain-displacement matrix (1x12)                          (out)
//
void cFrame3DCRTL :: BLvMat(double L, cVector &u, cMatrix &BLv)
{
  // Auxiliary variables.
  
  double t1 = L * L;
  double t2 = 0.1e1 / t1;
  double t5 = 0.1e1 / L;
  double t12 = 0.6e1 / 0.5e1 * u[1] * t2 + u[5] * t5 / 0.10e2 - 0.6e1 / 0.5e1 * u[7] * t2 + u[11] * t5 / 0.10e2;
  double t14 = u[1] * t5 / 0.10e2;
  double t17 = u[7] * t5 / 0.10e2;
  
  // Assembly BLv matrix.
  
  BLv[0][0] = 0.0e0;
  BLv[0][1] = t12;
  BLv[0][2] = 0.0e0;
  BLv[0][3] = 0.0e0;
  BLv[0][4] = 0.0e0;
  BLv[0][5] = t14 + 0.2e1 / 0.15e2 * u[5] - t17 - u[11] / 0.30e2;
  BLv[0][6] = 0.0e0;
  BLv[0][7] = -t12;
  BLv[0][8] = 0.0e0;
  BLv[0][9] = 0.0e0;
  BLv[0][10] = 0.0e0;
  BLv[0][11] = t14 - u[5] / 0.30e2 - t17 + 0.2e1 / 0.15e2 * u[11];
}

// ================================== BLwMat ===============================
//
// This method computes the [BLw] matrix (nonlinear) related to membrane
// strains. Obtained from transverse displacement effects.
//
//   L   - element length                                              (in)
//   u   - nodal displacement vector                                   (in)
//   BLw - strain-displacement matrix (1x12)                          (out)
//
void cFrame3DCRTL :: BLwMat(double L, cVector &u, cMatrix &BLw)
{
  // Auxiliary variables.
  
  double t1 = L * L;
  double t2 = 0.1e1 / t1;
  double t5 = 0.1e1 / L;
  double t12 = 0.6e1 / 0.5e1 * u[2] * t2 - u[4] * t5 / 0.10e2 - 0.6e1 / 0.5e1 * u[8] * t2 - u[10] * t5 / 0.10e2;
  double t14 = u[2] * t5 / 0.10e2;
  double t17 = u[8] * t5 / 0.10e2;
  
  // Assembly the BLw matrix.
  
  BLw[0][0] = 0.0e0;
  BLw[0][1] = 0.0e0;
  BLw[0][2] = t12;
  BLw[0][3] = 0.0e0;
  BLw[0][4] = -t14 + 0.2e1 / 0.15e2 * u[4] + t17 - u[10] / 0.30e2;
  BLw[0][5] = 0.0e0;
  BLw[0][6] = 0.0e0;
  BLw[0][7] = 0.0e0;
  BLw[0][8] = -t12;
  BLw[0][9] = 0.0e0;
  BLw[0][10] = -t14 - u[4] / 0.30e2 + t17 + 0.2e1 / 0.15e2 * u[10];
  BLw[0][11] = 0.0e0;
}

// ================================= BLthetaMat ============================
//
// This method computes the [BLtheta] matrix (nonlinear) related to membrane
// strains. Obtained from angle of twist.
//
//   L       - element length                                          (in)
//   Ip      - polar moment of inertia of the cross section            (in)
//   A       - area of the cross section                               (in)
//   u       - nodal displacement vector                               (in)
//   BLtheta - strain-displacement matrix (1x12)                      (out)
//
void cFrame3DCRTL :: BLthetaMat(double L, double Ip, double A,
                                cVector &u, cMatrix &BLtheta)
{
  // Auxiliary variables.
  
  double t1 = L * L;
  double t2 = 0.1e1 / t1;
  double t5 = Ip / A;
  double t9 = u[3] * t2 * t5 - u[9] * t2 * t5;
  
  // Assembly the BLtheta matrix.
  
  BLtheta[0][0] = 0.0e0;
  BLtheta[0][1] = 0.0e0;
  BLtheta[0][2] = 0.0e0;
  BLtheta[0][3] = t9;
  BLtheta[0][4] = 0.0e0;
  BLtheta[0][5] = 0.0e0;
  BLtheta[0][6] = 0.0e0;
  BLtheta[0][7] = 0.0e0;
  BLtheta[0][8] = 0.0e0;
  BLtheta[0][9] = -t9;
  BLtheta[0][10] = 0.0e0;
  BLtheta[0][11] = 0.0e0;
}

// =================================== BLMat ===============================
//
// This method computes the nonlinear strain-displacement matrix [BL].
//
//   L  - element length                                               (in)
//   Ip - polar moment of inertia of the cross section                 (in)
//   A  - area of the cross section                                    (in)
//   u  - nodal displacement vector                                    (in)
//   Bm - strain-displacement matrix (1x12)                           (out)
//
void cFrame3DCRTL :: BLMat(double L, double Ip, double A,
                           cVector &u, cMatrix &BL)
{
  // Auxiliary variables.
  
  double t1 = L * L;
  double t2 = 0.1e1 / t1;
  double t5 = 0.1e1 / L;
  double t12 = 0.6e1 / 0.5e1 * u[1] * t2 + u[5] * t5 / 0.10e2 - 0.6e1 / 0.5e1 * u[7] * t2 + u[11] * t5 / 0.10e2;
  double t21 = 0.6e1 / 0.5e1 * u[2] * t2 - u[4] * t5 / 0.10e2 - 0.6e1 / 0.5e1 * u[8] * t2 - u[10] * t5 / 0.10e2;
  double t24 = Ip / A;
  double t28 = u[3] * t2 * t24 - u[9] * t2 * t24;
  double t30 = u[2] * t5 / 0.10e2;
  double t33 = u[8] * t5 / 0.10e2;
  double t37 = u[1] * t5 / 0.10e2;
  double t40 = u[7] * t5 / 0.10e2;
  
  // Assembly the [BL] matrix.
  
  BL[0][0] = 0.0e0;
  BL[0][1] = t12;
  BL[0][2] = t21;
  BL[0][3] = t28;
  BL[0][4] = -t30 + 0.2e1 / 0.15e2 * u[4] + t33 - u[10] / 0.30e2;
  BL[0][5] = t37 + 0.2e1 / 0.15e2 * u[5] - t40 - u[11] / 0.30e2;
  BL[0][6] = 0.0e0;
  BL[0][7] = -t12;
  BL[0][8] = -t21;
  BL[0][9] = -t28;
  BL[0][10] = -t30 - u[4] / 0.30e2 + t33 + 0.2e1 / 0.15e2 * u[10];
  BL[0][11] = t37 - u[5] / 0.30e2 - t40 + 0.2e1 / 0.15e2 * u[11];
}

// =============================== ConvStiffMat ============================
//
// This method computes the element nonlinear matrix stiffness obtained
// from integration [kl] = int([Bbar]t[C][Bbar]), where section constitutive
// matrix [C] is diagonal.
//
//   L  - element length                                               (in)
//   u  - nodal displacement vector                                    (in)
//   kl - nonlinear stiffness matrix (12x12)                          (out)
//
void cFrame3DCRTL :: ConvStiffMat(double L, cVector &u, cMatrix &kl)
{
  // Compute the stiffness and geometric properties.
  
  cMatMec *mat = Section->GetMaterial( )->GetMecProp( );
  double matpar[2];
  mat->GetParam(matpar);
  double E  = matpar[0];
  double nu = matpar[1];
  double G  = 0.5*E/(1 + nu);
  
  double A  = Section->GetA( );
  double Ip = Section->GetJt( );   // ???
  double EA  = E*A;
  double EIz = E*Section->GetIz( );
  double EIy = E*Section->GetIy( );
  double GJ  = G*Ip;
    
  // Auxiliary variables.
  
  double t1 = 0.1e1 / L;
  double t2 = t1 * EA;
  double t3 = L * L;
  double t4 = 0.1e1 / t3;
  double t13 = 0.6e1 / 0.5e1 * u[1] * t4 + u[5] * t1 / 0.10e2 - 0.6e1 / 0.5e1 * u[7] * t4 + u[11] * t1 / 0.10e2;
  double t14 = EA * t13;
  double t23 = 0.6e1 / 0.5e1 * u[2] * t4 - u[4] * t1 / 0.10e2 - 0.6e1 / 0.5e1 * u[8] * t4 - u[10] * t1 / 0.10e2;
  double t24 = EA * t23;
  double t27 = Ip / A;
  double t31 = u[3] * t4 * t27 - u[9] * t4 * t27;
  double t32 = EA * t31;
  double t34 = u[2] * t1 / 0.10e2;
  double t37 = u[8] * t1 / 0.10e2;
  double t39 = -t34 + 0.2e1 / 0.15e2 * u[4] + t37 - u[10] / 0.30e2;
  double t40 = t39 * EA;
  double t42 = u[1] * t1 / 0.10e2;
  double t45 = u[7] * t1 / 0.10e2;
  double t47 = t42 + 0.2e1 / 0.15e2 * u[5] - t45 - u[11] / 0.30e2;
  double t48 = EA * t47;
  double t49 = -EA * t13;
  double t50 = -EA * t23;
  double t51 = -EA * t31;
  double t54 = -t34 - u[4] / 0.30e2 + t37 + 0.2e1 / 0.15e2 * u[10];
  double t55 = EA * t54;
  double t58 = t42 - u[5] / 0.30e2 - t45 + 0.2e1 / 0.15e2 * u[11];
  double t59 = t58 * EA;
  double t61 = 0.1e1 / t3 / L;
  double t63 = 0.12e2 * t61 * EIz;
  double t64 = t13 * t13;
  double t67 = t63 + t64 * EA * L;
  double t69 = t14 * t23 * L;
  double t70 = t31 * L;
  double t71 = t14 * t70;
  double t72 = t39 * L;
  double t73 = t14 * t72;
  double t75 = 0.6e1 * t4 * EIz;
  double t76 = t47 * L;
  double t78 = t75 + t14 * t76;
  double t79 = -t13 * L;
  double t81 = -t63 + t14 * t79;
  double t82 = -t23 * L;
  double t83 = t14 * t82;
  double t84 = -t31 * L;
  double t85 = t14 * t84;
  double t86 = t54 * L;
  double t87 = t14 * t86;
  double t88 = t58 * L;
  double t90 = t75 + t14 * t88;
  double t92 = 0.12e2 * t61 * EIy;
  double t93 = t23 * t23;
  double t96 = t92 + t93 * EA * L;
  double t97 = t24 * t70;
  double t99 = 0.6e1 * t4 * EIy;
  double t101 = -t99 + t24 * t72;
  double t102 = t24 * t76;
  double t103 = t24 * t79;
  double t105 = -t92 + t24 * t82;
  double t106 = t24 * t84;
  double t108 = -t99 + t24 * t86;
  double t109 = t24 * t88;
  double t110 = t31 * t31;
  double t113 = t1 * GJ;
  double t114 = t110 * EA * L + t113;
  double t115 = t32 * t72;
  double t116 = t32 * t76;
  double t117 = t32 * t79;
  double t118 = t32 * t82;
  double t120 = t32 * t84 - t113;
  double t121 = t32 * t86;
  double t122 = t32 * t88;
  double t123 = t1 * EIy;
  double t124 = 0.4e1 * t123;
  double t125 = t39 * t39;
  double t129 = t40 * t76;
  double t130 = t40 * t79;
  double t132 = t99 + t40 * t82;
  double t133 = t40 * t84;
  double t136 = 0.2e1 * t123 + t40 * t86;
  double t137 = t40 * t88;
  double t138 = t1 * EIz;
  double t139 = 0.4e1 * t138;
  double t140 = t47 * t47;
  double t145 = -t75 + t48 * t79;
  double t146 = t48 * t82;
  double t147 = t48 * t84;
  double t148 = t48 * t86;
  double t151 = 0.2e1 * t138 + t48 * t88;
  double t152 = t49 * t82;
  double t153 = t49 * t84;
  double t154 = t49 * t86;
  double t156 = -t75 + t49 * t88;
  double t157 = t50 * t84;
  double t159 = t99 + t50 * t86;
  double t160 = t50 * t88;
  double t161 = t51 * t86;
  double t162 = t51 * t88;
  double t163 = t54 * t54;
  double t167 = t55 * t88;
  double t168 = t58 * t58;
  
  // Assembly the [kl] matrix.
  
  kl[0][0] = t2;
  kl[0][1] = -t14;
  kl[0][2] = -t24;
  kl[0][3] = -t32;
  kl[0][4] = -t40;
  kl[0][5] = -t48;
  kl[0][6] = -t2;
  kl[0][7] = -t49;
  kl[0][8] = -t50;
  kl[0][9] = -t51;
  kl[0][10] = -t55;
  kl[0][11] = -t59;
  kl[1][0] = -t14;
  kl[1][1] = t67;
  kl[1][2] = t69;
  kl[1][3] = t71;
  kl[1][4] = t73;
  kl[1][5] = t78;
  kl[1][6] = t14;
  kl[1][7] = t81;
  kl[1][8] = t83;
  kl[1][9] = t85;
  kl[1][10] = t87;
  kl[1][11] = t90;
  kl[2][0] = -t24;
  kl[2][1] = t69;
  kl[2][2] = t96;
  kl[2][3] = t97;
  kl[2][4] = t101;
  kl[2][5] = t102;
  kl[2][6] = t24;
  kl[2][7] = t103;
  kl[2][8] = t105;
  kl[2][9] = t106;
  kl[2][10] = t108;
  kl[2][11] = t109;
  kl[3][0] = -t32;
  kl[3][1] = t71;
  kl[3][2] = t97;
  kl[3][3] = t114;
  kl[3][4] = t115;
  kl[3][5] = t116;
  kl[3][6] = t32;
  kl[3][7] = t117;
  kl[3][8] = t118;
  kl[3][9] = t120;
  kl[3][10] = t121;
  kl[3][11] = t122;
  kl[4][0] = -t40;
  kl[4][1] = t73;
  kl[4][2] = t101;
  kl[4][3] = t115;
  kl[4][4] = t124 + t125 * EA * L;
  kl[4][5] = t129;
  kl[4][6] = t40;
  kl[4][7] = t130;
  kl[4][8] = t132;
  kl[4][9] = t133;
  kl[4][10] = t136;
  kl[4][11] = t137;
  kl[5][0] = -t48;
  kl[5][1] = t78;
  kl[5][2] = t102;
  kl[5][3] = t116;
  kl[5][4] = t129;
  kl[5][5] = t139 + t140 * EA * L;
  kl[5][6] = t48;
  kl[5][7] = t145;
  kl[5][8] = t146;
  kl[5][9] = t147;
  kl[5][10] = t148;
  kl[5][11] = t151;
  kl[6][0] = -t2;
  kl[6][1] = t14;
  kl[6][2] = t24;
  kl[6][3] = t32;
  kl[6][4] = t40;
  kl[6][5] = t48;
  kl[6][6] = t2;
  kl[6][7] = t49;
  kl[6][8] = t50;
  kl[6][9] = t51;
  kl[6][10] = t55;
  kl[6][11] = t59;
  kl[7][0] = -t49;
  kl[7][1] = t81;
  kl[7][2] = t103;
  kl[7][3] = t117;
  kl[7][4] = t130;
  kl[7][5] = t145;
  kl[7][6] = t49;
  kl[7][7] = t67;
  kl[7][8] = t152;
  kl[7][9] = t153;
  kl[7][10] = t154;
  kl[7][11] = t156;
  kl[8][0] = -t50;
  kl[8][1] = t83;
  kl[8][2] = t105;
  kl[8][3] = t118;
  kl[8][4] = t132;
  kl[8][5] = t146;
  kl[8][6] = t50;
  kl[8][7] = t152;
  kl[8][8] = t96;
  kl[8][9] = t157;
  kl[8][10] = t159;
  kl[8][11] = t160;
  kl[9][0] = -t51;
  kl[9][1] = t85;
  kl[9][2] = t106;
  kl[9][3] = t120;
  kl[9][4] = t133;
  kl[9][5] = t147;
  kl[9][6] = t51;
  kl[9][7] = t153;
  kl[9][8] = t157;
  kl[9][9] = t114;
  kl[9][10] = t161;
  kl[9][11] = t162;
  kl[10][0] = -t55;
  kl[10][1] = t87;
  kl[10][2] = t108;
  kl[10][3] = t121;
  kl[10][4] = t136;
  kl[10][5] = t148;
  kl[10][6] = t55;
  kl[10][7] = t154;
  kl[10][8] = t159;
  kl[10][9] = t161;
  kl[10][10] = t124 + t163 * EA * L;
  kl[10][11] = t167;
  kl[11][0] = -t59;
  kl[11][1] = t90;
  kl[11][2] = t109;
  kl[11][3] = t122;
  kl[11][4] = t137;
  kl[11][5] = t151;
  kl[11][6] = t59;
  kl[11][7] = t156;
  kl[11][8] = t160;
  kl[11][9] = t162;
  kl[11][10] = t167;
  kl[11][11] = t139 + t168 * EA * L;
}

// =============================== CoupStiffMat ============================
//
// This method computes the element nonlinear matrix stiffness obtained
// from integration [kl] = int([Bbar]t[C][Bbar]), where section constitutive
// matrix [C] is fully coupled.
//
//   L  - element length                                               (in)
//   u  - nodal displacement vector                                    (in)
//   kl - nonlinear stiffness matrix (12x12)                          (out)
//
void cFrame3DCRTL :: CoupStiffMat(double L, cVector &u, cMatrix &kl)
{
  // Geometric properties.

  double A  = Section->GetA( );
  double Ip = Section->GetJt( );   // ???
  
  // Section stiffness - fully coupled 4x4 matrix.
  
  cMatrix Cb(4, 4);
  Section->SecStiffness3D(Cb);
  double C11 = Cb[0][0]; double C12 = Cb[0][1]; double C13 = Cb[0][2]; double C14 = Cb[0][3];
  double C21 = Cb[1][0]; double C22 = Cb[1][1]; double C23 = Cb[1][2]; double C24 = Cb[1][3];
  double C31 = Cb[2][0]; double C32 = Cb[2][1]; double C33 = Cb[2][2]; double C34 = Cb[2][3];
  double C41 = Cb[3][0]; double C42 = Cb[3][1]; double C43 = Cb[3][2]; double C44 = Cb[3][3];

  // Auxiliary variables (Maple code).

  double t1 = 0.1e1 / L;
  double t2 = t1 * C11;
  double t3 = L * L;
  double t4 = 0.1e1 / t3;
  double t13 = 0.6e1 / 0.5e1 * u[1] * t4 + u[5] * t1 / 0.10e2 - 0.6e1 / 0.5e1 * u[7] * t4 + u[11] * t1 / 0.10e2;
  double t14 = t13 * C11;
  double t23 = 0.6e1 / 0.5e1 * u[2] * t4 - u[4] * t1 / 0.10e2 - 0.6e1 / 0.5e1 * u[8] * t4 - u[10] * t1 / 0.10e2;
  double t24 = t23 * C11;
  double t27 = Ip / A;
  double t31 = u[3] * t4 * t27 - u[9] * t4 * t27;
  double t32 = t31 * C11;
  double t33 = t1 * C14;
  double t34 = -t32 + t33;
  double t35 = t1 * C13;
  double t37 = u[2] * t1 / 0.10e2;
  double t40 = u[8] * t1 / 0.10e2;
  double t42 = -t37 + 0.2e1 / 0.15e2 * u[4] + t40 - u[10] / 0.30e2;
  double t43 = t42 * C11;
  double t44 = -t35 - t43;
  double t45 = t1 * C12;
  double t47 = u[1] * t1 / 0.10e2;
  double t50 = u[7] * t1 / 0.10e2;
  double t52 = t47 + 0.2e1 / 0.15e2 * u[5] - t50 - u[11] / 0.30e2;
  double t53 = t52 * C11;
  double t54 = t45 - t53;
  double t55 = -t13 * C11;
  double t56 = -t23 * C11;
  double t57 = -t31 * C11;
  double t58 = -t57 - t33;
  double t61 = -t37 - u[4] / 0.30e2 + t40 + 0.2e1 / 0.15e2 * u[10];
  double t62 = t61 * C11;
  double t63 = t35 - t62;
  double t66 = t47 - u[5] / 0.30e2 - t50 + 0.2e1 / 0.15e2 * u[11];
  double t67 = t66 * C11;
  double t68 = -t45 - t67;
  double t70 = 0.1e1 / t3 / L;
  double t71 = t70 * C22;
  double t72 = 0.48e2 * t71;
  double t73 = t70 * C21;
  double t75 = 0.12e2 * t73 * t13;
  double t76 = t13 * C12;
  double t77 = t4 * C22;
  double t78 = 0.6e1 * t77;
  double t79 = t76 - t78;
  double t81 = 0.12e2 * t79 * t70;
  double t82 = t3 * t3;
  double t84 = 0.1e1 / t82 / L;
  double t86 = 0.72e2 * t84 * C22;
  double t90 = t4 * C21;
  double t91 = 0.6e1 * t90;
  double t92 = t14 - t91;
  double t96 = 0.6e1 * t79 * t1;
  double t98 = t70 * C23;
  double t99 = 0.48e2 * t98;
  double t101 = 0.12e2 * t73 * t23;
  double t102 = t13 * C13;
  double t103 = t4 * C23;
  double t104 = 0.6e1 * t103;
  double t105 = t102 - t104;
  double t107 = 0.12e2 * t105 * t70;
  double t109 = 0.72e2 * t84 * C23;
  double t116 = 0.6e1 * t105 * t1;
  double t119 = 0.1e1 / t82;
  double t120 = t119 * C24;
  double t121 = t73 * t31 - t120;
  double t126 = t13 * C14;
  double t128 = 0.6e1 * t4 * C24;
  double t130 = 0.48e2 * t103;
  double t132 = 0.12e2 * t73 * t42;
  double t134 = 0.6e1 * t105 * t4;
  double t135 = t119 * C23;
  double t136 = 0.48e2 * t135;
  double t144 = 0.48e2 * t77;
  double t146 = 0.12e2 * t73 * t52;
  double t148 = 0.6e1 * t79 * t4;
  double t149 = t119 * C22;
  double t150 = 0.48e2 * t149;
  double t159 = -0.12e2 * t73 * t13;
  double t167 = -0.12e2 * t73 * t23;
  double t175 = -t73 * t31 + t120;
  double t181 = 0.36e2 * t103;
  double t183 = 0.12e2 * t73 * t61;
  double t184 = 0.24e2 * t135;
  double t192 = 0.36e2 * t77;
  double t194 = 0.12e2 * t73 * t66;
  double t195 = 0.24e2 * t149;
  double t203 = t70 * C32;
  double t204 = 0.48e2 * t203;
  double t205 = t70 * C31;
  double t207 = 0.12e2 * t205 * t13;
  double t208 = t23 * C12;
  double t209 = t4 * C32;
  double t210 = 0.6e1 * t209;
  double t211 = t208 - t210;
  double t213 = 0.12e2 * t211 * t70;
  double t215 = 0.72e2 * t84 * C32;
  double t219 = t4 * C31;
  double t220 = 0.6e1 * t219;
  double t221 = t24 - t220;
  double t225 = 0.6e1 * t211 * t1;
  double t227 = t70 * C33;
  double t228 = 0.48e2 * t227;
  double t230 = 0.12e2 * t205 * t23;
  double t231 = t23 * C13;
  double t232 = t4 * C33;
  double t233 = 0.6e1 * t232;
  double t234 = t231 - t233;
  double t236 = 0.12e2 * t234 * t70;
  double t238 = 0.72e2 * t84 * C33;
  double t245 = 0.6e1 * t234 * t1;
  double t248 = t119 * C34;
  double t249 = t205 * t31 - t248;
  double t254 = t23 * C14;
  double t256 = 0.6e1 * t4 * C34;
  double t258 = 0.48e2 * t232;
  double t260 = 0.12e2 * t205 * t42;
  double t262 = 0.6e1 * t234 * t4;
  double t263 = t119 * C33;
  double t264 = 0.48e2 * t263;
  double t272 = 0.48e2 * t209;
  double t274 = 0.12e2 * t205 * t52;
  double t276 = 0.6e1 * t211 * t4;
  double t277 = t119 * C32;
  double t278 = 0.48e2 * t277;
  double t287 = -0.12e2 * t205 * t13;
  double t295 = -0.12e2 * t205 * t23;
  double t303 = -t205 * t31 + t248;
  double t309 = 0.36e2 * t232;
  double t311 = 0.12e2 * t205 * t61;
  double t312 = 0.24e2 * t263;
  double t320 = 0.36e2 * t209;
  double t322 = 0.12e2 * t205 * t66;
  double t323 = 0.24e2 * t277;
  double t331 = t1 * C41;
  double t332 = -t32 + t331;
  double t339 = t31 * C14;
  double t340 = t1 * C44;
  double t342 = t31 * C13;
  double t343 = t1 * C43;
  double t347 = t31 * C12;
  double t348 = t1 * C42;
  double t365 = t1 * C31;
  double t366 = -t365 - t43;
  double t367 = 0.24e2 * t209;
  double t369 = 0.6e1 * t219 * t13;
  double t370 = t42 * C12;
  double t371 = t1 * C32;
  double t373 = t370 + 0.4e1 * t371;
  double t375 = 0.12e2 * t373 * t70;
  double t376 = 0.36e2 * t277;
  double t381 = t43 + 0.4e1 * t365;
  double t385 = 0.6e1 * t373 * t1;
  double t387 = 0.24e2 * t232;
  double t389 = 0.6e1 * t219 * t23;
  double t390 = t42 * C13;
  double t391 = t1 * C33;
  double t393 = t390 + 0.4e1 * t391;
  double t395 = 0.12e2 * t393 * t70;
  double t396 = 0.36e2 * t263;
  double t403 = 0.6e1 * t393 * t1;
  double t406 = t70 * C34;
  double t409 = 0.3e1 * (-t219 * t31 + t406) * t3;
  double t412 = t42 * C14;
  double t413 = t1 * C34;
  double t414 = 0.4e1 * t413;
  double t418 = 0.6e1 * t219 * t42;
  double t420 = 0.6e1 * t393 * t4;
  double t421 = 0.24e2 * t227;
  double t431 = 0.6e1 * t219 * t52;
  double t433 = 0.6e1 * t373 * t4;
  double t434 = 0.24e2 * t203;
  double t443 = -0.6e1 * t219 * t13;
  double t451 = -0.6e1 * t219 * t23;
  double t461 = 0.3e1 * (t219 * t31 - t406) * t3;
  double t465 = 0.20e2 * t391;
  double t467 = 0.6e1 * t219 * t61;
  double t468 = 0.12e2 * t227;
  double t476 = 0.20e2 * t371;
  double t478 = 0.6e1 * t219 * t66;
  double t479 = 0.12e2 * t203;
  double t487 = t1 * C21;
  double t488 = t487 - t53;
  double t489 = 0.24e2 * t77;
  double t491 = 0.6e1 * t90 * t13;
  double t492 = t52 * C12;
  double t493 = t1 * C22;
  double t495 = t492 - 0.4e1 * t493;
  double t497 = 0.12e2 * t495 * t70;
  double t498 = 0.36e2 * t149;
  double t503 = t53 - 0.4e1 * t487;
  double t507 = 0.6e1 * t495 * t1;
  double t509 = 0.24e2 * t103;
  double t511 = 0.6e1 * t90 * t23;
  double t512 = t52 * C13;
  double t513 = t1 * C23;
  double t515 = t512 - 0.4e1 * t513;
  double t517 = 0.12e2 * t515 * t70;
  double t518 = 0.36e2 * t135;
  double t525 = 0.6e1 * t515 * t1;
  double t528 = t70 * C24;
  double t531 = 0.3e1 * (t90 * t31 - t528) * t3;
  double t534 = t52 * C14;
  double t535 = t1 * C24;
  double t536 = 0.4e1 * t535;
  double t540 = 0.6e1 * t90 * t42;
  double t542 = 0.6e1 * t515 * t4;
  double t543 = 0.24e2 * t98;
  double t553 = 0.6e1 * t90 * t52;
  double t555 = 0.6e1 * t495 * t4;
  double t556 = 0.24e2 * t71;
  double t565 = -0.6e1 * t90 * t13;
  double t573 = -0.6e1 * t90 * t23;
  double t583 = 0.3e1 * (-t90 * t31 + t528) * t3;
  double t587 = 0.20e2 * t513;
  double t589 = 0.6e1 * t90 * t61;
  double t590 = 0.12e2 * t98;
  double t598 = 0.20e2 * t493;
  double t600 = 0.6e1 * t90 * t66;
  double t601 = 0.12e2 * t71;
  double t609 = -t13 * C12;
  double t610 = t609 + t78;
  double t612 = 0.12e2 * t610 * t70;
  double t616 = t55 + t91;
  double t620 = 0.6e1 * t610 * t1;
  double t622 = -t13 * C13;
  double t623 = t622 + t104;
  double t625 = 0.12e2 * t623 * t70;
  double t632 = 0.6e1 * t623 * t1;
  double t638 = -t13 * C14;
  double t641 = 0.6e1 * t623 * t4;
  double t650 = 0.6e1 * t610 * t4;
  double t689 = -t23 * C12;
  double t690 = t689 + t210;
  double t692 = 0.12e2 * t690 * t70;
  double t696 = t56 + t220;
  double t700 = 0.6e1 * t690 * t1;
  double t702 = -t23 * C13;
  double t703 = t702 + t233;
  double t705 = 0.12e2 * t703 * t70;
  double t712 = 0.6e1 * t703 * t1;
  double t718 = -t23 * C14;
  double t721 = 0.6e1 * t703 * t4;
  double t730 = 0.6e1 * t690 * t4;
  double t769 = -t57 - t331;
  double t776 = -t31 * C14;
  double t778 = -t31 * C13;
  double t782 = -t31 * C12;
  double t799 = t365 - t62;
  double t800 = t61 * C12;
  double t802 = t800 + 0.2e1 * t371;
  double t804 = 0.12e2 * t802 * t70;
  double t809 = t62 + 0.2e1 * t365;
  double t813 = 0.6e1 * t802 * t1;
  double t815 = t61 * C13;
  double t817 = t815 + 0.2e1 * t391;
  double t819 = 0.12e2 * t817 * t70;
  double t826 = 0.6e1 * t817 * t1;
  double t830 = t61 * C14;
  double t831 = 0.2e1 * t413;
  double t834 = 0.6e1 * t817 * t4;
  double t843 = 0.6e1 * t802 * t4;
  double t882 = -t487 - t67;
  double t883 = t66 * C12;
  double t885 = t883 - 0.2e1 * t493;
  double t887 = 0.12e2 * t885 * t70;
  double t892 = t67 - 0.2e1 * t487;
  double t896 = 0.6e1 * t885 * t1;
  double t898 = t66 * C13;
  double t900 = t898 - 0.2e1 * t513;
  double t902 = 0.12e2 * t900 * t70;
  double t909 = 0.6e1 * t900 * t1;
  double t913 = t66 * C14;
  double t914 = 0.2e1 * t535;
  double t917 = 0.6e1 * t900 * t4;
  double t926 = 0.6e1 * t885 * t4;
  
  // Assembly the [kl] matrix.
  
  kl[0][0] = t2;
  kl[0][1] = -t14;
  kl[0][2] = -t24;
  kl[0][3] = t34;
  kl[0][4] = t44;
  kl[0][5] = t54;
  kl[0][6] = -t2;
  kl[0][7] = -t55;
  kl[0][8] = -t56;
  kl[0][9] = t58;
  kl[0][10] = t63;
  kl[0][11] = t68;
  kl[1][0] = -t14;
  kl[1][1] = t72 + (t75 + t81 - t86) * t3 / 0.2e1 + t92 * t13 * L - t96;
  kl[1][2] = t99 + (t101 + t107 - t109) * t3 / 0.2e1 + t92 * t23 * L - t116;
  kl[1][3] = 0.6e1 * t121 * t3 + t92 * t31 * L - t126 + t128;
  kl[1][4] = -t130 + (t132 - t134 + t136) * t3 / 0.2e1 + t92 * t42 * L + 0.4e1 * t102;
  kl[1][5] = t144 + (t146 + t148 - t150) * t3 / 0.2e1 + t92 * t52 * L - 0.4e1 * t76;
  kl[1][6] = t14;
  kl[1][7] = -t72 + (t159 - t81 + t86) * t3 / 0.2e1 - t92 * t13 * L + t96;
  kl[1][8] = -t99 + (t167 - t107 + t109) * t3 / 0.2e1 - t92 * t23 * L + t116;
  kl[1][9] = 0.6e1 * t175 * t3 - t92 * t31 * L + t126 - t128;
  kl[1][10] = -t181 + (t183 - t134 + t184) * t3 / 0.2e1 + t92 * t61 * L + 0.2e1 * t102;
  kl[1][11] = t192 + (t194 + t148 - t195) * t3 / 0.2e1 + t92 * t66 * L - 0.2e1 * t76;
  kl[2][0] = -t24;
  kl[2][1] = t204 + (t207 + t213 - t215) * t3 / 0.2e1 + t221 * t13 * L - t225;
  kl[2][2] = t228 + (t230 + t236 - t238) * t3 / 0.2e1 + t221 * t23 * L - t245;
  kl[2][3] = 0.6e1 * t249 * t3 + t221 * t31 * L - t254 + t256;
  kl[2][4] = -t258 + (t260 - t262 + t264) * t3 / 0.2e1 + t221 * t42 * L + 0.4e1 * t231;
  kl[2][5] = t272 + (t274 + t276 - t278) * t3 / 0.2e1 + t221 * t52 * L - 0.4e1 * t208;
  kl[2][6] = t24;
  kl[2][7] = -t204 + (t287 - t213 + t215) * t3 / 0.2e1 - t221 * t13 * L + t225;
  kl[2][8] = -t228 + (t295 - t236 + t238) * t3 / 0.2e1 - t221 * t23 * L + t245;
  kl[2][9] = 0.6e1 * t303 * t3 - t221 * t31 * L + t254 - t256;
  kl[2][10] = -t309 + (t311 - t262 + t312) * t3 / 0.2e1 + t221 * t61 * L + 0.2e1 * t231;
  kl[2][11] = t320 + (t322 + t276 - t323) * t3 / 0.2e1 + t221 * t66 * L - 0.2e1 * t208;
  kl[3][0] = t332;
  kl[3][1] = -t332 * t13 * L;
  kl[3][2] = -t332 * t23 * L;
  kl[3][3] = -t332 * t31 * L - t339 + t340;
  kl[3][4] = t342 - t343 - t332 * t42 * L;
  kl[3][5] = -t347 + t348 - t332 * t52 * L;
  kl[3][6] = -t332;
  kl[3][7] = t332 * t13 * L;
  kl[3][8] = t332 * t23 * L;
  kl[3][9] = t332 * t31 * L + t339 - t340;
  kl[3][10] = -t342 + t343 - t332 * t61 * L;
  kl[3][11] = t347 - t348 - t332 * t66 * L;
  kl[4][0] = t366;
  kl[4][1] = -t367 + (-t369 + t375 + t376) * t3 / 0.2e1 + t381 * t13 * L - t385;
  kl[4][2] = -t387 + (-t389 + t395 + t396) * t3 / 0.2e1 + t381 * t23 * L - t403;
  kl[4][3] = t409 + t381 * t31 * L - t412 - t414;
  kl[4][4] = 0.28e2 * t391 + (-t418 - t420 - t421) * t3 / 0.2e1 + t381 * t42 * L + 0.4e1 * t390;
  kl[4][5] = -0.28e2 * t371 + (-t431 + t433 + t434) * t3 / 0.2e1 + t381 * t52 * L - 0.4e1 * t370;
  kl[4][6] = -t366;
  kl[4][7] = t367 + (-t443 - t375 - t376) * t3 / 0.2e1 - t381 * t13 * L + t385;
  kl[4][8] = t387 + (-t451 - t395 - t396) * t3 / 0.2e1 - t381 * t23 * L + t403;
  kl[4][9] = t461 - t381 * t31 * L + t412 + t414;
  kl[4][10] = t465 + (-t467 - t420 - t468) * t3 / 0.2e1 + t381 * t61 * L + 0.2e1 * t390;
  kl[4][11] = -t476 + (-t478 + t433 + t479) * t3 / 0.2e1 + t381 * t66 * L - 0.2e1 * t370;
  kl[5][0] = t488;
  kl[5][1] = t489 + (t491 + t497 - t498) * t3 / 0.2e1 + t503 * t13 * L - t507;
  kl[5][2] = t509 + (t511 + t517 - t518) * t3 / 0.2e1 + t503 * t23 * L - t525;
  kl[5][3] = t531 + t503 * t31 * L - t534 + t536;
  kl[5][4] = -0.28e2 * t513 + (t540 - t542 + t543) * t3 / 0.2e1 + t503 * t42 * L + 0.4e1 * t512;
  kl[5][5] = 0.28e2 * t493 + (t553 + t555 - t556) * t3 / 0.2e1 + t503 * t52 * L - 0.4e1 * t492;
  kl[5][6] = -t488;
  kl[5][7] = -t489 + (t565 - t497 + t498) * t3 / 0.2e1 - t503 * t13 * L + t507;
  kl[5][8] = -t509 + (t573 - t517 + t518) * t3 / 0.2e1 - t503 * t23 * L + t525;
  kl[5][9] = t583 - t503 * t31 * L + t534 - t536;
  kl[5][10] = -t587 + (t589 - t542 + t590) * t3 / 0.2e1 + t503 * t61 * L + 0.2e1 * t512;
  kl[5][11] = t598 + (t600 + t555 - t601) * t3 / 0.2e1 + t503 * t66 * L - 0.2e1 * t492;
  kl[6][0] = -t2;
  kl[6][1] = t14;
  kl[6][2] = t24;
  kl[6][3] = -t34;
  kl[6][4] = -t44;
  kl[6][5] = -t54;
  kl[6][6] = t2;
  kl[6][7] = t55;
  kl[6][8] = t56;
  kl[6][9] = -t58;
  kl[6][10] = -t63;
  kl[6][11] = -t68;
  kl[7][0] = -t55;
  kl[7][1] = -t72 + (-t75 + t612 + t86) * t3 / 0.2e1 + t616 * t13 * L - t620;
  kl[7][2] = -t99 + (-t101 + t625 + t109) * t3 / 0.2e1 + t616 * t23 * L - t632;
  kl[7][3] = -0.6e1 * t121 * t3 + t616 * t31 * L - t638 - t128;
  kl[7][4] = t130 + (-t132 - t641 - t136) * t3 / 0.2e1 + t616 * t42 * L + 0.4e1 * t622;
  kl[7][5] = -t144 + (-t146 + t650 + t150) * t3 / 0.2e1 + t616 * t52 * L - 0.4e1 * t609;
  kl[7][6] = t55;
  kl[7][7] = t72 + (-t159 - t612 - t86) * t3 / 0.2e1 - t616 * t13 * L + t620;
  kl[7][8] = t99 + (-t167 - t625 - t109) * t3 / 0.2e1 - t616 * t23 * L + t632;
  kl[7][9] = -0.6e1 * t175 * t3 - t616 * t31 * L + t638 + t128;
  kl[7][10] = t181 + (-t183 - t641 - t184) * t3 / 0.2e1 + t616 * t61 * L + 0.2e1 * t622;
  kl[7][11] = -t192 + (-t194 + t650 + t195) * t3 / 0.2e1 + t616 * t66 * L - 0.2e1 * t609;
  kl[8][0] = -t56;
  kl[8][1] = -t204 + (-t207 + t692 + t215) * t3 / 0.2e1 + t696 * t13 * L - t700;
  kl[8][2] = -t228 + (-t230 + t705 + t238) * t3 / 0.2e1 + t696 * t23 * L - t712;
  kl[8][3] = -0.6e1 * t249 * t3 + t696 * t31 * L - t718 - t256;
  kl[8][4] = t258 + (-t260 - t721 - t264) * t3 / 0.2e1 + t696 * t42 * L + 0.4e1 * t702;
  kl[8][5] = -t272 + (-t274 + t730 + t278) * t3 / 0.2e1 + t696 * t52 * L - 0.4e1 * t689;
  kl[8][6] = t56;
  kl[8][7] = t204 + (-t287 - t692 - t215) * t3 / 0.2e1 - t696 * t13 * L + t700;
  kl[8][8] = t228 + (-t295 - t705 - t238) * t3 / 0.2e1 - t696 * t23 * L + t712;
  kl[8][9] = -0.6e1 * t303 * t3 - t696 * t31 * L + t718 + t256;
  kl[8][10] = t309 + (-t311 - t721 - t312) * t3 / 0.2e1 + t696 * t61 * L + 0.2e1 * t702;
  kl[8][11] = -t320 + (-t322 + t730 + t323) * t3 / 0.2e1 + t696 * t66 * L - 0.2e1 * t689;
  kl[9][0] = t769;
  kl[9][1] = -t769 * t13 * L;
  kl[9][2] = -t769 * t23 * L;
  kl[9][3] = -t769 * t31 * L - t776 - t340;
  kl[9][4] = t778 + t343 - t769 * t42 * L;
  kl[9][5] = -t782 - t348 - t769 * t52 * L;
  kl[9][6] = -t769;
  kl[9][7] = t769 * t13 * L;
  kl[9][8] = t769 * t23 * L;
  kl[9][9] = t769 * t31 * L + t776 + t340;
  kl[9][10] = -t778 - t343 - t769 * t61 * L;
  kl[9][11] = t782 + t348 - t769 * t66 * L;
  kl[10][0] = t799;
  kl[10][1] = -t367 + (-t369 + t804 + t376) * t3 / 0.2e1 + t809 * t13 * L - t813;
  kl[10][2] = -t387 + (-t389 + t819 + t396) * t3 / 0.2e1 + t809 * t23 * L - t826;
  kl[10][3] = t409 + t809 * t31 * L - t830 - t831;
  kl[10][4] = t465 + (-t418 - t834 - t421) * t3 / 0.2e1 + t809 * t42 * L + 0.4e1 * t815;
  kl[10][5] = -t476 + (-t431 + t843 + t434) * t3 / 0.2e1 + t809 * t52 * L - 0.4e1 * t800;
  kl[10][6] = -t799;
  kl[10][7] = t367 + (-t443 - t804 - t376) * t3 / 0.2e1 - t809 * t13 * L + t813;
  kl[10][8] = t387 + (-t451 - t819 - t396) * t3 / 0.2e1 - t809 * t23 * L + t826;
  kl[10][9] = t461 - t809 * t31 * L + t830 + t831;
  kl[10][10] = 0.16e2 * t391 + (-t467 - t834 - t468) * t3 / 0.2e1 + t809 * t61 * L + 0.2e1 * t815;
  kl[10][11] = -0.16e2 * t371 + (-t478 + t843 + t479) * t3 / 0.2e1 + t809 * t66 * L - 0.2e1 * t800;
  kl[11][0] = t882;
  kl[11][1] = t489 + (t491 + t887 - t498) * t3 / 0.2e1 + t892 * t13 * L - t896;
  kl[11][2] = t509 + (t511 + t902 - t518) * t3 / 0.2e1 + t892 * t23 * L - t909;
  kl[11][3] = t531 + t892 * t31 * L - t913 + t914;
  kl[11][4] = -t587 + (t540 - t917 + t543) * t3 / 0.2e1 + t892 * t42 * L + 0.4e1 * t898;
  kl[11][5] = t598 + (t553 + t926 - t556) * t3 / 0.2e1 + t892 * t52 * L - 0.4e1 * t883;
  kl[11][6] = -t882;
  kl[11][7] = -t489 + (t565 - t887 + t498) * t3 / 0.2e1 - t892 * t13 * L + t896;
  kl[11][8] = -t509 + (t573 - t902 + t518) * t3 / 0.2e1 - t892 * t23 * L + t909;
  kl[11][9] = t583 - t892 * t31 * L + t913 - t914;
  kl[11][10] = -0.16e2 * t513 + (t589 - t917 + t590) * t3 / 0.2e1 + t892 * t61 * L + 0.2e1 * t898;
  kl[11][11] = 0.16e2 * t493 + (t600 + t926 - t601) * t3 / 0.2e1 + t892 * t66 * L - 0.2e1 * t883;
}

// ================================ ConvNormal =============================
//
// This method computes the element normal force (N) for sections with
// conventional (i.e. diagonal) constitutive matrix [C].
//
//   L0 - element undeformed length                                    (in)
//   u  - nodal displacement vector                                    (in)
//
double cFrame3DCRTL :: ConvNormal(double L0, cVector &u)
{
  // Get material and geometric properties.
  
  cMatMec *mat = Section->GetMaterial( )->GetMecProp( );
  double matpar[2];
  mat->GetParam(matpar);
  double E  = matpar[0];
  double A  = Section->GetA( );
  double Ip = Section->GetJt( ); // ??

  // Compute [Bm] = [B0] + 1/2[BL].

  cMatrix B0(1, 12);
  B0Mat(L0, B0);

  cMatrix BL(1, 12);
  BLMat(L0, Ip, A, u, BL);

  cMatrix Bm(1, 12);
  Bm  = B0;
  Bm += 0.5*BL;

  // Compute the membrane strain => em = [Bm]{u}.

  double em = 0.0;
  for (int i = 0; i < 12; i++)
    em += Bm[0][i]*u[i];

  // Compute the normal force.

  double N = E*A*em;
  return(N);
}

// ================================ CoupNormal =============================
//
// This method computes the element normal force (N) for sections with
// fully coupled 4x4 constitutive matrix [C].
//
//   L  - element length                                               (in)
//   u  - nodal displacement vector                                    (in)
//
double cFrame3DCRTL :: CoupNormal(double L, cVector &u)
{
  // Geometric properties.
  
  double A  = Section->GetA( );
  double Ip = Section->GetJt( );   // ???

  // First line of the fully coupled section stiffness matrix.
  
  cMatrix Cb(4, 4);
  Section->SecStiffness3D(Cb);
  double C11 = Cb[0][0]; double C12 = Cb[0][1]; double C13 = Cb[0][2]; double C14 = Cb[0][3];

  // Auxiliary variables (Maple code).

  double t1 = L * L;
  double t3 = 0.1e1 / t1 / L;
  double t6 = 0.1e1 / t1;
  double t28 = 0.1e1 / L;
  double t30 = u[1] * t6;
  double t32 = u[5] * t28;
  double t34 = u[7] * t6;
  double t36 = u[11] * t28;
  double t38 = 0.3e1 / 0.5e1 * t30 + t32 / 0.20e2 - 0.3e1 / 0.5e1 * t34 + t36 / 0.20e2;
  double t40 = u[2] * t6;
  double t42 = u[4] * t28;
  double t44 = u[8] * t6;
  double t46 = u[10] * t28;
  double t48 = 0.3e1 / 0.5e1 * t40 - t42 / 0.20e2 - 0.3e1 / 0.5e1 * t44 - t46 / 0.20e2;
  double t52 = Ip / A;
  double t56 = u[3] * t6 * t52 - u[9] * t6 * t52;
  double t59 = u[2] * t28 / 0.20e2;
  double t62 = u[8] * t28 / 0.20e2;
  double t67 = u[1] * t28 / 0.20e2;
  double t70 = u[7] * t28 / 0.20e2;
  double t86 = -t28 * u[0] + t38 * u[1] + t48 * u[2] + t56 * u[3] / 0.2e1 + (-t59 + u[4] / 0.15e2 + t62 - u[10] / 0.60e2) * u[4] + (t67 + u[5] / 0.15e2 - t70 - u[11] / 0.60e2) * u[5] + t28 * u[6] - t38 * u[7] - t48 * u[8] - t56 * u[9] / 0.2e1 + (-t59 - u[4] / 0.60e2 + t62 + u[10] / 0.15e2) * u[10] + (t67 - u[5] / 0.60e2 - t70 + u[11] / 0.15e2) * u[11];
  
  // Compute the normal force.
  
  double N = (C12 * (0.12e2 * t3 * u[1] + 0.6e1 * t6 * u[5] - 0.12e2 * t3 * u[7] + 0.6e1 * t6 * u[11]) + C13 * (0.12e2 * t3 * u[2] - 0.6e1 * t6 * u[4] - 0.12e2 * t3 * u[8] - 0.6e1 * t6 * u[10])) * t1 / 0.2e1 + C11 * t86 * L + C12 * (-0.6e1 * t30 - 0.4e1 * t32 + 0.6e1 * t34 - 0.2e1 * t36) * L + C13 * (-0.6e1 * t40 + 0.4e1 * t42 + 0.6e1 * t44 + 0.2e1 * t46) * L + C14 * (-t28 * u[3] + t28 * u[9]) * L;

  return(N);
}

// =============================== ConvIntForce ============================
//
// This method computes the element internal force obtained from the
// integration {gl} = int([Bbar]t{sig}), where {sig} = [C]{u} and the
// section constitutive matrix [C] is diagonal.
//
//   L  - element length                                               (in)
//   u  - nodal displacement vector                                    (in)
//   gl - internal force vector (12)                                  (out)
//
void cFrame3DCRTL :: ConvIntForce(double L, cVector &u, cVector &gl)
{
  // Material properties.
  
  cMatMec *mat = Section->GetMaterial( )->GetMecProp( );
  double matpar[2];
  mat->GetParam(matpar);
  double E  = matpar[0];
  double nu = matpar[1];
  double G  = 0.5*E/(1 + nu);

  // Section properties.
  
  double A   = Section->GetA( );
  double Ip  = Section->GetJt( );  // ??
  double EA  = E*A;
  double EIz = E*Section->GetIz( );
  double EIy = E*Section->GetIy( );
  double GJ  = G*Ip;

  // Auxiliary variables (Maple code).

  double t1 = 0.1e1 / L;
  double t2 = t1 * EA;
  double t4 = L * L;
  double t5 = 0.1e1 / t4;
  double t6 = u[1] * t5;
  double t8 = u[5] * t1;
  double t10 = u[7] * t5;
  double t12 = u[11] * t1;
  double t14 = 0.3e1 / 0.5e1 * t6 + t8 / 0.20e2 - 0.3e1 / 0.5e1 * t10 + t12 / 0.20e2;
  double t17 = u[2] * t5;
  double t19 = u[4] * t1;
  double t21 = u[8] * t5;
  double t23 = u[10] * t1;
  double t25 = 0.3e1 / 0.5e1 * t17 - t19 / 0.20e2 - 0.3e1 / 0.5e1 * t21 - t23 / 0.20e2;
  double t30 = Ip / A;
  double t34 = u[3] * t5 * t30 - u[9] * t5 * t30;
  double t37 = u[2] * t1;
  double t38 = t37 / 0.20e2;
  double t40 = u[8] * t1;
  double t41 = t40 / 0.20e2;
  double t43 = -t38 + u[4] / 0.15e2 + t41 - u[10] / 0.60e2;
  double t46 = u[1] * t1;
  double t47 = t46 / 0.20e2;
  double t49 = u[7] * t1;
  double t50 = t49 / 0.20e2;
  double t52 = t47 + u[5] / 0.15e2 - t50 - u[11] / 0.60e2;
  double t64 = -t38 - u[4] / 0.60e2 + t41 + u[10] / 0.15e2;
  double t69 = t47 - u[5] / 0.60e2 - t50 + u[11] / 0.15e2;
  double t72 = t2 * u[0] - EA * t14 * u[1] - EA * t25 * u[2] - EA * t34 * u[3] / 0.2e1 - EA * t43 * u[4] - EA * t52 * u[5] - t2 * u[6] + EA * t14 * u[7] + EA * t25 * u[8] + EA * t34 * u[9] / 0.2e1 - EA * t64 * u[10] - EA * t69 * u[11];
  double t73 = t4 * t4;
  double t75 = 0.1e1 / t73 / L;
  double t76 = t75 * EIz;
  double t80 = 0.1e1 / t73 / t4;
  double t81 = t80 * EIz;
  double t88 = 0.72e2 * t76 * u[5] + 0.144e3 * t81 * u[1] - 0.144e3 * t81 * u[7] + 0.72e2 * t76 * u[11];
  double t89 = t4 * L;
  double t92 = 0.1e1 / t73;
  double t93 = t92 * EIz;
  double t94 = t93 * u[5];
  double t96 = t76 * u[7];
  double t98 = t76 * u[1];
  double t100 = t93 * u[11];
  double t102 = -0.84e2 * t94 + 0.144e3 * t96 - 0.144e3 * t98 - 0.60e2 * t100;
  double t109 = 0.6e1 / 0.5e1 * t6 + t8 / 0.10e2 - 0.6e1 / 0.5e1 * t10 + t12 / 0.10e2;
  double t110 = EA * t109;
  double t113 = 0.36e2 * t93;
  double t118 = -t34 * u[9] * L / 0.2e1;
  double t121 = t64 * u[10] * L;
  double t124 = t25 * u[2] * L;
  double t127 = t34 * u[3] * L / 0.2e1;
  double t130 = t43 * u[4] * L;
  double t133 = 0.1e1 / t89;
  double t134 = t133 * EIz;
  double t135 = 0.24e2 * t134;
  double t140 = 0.12e2 * t134;
  double t150 = -t25 * u[8] * L;
  double t152 = t88 * t89 / 0.3e1 + t102 * t4 / 0.2e1 - t110 * u[0] + (t110 * t14 + t113) * u[1] * L + t110 * t118 + t110 * t121 + t110 * t124 + t110 * t127 + t110 * t130 + (t110 * t52 + t135) * u[5] * L + (t110 * t69 + t140) * u[11] * L + t110 * u[6] + (-t110 * t14 - t113) * u[7] * L + t110 * t150;
  double t153 = t80 * EIy;
  double t156 = t75 * EIy;
  double t163 = -0.144e3 * t153 * u[8] - 0.72e2 * t156 * u[4] + 0.144e3 * t153 * u[2] - 0.72e2 * t156 * u[10];
  double t166 = t92 * EIy;
  double t167 = t166 * u[10];
  double t169 = t156 * u[2];
  double t171 = t156 * u[8];
  double t173 = t166 * u[4];
  double t175 = 0.60e2 * t167 - 0.144e3 * t169 + 0.144e3 * t171 + 0.84e2 * t173;
  double t182 = 0.6e1 / 0.5e1 * t17 - t19 / 0.10e2 - 0.6e1 / 0.5e1 * t21 - t23 / 0.10e2;
  double t183 = EA * t182;
  double t186 = t14 * u[1] * L;
  double t189 = 0.36e2 * t166;
  double t194 = t133 * EIy;
  double t195 = 0.12e2 * t194;
  double t205 = 0.24e2 * t194;
  double t210 = -t14 * u[7] * L;
  double t213 = t52 * u[5] * L;
  double t218 = t69 * u[11] * L;
  double t220 = t163 * t89 / 0.3e1 + t175 * t4 / 0.2e1 - t183 * u[0] + t183 * t186 + (t183 * t25 + t189) * u[2] * L + (t183 * t64 - t195) * u[10] * L + (-t183 * t25 - t189) * u[8] * L + t183 * t127 + (t183 * t43 - t205) * u[4] * L + t183 * t210 + t183 * t213 + t183 * u[6] + t183 * t118 + t183 * t218;
  double t221 = EA * t34;
  double t226 = t5 * GJ;
  double t241 = -t221 * u[0] + t221 * t186 + t221 * t124 + (t221 * t34 / 0.2e1 + t226) * u[3] * L + t221 * t130 + t221 * t213 + t221 * u[6] + t221 * t210 + t221 * t150 + (-t221 * t34 / 0.2e1 - t226) * u[9] * L + t221 * t121 + t221 * t218;
  double t248 = (0.72e2 * t171 + 0.36e2 * t173 - 0.72e2 * t169 + 0.36e2 * t167) * t89 / 0.3e1;
  double t249 = t194 * u[10];
  double t251 = t166 * u[2];
  double t253 = t166 * u[8];
  double t255 = t194 * u[4];
  double t260 = t37 / 0.10e2;
  double t262 = t40 / 0.10e2;
  double t265 = EA * (-t260 + 0.2e1 / 0.15e2 * u[4] + t262 - u[10] / 0.30e2);
  double t273 = t5 * EIy;
  double t274 = 0.8e1 * t273;
  double t293 = t248 + (-0.36e2 * t249 + 0.84e2 * t251 - 0.84e2 * t253 - 0.48e2 * t255) * t4 / 0.2e1 - t265 * u[0] + t265 * t186 + (t265 * t25 - t205) * u[2] * L + (t265 * t64 + t274) * u[10] * L + (-t265 * t25 + t205) * u[8] * L + t265 * t127 + (t265 * t43 + 0.16e2 * t273) * u[4] * L + t265 * t210 + t265 * t213 + t265 * u[6] + t265 * t118 + t265 * t218;
  double t300 = (0.36e2 * t94 + 0.72e2 * t98 - 0.72e2 * t96 + 0.36e2 * t100) * t89 / 0.3e1;
  double t301 = t134 * u[5];
  double t303 = t93 * u[7];
  double t305 = t93 * u[1];
  double t307 = t134 * u[11];
  double t312 = t46 / 0.10e2;
  double t314 = t49 / 0.10e2;
  double t317 = EA * (t312 + 0.2e1 / 0.15e2 * u[5] - t314 - u[11] / 0.30e2);
  double t329 = t5 * EIz;
  double t335 = 0.8e1 * t329;
  double t345 = t300 + (-0.48e2 * t301 + 0.84e2 * t303 - 0.84e2 * t305 - 0.36e2 * t307) * t4 / 0.2e1 - t317 * u[0] + (t317 * t14 + t135) * u[1] * L + t317 * t118 + t317 * t121 + t317 * t124 + t317 * t127 + t317 * t130 + (t317 * t52 + 0.16e2 * t329) * u[5] * L + (t317 * t69 + t335) * u[11] * L + t317 * u[6] + (-t317 * t14 - t135) * u[7] * L + t317 * t150;
  double t350 = -EA * t109;
  double t375 = -t88 * t89 / 0.3e1 - t102 * t4 / 0.2e1 - t350 * u[0] + (t350 * t14 - t113) * u[1] * L + t350 * t118 + t350 * t121 + t350 * t124 + t350 * t127 + t350 * t130 + (t350 * t52 - t135) * u[5] * L + (t350 * t69 - t140) * u[11] * L + t350 * u[6] + (-t350 * t14 + t113) * u[7] * L + t350 * t150;
  double t380 = -EA * t182;
  double t405 = -t163 * t89 / 0.3e1 - t175 * t4 / 0.2e1 - t380 * u[0] + t380 * t186 + (t380 * t25 - t189) * u[2] * L + (t380 * t64 + t195) * u[10] * L + (-t380 * t25 + t189) * u[8] * L + t380 * t127 + (t380 * t43 + t205) * u[4] * L + t380 * t210 + t380 * t213 + t380 * u[6] + t380 * t118 + t380 * t218;
  double t406 = -EA * t34;
  double t425 = -t406 * u[0] + t406 * t186 + t406 * t124 + (t406 * t34 / 0.2e1 - t226) * u[3] * L + t406 * t130 + t406 * t213 + t406 * u[6] + t406 * t210 + t406 * t150 + (-t406 * t34 / 0.2e1 + t226) * u[9] * L + t406 * t121 + t406 * t218;
  double t436 = EA * (-t260 - u[4] / 0.30e2 + t262 + 0.2e1 / 0.15e2 * u[10]);
  double t462 = t248 + (-0.24e2 * t249 + 0.60e2 * t251 - 0.60e2 * t253 - 0.36e2 * t255) * t4 / 0.2e1 - t436 * u[0] + t436 * t186 + (t436 * t25 - t195) * u[2] * L + (t436 * t64 + 0.4e1 * t273) * u[10] * L + (-t436 * t25 + t195) * u[8] * L + t436 * t127 + (t436 * t43 + t274) * u[4] * L + t436 * t210 + t436 * t213 + t436 * u[6] + t436 * t118 + t436 * t218;
  double t473 = EA * (t312 - u[5] / 0.30e2 - t314 + 0.2e1 / 0.15e2 * u[11]);
  double t499 = t300 + (-0.36e2 * t301 + 0.60e2 * t303 - 0.60e2 * t305 - 0.24e2 * t307) * t4 / 0.2e1 - t473 * u[0] + (t473 * t14 + t140) * u[1] * L + t473 * t118 + t473 * t121 + t473 * t124 + t473 * t127 + t473 * t130 + (t473 * t52 + t335) * u[5] * L + (t473 * t69 + 0.4e1 * t329) * u[11] * L + t473 * u[6] + (-t473 * t14 - t140) * u[7] * L + t473 * t150;

  // Internal force vector.
  
  gl[0] = t72;
  gl[1] = t152;
  gl[2] = t220;
  gl[3] = t241;
  gl[4] = t293;
  gl[5] = t345;
  gl[6] = -t72;
  gl[7] = t375;
  gl[8] = t405;
  gl[9] = t425;
  gl[10] = t462;
  gl[11] = t499;
}

// =============================== CoupIntForce ============================
//
// This method computes the element internal force obtained from the
// integration {gl} = int([Bbar]t{sig}), where {sig} = [C]{u} and the
// section constitutive matrix [C] is fully coupled.
//
//   L  - element length                                               (in)
//   u  - nodal displacement vector                                    (in)
//   gl - internal force vector (12)                                  (out)
//
void cFrame3DCRTL :: CoupIntForce(double L, cVector &u, cVector &gl)
{
  // Geometric properties.

  double A  = Section->GetA( );
  double Ip = Section->GetJt( );

  // Section stiffness - fully coupled 4x4 matrix.
  
  cMatrix Cb(4, 4);
  Cb.Zero( );
  double C11 = Cb[0][0]; double C12 = Cb[0][1]; double C13 = Cb[0][2]; double C14 = Cb[0][3];
  double C21 = Cb[1][0]; double C22 = Cb[1][1]; double C23 = Cb[1][2]; double C24 = Cb[1][3];
  double C31 = Cb[2][0]; double C32 = Cb[2][1]; double C33 = Cb[2][2]; double C34 = Cb[2][3];
  double C41 = Cb[3][0]; double C42 = Cb[3][1]; double C43 = Cb[3][2]; double C44 = Cb[3][3];

  // Auxiliary variables (Maple code).

  double t1 = L * L;
  double t2 = t1 * L;
  double t3 = 0.1e1 / t2;
  double t4 = t3 * C13;
  double t7 = t1 * t1;
  double t8 = 0.1e1 / t7;
  double t9 = t8 * C13;
  double t12 = t3 * C12;
  double t15 = t8 * C12;
  double t26 = 0.6e1 * t4 * u[4] - 0.12e2 * t9 * u[2] - 0.6e1 * t12 * u[5] - 0.12e2 * t15 * u[1] + 0.12e2 * t15 * u[7] + 0.6e1 * t4 * u[10] + 0.12e2 * t9 * u[8] - 0.6e1 * t12 * u[11];
  double t29 = 0.1e1 / L;
  double t30 = t29 * C11;
  double t31 = t30 * u[0];
  double t32 = 0.1e1 / t1;
  double t33 = u[1] * t32;
  double t35 = u[5] * t29;
  double t37 = u[7] * t32;
  double t39 = u[11] * t29;
  double t41 = 0.3e1 / 0.5e1 * t33 + t35 / 0.20e2 - 0.3e1 / 0.5e1 * t37 + t39 / 0.20e2;
  double t43 = 0.6e1 * t12;
  double t44 = -t30 * t41 + t43;
  double t49 = Ip / A;
  double t53 = -u[3] * t32 * t49 + u[9] * t32 * t49;
  double t55 = t32 * C14;
  double t56 = -t30 * t53 / 0.2e1 - t55;
  double t59 = u[2] * t32;
  double t61 = u[4] * t29;
  double t63 = u[8] * t32;
  double t65 = u[10] * t29;
  double t67 = 0.3e1 / 0.5e1 * t59 - t61 / 0.20e2 - 0.3e1 / 0.5e1 * t63 - t65 / 0.20e2;
  double t69 = 0.6e1 * t4;
  double t70 = -t30 * t67 + t69;
  double t73 = u[1] * t29;
  double t74 = t73 / 0.20e2;
  double t76 = u[7] * t29;
  double t77 = t76 / 0.20e2;
  double t79 = t74 + u[5] / 0.15e2 - t77 - u[11] / 0.60e2;
  double t81 = t32 * C12;
  double t83 = -t30 * t79 + 0.4e1 * t81;
  double t87 = t30 * t53 / 0.2e1 + t55;
  double t90 = u[2] * t29;
  double t91 = t90 / 0.20e2;
  double t93 = u[8] * t29;
  double t94 = t93 / 0.20e2;
  double t96 = -t91 + u[4] / 0.15e2 + t94 - u[10] / 0.60e2;
  double t98 = t32 * C13;
  double t100 = -t30 * t96 - 0.4e1 * t98;
  double t104 = t30 * t41 - t43;
  double t109 = -t91 - u[4] / 0.60e2 + t94 + u[10] / 0.15e2;
  double t112 = -t30 * t109 - 0.2e1 * t98;
  double t116 = t30 * t67 - t69;
  double t119 = t30 * u[6];
  double t122 = t74 - u[5] / 0.60e2 - t77 + u[11] / 0.15e2;
  double t125 = -t30 * t122 + 0.2e1 * t81;
  double t128 = t26 * t1 / 0.2e1 + t31 + t44 * u[1] * L + t56 * u[9] * L + t70 * u[2] * L + t83 * u[5] * L + t87 * u[3] * L + t100 * u[4] * L + t104 * u[7] * L + t112 * u[10] * L + t116 * u[8] * L - t119 + t125 * u[11] * L;
  double t130 = 0.1e1 / t7 / t1;
  double t131 = t130 * C23;
  double t135 = 0.1e1 / t7 / L;
  double t136 = t135 * C22;
  double t139 = t130 * C22;
  double t144 = t135 * C23;
  double t153 = 0.144e3 * t131 * u[2] + 0.72e2 * t136 * u[5] - 0.144e3 * t139 * u[7] - 0.144e3 * t131 * u[8] - 0.72e2 * t144 * u[10] - 0.72e2 * t144 * u[4] + 0.144e3 * t139 * u[1] + 0.72e2 * t136 * u[11];
  double t156 = t3 * C21;
  double t158 = -0.12e2 * t156 * t41;
  double t163 = 0.6e1 / 0.5e1 * t33 + t35 / 0.10e2 - 0.6e1 / 0.5e1 * t37 + t39 / 0.10e2;
  double t166 = 0.6e1 * t32 * C22;
  double t167 = t163 * C12 - t166;
  double t169 = 0.12e2 * t167 * t3;
  double t170 = 0.72e2 * t136;
  double t173 = t8 * C21;
  double t175 = 0.12e2 * t173 * u[6];
  double t177 = 0.12e2 * t173 * u[0];
  double t179 = 0.12e2 * t156 * t109;
  double t182 = 0.6e1 * t32 * C23;
  double t183 = t163 * C13 - t182;
  double t185 = 0.6e1 * t183 * t32;
  double t186 = t8 * C23;
  double t187 = 0.24e2 * t186;
  double t191 = t8 * C24;
  double t192 = t156 * t53 / 0.2e1 + t191;
  double t195 = -t156 * t53 / 0.2e1 - t191;
  double t198 = 0.12e2 * t156 * t67;
  double t200 = 0.12e2 * t183 * t3;
  double t201 = 0.72e2 * t144;
  double t205 = 0.12e2 * t156 * t96;
  double t206 = 0.48e2 * t186;
  double t210 = 0.12e2 * t156 * t122;
  double t212 = 0.6e1 * t167 * t32;
  double t213 = t8 * C22;
  double t214 = 0.24e2 * t213;
  double t218 = 0.12e2 * t156 * t79;
  double t219 = 0.48e2 * t213;
  double t223 = 0.12e2 * t156 * t41;
  double t227 = -0.12e2 * t156 * t67;
  double t230 = (t158 - t169 + t170) * u[7] + t175 - t177 + (t179 - t185 + t187) * u[10] + 0.12e2 * t192 * u[9] + 0.12e2 * t195 * u[3] + (t198 + t200 - t201) * u[2] + (t205 - t185 + t206) * u[4] + (t210 + t212 - t214) * u[11] + (t218 + t212 - t219) * u[5] + (t223 + t169 - t170) * u[1] + (t227 - t200 + t201) * u[8];
  double t234 = t32 * C21;
  double t235 = 0.6e1 * t234;
  double t236 = t163 * C11 - t235;
  double t250 = 0.6e1 * t32 * C24;
  double t252 = (t163 * C14 - t250) * t29;
  double t257 = t167 * t29;
  double t263 = t183 * t29;
  double t290 = t153 * t2 / 0.3e1 + t230 * t1 / 0.2e1 - t236 * u[0] + (t236 * t67 - t185) * u[2] * L + (t236 * t41 - t212) * u[1] * L + t236 * u[6] + (-t236 * t53 / 0.2e1 - t252) * u[3] * L + (t236 * t79 - 0.4e1 * t257) * u[5] * L + (t236 * t96 + 0.4e1 * t263) * u[4] * L + (t236 * t53 / 0.2e1 + t252) * u[9] * L + (t236 * t122 - 0.2e1 * t257) * u[11] * L + (-t236 * t67 + t185) * u[8] * L + (-t236 * t41 + t212) * u[7] * L + (t236 * t109 + 0.2e1 * t263) * u[10] * L;
  double t291 = t130 * C33;
  double t294 = t135 * C32;
  double t297 = t130 * C32;
  double t302 = t135 * C33;
  double t311 = 0.144e3 * t291 * u[2] + 0.72e2 * t294 * u[5] - 0.144e3 * t297 * u[7] - 0.144e3 * t291 * u[8] - 0.72e2 * t302 * u[10] - 0.72e2 * t302 * u[4] + 0.144e3 * t297 * u[1] + 0.72e2 * t294 * u[11];
  double t314 = t3 * C31;
  double t316 = -0.12e2 * t314 * t41;
  double t321 = 0.6e1 / 0.5e1 * t59 - t61 / 0.10e2 - 0.6e1 / 0.5e1 * t63 - t65 / 0.10e2;
  double t324 = 0.6e1 * t32 * C32;
  double t325 = t321 * C12 - t324;
  double t327 = 0.12e2 * t325 * t3;
  double t328 = 0.72e2 * t294;
  double t331 = t8 * C31;
  double t333 = 0.12e2 * t331 * u[6];
  double t335 = 0.12e2 * t331 * u[0];
  double t337 = 0.12e2 * t314 * t109;
  double t340 = 0.6e1 * t32 * C33;
  double t341 = t321 * C13 - t340;
  double t343 = 0.6e1 * t341 * t32;
  double t344 = t8 * C33;
  double t345 = 0.24e2 * t344;
  double t349 = t8 * C34;
  double t350 = t314 * t53 / 0.2e1 + t349;
  double t353 = -t314 * t53 / 0.2e1 - t349;
  double t356 = 0.12e2 * t314 * t67;
  double t358 = 0.12e2 * t341 * t3;
  double t359 = 0.72e2 * t302;
  double t363 = 0.12e2 * t314 * t96;
  double t364 = 0.48e2 * t344;
  double t368 = 0.12e2 * t314 * t122;
  double t370 = 0.6e1 * t325 * t32;
  double t371 = t8 * C32;
  double t372 = 0.24e2 * t371;
  double t376 = 0.12e2 * t314 * t79;
  double t377 = 0.48e2 * t371;
  double t381 = 0.12e2 * t314 * t41;
  double t385 = -0.12e2 * t314 * t67;
  double t388 = (t316 - t327 + t328) * u[7] + t333 - t335 + (t337 - t343 + t345) * u[10] + 0.12e2 * t350 * u[9] + 0.12e2 * t353 * u[3] + (t356 + t358 - t359) * u[2] + (t363 - t343 + t364) * u[4] + (t368 + t370 - t372) * u[11] + (t376 + t370 - t377) * u[5] + (t381 + t327 - t328) * u[1] + (t385 - t358 + t359) * u[8];
  double t392 = t32 * C31;
  double t393 = 0.6e1 * t392;
  double t394 = t321 * C11 - t393;
  double t408 = 0.6e1 * t32 * C34;
  double t410 = (t321 * C14 - t408) * t29;
  double t415 = t325 * t29;
  double t421 = t341 * t29;
  double t448 = t311 * t2 / 0.3e1 + t388 * t1 / 0.2e1 - t394 * u[0] + (t394 * t67 - t343) * u[2] * L + (t394 * t41 - t370) * u[1] * L + t394 * u[6] + (-t394 * t53 / 0.2e1 - t410) * u[3] * L + (t394 * t79 - 0.4e1 * t415) * u[5] * L + (t394 * t96 + 0.4e1 * t421) * u[4] * L + (t394 * t53 / 0.2e1 + t410) * u[9] * L + (t394 * t122 - 0.2e1 * t415) * u[11] * L + (-t394 * t67 + t343) * u[8] * L + (-t394 * t41 + t370) * u[7] * L + (t394 * t109 + 0.2e1 * t421) * u[10] * L;
  double t450 = t29 * C43;
  double t451 = -t53 * C13 - t450;
  double t452 = t451 * t32;
  double t455 = t451 * t3;
  double t459 = t29 * C42;
  double t460 = -t53 * C12 - t459;
  double t461 = t460 * t32;
  double t464 = t460 * t3;
  double t479 = t29 * C41;
  double t480 = -t53 * C11 - t479;
  double t483 = 0.6e1 * t461;
  double t489 = t29 * C44;
  double t491 = (-t53 * C14 - t489) * t29;
  double t496 = 0.6e1 * t452;
  double t501 = t460 * t29;
  double t511 = t451 * t29;
  double t535 = (-0.6e1 * t452 * u[4] + 0.12e2 * t455 * u[2] + 0.6e1 * t461 * u[5] + 0.12e2 * t464 * u[1] - 0.12e2 * t464 * u[7] - 0.6e1 * t452 * u[10] - 0.12e2 * t455 * u[8] + 0.6e1 * t461 * u[11]) * t1 / 0.2e1 - t480 * u[0] + (t480 * t41 - t483) * u[1] * L + (t480 * t53 / 0.2e1 + t491) * u[9] * L + (t480 * t67 - t496) * u[2] * L + (t480 * t79 - 0.4e1 * t501) * u[5] * L + (-t480 * t53 / 0.2e1 - t491) * u[3] * L + (t480 * t96 + 0.4e1 * t511) * u[4] * L + (-t480 * t41 + t483) * u[7] * L + (t480 * t109 + 0.2e1 * t511) * u[10] * L + (-t480 * t67 + t496) * u[8] * L + t480 * u[6] + (t480 * t122 - 0.2e1 * t501) * u[11] * L;
  double t554 = (-0.72e2 * t302 * u[2] - 0.36e2 * t371 * u[5] + 0.72e2 * t294 * u[7] + 0.72e2 * t302 * u[8] + 0.36e2 * t344 * u[10] + 0.36e2 * t344 * u[4] - 0.72e2 * t294 * u[1] - 0.36e2 * t371 * u[11]) * t2 / 0.3e1;
  double t556 = -0.6e1 * t392 * t41;
  double t557 = t90 / 0.10e2;
  double t559 = t93 / 0.10e2;
  double t561 = -t557 + 0.2e1 / 0.15e2 * u[4] + t559 - u[10] / 0.30e2;
  double t563 = t29 * C32;
  double t565 = t561 * C12 + 0.4e1 * t563;
  double t567 = 0.12e2 * t565 * t3;
  double t568 = 0.36e2 * t371;
  double t572 = 0.6e1 * t314 * u[6];
  double t574 = 0.6e1 * t314 * u[0];
  double t576 = 0.6e1 * t392 * t109;
  double t578 = t29 * C33;
  double t580 = t561 * C13 + 0.4e1 * t578;
  double t582 = 0.6e1 * t580 * t32;
  double t583 = t3 * C33;
  double t584 = 0.12e2 * t583;
  double t588 = t3 * C34;
  double t590 = 0.6e1 * (-t392 * t53 / 0.2e1 - t588) * u[9];
  double t593 = 0.6e1 * (t392 * t53 / 0.2e1 + t588) * u[3];
  double t595 = 0.6e1 * t392 * t67;
  double t597 = 0.12e2 * t580 * t3;
  double t598 = 0.36e2 * t344;
  double t602 = 0.6e1 * t392 * t96;
  double t603 = 0.24e2 * t583;
  double t607 = 0.6e1 * t392 * t122;
  double t609 = 0.6e1 * t565 * t32;
  double t610 = t3 * C32;
  double t611 = 0.12e2 * t610;
  double t615 = 0.6e1 * t392 * t79;
  double t616 = 0.24e2 * t610;
  double t620 = 0.6e1 * t392 * t41;
  double t624 = -0.6e1 * t392 * t67;
  double t627 = (-t556 - t567 - t568) * u[7] - t572 + t574 + (-t576 - t582 - t584) * u[10] + t590 + t593 + (-t595 + t597 + t598) * u[2] + (-t602 - t582 - t603) * u[4] + (-t607 + t609 + t611) * u[11] + (-t615 + t609 + t616) * u[5] + (-t620 + t567 + t568) * u[1] + (-t624 - t597 - t598) * u[8];
  double t631 = t29 * C31;
  double t633 = t561 * C11 + 0.4e1 * t631;
  double t646 = t29 * C34;
  double t649 = (t561 * C14 + 0.4e1 * t646) * t29;
  double t654 = t565 * t29;
  double t660 = t580 * t29;
  double t687 = t554 + t627 * t1 / 0.2e1 - t633 * u[0] + (t633 * t67 - t582) * u[2] * L + (t633 * t41 - t609) * u[1] * L + t633 * u[6] + (-t633 * t53 / 0.2e1 - t649) * u[3] * L + (t633 * t79 - 0.4e1 * t654) * u[5] * L + (t633 * t96 + 0.4e1 * t660) * u[4] * L + (t633 * t53 / 0.2e1 + t649) * u[9] * L + (t633 * t122 - 0.2e1 * t654) * u[11] * L + (-t633 * t67 + t582) * u[8] * L + (-t633 * t41 + t609) * u[7] * L + (t633 * t109 + 0.2e1 * t660) * u[10] * L;
  double t706 = (0.72e2 * t144 * u[2] + 0.36e2 * t213 * u[5] - 0.72e2 * t136 * u[7] - 0.72e2 * t144 * u[8] - 0.36e2 * t186 * u[10] - 0.36e2 * t186 * u[4] + 0.72e2 * t136 * u[1] + 0.36e2 * t213 * u[11]) * t2 / 0.3e1;
  double t708 = -0.6e1 * t234 * t41;
  double t709 = t73 / 0.10e2;
  double t711 = t76 / 0.10e2;
  double t713 = t709 + 0.2e1 / 0.15e2 * u[5] - t711 - u[11] / 0.30e2;
  double t715 = t29 * C22;
  double t717 = t713 * C12 - 0.4e1 * t715;
  double t719 = 0.12e2 * t717 * t3;
  double t720 = 0.36e2 * t213;
  double t724 = 0.6e1 * t156 * u[6];
  double t726 = 0.6e1 * t156 * u[0];
  double t728 = 0.6e1 * t234 * t109;
  double t730 = t29 * C23;
  double t732 = t713 * C13 - 0.4e1 * t730;
  double t734 = 0.6e1 * t732 * t32;
  double t735 = t3 * C23;
  double t736 = 0.12e2 * t735;
  double t740 = t3 * C24;
  double t742 = 0.6e1 * (t234 * t53 / 0.2e1 + t740) * u[9];
  double t745 = 0.6e1 * (-t234 * t53 / 0.2e1 - t740) * u[3];
  double t747 = 0.6e1 * t234 * t67;
  double t749 = 0.12e2 * t732 * t3;
  double t750 = 0.36e2 * t186;
  double t754 = 0.6e1 * t234 * t96;
  double t755 = 0.24e2 * t735;
  double t759 = 0.6e1 * t234 * t122;
  double t761 = 0.6e1 * t717 * t32;
  double t762 = t3 * C22;
  double t763 = 0.12e2 * t762;
  double t767 = 0.6e1 * t234 * t79;
  double t768 = 0.24e2 * t762;
  double t772 = 0.6e1 * t234 * t41;
  double t776 = -0.6e1 * t234 * t67;
  double t779 = (t708 - t719 + t720) * u[7] + t724 - t726 + (t728 - t734 + t736) * u[10] + t742 + t745 + (t747 + t749 - t750) * u[2] + (t754 - t734 + t755) * u[4] + (t759 + t761 - t763) * u[11] + (t767 + t761 - t768) * u[5] + (t772 + t719 - t720) * u[1] + (t776 - t749 + t750) * u[8];
  double t783 = t29 * C21;
  double t785 = t713 * C11 - 0.4e1 * t783;
  double t798 = t29 * C24;
  double t801 = (t713 * C14 - 0.4e1 * t798) * t29;
  double t806 = t717 * t29;
  double t812 = t732 * t29;
  double t839 = t706 + t779 * t1 / 0.2e1 - t785 * u[0] + (t785 * t67 - t734) * u[2] * L + (t785 * t41 - t761) * u[1] * L + t785 * u[6] + (-t785 * t53 / 0.2e1 - t801) * u[3] * L + (t785 * t79 - 0.4e1 * t806) * u[5] * L + (t785 * t96 + 0.4e1 * t812) * u[4] * L + (t785 * t53 / 0.2e1 + t801) * u[9] * L + (t785 * t122 - 0.2e1 * t806) * u[11] * L + (-t785 * t67 + t734) * u[8] * L + (-t785 * t41 + t761) * u[7] * L + (t785 * t109 + 0.2e1 * t812) * u[10] * L;
  double t862 = -t26 * t1 / 0.2e1 - t31 - t44 * u[1] * L - t56 * u[9] * L - t70 * u[2] * L - t83 * u[5] * L - t87 * u[3] * L - t100 * u[4] * L - t104 * u[7] * L - t112 * u[10] * L - t116 * u[8] * L + t119 - t125 * u[11] * L;
  double t866 = -t163 * C12 + t166;
  double t868 = 0.12e2 * t866 * t3;
  double t872 = -t163 * C13 + t182;
  double t874 = 0.6e1 * t872 * t32;
  double t880 = 0.12e2 * t872 * t3;
  double t886 = 0.6e1 * t866 * t32;
  double t895 = (-t158 - t868 - t170) * u[7] - t175 + t177 + (-t179 - t874 - t187) * u[10] - 0.12e2 * t192 * u[9] - 0.12e2 * t195 * u[3] + (-t198 + t880 + t201) * u[2] + (-t205 - t874 - t206) * u[4] + (-t210 + t886 + t214) * u[11] + (-t218 + t886 + t219) * u[5] + (-t223 + t868 + t170) * u[1] + (-t227 - t880 - t201) * u[8];
  double t899 = -t163 * C11 + t235;
  double t913 = (-t163 * C14 + t250) * t29;
  double t918 = t866 * t29;
  double t924 = t872 * t29;
  double t951 = -t153 * t2 / 0.3e1 + t895 * t1 / 0.2e1 - t899 * u[0] + (t899 * t67 - t874) * u[2] * L + (t899 * t41 - t886) * u[1] * L + t899 * u[6] + (-t899 * t53 / 0.2e1 - t913) * u[3] * L + (t899 * t79 - 0.4e1 * t918) * u[5] * L + (t899 * t96 + 0.4e1 * t924) * u[4] * L + (t899 * t53 / 0.2e1 + t913) * u[9] * L + (t899 * t122 - 0.2e1 * t918) * u[11] * L + (-t899 * t67 + t874) * u[8] * L + (-t899 * t41 + t886) * u[7] * L + (t899 * t109 + 0.2e1 * t924) * u[10] * L;
  double t955 = -t321 * C12 + t324;
  double t957 = 0.12e2 * t955 * t3;
  double t961 = -t321 * C13 + t340;
  double t963 = 0.6e1 * t961 * t32;
  double t969 = 0.12e2 * t961 * t3;
  double t975 = 0.6e1 * t955 * t32;
  double t984 = (-t316 - t957 - t328) * u[7] - t333 + t335 + (-t337 - t963 - t345) * u[10] - 0.12e2 * t350 * u[9] - 0.12e2 * t353 * u[3] + (-t356 + t969 + t359) * u[2] + (-t363 - t963 - t364) * u[4] + (-t368 + t975 + t372) * u[11] + (-t376 + t975 + t377) * u[5] + (-t381 + t957 + t328) * u[1] + (-t385 - t969 - t359) * u[8];
  double t988 = -t321 * C11 + t393;
  double t1002 = (-t321 * C14 + t408) * t29;
  double t1007 = t955 * t29;
  double t1013 = t961 * t29;
  double t1040 = -t311 * t2 / 0.3e1 + t984 * t1 / 0.2e1 - t988 * u[0] + (t988 * t67 - t963) * u[2] * L + (t988 * t41 - t975) * u[1] * L + t988 * u[6] + (-t988 * t53 / 0.2e1 - t1002) * u[3] * L + (t988 * t79 - 0.4e1 * t1007) * u[5] * L + (t988 * t96 + 0.4e1 * t1013) * u[4] * L + (t988 * t53 / 0.2e1 + t1002) * u[9] * L + (t988 * t122 - 0.2e1 * t1007) * u[11] * L + (-t988 * t67 + t963) * u[8] * L + (-t988 * t41 + t975) * u[7] * L + (t988 * t109 + 0.2e1 * t1013) * u[10] * L;
  double t1042 = t53 * C13 + t450;
  double t1043 = t1042 * t32;
  double t1046 = t1042 * t3;
  double t1050 = t53 * C12 + t459;
  double t1051 = t1050 * t32;
  double t1054 = t1050 * t3;
  double t1069 = t53 * C11 + t479;
  double t1072 = 0.6e1 * t1051;
  double t1079 = (t53 * C14 + t489) * t29;
  double t1084 = 0.6e1 * t1043;
  double t1089 = t1050 * t29;
  double t1099 = t1042 * t29;
  double t1123 = (-0.6e1 * t1043 * u[4] + 0.12e2 * t1046 * u[2] + 0.6e1 * t1051 * u[5] + 0.12e2 * t1054 * u[1] - 0.12e2 * t1054 * u[7] - 0.6e1 * t1043 * u[10] - 0.12e2 * t1046 * u[8] + 0.6e1 * t1051 * u[11]) * t1 / 0.2e1 - t1069 * u[0] + (t1069 * t41 - t1072) * u[1] * L + (t1069 * t53 / 0.2e1 + t1079) * u[9] * L + (t1069 * t67 - t1084) * u[2] * L + (t1069 * t79 - 0.4e1 * t1089) * u[5] * L + (-t1069 * t53 / 0.2e1 - t1079) * u[3] * L + (t1069 * t96 + 0.4e1 * t1099) * u[4] * L + (-t1069 * t41 + t1072) * u[7] * L + (t1069 * t109 + 0.2e1 * t1099) * u[10] * L + (-t1069 * t67 + t1084) * u[8] * L + t1069 * u[6] + (t1069 * t122 - 0.2e1 * t1089) * u[11] * L;
  double t1126 = -t557 - u[4] / 0.30e2 + t559 + 0.2e1 / 0.15e2 * u[10];
  double t1129 = t1126 * C12 + 0.2e1 * t563;
  double t1131 = 0.12e2 * t1129 * t3;
  double t1136 = t1126 * C13 + 0.2e1 * t578;
  double t1138 = 0.6e1 * t1136 * t32;
  double t1142 = 0.12e2 * t1136 * t3;
  double t1148 = 0.6e1 * t1129 * t32;
  double t1157 = (-t556 - t1131 - t568) * u[7] - t572 + t574 + (-t576 - t1138 - t584) * u[10] + t590 + t593 + (-t595 + t1142 + t598) * u[2] + (-t602 - t1138 - t603) * u[4] + (-t607 + t1148 + t611) * u[11] + (-t615 + t1148 + t616) * u[5] + (-t620 + t1131 + t568) * u[1] + (-t624 - t1142 - t598) * u[8];
  double t1162 = t1126 * C11 + 0.2e1 * t631;
  double t1177 = (t1126 * C14 + 0.2e1 * t646) * t29;
  double t1182 = t1129 * t29;
  double t1188 = t1136 * t29;
  double t1215 = t554 + t1157 * t1 / 0.2e1 - t1162 * u[0] + (t1162 * t67 - t1138) * u[2] * L + (t1162 * t41 - t1148) * u[1] * L + t1162 * u[6] + (-t1162 * t53 / 0.2e1 - t1177) * u[3] * L + (t1162 * t79 - 0.4e1 * t1182) * u[5] * L + (t1162 * t96 + 0.4e1 * t1188) * u[4] * L + (t1162 * t53 / 0.2e1 + t1177) * u[9] * L + (t1162 * t122 - 0.2e1 * t1182) * u[11] * L + (-t1162 * t67 + t1138) * u[8] * L + (-t1162 * t41 + t1148) * u[7] * L + (t1162 * t109 + 0.2e1 * t1188) * u[10] * L;
  double t1218 = t709 - u[5] / 0.30e2 - t711 + 0.2e1 / 0.15e2 * u[11];
  double t1221 = t1218 * C12 - 0.2e1 * t715;
  double t1223 = 0.12e2 * t1221 * t3;
  double t1228 = t1218 * C13 - 0.2e1 * t730;
  double t1230 = 0.6e1 * t1228 * t32;
  double t1234 = 0.12e2 * t1228 * t3;
  double t1240 = 0.6e1 * t1221 * t32;
  double t1249 = (t708 - t1223 + t720) * u[7] + t724 - t726 + (t728 - t1230 + t736) * u[10] + t742 + t745 + (t747 + t1234 - t750) * u[2] + (t754 - t1230 + t755) * u[4] + (t759 + t1240 - t763) * u[11] + (t767 + t1240 - t768) * u[5] + (t772 + t1223 - t720) * u[1] + (t776 - t1234 + t750) * u[8];
  double t1254 = t1218 * C11 - 0.2e1 * t783;
  double t1269 = (t1218 * C14 - 0.2e1 * t798) * t29;
  double t1274 = t1221 * t29;
  double t1280 = t1228 * t29;
  double t1307 = t706 + t1249 * t1 / 0.2e1 - t1254 * u[0] + (t1254 * t67 - t1230) * u[2] * L + (t1254 * t41 - t1240) * u[1] * L + t1254 * u[6] + (-t1254 * t53 / 0.2e1 - t1269) * u[3] * L + (t1254 * t79 - 0.4e1 * t1274) * u[5] * L + (t1254 * t96 + 0.4e1 * t1280) * u[4] * L + (t1254 * t53 / 0.2e1 + t1269) * u[9] * L + (t1254 * t122 - 0.2e1 * t1274) * u[11] * L + (-t1254 * t67 + t1230) * u[8] * L + (-t1254 * t41 + t1240) * u[7] * L + (t1254 * t109 + 0.2e1 * t1280) * u[10] * L;

  // Internal force vector.

  gl[0] = t128;
  gl[1] = t290;
  gl[2] = t448;
  gl[3] = t535;
  gl[4] = t687;
  gl[5] = t839;
  gl[6] = t862;
  gl[7] = t951;
  gl[8] = t1040;
  gl[9] = t1123;
  gl[10] = t1215;
  gl[11] = t1307;
}

// ======================================================= End of file =====
