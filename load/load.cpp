// -------------------------------------------------------------------------
// load.cpp - implementation of load class.
// -------------------------------------------------------------------------
// Created:      20-Jun-2005     Evandro Parente Junior
//
// Modified:     14-May-2011     Evandro Parente Junior
//               Generalization to handle non-isoparametric elements.
//
// Modified:     08-Apr-2015     Elias Saraiva Barroso
//               Input/output operations using c++ streams.
//
// Modified:     13-Oct-2015     Evandro Parente Junior
//               Evaluation of support reactions.
//
// Modified:     06-Feb-2019     Evandro Parente Junior
//               Implementation of Donnell uniform load.
// -------------------------------------------------------------------------

#include <iostream>
#include <sstream>
#include <math.h>
#include <stdlib.h>

using namespace std;

#include "load.h"
#include "loadiga.h"
#include "timefunc.h"
#include "node.h"
#include "element.h"
#include "anmodel.h"
#include "shape.h"
#include "intpoint.h"
#include "vec.h"
#include "mat.h"
#include "gblvar.h"
#include "gbldef.h"
#include "utl.h"
#include "field.h"

// -------------------------------------------------------------------------
// Static variables:
//
cLoad* cLoad :: Head = 0;
cLoad* cLoad :: Tail = 0;

// -------------------------------------------------------------------------
// Public methods:
//

// ============================= ReadNodalLoads ============================

void cLoad :: ReadNodalLoads(void)
{
  int nload;

  // Read the number of nodal loads.

  if (!(in >> nload) || (nload < 1))
  {
    cout << "Error in the input of the number of nodal loads!\n";
    exit(0);
  }

  // Read each load.

  int id;
  double fx,fy,fz,mx,my,mz;
  for (int i = 0; i < nload; i++)
  {
    if (!(in >> id) || !(in >> fx) || !(in >> fy) || !(in >> fz) ||
        !(in >> mx) || !(in >> my) || !(in >> mz))
    {
      cout << "Error in the input of the nodal load " << i+1 << "!\n";
      exit(0);
    }

    cNode *node = cNode :: FindNode(id);
    if (!node)
    {
      cout << "Error in the input of the nodal load " << i+1 << "(invalid node)!\n";
      exit(0);
    }

    new cNodalLoad(node, fx, fy, fz, mx, my, mz); // Create the nodal load
  }
}

// ============================= ReadParConcForces ============================

void cLoad :: ReadParConcForces(void)
{
  int nload;

  // Read the number of body loads.

  if (!(in >> nload) || (nload < 1))
  {
    cout << "Error in the input of the number of parametric concentrate forces!\n";
    exit(0);
  }

  // Read each load.

  for (int i = 0; i < nload; i++)
  {
    cParConcLoad *load = new cParConcLoad( );
    load->Read( );
  }
}

// ============================= ReadPatConcForces ============================

void cLoad :: ReadPatConcForces(void)
{
  int nload;

  // Read the number of body loads.

  if (!(in >> nload) || (nload < 1))
  {
    cout << "Error in the input of the number of parametric concentrate forces!\n";
    exit(0);
  }

  // Read each load.

  for (int i = 0; i < nload; i++)
  {
    cPatConcLoad *load = new cPatConcLoad( );
    load->Read( );
  }
}

// ============================= ReadBodyForces ============================

void cLoad :: ReadBodyForces(void)
{
  int nload;

  // Read the number of body loads.

  if (!(in >> nload) || (nload < 1))
  {
    cout << "Error in the input of the number of body forces!\n";
    exit(0);
  }

  // Read each load.

  for (int i = 0; i < nload; i++)
  {
    cBodyLoad *load = new cBodyLoad( );
    load->Read( );
  }
}

// ============================ ReadLineUnifLoads ==========================

void cLoad :: ReadLineUnifLoads(void)
{
  int nload;

  // Read the number of line loads.

  if (!(in >> nload) || (nload < 1))
  {
    cout << "Error in the input of the number of uniform line loads!\n";
    exit(0);
  }

  // Read each load.

  for (int i = 0; i < nload; i++)
  {
    cEdgeUnifLoad *load = new cEdgeUnifLoad( );
    load->Read( );
  }
}

// ============================ ReadLineUnifMoments ========================

void cLoad :: ReadLineUnifMoments(void)
{
  int nload;

  // Read the number of line loads.

  if (!(in >> nload) || (nload < 1))
  {
    cout << "Error in the input of the number of uniform line loads!\n";
    exit(0);
  }

  // Read each load.

  for (int i = 0; i < nload; i++)
  {
    cEdgeUnifMoment *load = new cEdgeUnifMoment( );
    load->Read( );
  }
}

// ============================ ReadLineGenLoads ===========================

void cLoad :: ReadLineGenLoads(void)
{
  int nload;

  // Read the number of line loads.

  if (!(in >> nload) || (nload < 1))
  {
    cout << "Error in the input of the number of uniform line loads!\n";
    exit(0);
  }

  // Read each load.

  for (int i = 0; i < nload; i++)
  {
    cEdgeGenLoad *load = new cEdgeGenLoad( );
    load->Read( );
  }
}

// ============================ ReadFaceUnifLoads ==========================

void cLoad :: ReadFaceUnifLoads(void)
{
  int nload;

  // Read the number of face loads.

  if (!(in >> nload) || (nload < 1))
  {
    cout << "Error in the input of the number of uniform face loads!\n";
    exit(0);
  }

  // Read each load.

  for (int i = 0; i < nload; i++)
  {
    cFaceUnifLoad *load = new cFaceUnifLoad( );
    load->Read( );
  }
}

// ============================ ReadShellUnifLoads =========================

void cLoad :: ReadShellUnifLoads(void)
{
  int nload;

  // Read the number of face loads.

  if (!(in >> nload) || (nload < 1))
  {
    cout << "Error in the input of the number of uniform shell loads!\n";
    exit(0);
  }

  // Read each load.

  for (int i = 0; i < nload; i++)
  {
    cShellUnifLoad *load = new cShellUnifLoad( );
    load->Read( );
  }
}


// ============================ ReadLineUnifIgaLoads =======================

void cLoad :: ReadLineUnifIgaLoads(void)
{
  int nload;

  // Read the number of line loads.

  if (!(in >> nload) || (nload < 1))
  {
    cout << "Error in the input of the number of IGA uniform line loads!\n";
    exit(0);
  }

  // Read each load.

  for (int i = 0; i < nload; i++)
  {
    cEdgeUnifLoad *load = new cEdgeUnifLoadIGA( );
    load->Read( );
  }
}

// ============================ ReadLineUnifIgaMoments ====================

void cLoad :: ReadLineUnifIgaMoments(void)
{
  int nload;

  // Read the number of line loads.

  if (!(in >> nload) || (nload < 1))
  {
    cout << "Error in the input of the number of IGA uniform line loads!\n";
    exit(0);
  }

  // Read each load.

  for (int i = 0; i < nload; i++)
  {
    cEdgeUnifMoment *load = new cEdgeUnifMomentIGA( );
    load->Read( );
  }
}

// ============================ ReadFaceUnifIgaLoads =======================

void cLoad :: ReadFaceUnifIgaLoads(void)
{
  int nload;

  // Read the number of face loads.

  if (!(in >> nload) || (nload < 1))
  {
    cout << "Error in the input of the number of IGA uniform face loads!\n";
    exit(0);
  }

  // Read each load.

  for (int i = 0; i < nload; i++)
  {
    cFaceUnifLoad *load = new cFaceUnifLoadIGA( );
    load->Read( );
  }
}

// ========================== ReadPlFrameUnifLoads =========================

void cLoad :: ReadPlFrameUnifLoads(void)
{
  int nload;

  // Read the number of plframe loads.

  if (!(in >> nload) || (nload < 1))
  {
    cout << "Error in the input of the number of plane frame uniform loads!\n";
    exit(0);
  }

  // Read each load.

  for (int i = 0; i < nload; i++)
  {
    cPlFrameUnifLoad *load = new cPlFrameUnifLoad( );
    load->Read( );
  }
}

// ========================== ReadGridUnifLoads =========================

void cLoad :: ReadGridUnifLoads(void)
{
  int nload;
    
  // Read the number of grid loads.
  
  if (!(in >> nload) || (nload < 1))
  {
    cout << "Erro in the input of the number of grid uniform loads!\n";
    exit(0);
  }
    
  // Read each load.
  
  for (int i = 0; i < nload; i++)
  {
    cGridUnifLoad *load = new cGridUnifLoad( );
    load->Read( );
  }
}

// ============================ ReadLineGenIgaLoads ========================

void cLoad :: ReadLineGenIgaLoads(void)
{
  int nload;

  // Read the number of line loads.

  if (!(in >> nload) || (nload < 1))
  {
    cout << "Error in the input of the number of uniform line loads!\n";
    exit(0);
  }

  // Read each load.

  for (int i = 0; i < nload; i++)
  {
    cEdgeGenLoadIGA *load = new cEdgeGenLoadIGA( );
    load->Read( );
  }
}

// ========================== ReadFrame3DUnifLoads =========================

void cLoad :: ReadFrame3DUnifLoads(void)
{
  int nload;

  // Read the number of three-dimensional frame loads.

  if (!(in >> nload) || (nload < 1))
  {
    cout << "Error in the input of the number of 3D frame uniform loads!\n";
    exit(0);
  }

  // Read each load.

  for (int i = 0; i < nload; i++)
  {
    cFrame3DUnifLoad *load = new cFrame3DUnifLoad( );
    load->Read( );
  }
}

// ============================ ReadLineFieldLoads =========================

void cLoad :: ReadLineFieldLoads(void)
{
  int nload;

  // Read the number of line loads.

  if (!(in >> nload) || (nload < 1))
  {
    cout << "Error in the input of the number of uniform line loads!\n";
    exit(0);
  }

  // Read each load.

  for (int i = 0; i < nload; i++)
  {
    cEdgeFieldLoad *load = new cEdgeFieldLoad( );
    load->Read( );
  }
}

// ========================== ReadDonnellUnifLoads =========================

void cLoad :: ReadDonnellUnifLoads(void)
{
  int nload;
    
  // Read the number of shell loads.
  
  if (!(in >> nload) || (nload < 1))
  {
    printf("Erro in the input of the number of shell uniform loads!\n");
    exit(0);
  }
    
  // Read each load.
  
  for (int i = 0; i < nload; i++)
  {
    cDonnellUnifLoad *load = new cDonnellUnifLoad( );
    load->Read( );
  }
}

// ============================== EqvForces ================================

void cLoad :: EqvForces(cElement *elm, double t, cVector &f)
{
  int n = f.Dim( );
  cVector fe(n);

  // Compute the total equivalent force applied in the given element.

  f.Zero( );
  for (cLoad *load = Head; load != 0; load = load->Next)
  {
    if (load->Type != NODAL_LOAD && ((cElemLoad*)load)->Elem == elm)
    {
      load->ExtForce(t, fe);
      f += fe;
    }
  }
}

// ================================ Destroy ================================

void cLoad :: Destroy(void)
{
  cLoad *load = Head;

  while (load)
  {
    cLoad *l = load->Next;  // Save the next element
    delete load;            // Delete the current element
    load = l;
  }
}

// ================================ cLoad ==================================

cLoad :: cLoad(void)
{
  this->Next = 0;
  if (Head == 0)    // First element
    Head = this;
  else                  // Other elements
    Tail->Next = this;
  Tail = this;

  // Assign the current time function or use the default f(t) = 1.

  Func = cTimeFunc :: GetCurr( );
  if (!Func) Func = new cFuncConst(1.0);
}

// =============================== ~cLoad ==================================

cLoad :: ~cLoad(void)
{
}

// -------------------------------------------------------------------------
// Class cNodalLoad:
// -------------------------------------------------------------------------


// ============================== cNodalLoad ===============================

cNodalLoad :: cNodalLoad(cNode *node, double fx, double fy, double fz,
                         double mx, double my, double mz)
{
  Type = NODAL_LOAD;
  Node = node;

  LoadVal[0] = fx;
  LoadVal[1] = fy;
  LoadVal[2] = fz;
  LoadVal[3] = mx;
  LoadVal[4] = my;
  LoadVal[5] = mz;
}

// ============================= ~cNodalLoad ===============================

cNodalLoad :: ~cNodalLoad(void)
{
}

// =============================== ExtForce ================================

void cNodalLoad :: ExtForce(double t, cVector &f)
{
  double tf = Func->GetVal(t);
  for (int i = 0; i < 6; i++) f[i] = tf*LoadVal[i];
}

// ============================== AddGlobVec ===============================

void cNodalLoad :: AddGlobVec(cVector &felm, cVector &fglob)
{
  for (int i = 0; i < 6; i++)
  {
    int dof = Node->GetDof(i);
    if (dof) fglob[dof-1] += felm[i];
  }
}

// ============================= AddReactions ==============================

void cNodalLoad :: AddReactions(cVector &felm, cMatrix &R)
{
  int idx = Node->GetLabel( ) - 1;
  for (int i = 0; i < 6; i++)
  {
    R(idx, i) -= felm[i];
  }
}

// -------------------------------------------------------------------------
// Class cElemLoad:
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Public methods:
//

// =============================== cElemLoad ===============================

cElemLoad :: cElemLoad(void)
{
  IsFollower = false;
}

// ============================== ~cElemLoad ===============================

cElemLoad :: ~cElemLoad(void)
{
}

// =============================== ExtForce ================================

void cElemLoad :: ExtForce(double t, cVector &f)
{
  int i,j,k;

  // Get active directions.

  int dir[8] = {0,0,0,0,0,0,0,0};
  Elem->GetActDir(dir);
  for (i = k = 0; i < 8; i++) if (dir[i]) k++;
  double loadval[8];
  cVector actfor(k);

  // Create the integration points.

  int nipt;
  int type = Elem->GetIntType( );
  int ord[] = {Elem->GetIntOrd( ), Elem->GetIntOrd( ), Elem->GetIntOrd( )};
  cIntPoint *intpnt = cIntPoint::CreateIntPoints(Shape->GetTopologyType(), type, ord, &nipt);

  // Get the nodal coordinates.

  int nnode = Shape->GetNumNode( );
  sNodeCoord *coord = new sNodeCoord[nnode];
  if (IsFollower)
    Shape->UpdatedNodalCoord(coord);
  else
    Shape->NodalCoord(coord);

  // Create N matrix.

  int ndofn = Elem->GetNumDofNode( );
  int ndof  = ndofn*nnode;
  cMatrix N(ndofn, ndof);

  // Get memory for shape functions and derivatives.

  int     nshp   = Shape->GetNumShp( );
  int     nmap   = Shape->GetNumMap( );
  double *shpfnc = new double[nshp];
  double *mapfnc = new double[nmap];
  sNodeCoord *shpdrv = new sNodeCoord[nshp];

  // Loop over integration points

  f.Zero( );
  for (i = 0; i < nipt; i++)
  {
    // Get integration point data.

    sNatCoord p = intpnt[i].GetCoord( );
    double wgt  = intpnt[i].GetWeight( );

    // Get the external loads and extract the forces in the active directions.

    GetVal(p, coord, loadval);
    for (j = k = 0; j < 8; j++) if (dir[j]) actfor[k++] = loadval[j];

    // Compute [N] matrix.

    Shape->ShpFunc(p, shpfnc);
    Shape->MapFunc(p, mapfnc);
    Shape->NMatrix(ndofn,shpfnc,N);

    // Compute {f} += coeff*[N]t{force}.

    double detJ;
    Shape->DrvShpXYZ(p, coord, &detJ, shpdrv);
    double volc  = Elem->VolCoeff(1.0, nnode, mapfnc, coord);
    double coeff = wgt*detJ*volc;
    MultTAcc(coeff, N, actfor, f);

   //cout << "Get val " << volc << "," << loadval[0] << ", " << loadval[1] << ", " << loadval[2] << endl; 
  }

  // Multiply by the time function value.

  double tf = Func->GetVal(t);
  f *= tf;

  // Release memory.

  delete []coord;
  delete []intpnt;
  delete []shpfnc;
  delete []mapfnc;
  delete []shpdrv;
}

// ============================== AddGlobVec ===============================

void cElemLoad :: AddGlobVec(cVector &felm, cVector &fglob)
{
  // Get the associated dofs.

  int ndof = GetNumDofs( );
  int *dof = new int[ndof];
  GetDofs(dof);

  // Add local vector to the correct positions of the global one.

  for (int i = 0; i < ndof; i++)
    if (dof[i]) fglob[dof[i]-1] += felm[i];

  delete []dof;
}

// ================================ HasSupp ================================

bool cElemLoad :: HasSupp(void)
{
  int nn = Shape->GetNumNode( );
  for (int i = 0; i < nn; i++)
  {
    cNode *node = Shape->GetNode(i);
    if (node->HasSupp( )) return(true);
  }

  return(false);
}

// ============================= AddReactions ==============================

void cElemLoad :: AddReactions(cVector &felm, cMatrix &R)
{
   // Get element active directions.

  int dir[6];
  Elem->GetActDir(dir);

  // Add external forces to nodal reactions.

  int k = 0;
  int nn = Shape->GetNumNode( );
  for (int i = 0; i < nn; i++)
  {
    cNode *node = Shape->GetNode(i);
    int idx = node->GetLabel( ) - 1;
    for (int j = 0; j < 6; j++) if (dir[j])
    {
      R(idx, j) -= felm[k++];
    }
  }
}

// -------------------------------------------------------------------------
// Protected methods:
//

// =============================== GetNumDofs ==============================

int cElemLoad :: GetNumDofs(void)
{
  return Shape->GetNumNode( )*Elem->GetNumDofNode( );
}

// ================================ GetDofs ================================

void cElemLoad :: GetDofs(int *dof)
{
  // Get element active directions.

  int dir[8] = {0,0,0,0,0,0,0,0};
  Elem->GetActDir(dir);

  // Get element equations from connected nodes.

  int k = 0;
  int nn = Shape->GetNumNode( );
  for (int i = 0; i < nn; i++)
  {
    cNode *node = Shape->GetNode(i);
    for (int j = 0; j < 6; j++) if (dir[j])
    {
      dof[k++] = node->GetDof(j);
    }
  }
}

// -------------------------------------------------------------------------
// Class cParConcLoad:
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Public methods:
//

// =============================== cParConcLoad ===============================

cParConcLoad :: cParConcLoad(void)
{
  Type  = PARCONC_LOAD;
  Elem  = 0;
  Shape = 0;
  LoadVal.Resize(3);
  LoadVal.Zero( );
}

// ============================== ~cParConcLoad ===============================

cParConcLoad :: ~cParConcLoad(void)
{
}

// ================================= Read ==================================

void cParConcLoad :: Read(void)
{
  // Get loaded element.

  int elmid;
  if (in >> elmid)
    Elem = cElement :: FindElm(elmid);
  if (!Elem)
  {
    cout << "Error in the input of body forces - element " << elmid << " not found!\n";
    exit(0);
  }
  Shape = Elem->GetShape( );

  // Get parametric coord

  if (!(in >> ParVal[0]) || !(in >> ParVal[1]) || !(in >> ParVal[2]))
  {
    cout << "Error in the input of parametric concentrate force!\n";
    exit(0);
  }

  // Check if the given parametric parameter within [-1,1]

  if ((ParVal[0] < -1 || ParVal[0] > 1) || (ParVal[1] < -1 || ParVal[1] > 1) ||
      (ParVal[2] < -1 || ParVal[2] > 1))
  {
    cout << "Error in the input of parametric concentrate force!\n";
    cout << "Parametric coordenates should be defined in [-1,1]!\n";
    exit(0);
  }

  // Read load data.

  if (!(in >> LoadVal[0]) || !(in >> LoadVal[1]) || !(in >> LoadVal[2]))
  {
    cout << "Error in the input of parametric concentrate force!\n";
    exit(0);
  }
}

// =============================== ExtForce ================================

void cParConcLoad :: ExtForce(double t, cVector &f)
{
  int i,j,k;

  // Get active directions.

  int dir[6];
  Elem->GetActDir(dir);
  for (i = k = 0; i < 6; i++) if (dir[i]) k++;
  cVector actfor(k);

  // Create N matrix.

  int nnode = Shape->GetNumNode( );
  int ndofn = Elem->GetNumDofNode( );
  int ndof  = ndofn*nnode;
  cMatrix N(ndofn, ndof);

  // Get memory for shape functions.

  double *shpfnc = new double[nnode];

  f.Zero( ); 
  sNatCoord p = {ParVal[0],ParVal[1],ParVal[2]};

  // Get the external loads and extract the forces in the active directions.

  double loadval[6] = {LoadVal[0],LoadVal[1],LoadVal[2],0,0,0};
  for (j = k = 0; j < 6; j++) if (dir[j]) actfor[k++] = loadval[j];

  // Compute [N] matrix.

  Shape->ShpFunc(p, shpfnc);
  N.Zero( );
  for (j = 0; j < ndofn; j++)
    for (k = 0; k < nnode; k++) N[j][ndofn*k+j] = shpfnc[k];

  // Compute {f} = [N]t{force}.

  MultTAcc(1.0, N, actfor, f);

  // Multiply by the time function value.

  double tf = Func->GetVal(t);
  f *= tf;

  // Release memory.

  delete []shpfnc;
}

// -------------------------------------------------------------------------
// Class cBodyLoad:
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Public methods:
//

// =============================== cBodyLoad ===============================

cBodyLoad :: cBodyLoad(void)
{
  Type  = BODY_LOAD;
  Elem  = 0;
  Shape = 0;
  LoadVal.Resize(3);
  LoadVal.Zero( );
}

// ============================== ~cBodyLoad ===============================

cBodyLoad :: ~cBodyLoad(void)
{
}

// ================================= Read ==================================

void cBodyLoad :: Read(void)
{
  // Get loaded element.

  int elmid;
  if (in >> elmid)
    Elem = cElement :: FindElm(elmid);
  if (!Elem)
  {
    cout << "Error in the input of body forces - element " << elmid << " not found!\n";
    exit(0);
  }
  Shape = Elem->GetShape( );

  // Read load data.

  if (!(in >> LoadVal[0]) || !(in >> LoadVal[1]) || !(in >> LoadVal[2]))
  {
    cout << "Error in the input of body forces!\n";
    exit(0);
  }
}

// ================================ GetVal =================================

void cBodyLoad :: GetVal(sNatCoord p, sNodeCoord *coord, double *val)
{
  int i;
  for (i = 0; i < 3; i++) val[i] = LoadVal[i];
  for (i = 3; i < 6; i++) val[i] = 0.0;
}

// -------------------------------------------------------------------------
// Class cShellUnifLoad:
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Public methods:
//

// =============================== cShellUnifLoad ==========================

cShellUnifLoad :: cShellUnifLoad(void)
{
  Type  = BODY_LOAD;
  Elem  = 0;
  Shape = 0;
  Local = 1;
  LoadVal.Resize(3);
  LoadVal.Zero( );
}

// ============================== ~cShellUnifLoad ==========================

cShellUnifLoad :: ~cShellUnifLoad(void)
{
}

// ================================= Read ==================================

void cShellUnifLoad :: Read(void)
{
  // Get loaded element.

  int elmid;
  if (in >> elmid)
    Elem = cElement :: FindElm(elmid);
  if (!Elem)
  {
    cout << "Error in the input of shell uniform loads - element " << elmid << " not found!\n";
    exit(0);
  }
  Shape = Elem->GetShape( );

  // Read load data.

  if (!(in >> Local) || !(in >> LoadVal[0]) || !(in >> LoadVal[1]) || !(in >> LoadVal[2]))
  {
    cout << "Error in the input of shell uniform loads!\n";
    exit(0);
  }
}

// ================================ GetVal =================================

void cShellUnifLoad :: GetVal(sNatCoord p, sNodeCoord *coord, double *val)
{
  // Evaluate thickness.
  int     nmap    = Shape->GetNumMap( );
  double *mapfunc = new double [nmap];
  Shape->MapFunc(p,mapfunc);
  double thk = 0.0;

  for(int i = 0; i < nmap-1; ++i)
    thk += mapfunc[i] * (*coord[i].thk);

  // Get load values
  if (Local)
  {
    cMatrix R(3, 3);
    Shape->LocalSys(p, coord, R);  // Get rotation matrix
    cVector qr(3);
    qr = R*LoadVal;                // Rotate to global system
    for (int i = 0; i < 3; i++) val[i] = qr[i];
  }
  else
  {
    for (int i = 0; i < 3; i++) val[i] = LoadVal[i];
  }
  for (int i = 3; i < 6; i++) val[i] = 0.0;

  // Divide by thickness
  for (int i = 0; i < 6; i++) val[i] = val[i]/thk;

  delete [] mapfunc;
}

// -------------------------------------------------------------------------
// Class cEdgeUnifLoad:
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Public methods:
//

// ============================= cEdgeUnifLoad =============================

cEdgeUnifLoad :: cEdgeUnifLoad(void)
{
  Type  = LINE_LOAD;
  Elem  = 0;
  Shape = 0;
  Local = 1;
  LoadVal.Resize(3);
  LoadVal.Zero( );
}

// ============================ ~cEdgeUnifLoad =============================

cEdgeUnifLoad :: ~cEdgeUnifLoad(void)
{
}

// ================================= Read ==================================

void cEdgeUnifLoad :: Read(void)
{
  // Get loaded element.

  int elmid;
  if (in >> elmid)
    Elem = cElement :: FindElm(elmid);
  if (!Elem)
  {
    cout << "Error in the input of uniform line loads !\n";
    cout << "Invalid element identify: " << elmid << "!\n";
    exit(0);
  }
  cShape *elmshp = Elem->GetShape( );

  // Get loaded edge.

  int corners[2];
  int nnode;
  eShpType type;
  cNode **conn = new cNode*[elmshp->GetNumNode( )];
  int code = 0;
  if ((in >> corners[0]) && (in >> corners[1]))
    code = elmshp->GetEdge(corners, &type, &nnode, conn);
  if (!code)
  {
    cout << "Error in the input of uniform line loads! \n";
    cout << "Invalid element corners input! \n";
    exit(0);
  }
  Shape = cShape :: CreateShape(type);
  Shape->SetNodes(conn, nnode);

  // Read load data.

  if (!(in >> Local) || !(in >> LoadVal[0]) || !(in >> LoadVal[1]) || !(in >> LoadVal[2]))
  {
    cout << "Error in the input of edge uniform loads!\n";
    exit(0);
  }

  delete [] conn;
}

// ================================ GetVal =================================

void cEdgeUnifLoad :: GetVal(sNatCoord p, sNodeCoord *coord, double *val)
{
  if (Local)
  {
    cMatrix R(3, 3);
    Shape->LocalSys(p, coord, R);  // Get rotation matrix
    cVector qr(3);
    qr = R*LoadVal;                // Rotate to global system
    for (int i = 0; i < 3; i++) val[i] = qr[i];
  }
  else
  {
    for (int i = 0; i < 3; i++) val[i] = LoadVal[i];
  }
  for (int i = 3; i < 6; i++) val[i] = 0.0;
}

// -------------------------------------------------------------------------
// Class cEdgeUnifMoment:
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Public methods:
//

// ============================= cEdgeUnifMoment ===========================

cEdgeUnifMoment :: cEdgeUnifMoment(void)
{
}

// ============================ ~cEdgeUnifMoment =============================

cEdgeUnifMoment :: ~cEdgeUnifMoment(void)
{
}

// ================================ GetVal =================================

void cEdgeUnifMoment :: GetVal(sNatCoord p, sNodeCoord *coord, double *val)
{
  if (Local)
  {
    cMatrix R(3, 3);
    Shape->LocalSys(p, coord, R);  // Get rotation matrix
    cVector qr(3);
    qr = R*LoadVal;                // Rotate to global system
    for (int i = 3; i < 6; i++) val[i] = qr[i-3];
  }
  else
  {
    for (int i = 3; i < 6; i++) val[i] = LoadVal[i-3];
  }
  for (int i = 0; i < 3; i++) val[i] = 0.0;
}

// -------------------------------------------------------------------------
// Class cEdgeGenLoad:
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Public methods:
//

// ============================= cEdgeGenLoad ==============================

cEdgeGenLoad :: cEdgeGenLoad(void)
{
  Type  = LINE_LOAD;
  Elem  = 0;
  Shape = 0;
}

// ============================ ~cEdgeGenLoad ==============================

cEdgeGenLoad :: ~cEdgeGenLoad(void)
{
}

// ================================= Read ==================================

void cEdgeGenLoad :: Read(void)
{
  // Get loaded element.

  int elmid;
  if (in >> elmid)
    Elem = cElement :: FindElm(elmid);
  if (!Elem)
  {
    cout << "Error in the input of uniform line loads !\n";
    cout << "Invalid element identify: " << elmid << "!\n";
    exit(0);
  }
  cShape *elmshp = Elem->GetShape( );

  // Get loaded edge.

  int corners[2];
  int nnode;
  eShpType type;
  cNode **conn = new cNode*[elmshp->GetNumNode( )];
  int code = 0;
  if ((in >> corners[0]) && (in >> corners[1]))
    code = elmshp->GetEdge(corners, &type, &nnode, conn);
  if (!code)
  {
    cout << "Error in the input of uniform line loads! \n";
    cout << "Invalid element corners input! \n";
    exit(0);
  }
  Shape = cShape :: CreateShape(type);
  Shape->SetNodes(conn, nnode);

  // Get input expression label.
  int exprid;
  if (!(in >> exprid) || (exprid <= 0) || (exprid > Utl :: GetNumExp( )))
  {
    cout << "Invalid expression label: " << exprid << "!\n";
    exit(0);
  }

  // Setup variables table.
  expr.AddVar("x",x);
  expr.AddVar("y",y);
  expr.AddVar("z",z);
  expr.AddVar("fx",fx);
  expr.AddVar("fy",fy);
  expr.AddVar("fz",fz);

  // Compile input expression.
  expr.register_symbol_table( );
  int c = expr.compile(Utl :: GetExp(exprid));
  if (!c)
  {
    cout << "Error in compilation of the expression " << exprid << "!" << endl;
    exit(0);
  }

  delete [] conn;
}

// ================================ GetVal =================================

void cEdgeGenLoad :: GetVal(sNatCoord p, sNodeCoord *coord, double *val)
{
  // Evaluate cartesian coordinates.
  int nnode = Shape->GetNumNode( );
  double *shpfnc = new double[nnode];
  Shape->ShpFunc(p, shpfnc);

  x  = y  = z  = 0.0;
  fx = fy = fz = 0.0;
  for(int n = 0; n < nnode; ++n)
  {
    x += shpfnc[n] * coord[n].x;   
    y += shpfnc[n] * coord[n].y;   
    z += shpfnc[n] * coord[n].z;   
  }

  // Evaluate expression.
  expr.evaluate( );

  // Store evaluated forces.
  val[0] = fx;
  val[1] = fy;
  val[2] = fz;

  delete [] shpfnc;
}

// -------------------------------------------------------------------------
// Class cEdgeFieldLoad:
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Public methods:
//

// ============================= cEdgeFieldLoad ============================

cEdgeFieldLoad :: cEdgeFieldLoad(void)
{
  Type  = LINE_LOAD;
  Elem  = 0;
  Shape = 0;
}

// ============================ ~cEdgeUnifLoad =============================

cEdgeFieldLoad :: ~cEdgeFieldLoad(void)
{
}

// ================================= Read ==================================

void cEdgeFieldLoad :: Read(void)
{
  // Get loaded element.

  int elmid;
  if (in >> elmid)
    Elem = cElement :: FindElm(elmid);
  if (!Elem)
  {
    cout << "Error in the input of uniform line loads !\n";
    cout << "Invalid element identify: " << elmid << "!\n";
    exit(0);
  }
  cShape *elmshp = Elem->GetShape( );

  // Get loaded edge.

  int corners[2];
  int nnode;
  eShpType type;
  cNode **conn = new cNode*[elmshp->GetNumNode( )];
  int code = 0;
  if ((in >> corners[0]) && (in >> corners[1]))
    code = elmshp->GetEdge(corners, &type, &nnode, conn);
  if (!code)
  {
    cout << "Error in the input of uniform line loads! \n";
    cout << "Invalid element corners input! \n";
    exit(0);
  }
  Shape = cShape :: CreateShape(type);
  Shape->SetNodes(conn, nnode);

  delete [] conn;
}

// ================================ GetVal =================================

void cEdgeFieldLoad :: GetVal(sNatCoord p, sNodeCoord *coord, double *val)
{
  // Evaluate cartesian coordinates.
  int nnode = Shape->GetNumNode( );
  double *shpfnc = new double[nnode];
  Shape->ShpFunc(p, shpfnc);

  double x = 0.0, y = 0.0, z = 0.0;
  for(int n = 0; n < nnode; ++n)
  {
    x += shpfnc[n] * coord[n].x;   
    y += shpfnc[n] * coord[n].y;   
    z += shpfnc[n] * coord[n].z;   
  }

  // Get stress tensor.
  double str[9];
  cField :: GetInstance( )->Sig(x,y,z,str);

  // Evaluate rotation matrix.
  cMatrix R(3,3);
  Shape->LocalSys(p,coord,R);

  // Compute normal vector.
  cVector nx(3), n(3);
  nx.Zero( );
  nx[1] = 1.0;
  n     = R * nx;
  n.Normalize( );
  n    *= -1.0;
  
  // Compute Forces.
  val[0] = str[0] * n[0] + str[1] * n[1] + str[2] * n[2];
  val[1] = str[3] * n[0] + str[4] * n[1] + str[5] * n[2];
  val[2] = str[6] * n[0] + str[7] * n[1] + str[8] * n[2];

  cout << "val[0] = " << val[0] << endl;
  cout << "val[1] = " << val[1] << endl;
  cout << "val[2] = " << val[2] << endl;

  delete [] shpfnc;
}

// -------------------------------------------------------------------------
// Class cFaceUnifLoad:
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Public methods:
//

// ============================= cFaceUnifLoad =============================

cFaceUnifLoad :: cFaceUnifLoad(void)
{
  Type  = FACE_LOAD;
  Elem  = 0;
  Shape = 0;
  Local = 1;
  LoadVal.Resize(3);
  LoadVal.Zero( );
}

// ============================ ~cFaceUnifLoad =============================

cFaceUnifLoad :: ~cFaceUnifLoad(void)
{
}

// ================================= Read ==================================

void cFaceUnifLoad :: Read(void)
{
  // Get loaded element.

  int elmid;
  if (in >> elmid)
    Elem = cElement :: FindElm(elmid);
  if (!Elem)
  {
    cout << "Error in the input of uniform area loads !\n";
    exit(0);
  }
  cShape *elmshp = Elem->GetShape( );

  // Get loaded face.

  int corners[3];
  int nnode;
  eShpType type;
  cNode **conn = new cNode*[elmshp->GetNumNode( )];
  int code = 0;
  if ((in >> corners[0]) && (in >> corners[1]) && (in >> corners[2]))
    code = elmshp->GetFace(corners, &type, &nnode, conn);
  if (!code)
  {
    cout << "Error in the input of uniform face (area) loads !\n";
    exit(0);
  }
  Shape = cShape :: CreateShape(type);
  Shape->SetNodes(conn);

  // Read load data.

  if (!(in >> Local) || !(in >> LoadVal[0]) || !(in >> LoadVal[1]) || !(in >> LoadVal[2]))
  {
    cout << "Error in the input of face uniform loads!\n";
    exit(0);
  }
}

// ================================ GetVal =================================

void cFaceUnifLoad :: GetVal(sNatCoord p, sNodeCoord *coord, double *val)
{
  if (Local)
  {
    cMatrix R(3, 3);
    Shape->LocalSys(p, coord, R);  // Get rotation matrix
    cVector qr(3);
    qr = R*LoadVal;                // Rotate to global system
    for (int i = 0; i < 3; i++) val[i] = qr[i];
  }
  else
  {
    for (int i = 0; i < 3; i++) val[i] = LoadVal[i];
  }
  for (int i = 3; i < 6; i++) val[i] = 0.0;
  
  return;
}

// -------------------------------------------------------------------------
// Class cFaceUnifPressLoad:
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Public methods:
//

// ============================= cFaceUnifPressLoad ========================

cFaceUnifPressLoad :: cFaceUnifPressLoad(void)
{
  Type  = FACE_LOAD;
  Elem  = 0;
  Shape = 0;
  IsFollower = true;
}

// ============================ ~cFaceUnifPressLoad ========================

cFaceUnifPressLoad :: ~cFaceUnifPressLoad(void)
{
}

// ================================= Read ==================================

void cFaceUnifPressLoad :: Read(void)
{
  // Get loaded element.

  int elmid;
  if (in >> elmid)
    Elem = cElement :: FindElm(elmid);
  if (!Elem)
  {
    cout << "Error in the input of uniform area pressure loads !\n";
    exit(0);
  }
  cShape *elmshp = Elem->GetShape( );

  // Get loaded face.

  int corners[3];
  int nnode;
  eShpType type;
  cNode **conn = new cNode*[elmshp->GetNumNode( )];
  int code = 0;
  if ((in >> corners[0]) && (in >> corners[1]) && (in >> corners[2]))
    code = elmshp->GetFace(corners, &type, &nnode, conn);
  if (!code)
  {
    cout << "Error in the input of uniform face (area) pressure loads !\n";
    exit(0);
  }
  Shape = cShape :: CreateShape(type);
  Shape->SetNodes(conn);

  // Read pressure value.

  if (!(in >> PressVal))
  {
    cout << "Error in the input of face uniform pressure loads!\n";
    exit(0);
  }
}

// ================================ GetVal =================================

void cFaceUnifPressLoad :: GetVal(sNatCoord p, sNodeCoord *coord, double *val)
{
  cVector tvec(3); 
  Shape->tVector(p,coord,tvec);
  tvec.Normalize( );
  //cout << "T vec: ";
  //tvec.Print( );

  tvec *= PressVal;

  for (int i = 0; i < 3; i++) val[i] = tvec[i];
  for (int i = 3; i < 6; i++) val[i] = 0.0;
}

// -------------------------------------------------------------------------
// Class cPlFrameUnifLoad:
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Public methods:
//

// ============================ cPlFrameUnifLoad ===========================

cPlFrameUnifLoad :: cPlFrameUnifLoad(void)
{
  Type  = PLFRAME_LOAD;
  Elem  = 0;
  Shape = 0;
  Local = 1;
  LoadVal.Resize(3);
  LoadVal.Zero( );
}

// =========================== ~cPlFrameUnifLoad ===========================

cPlFrameUnifLoad :: ~cPlFrameUnifLoad(void)
{
}

// ================================= Read ==================================

void cPlFrameUnifLoad :: Read(void)
{
  // Get loaded element.

  int elmid;
  if (in >> elmid)
    Elem = cElement :: FindElm(elmid);
  if (!Elem)
  {
    cout << "Error in the input of uniform plane frame loads !\n";
    exit(0);
  }

  // Get load shape.

  Shape = Elem->GetShape( );
  if (!Shape)
  {
    cout << "Error in the input of uniform plane frame loads !\n";
    exit(0);
  }

  // Read load data.

  if (!(in >> Local) || !(in >> LoadVal[0]) || !(in >> LoadVal[1]) || !(in >> LoadVal[2]))
  {
    cout << "Error in the input of face uniform loads!\n";
    exit(0);
  }
}

// =============================== ExtForce ================================

void cPlFrameUnifLoad :: ExtForce(double t, cVector &f)
{
  if (!Shape)
  {
    cout << "Error in the input of uniform plane frame loads !\n";
    exit(0);
  }

  // Evaluate the length and orientation of the element.

  sNodeCoord coord[2];
  Shape->NodalCoord(coord);
  double dx = coord[1].x - coord[0].x;
  double dy = coord[1].y - coord[0].y;
  double L  = sqrt(dx*dx + dy*dy);
  double c  = dx/L;
  double s  = dy/L;

  // Transforma as cargas para o sistema local.

  double qx,qy;
  if (Local)
  {
    qx = LoadVal[0];
    qy = LoadVal[1];
  }
  else
  {
    qx =  LoadVal[0]*c + LoadVal[1]*s;
    qy = -LoadVal[0]*s + LoadVal[1]*c;
  }
  double mz = LoadVal[2];

  // Evaluate the end forces in the local system.

  cVector fl(2);
  fl[0] = qx*L/2.0;
  fl[1] = qy*L/2.0;

  // Evaluate the end forces in the global system.

  f.Zero( );
  f[0] = fl[0]*c - fl[1]*s;
  f[1] = fl[0]*s + fl[1]*c;
  f[2] = mz*L/2.0 + qy*L*L/12.0;
  f[3] = f[0];
  f[4] = f[1];
  f[5] = mz*L/2.0 - qy*L*L/12.0;

  // Multiply by the time function value.

  double tf = Func->GetVal(t);
  f *= tf;
}

// ================================ GetVal =================================

void cPlFrameUnifLoad :: GetVal(sNatCoord p, sNodeCoord *coord, double *val)
{
  val[0] = LoadVal[0];
  val[1] = LoadVal[1];
  val[2] = 0.0;
  val[3] = 0.0;
  val[4] = 0.0;
  val[5] = LoadVal[2];
}

// -------------------------------------------------------------------------
// Class cGridunifLoad:
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Public methods:
//

// ============================ cGridUnifLoad ==========================

cGridUnifLoad :: cGridUnifLoad(void)
{
    Type  = GRID_LOAD;
    Elem  = 0;
    Shape = 0;
    Local = 1;
    LoadVal.Resize(1);
    LoadVal.Zero( );
}

// =========================== ~cGridUnifLoad ===========================

cGridUnifLoad :: ~cGridUnifLoad(void)
{
}

// ================================= Read ==================================

void cGridUnifLoad :: Read(void)
{
  // Get loaded element.

  int elmid;
  if (in >> elmid)
    Elem = cElement :: FindElm(elmid);
  if (!Elem)
  {
    printf("Erro in the input of uniform grid loads!\n");
    exit(0);
  }
  
  // Get load shape.
  
  Shape = Elem->GetShape( );
  if (!Shape)
  {
    printf("Erro in the input of uniform grid loads!\n");
    exit(0);
  }
    
  // Read load datas.
  
  if (!(in >> LoadVal[0]))
  {
    printf("Erro in the input of uniform grid loads!\n");
    exit(0);
  }
}
// =============================== ExtForce ================================

void cGridUnifLoad :: ExtForce(double t, cVector &f)
{
  if (!Shape)
  {
    printf("Erro in the unput of uniform grid loads!\n");
    exit(0);
  }
    
  // Computer the element length.

  sNodeCoord coord[2];
  Shape->NodalCoord(coord);
  
  double dx = coord[1].x - coord[0].x;
  double dy = coord[1].y - coord[0].y;
  double dz = coord[1].z - coord[0].z;
  double l  = sqrt(dx*dx + dy*dy + dz*dz);
  double l2 =l*l;
  double qz;
   
  qz = LoadVal[0];
  double s=dy/l;
  double c=dx/l;

   // Evaluate the end forces in the local system.
  
  cVector fl(6);
  fl.Zero( );
  
  fl[ 0] =   l * qz /  2.0;  // fz1
  fl[ 1] =   0.0;            // mx1
  fl[ 2] =   -l2 * qz / 12.0;  // my1
  fl[ 3] =   l * qz /  2.0;  // fz2
  fl[ 4] =   0.0;            // mx2
  fl[ 5] =   l2 * qz / 12.0;  // my2
  
  // Evaluate the end forces in the global system.
  
  f.Zero( );
  f[ 0] =   fl[0];  // fz1
  f[ 1] =   -fl[2]*s;  // mx1
  f[ 2] =   fl[2]*c;  // my1
  f[ 3] =   fl[3];  // fz2
  f[ 4] =   -fl[5]*s;  // mx2
  f[ 5] =   fl[5]*c;  // my2
 
  // Multiply by the time function value.
  double tf = Func->GetVal(t);
  f *= tf;
}

// ================================ GetVal =================================

void cGridUnifLoad :: GetVal(sNatCoord p, sNodeCoord *coord, double *val)
{
  val[0] = 0.0;
  val[1] = 0.0;
  val[2] = LoadVal[0];
  val[3] = 0.0;
  val[4] = 0.0;
  val[5] = 0.0;
}

// -------------------------------------------------------------------------
// Class cFrame3DUnifLoad:
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Public methods:
//

// ============================ cFrame3DUnifLoad ===========================

cFrame3DUnifLoad :: cFrame3DUnifLoad(void)
{
  Type  = FRAME3D_LOAD;
  Elem  = 0;
  Shape = 0;
  Local = 1;
  LoadVal.Resize(3);
  LoadVal.Zero( );
}

// =========================== ~cFrame3DUnifLoad ===========================

cFrame3DUnifLoad :: ~cFrame3DUnifLoad(void)
{
}

// ================================= Read ==================================

void cFrame3DUnifLoad :: Read(void)
{
  // Get loaded element.

  int elmid;
  if (in >> elmid)
    Elem = cElement :: FindElm(elmid);
  if (!Elem)
  {
    cout << "Erro in the input of uniform three-dimensional frame loads!\n";
    exit(0);
  }

  // Get load shape.

  Shape = Elem->GetShape( );
  if (!Shape)
  {
    cout << "Erro in the input of uniform three-dimensional frame loads!\n";
    exit(0);
  }

  // Read load datas.

  if (!(in >> Local) || !(in >> LoadVal[0]) || !(in >> LoadVal[1]) || !(in >> LoadVal[2]))
  {
    cout << "Erro in the input of face uniform loads!\n";
    exit(0);
  }
}

// =============================== ExtForce ================================

void cFrame3DUnifLoad :: ExtForce(double t, cVector &f)
{
  if (!Shape)
  {
    cout << "Erro in the unput of uniform three-dimensional frame loads!\n";
    exit(0);
  }

  // Compute the element length.

  sNodeCoord coord[2];
  Shape->NodalCoord(coord);

  double dx = coord[1].x - coord[0].x;
  double dy = coord[1].y - coord[0].y;
  double dz = coord[1].z - coord[0].z;
  double L  = sqrt(dx*dx + dy*dy + dz*dz);
  double L2 = L*L;

  // Trasform the external loadas to the local system.

  cMatrix R(3, 3);
  Elem->CalcRotMat(R);
  double lx = R[0][0], mx = R[0][1], nx = R[0][2];
  double ly = R[1][0], my = R[1][1], ny = R[1][2];
  double lz = R[2][0], mz = R[2][1], nz = R[2][2];

  double qx,qy,qz;
  if (Local)
  {
    qx = LoadVal[0];
    qy = LoadVal[1];
    qz = LoadVal[2];
  }
  else
  {
    qx = LoadVal[0]*lx + LoadVal[1]*mx + LoadVal[2]*nx;
    qy = LoadVal[0]*ly + LoadVal[1]*my + LoadVal[2]*ny;
    qz = LoadVal[0]*lz + LoadVal[1]*mz + LoadVal[2]*nz;
  }

  // Evaluate the end forces in the local system.

  cVector fl(12);
  fl.Zero( );

  fl[ 0] =  qx*L/2.0;   // fx1
  fl[ 1] =  qy*L/2.0;   // fy1
  fl[ 2] =  qz*L/2.0;   // fz1
  fl[ 3] =  0.0;        // mx1
  fl[ 4] = -qz*L2/12.0; // my1
  fl[ 5] =  qy*L2/12.0; // mz1
  fl[ 6] =  qx*L/2.0;   // fx2
  fl[ 7] =  qy*L/2.0;   // fy2
  fl[ 8] =  qz*L/2.0;   // fz2
  fl[ 9] =  0.0;        // mx2
  fl[10] =  qz*L2/12.0; // my2
  fl[11] = -qy*L2/12.0; // mz2

  // Evaluate the end forces in the global system.

  f.Zero( );
  f[ 0] = fl[ 0]*lx + fl[ 1]*ly + fl[ 2]*lz;
  f[ 1] = fl[ 0]*mx + fl[ 1]*my + fl[ 2]*mz;
  f[ 2] = fl[ 0]*nx + fl[ 1]*ny + fl[ 2]*nz;

  f[ 3] = fl[ 3]*lx + fl[ 4]*ly + fl[ 5]*lz;
  f[ 4] = fl[ 3]*mx + fl[ 4]*my + fl[ 5]*mz;
  f[ 5] = fl[ 3]*nx + fl[ 4]*ny + fl[ 5]*nz;

  f[ 6] = fl[ 6]*lx + fl[ 7]*ly + fl[ 8]*lz;
  f[ 7] = fl[ 6]*mx + fl[ 7]*my + fl[ 8]*mz;
  f[ 8] = fl[ 6]*nx + fl[ 7]*ny + fl[ 8]*nz;

  f[ 9] = fl[ 9]*lx + fl[10]*ly + fl[11]*lz;
  f[10] = fl[ 9]*mx + fl[10]*my + fl[11]*mz;
  f[11] = fl[ 9]*nx + fl[10]*ny + fl[11]*nz;

  // Multiply by the time function value.

  double tf = Func->GetVal(t);
  f *= tf;
}

// ================================ GetVal =================================

void cFrame3DUnifLoad :: GetVal(sNatCoord p, sNodeCoord *coord, double *val)
{
  val[0] = LoadVal[0];
  val[1] = LoadVal[1];
  val[2] = LoadVal[2];
  val[3] = 0.0;
  val[4] = 0.0;
  val[5] = 0.0;
}

// -------------------------------------------------------------------------
// Class cDonnellUnifLoad:
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Public methods:
//

// ============================ cDonnellUnifLoad ===========================

cDonnellUnifLoad :: cDonnellUnifLoad(void)
{
  Type  = DONNELL_LOAD;
  Elem  = 0;
  Shape = 0;
  Cartesian = 0;
  LoadVal.Resize(3);
  LoadVal.Zero( );
}

// =========================== ~cDonnellUnifLoad ===========================

cDonnellUnifLoad :: ~cDonnellUnifLoad(void)
{
}

// ================================= Read ==================================

void cDonnellUnifLoad :: Read(void)
{
  // Get loaded element.

  int elmid;
  if (in >> elmid)
    Elem = cElement :: FindElm(elmid);
  if (!Elem)
  {
    cout << "Error in the input of Donnell loads - element " << elmid << " not found!\n";
    exit(0);
  }

  // Check if the element is valid.

  cAnModel *model = Elem->GetAnModel( );
  if (!model || (model->GetType( ) != DONNELL_SHELL))
  {
    cout << "Error in the input of loads - element " << elmid << " is not Donnell shell!\n";
    exit(0);
  }
  Shape = Elem->GetShape( );

  // Read load data.

  if (!(in >> Cartesian) || !(in >> LoadVal[0]) || !(in >> LoadVal[1]) || 
      !(in >> LoadVal[2]))
  {
    cout << "Error in the input of Donnell loads!\n";
    exit(0);
  }
}

// ================================ GetVal =================================

void cDonnellUnifLoad :: GetVal(sNatCoord p, sNodeCoord *coord, double *val)
{
  // Surface forces.

  if (Cartesian) // Loads in the cartesian system
  {
    int nnode = Shape->GetNumNode( );
    double *mapfnc = new double[nnode];
    Shape->MapFunc(p, mapfnc);
    double s = 0.0;
    double R = 0.0;
    for (int i = 0; i < nnode; i++)
    {      
      s += mapfnc[i]*coord[i].y;
      R += mapfnc[i]*coord[i].z;
    }      
    double ang = s/R;
    double qx = LoadVal[0];
    double qy = LoadVal[1];
    double qz = LoadVal[2];
    val[0] = qx;                         // qx
    val[1] = qy*cos(ang) - qz*sin(ang);  // qs
    val[2] = qy*sin(ang) + qz*cos(ang);  // qr
    delete []mapfnc;
  }
  else  // Loads in the cylindrical system
  {  
    for (int i = 0; i < 3; i++) val[i] = LoadVal[i];
  }

  // Moments.

  val[3] = val[4] = val[5] = 0.0;
}

// ======================================================= End of file =====
