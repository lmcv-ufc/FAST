// -------------------------------------------------------------------------
// shpiga.cpp - implementation of the isogeometric shape class.
// -------------------------------------------------------------------------
// Created:      25-May-2015     Elias Saraiva Barroso
//
// Modified:     07-Feb-2018     Elias Saraiva Barroso
//               Implementation of equivalent mesh routines
//               
// -------------------------------------------------------------------------

#include <iomanip>

using namespace std;

#include "patch.h"
#include "shpiga.h"
#include "element.h"
#include "gblvar.h"
#include "bernbasis.h"
#include "utl.h"
#include "ctrlpnt.h"
#include "cpdata.h"
#include "section.h"

// -------------------------------------------------------------------------
// Static variables:
//

PatElmStdMap  cShapeIGA :: PatElmMap;
IsoPatStdMap  cShapeIGA :: IsoPatMap;
PatPntStdVec  cShapeIGA :: PatPntVec;
ElmStdVec     cShapeIGA :: ParElmVec;
EqvElmStdMap  cShapeIGA :: EqvElmMap;
int           cShapeIGA :: EqvMshSL   = 0;
int           cShapeIGA :: NumEqvNode = 0;
int           cShapeIGA :: NumEqvElm  = 0;
int           sEqvElm   :: buffsize   = 2;
int           cShapeIGA :: EqvMshDivRST[3] = {0,0,0};
double        cShapeIGA :: EqvMshTol = 1e-6;
cEqvNodeKDT   cShapeIGA :: EqvNodeTree;
static int    TargetVecSize = 150;
static int    dim;
static double val;

// -------------------------------------------------------------------------
// Auxiliary methods:
//
bool CompNodeCoordX(cNode *pnt){ return (pnt->GetCoord( ).x < val); }
bool CompNodeCoordY(cNode *pnt){ return (pnt->GetCoord( ).y < val); }
bool CompNodeCoordZ(cNode *pnt){ return (pnt->GetCoord( ).z < val); }

// ============================= SetKDTree =================================

void SetKDTree(int i, cNode **low, cNode **upp, cEqvNodeKDT &kdt)
{
  double Mean[3]     = {0.0,0.0,0.0};
  double Variance[3] = {0.0,0.0,0.0};
  int vecsize        = upp - low;

  // Evaluate mean.
  sNodeCoord coord;
  for(cNode** pnt = low; pnt != upp; ++pnt)
  {
    coord = (*pnt)->GetCoord( );
    Mean[0] += coord.x;
    Mean[1] += coord.y;
    Mean[2] += coord.z;
  }

  Mean[0] /= vecsize;
  Mean[1] /= vecsize;
  Mean[2] /= vecsize;

  // Evaluate variance.
  for(cNode** pnt = low; pnt != upp; ++pnt)
  {
    coord = (*pnt)->GetCoord( );
    Variance[0] += (Mean[0] - coord.x) * (Mean[0] - coord.x);
    Variance[1] += (Mean[1] - coord.y) * (Mean[1] - coord.y);
    Variance[2] += (Mean[2] - coord.z) * (Mean[2] - coord.z);
  }

  Variance[0] /= vecsize + 1;
  Variance[1] /= vecsize + 1;
  Variance[2] /= vecsize + 1;

  // Choose dimension based in veriance modulus.
  dim = 0;
  if (Variance[dim] < Variance[1])
    dim = 1;

  if (Variance[dim] < Variance[2])
    dim = 2;

  // Set kst current node dimension and value.
  val = Mean[dim];
  kdt.SetNode(i,dim,val);

  // shrink input vector based in the given dimension and mean value.
  cNode **mid;

  if (dim == 0)
    mid = std::partition(low,upp,CompNodeCoordX);
  else if (dim == 1)
    mid = std::partition(low,upp,CompNodeCoordY);
  else if (dim == 2)
    mid = std::partition(low,upp,CompNodeCoordZ);

  // Get sons indexes of current the node.
  int left  = kdt.GetLeft(i);
  int right = kdt.GetRight(i);

  if (right  < kdt.GetTreeSize( ))
  {
    SetKDTree(left,low,mid,kdt);
    SetKDTree(right,mid,upp,kdt);
  }
}

// ============================= FindEqvNode ===============================

sEqvNode* FindEqvNode(vector<sEqvNode*> &vec, double pos[])
{
  double tol = cShapeIGA :: GetEqvMshTol( );
  vector<sEqvNode*>::iterator it;

  for(it = vec.begin( ); it != vec.end( ); ++it)
    if (fabs((*it)->coord.x - pos[0]) < tol)
      if (fabs((*it)->coord.y - pos[1]) < tol)
        if (fabs((*it)->coord.z - pos[2]) < tol)
          return *it;

  return 0;
}


// ======================== GetEqvElmTag ===================================

string GetEqvElmTag(const eTopType &type)
{
  string tag;

  switch (type)
  {
    case LINE_TOPOLOGY:
    break;

    case TRIANGULAR_TOPOLOGY:
      tag = "%ELEMENT.T6";
    break;

    case TETRAHEDRAL_TOPOLOGY:
    break;

    case QUADRILATERAL_TOPOLOGY:
      tag = "%ELEMENT.Q8";
    break;

    case WEDGE_TOPOLOGY:
      tag = "%ELEMENT.T6";
    break;

    case HEXAHEDRAL_TOPOLOGY:
      tag = "%ELEMENT.BRICK20";
      tag = "%ELEMENT.Q8";
    break;
  }

  return tag;
}

// ======================== ReadCP =========================================

std::istream& ReadCP(std::istream &in, cPatch &pat)
{
  // Auxiliary data.
  int nl;                    // Node label.
  cNode *node;               // Node pointer.

  // Get the number of patch control points.
  int ncp = pat.getNumBas( ); 

  // Read each control point.
  sCtrlPnt** InpCP = new sCtrlPnt* [ncp];

  for(int i = 0; i < ncp; ++i)
  {
    // Read node label.
    if (!(in >> nl))
    {
      cout << "Error in read of node label!" << endl;
      in.setstate(ios::failbit);
      return in;
    }

    // Get node by its label.
    node = cNode :: FindNode(nl);
    if (!node)
    {
      cout << "Error in the read of control points!" << endl;
      cout << "Invalid node label." << endl;
      in.setstate(ios::failbit);
      return in;
    }

    // Store control point.
    InpCP[i] = new sCtrlPnt(sCtrlPntData(node));
  }
  pat.setCtrlPnt(InpCP,ncp);

  // Release memory
  delete [] InpCP;

  return in;
}

// ============================= operator>> (ePatType)  ====================

ifstream& operator>> (ifstream &in, eEdgeType &type)
{
  // Read the patch type label.

  char label[100];
  if (!Utl::ReadString(in, label))
    Utl::Exit("Error in the input of the pat type label");

  // Set the appropriate edge type.

  if (string(label) == "Top")
    type = EDGE_TOP;
  else if (string(label) == "Bottom")
    type = EDGE_BOTTOM;
  else if (string(label) == "Left")
    type = EDGE_LEFT;
  else if (string(label) == "Right")
    type = EDGE_RIGHT;
  else
    Utl::Exit("Unknown edge type: " + string(label));

  return in;
}

// ============================= operator>> (ePatType)  ====================

istream& operator>> (istream &in, eFaceType &type)
{
  // Read the patch type label.

  char label[100];
  if (!Utl::ReadString(in, label))
    Utl::Exit("Error in the input of the pat type label");

  // Set the appropriate edge type.

  if (string(label) == "Front")
    type = FACE_FRONT;
  else if (string(label) == "Back")
    type = FACE_BACK;
  else if (string(label) == "Left")
    type = FACE_LEFT;
  else if (string(label) == "Right")
    type = FACE_RIGHT;
  else if (string(label) == "Above")
    type = FACE_ABOVE;
  else if (string(label) == "Below")
    type = FACE_BELOW;
  else
    Utl::Exit("Unknown edge type: " + string(label));

  return in;
}

// -------------------------------------------------------------------------
// Wrapper methods:
// -------------------------------------------------------------------------

void ReadPatch( ) { cPatch :: ReadPatch(in);              }
void InitIGA( )   { cPatch :: setCtrlPntReadFunc(ReadCP); }

// -------------------------------------------------------------------------
// Equivalent Mesh
// -------------------------------------------------------------------------

// ============================= sEqvNode ==================================

sEqvNode :: sEqvNode(const sEqvNode &n)
{
  label = n.label;
  elm   = n.elm;
  p     = n.p;
  coord = n.coord;
}

// ============================= sEqvElm ===================================

sEqvElm :: sEqvElm(const sEqvElm &elm)
{
  // Alloc memory for node vector.
  nodes = new sEqvNode* [buffsize];

  // Copy data.
  (*this) = elm;
}


sEqvElm :: sEqvElm( )
{
  // Alloc memory for node vector.
  nodes = new sEqvNode* [buffsize];
}

// ============================= ~sEqvElm ==================================

sEqvElm :: ~sEqvElm( )
{
  // Release node vector.
  delete [] nodes;
}

// ============================= operator= =================================

void sEqvElm :: operator=(const sEqvElm &e)
{
  nn        = e.nn;
  label     = e.label;
  parentelm = e.parentelm;
  
  for(int n = 0; n < nn; ++n)
    nodes[n] = e.nodes[n];
}

// -------------------------------------------------------------------------
// Class cShapeIGA:
// -------------------------------------------------------------------------

// ============================== ReadEqvMshDiv ============================

void cShapeIGA :: ReadEqvMshDiv( )
{
  for(int i = 0; i < 3; ++i)
    if (!(in >> EqvMshDivRST[i]))
    {
      cout << "Error in the input of the equivalent mesh division! (dim = "; 
      cout << i+1 << ").\n"; 
      exit(0);
    }
}

// =========================== ReadPatPnt ==================================

void cShapeIGA :: ReadPatPnt( )
{
  // Auxiliary data.
  int     numpnt;
  int     patid;
  double  inf;
  double  sup;
  cPatch* pat;
  sPatPnt pnt;

  // Read the number of input displacement points. 
  if (!(in >> numpnt))
  {
    cout << "Error in read of number of displacement points!\n";
    exit(0);
  }

  PatPntVec.reserve(PatPntVec.size( ) + numpnt);

  // Read each point.

  for(int i = 0; i < numpnt; ++i)
  {
    if (!(in >> patid))
    {
      cout << "Error in read of patch label of patch displacement point!\n";
      exit(0);
    }
    
    pat = cPatch :: FindPatch(patid);

    if (!pat)
    {
      cout << "Invalid patch label:" << patid << "!\n";
      exit(0);
    }

    pnt.patid = patid;

    // Read parametric point coordinates.
    pnt.parcoord[0] = pnt.parcoord[1] = pnt.parcoord[2] = 0;
    for(int pvar = 0; pvar < pat->getNumParVar( ); ++pvar)
    {
      if (!(in >> pnt.parcoord[pvar]))
      {
        cout << "Error in read of patch displacement ";
	cout << "point parametric coordinate!\n";
        exit(0);
      } 

      // Check if the input parameter is located within parametric limits.
      pat->getParVarLimits(pvar,inf,sup);

      if ((pnt.parcoord[pvar] < inf) || (pnt.parcoord[pvar] > sup))
      {
        cout << "Invalid parametric location of patch displacement!\n;";
	cout << "Parametric coordinate must be located within patch bounds!\n";
	cout << "Patch parametric axis: " << pvar+1 << "\n";
	cout << "Input: " << pnt.parcoord[pvar] << "\n";
	cout << "Patch parametric limits, inf: " << inf << ", sup:" << sup << "\n";
	exit(0);
      }
    }
   
    PatPntVec.push_back(pnt);
  }
}

// ============================= Destroy =================================

void cShapeIGA :: Destroy( )
{
  // Delete Patch data.

  // Delete equivalent mesh data.
}

// ============================= InitEqvMshQuad =================================

void cShapeIGA :: InitEqvMshQuad(eTopType type)
{
  // Get quadrature object.
  sEqvElmQuad *Quad = &(GetEqvElmPnt(type));

  // Stop if quadrature is loaded.
  if (Quad->coords.size( ) > 0)
    return;
 
  // Load appropriate equivalent element mesh.
  if (type == TRIANGULAR_TOPOLOGY || type == WEDGE_TOPOLOGY )
  {
    // TODO Implement discretization here!
    
    // Compute quadrature data.
    Quad->NumElmPnt = 6;
    Quad->NumElm    = 1;
    Quad->NumPnt    = 6;

    // Alloc Memory for quadrature points.
    Quad->coords.resize(Quad->NumPnt);
    
    // Auxiliar data.
    double r,s;

    // Loop over each sub-element and evaluate its parametric
    // coords.
    sNatCoord *p = &(Quad->coords[0]);

    p->r = 0.0;     p->s = 0.0;      ++p;
    p->r = 0.5;     p->s = 0.0;      ++p;
    p->r = 1.0;     p->s = 0.0;      ++p;
    p->r = 0.5;     p->s = 0.5;      ++p;
    p->r = 0.0;     p->s = 1.0;      ++p;
    p->r = 0.0;     p->s = 0.5;      ++p;
  }
  //else if (type == QUADRILATERAL_TOPOLOGY)
  else if (type == HEXAHEDRAL_TOPOLOGY)
  {
    // Get discretization subdivision.
    int divR = EqvMshDivRST[0];
    int divS = EqvMshDivRST[1];
    
    // Compute increments in each direction.
    double dr = 2.0/(divR+2);
    double ds = 2.0/(divS+2);
    
    // Compute quadrature data.
    Quad->NumElmPnt = 8;
    Quad->NumElm    = (divR + 1) * (divS + 1);
    Quad->NumPnt    = 8 * Quad->NumElm;
    
    // Alloc Memory for quadrature points.
    Quad->coords.resize(Quad->NumPnt);
    
    // Auxiliar data.
    double r,s;

    // Loop over each sub-element and evaluate its parametric
    // coords.
    sNatCoord *p = &(Quad->coords[0]);
    for(int j = 0; j < divS+1; ++j)
      for(int i = 0; i < divR+1; ++i)
      {
        r = -1.0 + dr * i; 
        s = -1.0 + ds * j;
    
        p->r = r;       p->s = s;      ++p;
        p->r = r+dr;    p->s = s;      ++p;
        p->r = r+2*dr;  p->s = s;      ++p;
        p->r = r+2*dr;  p->s = s+ds;   ++p;
        p->r = r+2*dr;  p->s = s+2*ds; ++p;
        p->r = r+dr;    p->s = s+2*ds; ++p;
        p->r = r;       p->s = s+2*ds; ++p;
        p->r = r;       p->s = s+ds;   ++p;
      }
  }  
  //else if (type == HEXAHEDRAL_TOPOLOGY)
  else if (0)
  {
    // Get discretization subdivision.
    int divR = EqvMshDivRST[0];
    int divS = EqvMshDivRST[1];
    int divT = EqvMshDivRST[2];
    
    // Compute increments in each direction.
    double dr = 2.0/(divR+2);
    double ds = 2.0/(divS+2);
    double dt = 2.0/(divT+2);
    
    // Compute quadrature data.
    Quad->NumElmPnt = 20;
    Quad->NumElm    = (divR + 1) * (divS + 1) * (divT + 1);
    Quad->NumPnt    = 20 * Quad->NumElm;
    
    // Alloc Memory for quadrature points.
    Quad->coords.resize(Quad->NumPnt);
    
    // Auxiliar data.
    double r,s,t;

    // Loop over each sub-element and evaluate its parametric
    // coords.
    sNatCoord *p = &(Quad->coords[0]);
    for(int k = 0; k < divT+1; ++k)
      for(int j = 0; j < divS+1; ++j)
        for(int i = 0; i < divR+1; ++i)
        {
          r = -1.0 + dr * i; 
          s = -1.0 + ds * j;
          t = -1.0 + dt * k;
        
          p->r = r;       p->s = s;       p->t = t;       ++p;
          p->r = r+dr;    p->s = s;       p->t = t;       ++p;
          p->r = r+2*dr;  p->s = s;       p->t = t;       ++p;
          p->r = r+2*dr;  p->s = s+ds;    p->t = t;       ++p;
          p->r = r+2*dr;  p->s = s+2*ds;  p->t = t;       ++p;
          p->r = r+dr;    p->s = s+2*ds;  p->t = t;       ++p;
          p->r = r;       p->s = s+2*ds;  p->t = t;       ++p;
          p->r = r;       p->s = s+ds;    p->t = t;       ++p;

          p->r = r;       p->s = s;       p->t = t+dt;    ++p;
          p->r = r+2*dr;  p->s = s;       p->t = t+dt;    ++p;
          p->r = r+2*dr;  p->s = s+2*ds;  p->t = t+dt;    ++p;
          p->r = r;       p->s = s+2*ds;  p->t = t+dt;    ++p;

          p->r = r;       p->s = s;       p->t = t+2*dt;  ++p;
          p->r = r+dr;    p->s = s;       p->t = t+2*dt;  ++p;
          p->r = r+2*dr;  p->s = s;       p->t = t+2*dt;  ++p;
          p->r = r+2*dr;  p->s = s+ds;    p->t = t+2*dt;  ++p;
          p->r = r+2*dr;  p->s = s+2*ds;  p->t = t+2*dt;  ++p;
          p->r = r+dr;    p->s = s+2*ds;  p->t = t+2*dt;  ++p;
          p->r = r;       p->s = s+2*ds;  p->t = t+2*dt;  ++p;
          p->r = r;       p->s = s+ds;    p->t = t+2*dt;  ++p;
        }   
  }
}

// ============================= GetEqvElmPnt ==============================

sEqvElmQuad& cShapeIGA :: GetEqvElmPnt(eTopType type)
{
  if (type == LINE_TOPOLOGY)
  {
    static sEqvElmQuad Line;
    return Line;
  }
  else if (type == TRIANGULAR_TOPOLOGY || type == WEDGE_TOPOLOGY)
  {
    static sEqvElmQuad Trian;
    return Trian;
  }
  else if (type == QUADRILATERAL_TOPOLOGY)
  {
    static sEqvElmQuad Quad;
    return Quad;
  }
  else if (type == TETRAHEDRAL_TOPOLOGY)
  {
    static sEqvElmQuad Tet;
    return Tet;
  }
  else //if (type == HEXAHEDRAL_TOPOLOGY)
  {
    static sEqvElmQuad Hex;
    return Hex;
  }
}

// ============================= LoadEqvMsh ================================

void cShapeIGA :: LoadEqvMsh( )
{
  // Auxiliary data.
  int                nn;         // Total number of nodes in model.
  int                nen;        // Max number of nodes in element.
  double            *shpfunc;    // Element shape functions.
  sNodeCoord        *nodes;      // Element nodes.
  cNode            **modelnodes; // Model nodes.
  double             currpos[3]; // Eqv. node cartesian coords.
  sNatCoord          p;          // Eqv. node parametric coords.
  sEqvElm           *eqvelm;     // Eqv. elem data.
  sEqvNode          *eqvnode;    // Eqv. node data.
  vector<sEqvNode*> *vec;        // Pointer to eqv. node vector.
  cElement          *elm;        // Parent element.
  cShape            *shp;        // Shape object of the parrent element.
  sEqvElmQuad       *Quad;       // Quadrature of parent element.

  // Define element counters by type.
  vector<int> EqvElmByTypeCount(20,0);

  // Get a vector of model nodes (control points). 
  nn         = cNode :: GetNumNode( );
  modelnodes = new cNode* [nn];

  for(int i = 0; i < nn; ++i)
    modelnodes[i] = cNode :: GetNode(i);
  
  // Compute maximum depth considering defined target vector size. 
  int md = round(log2( double (nn)/TargetVecSize));
  md = md - 1;

  if (md < 0) md = 0;

  // Construct KDT.
  EqvNodeTree.Load(md);
  SetKDTree(0,modelnodes,modelnodes+nn,EqvNodeTree);

  // Check if KDT was constructed correctly.
  if (!EqvNodeTree.IsReady( ))
    cout << "Error in the building of KDTree!"<< endl;
   
  // Reserve the target size in each tree node data.
  for(int i = 0; i < EqvNodeTree.GetLeafSize( ); ++i)
    EqvNodeTree.GetNodeData(i).reserve(TargetVecSize);
 
  // Compute the number of elements of each type.
  int nelm = cElement :: GetNumElm( );
  eTopType toptype;
  for(int i = 0; i < nelm; ++i)
  {
    // Get current element topology type.
    toptype = cElement :: GetElm(i)->GetShape( )->GetTopologyType( );

    // Increment corresponding counters.
    EqvElmByTypeCount[(int) toptype] += 1;//Quad->NumElm;
  }

  // Alloc memory for equivalent mesh element data.
  ParElmVec.reserve(cElement :: GetNumElm( ));
  for(unsigned int t = 0; t < EqvElmByTypeCount.size( ); ++t)
    if (EqvElmByTypeCount[t] > 0)
    {
      // Load quadrature object.
      InitEqvMshQuad( (eTopType) t);

      // Get equivalent element quadrature points.
      Quad = &(GetEqvElmPnt( (eTopType) t));

      // Update eqv. element node vector buff size.
      if (sEqvElm :: buffsize < Quad->NumElmPnt)
        sEqvElm :: buffsize = Quad->NumElmPnt;
      
      // Compute total number of equivalent elements.
      EqvElmByTypeCount[t] *= Quad->NumElm;
      
      // Alloc memory for the specific type.
      EqvElmMap[(eTopType) t].reserve(EqvElmByTypeCount[t]);
    }

  // Create eqv. element.
  eqvelm = new sEqvElm( );

  // Alloc memory for shape functions.
  nen     = cElement :: GetMaxNode( );
  shpfunc = new double [nen]; 
  nodes   = new sNodeCoord[nen];
 
  // Insert all element nodes into KDT and create element lists.
  for(int i = 0; i < nelm; ++i)
  {
    // Get currnt element.
    elm = cElement :: GetElm(i);

    // Get element shape.
    shp = elm->GetShape( );

    // Get equivalent element quadrature points.
    Quad = &(GetEqvElmPnt(shp->GetTopologyType( )));

    // Insert current element in parent vector list.
    ParElmVec.push_back(elm);
    
    // Update equivalent elements counter.
    NumEqvElm  += Quad->NumElm;

    // Loop over each equivalent element of current element.
    int ncount = 0;
    for(int e = 0; e < Quad->NumElm; ++e)
    {
      eqvelm->nn = Quad->NumElmPnt;

      // Loop over each equivalent nodes defined by quadrature.
      for(int n = 0; n < Quad->NumElmPnt; ++n)
      {
        // Get current nodal parametric coords.
        p = Quad->coords[ncount++];
    
        // Evaluate element shape functions.
        elm->GetShape( )->ShpFunc(p,shpfunc);
    
        // Get element nodes.
        elm->GetShape( )->NodalCoord(nodes);
    
        // Evaluate current nodal cartesian coordinates.
        nn = elm->GetShape( )->GetNumNode( );
        currpos[0] = currpos[1] = currpos[2] = 0.0;
        for(int k = 0; k < nn; ++k)
        {
          currpos[0] += shpfunc[k] * nodes[k].x; 
          currpos[1] += shpfunc[k] * nodes[k].y; 
          currpos[2] += shpfunc[k] * nodes[k].z; 
        }

        // Get corresponding vector of eqv. nodes in kdt.
        vec = &(EqvNodeTree.GetNodeData(currpos));
    
        // Search for current node in vector.
        eqvnode = FindEqvNode(*vec,currpos); 
        if (!eqvnode)
        {
	  // Initialize eqv. node data.
	  eqvnode          = new sEqvNode( );
          eqvnode->label   = -1;
          eqvnode->elm     = elm->GetLabel( );
          eqvnode->p       = p;
          eqvnode->coord.x = currpos[0]; 
          eqvnode->coord.y = currpos[1]; 
          eqvnode->coord.z = currpos[2];
    
          // Insert eqv. node.
          vec->push_back(eqvnode);
    
          // Store created node reference in current eqv. element.
          eqvelm->nodes[n] = vec->back( );
        }
        else
          eqvelm->nodes[n] = eqvnode;
      }
      eqvelm->parentelm = elm;
      
      // Insert eqv. element.
      EqvElmMap[shp->GetTopologyType( )].push_back(*eqvelm);
    }
  }

  // Set label of each equivalent node.
  int GblNodeLabel = 0;
  NumEqvNode       = 0;
  for(int i = 0; i < EqvNodeTree.GetLeafSize( ); ++i)
  {
    vec = &(EqvNodeTree.GetNodeData(i));

    for(unsigned int n = 0; n < vec->size( ); ++n)
      (*vec)[n]->label = ++GblNodeLabel;

    NumEqvNode += vec->size( );
  }

  // Set the label of each equivalent element.
  int GblElmLabel = 0;
  std::map<eTopType, std::vector<sEqvElm> >::iterator it;
  for(it = EqvElmMap.begin( ); it != EqvElmMap.end( ); ++it)
    for(unsigned int i = 0; i < it->second.size( ); ++i)
      it->second[i].label = ++GblElmLabel;

      cout << "NumEqvNode: " << NumEqvNode << endl;
      cout << "NumEqvElm: " << NumEqvElm << endl;

  // Release memory.
  delete [] modelnodes;
  delete [] shpfunc;
  delete [] nodes;
}

// ============================= PrintEqvMshData ===========================

void cShapeIGA :: PrintEqvMshData(std::ostream &out)
{
  // Print Eqv node data.
  vector<sEqvNode*> *vec;
  out << "%NODE\n" << NumEqvNode << "\n\n";

  // Print each node.
  out << "%NODE.COORD\n" << NumEqvNode << "\n";
  out << scientific << setprecision(cControl :: GetOutPrec( ));
  for(int i = 0; i < EqvNodeTree.GetLeafSize( ); ++i)
  {
    vec = &(EqvNodeTree.GetNodeData(i));

    for(unsigned int n = 0; n < vec->size( ); ++n)
    {
      out << setw(5) << left << (*vec)[n]->label; 
      out << showpos;
      out << "     " << (*vec)[n]->coord.x;
      out << "     " << (*vec)[n]->coord.y;
      out << "     " << (*vec)[n]->coord.z;
      out << noshowpos << "\n";
    }
  }
  out << "\n";

  // Print section and integration data.
  out << "%SECTION\n0\n\n";
  out << "%INTEGRATION.ORDER\n1 1 2 2 2 2 2 2\n\n";

  // Print Eqv element data.
  out << "%ELEMENT\n" << NumEqvElm << "\n\n";

  std::map<eTopType, std::vector<sEqvElm> >::iterator it;
  for(it = EqvElmMap.begin( ); it != EqvElmMap.end( ); ++it)
  {
    out << GetEqvElmTag(it->first) << "\n";
    out << it->second.size( ) << "\n";

    for(unsigned int e = 0; e < it->second.size( ); ++e)
    { 
      // Print element label.
      out << setw(5) << left << it->second[e].label;

      // Print element section
      out << " ";
      out << setw(5) << left << it->second[e].parentelm->GetSection( )->GetLabel( ); 

      // Print element quadrature.
      if (it->second[e].parentelm->GetIntOrdIdx( ) > -1)   
      {
        out << " " << setw(5) << left;
        out << it->second[e].parentelm->GetIntOrdIdx( ) + 1; 
      }

      // Print element nodes.
      for(int n = 0; n < it->second[e].nn; ++n)
        out << " " << setw(5) << left << it->second[e].nodes[n]->label;
      out << "\n";	
    }
    out << "\n";
  }
}

// ============================= PrintEqvMshDispl ==========================

void cShapeIGA :: PrintEqvMshDispl(ostream &out)
{
  // Auxiliary data.
  int          nen;        // Max number of nodes in element.
  double      *shpfunc;    // Element shape functions.
  double       displ[6];   // Displacements.
  int          nn;         // number of nodes in element;
  cElement    *elm;        // Eqv. node element.
  cShape      *shp;        // Element shape.

  // Alloc memory for shape functions.
  nen     = cElement :: GetMaxNode( );
  shpfunc = new double [nen]; 

  // Print Eqv node data.
  vector<sEqvNode*> *vec;
  out << "%RESULT.CASE.STEP.NODAL.DISPLACEMENT\n"; 
  out << NumEqvNode << "  'displacement'\n";

  // Print each node.
  for(int i = 0; i < EqvNodeTree.GetLeafSize( ); ++i)
  {
    vec = &(EqvNodeTree.GetNodeData(i));

    for(unsigned int n = 0; n < vec->size( ); ++n)
    {
      // Get parent element of current eqv. node and its shape.
      elm = cElement :: GetElm((*vec)[n]->elm - 1);
      shp = elm->GetShape( );

      // Get number of nodes.
      nn = shp->GetNumNode( );

      // Evaluate node shape functions.
      shp->ShpFunc((*vec)[n]->p,shpfunc);

      // Evaluate displacements.
      for(int i = 0; i < 6; ++i) displ[i] = 0;
      for(int j = 0; j < nn; ++j)
        for(int k = 0; k < 6; ++k)
	  displ[k] += shpfunc[j] * shp->GetNode(j)->GetDispl(k); 

      // Print displacements.
      out << scientific << setprecision(cControl::GetOutPrec( ));
      out << setw(5) << left << (*vec)[n]->label;
      out << showpos; 
      for(int i = 0; i < 6; ++i) out << " " << displ[i];
      out << noshowpos;
      out << "\n";
    }
  }
  out << "\n";
  
  // Release memory.
  delete [] shpfunc;
}

// ============================= PrintEqvMshTemp ===========================

void cShapeIGA :: PrintEqvMshTemp(ostream &out)
{
  // Auxiliary data.
  int          nen;        // Max number of nodes in element.
  double      *shpfunc;    // Element shape functions.
  double       temp;       // Displacements.
  int          nn;         // number of nodes in element;
  cElement    *elm;        // Eqv. node element.
  cShape      *shp;        // Element shape.

  // Alloc memory for shape functions.
  nen     = cElement :: GetMaxNode( );
  shpfunc = new double [nen]; 

  // Print Eqv node data.
  vector<sEqvNode*> *vec;
  out << "%RESULT.CASE.STEP.NODAL.SCALAR\n"; 
  out << 1 << "  'temperature'\n";

  out << "\n%RESULT.CASE.STEP.NODAL.SCALAR.DATA\n"; 
  out << NumEqvNode << "\n";

  // Print each node.
  for(int i = 0; i < EqvNodeTree.GetLeafSize( ); ++i)
  {
    vec = &(EqvNodeTree.GetNodeData(i));

    for(unsigned int n = 0; n < vec->size( ); ++n)
    {
      // Get parent element of current eqv. node and its shape.
      elm = cElement :: GetElm((*vec)[n]->elm - 1);
      shp = elm->GetShape( );

      // Get number of nodes.
      nn = shp->GetNumNode( );

      // Evaluate node shape functions.
      shp->ShpFunc((*vec)[n]->p,shpfunc);

      // Evaluate displacements.
      temp = 0;
      for(int j = 0; j < nn; ++j)
        temp += shpfunc[j] * shp->GetNode(j)->GetTemp( );

      // Print temperature.
      out << scientific << setprecision(cControl::GetOutPrec( ));
      out << setw(5) << left << (*vec)[n]->label;
      out << showpos << " " << temp << noshowpos << "\n";
    }
  }
  out << "\n";
  
  // Release memory.
  delete [] shpfunc;
}


// =============================== PrintMode ===============================

void cShapeIGA :: PrintMode(int mode, double lbd, cVector &v,ostream &out,string tag1, string tag2)
{
  // Auxiliary data.
  int          nen;        // Max number of nodes in element.
  double      *shpfunc;    // Element shape functions.
  double       displ[6];   // Displacements.
  int          nn;         // number of nodes in element;
  cElement    *elm;        // Eqv. node element.
  cShape      *shp;        // Element shape.

  // Alloc memory for shape functions.
  nen     = cElement :: GetMaxNode( );
  shpfunc = new double [nen]; 

  // Print eigen value.
  out << "\n%RESULT.CASE.STEP.";
  out << tag1 << "\n";
  out << mode + 1 << "  ";
  out << showpos << right << scientific << setw(14) << setprecision(6) << lbd << "\n";

  // Print Eqv node data.
  vector<sEqvNode*> *vec;
  out << "\n%RESULT.CASE.STEP.NODAL.EIGENVECTOR";
  out << noshowpos << left << "\n" << NumEqvNode << "  '" << tag2 <<"'\n";

  // Print each node.
  for(int i = 0; i < EqvNodeTree.GetLeafSize( ); ++i)
  {
    vec = &(EqvNodeTree.GetNodeData(i));

    for(unsigned int n = 0; n < vec->size( ); ++n)
    {
      // Get parent element of current eqv. node and its shape.
      elm = cElement :: GetElm((*vec)[n]->elm - 1);
      shp = elm->GetShape( );

      // Get number of nodes.
      nn = shp->GetNumNode( );

      // Evaluate node shape functions.
      shp->ShpFunc((*vec)[n]->p,shpfunc);

      // Evaluate displacements.
      for(int i = 0; i < 6; ++i) displ[i] = 0;
      for(int j = 0; j < nn; ++j)
        for(int k = 0; k < 6; ++k)
	{
          int dof = shp->GetNode(j)->GetDof(k);
          double val = 0.0;
          if (dof > 0) val = v[dof-1];
	  displ[k] += shpfunc[j] * val; 
	}

      // Print displacements.
      out << scientific << setprecision(cControl::GetOutPrec( ));
      out << setw(5) << left << (*vec)[n]->label;
      out << showpos; 
      for(int i = 0; i < 6; ++i) out << " " << displ[i];
      out << noshowpos;
      out << "\n";
    }
  }
  out << "\n";
  
  // Release memory.
  delete [] shpfunc;
}

// ======================================================= End of file =====
