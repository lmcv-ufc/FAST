#include "cpdata.h"
#include "node.h"

sCtrlPntData :: sCtrlPntData(double x, double y, double z, double w)
{
  coord[0] = x;
  coord[1] = y;
  coord[2] = z;
  coord[3] = w;
  label    = -1;
}

sCtrlPntData :: sCtrlPntData(const sCtrlPntData &data)
{
  coord[0] = data.coord[0];
  coord[1] = data.coord[1];
  coord[2] = data.coord[2];
  coord[3] = data.coord[3];
  label    = data.label;
}

sCtrlPntData :: sCtrlPntData(cNode *node)
{
  coord[0] = node->GetCoord( ).x;
  coord[1] = node->GetCoord( ).y;
  coord[2] = node->GetCoord( ).z;
  coord[3] = node->GetW( );

  label   = node->GetLabel( );
}

// ======================================================= End of file =====
