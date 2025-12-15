#ifndef CTRLPNTDATA_H
#define CTRLPNTDATA_H

class cNode;

struct sCtrlPntData
{
  int label;
  double coord[4];


  sCtrlPntData(double x = 0.0, double y = 0.0, double z = 0.0, double w = 1.0);
  sCtrlPntData(const sCtrlPntData&);
  sCtrlPntData(cNode*);

  double getX(void) const      { return coord[0]; }
  double getY(void) const      { return coord[1]; }
  double getZ(void) const      { return coord[2]; }
  double getW(void) const      { return coord[3]; }
  int    getLabel(void) const  { return label;    }
  void setX(const double &val) { coord[0] = val;  }
  void setY(const double &val) { coord[1] = val;  }
  void setZ(const double &val) { coord[2] = val;  }
  void setW(const double &val) { coord[3] = val;  }
  void setLabel(const int &l)  { label = l;       }
};



#endif
