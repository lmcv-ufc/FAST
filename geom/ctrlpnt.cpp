#include "ctrlpnt.h"
#include <iostream>

// -------------------------------------------------------------------------
// Struct sCtrlPnt:
//

// -------------------------------------------------------------------------
// Public methods:
//

// ============================== sCtrlPnt =================================

sCtrlPnt :: sCtrlPnt(double x, double y, double z, double w) : data(x,y,z,w)
{
}

sCtrlPnt :: sCtrlPnt(const sCtrlPntData &cpdata) : data(cpdata)
{
}

sCtrlPnt :: sCtrlPnt(const sCtrlPnt &p) : data(p.data)
{
}

// ============================== ~sCtrlPnt ==================================

sCtrlPnt :: ~sCtrlPnt(void)
{
}

// ============================== copy =====================================

void sCtrlPnt :: copy(const sCtrlPnt &p)
{
  setX(p.getX( ));
  setY(p.getY( ));
  setZ(p.getZ( ));
  setW(p.getW( ));
}

// ============================== getX =====================================

double sCtrlPnt :: getX( ) const
{
  return data.getX( );
}

// ============================== getY =====================================

double sCtrlPnt :: getY( ) const
{
  return data.getY( );
}

// ============================== getZ =====================================

double sCtrlPnt :: getZ( ) const
{
  return data.getZ( );
}

// ============================== getW =====================================

double sCtrlPnt :: getW( ) const
{
  return data.getW( );
}

// ============================== getLabel =================================

int sCtrlPnt :: getLabel( ) const
{
  return data.getLabel( );
}

// ============================== setX =====================================

void sCtrlPnt :: setX(const double &val)
{
  data.setX(val);
}

// ============================== setY =====================================

void sCtrlPnt :: setY(const double &val)
{
  data.setY(val);
}

// ============================== setZ =====================================

void sCtrlPnt :: setZ(const double &val)
{
  data.setZ(val);
}

// ============================== setW =====================================

void sCtrlPnt :: setW(const double &val)
{
  data.setW(val);
}

// ============================== setLabel =================================

void sCtrlPnt :: setLabel(const int &l)
{
  data.setLabel(l);
}

// ============================== zero =====================================

void sCtrlPnt :: zero(void)
{
  data.setX(0.0);
  data.setY(0.0);
  data.setZ(0.0);
  data.setW(0.0);
}

// ============================== operator= ================================

sCtrlPnt& sCtrlPnt :: operator=(const sCtrlPnt &p)
{
  setX(p.getX( ));
  setY(p.getY( ));
  setZ(p.getZ( ));

  return *this;
}

// ============================== operator+= ===============================

sCtrlPnt& sCtrlPnt :: operator+=(const sCtrlPnt &p)
{
  setX(getX( ) + p.getX( ));
  setY(getY( ) + p.getY( ));
  setZ(getZ( ) + p.getZ( ));

  return *this;
}

// ============================== operator-= ===============================

sCtrlPnt& sCtrlPnt :: operator-=(const sCtrlPnt &p)
{
  setX(getX( ) - p.getX( ));
  setY(getY( ) - p.getY( ));
  setZ(getZ( ) - p.getZ( ));

  return *this;
}

// ============================== operator/= ===============================

sCtrlPnt& sCtrlPnt :: operator/=(const double &scalar)
{
  setX(getX( )/scalar);
  setY(getY( )/scalar);
  setZ(getZ( )/scalar);

  return *this;
}

// ============================== operator*= ===============================

sCtrlPnt& sCtrlPnt :: operator*=(const double &scalar)
{
  setX(getX( )*scalar);
  setY(getY( )*scalar);
  setZ(getZ( )*scalar);

  return *this;
}

// ============================== operator+ ===============================

sCtrlPnt sCtrlPnt :: operator+(const sCtrlPnt &p) const
{
  return sCtrlPnt(*this) += p;
}

// ============================== operator- ===============================

sCtrlPnt sCtrlPnt :: operator-(const sCtrlPnt &p) const
{
  return sCtrlPnt(*this) -= p;
}

// ============================== operator* ===============================

sCtrlPnt sCtrlPnt :: operator*(const double &s) const
{
  return sCtrlPnt(*this) *= s;
}

// ============================== operator/ ===============================

sCtrlPnt sCtrlPnt :: operator/(const double &s) const
{
  return sCtrlPnt(*this) /= s;
}

// -------------------------------------------------------------------------
// Auxiliary methods:
//

// ============================== operator* ================================

sCtrlPnt operator*(const double &scalar, const sCtrlPnt &v)
{
  return sCtrlPnt(v) *= scalar;
}

// ============================== operator<< ===============================

std::ostream& operator<<(std::ostream &out, const sCtrlPnt &p)
{
  out << p.getX( ) << ' ' << p.getY( ) << " " << p.getZ( ) << " " << p.getW( );
  return out;
}

// ============================== operator>> ===============================

std::istream& operator>>(std::istream &in, sCtrlPnt &p)
{
  double x,y,z,w;
  if (!(in >> x) || !(in >> y) || !(in >> z) || !(in >> w))
  {
    std::cout << "Error in the input of control point coordinates!" << std::endl;
    in.setstate(std::ios::failbit);
    return in;
  }

  p.setX(x);
  p.setY(y);
  p.setZ(z);
  p.setW(w);

  return in;
}

// ======================================================= End of file =====
