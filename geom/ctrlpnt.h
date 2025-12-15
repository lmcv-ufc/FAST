#ifndef CTRLPNT_H
#define CTRLPNT_H

#include "cpdata.h"

// -------------------------------------------------------------------------
// Foward declaration:
//
#include <iosfwd>

// -------------------------------------------------------------------------
// Struct sCtrlPnt:
//
struct sCtrlPnt
{
  protected:
   sCtrlPntData data;

  public:
   sCtrlPnt(double x = 0.0, double y = 0.0, double z = 0.0, double w = 1.0);
   sCtrlPnt(const sCtrlPntData&);
   sCtrlPnt(const sCtrlPnt&);
  ~sCtrlPnt(void);

   double getX(void) const;
   double getY(void) const;
   double getZ(void) const;
   double getW(void) const;
   int    getLabel(void) const;

   void setX(const double&);
   void setY(const double&);
   void setZ(const double&);
   void setW(const double&);
   void setLabel(const int&);

   void copy(const sCtrlPnt&);
   void zero(void);

   sCtrlPnt& operator=(const sCtrlPnt&);
   sCtrlPnt& operator+=(const sCtrlPnt&);
   sCtrlPnt& operator-=(const sCtrlPnt&);
   sCtrlPnt& operator*=(const double&);
   sCtrlPnt& operator/=(const double&);

   sCtrlPnt operator+(const sCtrlPnt &p) const;
   sCtrlPnt operator-(const sCtrlPnt &p) const;
   sCtrlPnt operator*(const double &s)   const;
   sCtrlPnt operator/(const double &s)   const;
};

std::ostream& operator<<(std::ostream&,const sCtrlPnt&);
std::istream& operator>>(std::istream&,sCtrlPnt&);
sCtrlPnt operator*(const double&,const sCtrlPnt&);
sCtrlPnt operator/(const double&,const sCtrlPnt&);

#endif // CTRLPNT_H
