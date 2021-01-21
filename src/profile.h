#ifndef PROFILE_H
#define PROFILE_H
#include "main.h"

#include "slice.h"
#include "grain.h"

class Cprofile
{
public:
  QList <Cslice> slice;
  double dy;  /**<width of the slices */
  double H,L;
  Cprofile(double ,double , double );
  void make_profile(QList <Cgrain> & grain);

  void average_y(Ccell*);
  void save(string file);

  void operator+=(Cprofile a);
  void operator/=(double d);

};

#endif // PROFILE_H
