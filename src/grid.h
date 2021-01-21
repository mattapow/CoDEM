#ifndef GRID_H
#define GRID_H
#include "main.h"
#include "box.h"
class Cgrid
{
public:

  double dx;  //lenght of one box
  double Lx;  //total lenght of the system
  int Nx;     //number of boxes
  QList <Cbox> box;  /*<< List of boxes*/


  Cgrid();
  void init(Cvector L, double step);
  void set_grain(Cgrain *);
  void clear_grain();
  void PRINT();

};

#endif // GRID_H
