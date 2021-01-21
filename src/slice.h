#ifndef SLICE_H
#define SLICE_H

#include "main.h"
#include "vector.h"
#include "grain.h"


class Cslice
{
public:
  Cslice();

double y,ymin,ymax,y_H;
double L,dy;
double phi;
double Ome;
double shear_rate_normalized;
Cvector V;
Cvector V_H;
double dV2;

Cmatrix stress;
Cmatrix gradV;


double surfcace_area(Cgrain &g);
void add_grain(Cgrain &g);
void average();

void operator+=(Cslice a);
void operator/=(double d);
Cslice operator *(double d);


void PRINT(FILE *fp);
void READ(FILE *fp);
};

#endif // SLICE_H
