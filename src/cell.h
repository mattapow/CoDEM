#ifndef CELL_H
#define CELL_H

#include "vector.h"
#include "parameter.h"
#include "grain.h"

class Cgrain;


class Ccell
{
public:
  Ccell(){ mass=0;shift=0;};

  Cvector L,V,A;   //size of the cell velocity and acceleration of top wall (inverse for bottom wall)
  Cvector Fsum;    //sum of forces that act of the two walls - positive force tend to expend the cell
  double mass;

  double  volume_grain,phi;


  double shift;
  double fixh; // If negative value, height will be varied, if positive value, height of system will be fixed at this value
  double eta; // viscous damping used to control system height in order to maintain pressure

  Cmatrix stress;       /*<< Average stress in the cell*/
//  Cmatrix stress_boundary;       /*<< Average stress near boundary the cell*/
  double boundary_normal_stress;
  double boundary_shear_stress;

  Cparameter *para;

  void CLP(Cvector &X);                 //
  void CLP(Cvector &,Cvector &) ;       //reset the vector first component within +-L.x[0]/2
  void CLP_grain(Cvector &,Cvector &);  //set the grain x position between 0 and L.x[0]
  void CLP_grain_LVL(Cvector &X, Cvector &Xtrue);

  void evale_volume_grain(QList <Cgrain > &);
  void evale_solid_fraction();
  void predictor(const double dt);  //change the velocity of walls
  void corrector(double dt);

  void PRINT(FILE *,double t);   /*<< prints cell data in file with the time*/
  void READ(FILE *,double &t);  /*<< reads cell data from file*/

};

#endif // CELL_H
