#ifndef PARAMETER_H
#define PARAMETER_H

#include "main.h"


class Cparameter
{
public:
  Cparameter(){ P=UNIT_STRESS; shear_rate=0; E=1000*UNIT_STRESS; e=0.5;   evale_stiff(); evale_visco();};


  double E;     /*<<young modulus*/
  double Stiff; /*<<stifness of contacts N/m = E*unit lenght*/
  double e;     /*<<coefficient of restitution for grains of mass 1*/
  double Visco; /*<<damping coeficient N/(m/s)*/
  double mu=0.5; // particle friction, default 0.5

  double Cohesion;/*<< Tensile strenght of contact */

  double P;           /*<< Imposed normal stress */
  double shear_rate;  /*<< Imposed shear rate */
  void evale_visco(void);
  void evale_stiff(void);


  void PRINT();
  void PRINT(FILE *fp);
  void READ(FILE *fp);
};

#endif // PARAMETER_H
