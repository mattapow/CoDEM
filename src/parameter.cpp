#include "parameter.h"


void Cparameter::evale_stiff(void)
{
Stiff=E*UNIT_LENGHT;
}


void Cparameter::evale_visco(void)
{
  Visco = sqrt(UNIT_MASS/2.*Stiff)*(-2.*log(e)/sqrt(PI*PI+log(e)*log(e)) );
}


void Cparameter::PRINT()
{
  printf("Parameters\n");
  printf("\tE = %f \n\tStiffness = %f \n\tRestitution coefficient %f \n\tViscous dampling coefficient %f\n\tNormal stress = %f \n\tshear rate = %f \n\tCohesion %f\n",E,Stiff,e,Visco,P, shear_rate, Cohesion);

}


void Cparameter::PRINT(FILE *fp)
{
  fprintf(fp,"%f\t%f\t%f\t%f\t%f\t%f\t%f",E,Stiff,e,Visco,P, shear_rate,Cohesion);
  fprintf(fp,"\n");
}


void Cparameter::READ(FILE *fp)
{
  fscanf(fp,"%lf %lf %lf %lf %lf %lf %lf",&E,&Stiff,&e,&Visco,&P, &shear_rate, &Cohesion);

  e=0.5;
  evale_stiff(); evale_visco();

}

