#include "node.h"

Cnode::Cnode()
{
    phi=0;
    weight=0;
    shear_rate=0;
}


Cvector X,V;        /**<POsition and CG velocity */
double phi;         /**<Solid fraction */
double shear_rate;  /**<Shear rate */

double weight;      /**<Total weight of grains */

double Cnode::get_weigth(Cgrain &g)
{
 Cvector dX=g.X-X;
 Cvector temp;
 g.cell->CLP(dX,temp);

 return ( exp(-(dX*dX)/(2*g.R*g.R))/sqrt(2*PI*g.R*g.R) );
}

void Cnode::add_grain(Cgrain &g)
{
    double w = get_weigth(g);

    weight+=w;
    V+=g.V*w;
    stress+=g.stress*w;
}

void Cnode::average()
{
    if(weight==0)return;
    V/=weight;
    stress/=weight;

    shear_rate/=weight;
}


void Cnode::operator +=(Cnode n)
{
    weight+=n.weight;
    V+=n.V*n.weight;
    shear_rate+=n.shear_rate*n.weight;
    stress+=n.stress*n.weight;
}

void Cnode::operator /=(double d)
{
    if(d==0)return;
    V/=d;
    shear_rate/=d;
    stress/=d;
}

void Cnode::PRINT(FILE *fp)
{
 // fprintf(fp,"%f\t%f\t",y,phi);
  X.PRINT(fp);
  V.PRINT(fp);
  fprintf(fp,"%f\t",shear_rate);//c 7
  stress.PRINT(fp);//c 8-16
  fprintf(fp,"\n");
}

void Cnode::READ(FILE *fp)
{
//  fscanf(fp,"%lf %lf",&y,&phi);
  X.READ(fp);
  V.READ(fp);
  fscanf(fp,"%lf",&shear_rate);
  stress.READ(fp);
}
