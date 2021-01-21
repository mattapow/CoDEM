#include "slice.h"

Cslice::Cslice()
{
  y=0;phi=0;Ome=0;
  L=dy=0;
  dV2=0;
  //shear_rate=0;
  //vorticity=0;
}




double Cslice::surfcace_area(Cgrain &g)
{

    //double YMIN =ymin;//lowest part of the grain within the slice
    //double YMAX =ymax;//highest part of the grain within the slice

    //if(YMIN<g.X.x[1] - g.X.R) YMIN = g.X.x[1] - g.X.R;
    //if(YMAX>g.X.x[1] + g.X.R) YMAX = g.X.x[1] + g.X.R;

    double Dy =y-g.X.x[1];
    if(Dy*Dy>g.R*g.R)return 0;
    double S;

    S = 2.*dy*sqrt(g.R*g.R-Dy*Dy);

    return S;
}


void Cslice::add_grain(Cgrain &g)
{
   // if(!g.FLOW)return;
  double Sgrain=surfcace_area(g);

  phi+=Sgrain;
  Ome+=g.Ome*Sgrain;
  V+=(g.V*Sgrain);
  dV2+=g.V.x[1]*g.V.x[1]*Sgrain;
  stress+=(g.stress*Sgrain);
  gradV+=(g.gradV*Sgrain);
  
}

void Cslice::average()
{
  if(phi==0)return;
  Ome/=phi;
  V/=phi;
  dV2/=phi;
  stress/=(phi);
  gradV/=phi;

  phi/=(dy*L);

  //stress/=(dy*L);

 // stress*=phi;
  //stress=phi;

}

void Cslice::operator+=(Cslice a)
{
  V+=a.V;
  Ome+=a.Ome;
  stress+=a.stress;
  phi+=a.phi;
  dV2+=a.dV2;
  gradV+=a.gradV;

}

void Cslice::operator/=(double d)
{
  V/=d;
  Ome/=d;
  phi/=d;
  stress/=d;
  dV2/=d;
  gradV/=d;
}

Cslice Cslice::operator *(double d)
{
    Cslice temp;

    temp.Ome=Ome*d;
    temp.V=V*d;
    temp.dV2=dV2*d;
    temp.phi=phi*d;
    temp.stress=stress*d;
    temp.gradV=gradV*d;
    return (temp);
}


void Cslice::PRINT(FILE *fp)
{
  fprintf(fp,"%f\t%f\t%f\t",y,y_H,phi);
  V.PRINT(fp);//4-6
  gradV.PRINT(fp);//7-15
  stress.PRINT(fp);//16-24
  V_H.PRINT(fp);//25-27
  fprintf(fp,"%f\t",shear_rate_normalized);//28
  fprintf(fp,"%f\t",Ome);//29
  fprintf(fp,"%f\t",dV2);//30
  fprintf(fp,"\n");
}

void Cslice::READ(FILE *fp)
{
  fscanf(fp,"%lf %lf",&y,&phi);
  V.READ(fp);
  gradV.READ(fp);
  stress.READ(fp);
    fscanf(fp,"%lf",&shear_rate_normalized);
    fscanf(fp,"%lf",&Ome);
    fscanf(fp,"%lf",&dV2); // was %e instead of %lf
}

