#include "profile.h"

Cprofile::Cprofile(double h, double l,double d )
{
  dy=d;
  H=h;
  L=l;
   int Ny=H/dy;

  for(int i=0;i<Ny+1;i++)
  {
     Cslice s;
     s.y=i*dy-H/2.;
     s.dy=dy;
     s.L=L;
     slice.push_back(s);
   }
}


void Cprofile::make_profile(QList <Cgrain> & grain)
{

//    for(int i=0;i<grain.size();i++)    grain[i].stress*=0.;
//    for(int i=0;i<grain.size();i++)    grain[i].evale_stress();

    for(int i=0;i<grain.size();i++)
    {
        for(int j=0;j<slice.size();j++)
            slice[j].add_grain(grain[i]);
    }

    for(int i=0;i<slice.size();i++)
        slice[i].average();

}

void Cprofile::average_y(Ccell *cell)
{
    Cprofile raw(H,L,dy);
    for(int i=0;i<raw.slice.size();i++)
        raw.slice[i]=slice[i];

    int Di=20;
    for(int i=0;i<raw.slice.size();i++)
    {
        double wT=0;
        for(int k=-Di;k<=Di;k++)
        {
            double j = i+k;
            if(j>=0&&j<raw.slice.size())
                if(fabs(raw.slice[j].y)<H/2.)
                {
                    double weight = 1;
                   // cout<<weight<<endl;
                    wT+=weight;
                    slice[i].phi+=(raw.slice[j].phi*weight);
                    slice[i].stress+=(raw.slice[j].stress*weight);
                }
         }
     //   cout<< wT<<endl<<endl;
        if(wT>0) slice[i].phi/=wT;
        if(wT>0) slice[i].stress/=wT;
    }

     for(int i=0;i<slice.size();i++)
     {
         slice[i].y_H=slice[i].y/(H/2.);
         slice[i].V_H=slice[i].V;
         slice[i].V_H[0]/=cell->V[0];
//         slice[i].shear_rate_normalized=slice[i].shear_rate/(cell->V[0]/(cell->L[1]/2.));
     }

}

void Cprofile::operator+=(Cprofile a)
{
 for(int i=0;i<slice.size();i++)
      slice[i]+=a.slice[i];
}

void Cprofile::operator/=(double d)
{
  for(int i=0;i<slice.size();i++)
  slice[i]/=d;
}


void Cprofile::save(string file)
{
  FILE *fp;

  fp=fopen(file.c_str(),"w");
  foreach(Cslice s,slice)if(fabs(s.y)<H/2. /*-UNIT_LENGHT*/)s.PRINT(fp);
  fclose(fp);
}

