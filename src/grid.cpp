#include "grid.h"


Cgrid::Cgrid(){};

void Cgrid::init(Cvector L, double step)
{
  Lx=L.x[0];
  if(step>L.x[0]) { Nx=1; dx = L.x[0];}
  else { Nx= (int) (L.x[0]/step); dx=L.x[0]/ ((double)Nx); }

  for(int i = 0;i<Nx;i++){Cbox b;  box.push_back(b);}

  for(int i = 1;i<Nx-1;i++)
    {
      box[i].left=&box[i-1];
      box[i].right=&box[i+1];
    }
  //extremes boxes and periodic boundary
  box[0].left=&box.last();
  box[0].right=&box[1];
  box[Nx-1].left=&box[Nx-2];
  box[Nx-1].right=&box[0];
}

void Cgrid::PRINT()
{
  cout<<"Grid properties"<<endl;
  cout<< "Number of boxes in x direction: "<<Nx<<" step dx "<<dx<< endl;

}



void Cgrid::clear_grain()
{
  for(int i = 0;i<Nx;i++)  box[i].grain.clear();
}


void Cgrid::set_grain(Cgrain *g)
{

  int i=(int) (g->X.x[0]/dx);
  if(i==Nx)i=Nx-1;//just in case one grain is on exacly on the right border

  if(i<0||i>Nx-1) {g->X.PRINT(); cout<<i<<" "<<Nx <<endl;}
  box[i].grain.push_back(g);
  g->box=&box[i];

}






