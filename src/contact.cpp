#include "contact.h"

void Ccontact::get_distance()
{
  Cvector dX = B->X-A->X;//center to center distance

  cell->CLP(dX,Vclp);

  double dx=dX.norm();
  nAB = dX/dx;
  nABt.x[0]=nAB.x[1];
  nABt.x[1]=-nAB.x[0];
  delta =dx-A->R-B->R;
}


void Ccontact::evale_force(double dt)
{
if(delta>=0) {Force*=0; ft*=0; return;}
  Cvector dV=B->V-A->V+Vclp;


  fn_el = para->Stiff*delta;
  fn_vis= para->Visco*(dV*nAB);
 // if(fabs(fn_vis)>0.3*fabs(fn_el))fn_vis=0;
  //friction force

  nABt.x[0]=nAB.x[1];
  nABt.x[1]=-nAB.x[0];


  double VAC,VBC,uDot_t;//relative velocity at the contact point along the tangent vector direction
  VAC = A->V*nABt + A->R*A->Ome;
  VBC = B->V*nABt - B->R*B->Ome;
  uDot_t=VBC-VAC+Vclp*nABt;

  ft+=(para->Stiff*uDot_t*dt);//force (positive or negative) along the tangent direction



  double mu=para->mu;


  double fc=fabs(mu*fn_el);
  if(ft>fc)ft=fc;
  else if(ft<-fc)ft=-fc;

  //sum force in the cell frame
  // remove cohesion:
  Force = nAB*(fn_el+fn_vis+para->Cohesion)+(nABt*ft);
}

bool Ccontact::AM_I_CROSS_BOUNDARY()
{
//    if((cell->L[1]/2-fabs(A->X[1])<UNIT_LENGHT) &&(cell->L[1]/2-fabs(B->X[1])<UNIT_LENGHT) ) return true;
//    return false;
    Cvector dX = B->X-A->X;
    if (fabs(dX.x[1]) > cell->L.x[1]/2) return true;
    return false;
}


bool Ccontact::not_in_contact()
{
  if(delta>=0)return true;//not in contact
  return false;  //in contact
}

void Ccontact::PRINT(FILE *fp)
{
  fprintf(fp,"%d\t%d\t",A->ID,B->ID);
  Force.PRINT2d(fp);
  fprintf(fp,"\n");
}

void Ccontact::READ(FILE *fp)
{
//  int aID, bID;//temporary
  fscanf(fp,"%d %d",&aID,&bID);
  Force.READ2d(fp);
//  get_distance();
//  ft = Force*nABt;

}



void Ccontact::PRINT()
{
  printf("%d\t%d\t",A->ID,B->ID);
  Force.PRINT();
  printf("\n");

}








