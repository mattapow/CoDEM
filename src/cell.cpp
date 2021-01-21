#include "cell.h"



void Ccell::CLP(Cvector &X)
{
    // Correct dX for each contact
    // Previously used X[i] less than or equal to L[i]/2
    if(X[0]>L[0]/2.){ X[0]-=L[0];}
    else if(X[0]<-L[0]/2.){ X[0]+=L[0];}

    if(X[1]>L[1]/2.){ X[1]-=L[1];X[0]-=shift; }
    else if(X[1]<-L[1]/2.){ X[1]+=L[1]; X[0]+=shift; }

    if(X[0]>L[0]/2.){ X[0]-=L[0];}
    else if(X[0]<-L[0]/2.){ X[0]+=L[0];}
}

void Ccell::CLP(Cvector &X,Cvector &Vclp)
{
  Vclp*=0;


//  if(X[0]>=L[0]/2.){ X[0]-=L[0];}
//  else if(X[0]<-L[0]/2.){ X[0]+=L[0];}


//  if(X[1]>=L[1]/2.){ X[1]-=L[1];X[0]-=shift;  Vclp= V*(-1.);}
//  else if(X[1]<-L[1]/2.){ X[1]+=L[1]; X[0]+=shift; Vclp=V;}


//  if(X[0]>=L[0]/2.){ X[0]-=L[0];}
//  else if(X[0]<-L[0]/2.){ X[0]+=L[0];}

  if(X[0]>L[0]/2.){ X[0]-=L[0];}
  else if(X[0]<-L[0]/2.){ X[0]+=L[0];}


//  // using only x velocity
//  if(X[1]>L[1]/2.){ X[1]-=L[1];X[0]-=shift;  Vclp.x[0]= V.x[0]*(-1.);}
//  else if(X[1]<-L[1]/2.){ X[1]+=L[1]; X[0]+=shift; Vclp.x[0]=V.x[0];}

  // using both x and y velocity
  if(X[1]>L[1]/2.){ X[1]-=L[1];X[0]-=shift;  Vclp-=V;}
  else if(X[1]<-L[1]/2.){ X[1]+=L[1]; X[0]+=shift; Vclp+=V;}


  if(X[0]>L[0]/2.){ X[0]-=L[0];}
  else if(X[0]<-L[0]/2.){ X[0]+=L[0];}
}


void Ccell::CLP_grain(Cvector &X,Cvector &Vgrain)
{
//    if(X[0]>=L[0]){X[0]-=L[0];}
//    else if(X[0]<0){X[0]+=L[0];}


//  if(X[1]>=L[1]/2){X[1]-=L[1];X[0]-=shift;Vgrain-=V;}
//  else if(X[1]<-L[1]/2){X[1]+=L[1];X[0]+=shift;Vgrain+=V;}

//  if(X[0]>=L[0]){X[0]-=L[0];}
//  else if(X[0]<0){X[0]+=L[0];}


    if(X[0]>L[0]){X[0]-=L[0];}
    else if(X[0]<0){X[0]+=L[0];}

  // Don't add y velocity while exchanging cell vertically
  if(X[1]>L[1]/2){X[1]-=L[1];X[0]-=shift;Vgrain.x[0]-=V.x[0];}
  else if(X[1]<-L[1]/2){X[1]+=L[1];X[0]+=shift;Vgrain.x[0]+=V.x[0];}

//      // add both x and y velocity while exchanging cell vertically
//      if(X[1]>L[1]/2){X[1]-=L[1];X[0]-=shift;Vgrain-=V;}
//      else if(X[1]<-L[1]/2){X[1]+=L[1];X[0]+=shift;Vgrain+=V;}

  if(X[0]>L[0]){X[0]-=L[0];}
  else if(X[0]<0){X[0]+=L[0];}
}

void Ccell::CLP_grain_LVL(Cvector &X, Cvector &Xtrue)
{
    // Note contents of Xtrue:
    // Xtrue.x[0] = no. of level of cell (x direction), i.e. dx
    // Xtrue.x[1] = no. of level of cell (x direction), i.e. dx
    // Xtrue.x[2] = nothing

    // Note when grain changes level along x direction
    if(X[0]>L[0]){Xtrue[0]+=1;}
    else if(X[0]<0){Xtrue[0]-=1;}

    // Note when grain changes level along y direction
    if(X[1]>L[1]/2){Xtrue[1]+=1;}
    else if(X[1]<-L[1]/2){Xtrue[1]-=1;}
}



void Ccell::predictor(const double dt)
{

    if (fixh < 0) {
        //i.e. no fixed height
        //V[1]+=A[1]*dt;
//        double eta = 1;
        V[1] = -(para->P+stress.x[1][1])/eta;
//        V[1] = -(para->P+boundary_normal_stress)/eta;

        double Vmax = 1;

        if(V[1]>Vmax)V[1]=Vmax;
        else if(V[1]<-Vmax)V[1]=-Vmax;
//        if (L.x[1]>(L.x[0]+10)) V[1] = -1;
//        cout<<"Lbefore: "<<L[1];
        L[1]+=2.*V[1]*dt;//+A[1]*dt*dt/2.;
//        cout<<" v1="<<V[1]<<" byy=:"<<boundary_normal_stress<<" Lafter: "<<L[1]<<endl;

    } else if (abs(fixh - L.x[1]) < 0.001){ //para->shear_rate*dt) {
        L.x[1] = fixh;
        V.x[1] = 0;
    } else {
        // Slowly converge towards fixed height
//        V.x[1] = 0.1 * (fixh - L.x[1]);
        // converse so that it will close
        V.x[1] = (fixh - L.x[1])/abs((fixh - L.x[1]))*1;//para->shear_rate/2;
        L.x[1]+=2.*V[1]*dt;
    }


  V[0]=L[1]*para->shear_rate;
  shift+=  V[0]*dt;
  if(shift>L[0]/2.)shift-=L.x[0];
  else if(shift<-L[0]/2.)shift+=L.x[0];
}


void Ccell::corrector(double dt)
{
//  double ap1 = -(para->P+stress.x[1][1])*L[0]*UNIT_LENGHT/(100*mass);

//  cout<<para->P<<"\t"<<stress.x[1][1]<<"\t"<<ap1<<endl;
// V[1]+=(ap1-A[1])*dt/2;
// A[1]=ap1;
}

void Ccell::evale_volume_grain(QList<Cgrain > &grain)
{
    volume_grain=0;
    foreach(Cgrain g, grain)
        volume_grain+=g.R*g.R*PI*UNIT_LENGHT;
}

void Ccell::evale_solid_fraction()
{
    phi=volume_grain/(L[0]*L[1]);
}


void Ccell::PRINT(FILE *fp,  double t)
{
  fprintf(fp,"%f\t",t);//time 1
  L.PRINT2d(fp);//2-3
  V.PRINT2d(fp);//4-5

  fprintf(fp,"%f\t",mass);//mass of the top wall 6
  fprintf(fp,"%f\t",shift);//shift of the top wall 7

  evale_solid_fraction();
  fprintf(fp,"%f\t",phi);//phi 8
  stress.PRINT2d(fp); // 9-12
  fprintf(fp,"%.8e\t%.8e\t",boundary_normal_stress,boundary_shear_stress);//boundary stresses 13,14

  fprintf(fp,"\n");
}


void Ccell::READ(FILE *fp,  double &t)
{
  fscanf(fp,"%lf",&t);//time
  L.READ2d(fp);
  V.READ2d(fp);

  fscanf(fp,"%lf",&mass); //mass of the top wall
  fscanf(fp,"%lf",&shift); //shift of the top wall
  fscanf(fp,"%lf",&phi); //phi of the top wal
  stress.READ2d(fp);
  fscanf(fp,"%lf%lf",&boundary_normal_stress,&boundary_shear_stress);

}
