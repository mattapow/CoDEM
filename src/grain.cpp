#include "grain.h"
#include <math.h>
#include "voro/voro++.hh"
using namespace voro;


void Cgrain::evale_mass(void)
{
  mass = R*R*UNIT_LENGHT*4.; //mass corresponding to a cylinder of radius R, length unity, and density 4/PI (->mass=1 for R=0.5)
  Imass =  R*R*mass;
}



/*void Cgrain::get_neighbour() //add grains from the box, the right box and the left box ; only add grains of higher Id not to count the neighbour pair twice
{
//  if(box==null)return;
  neighbour.clear();
  foreach(Cgrain *B,box->grain)         if(B->ID > ID) if(fabs(X[1]-B->X[1])<R+B->R+0.2*UNIT_LENGHT) neighbour.push_back(B);
  foreach(Cgrain *B,box->right->grain)  if(B->ID > ID) if(fabs(X[1]-B->X[1])<R+B->R+0.2*UNIT_LENGHT) neighbour.push_back(B);
  foreach(Cgrain *B,box->left->grain)   if(B->ID > ID) if(fabs(X[1]-B->X[1])<R+B->R+0.2*UNIT_LENGHT) neighbour.push_back(B);
}*/

void Cgrain::predictor(const double dt)
{
  V+=A*dt;
  double old_y = X.x[1];
  X+=V*dt+A*(dt*dt/2);
//  Xtrue+=V*dt+A*(dt*dt/2);



  Ome+=Ome_dot*dt;
  theta+=Ome*dt+Ome_dot*(dt*dt/2);


  // Note contents of Xtrue:
  // Xtrue.x[0] = no. of level of cell (x direction), i.e. dx [+1 if went right]
  // Xtrue.x[1] = no. of level of cell (x direction), i.e. dx [+1 if went up]
  // Xtrue.x[2] = nothing
  //
  // update level change
   cell->CLP_grain_LVL(X,Xtrue);
  cell->CLP_grain(X,V);

}
void Cgrain::corrector(const double dt)
{
  Cvector Ap = Fsum/mass;
  double Ome_dotp=Msum/Imass;

  Ome+=(Ome_dotp-Ome_dot)*dt/2.;
  Ome_dot=Ome_dotp;
  V+=(Ap-A)*dt/2.;
  A=Ap;

}

bool is_the_same(Ccontact C,Ccontact D)
{
    if(C.A->ID==D.A->ID)
        if(C.B->ID==D.B->ID) return true;
    if(C.A->ID==D.B->ID)
        if(C.B->ID==D.A->ID) return true;
    return false;
}

void Cgrain::update_contact(const double dt)
{

  QList <Ccontact> contact_old;//copy the list of contact
  for(int c=0;c<contact.size();c++)contact_old.push_back(contact[c]);

  contact.clear();//clear the old list of contact

 for(int i=0;i<neighbour.size();i++)//for all the grain's neighbours
  {
    Ccontact C(this,neighbour[i],cell,para);//create a new contact
    C.get_distance();
    if(!C.not_in_contact()) contact.push_back(C);
   }

   for(int c=0;c<contact.size();c++)
    for(int c1=0;c1<contact_old.size();c1++)
        if(is_the_same(contact[c],contact_old[c1]))        {contact[c].ft=contact_old[c1].ft; break;}

    //once the contacts have been detected, find their force
  for(int c=0;c<contact.size();c++)  contact[c].evale_force(dt);

}

void Cgrain::evale_shear_rate(QList <Cgrain> grain)
{
    shear_rate=0;

    //  double w_tot=0,w;
 //   cout<<R<<endl;
    double Dist=1.5*UNIT_LENGHT; //center to center distance max for getting near grains
    QList <Cvector> DX;//list of grains not far
    QList <Cvector> DV;
    for(int i=0;i<grain.size();i++)
        //if(grain[i].FLOW)
        if(ID!=grain[i].ID)
        {
            Cvector DX_temp;
            Cvector DV_temp;
            DX_temp=grain[i].X-X;
            cell->CLP(DX_temp,DV_temp);
            if(DX_temp.norm()<Dist)
                {
                DX.push_back(DX_temp);
                DV_temp+=grain[i].V-V;
                DV.push_back(DV_temp);
                }
        }

    if(DX.size()==0) return;
    Cmatrix LL,VL;

    for(int i=0;i<DX.size();i++)
    {
        LL+= DX[i]^DX[i];
        VL+= DX[i]^DV[i];
    }

        LL/=DX.size();
        VL/=DX.size();


    Cmatrix LL_inv;
    LL_inv=LL.inverse();

    gradV = LL_inv*VL;
    shear_rate= gradV.sym();
    vorticity=gradV.asym();
}

void Cgrain::evale_stress()
{
     foreach(Ccontact c, contact)
        {
         c.get_distance();
         stress+=c.Force^(c.nAB*R/(PI*R*R));
         c.B->stress+=c.Force^(c.nAB*c.B->R/(PI*c.B->R*c.B->R))*(R)/(c.A->R+c.B->R);
        }
}

bool is_this_g_ok(Cgrain * g_check, QList <Cgrain *> grns)
{
    bool this_g_ok = true;
    foreach(Cgrain *nG, grns)
        if(g_check->ID==nG->ID ) {this_g_ok = this_g_ok && false; break;}
    return this_g_ok;
}

void Cgrain::evale_from_voro(int lvl)
{
    // For stress and gradV
    // list of grains to consider
    QList <Cgrain *> grain_for_stress;
    QList <Cgrain *> boundary_grain_1;
    QList <Cgrain *> boundary_grain_2;

    // level 0
    // insert
    grain_for_stress.push_back(this);
    double net_vol = voro_area; // start with this vol
    double net_s_vol = PI*R*R;

    // level 1
    foreach(Cgrain *nG1, voro_neighbour) {
        // insert
        grain_for_stress.push_back(nG1);
        boundary_grain_1.push_back(nG1);
        net_vol+=nG1->voro_area;
        net_s_vol+=PI*nG1->R*nG1->R;
    }

    if (lvl >= 2) {
        // level 2
        foreach(Cgrain *bG1,boundary_grain_1)
        foreach(Cgrain *nG2, bG1->voro_neighbour) {
            if (is_this_g_ok(nG2,grain_for_stress)) {
                // insert
                grain_for_stress.push_back(nG2);
                boundary_grain_2.push_back(nG2);
                net_vol+=nG2->voro_area;
                net_s_vol+=PI*nG2->R*nG2->R;
            }
        }
    }

    if (lvl >= 3) {
        // level 3
        foreach(Cgrain *nG3, boundary_grain_2) {
            if (is_this_g_ok(nG3,grain_for_stress)) {
                // insert
                grain_for_stress.push_back(nG3);
                net_vol+=nG3->voro_area;
            }
        }
    }

    // calulates gradV for each grain_for_stress
    QList <Cvector> DX;//list of grains not far
    QList <Cvector> DV;

    // loop through each neighbouring grains
    foreach(Cgrain *nG, grain_for_stress)
        if(ID!= nG->ID)
        {
            // for gradV
            Cvector DX_temp;
            Cvector DV_temp;
            DX_temp=nG->X-X;
            cell->CLP(DX_temp,DV_temp);
            DX.push_back(DX_temp);
            DV_temp+=nG->V-V;
            DV.push_back(DV_temp);
        }

    if(DX.size()==0) return;
    Cmatrix LL,VL;

    for(int i=0;i<DX.size();i++)
    {
        LL+= DX[i]^DX[i];
        VL+= DX[i]^DV[i];
    }
    LL/=DX.size();
    VL/=DX.size();

    Cmatrix LL_inv;
    LL_inv=LL.inverse();

    gradV = LL_inv*VL;
    shear_rate= gradV.sym();
    vorticity=gradV.asym();

//    evale_stress();

    // Compute stress for each grain_for_stress
    foreach(Cgrain *sG, grain_for_stress) {
//        if (ID==47) cout<<"contact size of "<<sG->ID<< "is: "<<sG->contact.size()<<endl;
        foreach(Ccontact c, sG->contact)
        {
            c.get_distance();
            stress-=c.Force^(c.nAB*(c.A->R+c.B->R+c.delta))*(sG->R)/(c.A->R+c.B->R);
            // add to list
//            if (ID==47) cout<<"adding stress from A"<<endl;
         }
    }

    // Compute stress from boundary grains
//    if (ID==88) cout<<"boundary size "<<boundary_grain_2.size()<<endl;
//    if (ID==88) cout<<"Msum: "<<Msum<<endl;

//    Cmatrix net_body_torque;
//    foreach(Cgrain *nG, grain_for_stress){
//        net_body_torque.x[0][1]-=nG->Msum/2;
//        net_body_torque.x[1][0]+=nG->Msum/2;
//    }


//    foreach(Cgrain *bG, boundary_grain_2) {
////        if (ID==88) cout<<"contact size "<<bG->contact.size()<<endl;
//        foreach (Ccontact c, bG->contact) {
//            if (c.A->ID == bG->ID && is_this_g_ok(c.B,grain_for_stress)) {
//                // this means, the force is on bG
//                c.get_distance();
//                Cvector ab = c.B->X - bG->X;
//                cell->CLP(ab);
//                // postion of contact point, con_x
//                Cvector con_x = (bG->X + ab*(c.B->R/(bG->R+c.B->R)));
//                Cvector del_x = con_x-X;
//                cell->CLP(del_x);
//                // sum to stress
////                if (ID==88) cout<<"adding stress from A"<<endl;
//                stress-=(c.Force^del_x);
//            } else if (c.B->ID == bG->ID && is_this_g_ok(c.A,grain_for_stress)){
//                // this means, the force is from bG
//                c.get_distance();
//                Cvector ab = c.A->X - bG->X;
//                cell->CLP(ab);
//                // postion of contact point, con_x
//                Cvector con_x = (bG->X + ab*(c.A->R/(c.A->R+bG->R)));
//                Cvector del_x = con_x-X;
//                cell->CLP(del_x);
//                // sum to stress
//                stress+=(c.Force^del_x);
////                if (ID==88) cout<<"adding stress from B"<<endl;
//            }
//        }
//    }

    // add net body torque
//    stress-=net_body_torque;
    // divide by net volume
    stress/=net_vol;
    voro_area_n = net_vol;
    solid_area_n = net_s_vol;

}



Cvector Cgrain::interp_vel(Cvector pos, QList <Cgrain> interp_grains)
{
    //QList <double> dist; // stores dist to each grain center from pos
    QList <double> dist_mr; // dist minus radius
    QList <Cvector> dir_v; // direction vector towards grain

    // compute dist to each grain center from pos
    for (int i=0; i<interp_grains.size();i++)
    {
        Cvector dx = interp_grains[i].X-pos; // normal vector
        Cvector dx_tan; // tangent vector
        dx_tan.x[0] = dx.x[1];
        dx_tan.x[1] = -dx.x[0];
        //dist.push_back(dx.norm());
        dist_mr.push_back(dx.norm()-interp_grains[i].R);
        dir_v.push_back(dx_tan/dx.norm());
    }

    Cvector final_vel;
    for (int i=0; i<interp_grains.size();i++)
    {
        double this_factor = 1;
        for (int j=0; j<interp_grains.size();j++)
            if (i!=j) this_factor *= dist_mr[j]/(dist_mr[j]+dist_mr[i]);
        final_vel += (interp_grains[i].V + dir_v[i]*(interp_grains[i].Ome*interp_grains[i].R))*this_factor;
    }
    return final_vel;
}


void Cgrain::PRINT(FILE *fp)
{
  fprintf(fp,"%d\t",ID); //1
  X.PRINT2d(fp); //2-3
  V.PRINT2d(fp); //4-5
  Xtrue.PRINT2d(fp); //6-7


  fprintf(fp,"%.6e\t%.6e\t%.6e\t%.6e\t",Ome,mass,Imass,R); //8-11
  //A.PRINT(fp);//15-17

  fprintf(fp,"\n");

}

void Cgrain::READ(FILE *fp)
{
  fscanf(fp,"%d",&ID);
  X.READ2d(fp);
  V.READ2d(fp);
  Xtrue.READ2d(fp);
//  Xtrue*=0;
  fscanf(fp,"%lf %lf %lf %lf",&Ome, &mass,&Imass,&R);
//  A.READ(fp);

}

void Cgrain::PRINT_pp(FILE *fp)
{

    stress.PRINT2d(fp); //1-4
    gradV.PRINT2d(fp); //5-8
    V.PRINT(fp); // 9-11

    fprintf(fp,"%f\t%f\t%f\t",voro_area,voro_area_n,solid_area_n); //12-14
    for (int i=0; i<12; i++) { //15-26
        if (i<voro_neighbour.size()) fprintf(fp,"%d\t", (int) (voro_neighbour[i]->ID));
        else fprintf(fp,"%d\t",-1);
    }

    fprintf(fp,"\n");

}

void Cgrain::READ_pp(FILE *fp)
{

    stress.READ2d(fp); //1-4
    gradV.READ2d(fp); //5-8
    V.READ2d(fp); // 9-12

    fscanf(fp,"%lf %lf %lf",voro_area,voro_area_n,solid_area_n); // 12-14
    for (int i=0; i<12; i++) { //15-26
        int neigh_id=0;
        fscanf(fp,"%d",neigh_id);
        cout<<neigh_id<<endl;
//        if (neigh_id > 0)
//            voro_neighbour.push_back(neigh_id);
    }

}

void Cgrain::convert_unit(double pix_per_d,double fps)
{
    X/=pix_per_d;
    V/=(pix_per_d*100);
    V*=fps;
    R/=pix_per_d;
}


