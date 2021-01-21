#include "config.h"


Cconfig::Cconfig()
{
  cell.para=&parameter;
}

void Cconfig::evolve()
{
//  Ctimer time_iter;
 // time_iter.start();

//    for (int i=0;i<grain.size();i++) {
//        cout << "time: "<<t<<endl;
//        cout << "Velocity x: "<<grain[i].V.x[0]<<endl;
//        cout << "Velocity y: "<<grain[i].V.x[1]<<endl;
//    }


  predictor();
  refresh_contact(false);
  sum_force();
  corrector();


 // time_iter.end();
  //cout<<endl<<"Corrector duration:\t "<<time_iter.duration() <<endl<<endl;
}


void Cconfig::predictor( )
{
  cell.predictor(dt);
  for(int i=0;i<grain.size();i++)
    {
      grain[i].predictor(dt);

    }
}


void Cconfig::refresh_contact(bool force_check_neighbours)
{
  /*<< Use the grid to find the neighbour*/
  if(check_neighbours.should_do(t)||force_check_neighbours)
    {
  for(int i=0;i<grain.size();i++)       grain[i].neighbour.clear();

    for(int i=0;i<grain.size()-1;i++)
        for(int j=i+1;j<grain.size();j++)
            {
            Cvector dX = grain[j].X-grain[i].X;
            cell.CLP(dX);
            if(dX.norm()<3.5*UNIT_LENGHT)
                grain[i].neighbour.push_back(&grain[j]);
            }
    }
  for(int i=0;i<grain.size();i++) grain[i].update_contact(dt);

//  cout<<"before refresh"<<endl;
//  cell.stress.PRINT();
//  evale_stress();
//  cell.stress.PRINT();

};



void Cconfig::sum_force()
{
  for(int i=0;i<grain.size();i++)
  {
      grain[i].Fsum*=0;//set to zero
      grain[i].Msum=0;
  }

  // CELL
  cell.stress*=0;
  cell.boundary_normal_stress*=0;
  cell.boundary_shear_stress*=0;

  foreach(Cgrain g, grain)
    foreach(Ccontact c, g.contact)
    {
      c.A->Fsum +=c.Force;
      c.B->Fsum -=c.Force;
      c.A->Msum += c.ft*c.A->R;
      c.B->Msum += c.ft*c.B->R;
      cell.stress+=(c.Force^c.nAB)*(c.A->R+c.B->R+c.delta);
      if(c.AM_I_CROSS_BOUNDARY()){
          if (c.B->X.x[1]-c.A->X.x[1]<-cell.L.x[1]/2) {
              cell.boundary_normal_stress+=c.Force.x[1];
              cell.boundary_shear_stress+=c.Force.x[0];
          }else if(c.B->X.x[1]-c.A->X.x[1]>cell.L.x[1]/2){
              cell.boundary_normal_stress-=c.Force.x[1];
              cell.boundary_shear_stress-=c.Force.x[0];
          }

      }
    }

  cell.stress/=cell.L[0]*cell.L[1]*UNIT_LENGHT;

//  cell.stress.PRINT();
  //cell.stress_boundary/=cell.L[0]*2*UNIT_LENGHT*UNIT_LENGHT;
  cell.boundary_normal_stress/=cell.L[0]; // force per unit length
  cell.boundary_shear_stress/=cell.L[0]; // force per unit length
}

void Cconfig::corrector( )
{
  cell.corrector(dt);
  for(int i=0;i<grain.size();i++) grain[i].corrector(dt);
};

void Cconfig::Vy_sum_corrector()
{
    // Calculate yvelocity error
    double vy_sum = 0;
    double vx_sum = 0;
    double mass_sum = 0;
    for (int i=0;i<grain.size();i++) {
//        vx_sum += grain[i].V.x[0]*grain[i].mass;
        vy_sum += grain[i].V.x[1]*grain[i].mass;
        mass_sum += grain[i].mass;
    }

    // Redistribute error to all grains
    for (int i=0;i<grain.size();i++) {
//        grain[i].V.x[0]-=vx_sum/mass_sum;
        grain[i].V.x[1]-=vy_sum/mass_sum;
    }

    if (abs(t - (int) t) < dt) {
        cout << "Vy_sum = " << vy_sum << endl;
    }
}


void Cconfig::evale_stress()
{
 cell.stress*=0;

foreach(Cgrain g, grain)
    foreach(Ccontact c, g.contact)
   {
    c.get_distance();
//    cell.stress+=(c.Force^c.nAB)*(c.A->R+c.B->R)*(g.R)/(c.A->R+c.B->R);
    cell.stress+=(c.Force^c.nAB)*(c.A->R+c.B->R+c.delta);
}
cell.stress/=cell.L.x[0]*cell.L.x[1]*cell.L.x[2];

foreach(Cgrain g, grain) g.evale_stress();

//cout<<stress.x[1][1]<<"\t"<<stress.x[0][1]<<endl;

}


void Cconfig::set_linear_velocity_profile()
{
    for(int i=0;i<grain.size();i++) {
        grain[i].V[0]=grain[i].X[1]*parameter.shear_rate;
        grain[i].V[1]=0;
    }
}





void Cconfig::save(Cin_out out_file)
{
  out_file.set_file_name();

  FILE *fp;
  fp=fopen(out_file.file_grain.c_str(),"w");
  if(fp==NULL){string err="Can not write in '" + out_file.file_grain+"' \n this file does not exist\n";  out_file.error("Config.cpp","save",err);}
  for(int i=0;i<grain.size();i++) grain[i].PRINT(fp);
  fclose(fp);

  fp=fopen(out_file.file_contact.c_str(),"w");
  if(fp==NULL){string err="Can not write in '" + out_file.file_contact+"' \n this file does not exist\n";  out_file.error("Config.cpp","save",err);}
  for(int i=0;i<grain.size();i++)
      for(int c=0;c<grain[i].contact.size();c++)
          grain[i].contact[c].PRINT(fp);

  fclose(fp);

  fp=fopen(out_file.file_cell.c_str(),"w");
  if(fp==NULL){string err="Can not write in '" + out_file.file_cell+"' \n this file does not exist\n";  out_file.error("Config.cpp","save",err);}
  cell.PRINT(fp,t);
  fclose(fp);

  fp=fopen(out_file.file_parameter.c_str(),"w");
  if(fp==NULL){string err="Can not write in '" + out_file.file_parameter+"' \n this file does not exist\n";  out_file.error("Config.cpp","save",err);}
  parameter.PRINT(fp);
  fclose(fp);

}

void Cconfig::save_dt_grains(Cin_out out_file)
{
    // Save each grain time,X,Y,dx,dy,H,shift only at finer time interval
    //out_file.set_file_name();
    // Save grain info in mesh folder ,i.e. 5

    // Save x position file
    string curr_file_nm = out_file.path[5]+"/grain_t_xs_"+out_file.to_string(out_file.iSave-1);
    FILE *fp;
    fp=fopen(curr_file_nm.c_str(),"a");
    fprintf(fp,"%.10e\t",t); //1: time
    for(int i=0;i<grain.size();i++) {
        fprintf(fp,"%.10e\t",grain[i].X.x[0]);
    }
    fprintf(fp,"\n");
    fclose(fp);

    // Save y position file
    curr_file_nm = out_file.path[5]+"/grain_t_ys_"+out_file.to_string(out_file.iSave-1);
    fp=fopen(curr_file_nm.c_str(),"a");
    fprintf(fp,"%.10e\t",t); //1: time
    for(int i=0;i<grain.size();i++) {
        fprintf(fp,"%.10e\t",grain[i].X.x[1]);
    }
    fprintf(fp,"\n");
    fclose(fp);

    // Save x vel file
    curr_file_nm = out_file.path[5]+"/grain_t_vxs_"+out_file.to_string(out_file.iSave-1);
    fp=fopen(curr_file_nm.c_str(),"a");
    fprintf(fp,"%.10e\t",t); //1: time
    for(int i=0;i<grain.size();i++) {
        fprintf(fp,"%.10e\t",grain[i].V.x[0]);
    }
    fprintf(fp,"\n");
    fclose(fp);

    // Save y vel file
    curr_file_nm = out_file.path[5]+"/grain_t_vys_"+out_file.to_string(out_file.iSave-1);
    fp=fopen(curr_file_nm.c_str(),"a");
    fprintf(fp,"%.10e\t",t); //1: time
    for(int i=0;i<grain.size();i++) {
        fprintf(fp,"%.10e\t",grain[i].V.x[1]);
    }
    fprintf(fp,"\n");
    fclose(fp);

    // Save xx stress file
    curr_file_nm = out_file.path[5]+"/grain_dt_stress_xx_"+out_file.to_string(out_file.iSave-1);
//    FILE *fp;
    fp=fopen(curr_file_nm.c_str(),"a");
    fprintf(fp,"%.10e\t",t); //1: time
    for(int i=0;i<(0.1*grain.size());i++) {
        fprintf(fp,"%.10e\t",grain[i].stress.x[0][0]);
    }
    fprintf(fp,"\n");
    fclose(fp);

    // Save xy stress file
    curr_file_nm = out_file.path[5]+"/grain_dt_stress_xy_"+out_file.to_string(out_file.iSave-1);
//    FILE *fp;
    fp=fopen(curr_file_nm.c_str(),"a");
    fprintf(fp,"%.10e\t",t); //1: time
    for(int i=0;i<(0.1*grain.size());i++) {
        fprintf(fp,"%.10e\t",grain[i].stress.x[0][1]);
    }
    fprintf(fp,"\n");
    fclose(fp);

    // Save yx stress file
    curr_file_nm = out_file.path[5]+"/grain_dt_stress_yx_"+out_file.to_string(out_file.iSave-1);
//    FILE *fp;
    fp=fopen(curr_file_nm.c_str(),"a");
    fprintf(fp,"%.10e\t",t); //1: time
    for(int i=0;i<(0.1*grain.size());i++) {
        fprintf(fp,"%.10e\t",grain[i].stress.x[1][0]);
    }
    fprintf(fp,"\n");
    fclose(fp);

    // Save yy stress file
    curr_file_nm = out_file.path[5]+"/grain_dt_stress_yy_"+out_file.to_string(out_file.iSave-1);
//    FILE *fp;
    fp=fopen(curr_file_nm.c_str(),"a");
    fprintf(fp,"%.10e\t",t); //1: time
    for(int i=0;i<(0.1*grain.size());i++) {
        fprintf(fp,"%.10e\t",grain[i].stress.x[1][1]);
    }
    fprintf(fp,"\n");
    fclose(fp);

    // dx and dy position not required
//    // Save dx position file
//    curr_file_nm = out_file.path[5]+"/grain_t_dxs_"+out_file.to_string(out_file.iSave-1);

//    fp=fopen(curr_file_nm.c_str(),"a");
//    fprintf(fp,"%.6e\t",t); //1: time
//    for(int i=0;i<grain.size();i++) {
//        fprintf(fp,"%d\t",(int) grain[i].Xtrue.x[0]);
//    }
//    fprintf(fp,"\n");
//    fclose(fp);

//    // Save dy position file
//     curr_file_nm = out_file.path[5]+"/grain_t_dys_"+out_file.to_string(out_file.iSave-1);

//    fp=fopen(curr_file_nm.c_str(),"a");
//    fprintf(fp,"%.6e\t",t); //1: time
//    for(int i=0;i<grain.size();i++) {
//        fprintf(fp,"%d\t",(int) grain[i].Xtrue.x[1]);
//    }
//    fprintf(fp,"\n");
//    fclose(fp);

    // Save H and shift position file
    curr_file_nm = out_file.path[5]+"/cell_t_H_shift_"+out_file.to_string(out_file.iSave-1);

    fp=fopen(curr_file_nm.c_str(),"a");
    fprintf(fp,"%.10e\t%.10e\t%.10e\t",t,cell.L.x[1],cell.shift); //1: time
    fprintf(fp,"\n");
    fclose(fp);
}


void Cconfig::save_dt_forces(Cin_out out_file)
{
    // Save each grain time,X,Y,dx,dy,H,shift only at finer time interval
    //out_file.set_file_name();
    // Save grain info in mesh folder ,i.e. 5

    // Save x force file
    string curr_file_nm = out_file.path[5]+"/grain_dt_forces_xs";
    FILE *fp;
    fp=fopen(curr_file_nm.c_str(),"a");
    fprintf(fp,"%.10e\t",t); //1: time
    for(int i=0;i<(0.1*grain.size());i++) {
        fprintf(fp,"%.10e\t",grain[i].Fsum.x[0]);
    }
    fprintf(fp,"\n");
    fclose(fp);

    // Save y force file
    curr_file_nm = out_file.path[5]+"/grain_dt_forces_ys";
//    FILE *fp;
    fp=fopen(curr_file_nm.c_str(),"a");
    fprintf(fp,"%.10e\t",t); //1: time
    for(int i=0;i<(0.1*grain.size());i++) {
        fprintf(fp,"%.10e\t",grain[i].Fsum.x[1]);
    }
    fprintf(fp,"\n");
    fclose(fp);

    // Save x position file
    curr_file_nm = out_file.path[5]+"/grain_dt_xs";
//    FILE *fp;
    fp=fopen(curr_file_nm.c_str(),"a");
    fprintf(fp,"%.10e\t",t); //1: time
    for(int i=0;i<(0.1*grain.size());i++) {
        fprintf(fp,"%.10e\t",grain[i].X.x[0]);
    }
    fprintf(fp,"\n");
    fclose(fp);

    // Save y position file
    curr_file_nm = out_file.path[5]+"/grain_dt_ys";
//    FILE *fp;
    fp=fopen(curr_file_nm.c_str(),"a");
    fprintf(fp,"%.10e\t",t); //1: time
    for(int i=0;i<(0.1*grain.size());i++) {
        fprintf(fp,"%.10e\t",grain[i].X.x[1]);
    }
    fprintf(fp,"\n");
    fclose(fp);

    // Save x vel file
    curr_file_nm = out_file.path[5]+"/grain_dt_vel_xs";
//    FILE *fp;
    fp=fopen(curr_file_nm.c_str(),"a");
    fprintf(fp,"%.10e\t",t); //1: time
    for(int i=0;i<(0.1*grain.size());i++) {
        fprintf(fp,"%.10e\t",grain[i].V.x[0]);
    }
    fprintf(fp,"\n");
    fclose(fp);

    // Save y vel file
    curr_file_nm = out_file.path[5]+"/grain_dt_vel_ys";
//    FILE *fp;
    fp=fopen(curr_file_nm.c_str(),"a");
    fprintf(fp,"%.10e\t",t); //1: time
    for(int i=0;i<(0.1*grain.size());i++) {
        fprintf(fp,"%.10e\t",grain[i].V.x[1]);
    }
    fprintf(fp,"\n");
    fclose(fp);

    // Save xx stress file
    curr_file_nm = out_file.path[5]+"/grain_dt_stress_xx";
//    FILE *fp;
    fp=fopen(curr_file_nm.c_str(),"a");
    fprintf(fp,"%.10e\t",t); //1: time
    for(int i=0;i<(0.1*grain.size());i++) {
        fprintf(fp,"%.10e\t",grain[i].stress.x[0][0]);
    }
    fprintf(fp,"\n");
    fclose(fp);

    // Save xy stress file
    curr_file_nm = out_file.path[5]+"/grain_dt_stress_xy";
//    FILE *fp;
    fp=fopen(curr_file_nm.c_str(),"a");
    fprintf(fp,"%.10e\t",t); //1: time
    for(int i=0;i<(0.1*grain.size());i++) {
        fprintf(fp,"%.10e\t",grain[i].stress.x[0][1]);
    }
    fprintf(fp,"\n");
    fclose(fp);

    // Save yx stress file
    curr_file_nm = out_file.path[5]+"/grain_dt_stress_yx";
//    FILE *fp;
    fp=fopen(curr_file_nm.c_str(),"a");
    fprintf(fp,"%.10e\t",t); //1: time
    for(int i=0;i<(0.1*grain.size());i++) {
        fprintf(fp,"%.10e\t",grain[i].stress.x[1][0]);
    }
    fprintf(fp,"\n");
    fclose(fp);

    // Save yy stress file
    curr_file_nm = out_file.path[5]+"/grain_dt_stress_yy";
//    FILE *fp;
    fp=fopen(curr_file_nm.c_str(),"a");
    fprintf(fp,"%.10e\t",t); //1: time
    for(int i=0;i<(0.1*grain.size());i++) {
        fprintf(fp,"%.10e\t",grain[i].stress.x[1][1]);
    }
    fprintf(fp,"\n");
    fclose(fp);


}


void Cconfig::save_grain_postprocess(Cin_out in_file )
{
    FILE *fp;
    fp=fopen(in_file.file_grain_pp.c_str(),"w");
    if(fp==NULL) return;

    for(int i=0;i<grain.size();i++)
        grain[i].PRINT_pp(fp);
    fclose(fp);
}

void Cconfig::read_grain_postprocess(Cin_out in_file )
{
    in_file.set_file_name();
    FILE *fp;
    cout<<"pp path: "<<in_file.file_grain_pp.c_str()<<endl;
    fp=fopen(in_file.file_grain_pp.c_str(),"r");
    if(fp==NULL) return;

    cout<<"Before pp loop"<<endl;

    for(int i=0;i<grain.size();i++){
//        grain[i].READ_pp(fp);
        grain[i].stress.READ(fp); // 1-9
        grain[i].gradV.READ(fp); // 10-18

        fscanf(fp,"%lf %lf %lf",&grain[i].voro_area,&grain[i].voro_area_n,&grain[i].solid_area_n); // 19-21
        for (int i=0; i<12; i++) { //22-33
            int neigh_id;
            fscanf(fp,"%d ",&neigh_id);
            if (neigh_id > 0) grain[i].voro_neighbour.push_back(&grain[neigh_id]);
        }
    }
    cout<<"Before cloasing"<<endl;
    fclose(fp);
}


void Cconfig::read(Cin_out in_file)
{
  in_file.set_file_name();

  cout<<"Start read config from "<<in_file.file_grain<<"\t";

  FILE *fp;

  // READ GRAINS
  fp=fopen(in_file.file_grain.c_str(),"r");
  if(fp==NULL){string err="Can not read in '" + in_file.file_grain+"' \n this file does not exist\n";  in_file.error("Config.cpp","read",err);}
  int Ngrain = in_file.get_Nline(in_file.file_grain);

  grain.clear();//clear pre existing grains if any
  for(int i=0;i<Ngrain;i++)
    {
      Cgrain g(&cell, &parameter);//set the pointers to the parameter and to the cell grains
      g.READ(fp);
      grain.push_back(g);

    }
  fclose(fp);

  // READ CONTACTS
  //read contact ; must be after reading grains
  fp=fopen( in_file.file_contact.c_str(),"r");
  if(fp==NULL){string err="Can not read in '" +  in_file.file_contact+"' \n this file does not exist\n";   in_file.error("Config.cpp","read",err);}
  int Ncontact =  in_file.get_Nline( in_file.file_contact);

  for(int c=0;c<Ncontact;c++)
    {
      Ccontact C;
      C.READ(fp);
      C.A = &grain[C.aID]; C.B = &grain[C.bID];//set the pointers to the two grains
      C.cell = &cell;      C.para = &parameter ;//set the pointers to the parameter and to the cell grains
//      grain[C.bID].contact.push_back(C);
      C.get_distance();
      C.ft = C.nABt*C.Force;
      grain[C.aID].contact.push_back(C);//store the contact in the A grain's list

  }
  fclose(fp);

   // READ CELL
  fp=fopen( in_file.file_cell.c_str(),"r");
  if(fp==NULL){string err="Can not read in '" +  in_file.file_cell+"' \n this file does not exist\n";   in_file.error("Config.cpp","read",err);}

  cell.READ(fp,t);
  fclose(fp);

  cell.evale_volume_grain(grain);

//  This one..................................
//  cell.stress.PRINT();
//  evale_stress();
//  cell.stress.PRINT();

  // READ PARAMETER
  fp=fopen( in_file.file_parameter.c_str(),"r");
  if(fp==NULL){string err="Can not read in '" +  in_file.file_parameter+"' \n this file does not exist\n";   in_file.error("Config.cpp","read",err);}

  parameter.READ(fp);
  fclose(fp);

//  // READ GRAIN POSTPROCESS
//  // read_grain_postprocess(in_file ) into grain, must be read after grains
//  fp=fopen(in_file.file_grain_pp.c_str(),"r");
//  if(fp==NULL){string err="Can not read in '" +  in_file.file_grain_pp+"' \n this file does not exist\n";   in_file.error("Config.cpp","read",err);}

//  for(int i=0;i<Ngrain;i++)
//  {
//    Cgrain g;
//    g.READ_pp(fp);
//    grain[i].gradV = g.gradV;
//    grain[i].stress = g.stress;

//  }
//  fclose(fp);

  cout<<"Success "<<endl;

  //grid.init(cell.L,UNIT_LENGHT*1.6);
  check_neighbours.next=t;
//  evolve();
  refresh_contact(true);
  sum_force();

}


