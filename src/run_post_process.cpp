#include "run_post_process.h"
#include "voro/voro++.hh"

using namespace voro;



Crun_post_process::Crun_post_process()
{
}



void Crun_post_process::evale_profile( )
{

  in_files.set_path("READ_WRITE");
  int iLast;
  int iFirst=in_files.iSave;
  cout<<"Enter the last file to be read (LAST_FILE)"<<endl;  in_files.get_secure("LAST_FILE",iLast);

  int lvl;
  cout<<"Enter no. of levels to compute Vel Grad (Eg. 1 = ~3d, 2 = ~5d, etc.) (key LEVEL)"<<endl; in_files.get_secure("LEVEL",lvl);


  config.read(in_files);  //init the profile average in time
  //Cvector L=config.cell.L;
/*
  Cmesh mesh_time_average(L,UNIT_LENGHT/1.);
  int iFirst=in_files.iSave;

  for( int i=iFirst; i<=iLast;i++)
    {
         config.read(in_files);


         for(int i=0;i<config.grain.size();i++)config.grain[i].evale_stress();


         Cmesh mesh(L,UNIT_LENGHT/1.);
         mesh.get_avergage(config.grain);

         mesh_time_average+=mesh;

         in_files.set_file_name();
         mesh.save(in_files.file_profile,in_files.file_mesh);
         in_files.iSave++;
  }

  //mesh_time_average/=(iLast-iFirst);
  for(int y=0;y<mesh_time_average.Yprofile.size();y++)
      mesh_time_average.Yprofile[y]/=mesh_time_average.Yprofile[y].weight;

  mesh_time_average.save(in_files.path[4]+"/profile_average",in_files.path[5]+"/mesh_average");
*/
/*
   QList <Cmatrix > stress;
   QList <Cvector> L;
    QList <double> phi;



   for( int i=iFirst; i<=iLast;i++)
    {
       config.read(in_files);
       config.evale_stress();

       stress.push_back(config.cell.stress);
       L.push_back(config.cell.L);
       phi.push_back(config.evale_solid_fraction());

       in_files.iSave++;
   }

    //temporal average
   Cmatrix stress_average;
   Cvector L_average;
   double phi_average=0;
   for(int i=0;i<stress.size();i++)
   {
       stress_average+=stress[i];
       L_average+=L[i];
       phi_average+=phi[i];
   }

   stress_average/=stress.size();
   L_average/=L.size();
   phi_average/=phi.size();

   FILE *fp;

   string file_name= "mu_I_H";
   fp=fopen(file_name.c_str(),"a");

   fprintf(fp,"%lf \t %lf \t %lf \t %lf \n",config.parameter.shear_rate/sqrt(fabs(stress_average.x[1][1])),stress_average.x[0][1]/fabs(stress_average.x[1][1]),phi_average,L_average.x[1]);
   fclose(fp);

   exit(0);

*/


  double H=config.cell.L[1];
  double num_slice_per_grain = 2.;
  Cprofile profile_average(H,config.cell.L[0], UNIT_LENGHT/num_slice_per_grain);

 // int iFirst=in_files.iSave;

  for( int i=iFirst; i<=iLast;i++)
  {
      config.read(in_files);

      // for dt
      Cconfig conf_next;
      in_files.iSave++;
      in_files.set_file_name();
      conf_next.read(in_files);
      in_files.iSave--;
      in_files.set_file_name();

       for(int i=0;i<config.grain.size();i++)
       {
          Cvector dX = conf_next.grain[i].X-config.grain[i].X;
          config.cell.CLP(dX);
//          config.cell.CLP(dX);
          config.grain[i].V=dX/(conf_next.t-config.t);
      }


      //// VORONOI
      // create a container for voronoi. Periodic in x
      container_poly con(0,config.cell.L.x[0],-2*config.cell.L.x[1]/2,2*config.cell.L.x[1]/2,-0.5,0.5,
              (int) (config.cell.L.x[0]/2),(int) (config.cell.L.x[1]/2),1,true,false,false,8);
      particle_order po;

      // add all grains
      for(int i=0;i<config.grain.size();i++)
      {
          // add ith grain to the container and to voroids
          con.put(po,config.grain[i].ID,config.grain[i].X.x[0],config.grain[i].X.x[1],config.grain[i].X.x[2],config.grain[i].R);

          // add top periodic grain if 5 grain thick, i.e. y < (-H/2+5)
          if (config.grain[i].X.x[1]<(-(config.cell.L.x[1]/2-6)))
          {
              // add id by 1*10^6 for top periodic
              double new_x = config.grain[i].X.x[0]+config.cell.shift;
              if (new_x >= config.cell.L.x[0]) new_x-=config.cell.L.x[0];
              if (new_x < 0) new_x+=config.cell.L.x[0];
              con.put(po,config.grain[i].ID+1000000,new_x,config.grain[i].X.x[1]+config.cell.L.x[1],config.grain[i].X.x[2],config.grain[i].R);
          }

          // add bottom periodic grain if 5 grain thick, i.e. y > (H/2-5)
          if (config.grain[i].X.x[1]>(config.cell.L.x[1]/2-6))
          {
              // add id by 2*10^6 for top periodic
              double new_x = config.grain[i].X.x[0]-config.cell.shift;
              if (new_x >= config.cell.L.x[0]) new_x-=config.cell.L.x[0];
              if (new_x < 0) new_x+=config.cell.L.x[0];
              con.put(po,config.grain[i].ID+2000000,new_x,config.grain[i].X.x[1]-config.cell.L.x[1],config.grain[i].X.x[2],config.grain[i].R);
          }
      }

      cout << "all grains read" <<endl;


      // create loop and cell
      c_loop_order clo(con,po);
      voronoicell_neighbor c;

      // loop through each voronoicells c
      if(clo.start()) do if(con.compute_cell(c,clo)) {
          int this_grain_id = clo.pid();

          // to avoid repetation on periodic grains
          if (this_grain_id < 1000000)
          {
              // calculate and save the volume of this grain's voronoi cell

              config.grain[this_grain_id].voro_area = c.volume();

              // calculate and save the neighbouring grains as per voronoi
              std::vector <int> nvector; //vector of neighoubrs
              c.neighbors(nvector);

              // calculate face areas
              std::vector <double> nfarea; // neighbouring face area
              c.face_areas(nfarea);

              for (int i=0;i<nvector.size();i++)
              {
                  // do modulus of 10^6, in case the neibour is periodic grain
                  if ((int) nvector[i]>=0){
                      config.grain[this_grain_id].voro_neighbour.push_back(&config.grain[(int)(nvector[i])%1000000]);

                  }
              }
          }
      } while (clo.inc());


//        con.draw_cells_gnuplot("/users/prashidhakharel/pack_ten_cube_pp.gnu");

      // loop through each grain to compute values from voro
      for(int i=0;i<config.grain.size();i++)config.grain[i].evale_from_voro(lvl);

      ///////////////////////////


//       for(int i=0;i<config.grain.size();i++)config.grain[i].evale_stress();
//       for(int i=0;i<config.grain.size();i++)config.grain[i].evale_shear_rate(config.grain);

      config.save_grain_postprocess(in_files);
//      in_files.iSave++;


      Cprofile profile(H,config.cell.L[0], UNIT_LENGHT/num_slice_per_grain);

      profile.make_profile(config.grain);
      profile_average+=profile;
      in_files.set_file_name();
      profile.save(in_files.file_profile);
      in_files.iSave++;
      config.read_grain_postprocess(in_files);

     config.evale_stress();

    }

  profile_average/=(iLast-iFirst);

  profile_average.average_y(&config.cell);

//  for(int i=1;i<profile_average.slice.size()-1;i++)profile_average.slice[i].gradV.x[0][1] = (profile_average.slice[i+1].V[0]-profile_average.slice[i-1].V[0])/(2.*profile_average.dy);


  cout<<in_files.path[4];
  profile_average.save(in_files.path[4]+"/profile_average");
  profile_average.save("./profile");

}

void Crun_post_process::evale_voro()
{
    // inputs
    in_files.set_path("READ_WRITE");
    int iLast;
      int iFirst=in_files.iSave;
    cout<<"Enter the last file to be read (LAST_FILE)"<<endl;  in_files.get_secure("LAST_FILE",iLast);
    int sav_int;
    cout<<"Enter saving interval for files. For default enter 1 (key SAV_INT)"<<endl; in_files.get_secure("SAV_INT",sav_int);

    int lvl;
    cout<<"Enter no. of levels to compute Vel Grad (Eg. 1 = ~3d, 2 = ~5d, etc.) (key LEVEL)"<<endl; in_files.get_secure("LEVEL",lvl);


    int nfram;
    cout<<"No of frames used to calculate velocity [-1 for instantanious; n for number of fram] (key NFRAM)"<<endl; in_files.get_secure("NFRAM",nfram);

    for( int cur_file=iFirst; cur_file<=iLast;cur_file+=sav_int)
    {

        config.read(in_files);  //init the profile average in time

        cout<<"first grain: "<<config.grain[10].ID<<endl;

        // if not instantanious then calc vel from two frames
        if (nfram > 0) {
            Cconfig conf_prev;
            Cconfig conf_next;
            in_files.iSave+=nfram;
            in_files.set_file_name();
            conf_next.read(in_files);
            in_files.iSave-=nfram;
            in_files.set_file_name();

            in_files.iSave-=nfram;
            in_files.set_file_name();
            conf_prev.read(in_files);
            in_files.iSave+=nfram;
            in_files.set_file_name();

             for(int i=0;i<config.grain.size();i++)
             {
                Cvector dX = conf_next.grain[i].X-conf_prev.grain[i].X;

                //config.cell.CLP(dX);
                // use next config cell's CLP

                conf_next.cell.CLP(dX);




                // config.cell.CLP(dX);

                config.grain[i].V = dX/(conf_next.t-conf_prev.t);
                if (conf_prev.grain[i].X.x[1]-config.grain[i].X.x[1]>config.cell.L.x[1]/2) {
                    config.grain[i].V -= config.cell.V;
                } else if (conf_prev.grain[i].X.x[1]-config.grain[i].X.x[1]<-(1*config.cell.L.x[1]/2)) {
                    config.grain[i].V += config.cell.V;
                }

            }
        }


        cout << "first file read"<<endl;

        in_files.set_file_name();


        // Divides each grain to corresponding voronoi cells
        // then add volume of each cell and voronoi neighbours to each grain

        // create a container for voronoi. Periodic in x
        container_poly con(0,config.cell.L.x[0],-2*config.cell.L.x[1]/2,2*config.cell.L.x[1]/2,-0.5,0.5,
                (int) (config.cell.L.x[0]/2),(int) (config.cell.L.x[1]/2),1,true,false,false,8);
        particle_order po;

        // add all grains
        for(int i=0;i<config.grain.size();i++)
        {
            // add ith grain to the container and to voroids
            con.put(po,config.grain[i].ID,config.grain[i].X.x[0],config.grain[i].X.x[1],config.grain[i].X.x[2],config.grain[i].R);

            // add top periodic grain if 5 grain thick, i.e. y < (-H/2+5)
            if (config.grain[i].X.x[1]<(-(config.cell.L.x[1]/2-6)))
            {
                // add id by 1*10^6 for top periodic
                double new_x = config.grain[i].X.x[0]+config.cell.shift;
                if (new_x >= config.cell.L.x[0]) new_x-=config.cell.L.x[0];
                if (new_x < 0) new_x+=config.cell.L.x[0];
                con.put(po,config.grain[i].ID+1000000,new_x,config.grain[i].X.x[1]+config.cell.L.x[1],config.grain[i].X.x[2],config.grain[i].R);
            }

            // add bottom periodic grain if 5 grain thick, i.e. y > (H/2-5)
            if (config.grain[i].X.x[1]>(config.cell.L.x[1]/2-6))
            {
                // add id by 2*10^6 for top periodic
                double new_x = config.grain[i].X.x[0]-config.cell.shift;
                if (new_x >= config.cell.L.x[0]) new_x-=config.cell.L.x[0];
                if (new_x < 0) new_x+=config.cell.L.x[0];
                con.put(po,config.grain[i].ID+2000000,new_x,config.grain[i].X.x[1]-config.cell.L.x[1],config.grain[i].X.x[2],config.grain[i].R);
            }
        }

        cout << "all grains read" <<endl;


        // create loop and cell
        c_loop_order clo(con,po);
        voronoicell_neighbor c;

        // loop through each voronoicells c
        if(clo.start()) do if(con.compute_cell(c,clo)) {
            int this_grain_id = clo.pid();

            // to avoid repetation on periodic grains
            if (this_grain_id < 1000000)
            {
                // calculate and save the volume of this grain's voronoi cell

                config.grain[this_grain_id].voro_area = c.volume();

                // calculate and save the neighbouring grains as per voronoi
                std::vector <int> nvector; //vector of neighoubrs
                c.neighbors(nvector);

                // calculate face areas
                std::vector <double> nfarea; // neighbouring face area
                c.face_areas(nfarea);

                for (int i=0;i<nvector.size();i++)
                {
                    // do modulus of 10^6, in case the neibour is periodic grain
                    if ((int) nvector[i]>=0){
                        config.grain[this_grain_id].voro_neighbour.push_back(&config.grain[(int)(nvector[i])%1000000]);

                    }
                }
            }
        } while (clo.inc());

//        con.draw_cells_gnuplot("/users/prashidhakharel/pack_ten_cube_pp.gnu");

        // loop through each grain to compute values from voro
        for(int i=0;i<config.grain.size();i++)config.grain[i].evale_from_voro(lvl);
        config.save_grain_postprocess(in_files);

        in_files.set_file_name();

        in_files.iSave+=sav_int;// increase infiles

    }
}

double gaussian(double x, double wdth){
    return 1.0/sqrt(2*PI*wdth*wdth)*exp(-(x*x)/(2*wdth*wdth));
}

void Crun_post_process::evale_corr() {
    // inputs
    in_files.set_path("READ_WRITE");
    int iLast;
    int iFirst=in_files.iSave;
    cout<<"Enter the last file to be read (LAST_FILE)"<<endl;  in_files.get_secure("LAST_FILE",iLast);

    double max_corr_len;
    cout<<"Enter maximum correlation lengh (MAX_CORR_LEN)"<<endl; in_files.get_secure("MAX_CORR_LEN",max_corr_len);

    double corr_len_intvl;
    cout<<"Enter interval for correlation length (CORR_LEN_INTVL)"<<endl; in_files.get_secure("CORR_LEN_INTVL",corr_len_intvl);

    // loop to get corr for various lengths
    QList <double> corr_len;
    QList <double> corr_Fxy_x;
    QList <double> corr_Fxy_y;
    QList <double> corr_Fyx_x;
    QList <double> corr_Fyx_y;

    for (double c_len = 0; c_len <= max_corr_len; c_len+=corr_len_intvl) {
        double num_Fxy_x = 0; // numerator
        double den_Fxy_x = 0; // denominator
        double num_Fxy_y = 0;
        double den_Fxy_y = 0;
        double num_Fyx_x = 0;
        double den_Fyx_x = 0;
        double num_Fyx_y = 0;
        double den_Fyx_y = 0;

        // loop through different time
        for( int cur_file=iFirst; cur_file<=iLast;cur_file++)
        {
            in_files.iSave = cur_file;
            in_files.set_file_name();

            cout<<"Read config"<<endl;
            config.read(in_files);
            cout<<"Read pp config"<<endl;
            config.read_grain_postprocess(in_files);

            cout << "This file read"<<endl;

            in_files.set_file_name();

            cout << "File name set, about to loop" <<endl;

            // loop through each grain twice
            foreach (Cgrain gi, config.grain)
            foreach (Cgrain gj, config.grain) {
//                cout<<"inside loop"<<endl;
                Cvector dx = gi.X - gj.X;
                config.cell.CLP(dx);

                double g_width = 0.4; //gaussian width;

                num_Fxy_x += gi.gradV.x[1][0]*gj.gradV.x[1][0]\
                        *gaussian(dx[0] + c_len, g_width)\
                        *gaussian(dx[1], g_width);
                den_Fxy_x += gaussian(dx[0] + c_len,g_width)\
                        *gaussian(dx[1],g_width);

                num_Fxy_y += gi.gradV.x[1][0]*gj.gradV.x[1][0]\
                        *gaussian(dx[0],g_width)\
                        *gaussian(dx[1] + c_len,g_width);
                den_Fxy_y += gaussian(dx[0],g_width)\
                        *gaussian(dx[1] + c_len,g_width);

                num_Fyx_x += gi.gradV.x[0][1]*gj.gradV.x[0][1]\
                        *gaussian(dx[0] + c_len,g_width)\
                        *gaussian(dx[1],g_width);
                den_Fyx_x += gaussian(dx[0] + c_len,g_width)\
                        *gaussian(dx[1],g_width);

                num_Fyx_y += gi.gradV.x[0][1]*gj.gradV.x[0][1]\
                        *gaussian(dx[0],g_width)\
                        *gaussian(dx[1] + c_len,g_width);
                den_Fyx_y += gaussian(dx[0],g_width)\
                        *gaussian(dx[1] + c_len,g_width);
            } // each grains

            cout<<"curr file:"<<cur_file<<endl;

            cout<<"set file name"<<endl;

            cout<< "increase iSave"<<endl;
        }// each time

        cout<<"Done c_len: "<<c_len<<endl;

        corr_len.push_back(c_len);
        corr_Fxy_x.push_back(num_Fxy_x/den_Fxy_x);
        corr_Fxy_y.push_back(num_Fxy_y/den_Fxy_y);
        corr_Fyx_x.push_back(num_Fyx_x/den_Fyx_x);
        corr_Fyx_y.push_back(num_Fyx_y/den_Fyx_y);
    } // each corr lengths


    // Write
    cout<<"Before writing"<<endl;
    FILE *fp;
    fp = fopen((in_files.path_root+"corr").c_str(),"w");

    cout<<"Write loc: "<<(in_files.path_root+"corr").c_str()<<endl;

    for (int i = 0; i<corr_len.size(); i++) {
        fprintf(fp,"%f\t%f\t%f\t%f\t%f\n",corr_len[i],corr_Fxy_x[i],corr_Fxy_y[i],corr_Fyx_x[i],corr_Fyx_y[i]);
    }
    cout<<"After writing"<<endl;
    fclose(fp);
}

void Crun_post_process::evale_length() {

    // Read/write file paths
    in_files.set_path("READ_WRITE");
    int iFirst=in_files.iSave;
    int iLast;
    cout<<"Enter the last file to be read (LAST_FILE)"<<endl;  in_files.get_secure("LAST_FILE",iLast);

    // Initialise squared velocity and vorticity fluctuations
    double dv = 0, dw = 0;

    // For each file number
    for(int i=iFirst; i<=iLast; i++) {

        // Read config
        config.read(in_files);
        in_files.set_file_name();

        // Compute velocity fluctuation
        for(int i=0;i<config.grain.size();i++) {dv += pow(config.grain[i].V[0] - config.parameter.shear_rate*config.grain[i].X[1], 2);}
        dv /= static_cast<double>(config.grain.size());

        // Compute vorticity fluctuation
        for(int i=0;i<config.grain.size();i++) {dw += pow(config.grain[i].vorticity - config.parameter.shear_rate/2, 2);}

        in_files.iSave++;
    }

    // Average of the number of files
//    dv /= (iLast-iFirst+1); dw /= (iLast-iFirst+1);

    // Extract length scale
//    double l = sqrt(dv / dw);
//    cout << l << endl;

}
