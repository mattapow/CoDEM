#include "run.h"
//#include <math.h>       /* ceil */ worked instead of unistd for Prashidha
#include <unistd.h>

Crun::Crun()
{
}



void Crun::simule()
{
  in_files.set_path("READ");
  out_files.set_path("WRITE");

  // read config
  config.read(in_files);

// initial save, not needed
//  config.save(out_files);

// cout<<"Enter the begining time (key TSTART)"<<endl;  in_files.get_secure("TSTART",tstart);
// cout<<"Enter the end time (key TEND)"<<endl;  in_files.get_secure("TEND",tend);


//  cout<<"Enter the begining time of saving (key SAVE_START)"<<endl;in_files.get_secure("SAVE_START",save_event.next);
// cout<<"Enter the period of time between saving (key SAVE_PERIOD)"<<endl;in_files.get_secure("SAVE_PERIOD",save_event.period);

  cout<<"Enter the shear rate (key SHEAR_RATE)"<<endl;in_files.get_secure("SHEAR_RATE",config.parameter.shear_rate);
  cout<<"Enter the cohesion (key COHESION)"<<endl;  in_files.get_secure("COHESION",config.parameter.Cohesion);

  double in_tend;
  cout<<"Enter no. of tend (key TEND)"<<endl;  in_files.get_secure("TEND",in_tend);
  double in_period;
  cout<<"Enter no. of delta t (key PERIOD)"<<endl; in_files.get_secure("PERIOD",in_period);

  int ndt;
  cout<<"Enter no. 1 if you want to save finer pos. data at mesh folder, no. 2 to save net force on 100 grains, 3 for both (or -1 otherwise)] (key NDT)"<<endl; in_files.get_secure("NDT",ndt);

  double fixh;
  cout<<"Fixed height? [if no enter -1 (or any negative value), if yes enter the height] (key FIXH)"<<endl; in_files.get_secure("FIXH",fixh);

  bool need_reset;
  if (fixh>0)  {
      need_reset = true;
  } else {
      need_reset = false;
  }
  double mu_p;
  cout<<"Friction? (key MU)"<<endl; in_files.get_secure("MU",mu_p);

  double K_value;
  cout<<"K value? (for default input -1) (key KVALUE)"<<endl; in_files.get_secure("KVALUE",K_value);

  double eta_value;
  cout<<"Eta value, i.e. viscous damping used to control system height in order to maintain pressure? (default 1, key ETA_VALUE)"<<endl; in_files.get_secure("ETA_VALUE",eta_value);

  config.parameter.mu = mu_p;

  if (K_value > 0){
      config.parameter.E = K_value;
      config.parameter.Stiff = K_value*UNIT_LENGHT;
  }

  config.cell.fixh = fixh;

  config.cell.eta = eta_value;

//  if (fixh > 0) {
//      config.cell.L.x[1] = fixh;
//      config.cell.V.x[1] = 0;
//  }


   Cevent save_event;
   tstart=0;
   tend=in_tend;
   save_event.next=0;

// **************************************************************


//   config.dt = sqrt(UNIT_MASS/(config.parameter.E*UNIT_LENGHT)) /20; //fraction of the collision time
//   config.dt = 1./1000;

   cout<<"------E="<<config.parameter.E<<endl;
   double no_collisions = 20;
//   config.dt = 1./pow(10,ceil(abs(log10(sqrt(UNIT_MASS/(config.parameter.E*UNIT_LENGHT))/20))));
   config.dt = floor(sqrt(UNIT_MASS/(config.parameter.E*UNIT_LENGHT)) /no_collisions *1000)/1000;
   if (config.dt == 0) {
       config.dt = floor(sqrt(UNIT_MASS/(config.parameter.E*UNIT_LENGHT)) /no_collisions *10000 )/10000;
   }
   if (config.dt == 0) {
       config.dt = sqrt(UNIT_MASS/(config.parameter.E*UNIT_LENGHT)) /no_collisions;
   }

   cout<<"selected dt: "<<config.dt<<endl;

   save_event.next=0;
   save_event.period=in_period;
   save_event.config_dt = config.dt;

   Cevent save_dt_grain_event;
   save_dt_grain_event.next=0;
   save_dt_grain_event.period = config.dt * ndt;
   save_dt_grain_event.config_dt = config.dt;

   // create event to track when to correct Vy sum
   Cevent Vy_sum_corrector_event;
   Vy_sum_corrector_event.next=0;
   Vy_sum_corrector_event.period = config.dt * 50;
   Vy_sum_corrector_event.config_dt = config.dt;

//  config.grid.init(config.cell.L,UNIT_LENGHT*1.6);

  config.check_neighbours.next=tstart;

    if(config.parameter.shear_rate==0)config.check_neighbours.period=0.01;
    else  config.check_neighbours.period= 0.02/config.parameter.shear_rate;
    config.check_neighbours.config_dt=0;

    // LINEAR VELOCITY
    config.set_linear_velocity_profile();


//config.check_neighbours.period = 5*config.dt;
  Ctimer run_time;
  run_time.start();

  // FINER SAVE TIME array
  double fine_dts[40] = {0,0.001,0.002,0.005,
                       0.01,0.015,0.02,0.03,0.05,0.07,
                       0.1,0.15,0.2,0.3,0.5,0.7,
                       1,1.5,2,3,5,7,
                       10,15,20,30,50,70,
                       100,150,200,300,500,700,
                       1000,1500,2000,3000,5000,7000};
  double prev_saved_t = 0;
  double next_save_t_force = 0;
  int fine_count = 0;

  // Reset XTrue values
  for(int i=0;i<config.grain.size();i++){
      config.grain[i].Xtrue*=0;
  }

  cout<<"before loop"<<endl;
  for(config.t=0;config.t<tend;config.t+=config.dt)
    {
      config.evolve();

      if(save_event.should_do(config.t))
      {
        config.cell.evale_volume_grain(config.grain);

        config.save(out_files);
        out_files.iSave++;

        cout<<out_files.path_root<<" ";
        cout << "dt ="<<config.dt<<" time: "<<config.t<<" finishes at " <<tend<<endl;
        cout <<"Cell dimension "; config.cell.L.PRINT();
        cout <<"Vcell ";  config.cell.V.PRINT();

        cout <<"internal stress_yy,xy: "<< config.cell.stress.x[1][1]<<", "<<config.cell.stress.x[0][1]<<endl;
        cout <<"boundary stress_yy,xy: "<< config.cell.boundary_normal_stress<<", "<<config.cell.boundary_shear_stress<<endl;

        run_time.end();
        cout<< "run time (s)"<<0.001*run_time.duration()<<endl;
        run_time.start();
        cout<< "need reset?: "<< need_reset <<endl;

        if (need_reset & (config.cell.V.x[1]==0 | config.cell.L.x[1] == config.cell.fixh)) {
            config.set_linear_velocity_profile();
            need_reset = false;
        }

        // reset fine count
        prev_saved_t = config.t;
        fine_count = 0;

      }

      if((ndt ==1 | ndt>=3) & config.t > (prev_saved_t+fine_dts[fine_count]-config.dt/2)){
          config.save_dt_grains(out_files);
          fine_count+=1;
      }

      if (ndt >= 2 & config.t > (next_save_t_force-config.dt/2)) {
          foreach(Cgrain g, config.grain) g.evale_stress();
          config.save_dt_forces(out_files);
          next_save_t_force += 0.01;
      }

//      if(Vy_sum_corrector_event.should_do(config.t)){

//          config.Vy_sum_corrector();

//      }

    }
  cout<<"after loop"<<endl;

    run_time.end();
    cout <<"time step"<<config.dt<<endl;
    cout<< "run time (s)"<<0.001*run_time.duration()<<endl;
    cout<< "time for 100 000 (h) "<<0.001*run_time.duration()*100000/3600.<<endl;

  config.parameter.PRINT();

}




// SIMULE_CON
void Crun::simule_con()
{
  in_files.set_path("READ");
  out_files.set_path("READ_WRITE");

  config.read(in_files);
  cout<<config.grain.size()<<endl;
  cout<<"contact size: "<<config.grain[0].contact.size()<<endl;



//  config.save(out_files);

 // cout<<"Enter the begining time (key TSTART)"<<endl;  in_files.get_secure("TSTART",tstart);
  //cout<<"Enter the end time (key TEND)"<<endl;  in_files.get_secure("TEND",tend);


//  cout<<"Enter the begining time of saving (key SAVE_START)"<<endl;in_files.get_secure("SAVE_START",save_event.next);
 // cout<<"Enter the period of time between saving (key SAVE_PERIOD)"<<endl;in_files.get_secure("SAVE_PERIOD",save_event.period);


  cout<<"Enter the shear rate (key SHEAR_RATE)"<<endl;in_files.get_secure("SHEAR_RATE",config.parameter.shear_rate);
  cout<<"Enter the cohesion (key COHESION)"<<endl;  in_files.get_secure("COHESION",config.parameter.Cohesion);

  double swips;
  cout<<"Enter no. of swips (key SWIPS)"<<endl;  in_files.get_secure("SWIPS",swips);
  double fps;
  cout<<"Enter no. of frames per swip (key FPS)"<<endl; in_files.get_secure("FPS",fps);

  int ndt;
  cout<<"Enter no. n for finer grain frames at n*dt interval (or -1 otherwise)] (key NDT)"<<endl; in_files.get_secure("NDT",ndt);

  double fixh;
  cout<<"Fixed height? [if no enter -1 (or any negative value), if yes enter the height] (key FIXH)"<<endl; in_files.get_secure("FIXH",fixh);


  double mu_p;
  cout<<"Friction? (key FIXH)"<<endl; in_files.get_secure("MU",mu_p);

  double K_value;
  cout<<"K value? (for default input -1) (key KVALUE)"<<endl; in_files.get_secure("KVALUE",K_value);

  config.parameter.mu = mu_p;

  if (K_value > 0){
      config.parameter.E = K_value;
      config.parameter.Stiff = K_value*UNIT_LENGHT;
  }

  config.cell.fixh = fixh;
//  if (fixh > 0) {
//      config.cell.L.x[1] = fixh;
//      config.cell.V.x[1] = 0;
//  }

//  double gettstart;
//  cout<<"Enter tstart, (key GETTSTART)"<<endl; in_files.get_secure("GETTSTART",gettstart);



   // **************************************************************


//   config.dt = sqrt(UNIT_MASS/(config.parameter.E*UNIT_LENGHT)) /20; //fraction of the collision time

//   config.dt = 1./1000;

  //Round to the nearest power of ten
   config.dt = 1./pow(10,ceil(abs(log10(sqrt(UNIT_MASS/(config.parameter.E*UNIT_LENGHT))/20))));


   Cevent save_event;
   tstart=config.t;
   cout<<"tStart: "<<tstart<<endl;
   tend=tstart+swips/config.parameter.shear_rate;


   save_event.period=1/config.parameter.shear_rate/fps;
   save_event.next=tstart + save_event.period;
    save_event.config_dt = config.dt;

   Cevent save_dt_grain_event;
   save_dt_grain_event.next=0;
   save_dt_grain_event.period = config.dt * ndt;
   save_dt_grain_event.config_dt = config.dt;

   // create event to track when to correct Vy sum
   Cevent Vy_sum_corrector_event;
   Vy_sum_corrector_event.next=0;
   Vy_sum_corrector_event.period = config.dt * 50;
   Vy_sum_corrector_event.config_dt = config.dt;

//  config.grid.init(config.cell.L,UNIT_LENGHT*1.6);

  config.check_neighbours.next=tstart;

    if(config.parameter.shear_rate==0)config.check_neighbours.period=0.2;
    else  config.check_neighbours.period= 0.02/config.parameter.shear_rate;

    // LINEAR VELOCITY
//       config.set_linear_velocity_profile();


//config.check_neighbours.period = 5*config.dt;
  Ctimer run_time;
  run_time.start();


  for(config.t=tstart+config.dt;config.t<tend;config.t+=config.dt)
    {
      config.evolve();

      if(save_event.should_do(config.t))
      {
        config.save(out_files);
        out_files.iSave++;

        cout<<out_files.path_root<<" ";
        cout << "dt ="<<config.dt<<" time: "<<config.t<<" finishes at " <<tend<<endl;
        cout <<"Cell dimension "; config.cell.L.PRINT();
        cout <<"Vcell ";  config.cell.V.PRINT();

        cout <<"internal stress_yy: "<< config.cell.stress.x[1][1]<<endl;
     //   cout <<"boundary stress_yy: "<< config.cell.stress_boundary.x[1][1]<<endl;
        run_time.end();
        cout<< "run time (s)"<<0.001*run_time.duration()<<endl;
        run_time.start();


        }

      if(save_dt_grain_event.should_do(config.t)){
          config.save_dt_grains(out_files);
      }

//      if(Vy_sum_corrector_event.should_do(config.t)){

//          config.Vy_sum_corrector();

//      }


    }

//    run_time.end();
  //  cout <<"time step"<<config.dt<<endl;
   // cout<< "run time (s)"<<0.001*run_time.duration()<<endl;
   // cout<< "time for 100 000 (h) "<<0.001*run_time.duration()*100000/3600.<<endl;

  config.parameter.PRINT();

}


double Crun::uniform_random(double min, double max)
{
     return min+(max-min)*(rand()/((double)RAND_MAX+1));
}

void Crun::create_random_config()
{
  cout<<"Start creating a random configuration"<<endl;
  srand( (unsigned int) time(NULL)+getpid());//Init of random function


  out_files.set_path("WRITE");

  cout<<"Enter the cell dimension (x,y,z) (key CELL)"<<endl; in_files.get_secure("CELL",config.cell.L);

  int Ngrain;
  cout<<"Enter the number of grains (key GRAIN)"<<endl; in_files.get_secure("GRAIN",Ngrain);

  double Rmin,Rmax;
  cout<<"Enter Rmin (key RMIN)"<<endl;  in_files.get_secure("RMIN",Rmin);
  cout<<"Enter Rmax (key RMAX)"<<endl;  in_files.get_secure("RMAX",Rmax);


  config.cell.mass= config.cell.L.x[0]*4/PI*UNIT_LENGHT*UNIT_LENGHT;


  //======== Set the flowing grains
 cout<<" - set the flowing grains"<<endl;

  QList <double> radius;
  for(int i=0;i<Ngrain;i++) radius.push_back(uniform_random(Rmin,Rmax));
  qSort(radius.begin(), radius.end(), qGreater<double>());// sort the list from larger particle to smaller ones, so that it's would be easier to set them randomly (set bigger first)


  for(int i=0;i<Ngrain;i++)//flowing grains
    {
      Cgrain g;
      g.R=radius[i];
      g.evale_mass();
      g.ID=config.grain.size();

      bool overlap=false;
      do
        {
          g.X.x[0]=uniform_random(0,config.cell.L.x[0]);
          g.X.x[1]=uniform_random(-config.cell.L.x[1]/2.,config.cell.L.x[1]/2.);

          overlap=false;
          for(int i=0;i<config.grain.size();i++)
            {
              Ccontact C(&g,&config.grain[i],&config.cell,&config.parameter);
              C.get_distance();
              if((C.delta < -0.00*UNIT_LENGHT)){overlap=true;break;} //at the first contact, stop checking and rember that there is an overlap
            }
        } while(overlap);
      config.grain.push_back(g);
      if(config.grain.size()%20==0)cout<<config.grain.size() << " grains set in cell"<<endl;
    }


  out_files.set_file_name();
  config.save(out_files);

  cout<<"Random configuration created"<<endl;
}





