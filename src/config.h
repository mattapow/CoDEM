#ifndef CONFIG_H
#define CONFIG_H

#include "grain.h"
#include "cell.h"
#include "parameter.h"
#include "grid.h"
#include "in_out.h"
#include "event.h"




class Cconfig
{
public:
  Cconfig();

  QList < Cgrain > grain;
  QList <Cgrain *> top_wall;    /**<List of pointer to the grains of the top wall - implemented in the config.read function */
  QList <Cgrain *> bot_wall;    /**<List of pointer to the grains of the top wall - implemented in the config.read function */

  //  QList < Ccontact > contact;
  Ccell cell;
  Cparameter parameter;

  Cgrid grid;
  double dt,t;


  Cevent check_neighbours;

  void evolve();
  void predictor( );
  void refresh_contact(bool force_check_neighbours);
  void sum_force();
  void corrector( );
  void Vy_sum_corrector();

  void evale_stress();

  void save(Cin_out );
  void save_dt_grains(Cin_out out_file);
  void save_dt_forces(Cin_out out_file);
  void read(Cin_out );

  void save_grain_postprocess(Cin_out);
  void read_grain_postprocess(Cin_out);

  void set_linear_velocity_profile();

};


class Ctimer
{
public:
  clock_t tstart,tend;
  void start(){tstart=clock(); };
  void end(){tend=clock(); };
  double duration(){ return (1000*((double)tend-(double)tstart)/CLOCKS_PER_SEC);}/**< Return time in ms*/
};


#endif // CONFIG_H
