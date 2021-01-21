#ifndef RUN_H
#define RUN_H

#include "main.h"
#include "in_out.h"
#include "event.h"
#include "config.h"

class Crun
{
public:
  Crun();
  Cconfig config;
  Cin_out in_files,out_files;
  double tstart,tend;

  void simule();
  void simule_con();
  void create_random_config();

  double uniform_random(double min, double max);//return a random number between min and max

};




#endif // RUN_H
