#ifndef RUN_POST_PROCESS_H
#define RUN_POST_PROCESS_H
#include "main.h"
#include "run.h"
#include "profile.h"
#include "mesh.h"

class Crun_post_process : public Crun
{
public:
  Crun_post_process();


  void evale_profile();

  void evale_voro();

  void evale_corr(); // Evaluate auto correlation

  void evale_length(); // Evaluate the typical length

};

#endif // RUN_POST_PROCESS_H
