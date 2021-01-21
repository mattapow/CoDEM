
#include "main.h"

#include "in_out.h"
#include "run.h"
// #include "run_visu.h"
#include "run_post_process.h"



int main(int argc, char *argv[])
{
  cout << "=============== DEM 2d ================" << endl;
  cout << "========== P. Rognon 03/2013 ========== " << endl;
  cout << "= Modified by P. Kharel & M. Macaulay =" << endl << endl;

  Cin_out in_out;

  string  action;
  // cout<<"Enter the action to perform (key: 'ACTION'', choices: 'CREATE_RANDOM' 'SIMULE' 'SIMULE_CON' 'VISU' 'POST_PROCESS' 'POST_PROCESS_EXP' 'VORO_PP' 'CORR', 'LENGTH') "<<endl;
  cout<<"Enter the action to perform (key: 'ACTION'', choices: 'CREATE_RANDOM' 'SIMULE' 'POST_PROCESS') "<<endl;
  in_out.get_secure("ACTION",action);

    if(action == "CREATE_RANDOM")    {Crun run; run.create_random_config();}
    if(action == "SIMULE") {Crun run; run.simule();}
    // if(action == "SIMULE_CON") {Crun run; run.simule_con();}
    if(action == "POST_PROCESS") {Crun_post_process run_pp; run_pp.evale_profile();}
    // if(action == "VORO_PP"){Crun_post_process run_pp; run_pp.evale_voro();}
    // if(action == "CORR"){Crun_post_process run_pp; run_pp.evale_corr();}
    // if(action == "LENGTH"){Crun_post_process run_pp; run_pp.evale_length();}
    // if(action == "VISU")    { Crun_visu run_visu; run_visu.visu(argc,argv);}

    return 0;
}

