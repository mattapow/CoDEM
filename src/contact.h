#ifndef CONTACT_H
#define CONTACT_H

#include "parameter.h"
#include "cell.h"
#include "grain.h"

class Cgrain;
class Ccell;


class Ccontact
{
public:

  Cparameter *para;
  Ccell *cell;

  Ccontact(){init();};                                /*<< Default constructor, set data to zero but there is no pointer to grains*/
//  Ccontact(Cgrain *a,Cgrain *b){A=a;B=b; init();};    /*<< Constructor setting data to zero and recording the pointers to grains*/
  Ccontact(Cgrain *a,Cgrain *b,Ccell *c,Cparameter *p){A=a;B=b; init(); cell=c;para=p;};    /*<< Constructor setting data to zero, recording the pointers to grains, and recording the cell pointer*/


  void init(){Force*=0;nAB*=0;delta=0;age=0;fn_el=0;fn_vis=0;};        /*<<set all the data to zero*/

  Cgrain *A,*B;                                       /*<<Pointers to the two grains involved*/

  Cvector Force;  /*<< Total contact force */

  Cvector Vclp;  /*<<Contain the velocity related to periodic condiction, may be null, or +-Vcell */
  double fn_el,fn_vis, ft;   /*<<Normal elastic and viscous forces, tangential force */
//    Cvector Ft;
  Cvector nAB,nABt;    /*<< Unit vector from A to B centers */
  double delta;   /*<< Interpenetration, negative if in contact*/
  double age;     /*<< Age of the contact */

  int aID,bID;  /*<< Id of involved grains, only needed for read function */
  void get_distance();
  bool not_in_contact();
  void evale_force(double dt);
  void PRINT();
  bool AM_I_CROSS_BOUNDARY();
  void PRINT(FILE *fp);
  void READ(FILE *fp);

};

#endif // CONTACT_H
