#ifndef GRAIN_H
#define GRAIN_H

#include "vector.h"
#include "parameter.h"
#include "cell.h"
#include "box.h"
#include "contact.h"

class Cbox;
class Ccontact;
class Ccell;

class Cgrain
{
public:
  Cvector X,V,A,Xtrue;
  // Note contents of Xtrue:
  // Xtrue.x[0] = no. of level of cell (x direction), i.e. dx [+1 if went right, -1 if left]
  // Xtrue.x[1] = no. of level of cell (y direction), i.e. dy [+1 if went up, -1 if down]
  // Xtrue.x[2] = nothing

  double Ome,Ome_dot,theta;
  Cvector Fsum;                        //sum of forces
  double Msum;                         //sum of moments
  double R;                            //radius
  double mass,Imass;                   //mass - moment inertia
  int ID;                              //grain number

  // added for voro
  double voro_area; // voronoi area for this grain
  double voro_area_n; // voronoi area of collection of grains
  double solid_area_n; // solid area of collection of grains
  QList <Cgrain *> voro_neighbour; // list of neighbouring grains as per voro

  Cmatrix gradV;
  double shear_rate,vorticity;         //averaged shear rate dVx/dy around the grain
  Cmatrix stress;

  Cparameter *para;
  Ccell *cell;
  Cbox *box;                        //pointer to the cell the grain is in
  QList <Cgrain *> neighbour;       //list of neighbouring grains
  QList <Ccontact> contact;         //list of contacting grains

  Cgrain(Ccell *c,Cparameter*p){para=p;cell=c;Ome=theta=0;voro_area=0;}
 // Cgrain(Cparameter *p){para=p;};

  Cgrain(){};

  void evale_mass(void);
  void evale_shear_rate(QList <Cgrain> grain);
  void evale_stress(void);

  void evale_from_voro(int lvl);
  Cvector interp_vel(Cvector pos, QList<Cgrain> interp_grains);

  void predictor(const double);
  void corrector(const double);

  //void get_neighbour();
  void update_contact(const double dt);          //update the list of contact from the list of neighbours


  void convert_unit(double pix_per_d,double fps);


  void PRINT(FILE *fp);
  void READ(FILE *fp);

  void PRINT_pp(FILE *fp);
  void READ_pp(FILE *fp);

};

#endif // GRAIN_H
