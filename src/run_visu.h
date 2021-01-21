#ifndef RUN_VISU_H
#define RUN_VISU_H

#include "main.h"
#include "run.h"

#include <GLUT/glut.h>


#if defined(__APPLE__) || defined(MACOSX)
#include <GLUT/glut.h>
#else
#include <gl/glut.h>
#endif





class Crun_visu : public Crun
{
public:
  Crun_visu();
 void visu( int argc, char *argv[ ]);

 static void Keyboard_input(unsigned char key, int x, int y);
 static void Reshape(int width, int height);

 static void show_cell();
 static void draw_grain(Cgrain &);
 static void draw_contact(Ccontact &);
 static void draw_velocity(Cgrain &g);
 static void save_image();


 static void Draw();
 static void set_camera();
 static void set_light();

 static void InitGL();

  double Rad(double deg){return deg*PI/180.;};
};

#endif // RUN_VISU_H
