#include "run_visu.h"

//#include "/usr/local/include/tiffio.h"


double R_cam;//disance from the center
double theta_cam;//angle from x axis
double phi_cam;//angle from xz plane
double step; //angle step for rotation
double x_look,y_look,z_look;//position where the camera look at
double l_step;//distance of translation in x y z direction

double zoom;

double magnify ;//magnification of velocities arrow
bool show_grain=true;
bool show_contact=true;
bool show_velocity=false;

QList <int >ID_to_show;

Crun_visu *run_v;

Crun_visu::Crun_visu()
{

  run_v=this;
  run_v->in_files.set_path("READ");

  run_v->config.read(in_files);
  //run_v->config.read_grain_exp(in_files,40.,1.);

   R_cam=200.;//disance from the center
   theta_cam=Rad(90);//angle from x axis
   phi_cam=Rad(0);//angle from xz plane
   step = Rad(5); //angle step for rotation
    zoom=0.5;

   double Vmax=0;
   foreach(Cgrain g,run_v->config.grain)
   {if(g.V.norm()>Vmax)Vmax=g.V.norm();}
   magnify=2./Vmax;


   x_look= 0;
   y_look=z_look=0;//position where the camera look at
   l_step=0.2;//distance of translation in x y z direction

   foreach(Cgrain g,run_v->config.grain)
    {
       if(fabs(g.X[1])<5)
       ID_to_show.push_back(g.ID);
   }
}

void Crun_visu::save_image()
{
//    TIFF *file;
//    GLubyte *image, *p;
//    int i;
//    char filename[255];
//    int x, y, width, height,  compression;

//    cout <<"test"<<endl;
//    sprintf(filename,"/Users/mattmacaulay/Desktop/test.tif");
//    x=0;
//    y=0;
//    //compression=COMPRESSION_PACKBITS;
//    compression= COMPRESSION_NONE;
//    width= 1000;//WIDTH; //taille de la fenètre
//    height=800;//HEIGHT;

//    file = TIFFOpen(filename, "w");
//    if (file == NULL)
//    {
//        cout<<"Can't open "<<endl;
//        return;
//    }
//    image = (GLubyte *) malloc(width * height * sizeof(GLubyte) * 3);

//    /* OpenGL's default 4 byte pack alignment would leave extra bytes at the
//       end of each image row so that each full row contained a number of bytes
//       divisible by 4.  Ie, an RGB row with 3 pixels and 8-bit componets would
//       be laid out like "RGBRGBRGBxxx" where the last three "xxx" bytes exist
//       just to pad the row out to 12 bytes (12 is divisible by 4). To make sure
//       the rows are packed as tight as possible (no row padding), set the pack
//       alignment to 1. */

//    glPixelStorei(GL_PACK_ALIGNMENT, 1);

//    glReadPixels(x, y, width, height, GL_RGB, GL_UNSIGNED_BYTE, image);
//    TIFFSetField(file, TIFFTAG_IMAGEWIDTH, (uint32) width);
//    TIFFSetField(file, TIFFTAG_IMAGELENGTH, (uint32) height);
//    TIFFSetField(file, TIFFTAG_BITSPERSAMPLE, 8);
//    TIFFSetField(file, TIFFTAG_COMPRESSION, compression);
//    TIFFSetField(file, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_RGB);
//    TIFFSetField(file, TIFFTAG_SAMPLESPERPIXEL, 3);
//    TIFFSetField(file, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
//    TIFFSetField(file, TIFFTAG_ROWSPERSTRIP, 1);
//    TIFFSetField(file, TIFFTAG_IMAGEDESCRIPTION, "OpenGL-rendered sauv");
//    p = image;
//    for (i = height - 1; i >= 0; i--)
//    {
//        if (TIFFWriteScanline(file, p, i, 0) < 0)
//        {
//            free(image);
//            TIFFClose(file);
//            return;
//        }
//        p += width * sizeof(GLubyte) * 3;
//    }
//    TIFFClose(file);




}




void Crun_visu::set_light()
{
  //first light
  //int LightPos0[4] = {3000,3000,3000,1};
 // glLightiv(GL_LIGHT0,GL_DIFFUSE,LightPos0);

  //float LightDif0[4] = {1.f,1.f,-1.f,.5f};
  //glLightfv(GL_LIGHT0,GL_DIFFUSE ,LightDif0);

 /* //second light
  int LightPos1[4] = {-3000,0,3000,1};
  glLightiv(GL_LIGHT1,GL_POSITION,LightPos1);
  float LightDif1[4] = {1.f,.5f,.5f,1.f};
  glLightfv(GL_LIGHT1,GL_DIFFUSE,LightDif1);


  //second light
  int LightPos2[4] = {0,-3000,+3000,1};
  glLightiv(GL_LIGHT1,GL_POSITION,LightPos2);
  float LightDif2[4] = {.5f,1.f,.5f,1.f};
  glLightfv(GL_LIGHT2,GL_DIFFUSE,LightDif2);
*/
}

void Crun_visu::set_camera()
{

  float x,y,z;//position of the camera

  x= R_cam*cos(theta_cam)*cos(phi_cam);
  y= R_cam*sin(phi_cam);
  z= R_cam*sin(theta_cam)*cos(phi_cam);

  //float nx=0,ny=1,nz=0;
  float nx,ny,nz;

  nx= cos(theta_cam+PI)*cos(phi_cam+PI/2);
  ny= sin(phi_cam+PI/2);
  nz= sin(theta_cam+PI)*cos(phi_cam+PI/2);

  gluLookAt(x,y,z,x_look,y_look,z_look,nx,ny,nz);//look at 000
}



void Crun_visu::show_cell()
{

}

void Crun_visu::draw_grain(Cgrain &g)
{
  if(!show_grain) return;

  glPushMatrix();
  glTranslatef(g.X.x[0],g.X.x[1],g.X.x[2]);
  glScalef(g.R,g.R,g.R);


  double dV =0;
  dV = g.V[0] -g.X[1] * g.cell->V[0]/g.cell->L[1];



//    if(dV<0)
     // GLfloat color[] = {0.44f, 0.16f, 0.78f, 1.0f};
    //  glMaterialfv(GL_FRONT_AND_BACK, GL_EMISSION/*GL_AMBIENT_AND_DIFFUSE*/, color);

//    g.evale_shear_rate(run_v->config.grain);

    double c= (g.gradV.x[0][1]+g.gradV.x[1][0])/run_v->config.parameter.shear_rate - 1;
    c*=2;

    glColor4f(0.369*c+0.55,0.314,-0.378*c+0.484,1.   );
   // else   glColor4f(0,/*111./255.*/,0.,1.0);



  //    glPushMatrix();
  glutSolidSphere(1,50,50);
//  glPopMatrix();
  glPopMatrix();
}

void Crun_visu::draw_contact(Ccontact &c)
{
  if(!show_contact) return;

  double magn = c.Force.norm()*2;
      //log(10*c.Force.norm());
  if (magn<1) magn=1;
  glLineWidth(magn);

//  glColor4f(1.0, 0.0, 0.0,0.5);
//glMateriali(GL_FRONT_AND_BACK,GL_SHININESS,50);

//int MatSpec [4] = {1,0,0,1};

//glMaterialiv(GL_FRONT_AND_BACK,GL_EMISSION,MatSpec);

  glBegin(GL_LINES);

  glColor4f(/*1.0*/magn/10., 0.0, 0.0,0.5);
  glVertex3f(c.A->X[0], c.A->X[1], c.A->X[2]);

  Cvector dX = c.B->X-c.A->X;
  Cvector dv;
  run_v->config.cell.CLP(dX,dv);
  dX+=c.A->X;
  glColor4f(1.0, 0.0, 0.0,0.5);
  glVertex3f(dX[0], dX[1], dX[2]);
  glEnd();


}

void Crun_visu::draw_velocity(Cgrain &g)
{
  if(!show_velocity)return;

  Cvector dV = g.V;
    g.V.PRINT();
  dV[0]-=g.X[1]*run_v->config.parameter.shear_rate;
  dV*=magnify;

  glLineWidth(1.5);

  glBegin(GL_LINES);
  glColor3f(0.0, 0.0, 1.0);
  glVertex3f(g.X[0], g.X[1], g.X[2]);
  glColor3f(1.0, 0.0, 0.0);

  glVertex3f(g.X[0]+dV[0], g.X[1]+dV[1], g.X[2]+dV[2]);
  glEnd();

}


void Crun_visu::Draw()
{
  glutSetWindowTitle(run_v->in_files.file_grain.c_str());

  //glMateriali(GL_FRONT_AND_BACK,GL_SHININESS,50);

  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glMatrixMode(GL_MODELVIEW);

  glLoadIdentity();

  set_camera();
  set_light();

  //show_cell();
  glScalef(zoom, zoom, 1.0f);

  glTranslatef(-run_v->config.cell.L[0]/2,0,0 );


  foreach(Cgrain g, run_v->config.grain)

      //for(int i=0;i<ID_to_show.size();i++)
        //if(g.ID==ID_to_show[i])
        {
       draw_grain(g);
      if(show_contact)      foreach(Ccontact c, g.contact)draw_contact(c);
      if(show_velocity)   draw_velocity(g);
    //        break;
        }

//  if(c>=0)   glColor4f(0,+111./255.,c/2,1.0);
 // else   glColor4f(0.,-c/2.+111./255.,0.,1.0);

  glLineWidth(3);
  glBegin(GL_LINE_LOOP);
 glColor3f(0, 0, 0);
  glVertex3f( 0, -run_v->config.cell.L[1]/2,1 );
    glVertex3f( 0, +run_v->config.cell.L[1]/2,1 );
    glVertex3f(run_v->config.cell.L[0], +run_v->config.cell.L[1]/2,1  );
    glVertex3f(run_v->config.cell.L[0], -run_v->config.cell.L[1]/2,1);
    glEnd();

//color scale
/*  glLineWidth(15);
  glBegin(GL_LINES);
  glColor3f(46./255, 80./255, 220./250);
  glVertex3f(run_v->config.cell.L[0]+2, - run_v->config.cell.L[1]/2, 1);
  glColor3f(234./255, 80./255, 27./250);
  glVertex3f(run_v->config.cell.L[0]+2, run_v->config.cell.L[1]/2, 1);
  glEnd();
*/


  glutSwapBuffers();// pour glut
  }



void Crun_visu::InitGL(){
  glEnable(GL_DEPTH_TEST);    //make sure faces that are behind other faces are not shown
 // glEnable(GL_LIGHTING);	// Switch on lightening
 // glEnable(GL_LIGHT0);          // Switch on first light
  //glEnable(GL_LIGHT1);          // Switch on second light
  //glEnable(GL_LIGHT2);          // Switch on second light

  glEnable(GL_COLOR_MATERIAL);   //enable the coloring of material

//  int MatSpec [4] = {1,0,0,1};

   // glMaterialiv(GL_FRONT_AND_BACK,GL_,MatSpec);
   glMateriali(GL_FRONT_AND_BACK,GL_SHININESS,50);
  glClearColor(1,1,1,1);    //background color: white
};


void Crun_visu::Reshape(int width, int height)
{
  glViewport(0,0,width,height);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluPerspective(20,float(width)/float(height),0.1,400);	//Pour les explications, lire le tutorial sur OGL et win
  glMatrixMode(GL_MODELVIEW);	//Optionnel
  glutPostRedisplay();

}


 void Crun_visu::Keyboard_input(unsigned char key, int x, int y)
{
   x=y;
  switch (key)
    {
    case 'r' :  theta_cam-=step;  break;
    case 'R' :  theta_cam+=step;  break;
    case 's' :  phi_cam-=step;  break;
    case 'S' :  phi_cam+=step;  break;
    case  'x':  x_look-=l_step;  break;
    case  'X':  x_look+=l_step;  break;
    case  'y':  y_look-=l_step;  break;
    case  'Y':  y_look+=l_step;  break;
    case  'z':  z_look-=l_step;  break;
    case  'Z':  z_look+=l_step;  break;

    case '-'  : zoom/=1.3;break;//R_cam*=1.3; break;
    case '+'  : zoom*=1.3;break; //R_cam/=1.3;break;


    case '4': save_image();break;

    case 'n'  : run_v->in_files.iSave++;
      if(run_v->in_files.exist_files())
      {
          run_v->config.read(run_v->in_files);
          Cconfig conf_next;
           run_v->in_files.iSave++;
           run_v->in_files.set_file_name();
          conf_next.read( run_v->in_files);
           run_v->in_files.iSave--;
           run_v->in_files.set_file_name();

      }
      else  run_v->in_files.iSave--; break;

    case 'p'  : run_v->in_files.iSave--;
      if(run_v->in_files.exist_files())
      {
          run_v->config.read(run_v->in_files);
          Cconfig conf_next;
           run_v->in_files.iSave++;
           run_v->in_files.set_file_name();
          conf_next.read( run_v->in_files);
           run_v->in_files.iSave--;
           run_v->in_files.set_file_name();
      }
      else  run_v->in_files.iSave++; break;


    case 'g'  : if(show_grain) show_grain=false; else show_grain=true ;break;
    case 'c'  : if(show_contact) show_contact=false; else show_contact=true; break;
    case 'v'  : if(show_velocity) show_velocity=false;else show_velocity=true;break;

    case 'M'  : magnify*=1.5;break;
    case 'm'  : magnify/=1.5;break;

    }
  glutPostRedisplay();

 }


void Crun_visu::visu(int argc, char *argv[])
{
  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH);
  glutInitWindowSize(1000,800);
  glutInitWindowPosition(100,100);

  glutCreateWindow("");
  InitGL();

  glutReshapeFunc(Reshape);
  glutDisplayFunc(Draw);
  glutKeyboardFunc(Keyboard_input);
  glutMainLoop();


}
