#ifndef IN_OUT_H
#define IN_OUT_H

#include "main.h"
#include "vector.h"
#include <sys/stat.h> /**<used for mkdir();*/
#include <dirent.h>   /**<used for opendir();*/


class Cin_out
{
  public:
  int iSave;

  string path_root;
  QList <string> path;  //contain in this order: path root + /grain/ /contact/ /cell/ /parameter/ /profile/ /mesh/ /grain_post_process/

  string file_grain;
  string file_grain_pp;
  string file_contact;
  string file_cell;
  string file_parameter;
  string file_profile;
  string file_mesh;

  Cin_out();

  void set_path( string);
  void set_file_name();
  int get_Nline(string);    /*<< returns the number of line in the file*/

  void error(string, string ,string );   /*<<  writes an error message and quit the program*/
  bool exist_files();                     /*<<Return false if one on the file can not be opened, true if they all can be open*/

  template<typename T> string to_string( const T & Value ); //convert anything into string

  void get_secure(string lock, double & value);
  void get_secure(string lock, int & value);
  void get_secure(string lock, Cvector & v);
  void get_secure(string lock, string & s);

};

#endif // IN_OUT_H
