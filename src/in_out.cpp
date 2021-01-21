#include "in_out.h"


template<typename T>
string Cin_out::to_string( const T & Value )//convert anything into string
{
  ostringstream oss;
  oss << Value;
  return oss.str();
}

void Cin_out::get_secure(string lock,  double & value )
{
  string key;
  cin>>key;
  if(key!=lock){string err="The entry'"+key+"' is not valid"+"\n"+"The entry'"+lock+"' was expected"; error("in_out.cpp","Cinout()",err);}
  cin>>value;
}

void Cin_out::get_secure(string lock,  int & value )
{
  string key;
  cin>>key;
  if(key!=lock){string err="The entry'"+key+"' is not valid"+"\n"+"The entry'"+lock+"' was expected"; error("in_out.cpp","Cinout()",err);}
  cin>>value;
}


void Cin_out::get_secure(string lock,  Cvector & v )
{
  string key;
  cin>>key;
  if(key!=lock){string err="The entry'"+key+"' is not valid"+"\n"+"The entry'"+lock+"' was expected"; error("in_out.cpp","Cinout()",err);}
  v.READ();
}

void Cin_out::get_secure(string lock,  string & v )
{
  string key;
  cin>>key;
  if(key!=lock){string err="The entry'"+key+"' is not valid"+"\n"+"The entry'"+lock+"' was expected"; error("in_out.cpp","Cinout()",err);}
  cin>>v;
}


void Cin_out::error(string file,string fct,string message){
  cout<<endl<<endl;
  cout<<"\tError in file '"<< file<<"' function '" << fct<<"':"<<endl;
  cout<<"====>>"<<message<<endl;
  cout<<endl<<endl;
  exit(0);
};

Cin_out::Cin_out(){};

bool Cin_out::exist_files()
{
  set_file_name();

  ifstream fp;

  fp.open(file_grain.c_str());
  if(!fp.is_open()) return false;  // can not open the file, return 0
  fp.close();
  fp.open(file_contact.c_str());
  if(!fp.is_open()) return false;  // can not open the file, return 0
  fp.close();
  fp.open(file_cell.c_str());
  if(!fp.is_open()) return false;  // can not open the file, return 0
  fp.close();
  fp.open(file_parameter.c_str());
  if(!fp.is_open()) return false;  // can not open the file, return 0
  fp.close();return true;
}

void Cin_out::set_path(string mode)
{
 if(mode!="WRITE"&&mode!="READ"&&mode!="READ_WRITE"){string err="Path opening mode should be either WRITE, READ or READ_WRITE, the mode '"+mode+"' is not valid"; error("in_out.cpp","Cinout()",err);}


  string key = "PATH_"+mode;
  cout<<"Enter the path where to "<<mode<<" (key "<< key<<")"<<endl;
  get_secure(key,path_root);

  cout<<"Enter the number of the first file (key 'FIRST_FILE')"<<endl;
  get_secure("FIRST_FILE",iSave);

  path.clear();                       //clear previous list of path
  path.push_back(path_root+"/grain");
  path.push_back(path_root+"/contact");
  path.push_back(path_root+"/para");
  path.push_back(path_root+"/cell");
  path.push_back(path_root+"/profile");
  path.push_back(path_root+"/mesh");
  path.push_back(path_root+"/grain_post_process");

  if(mode=="READ_WRITE"){mkdir(path.last().c_str(),0755); return;}

  bool exist = opendir(path_root.c_str());

  if(exist&&mode=="WRITE"){string err="Can not write on: '"+path_root+"' the directory already exists."; error("in_out.cpp","Cinout()",err);}
  if(!exist&&mode=="WRITE")
    {
      if(mkdir(path_root.c_str(),0755)==1) {string err="Can not creat the path: '"+path_root+"'"; error("in_out.cpp","Cinout()",err);}  //create the root path

      for(int p=0;p<path.size();p++)
        if(mkdir(path[p].c_str(),0755)==1) {string err="Can not creat the path: '"+path[p]+"'"; error("in_out.cpp","Cinout()",err);}//create the sub path
    }

  if(!exist&&mode=="READ"){string err="Can not read from: '"+path_root+"' the directory do not exists."; error("in_out.cpp","Cinout()",err);}



}

void Cin_out::set_file_name()
{
  file_grain=path[0]+to_string("/grain_")+to_string(iSave);
  file_contact=path[1]+to_string("/contact_")+to_string(iSave);
  file_parameter=path[2]+to_string("/parameter_")+to_string(iSave);
  file_cell=path[3]+to_string("/cell_")+to_string(iSave);
  file_profile=path[4]+to_string("/profile_")+to_string(iSave);
  file_mesh=path[5]+to_string("/mesh_")+to_string(iSave);
  file_grain_pp=path[6]+to_string("/grain_pp")+to_string(iSave);

}

int Cin_out::get_Nline(string file)
{
  int Nline=0;

  ifstream fp;

  fp.open(file.c_str());
  if(!fp.is_open()) return Nline;  // can not open the file, return 0 line

  while(1<2)
    {
      string line;
      getline(fp,line);
      if(!fp.eof()) Nline++;
      else break;
    }
  fp.close();


  return Nline;
}
