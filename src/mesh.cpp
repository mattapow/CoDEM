#include "mesh.h"

Cmesh::Cmesh(Cvector L, double step)
{
  // initialise the 2D array list
   cout <<"start mesh init"<<endl;

   height=L[1];

    Cnode n;
    QList <Cnode> col;


    for(double y=-L[1]/2.;y<L[1]/2.;y+=step)    col.push_back(n);
    for(double x=-L[0]/2.;x<L[0]/2.;x+=step)   node.push_back(col);


    for(int i=0;i<node.size();i++)
        for(int j=0;j<node[i].size();j++)
    {
        node[i][j].X[0] = i*step-L[0]/2;
        node[i][j].X[1] = j*step-L[1]/2;
    }
    for(int i=1;i<node.size()-1;i++)
        for(int j=1;j<node[i].size()-1;j++)
        {
          node[i][j].left = &node[i-1][j];
            node[i][j].right = &node[i+1][j];
              node[i][j].bot = &node[i][j-1];
                node[i][j].top = &node[i][j+1];

        }


    // initialise the Yprofile list
    for(double y=-L[1]/2.;y<L[1]/2.;y+=step)
    {
        Cnode m;
        m.X[1] = y;
        Yprofile.push_back(m);
    }
}

void Cmesh::get_avergage(QList <Cgrain> grain)
{

    for(int i=0;i<grain.size();i++)
    {

        for(int a=0;a<node.size();a++)
            for(int b=0;b<node[a].size();b++)
                node[a][b].add_grain(grain[i]);
    }

    for(int a=0;a<node.size();a++)
        for(int b=0;b<node[a].size();b++)
            node[a][b].average();

    for(int i=1;i<node.size()-1;i++)
        for(int j=1;j<node[i].size()-1;j++)
        {
            node[i][j].shear_rate = (node[i][j+1].V.x[0] - node[i][j-1].V.x[0])/ (node[i][j+1].X[1]-node[i][j-1].X[1]);

        }




    //perform the averge along the y axis
    for(int a=0;a<node.size();a++)
        for(int b=0;b<node[a].size();b++)
            Yprofile[b]+=node[a][b];

    for(int y=0;y<Yprofile.size();y++)
        Yprofile[y]/=Yprofile[y].weight;

   // for(int y=1;y<Yprofile.size()-1;y++)
     //   Yprofile[y].shear_rate= (Yprofile[y+1].V[0]-Yprofile[y-1].V[0])/(Yprofile[y+1].X[1]-Yprofile[y-1].X[1]);


}

void Cmesh::operator += (Cmesh m)
{
    for(int y=0;y<Yprofile.size();y++)
        Yprofile[y]+=m.Yprofile[y];
}
void Cmesh::operator /= (double d)
{
    if(d==0)return;
    for(int y=0;y<Yprofile.size();y++)
        Yprofile[y]/=d;
}

void Cmesh::save(string file_profile, string file_mesh)
{
  FILE *fp;

  fp=fopen(file_profile.c_str(),"w");
  foreach(Cnode n,Yprofile) if( fabs(n.X[1])<=height/2.-UNIT_LENGHT) n.PRINT(fp);
  fclose(fp);

  fp=fopen(file_mesh.c_str(),"w");

  for(int a=0;a<node.size();a++)
      for(int b=0;b<node[a].size();b++)
          if(node[a][b].X[1]<=height/2.) node[a][b].PRINT(fp);
  fclose(fp);

}

