#ifndef VECTOR_H
#define VECTOR_H
#define DIM 3
#define DIMuse 3
#include "main.h"

class Cmatrix;
class Cvector;


class Cmatrix
{
public:
    double x[DIM][DIM];
    Cmatrix(){for(int i=0;i<DIM;i++) for(int j=0;j<DIM;j++) x[i][j]=0;}

    //double norm(){double n=0; for(int i=0;i<DIMuse;i++) n+=x[i]*x[i];return sqrt(n);}//return the norm of the vector
    Cmatrix operator+(Cmatrix a){ Cmatrix res; for(int i=0;i<DIMuse;i++)for(int j=0;j<DIMuse;j++)res.x[i][j]=x[i][j]+a.x[i][j];return res; };
    Cmatrix operator-(Cmatrix a){ Cmatrix res; for(int i=0;i<DIMuse;i++)for(int j=0;j<DIMuse;j++)res.x[i][j]=x[i][j]-a.x[i][j];return res; };

    void operator+=(Cmatrix a){ for(int i=0;i<DIMuse;i++) for(int j=0;j<DIMuse;j++)x[i][j]+= a.x[i][j]; };
    void operator-=(Cmatrix a){ for(int i=0;i<DIMuse;i++) for(int j=0;j<DIMuse;j++)x[i][j]-= a.x[i][j]; };
    void operator*=(double d){ for(int i=0;i<DIMuse;i++) for(int j=0;j<DIMuse;j++) x[i][j]*= d; };
    void operator/=(double d){ for(int i=0;i<DIMuse;i++) for(int j=0;j<DIMuse;j++) x[i][j]/= d; };

    Cmatrix operator*(double d){Cmatrix res; for(int i=0;i<DIMuse;i++) for(int j=0;j<DIMuse;j++) res.x[i][j]=x[i][j]*d;return res; };
    Cmatrix operator/(double d){Cmatrix res; for(int i=0;i<DIMuse;i++) for(int j=0;j<DIMuse;j++) res.x[i][j]=x[i][j]/d;return res; };

    Cmatrix operator*(Cmatrix a){Cmatrix res; for(int i=0;i<DIMuse;i++)  for(int j=0;j<DIMuse;j++)for(int k=0;k<DIMuse;k++)   res.x[i][j]+=x[i][k]*a.x[k][j];   return res; };



    double sym(){return( (x[0][1]+x[1][0])/2.);};
    double asym(){return( (x[0][1]-x[1][0])/2.);};

    Cmatrix inverse()
    {
        Cmatrix inv;
        double det = x[0][0]*x[1][1]-x[1][0]*x[0][1];
        if(det==0)return inv;

        inv.x[0][0]= x[1][1];
        inv.x[0][1]=-x[0][1];
        inv.x[1][0]=-x[1][0];
        inv.x[1][1]= x[0][0];

        return(inv/det);
    }
    void PRINT() {for(int i=0;i<DIM;i++)for(int j=0;j<DIM;j++)cout<<x[i][j]<<"\t";cout<<endl;} /*<< Print in to screen*/
    void PRINT(FILE *fp) {for(int i=0;i<DIM;i++)for(int j=0;j<DIM;j++)fprintf(fp,"%.6e\t",x[i][j]);} /*<< Print in file*/
    void READ(FILE *fp) {for(int i=0;i<DIM;i++)for(int j=0;j<DIM;j++)fscanf(fp,"%lf",&x[i][j]);} /*<< Read in file*/
    void PRINT2d(FILE *fp) {for(int i=0;i<2;i++)for(int j=0;j<2;j++)fprintf(fp,"%.8e\t",x[i][j]);} /*<< Print in file*/
    void READ2d(FILE *fp) {for(int i=0;i<2;i++)for(int j=0;j<2;j++)fscanf(fp,"%lf",&x[i][j]);} /*<< Read in file*/


};

class Cvector
{
public:
  double x[DIM];

  Cvector(){for(int i=0;i<DIM;i++) x[i]=0;}
  double norm(){double n=0; for(int i=0;i<DIMuse;i++) n+=x[i]*x[i];return sqrt(n);}//return the norm of the vector
  Cvector operator+(Cvector a){ Cvector res; for(int i=0;i<DIMuse;i++)res.x[i]=x[i]+a.x[i];return res; };
  Cvector operator-(Cvector a){ Cvector res; for(int i=0;i<DIMuse;i++)res.x[i]=x[i]-a.x[i];return res; };
  void operator+=(Cvector a){ for(int i=0;i<DIMuse;i++) x[i]+= a.x[i]; };
  void operator-=(Cvector a){ for(int i=0;i<DIMuse;i++) x[i]-= a.x[i]; };
  void operator*=(double d){ for(int i=0;i<DIMuse;i++) x[i]*= d; };
  void operator/=(double d){ for(int i=0;i<DIMuse;i++) x[i]/= d; };

  Cvector operator*(double d){Cvector res; for(int i=0;i<DIMuse;i++) res.x[i]=x[i]*d;return res; };
  Cvector operator/(double d){Cvector res; for(int i=0;i<DIMuse;i++) res.x[i]=x[i]/d;return res; };
  Cmatrix operator^(Cvector b){Cmatrix m;for(int i=0;i<DIMuse;i++)for(int j=0;j<DIMuse;j++)m.x[i][j]=x[i]*b.x[j];return(m); }
//  Cvector operator*(Cmatrix m){Cvector res; for(int i=0;i<DIMuse;i++)  for(int k=0;k<DIMuse;k++)   res.x[i]+=m.x[i][k]*x[k];   return(res); };


  double operator*(Cvector a){double scalar=0;for(int i=0;i<DIMuse;i++) scalar+=x[i]*a.x[i]; return scalar;}/*<< Scalar product*/
  double &operator[](int i){return x[i];}


  void READ(){for(int i=0;i<DIM;i++) cin>>x[i];}; /*<< Read from screen the 3 componnents*/
  void PRINT(){for(int i=0;i<DIM;i++) cout<<x[i]<<"\t"; cout<<endl;}; /*<< Print on screen the 3 componnents*/
  void PRINT(FILE *fp) {fprintf(fp,"%.6e\t%.6e\t%.6e\t",x[0], x[1], x[2]);} /*<< Print in file*/
  void READ(FILE *fp) {fscanf(fp,"%lf %lf %lf\t",&x[0],&x[1], &x[2]);} /*<< Read in file*/
  void PRINT2d(FILE *fp) {fprintf(fp,"%.8e\t%.8e\t",x[0], x[1]);} /*<< Print in file*/
  void READ2d(FILE *fp) {fscanf(fp,"%lf %lf\t",&x[0],&x[1]);} /*<< Read in file*/


};


#endif // VECTOR_H
