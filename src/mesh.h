#ifndef MESH_H
#define MESH_H

#include "main.h"
#include "node.h"

class Cmesh
{
public:


    QList <QList <Cnode> > node;/**<Two dimensional array of node */

    QList <Cnode> Yprofile;/**<One dimensional array of node for averaging on x direction*/
    double height;
    Cmesh(Cvector L, double step); /**< Generate the list of nodes */
    void get_avergage(QList <Cgrain> grain);

    void operator += (Cmesh m);
    void operator /=(double d);
    void save(string file_profile,string file_mesh);
};

#endif // MESH_H
