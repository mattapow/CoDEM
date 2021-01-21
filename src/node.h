#ifndef NODE_H
#define NODE_H
#include "main.h"
#include "grain.h"

class Cnode
{
public:
    Cnode();
    Cvector X,V;        /**<POsition and CG velocity */
    double phi;         /**<Solid fraction */
    double shear_rate;  /**<Shear rate */
    Cmatrix stress;     /**<Stress tensor */

    double weight;      /**<Total weight of grains */

    Cnode *top,*bot,*right,*left;

    double get_weigth(Cgrain &g);
    void add_grain(Cgrain &g);
    void average();

    void operator +=(Cnode n);
    void operator /=(double d);

    void PRINT(FILE *fp);
    void READ(FILE *fp);
};

#endif // NODE_H
