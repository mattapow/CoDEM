#ifndef BOX_H
#define BOX_H
#include "main.h"
#include "grain.h"
class Cgrain;

class Cbox
{
public:
  Cbox();
  Cbox *right;  /*<<Pointer to the box on the right */
  Cbox *left;   /*<<Pointer to the box on the left */
  QList <Cgrain *> grain; /*<<List of pointer to the grains within that box*/
};

#endif // BOX_H
