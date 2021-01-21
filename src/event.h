#ifndef EVENT_H
#define EVENT_H

class Cevent
{
public:
  double period;
  double next;
  double config_dt;

  bool should_do(double t){
      if (period>0 && t>(next-config_dt/2)) {
          next+=period;
          return true;
      }
      return false;
  }
};

#endif // EVENT_H
