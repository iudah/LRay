#ifndef RAY_REPR_H
#define RAY_REPR_H

#include "vec4.h"

struct ray {
  vec4 *origin, *direction;
  float distance;
};

#endif