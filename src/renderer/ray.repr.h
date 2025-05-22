#ifndef RAY_REPR_H
#define RAY_REPR_H

#include "math/vec4.h"

struct ray {
  vec4 *origin, *direction;
  float distance;
};

#endif