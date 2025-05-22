#ifndef RAY_H
#define RAY_H

#include "math/vec4.h"

typedef struct ray  ray;

ray *make_ray(ray *ray, vec4 *origin, vec4 *direction);

#endif