#ifndef RAY_H
#define RAY_H

#include "vec4.h"

typedef struct ray ray;

ray *make_ray(ray *ray, vec4 *origin, vec4 *direction);
vec4 *ray_hit_point(vec4*v, ray *ray);

#endif