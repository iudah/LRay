#ifndef OBJECT_H
#define OBJECT_H

#include "ray.h"
#include "vec4.h"

typedef struct object object;

vec4 *trace_color(vec4 *color, ray *ray, object *object);
object *make_sphere(object *object, float radius, vec4 *center, vec4 *color);
bool intersect_object(ray *ray, object *object);
vec4 *get_object_hit_normal(vec4 *v, object *hit_object, ray *primary_ray,
                            vec4 *hit_point);

#endif
