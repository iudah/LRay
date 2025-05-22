#ifndef OBJECT_H
#define OBJECT_H

#include "math/vec4.h"
#include "ray.h"

typedef struct object object;

vec4 *trace_color(vec4 *color, ray *ray, object *object);
object *make_sphere(object *object, float radius, vec4 *center, vec4 *color);
bool intersect_object(ray *ray, object *object);

#endif
