#ifndef LIGHT_H
#define LIGHT_H

#include "ray.h"
#include "vec4.h"
#include <stdint.h>

typedef struct light light;

enum light_type { POINT_LIGHT, SURFACE_LIGHT };

// vec4 *trace_color(vec4 *color, ray *ray, light *light);
light *make_light(light *light, uint8_t type, float brightness, vec4 *center,
                  vec4 *color);
ray *make_secondary_ray(ray *ray_res, ray *primary_ray, light *light,
                              vec4 *hit_point, vec4 *hit_normal);
vec4 *light_contribution(vec4 *color, ray *light_ray, light *light,
                         vec4 *hit_normal);

#endif
