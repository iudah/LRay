#include "light.h"
#include "ray.repr.h"
#include <math.h>
#include <stdint.h>
#include <string.h>
#include <zot.h>

struct light {
  ray *(*secondary_ray_fn)(ray *ray_res, ray *primary_ray, light *light,
                           vec4 *hit_point, vec4 *hit_normal);
  vec4 *color;
  vec4 *center;
  float brightness;
};

ray *make_point_secondary_ray(ray *ray_res, ray *primary_ray, light *light,
                              vec4 *hit_point, vec4 *hit_normal) {
  vec4 *offset_point = vadd(NULL, hit_point, vscale(NULL, EPSILON, hit_normal));

  vec4 *light_to_hitpoint = vsub(NULL, light->center, hit_point);
  // point light
  ray_res = make_ray(ray_res, offset_point, vnorm(NULL, light_to_hitpoint));
  ray_res->distance = vmag(NULL, light_to_hitpoint);

  return ray_res;
}

light *make_light(light *light, uint8_t type, float brightness, vec4 *center,
                  vec4 *color) {
  struct light *l = light ? light : zcalloc(1, sizeof(struct light));
  l->secondary_ray_fn = (void *)make_point_secondary_ray;
  l->color = color;
  l->center = center;
  l->brightness = brightness;

  return (struct light *)l;
}

ray *make_secondary_ray(ray *ray_res, ray *primary_ray, light *light,
                        vec4 *hit_point, vec4 *hit_normal) {
  return light->secondary_ray_fn(ray_res, primary_ray, light, hit_point,
                                 hit_normal);
}

vec4 *light_contribution(vec4 *color, ray *light_ray, light *light,
                         vec4 *hit_normal) {
  float light_factor = fmaxf(0, vdot(NULL, light_ray->direction, hit_normal));
  float r =light_ray->distance;// vmag(NULL, vsub(NULL, light_ray->origin, light->center));
// light->brightness

  return vscale(color,  light_factor / r / r, light->color);
}
