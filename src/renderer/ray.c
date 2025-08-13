#include "ray.h"
#include "../math/vec4.h"
#include "ray.repr.h"
#include <math.h>
#include <string.h>
#include <zot.h>

ray *make_ray(ray *ray, vec4 *origin, vec4 *direction) {
  ray = ray ? ray : zcalloc(1, sizeof(struct ray));

  ray->direction = direction;
  ray->origin = origin;
  ray->distance = INFINITY;
  return ray;
}

vec4 *ray_hit_point(vec4 *v, ray *ray) {
  vec4 *dist_dir = vscale(NULL, ray->distance, ray->direction);
  return vadd(v, dist_dir, ray->origin);
}
