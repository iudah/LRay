#include "object.h"
#include "ray.repr.h"
#include <math.h>
#include <zot.h>

struct object {
  bool (*intersect_fn)(ray *ray, object *sphere);
  vec4 *color;
};

struct sphere {
  struct object object;
  vec4 *center;
  float radius;
};

bool solve_quadratic(float *t1, float *t2, float a, float b, float c) {
  // https://scratchapixel.com/lessons/3d-basic-rendering/minimal-ray-tracer-rendering-simple-shapes/ray-sphere-intersection.html
  auto D = b * b - 4 * a * c;
  if (D < 0)
    return false;
  else if (D == 0) {
    *t1 = *t2 = -0.5 * b / a;
  } else {
    auto q = -0.5f * ((b > 0) ? b + sqrtf(D) : b - sqrtf(D));
    *t1 = q / a;
    *t2 = c / q;
    if (*t1 > *t2) {
      auto tmp = *t1;
      *t1 = *t2;
      *t2 = tmp;
    }
  }

  return true;
}

bool intersect_sphere(struct ray *ray, object *object) {
  struct sphere *sphere = (void *)object;

  auto a = 1;
  auto O_minus_C = vsub(NULL, ray->origin, sphere->center);
  auto b = 2 * vdot(NULL, ray->direction, O_minus_C);
  auto c = vdot(NULL, O_minus_C, O_minus_C) - sphere->radius * sphere->radius;

  float t_1, t_2;
  bool sphere_is_hit = solve_quadratic(&t_1, &t_2, a, b, c);
  if (sphere_is_hit) {
    if (t_1 < 0) {
      t_1 = t_2;
      if (t_2 < 0) {
        return false;
      }
    }

    if (t_1 < ray->distance) {
      ray->distance = t_1;

      return true;
    }
  }

  return false;
}

object *make_sphere(object *object, float radius, vec4 *center, vec4 *color) {
  struct sphere *sphere = object ? object : zcalloc(1, sizeof(struct sphere));
  sphere->object.intersect_fn = (void *)intersect_sphere;
  sphere->object.color = color;
  sphere->center = center;
  sphere->radius = radius;

  return (struct object *)sphere;
}

// bool intersect_triangle(ray *ray, void *triangle) { return false; }

bool intersect_object(ray *ray, object *object) {
  return object->intersect_fn(ray, object);
}

vec4 *trace_color(vec4 *color, ray *ray, object *object) {
  // auto hit_point = vscale(NULL, ray->distance, ray->direction);
  make_vec4(color, (float *)object->color);
  printf(".");
  return color;
}
