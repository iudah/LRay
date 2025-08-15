#include "object.h"
#include "ray.repr.h"
#include <math.h>
#include <stdint.h>
#include <zot.h>

struct object {
  bool (*intersect_fn)(ray *ray, object *object);
  vec4 *(*hit_normal)(vec4 *v, object *hit_object, ray *primary_ray,
                      vec4 *hit_point);
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
      if (t_1 < 0) {
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

vec4 *get_sphere_hit_normal(vec4 *vres, object *hit_object, ray *primary_ray,
                            vec4 *hit_point) {
  struct sphere *sphere = (void *)hit_object;

  return vnorm(vres, vsub(NULL, hit_point, sphere->center));
}

object *make_sphere(object *object, float radius, vec4 *center, vec4 *color) {
  struct sphere *sphere = object ? object : zcalloc(1, sizeof(struct sphere));
  sphere->object.intersect_fn = (void *)intersect_sphere;
  sphere->object.hit_normal = (void *)get_sphere_hit_normal;
  sphere->object.color = color;
  sphere->center = center;
  sphere->radius = radius;

  return (struct object *)sphere;
}

struct triangle {
  struct object object;
  vec4 *A;
  vec4 *B;
  vec4 *C;
  bool culling;
};

static inline bool moller_trumbore(struct ray *ray, object *object) {
  struct triangle *triangle = (void *)object;

  auto AB = vsub(NULL, triangle->B, triangle->A);
  auto AC = vsub(NULL, triangle->C, triangle->A);

  auto AO = vsub(NULL, ray->origin, triangle->A);
  auto D = ray->direction;

  float det = 0, t = 0, u = 0, v = 0;

  compute_moller_trumbore_unknowns(D, AB, AC, AO, &det, &t, &u, &v);

  if (triangle->culling) {
    if (det < EPSILON)
      return false;
  } else {
    if (fabs(det) < EPSILON)
      return false;
  }

  if (u < 0 || u > 1)
    return false;
  if (v < 0 || (u + v) > 1)
    return false;

  if (t >= ray->distance)
    return false;

  ray->distance = t;

  return true;
}

bool intersect_triangle(struct ray *ray, object *object) {

  return moller_trumbore(ray, object);
}

vec4 *get_triangle_hit_normal(vec4 *vres, object *hit_object, ray *primary_ray,
                              vec4 *hit_point) {
  struct triangle *triangle = (void *)hit_object;

  auto P = hit_point;
  auto AB = vsub(NULL, triangle->B, triangle->A);
  auto AP = vsub(NULL, P, triangle->A);

  return vnorm(vres, vcross(NULL, AB, AP));
}

object *make_triangle(object *object, vec4 *A, vec4 *B, vec4 *C, vec4 *color) {
  struct triangle *triangle =
      object ? object : zcalloc(1, sizeof(struct triangle));
  triangle->object.intersect_fn = (void *)intersect_triangle;
  triangle->object.hit_normal = (void *)get_triangle_hit_normal;
  triangle->object.color = color;
  triangle->A = A;
  triangle->B = B;
  triangle->C = C;

  triangle->culling = true;

  return (struct object *)triangle;
}

bool intersect_object(ray *ray, object *object) {
  return object->intersect_fn(ray, object);
}

vec4 *trace_color(vec4 *color, ray *ray, object *object) {
  // auto hit_point = vscale(NULL, ray->distance, ray->direction);
  color = make_vec4(color, (float *)object->color);
  printf(".");
  return color;
}

vec4 *get_object_hit_normal(vec4 *v, object *hit_object, ray *primary_ray,
                            vec4 *hit_point) {
  return hit_object->hit_normal(v, hit_object, primary_ray, hit_point);
}
