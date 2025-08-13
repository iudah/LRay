#include "scene.h"
#include "light.h"
#include "object.h"
#include "ray.h"
#include <stdint.h>
#include <string.h>
#include <zot.h>

#define CAPACITY 20

struct array_ptr {
  void **ptr;
  uint32_t n_ptr;
  uint32_t cap_ptr;
};

struct scene {
  struct array_ptr object;
  struct array_ptr light;
};

void make_array(struct array_ptr *array) {
  array->ptr = zmalloc(CAPACITY * sizeof(void *));
  array->cap_ptr = CAPACITY;
  array->n_ptr = 0;
}

bool array_append(struct array_ptr *array, void *ptr) {
  if (array->cap_ptr == array->n_ptr) {
    void *tmp = zrealloc(array->ptr, array->cap_ptr * 2);
    if (!tmp) {
      return false;
    }
    array->ptr = tmp;
    array->cap_ptr *= 2;
  }
  array->ptr[array->n_ptr++] = ptr;
  return true;
}

scene *make_scene(scene *scene) {
  scene = scene ? scene : zcalloc(1, sizeof(struct scene));
  make_array(&scene->object);
  make_array(&scene->light);
  return scene;
}

bool scene_add_object(scene *_Nonnull scene, object *_Nonnull object) {
  return array_append(&scene->object, object);
}

bool scene_add_light(scene *_Nonnull scene, light *_Nonnull light) {
  return array_append(&scene->light, light);
}

vec4 *scan_pixel(vec4 *color, scene *scene, ray *primary_ray) {

  object *hit_object = NULL;

  for (uint32_t i = 0; i < scene->object.n_ptr; i++) {
    bool ray_intersects = intersect_object(primary_ray, scene->object.ptr[i]);
    if (ray_intersects) {
      hit_object = scene->object.ptr[i];
    }
  }

  if (hit_object) {

    vec4 *light_contrib = make_vec4(NULL, NULL);
     vec4 *hit_point = ray_hit_point(NULL, primary_ray);
   vec4 *hit_normal = get_object_hit_normal(NULL, hit_object, primary_ray, hit_point);

    // for each light in the scene
    for (uint32_t j = 0; j < scene->light.n_ptr; j++) {
      light *light = scene->light.ptr[j];
      // get secondary ray
      struct ray *shadow_ray =
          make_secondary_ray(NULL, primary_ray, light, hit_point, hit_normal);

      bool in_shadow = false;
      for (uint32_t i = 0; i < scene->object.n_ptr; i++) {
        if (scene->object.ptr[i] == hit_object)
          continue;

        in_shadow = intersect_object(shadow_ray, scene->object.ptr[i]);

        if (in_shadow) {
          break;
        }
      }
      if (in_shadow) {
        continue;
      }
      light_contrib =
          vadd(NULL, light_contrib,
               light_contribution(NULL, shadow_ray, light, hit_normal));
    }

    vec4 *object_color = trace_color(NULL, primary_ray, hit_object);
    vmul(color, object_color, light_contrib);
  }

  return color;
}