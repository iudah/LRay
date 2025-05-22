#include "scene.h"
#include "object.h"
#include <stdint.h>
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

vec4 *scan_pixel(vec4 *color, scene *scene, ray *ray) {

  object *hit_object = NULL;

  for (uint32_t i = 0; i < scene->object.n_ptr; i++) {
    bool ray_intersects = intersect_object(ray, scene->object.ptr[i]);
    if (ray_intersects) {
      hit_object = scene->object.ptr[i];
    }
  }

  if (hit_object) {

    vec4 *new_color = make_vec4(NULL, NULL);
    trace_color(color, ray, hit_object);
  }

  return color;
}