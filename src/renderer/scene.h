#ifndef SCENE_H
#define SCENE_H

#include "math/vec4.h"
#include "object.h"
#include "ray.h"

typedef struct scene scene;

scene *_Nullable make_scene(scene *_Nullable scene);
bool scene_add_object(scene *_Nonnull scene, object *_Nonnull object);
vec4 *_Nullable scan_pixel(vec4 *_Nonnull color, scene *_Nonnull scene,
                           ray *_Nonnull ray);

#endif