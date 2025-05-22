#include "camera.h"
#include "object.h"
#include "scene.h"
#include <stdio.h>
#include <stdlib.h>

int main() {
  camera *camera = make_camera(nullptr, 240, 320, 51.52);
  printf("Camera(%p): %s\n", camera,
         camera_dim_string(camera, (char[24]){0}, 24));

  scene *scene = make_scene(nullptr);

  for (int i = 0; i < 5; i++) {
    object *sphere = make_sphere(
        nullptr, .21,
        make_vec4(nullptr, (float[]){1, 2 - 4 * (random() / (float)RAND_MAX),
                                     5 - 10 * (random() / (float)RAND_MAX), 1}),
        // make_vec4(nullptr, (float[]){0, 1, 1, 1}),
        make_vec4(nullptr, (float[]){1, 0, random() / (float)RAND_MAX, 1}));
    scene_add_object(scene, sphere);
  }

  object *sphere =
      make_sphere(nullptr, .51, make_vec4(nullptr, (float[]){0, 0, 0, 1}),
                  make_vec4(nullptr, (float[]){.7, .7, .7, 1}));
  scene_add_object(scene, sphere);

  render(camera, scene);

  return 0;
}
