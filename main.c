#include "camera.h"
#include "light.h"
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
        make_vec4(nullptr,
                  (float[]){-5. / 2 + i + .42, -5. / 2 + i + .42, 0, 1}),
        make_vec4(nullptr, (float[]){1, (random() / (float)RAND_MAX),
                                     (random() / (float)RAND_MAX), 1}));
    scene_add_object(scene, sphere);
   
    light *light =
        make_light(nullptr, POINT_LIGHT, 5,
                   make_vec4(nullptr, (float[]){-5. / 2 + i + .42,
                                                -5. / 2 + i + .42, .50, 1}),
                   make_vec4(nullptr, (float[]){1, 1, 1, 1}));
   scene_add_light(scene, light);
 }


  render(camera, scene);

  return 0;
}
