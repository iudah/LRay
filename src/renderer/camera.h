#ifndef CAMERA_H
#define CAMERA_H

#include "scene.h"
#include <inttypes.h>

typedef struct camera camera;

camera *make_camera(camera *cam, int image_width, int image_height, float AOV);
bool render(camera *cam, scene *scene);
char *camera_dim_string(camera *camera, char *str, uint64_t bufflen);

#endif