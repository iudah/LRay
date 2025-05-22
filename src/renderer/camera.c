
#include "camera.h"
#include "math/vec4.h"
#include "ray.h"
#include "scene.h"
#include <math.h>
#include <stdlib.h>
#include <zot.h>

#define STB_IMAGE_WRITE_IMPLEMENTATION 1
#include "stb_image/stb_image_write.h"

typedef enum { m, cm, mm } world_unit;
struct camera {
  mat4 *camera_to_world;
  float AOV; // Measured in degree
  float FOV;
  float aspect_ratio_x, aspect_ratio_y;
  int width, height; // Measured in pixel
};

// float unit_per_mm[] = {[m] = 0.001, [mm] = 1};

camera *make_camera(camera *cam, int image_width, int image_height, float AOV) {
  cam = cam ? cam : zcalloc(1, sizeof(struct camera));

  cam->width = image_width;
  cam->height = image_height;

  cam->aspect_ratio_x = (float)cam->width / cam->height;
  cam->aspect_ratio_y = (float)cam->height / cam->width;

  cam->AOV = AOV;
  cam->FOV = tanf(cam->AOV / 2. * M_PI / 180.);

  cam->camera_to_world = make_mat4(cam->camera_to_world, NULL);
  mtranslate(cam->camera_to_world, 0, 0, 5);
  mrotate_z(cam->camera_to_world, 29);

  return cam;
}

ray *make_primary_ray(ray *ray, camera *cam, int x, int y) {

  float NDC_x = (x + 0.5f) / cam->width;
  float NDC_y = (y + 0.5f) / cam->height;
  float cam_x = (2.f * NDC_x - 1.f) * cam->FOV * cam->aspect_ratio_x;
  float cam_y = (1.f - 2.f * NDC_y) * cam->FOV; //* cam->aspect_ratio_y;

  vec4 *dir = make_vec4(NULL, (float[]){cam_x, cam_y, -1, 0});
  dir = v3mdot(NULL, dir, cam->camera_to_world);

  vec4 *orig = make_vec4(NULL, (float[]){0, 0, 0, 1});
  orig = vmdot(NULL, orig, cam->camera_to_world);

  vnorm(dir, dir);

  return make_ray(ray, orig, dir);
}

#if 1
bool render(camera *cam, scene *scene) {
  char *image = zcalloc(cam->height * cam->width * 3, sizeof(char));
  for (int y = 0; y < cam->height; y++) {
    for (int x = 0; x < cam->width; x++) {
      ray *primary = make_primary_ray(NULL, cam, x, y);
      vec4 *color = make_vec4(NULL, (float[]){0, 0, 0, 1});

      scan_pixel(color, scene, primary);
      char *pix = &image[(y * cam->width + x) * 3];
      pix[0] = (char)(255 * ((float *)color)[0]);
      pix[1] = (char)(255 * ((float *)color)[1]);
      pix[2] = (char)(255 * ((float *)color)[2]);
    }
  }

  stbi_write_jpg("/mnt/sdcard/Jay/Projects/LRay/img7.jpg", cam->width,
                 cam->height, 3, image, 100);

  return true;
}
#else
#include "ray.repr.h"
struct v4 {
  float x, y, z, w;
};

bool render(camera *cam, scene *scene) {
  unsigned char *image = calloc(cam->width * cam->height * 3, 1);
  for (int y = 0; y < cam->height; y++) {
    for (int x = 0; x < cam->width; x++) {
      int idx = (y * cam->width + x) * 3;
      // light gray background
      image[idx + 0] = image[idx + 1] = image[idx + 2] = 200;

      ray *primary = make_primary_ray(NULL, cam, x, y);
      vec4 color = {0, 0, 0, 1};

      // debug: log first few rays
      if (x < 3 && y < 3) {
        printf("ray[%d,%d] dir = %f %f %f\n", x, y,
               ((float *)primary->direction)[0],
               ((float *)primary->direction)[1],
               ((float *)primary->direction)[2]);
        fflush(stdout);
      }
      scan_pixel(&color, scene, primary);
      if (primary->distance < INFINITY) {
        printf("hit pixel (%d,%d)\n", x, y);
        fflush(stdout);
      }
      // free(primary);
      // write color into image
      image[idx + 0] = (unsigned char)(255 * color.x);
      image[idx + 1] = (unsigned char)(255 * color.y);
      image[idx + 2] = (unsigned char)(255 * color.z);
    }
  }

  int ok = stbi_write_jpg("img7.jpg", cam->width, cam->height, 3, image, 100);
  printf("wrote img7.jpg? %d\n", ok);
  fflush(stdout);
  free(image);
  return ok;
}
#endif

char *camera_dim_string(camera *camera, char *str, uint64_t bufflen) {
  snprintf(str, bufflen, "%i by %i", camera->width, camera->height);
  return str;
}