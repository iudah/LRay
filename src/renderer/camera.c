
#include "camera.h"
#include "ray.h"
#include "scene.h"
#include "vec4.h"
#include <math.h>
#include <stdlib.h>
#include <zot.h>

#define STB_IMAGE_WRITE_IMPLEMENTATION 1
#include "stb_image_write.h"

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

ray *make_primary_ray(ray *ray, camera *cam, float x, float y) {

  float NDC_x = (x) / cam->width;
  float NDC_y = (y) / cam->height;
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
  uint8_t *image = zcalloc(cam->height * cam->width * 3, sizeof(char));
  uint8_t n_samples = 8;
  uint8_t j_max = (uint8_t)sqrtf(n_samples);
  uint8_t k_max = n_samples / j_max;

  for (int y = 0; y < cam->height; ++y) {
    for (int x = 0; x < cam->width; ++x) {

      uint8_t *pix = &image[(y * cam->width + x) * 3];

      uint32_t pxl[] = {0, 0, 0};

      for (int j = 0; j < j_max; ++j) {
        for (int k = 0; k < k_max; ++k) {
          ray *primary = make_primary_ray(NULL, cam, (float)x + (float)j / j_max,
                                          (float)y + (float)k / k_max);
          vec4 *color = make_vec4(NULL, (float[]){0, 0, 0, 1});

          scan_pixel(color, scene, primary);

#define CLAMP(x) ((x) < 0 ? 0 : (x) > 1 ? 1 : (x))

          pxl[0] += (uint8_t)(255 * CLAMP(((float *)color)[0]));
          pxl[1] += (uint8_t)(255 * CLAMP(((float *)color)[1]));
          pxl[2] += (uint8_t)(255 * CLAMP(((float *)color)[2]));
        }
      }
      pix[0] = (uint8_t)(pxl[0] / n_samples);
      pix[1] = (uint8_t)(pxl[1] / n_samples);
      pix[2] = (uint8_t)(pxl[2] / n_samples);
    }
  }

  stbi_write_jpg("/mnt/sdcard/Jay/lray/img7.jpg", cam->width, cam->height, 3,
                 image, 100);

  return true;
}
#endif

char *camera_dim_string(camera *camera, char *str, uint64_t bufflen) {
  snprintf(str, bufflen, "%i by %i", camera->width, camera->height);
  return str;
}