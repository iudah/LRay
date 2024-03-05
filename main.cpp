#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <utility>
#include <vector>

#include "Vec4.hpp"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

namespace LRay {
// A class declaration representing 3d object
class TDObject;

// A class declaration representing camera
class Camera;

// A class declaration to represent vector of 4 dimensions
class Vec4;

// A class declaration to represent vector of 3 dimensions
class Vec3;

// A ray
class Ray;
// the 3d scene

class Scene;

class Sphere;

std::ostream &operator<<(std::ostream &os, const LRay::Ray &ray);
std::ostream &operator<<(std::ostream &os, const LRay::TDObject &object);

}; // namespace LRay

class LRay::Vec3 : public LRay::Vec4 {
public:
  Vec3(float x = 0, float y = 0, float z = 0);
  Vec3(const LRay::Vec4 &vec4);
};

class LRay::TDObject {
public:
  virtual bool intersect(LRay::Ray &ray) const;
  virtual LRay::Vec3 trace_color(LRay::Ray &ray) const;
  friend std::ostream &operator<<(std::ostream &os,
                                  const LRay::TDObject &object);
};

class LRay::Sphere : public LRay::TDObject {

protected:
  LRay::Vec3 origin;
  float radius;

public:
  Sphere(const Vec3 &origin = Vec3(), const float radius = 1);
  bool intersect(LRay::Ray &ray) const override;
  LRay::Vec3 trace_color(LRay::Ray &ray) const override;
  friend std::ostream &operator<<(std::ostream &os, const LRay::Sphere &sphere);
};

class LRay::Ray {
  LRay::Vec3 origin;
  LRay::Vec3 direction;
  float distance = INFINITY;

public:
  Ray(LRay::Vec4 origin, LRay::Vec4 dir);
  float get_distance() const { return distance; };
  const LRay::Vec3 &get_origin() const { return origin; }
  const LRay::Vec3 &get_direction() const { return direction; }
  bool set_distance(float hit_distance);
  friend std::ostream &LRay::operator<<(std::ostream &os, const LRay::Ray &ray);
};

class LRay::Scene : public std::vector<LRay::TDObject *> {
public:
  LRay::Vec3 trace(LRay::Ray &ray) const;
};

class LRay::Camera {
private:
  LRay::Vec4 position;
  LRay::Scene scene;
  const int width;
  const int height;
  const float cm_per_px;
  const float focal_length;

public:
  Camera(int width = 720, int height = 960, float focal_length = 3.5 /*cm*/,
         int PPI = 300);

  // Render image
  int render() const;
  LRay::Ray make_camera_ray(int x, int y) const;
  void set_scene(LRay::Scene &scene);
};

LRay::Camera::Camera(int width, int height, float focal_length, int ppi)
    : width(width), height(height), cm_per_px(2.54 /*cm/in*/ / ppi),
      focal_length(focal_length) {
  std::cout << "Camera: " << width << "x" << height << " @ " << ppi << "PPI\n";
}

int LRay::Camera::render() const {
  std::cout << "Camera rendering scene\n";

  int x_bound[] = {-width / 2, width / 2};
  int y_bound[] = {-height / 2, height / 2};

  unsigned char *image_data =
      (unsigned char *)malloc(width * height * 4 * sizeof(unsigned char));

  unsigned char *current_color = image_data;

    for (int y = y_bound[0]; y < y_bound[1]; y++) {
  for (int x = x_bound[0]; x < x_bound[1]; x++) {
      LRay::Ray camera_ray = make_camera_ray(x, y);

      std::cout << "Camera ray: " << camera_ray << " ";

      LRay::Vec3 color = scene.trace(camera_ray);

      *current_color = (unsigned char)color.get_x();
      current_color++;
      *current_color = (unsigned char)color.get_y();
      current_color++;
      *current_color = (unsigned char)color.get_z();
      current_color++;
    }
  }


  stbi_write_jpg("./camera_test.jpg", width, height, 3, image_data, 100);

  return 0;
}

void LRay::Camera::set_scene(LRay::Scene &scene_of_interest) {
  scene = scene_of_interest;
}

LRay::Ray LRay::Camera::make_camera_ray(int x, int y) const {

  return LRay::Ray(Vec3(),
                   Vec3(x * cm_per_px, y * cm_per_px, -focal_length).unit());
}

LRay::Vec3::Vec3(float x, float y, float z) : LRay::Vec4(x, y, z) {}

LRay::Vec3::Vec3(const LRay::Vec4 &vec4) : LRay::Vec4(vec4) {}

LRay::Ray::Ray(const Vec4 origin, const Vec4 dir)
    : origin(origin), direction(dir) {}

std::ostream &LRay::operator<<(std::ostream &os, const LRay::Ray &ray) {
  os << ray.get_origin() << "->" << ray.get_direction();
  return os;
}

bool LRay::Ray::set_distance(float new_distance) {

  if (new_distance < distance && new_distance > 0) {

    distance = new_distance;
    return true;
  }

  return false;
}

// static bool intersect(LRay::Ray)

LRay::Vec3 LRay::Scene::trace(LRay::Ray &ray) const {

  const LRay::TDObject *hit_object = NULL;
  for (const LRay::TDObject *object : *this) {

    if (object->intersect(ray)) {
      hit_object = (LRay::TDObject *)object;
    }
    // else {
    //   std::cout << object << " isn't hit\n";
    // }
  }

  std::cout << (hit_object ? "hit\n" : "no hit\n");

  return hit_object ? hit_object->trace_color(ray) : LRay::Vec4();
}

static bool solve_quadratic_root(const float a, const float b, const float c,
                                 float &x_1, float &x_2) {
  float d = b * b - 4 * a * c;

  if (d < 0)
    return false;
  else if (d == 0) {
    x_1 = x_2 = -0.5 * b / a;
  } else {
    float q = -0.5 * (b + (b > 0 ? sqrt(d) : -sqrt(d)));
    x_1 = q / a;
    x_2 = c / q;

    if (x_2 < x_1) {
      std::swap(x_1, x_2);
    }

    // std::cout << a << "d^2 + " << b << "d + " << c << " = " << x_1 << ", "
    //           << x_2 << "\n";
  }
  return true;
}

bool LRay::TDObject::intersect(LRay::Ray &ray) const { return false; }
std::ostream &LRay::operator<<(std::ostream &os, const LRay::TDObject &object) {

  os << "LRay::TDObject @ " << &object << " ";
  return os;
}

LRay::Vec3 LRay::TDObject::trace_color(LRay::Ray &ray) const {
  std::cout << "TDObject::trace_color\n";
  return Vec3(100, 10, 20);
}

LRay::Vec3 LRay::Sphere::trace_color(LRay::Ray &ray) const {
  // std::cout << "Sphere::trace_color\n";
  return Vec3(200, 200, 200);
}

bool LRay::Sphere::intersect(LRay::Ray &ray) const {

  LRay::Vec4 D = ray.get_direction();
  LRay::Vec4 O_C = ray.get_origin() - origin;

  float a = 1;
  float b = 2 * O_C.dot(D);
  float c = O_C.dot(O_C) - radius * radius;

  float x_1;
  float x_2;

  if (solve_quadratic_root(a, b, c, x_1, x_2)) {
    if (x_1 > 0)
      return ray.set_distance(x_1);
    if (x_2 > 0)
      return ray.set_distance(x_2);
  }
  return false;
}

LRay::Sphere::Sphere(const LRay::Vec3 &origin, const float radius)
    : origin(origin), radius(radius) {}

std::ostream &LRay::operator<<(std::ostream &os, const LRay::Sphere &sphere) {

  os << "LRay::Sphere @ " << sphere.origin << " ";
  return os;
}
LRay::Scene sphere_scene() {
  //  std::vector<LRay::TDObject> scene;
  LRay::Scene scene;
  //  scene.reserve(2);
  scene.push_back(new LRay::Sphere(LRay::Vec3(0, 0, 5.4)));
  scene.push_back(new LRay::Sphere(LRay::Vec3(0, 0, -5.4)));
  // scene.push_back(new LRay::Sphere(LRay::Vec3(0, 0, 0)));

  // std::cout << "size " << scene.size() << "\n";
  // abort();

  return scene;
  //  return LRay::Scene();
}

int main() {
  LRay::Scene sample_scene = sphere_scene();

  // std::cout<<"sz "<< sample_scene.size()<<"\n";

  //   LRay::Ray ray = LRay::Ray(LRay::Vec3(), LRay::Vec3());

  //   LRay::Sphere sf = LRay::Sphere();
  //   sf.intersect(ray);

  //   LRay::TDObject tdo = sf;
  //   tdo.intersect(ray);

  //   LRay::TDObject *tdd = &tdo;
  //   tdd->intersect(ray);

  //   std::vector<LRay::TDObject> a = {LRay::Sphere(), LRay::TDObject()};
  //   a[0].intersect(ray);

  //   abort();

  LRay::Camera camera_0 = LRay::Camera(320, 400);
  camera_0.set_scene(sample_scene);
  camera_0.render();
}