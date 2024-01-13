#include <cmath>

namespace LRay {
// A class declaration representing 3d object
class _3DObject;

// A class declaration representing camera
class Camera;

// A class declaration to represent vector of 4 dimensions
class Vec4;

// A class declaration to represent vector of 3 dimensions
class Vec3;

// A ray
class Ray;
}; // namespace LRay

// Vec4
class LRay::Vec4 {
protected:
  float x;
  float y;
  float z;
  float w;

public:
  // create a 4d vector i.e. A =(x,y,z,w)
  Vec4(float x = 0, float y = 0, float z = 0, float w = 0);

  // create a copy of a 4d vector i.e. B = A
  Vec4(Vec4 &vec4);
  
  // Add two vectors i.e. C = A + B => C_x = A_x + B_x; C_y = A_y + B_y ...
  LRay::Vec4 add(LRay::Vec4 &b);
  
  // Minus two vectors i.e. C = A - B => C_x = A_x - B_x; C_y = A_y - B_y ...
  LRay::Vec4 minus(LRay::Vec4 &b);
  
  // Multiply two vectors i.e. C = A * B => C_x = A_x * B_x; C_y = A_y * B_y ...
  LRay::Vec4 times(LRay::Vec4 &b);
  
  // Multiply vector by scalar s i.e. C = A * s
  
  // => C_x = A_x * s; C_y = A_y * s...
  LRay::Vec4 times(float s);
  
  // Dot product of two vectors
  float dot(LRay::Vec4 &b);
  
  // Angle btwn two vectors
  float angle(LRay::Vec4 &b){
    float result = acos(this->dot(b)/(this->mag() * b.mag()));
    return result;
  };
  
  // magnitude of vector
  float mag();
  
  // compute the unit vector
  LRay::Vec4 unit();
  
  // Divide vector by scalar s i.e. C = A / s
  // => C_x = A_x / s; C_y = A_y / s...
  LRay::Vec4 divide(float s);
  
  // override division sign to divide by scalar
  LRay::Vec4 operator/(float s);
};

LRay::Vec4::Vec4(float x, float y, float z, float w) : x(x), y(y), z(z), w(w) {}

LRay::Vec4 LRay::Vec4::add(LRay::Vec4 &b) {
  return LRay::Vec4(this->x + b.x, this->y + b.y, this->z + b.z, this->w + b.w);
}
LRay::Vec4 LRay:: Vec4::minus(LRay::Vec4 &b){
  float x = this -> x - b.x;
  float y = this -> y - b.y;
  float z = this -> z - b.z;
  float w = this -> w - b.w;
  return LRay::Vec4(x,y,z,w);
}
LRay::Vec4 LRay::Vec4::unit() { return (*this) / this->mag(); }

class LRay::Vec3 : public LRay::Vec4 {
  Vec3(float x = 0, float y = 0, float z = 0);
  Vec3(Vec3 &vec3);
};

class LRay::_3DObject {};

class LRay::Ray {};

class LRay::Camera {
private:
  LRay::Vec4 position;
  const int width;
  const int height;

public:
  Camera(int width = 720, int height = 960);

  // Render image
  int render();
};

int main() {
  LRay::Camera camera_0 = LRay::Camera();
  camera_0.render();
}
