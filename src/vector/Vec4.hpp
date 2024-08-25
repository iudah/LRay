#ifndef VEC4_HPP
#define VEC4_HPP

#include <ostream>
namespace LRay {

// A class declaration to represent vector of 4 dimensions
// Vec4
class Vec4 {
protected:
  float x;
  float y;
  float z;
  float w;

public:
  // create a 4d vector i.e. A =(x,y,z,w)
  explicit Vec4(float x = 0, float y = 0, float z = 0, float w = 0);

  // create a copy of a 4d vector i.e. B = A
  Vec4(const Vec4 &vec4);

  // Add two vectors i.e. C = A + B => C_x = A_x + B_x; C_y = A_y + B_y ...
  LRay::Vec4 add(const LRay::Vec4 &b) const;

  // Minus two vectors i.e. C = A - B => C_x = A_x - B_x; C_y = A_y - B_y ...
  LRay::Vec4 minus(const LRay::Vec4 &b) const;

  // Multiply two vectors i.e. C = A * B => C_x = A_x * B_x; C_y = A_y * B_y ...
  LRay::Vec4 times(const LRay::Vec4 &b);

  // Multiply vector by scalar s i.e. C = A * s
  // => C_x = A_x * s; C_y = A_y * s...
  LRay::Vec4 times(const float s) const;

  // override division sign to divide by scalar
  LRay::Vec4 operator/(const float s) const;
  LRay::Vec4 operator+(const Vec4 b) const;
  LRay::Vec4 operator-(const Vec4 b) const;

  // Dot product of two vectors
  float dot(const LRay::Vec4 &b) const;

  // Angle btwn two vectors
  float angle(const LRay::Vec4 &b) const;

  // magnitude of vector
  float mag() const;

  // compute the unit vector
  LRay::Vec4 unit() const;

  // Divide vector by scalar s i.e. C = A / s
  // => C_x = A_x / s; C_y = A_y / s...
  LRay::Vec4 divide_by(const float s) const;

  // to string
  friend std::ostream &operator<<(std::ostream &os, const LRay::Vec4 &vec4);

  float get_x() const { return x; }
  float get_y() const { return y; }
  float get_z() const { return z; }
  float get_w() const { return w; }
};

std::ostream &operator<<(std::ostream &os, const LRay::Vec4 &vec4);

} // namespace LRay
#endif
