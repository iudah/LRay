#include <cmath>

#include "Vec4.hpp"

namespace LRay {

Vec4::Vec4(float x, float y, float z, float w) : x(x), y(y), z(z), w(w) {}

Vec4 Vec4::add(const Vec4 &b) const {
  return LRay::Vec4(x + b.x, y + b.y, z + b.z, w + b.w);
}

Vec4 Vec4::minus(const Vec4 &b) const {
  return Vec4(x - b.x, y - b.y, z - b.z, w - b.w);
}

Vec4 Vec4::unit() const { return (*this) / mag(); }

float Vec4::angle(const Vec4 &b) const {
  float result = std::acos(dot(b) / (mag() * b.mag()));
  return result;
}


// Dot product of two vectors
  float Vec4::dot(const Vec4 &b) const{
    return (x*b.x + y*b.y + z*b.z + w*b.w);
  }

} // namespace LRay