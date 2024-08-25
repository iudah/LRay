#include <cmath>

#include "Vec4.hpp"

// namespace LRay
namespace LRay {

    //vec4 constructor
    Vec4::Vec4(float x, float y, float z, float w): x(x),
    y(y),
    z(z),
    w(w) {}

    Vec4::Vec4(const Vec4 &vec4): x(vec4.x),
    y(vec4.y),
    z(vec4.z),
    w(vec4.w) {}

    //sum of vector
    Vec4 Vec4::add(const Vec4 &b) const {
        return LRay::Vec4(x + b.x, y + b.y, z + b.z, w + b.w);
    }

    //difference between two vectors.
    Vec4 Vec4::minus(const Vec4 &b) const {
        return Vec4(x - b.x, y - b.y, z - b.z, w - b.w);
    }

    Vec4 Vec4::unit() const {
        return (*this) / mag();
    }

    //function to calculate angle between vector
    float Vec4::angle(const Vec4 &b) const {
        float result = std::acos(dot(b) / (mag() * b.mag()));
        return result;
    }

    float Vec4::mag() const {
        return std::sqrt(x * x + y * y + z * z + w * w);
    }

    Vec4 Vec4::operator/(const float s) const {
        return (*this).divide_by(s);
    }

    // Dot product of two vectors
    float Vec4::dot(const Vec4 &b) const {
        return x * b.x + y * b.y + z * b.z + w * b.w;
    }

    Vec4 Vec4::divide_by(const float s) const {
        return LRay::Vec4(x / s, y / s, z / s, w / s);
    }

    std::ostream &operator<<(std::ostream &os, const LRay::Vec4 &vec4) {
        os << "(" << vec4.get_x() << ", " << vec4.get_y() << ", " << vec4.get_z()
        << ", " << vec4.get_w() << ")";
        return os;
    }
    Vec4 Vec4::operator-(const Vec4 b) const {
        return minus(b);
    }

} // namespace LRay