// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once
#include <cmath>
#include <cstdint>
#include <array>

namespace paacs {

using Index = std::int64_t;
constexpr Index kInvalidIndex = -1;

struct Vec3 {
    double x = 0, y = 0, z = 0;

    Vec3() = default;
    Vec3(double X, double Y, double Z) : x(X), y(Y), z(Z) {}

    Vec3 operator+(const Vec3& o) const { return {x + o.x, y + o.y, z + o.z}; }
    Vec3 operator-(const Vec3& o) const { return {x - o.x, y - o.y, z - o.z}; }
    Vec3 operator*(double s)      const { return {x * s, y * s, z * s}; }
    Vec3 operator/(double s)      const { return {x / s, y / s, z / s}; }

    double dot(const Vec3& o)   const { return x * o.x + y * o.y + z * o.z; }
    Vec3   cross(const Vec3& o) const {
        return {y * o.z - z * o.y, z * o.x - x * o.z, x * o.y - y * o.x};
    }
    double norm2() const { return dot(*this); }
    double norm()  const { return std::sqrt(norm2()); }
    Vec3   normalized() const { double n = norm(); return n > 0 ? *this / n : Vec3{}; }
};

inline Vec3 operator*(double s, const Vec3& v) { return v * s; }

struct Vec2 {
    double u = 0, v = 0;
    Vec2() = default;
    Vec2(double U, double V) : u(U), v(V) {}
    Vec2 operator+(const Vec2& o) const { return {u + o.u, v + o.v}; }
    Vec2 operator-(const Vec2& o) const { return {u - o.u, v - o.v}; }
    Vec2 operator*(double s)      const { return {u * s, v * s}; }
    double norm() const { return std::hypot(u, v); }
};

} // namespace paacs
