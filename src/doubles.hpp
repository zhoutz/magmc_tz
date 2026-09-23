#pragma once

#include <cmath>

struct double3 {
  double x, y, z;

  double length() const { return std::sqrt(x * x + y * y + z * z); }
};

inline double3 cross(double3 const &a, double3 const &b) {
  return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}

inline double dot(double3 const &a, double3 const &b) { return a.x * b.x + a.y * b.y + a.z * b.z; }

inline double3 operator+(double3 const &a, double3 const &b) {
  return {a.x + b.x, a.y + b.y, a.z + b.z};
}

inline double3 operator*(double a, double3 const &b) { return {a * b.x, a * b.y, a * b.z}; }
inline double3 operator/(double3 const &a, double b) { return {a.x / b, a.y / b, a.z / b}; }

inline double3 to_unit(double3 const &a) { return a / a.length(); }
