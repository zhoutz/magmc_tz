#pragma once

#include <cmath>

struct double3 {
  double x, y, z;
};

inline double3 cross(double3 const &a, double3 const &b) {
  return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}

inline double3 to_unit(double3 const &a) {
  double norm = std::sqrt(a.x * a.x + a.y * a.y + a.z * a.z);
  return {a.x / norm, a.y / norm, a.z / norm};
}

inline double dot(double3 const &a, double3 const &b) { return a.x * b.x + a.y * b.y + a.z * b.z; }

inline double3 operator+(double3 const &a, double3 const &b) {
  return {a.x + b.x, a.y + b.y, a.z + b.z};
}

inline double3 operator*(double a, double3 const &b) { return {a * b.x, a * b.y, a * b.z}; }
