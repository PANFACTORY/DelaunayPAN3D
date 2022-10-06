#pragma once
#include <array>
#include <cmath>

namespace DelaunayPAN3D {
template <class T>
class Vec3 {
   public:
    Vec3() {
        this->x[0] = T();
        this->x[1] = T();
        this->x[2] = T();
        this->type = -1;
    };
    ~Vec3() {}
    Vec3(T _x, T _y, T _z, int _type) {
        this->x[0] = _x;
        this->x[1] = _y;
        this->x[2] = _z;
        this->type = _type;
    }

    const Vec3<T> operator+(const Vec3<T>& _node) const {
        return Vec3<T>(this->x[0] + _node.x[0], this->x[1] + _node.x[1],
                       this->x[2] + _node.x[2], -1);
    }
    const Vec3<T> operator-(const Vec3<T>& _node) const {
        return Vec3<T>(this->x[0] - _node.x[0], this->x[1] - _node.x[1],
                       this->x[2] - _node.x[2], -1);
    }
    const Vec3<T> operator*(T _a) const {
        return Vec3<T>(this->x[0] * _a, this->x[1] * _a, this->x[2] * _a, -1);
    }
    const Vec3<T> operator/(T _a) const {
        return Vec3<T>(this->x[0] / _a, this->x[1] / _a, this->x[2] / _a, -1);
    }
    bool operator==(const Vec3<T>& _node) const {
        T EPS = 1e-15;
        if (fabs(this->x[0] - _node.x[0]) < EPS &&
            fabs(this->x[1] - _node.x[1]) < EPS &&
            fabs(this->x[2] - _node.x[2]) < EPS) {
            return true;
        }
        return false;
    }
    T& operator[](int i) { return this->x[i]; }
    T& operator()(int i) { return this->x[i]; }

    T dot(const Vec3<T>& node) const {
        return this->x[0] * node.x[0] + this->x[1] * node.x[1] +
               this->x[2] * node.x[2];
    }
    const Vec3<T> cross(const Vec3<T>& node) const {
        return Vec3<T>(this->x[1] * node.x[2] - this->x[2] * node.x[1],
                       this->x[2] * node.x[0] - this->x[0] * node.x[2],
                       this->x[0] * node.x[1] - this->x[1] * node.x[0], -1);
    }
    T norm() const { return sqrt(this->dot(*this)); }

    int type;

   private:
    std::array<T, 3> x;
};
}  // namespace DelaunayPAN3D
