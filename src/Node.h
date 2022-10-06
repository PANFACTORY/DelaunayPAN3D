//*********************************************************
// Title		:src/cpp/Node.h
// Author	:Tanabe Yuta
// Date		:2019/01/26
// Copyright	:(C)2019 TanabeYuta
//*********************************************************

#pragma once
#include <array>
#include <cmath>

namespace DelaunayPAN3D {
template <class T>
class Node {
   public:
    Node() {}
    ~Node() {}
    Node(T _x, T _y, T _z, int _type, int _id) {
        this->x[0] = _x;
        this->x[1] = _y;
        this->x[2] = _z;
        this->type = _type;
        this->id = _id;
    }

    const Node<T> operator+(const Node<T>& _node) const {
        return Node<T>(this->x[0] + _node.x[0], this->x[1] + _node.x[1],
                       this->x[2] + _node.x[2], -1, -1);
    }
    const Node<T> operator-(const Node<T>& _node) const {
        return Node<T>(this->x[0] - _node.x[0], this->x[1] - _node.x[1],
                       this->x[2] - _node.x[2], -1, -1);
    }
    const Node<T> operator*(T _a) const {
        return Node<T>(this->x[0] * _a, this->x[1] * _a, this->x[2] * _a, -1,
                       -1);
    }
    const Node<T> operator/(T _a) const {
        return Node<T>(this->x[0] / _a, this->x[1] / _a, this->x[2] / _a, -1,
                       -1);
    }
    bool operator==(const Node<T>& _node) const {
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

    T dot(const Node<T>& node) const {
        return this->x[0] * node.x[0] + this->x[1] * node.x[1] +
               this->x[2] * node.x[2];
    }
    const Node<T> cross(const Node<T>& node) const {
        return Node<T>(this->x[1] * node.x[2] - this->x[2] * node.x[1],
                       this->x[2] * node.x[0] - this->x[0] * node.x[2],
                       this->x[0] * node.x[1] - this->x[1] * node.x[0], -1, -1);
    }
    T norm() const { return sqrt(this->dot(*this)); }

    int type;
    int id;

   private:
    std::array<T, 3> x;
};
}  // namespace DelaunayPAN3D
