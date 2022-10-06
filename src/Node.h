//*********************************************************
// Title		:src/cpp/Node.h
// Author	:Tanabe Yuta
// Date		:2019/01/26
// Copyright	:(C)2019 TanabeYuta
//*********************************************************

#pragma once
#include <cmath>

namespace DelaunayPAN3D {
template <class T>
class Node {
   public:
    Node();
    ~Node();
    Node(T, T, T, int, int);

    T x, y, z;
    int type;
    int id;

    const Node<T> operator+(const Node<T>& _node) const {
        return Node<T>(this->x + _node.x, this->y + _node.y, this->z + _node.z,
                       -1, -1);
    }
    const Node<T> operator-(const Node<T>& _node) const {
        return Node<T>(this->x - _node.x, this->y - _node.y, this->z - _node.z,
                       -1, -1);
    }
    const Node<T> operator*(T _a) const {
        return Node<T>(this->x * _a, this->y * _a, this->z * _a, -1, -1);
    }
    const Node<T> operator/(T _a) const {
        return Node<T>(this->x / _a, this->y / _a, this->z / _a, -1, -1);
    }
    bool operator==(const Node<T>&);

    T dot(const Node<T>& node) const {
        return this->x * node.x + this->y * node.y + this->z * node.z;
    }
    const Node<T> cross(const Node<T>& node) const {
        return Node<T>(this->y * node.z - this->z * node.y,
                       this->z * node.x - this->x * node.z,
                       this->x * node.y - this->y * node.x, -1, -1);
    }
    T norm() const { return sqrt(this->dot(*this)); }
};

template <class T>
Node<T>::Node() {}

template <class T>
Node<T>::~Node() {}

template <class T>
Node<T>::Node(T _x, T _y, T _z, int _type, int _id) {
    this->x = _x;
    this->y = _y;
    this->z = _z;
    this->type = _type;
    this->id = _id;
}

template <class T>
bool Node<T>::operator==(const Node<T>& _node) {
    T EPS = 1e-15;
    if (fabs(this->x - _node.x) < EPS && fabs(this->y - _node.y) < EPS &&
        fabs(this->z - _node.z) < EPS) {
        return true;
    }
    return false;
}
}  // namespace DelaunayPAN3D
