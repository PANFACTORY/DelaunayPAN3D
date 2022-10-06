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

    Node<T> operator+(const Node<T>&);
    Node<T> operator-(const Node<T>&);
    Node<T> operator*(const Node<T>&);
    T operator^(const Node<T>&);
    Node<T> operator*(T);
    Node<T> operator/(T);
    bool operator==(const Node<T>&);

    T Norm();
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
Node<T> Node<T>::operator+(const Node<T>& _node) {
    return Node<T>(this->x + _node.x, this->y + _node.y, this->z + _node.z, -1,
                   -1);
}

template <class T>
Node<T> Node<T>::operator-(const Node<T>& _node) {
    return Node<T>(this->x - _node.x, this->y - _node.y, this->z - _node.z, -1,
                   -1);
}

template <class T>
Node<T> Node<T>::operator*(const Node<T>& _node) {
    return Node<T>(this->y * _node.z - this->z * _node.y,
                   this->z * _node.x - this->x * _node.z,
                   this->x * _node.y - this->y * _node.x, -1, -1);
}

template <class T>
T Node<T>::operator^(const Node<T>& _node) {
    return this->x * _node.x + this->y * _node.y + this->z * _node.z;
}

template <class T>
Node<T> Node<T>::operator*(T _a) {
    return Node<T>(this->x * _a, this->y * _a, this->z * _a, -1, -1);
}

template <class T>
Node<T> Node<T>::operator/(T _a) {
    return Node<T>(this->x / _a, this->y / _a, this->z / _a, -1, -1);
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

template <class T>
T Node<T>::Norm() {
    return sqrt(pow(this->x, 2.0) + pow(this->y, 2.0) + pow(this->z, 2.0));
}
}  // namespace DelaunayPAN3D
