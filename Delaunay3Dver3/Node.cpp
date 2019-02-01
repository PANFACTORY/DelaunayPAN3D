//*********************************************************
//Title		:Node.cpp
//Author	:Tanabe Yuta
//Date		:2019/01/26
//Copyright	:(C)2019 TanabeYuta
//*********************************************************


#include "pch.h"
#include <cmath>

#include "Node.h"
#include "Parameter.h"


Node::Node() {}


Node::~Node() {}


Node::Node(double _x, double _y, double _z, int _type, int _id) {
	this->x = _x;	this->y = _y;	this->z = _z;	this->type = _type;		this->id = _id;
}


Node Node::operator+(const Node &_node) {
	return Node(this->x + _node.x, this->y + _node.y, this->z + _node.z, -1, -1);
}


Node Node::operator-(const Node &_node) {
	return Node(this->x - _node.x, this->y - _node.y, this->z - _node.z, -1, -1);
}


Node Node::operator*(const Node &_node) {
	return Node(this->y*_node.z - this->z*_node.y, this->z*_node.x - this->x*_node.z, this->x*_node.y - this->y*_node.x, -1, -1);
}


double Node::operator^(const Node &_node) {
	return this->x*_node.x + this->y*_node.y + this->z*_node.z;
}


Node Node::operator*(double _a) {
	return Node(this->x *_a, this->y *_a, this->z *_a, -1, -1);
}


Node Node::operator/(double _a) {
	return Node(this->x / _a, this->y / _a, this->z / _a, -1, -1);
}


bool Node::operator==(const Node &_node) {
	if (fabs(this->x - _node.x) < EPS && fabs(this->y - _node.y) < EPS && fabs(this->z - _node.z) < EPS){
		return true;
	}
	return false;
}


double Node::Size() {
	return sqrt(pow(this->x, 2.0) + pow(this->y, 2.0) + pow(this->z, 2.0));
}