//*********************************************************
//Title		:Surface.cpp
//Author	:Tanabe Yuta
//Date		:2019/01/26
//Copyright	:(C)2019 TanabeYuta
//*********************************************************


#include "pch.h"
#include "Parameter.h"
#include "Surface.h"


Surface::Surface(){}


Surface::~Surface(){}


Surface::Surface(Node* _pnode0, Node* _pnode1, Node* _pnode2, Element* _pparent, Element* _pneighbor){
	this->pnodes[0] = _pnode0;	this->pnodes[1] = _pnode1;	this->pnodes[2] = _pnode2;	
	this->pparent = _pparent;		this->pneighbor = _pneighbor;	this->IsActive = true;
}


bool Surface::operator==(const Surface &_surface){
	if ((this->pnodes[0] == _surface.pnodes[0] && this->pnodes[1] == _surface.pnodes[2] && this->pnodes[2] == _surface.pnodes[1])
		|| (this->pnodes[0] == _surface.pnodes[1] && this->pnodes[1] == _surface.pnodes[0] && this->pnodes[2] == _surface.pnodes[2])
		|| (this->pnodes[0] == _surface.pnodes[2] && this->pnodes[1] == _surface.pnodes[1] && this->pnodes[2] == _surface.pnodes[0])) {
		return true;
	}
	return false;
}


bool Surface::IsRayCross(Node _sp, Node _dir){
	Node v01 = *(this->pnodes[1]) - *(this->pnodes[0]);
	Node v02 = *(this->pnodes[2]) - *(this->pnodes[0]);
	Node v0g = _sp - *(this->pnodes[0]);

	double det = v01 ^ (v02 * _dir);
	if (det > EPS) {
		double u = (v0g ^ (v02 * _dir)) / det;
		if (-EPS < u && u < 1.0 + EPS) {
			double v = (v01 ^ (v0g * _dir)) / det;
			if (-EPS < v && u + v < 1.0 + EPS) {
				double t = (v01 ^ (v02 * v0g)) / det;
				if (t > -EPS && t < 1.0 - EPS) {
					return true;
				}
			}
		}
	}
	return false;
}
