//*********************************************************
//Title		:Element.cpp
//Author	:Tanabe Yuta
//Date		:2019/01/26
//Copyright	:(C)2019 TanabeYuta
//*********************************************************


#include "pch.h"
#include "Surface.h"
#include "Element.h"
#include "Parameter.h"


Element::Element(){}


Element::~Element(){
	for (auto& psurface : this->psurfaces) {
		delete psurface;
	}
}


Element::Element(Node *_pnode0, Node *_pnode1, Node *_pnode2, Node *_pnode3){
	this->IsActive = true;

	//----------各面を設定----------
	this->psurfaces[0] = new Surface(_pnode1, _pnode3, _pnode2, this, nullptr);
	this->psurfaces[1] = new Surface(_pnode0, _pnode2, _pnode3, this, nullptr);
	this->psurfaces[2] = new Surface(_pnode0, _pnode3, _pnode1, this, nullptr);
	this->psurfaces[3] = new Surface(_pnode0, _pnode1, _pnode2, this, nullptr);

	//----------外接球の中心座標と半径を計算----------
	Node v0 = *_pnode1 - *_pnode0;
	Node v1 = *_pnode2 - *_pnode0;
	Node v2 = *_pnode3 - *_pnode0;

	Node ABC = Node(0.5*((*_pnode1 ^ *_pnode1) - (*_pnode0 ^ *_pnode0)), 0.5*((*_pnode2 ^ *_pnode2) - (*_pnode0 ^ *_pnode0)), 0.5*((*_pnode3 ^ *_pnode3) - (*_pnode0 ^ *_pnode0)), -1);

	double detP = v0 ^ (v1 * v2);
	Node P0 = v1 * v2;
	Node P1 = v2 * v0;
	Node P2 = v0 * v1;

	this->scenter = Node((ABC.x * P0.x + ABC.y * P1.x + ABC.z * P2.x) / detP, (ABC.x * P0.y + ABC.y * P1.y + ABC.z * P2.y) / detP, (ABC.x * P0.z + ABC.y * P1.z + ABC.z * P2.z) / detP, -1);

	this->sround = (this->scenter - *_pnode0).Size();

	//----------重心座標の計算----------
	this->gcenter = (*_pnode0 + *_pnode1 + *_pnode2 + *_pnode3) / 4.0;
}


Element* Element::GetLocateId(Node *_pnode){
	for (auto surface : this->psurfaces) {
		if (surface->IsRayCross(this->gcenter, this->gcenter - *_pnode) == true) {
			return surface->pneighbor;
		}
	}
	return this;
}


bool Element::IsInSphere(Node *_pnode){
	if (sround > (this->scenter - *_pnode).Size() - EPS) {
		return true;
	}
	return false;
}


Surface* Element::GetAdjacentSurface(Element* _pelement){
	for (auto& psurface : this->psurfaces) {
		if (psurface->pneighbor == _pelement) {
			return psurface;
		}
	}
	return nullptr;
}