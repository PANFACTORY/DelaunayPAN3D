//*********************************************************
//Title		:src/cpp/Element.h
//Author	:Tanabe Yuta
//Date		:2019/01/26
//Copyright	:(C)2019 TanabeYuta
//*********************************************************


#pragma once
#include <fenv.h>
#include <cmath>
#include <array>


#include "Node.h"
#include "Surface.h"
#include "Element.h"
#include "Parameter.h"


namespace DelaunayPAN3D {
	template<class T>
	class Surface;


	template<class T>
	class Element
	{
public:
		Element();
		~Element();
		Element(Node<T>*, Node<T>*, Node<T>*, Node<T>*);


		bool IsActive;
		std::array<Surface<T>*, 4> psurfaces;
		std::array<Node<T>*, 4> pnodes;		
		Node<T> scenter;					
		T sround;						
		Node<T> gcenter;							
		T volume;							
		T aspect;						


		Element<T>* GetLocateId(Node<T>*);			
		bool IsInSphere(Node<T>*);					
		Surface<T>* GetAdjacentSurface(Element<T>*);	
	};


	template<class T>
	Element<T>::Element(){}


	template<class T>
	Element<T>::~Element(){
		for (auto& psurface : this->psurfaces) {
			delete psurface;
		}
	}


	template<class T>
	Element<T>::Element(Node<T>* _pnode0, Node<T>* _pnode1, Node<T>* _pnode2, Node<T>* _pnode3){
		this->IsActive = true;

		//----------Set nodes----------
		this->pnodes[0] = _pnode0;	
		this->pnodes[1] = _pnode1;	
		this->pnodes[2] = _pnode2;	
		this->pnodes[3] = _pnode3;

		//----------Set surfaces----------
		this->psurfaces[0] = new Surface<T>(_pnode1, _pnode3, _pnode2, this, nullptr);
		this->psurfaces[1] = new Surface<T>(_pnode0, _pnode2, _pnode3, this, nullptr);
		this->psurfaces[2] = new Surface<T>(_pnode0, _pnode3, _pnode1, this, nullptr);
		this->psurfaces[3] = new Surface<T>(_pnode0, _pnode1, _pnode2, this, nullptr);

		//----------Get center and radius of external sphere----------
		Node<T> v0 = *_pnode1 - *_pnode0;
		Node<T> v1 = *_pnode2 - *_pnode0;
		Node<T> v2 = *_pnode3 - *_pnode0;

		Node<T> ABC = Node<T>(0.5*((*_pnode1 ^ *_pnode1) - (*_pnode0 ^ *_pnode0)), 0.5*((*_pnode2 ^ *_pnode2) - (*_pnode0 ^ *_pnode0)), 0.5*((*_pnode3 ^ *_pnode3) - (*_pnode0 ^ *_pnode0)), -1, -1);

		T detP = v0 ^ (v1 * v2);
		Node<T> P0 = v1 * v2;
		Node<T> P1 = v2 * v0;
		Node<T> P2 = v0 * v1;

		this->scenter = Node<T>((ABC.x * P0.x + ABC.y * P1.x + ABC.z * P2.x) / detP, (ABC.x * P0.y + ABC.y * P1.y + ABC.z * P2.y) / detP, (ABC.x * P0.z + ABC.y * P1.z + ABC.z * P2.z) / detP, -1, -1);
		this->sround = (this->scenter - *_pnode0).Norm();

		//----------Get center of gravity----------
		this->gcenter = (*_pnode0 + *_pnode1 + *_pnode2 + *_pnode3) / 4.0;

		//----------Get volume----------
		this->volume = ((*_pnode1 - *_pnode0) * (*_pnode2 - *_pnode0)) ^ (*_pnode3 - *_pnode0);
		
		//----------Get aspect ratio----------
		this->aspect = this->volume / pow(this->sround, 3.0) / ARTETRAHEDRON;
		if (fetestexcept(FE_DIVBYZERO)) {
			feclearexcept(FE_ALL_EXCEPT);
			this->aspect = T();
		}
	}


	template<class T>
	Element<T>* Element<T>::GetLocateId(Node<T>* _pnode){
		for (auto surface : this->psurfaces) {
			if (surface->IsRayCross(this->gcenter, this->gcenter - *_pnode)) {
				return surface->pneighbor;
			}
		}
		return this;
	}


	template<class T>
	bool Element<T>::IsInSphere(Node<T>* _pnode){
		if (sround + EPS > (this->scenter - *_pnode).Norm()) {
			return true;
		}
		return false;
	}


	template<class T>
	Surface<T>* Element<T>::GetAdjacentSurface(Element<T>* _pelement){
		for (auto& psurface : this->psurfaces) {
			if (psurface->pneighbor == _pelement) {
				return psurface;
			}
		}
		return nullptr;
	}
}