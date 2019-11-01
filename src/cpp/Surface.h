//*********************************************************
//Title		:src/cpp/Surface.h
//Author	:Tanabe Yuta
//Date		:2019/01/26
//Copyright	:(C)2019 TanabeYuta
//*********************************************************


#pragma once
#include <array>


#include "Parameter.h"
#include "Node.h"


namespace DelaunayPAN3DV3{
	template<class T>
	class Element;


	template<class T>
	class Surface
	{
	public:
		Surface();
		~Surface();
		Surface(Node<T>*, Node<T>*, Node<T>*, Element<T>*, Element<T>*);


		Element<T>* pneighbor;					
		bool IsActive;					
		std::array<Node<T>*, 3> pnodes;	
		Element<T>* pparent;				
		
		
		bool operator==(const Surface<T>&);	


		bool IsRayCross(Node<T>, Node<T>);	
	};


	template<class T>
	Surface<T>::Surface(){}


	template<class T>
	Surface<T>::~Surface(){}


	template<class T>
	Surface<T>::Surface(Node<T>* _pnode0, Node<T>* _pnode1, Node<T>* _pnode2, Element<T>* _pparent, Element<T>* _pneighbor){
		this->pnodes[0] = _pnode0;	
		this->pnodes[1] = _pnode1;	
		this->pnodes[2] = _pnode2;	
		this->pparent = _pparent;		
		this->pneighbor = _pneighbor;	
		this->IsActive = true;
	}


	template<class T>
	bool Surface<T>::operator==(const Surface<T>& _surface){
		if ((this->pnodes[0] == _surface.pnodes[0] && this->pnodes[1] == _surface.pnodes[2] && this->pnodes[2] == _surface.pnodes[1])
			|| (this->pnodes[0] == _surface.pnodes[1] && this->pnodes[1] == _surface.pnodes[0] && this->pnodes[2] == _surface.pnodes[2])
			|| (this->pnodes[0] == _surface.pnodes[2] && this->pnodes[1] == _surface.pnodes[1] && this->pnodes[2] == _surface.pnodes[0])) {
			return true;
		}
		return false;
	}


	template<class T>
	bool Surface<T>::IsRayCross(Node<T> _sp, Node<T> _dir){
		Node<T> v01 = *(this->pnodes[1]) - *(this->pnodes[0]);
		Node<T> v02 = *(this->pnodes[2]) - *(this->pnodes[0]);
		Node<T> v0g = _sp - *(this->pnodes[0]);

		T det = v01 ^ (v02 * _dir);
		if (det > EPS) {
			T u = (v0g ^ (v02 * _dir)) / det;
			if (-EPS < u && u < 1.0 + EPS) {
				T v = (v01 ^ (v0g * _dir)) / det;
				if (-EPS < v && u + v < 1.0 + EPS) {
					T t = (v01 ^ (v02 * v0g)) / det;
					if (t > -EPS && t < 1.0 - EPS) {
						return true;
					}
				}
			}
		}
		return false;
	}
}