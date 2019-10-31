//*********************************************************
//Title		:Delaunay.h
//Author	:Tanabe Yuta
//Date		:2019/01/08
//Copyright	:(C)2019 TanabeYuta
//*********************************************************


#pragma once
#define _USE_MATH_DEFINES
#include <math.h>
#include <iostream>
#include <vector>


#include "Parameter.h"
#include "Node.h"
#include "Surface.h"
#include "Element.h"


namespace DelaunayPAN3DV3{
	template<class T>
	void MakeSupertetrahedran(std::vector<Node<T>*> &_nlist, std::vector<Element<T>*> &_elist) {
		std::cout << "Make supertetraedron\n";
				
		//----------�ߓ_�̑��ݔ͈͂����߂�----------
		T xmax = 0.0, xmin = 0.0, ymax = 0.0, ymin = 0.0, zmax = 0.0, zmin = 0.0;
		for (auto pnode : _nlist) {
			if (pnode->x > xmax) { xmax = pnode->x; }
			if (pnode->x < xmin) { xmin = pnode->x; }
			if (pnode->y > ymax) { ymax = pnode->y; }
			if (pnode->y < ymin) { ymin = pnode->y; }
			if (pnode->z > zmax) { zmax = pnode->z; }
			if (pnode->z < zmin) { zmin = pnode->z; }
		}

		//----------���W�̐��K��----------
		T xrange = 0.5*(xmax - xmin), yrange = 0.5*(ymax - ymin), zrange = 0.5*(zmax - zmin);
		T dmax = xrange;
		if (dmax < yrange) {	dmax = yrange;	}
		if (dmax < zrange) {	dmax = zrange;	}
		for (auto& pnode : _nlist) {
			pnode->x = (pnode->x - xmin) / dmax + 0.5*(ALPHA - 1.0)*xrange / dmax;
			pnode->y = (pnode->y - ymin) / dmax + 0.5*(ALPHA - 1.0)*yrange / dmax;
			pnode->z = (pnode->z - zmin) / dmax + 0.5*(ALPHA - 1.0)*zrange / dmax;
		}

		//----------���z�l�ʑ̐ߓ_�̐���----------
		Node<T>* nst0 = new Node<T>(0.0, 0.0, 0.0, -1, _nlist.size());
		_nlist.push_back(nst0);
		Node<T>* nst1 = new Node<T>(ALPHA*xrange / dmax, 0.0, 0.0, -1, _nlist.size());
		_nlist.push_back(nst1);
		Node<T>* nst2 = new Node<T>(ALPHA*xrange / dmax, ALPHA*yrange / dmax, 0.0, -1, _nlist.size());
		_nlist.push_back(nst2);
		Node<T>* nst3 = new Node<T>(0.0, ALPHA*yrange / dmax, 0.0, -1, _nlist.size());
		_nlist.push_back(nst3);
		Node<T>* nst4 = new Node<T>(0.0, 0.0, ALPHA*zrange / dmax, -1, _nlist.size());
		_nlist.push_back(nst4);
		Node<T>* nst5 = new Node<T>(ALPHA*xrange / dmax, 0.0, ALPHA*zrange / dmax, -1, _nlist.size());
		_nlist.push_back(nst5);
		Node<T>* nst6 = new Node<T>(ALPHA*xrange / dmax, ALPHA*yrange / dmax, ALPHA*zrange / dmax, -1, _nlist.size());
		_nlist.push_back(nst6);
		Node<T>* nst7 = new Node<T>(0.0, ALPHA*yrange / dmax, ALPHA*zrange / dmax, -1, _nlist.size());
		_nlist.push_back(nst7);
		
		//----------���z�l�ʑ̌Q�̐���----------
		_elist.push_back(new Element<T>(nst1, nst3, nst0, nst7));
		_elist.push_back(new Element<T>(nst2, nst1, nst6, nst7));
		_elist.push_back(new Element<T>(nst2, nst3, nst1, nst7));
		_elist.push_back(new Element<T>(nst1, nst5, nst6, nst7));
		_elist.push_back(new Element<T>(nst1, nst0, nst5, nst7));
		_elist.push_back(new Element<T>(nst4, nst5, nst0, nst7));

		//----------�v�f���m�̗אڊ֌W���v�Z----------
		for (auto& pelement : _elist) {
			for (auto& psurface : pelement->psurfaces) {
				if (psurface->pneighbor == nullptr) {
					for (auto& pelement2 : _elist) {
						for (auto& psurface2 : pelement2->psurfaces) {
							if (*psurface == *psurface2) {
								psurface->pneighbor = pelement2;
								psurface2->pneighbor = pelement;
								break;
							}
						}
					}
				}
			}
		}
	}


	template<class T>
	void MeshLocal(Node<T> *_node, Element<T>* _pethis, std::vector<Element<T>*>& _elist) {
		std::vector<Element<T>*> stack, substack;
		std::vector<Surface<T>*> sstack;

		//----------�ǉ������_���O�ڋ����Ɏ��v�f���W�߂�----------
		substack.push_back(_pethis);
		while (substack.size()) {
			Element<T>* pend = *(substack.end() - 1);			
			substack.pop_back();

			if (pend->IsActive) {
				stack.push_back(pend);
				pend->IsActive = false;

				for (auto& psurface : pend->psurfaces) {
					Element<T>* pneighbor = psurface->pneighbor;
					if (pneighbor != nullptr && pneighbor->IsInSphere(_node)) {
						substack.push_back(pneighbor);
					}
					else {
						sstack.push_back(psurface);
					}
				}
			}
		}

		//----------���ʑ̂̉����𖄂߂�----------
		bool is_anysurface_invalid = true;
		while (is_anysurface_invalid) {
			is_anysurface_invalid = false;

			for (int i = 0; i < sstack.size(); i++) {
				if (sstack[i]->IsActive) {
					Element<T> D = Element<T>(sstack[i]->pnodes[0], sstack[i]->pnodes[1], sstack[i]->pnodes[2], _node);
					
					//----------�s�ǂȖʂ�����ꍇ----------
					if (D.volume < EPS) {
						Element<T>* peadd = sstack[i]->pneighbor;		

						//----------�V�����ǉ��ł���v�f������ꍇ----------
						if (peadd != nullptr) {
							if (peadd->IsActive) {
								is_anysurface_invalid = true;
								peadd->IsActive = false;
								stack.push_back(peadd);
								
								//----------�ǉ������v�f�̊e�ʂɂ��ċ��L�ʂ�False��----------
								for (auto& psurface : peadd->psurfaces) {
									Element<T>* pneighbor = psurface->pneighbor;
									if (pneighbor != nullptr && !pneighbor->IsActive) {
										pneighbor->GetAdjacentSurface(peadd)->IsActive = false;
									}
									else {
										sstack.push_back(psurface);
									}
								}
								break;
							}
						}

						//----------�V�����ǉ��ł���v�f�������ꍇ----------
						else if (fabs(D.volume) < EPS) {
							sstack[i]->IsActive = false;
						}
					}
				}
			}
		}

		//----------�V�����v�f�𐶐�----------
		std::vector<Element<T>*> penew;			
		for (auto& psurface : sstack) {
			if (psurface->IsActive) {
				Element<T>* tmp = new Element<T>(psurface->pnodes[0], psurface->pnodes[1], psurface->pnodes[2], _node);
				tmp->psurfaces[3]->pneighbor = psurface->pneighbor;
				if (psurface->pneighbor != nullptr) {
					psurface->pneighbor->GetAdjacentSurface(psurface->pparent)->pneighbor = tmp;
				}
				penew.push_back(tmp);
				_elist.push_back(tmp);
			}
		}

		//----------�V�����������ꂽ�v�f���m�̗אڊ֌W���v�Z----------
		for (auto& pelement : penew) {
			for (auto& psurface : pelement->psurfaces) {
				OUT:
				if (psurface->pneighbor == nullptr) {
					for (auto& pelement2 : penew) {
						for (auto& psurface2 : pelement2->psurfaces) {
							if (*psurface == *psurface2) {

								//----------�G���v�f�����o���ꂽ�Ƃ�----------
								if (psurface2->pneighbor != nullptr) {
									std::cout << "!!";
								}

								psurface->pneighbor = pelement2;
								psurface2->pneighbor = pelement;
								goto OUT;
							}
						}
					}
				}
			}
		}

		//----------stack���̌Â��v�f���폜----------
		for (auto pelement = _elist.begin(); pelement != _elist.end(); ) {
			if (!(*pelement)->IsActive) {
				delete *pelement;
				pelement = _elist.erase(pelement);
			}
			else {
				++pelement;
			}
		}
	}


	template<class T>
	void MakeRoughMesh(std::vector<Node<T>*> _nlist, std::vector<Element<T>*> &_elist) {
		std::cout << "Make rough mesh\n";

		Element<T>* pethis = _elist[0];									
		for (auto& pnode : _nlist) {
			if (pnode->type != -1) {
				int count = 0;
				while (1) {
					Element<T>* penext = pethis->GetLocateId(pnode);				
					//----------�v�f���ɓ_������Ƃ�----------
					if (penext == pethis) {
						//std::cout << "\tat\t" << pethis << "\n";
						MeshLocal(pnode, pethis, _elist);
						pethis = *(_elist.end() - 1);
						break;
					}
					//----------�Ȃ��Ƃ�----------
					else {
						pethis = penext;
					}
				}
			}
		}
	}


	template<class T>
	void DeleteSupertetrahedran(std::vector<Element<T>*> &_elist) {
		std::cout << "Delete supertetraedron\n";
		
		for (auto pelement = _elist.begin(); pelement != _elist.end(); ) {
			if ((*pelement)->pnodes[0]->type == -1 
				|| (*pelement)->pnodes[1]->type == -1 
				|| (*pelement)->pnodes[2]->type == -1 
				|| (*pelement)->pnodes[3]->type == -1) {
				for (auto& psurface : (*pelement)->psurfaces) {
					if (psurface->pneighbor != nullptr) {
						psurface->pneighbor->GetAdjacentSurface(psurface->pparent)->pneighbor = nullptr;
					}
				}				
				delete *pelement;
				pelement = _elist.erase(pelement);
			}
			else {
				++pelement;
			}
		}
	}


	template<class T>
	void DeleteCreviceElement(std::vector<Element<T>*> &_elist) {
		std::cout << "Delete Crevice Element\n";

		for (auto pelement = _elist.begin(); pelement != _elist.end(); ) {
			if ((*pelement)->pnodes[0]->type == (*pelement)->pnodes[1]->type
				&& (*pelement)->pnodes[1]->type == (*pelement)->pnodes[2]->type
				&& (*pelement)->pnodes[2]->type == (*pelement)->pnodes[3]->type) {
				for (auto& psurface : (*pelement)->psurfaces) {
					if (psurface->pneighbor != nullptr) {
						psurface->pneighbor->GetAdjacentSurface(psurface->pparent)->pneighbor = nullptr;
					}
				}
				delete *pelement;
				pelement = _elist.erase(pelement);
			}
			else {
				++pelement;
			}
		}
	}


	template<class T>
	void MakeFineMesh(std::vector<Node<T>*> &_nlist, std::vector<Element<T>*> &_elist) {
		std::cout << "Make fine mesh\n";
		
		for (int i = 0; i < ADDNODE; i++) {
			//----------�Œ��̕ӂ�T��----------
			T edgelengthmax = 0.0;
			Element<T>* pethis = nullptr;
			Node<T>* pnode0 = nullptr;
			Node<T>* pnode1 = nullptr;
					   
			for (auto pelement : _elist) {
				for (int j = 0; j < 3; j++) {
					for (int k = j + 1; k < 3; k++) {
						T edgelength = (*pelement->pnodes[k] - *pelement->pnodes[j]).Norm();
						if (edgelength > edgelengthmax) {
							edgelengthmax = edgelength;
							pethis = pelement;
							pnode0 = pelement->pnodes[j];
							pnode1 = pelement->pnodes[k];
						}
					}
				}
			}
					   
			//----------�Œ��̕ӂ̒��_��ߓ_�ɒǉ�----------
			Node<T>* tmp = new Node<T>((*pnode0 + *pnode1) / 2.0);
			tmp->type = 2;
			tmp->id = _nlist.size();
			_nlist.push_back(tmp);
			MeshLocal(tmp, pethis, _elist);
		}
	}
}