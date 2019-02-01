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


namespace Delaunay3D {
	//**********仮想四面体の生成**********
	void MakeSupertetrahedran(std::vector<Node*> &_nlist, std::vector<Element*> &_elist) {
		std::cout << "Make supertetraedron\n";
				
		//----------節点の存在範囲を求める----------
		double xmax = 0.0, xmin = 0.0, ymax = 0.0, ymin = 0.0, zmax = 0.0, zmin = 0.0;
		for (auto pnode : _nlist) {
			if (pnode->x > xmax) { xmax = pnode->x; }
			if (pnode->x < xmin) { xmin = pnode->x; }
			if (pnode->y > ymax) { ymax = pnode->y; }
			if (pnode->y < ymin) { ymin = pnode->y; }
			if (pnode->z > zmax) { zmax = pnode->z; }
			if (pnode->z < zmin) { zmin = pnode->z; }
		}

		//----------座標の正規化----------
		double xrange = 0.5*(xmax - xmin), yrange = 0.5*(ymax - ymin), zrange = 0.5*(zmax - zmin);
		double dmax = xrange;
		if (dmax < yrange) {	dmax = yrange;	}
		if (dmax < zrange) {	dmax = zrange;	}
		for (auto& pnode : _nlist) {
			pnode->x = (pnode->x - xmin) / dmax + 0.5*(ALPHA - 1.0)*xrange / dmax;
			pnode->y = (pnode->y - ymin) / dmax + 0.5*(ALPHA - 1.0)*yrange / dmax;
			pnode->z = (pnode->z - zmin) / dmax + 0.5*(ALPHA - 1.0)*zrange / dmax;
		}

		//----------仮想四面体節点の生成----------
		Node* nst0 = new Node(0.0, 0.0, 0.0, -1, _nlist.size());
		_nlist.push_back(nst0);
		Node* nst1 = new Node(ALPHA*xrange / dmax, 0.0, 0.0, -1, _nlist.size());
		_nlist.push_back(nst1);
		Node* nst2 = new Node(ALPHA*xrange / dmax, ALPHA*yrange / dmax, 0.0, -1, _nlist.size());
		_nlist.push_back(nst2);
		Node* nst3 = new Node(0.0, ALPHA*yrange / dmax, 0.0, -1, _nlist.size());
		_nlist.push_back(nst3);
		Node* nst4 = new Node(0.0, 0.0, ALPHA*zrange / dmax, -1, _nlist.size());
		_nlist.push_back(nst4);
		Node* nst5 = new Node(ALPHA*xrange / dmax, 0.0, ALPHA*zrange / dmax, -1, _nlist.size());
		_nlist.push_back(nst5);
		Node* nst6 = new Node(ALPHA*xrange / dmax, ALPHA*yrange / dmax, ALPHA*zrange / dmax, -1, _nlist.size());
		_nlist.push_back(nst6);
		Node* nst7 = new Node(0.0, ALPHA*yrange / dmax, ALPHA*zrange / dmax, -1, _nlist.size());
		_nlist.push_back(nst7);
		
		//----------仮想四面体群の生成----------
		_elist.push_back(new Element(nst1, nst3, nst0, nst7));
		_elist.push_back(new Element(nst2, nst1, nst6, nst7));
		_elist.push_back(new Element(nst2, nst3, nst1, nst7));
		_elist.push_back(new Element(nst1, nst5, nst6, nst7));
		_elist.push_back(new Element(nst1, nst0, nst5, nst7));
		_elist.push_back(new Element(nst4, nst5, nst0, nst7));

		//----------要素同士の隣接関係を計算----------
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


	//**********局所Delaunay分割**********
	void MeshLocal(Node *_node, Element *_pethis, std::vector<Element*> &_elist) {
		std::vector<Element*> stack, substack;
		std::vector<Surface*> sstack;

		//----------追加した点を外接球内に持つ要素を集める----------
		substack.push_back(_pethis);
		while (substack.size()) {
			Element* pend = *(substack.end() - 1);				//サブスタック末尾の要素を指すポインタ
			substack.pop_back();

			if (pend->IsActive) {
				stack.push_back(pend);
				pend->IsActive = false;

				for (auto& psurface : pend->psurfaces) {
					Element* pneighbor = psurface->pneighbor;
					if (pneighbor != nullptr && pneighbor->IsInSphere(_node)) {
						substack.push_back(pneighbor);
					}
					else {
						sstack.push_back(psurface);
					}
				}
			}
		}

		//----------多面体の凹部を埋める----------
		//「新しく追加できる要素が無い場合」の部分を経由したときバグが発生
		bool is_anysurface_invalid = true;
		while (is_anysurface_invalid) {
			is_anysurface_invalid = false;

			for (int i = 0; i < sstack.size(); i++) {
				if (sstack[i]->IsActive) {
					Element D = Element(sstack[i]->pnodes[0], sstack[i]->pnodes[1], sstack[i]->pnodes[2], _node);
					
					//----------不良な面がある場合----------
					if (D.volume < EPS) {
						Element* peadd = sstack[i]->pneighbor;			//不良な面に隣接する要素を指すポインタ

						//----------新しく追加できる要素がある場合----------
						if (peadd != nullptr) {
							if (peadd->IsActive) {
								is_anysurface_invalid = true;
								peadd->IsActive = false;
								stack.push_back(peadd);
								
								//----------追加した要素の各面について共有面をFalseに----------
								for (auto& psurface : peadd->psurfaces) {
									Element* pneighbor = psurface->pneighbor;
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

						//----------新しく追加できる要素が無い場合----------
						else if (fabs(D.volume) < EPS) {
							sstack[i]->IsActive = false;
						}
					}
				}
			}
		}

		//----------新しい要素を生成----------
		std::vector<Element*> penew;				//新しく生成される要素を指すポインタのスタック
		for (auto& psurface : sstack) {
			if (psurface->IsActive) {
				Element* tmp = new Element(psurface->pnodes[0], psurface->pnodes[1], psurface->pnodes[2], _node);
				tmp->psurfaces[3]->pneighbor = psurface->pneighbor;
				if (psurface->pneighbor != nullptr) {
					psurface->pneighbor->GetAdjacentSurface(psurface->pparent)->pneighbor = tmp;
				}
				penew.push_back(tmp);
				_elist.push_back(tmp);
			}
		}

		//----------新しく生成された要素同士の隣接関係を計算----------
		//扁平要素を検出したときに自己修正できるのがベスト
		for (auto& pelement : penew) {
			for (auto& psurface : pelement->psurfaces) {
				OUT:
				if (psurface->pneighbor == nullptr) {
					for (auto& pelement2 : penew) {
						for (auto& psurface2 : pelement2->psurfaces) {
							if (*psurface == *psurface2) {

								//----------扁平要素が検出されたとき----------
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

		//----------stack内の古い要素を削除----------
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


	//**********粗いDelaunay分割**********
	void MakeRoughMesh(std::vector<Node*> _nlist, std::vector<Element*> &_elist) {
		std::cout << "Make rough mesh\n";

		Element* pethis = _elist[0];										//現在調べている要素を指すポインタ
		for (auto& pnode : _nlist) {
			if (pnode->type != -1) {
				int count = 0;
				while (1) {
					Element* penext = pethis->GetLocateId(pnode);				//次に調べる要素を指すポインタ
					//----------要素内に点があるとき----------
					if (penext == pethis) {
						//std::cout << "\tat\t" << pethis << "\n";
						MeshLocal(pnode, pethis, _elist);
						pethis = *(_elist.end() - 1);
						break;
					}
					//----------ないとき----------
					else {
						pethis = penext;
					}
				}
			}
		}
	}


	//**********仮想四面体の削除**********
	void DeleteSupertetrahedran(std::vector<Element*> &_elist) {
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


	//**********凹部要素の削除**********
	void DeleteCreviceElement(std::vector<Element*> &_elist) {
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


	//**********細かいDelaunay分割**********
	void MakeFineMesh(std::vector<Node*> &_nlist, std::vector<Element*> &_elist) {
		std::cout << "Make fine mesh\n";
		
		for (int i = 0; i < ADDNODE; i++) {
			//----------最長の辺を探索----------
			double edgelengthmax = 0.0;
			Element* pethis = nullptr;
			Node* pnode0 = nullptr;
			Node* pnode1 = nullptr;
					   
			for (auto pelement : _elist) {
				for (int j = 0; j < 3; j++) {
					for (int k = j + 1; k < 3; k++) {
						double edgelength = (*pelement->pnodes[k] - *pelement->pnodes[j]).Size();
						if (edgelength > edgelengthmax) {
							edgelengthmax = edgelength;
							pethis = pelement;
							pnode0 = pelement->pnodes[j];
							pnode1 = pelement->pnodes[k];
						}
					}
				}
			}
					   
			//----------最長の辺の中点を節点に追加----------
			Node* tmp = new Node((*pnode0 + *pnode1) / 2.0);
			tmp->type = 2;
			tmp->id = _nlist.size();
			_nlist.push_back(tmp);
			MeshLocal(tmp, pethis, _elist);
		}
	}
}