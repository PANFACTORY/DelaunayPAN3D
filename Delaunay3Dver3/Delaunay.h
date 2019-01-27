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

		//----------節点のうち最遠方点を探索----------
		double rmax = 0.0;
		for (auto pnode : _nlist) {
			if (pnode->Size() > rmax) {
				rmax = pnode->Size();
			}
		}
		rmax *= 1.5 * 3.0;

		//----------仮想四面体の生成----------
		Node* nst0 = new Node(rmax * 2.0*sqrt(2.0) / 3.0 *cos(2.0*0.0*M_PI / 3.0), rmax * 2.0*sqrt(2.0) / 3.0 *sin(2.0*0.0*M_PI / 3.0), -rmax / 3.0, -1);
		Node* nst1 = new Node(rmax * 2.0*sqrt(2.0) / 3.0 *cos(2.0*1.0*M_PI / 3.0), rmax * 2.0*sqrt(2.0) / 3.0 *sin(2.0*1.0*M_PI / 3.0), -rmax / 3.0, -1);
		Node* nst2 = new Node(rmax * 2.0*sqrt(2.0) / 3.0 *cos(2.0*2.0*M_PI / 3.0), rmax * 2.0*sqrt(2.0) / 3.0 *sin(2.0*2.0*M_PI / 3.0), -rmax / 3.0, -1);
		Node* nst3 = new Node(0.0, 0.0, rmax, -1);

		_nlist.push_back(nst0);		_nlist.push_back(nst1);		_nlist.push_back(nst2);		_nlist.push_back(nst3);
		_elist.push_back(new Element(nst0, nst1, nst2, nst3));
	}


	//**********局所Delaunay分割**********
	//局所的にDelaunay分割した後末尾の要素のidを返す
	Element* MeshLocal(Node *_node, Element *_pethis, std::vector<Element*> &_elist) {
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
		bool is_anysurface_invalid = true;
		while (is_anysurface_invalid) {
			is_anysurface_invalid = false;

			for (int i = 0; i < sstack.size(); i++) {
				if (sstack[i]->IsActive) {
					Node vn = (*sstack[i]->pnodes[1] - *sstack[i]->pnodes[0]) * (*sstack[i]->pnodes[2] - *sstack[i]->pnodes[0]);
					Node vp = *_node - *sstack[i]->pnodes[0];
					double D = vn ^ vp;

					//----------不良な面がある場合----------
					if (D < EPS) {
						Element* peadd = sstack[i]->pneighbor;			//不良な面に隣接する要素を指すポインタ

						//----------新しく追加できる要素がある場合----------
						if (peadd != nullptr) {
							if (peadd->IsActive) {
								is_anysurface_invalid = true;
								stack.push_back(peadd);
								peadd->IsActive = false;

								//----------追加した要素の各面について隣接面をFalseに----------
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
						else if (fabs(D) < EPS) {
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
			}
		}

		//----------新しく生成された要素同士の隣接関係を計算----------
		for (auto pelement : penew) {
			for (auto psurface : pelement->psurfaces) {
				if (psurface->pneighbor == nullptr) {
					for (auto pelement2 : penew) {
						if (pelement != pelement2) {
							for (auto psurface2 : pelement2->psurfaces) {
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
		
		//----------新しく生成された要素をelistに追加----------
		for (auto pelement : penew) {
			_elist.push_back(pelement);
		}

		//----------stack内の古い要素を削除----------
		for (auto& pelement : _elist) {
			if (!pelement->IsActive) {
				auto tmp = pelement;
				pelement = *(_elist.end() - 1);
				_elist.pop_back();
				delete tmp;
			}
		}
		
		return *(_elist.end() - 1);
	}


	//**********粗いDelaunay分割**********
	void MakeRoughMesh(std::vector<Node*> _nlist, std::vector<Element*> &_elist) {
		std::cout << "Make rough mash\n";

		Element* pethis = _elist[0];										//現在調べている要素を指すポインタ
		for (auto& pnode : _nlist) {
			if (pnode->type != -1) {
				while (1) {
					int count = 0;
					Element* penext = pethis->GetLocateId(pnode);				//次に調べる要素を指すポインタ
					//----------要素内に点があるとき----------
					if (penext == pethis) {
						std::cout << "at\t" << pethis << "\n";
						pethis = MeshLocal(pnode, pethis, _elist);
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
		for (auto& pelement : _elist) {
			if (pelement->pnodes[0]->type == -1 || pelement->pnodes[1]->type == -1 || pelement->pnodes[2]->type == -1 || pelement->pnodes[3]->type == -1) {
				for (auto& psurface : pelement->psurfaces) {
					if (psurface->pneighbor != nullptr) {
						psurface->pneighbor->GetAdjacentSurface(psurface->pparent)->pneighbor = nullptr;
					}
				}
				auto tmp = pelement;
				pelement = *(_elist.end() - 1);
				_elist.pop_back();
				delete tmp;
			}
		}
	}
}