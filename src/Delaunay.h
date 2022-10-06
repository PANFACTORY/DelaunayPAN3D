//*********************************************************
// Title		:src/cpp/Delaunay.h
// Author	:Tanabe Yuta
// Date		:2019/01/08
// Copyright	:(C)2019 TanabeYuta
//*********************************************************

#pragma once
#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <vector>

#include "Element.h"
#include "Surface.h"

namespace DelaunayPAN3D {
//*************************************************************************
//	Make Mesh
//*************************************************************************
template <class N, class T>
void MakeMesh(std::vector<N*>& _pnodes, std::vector<Element<N, T>*>& _pelements,
              int _addnodenum, bool _iscopynodeexist, T ALPHA, T EPS) {
    //----------Get region which nodes exist----------
    T xmax = T(), xmin = T(), ymax = T(), ymin = T(), zmax = T(), zmin = T();
    for (auto pnode : _pnodes) {
        if ((*pnode)[0] > xmax) {
            xmax = (*pnode)[0];
        }
        if ((*pnode)[0] < xmin) {
            xmin = (*pnode)[0];
        }
        if ((*pnode)[1] > ymax) {
            ymax = (*pnode)[1];
        }
        if ((*pnode)[1] < ymin) {
            ymin = (*pnode)[1];
        }
        if ((*pnode)[2] > zmax) {
            zmax = (*pnode)[2];
        }
        if ((*pnode)[2] < zmin) {
            zmin = (*pnode)[2];
        }
    }

    //----------Normalize cordinate----------
    T xrange = 0.5 * (xmax - xmin), yrange = 0.5 * (ymax - ymin),
      zrange = 0.5 * (zmax - zmin);
    T dmax = xrange;
    if (dmax < yrange) {
        dmax = yrange;
    }
    if (dmax < zrange) {
        dmax = zrange;
    }
    for (auto& pnode : _pnodes) {
        (*pnode)[0] =
            ((*pnode)[0] - xmin) / dmax + 0.5 * (ALPHA - 1.0) * xrange / dmax;
        (*pnode)[1] =
            ((*pnode)[1] - ymin) / dmax + 0.5 * (ALPHA - 1.0) * yrange / dmax;
        (*pnode)[2] =
            ((*pnode)[2] - zmin) / dmax + 0.5 * (ALPHA - 1.0) * zrange / dmax;
    }

    //----------Make supertetrahedron----------
    T x = ALPHA * xrange / dmax, y = ALPHA * yrange / dmax,
      z = ALPHA * zrange / dmax;
    MakeSupertetrahedron(_pnodes, _pelements, x, y, z);

    //----------Make rough mesh----------
    MakeRoughMesh(_pnodes, _pelements, EPS);

    //----------Delete needless elements----------
    DeleteSupertetrahedron(_pelements);
    if (_iscopynodeexist) {
        DeleteCreviceElement(_pelements);
    }

    //----------Make fine mesh----------
    MakeFineMesh(_pnodes, _pelements, _addnodenum, EPS);

    //----------Renormalize cordinate----------
    for (auto& pnode : _pnodes) {
        (*pnode)[0] = (*pnode)[0] * dmax - 0.5 * (ALPHA - 1.0) * xrange + xmin;
        (*pnode)[1] = (*pnode)[1] * dmax - 0.5 * (ALPHA - 1.0) * yrange + ymin;
        (*pnode)[2] = (*pnode)[2] * dmax - 0.5 * (ALPHA - 1.0) * zrange + zmin;
    }
}

//*************************************************************************
//	Make supertetrahedron
//*************************************************************************
template <class N, class T>
void MakeSupertetrahedron(std::vector<N*>& _pnodes,
                          std::vector<Element<N, T>*>& _pelements, T _xmax,
                          T _ymax, T _zmax) {
    std::cout << "Make supertetraedron\n";

    //----------Make nodes of supertetrahedron----------
    N* nst0 = new N(T(), T(), T(), -1, _pnodes.size());
    _pnodes.push_back(nst0);
    N* nst1 = new N(_xmax, T(), T(), -1, _pnodes.size());
    _pnodes.push_back(nst1);
    N* nst2 = new N(_xmax, _ymax, T(), -1, _pnodes.size());
    _pnodes.push_back(nst2);
    N* nst3 = new N(T(), _ymax, T(), -1, _pnodes.size());
    _pnodes.push_back(nst3);
    N* nst4 = new N(T(), T(), _zmax, -1, _pnodes.size());
    _pnodes.push_back(nst4);
    N* nst5 = new N(_xmax, T(), _zmax, -1, _pnodes.size());
    _pnodes.push_back(nst5);
    N* nst6 = new N(_xmax, _ymax, _zmax, -1, _pnodes.size());
    _pnodes.push_back(nst6);
    N* nst7 = new N(T(), _ymax, _zmax, -1, _pnodes.size());
    _pnodes.push_back(nst7);

    //----------Make elements of supertetrahedron----------
    _pelements.push_back(new Element<N, T>(nst1, nst3, nst0, nst7));
    _pelements.push_back(new Element<N, T>(nst2, nst1, nst6, nst7));
    _pelements.push_back(new Element<N, T>(nst2, nst3, nst1, nst7));
    _pelements.push_back(new Element<N, T>(nst1, nst5, nst6, nst7));
    _pelements.push_back(new Element<N, T>(nst1, nst0, nst5, nst7));
    _pelements.push_back(new Element<N, T>(nst4, nst5, nst0, nst7));

    //----------Make connection of supertetrahedron----------
    for (auto& pelement : _pelements) {
        for (auto& psurface : pelement->psurfaces) {
            if (psurface->pneighbor == nullptr) {
                for (auto& pelement2 : _pelements) {
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

//*************************************************************************
//	Mesh local
//*************************************************************************
template <class N, class T>
void MeshLocal(N* _pnode, Element<N, T>* _pethis,
               std::vector<Element<N, T>*>& _pelements, T EPS) {
    std::vector<Element<N, T>*> stack, substack;
    std::vector<Surface<N, T>*> sstack;

    //----------Get elements which node is in----------
    substack.push_back(_pethis);
    while (substack.size()) {
        Element<N, T>* pend = *(substack.end() - 1);
        substack.pop_back();

        if (pend->IsActive) {
            stack.push_back(pend);
            pend->IsActive = false;

            for (auto& psurface : pend->psurfaces) {
                Element<N, T>* pneighbor = psurface->pneighbor;
                if (pneighbor != nullptr &&
                    pneighbor->IsInSphere(_pnode, EPS)) {
                    substack.push_back(pneighbor);
                } else {
                    sstack.push_back(psurface);
                }
            }
        }
    }

    //----------Modify crevice of polyhedron----------
    bool is_anysurface_invalid = true;
    while (is_anysurface_invalid) {
        is_anysurface_invalid = false;

        for (int i = 0; i < sstack.size(); i++) {
            if (sstack[i]->IsActive) {
                Element<N, T> D(sstack[i]->pnodes[0], sstack[i]->pnodes[1],
                                sstack[i]->pnodes[2], _pnode);

                //----------if there are crevices----------
                if (D.volume < EPS) {
                    Element<N, T>* peadd = sstack[i]->pneighbor;

                    //----------if able to add elements----------
                    if (peadd != nullptr) {
                        if (peadd->IsActive) {
                            is_anysurface_invalid = true;
                            peadd->IsActive = false;
                            stack.push_back(peadd);

                            //----------make surfaces isactive false----------
                            for (auto& psurface : peadd->psurfaces) {
                                Element<N, T>* pneighbor = psurface->pneighbor;
                                if (pneighbor != nullptr &&
                                    !pneighbor->IsActive) {
                                    pneighbor->GetAdjacentSurface(peadd)
                                        ->IsActive = false;
                                } else {
                                    sstack.push_back(psurface);
                                }
                            }
                            break;
                        }
                    } else if (fabs(D.volume) < EPS) {
                        sstack[i]->IsActive = false;
                    }
                }
            }
        }
    }

    //----------Make new elements----------
    std::vector<Element<N, T>*> penew;
    for (auto& psurface : sstack) {
        if (psurface->IsActive) {
            Element<N, T>* tmp =
                new Element<N, T>(psurface->pnodes[0], psurface->pnodes[1],
                                  psurface->pnodes[2], _pnode);
            tmp->psurfaces[3]->pneighbor = psurface->pneighbor;
            if (psurface->pneighbor != nullptr) {
                psurface->pneighbor->GetAdjacentSurface(psurface->pparent)
                    ->pneighbor = tmp;
            }
            penew.push_back(tmp);
            _pelements.push_back(tmp);
        }
    }

    //----------Make connection of new elements----------
    for (auto& pelement : penew) {
        for (auto& psurface : pelement->psurfaces) {
        OUT:
            if (psurface->pneighbor == nullptr) {
                for (auto& pelement2 : penew) {
                    for (auto& psurface2 : pelement2->psurfaces) {
                        if (*psurface == *psurface2) {
                            //----------if invalid element is made----------
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

    //----------Delete needless elements in stack----------
    for (auto pelement = _pelements.begin(); pelement != _pelements.end();) {
        if (!(*pelement)->IsActive) {
            delete *pelement;
            pelement = _pelements.erase(pelement);
        } else {
            ++pelement;
        }
    }
}

//*************************************************************************
//	Make rough mesh
//*************************************************************************
template <class N, class T>
void MakeRoughMesh(std::vector<N*> _pnodes,
                   std::vector<Element<N, T>*>& _pelements, T EPS) {
    std::cout << "Make rough mesh\n";

    Element<N, T>* pethis = _pelements[0];
    for (auto& pnode : _pnodes) {
        if (pnode->type != -1) {
            int count = 0;
            while (1) {
                Element<N, T>* penext = pethis->GetLocateId(pnode, EPS);
                //----------if node is in the element----------
                if (penext == pethis) {
                    MeshLocal(pnode, pethis, _pelements, EPS);
                    pethis = *(_pelements.end() - 1);
                    break;
                } else {
                    pethis = penext;
                }
            }
        }
    }
}

//*************************************************************************
//	Delete supertetrahedron
//*************************************************************************
template <class N, class T>
void DeleteSupertetrahedron(std::vector<Element<N, T>*>& _pelements) {
    std::cout << "Delete supertetraedron\n";

    for (auto pelement = _pelements.begin(); pelement != _pelements.end();) {
        if ((*pelement)->pnodes[0]->type == -1 ||
            (*pelement)->pnodes[1]->type == -1 ||
            (*pelement)->pnodes[2]->type == -1 ||
            (*pelement)->pnodes[3]->type == -1) {
            for (auto& psurface : (*pelement)->psurfaces) {
                if (psurface->pneighbor != nullptr) {
                    psurface->pneighbor->GetAdjacentSurface(psurface->pparent)
                        ->pneighbor = nullptr;
                }
            }
            delete *pelement;
            pelement = _pelements.erase(pelement);
        } else {
            ++pelement;
        }
    }
}

//*************************************************************************
//	Delete elements at crevice part
//*************************************************************************
template <class N, class T>
void DeleteCreviceElement(std::vector<Element<N, T>*>& _pelements) {
    std::cout << "Delete Crevice Element\n";

    for (auto pelement = _pelements.begin(); pelement != _pelements.end();) {
        if ((*pelement)->pnodes[0]->type == (*pelement)->pnodes[1]->type &&
            (*pelement)->pnodes[1]->type == (*pelement)->pnodes[2]->type &&
            (*pelement)->pnodes[2]->type == (*pelement)->pnodes[3]->type) {
            for (auto& psurface : (*pelement)->psurfaces) {
                if (psurface->pneighbor != nullptr) {
                    psurface->pneighbor->GetAdjacentSurface(psurface->pparent)
                        ->pneighbor = nullptr;
                }
            }
            delete *pelement;
            pelement = _pelements.erase(pelement);
        } else {
            ++pelement;
        }
    }
}

//*************************************************************************
//	Make fine mesh
//*************************************************************************
template <class N, class T>
void MakeFineMesh(std::vector<N*>& _pnodes,
                  std::vector<Element<N, T>*>& _pelements, int _addnodenum,
                  T EPS) {
    std::cout << "Make fine mesh\n";

    for (int i = 0; i < _addnodenum; i++) {
        //----------Find element which has longest edge----------
        T edgelengthmax = T();
        Element<N, T>* pethis = nullptr;
        N* pnode0 = nullptr;
        N* pnode1 = nullptr;

        for (auto pelement : _pelements) {
            for (int j = 0; j < 3; j++) {
                for (int k = j + 1; k < 3; k++) {
                    T edgelength =
                        (*pelement->pnodes[k] - *pelement->pnodes[j]).norm();
                    if (edgelength > edgelengthmax) {
                        edgelengthmax = edgelength;
                        pethis = pelement;
                        pnode0 = pelement->pnodes[j];
                        pnode1 = pelement->pnodes[k];
                    }
                }
            }
        }

        //----------Add Node----------
        N* tmp = new N((*pnode0 + *pnode1) / 2.0);
        tmp->type = 2;
        tmp->id = _pnodes.size();
        _pnodes.push_back(tmp);
        MeshLocal(tmp, pethis, _pelements, EPS);
    }
}
}  // namespace DelaunayPAN3D
