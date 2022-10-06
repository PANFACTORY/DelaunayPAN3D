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
#include "Node.h"
#include "Surface.h"

namespace DelaunayPAN3D {
//*************************************************************************
//	Make Mesh
//*************************************************************************
template <class T>
void MakeMesh(std::vector<Node<T>*>& _pnodes,
              std::vector<Element<T>*>& _pelements, int _addnodenum,
              bool _iscopynodeexist, T ALPHA, T EPS) {
    //----------Get region which nodes exist----------
    T xmax = T(), xmin = T(), ymax = T(), ymin = T(), zmax = T(), zmin = T();
    for (auto pnode : _pnodes) {
        if (pnode->x > xmax) {
            xmax = pnode->x;
        }
        if (pnode->x < xmin) {
            xmin = pnode->x;
        }
        if (pnode->y > ymax) {
            ymax = pnode->y;
        }
        if (pnode->y < ymin) {
            ymin = pnode->y;
        }
        if (pnode->z > zmax) {
            zmax = pnode->z;
        }
        if (pnode->z < zmin) {
            zmin = pnode->z;
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
        pnode->x =
            (pnode->x - xmin) / dmax + 0.5 * (ALPHA - 1.0) * xrange / dmax;
        pnode->y =
            (pnode->y - ymin) / dmax + 0.5 * (ALPHA - 1.0) * yrange / dmax;
        pnode->z =
            (pnode->z - zmin) / dmax + 0.5 * (ALPHA - 1.0) * zrange / dmax;
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
        pnode->x = pnode->x * dmax - 0.5 * (ALPHA - 1.0) * xrange + xmin;
        pnode->y = pnode->y * dmax - 0.5 * (ALPHA - 1.0) * yrange + ymin;
        pnode->z = pnode->z * dmax - 0.5 * (ALPHA - 1.0) * zrange + zmin;
    }
}

//*************************************************************************
//	Make supertetrahedron
//*************************************************************************
template <class T>
void MakeSupertetrahedron(std::vector<Node<T>*>& _pnodes,
                          std::vector<Element<T>*>& _pelements, T _xmax,
                          T _ymax, T _zmax) {
    std::cout << "Make supertetraedron\n";

    //----------Make nodes of supertetrahedron----------
    Node<T>* nst0 = new Node<T>(T(), T(), T(), -1, _pnodes.size());
    _pnodes.push_back(nst0);
    Node<T>* nst1 = new Node<T>(_xmax, T(), T(), -1, _pnodes.size());
    _pnodes.push_back(nst1);
    Node<T>* nst2 = new Node<T>(_xmax, _ymax, T(), -1, _pnodes.size());
    _pnodes.push_back(nst2);
    Node<T>* nst3 = new Node<T>(T(), _ymax, T(), -1, _pnodes.size());
    _pnodes.push_back(nst3);
    Node<T>* nst4 = new Node<T>(T(), T(), _zmax, -1, _pnodes.size());
    _pnodes.push_back(nst4);
    Node<T>* nst5 = new Node<T>(_xmax, T(), _zmax, -1, _pnodes.size());
    _pnodes.push_back(nst5);
    Node<T>* nst6 = new Node<T>(_xmax, _ymax, _zmax, -1, _pnodes.size());
    _pnodes.push_back(nst6);
    Node<T>* nst7 = new Node<T>(T(), _ymax, _zmax, -1, _pnodes.size());
    _pnodes.push_back(nst7);

    //----------Make elements of supertetrahedron----------
    _pelements.push_back(new Element<T>(nst1, nst3, nst0, nst7));
    _pelements.push_back(new Element<T>(nst2, nst1, nst6, nst7));
    _pelements.push_back(new Element<T>(nst2, nst3, nst1, nst7));
    _pelements.push_back(new Element<T>(nst1, nst5, nst6, nst7));
    _pelements.push_back(new Element<T>(nst1, nst0, nst5, nst7));
    _pelements.push_back(new Element<T>(nst4, nst5, nst0, nst7));

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
template <class T>
void MeshLocal(Node<T>* _pnode, Element<T>* _pethis,
               std::vector<Element<T>*>& _pelements, T EPS) {
    std::vector<Element<T>*> stack, substack;
    std::vector<Surface<T>*> sstack;

    //----------Get elements which node is in----------
    substack.push_back(_pethis);
    while (substack.size()) {
        Element<T>* pend = *(substack.end() - 1);
        substack.pop_back();

        if (pend->IsActive) {
            stack.push_back(pend);
            pend->IsActive = false;

            for (auto& psurface : pend->psurfaces) {
                Element<T>* pneighbor = psurface->pneighbor;
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
                Element<T> D =
                    Element<T>(sstack[i]->pnodes[0], sstack[i]->pnodes[1],
                               sstack[i]->pnodes[2], _pnode);

                //----------if there are crevices----------
                if (D.volume < EPS) {
                    Element<T>* peadd = sstack[i]->pneighbor;

                    //----------if able to add elements----------
                    if (peadd != nullptr) {
                        if (peadd->IsActive) {
                            is_anysurface_invalid = true;
                            peadd->IsActive = false;
                            stack.push_back(peadd);

                            //----------make surfaces isactive false----------
                            for (auto& psurface : peadd->psurfaces) {
                                Element<T>* pneighbor = psurface->pneighbor;
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
    std::vector<Element<T>*> penew;
    for (auto& psurface : sstack) {
        if (psurface->IsActive) {
            Element<T>* tmp =
                new Element<T>(psurface->pnodes[0], psurface->pnodes[1],
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
template <class T>
void MakeRoughMesh(std::vector<Node<T>*> _pnodes,
                   std::vector<Element<T>*>& _pelements, T EPS) {
    std::cout << "Make rough mesh\n";

    Element<T>* pethis = _pelements[0];
    for (auto& pnode : _pnodes) {
        if (pnode->type != -1) {
            int count = 0;
            while (1) {
                Element<T>* penext = pethis->GetLocateId(pnode, EPS);
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
template <class T>
void DeleteSupertetrahedron(std::vector<Element<T>*>& _pelements) {
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
template <class T>
void DeleteCreviceElement(std::vector<Element<T>*>& _pelements) {
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
template <class T>
void MakeFineMesh(std::vector<Node<T>*>& _pnodes,
                  std::vector<Element<T>*>& _pelements, int _addnodenum,
                  T EPS) {
    std::cout << "Make fine mesh\n";

    for (int i = 0; i < _addnodenum; i++) {
        //----------Find element which has longest edge----------
        T edgelengthmax = T();
        Element<T>* pethis = nullptr;
        Node<T>* pnode0 = nullptr;
        Node<T>* pnode1 = nullptr;

        for (auto pelement : _pelements) {
            for (int j = 0; j < 3; j++) {
                for (int k = j + 1; k < 3; k++) {
                    T edgelength =
                        (*pelement->pnodes[k] - *pelement->pnodes[j]).Norm();
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
        Node<T>* tmp = new Node<T>((*pnode0 + *pnode1) / 2.0);
        tmp->type = 2;
        tmp->id = _pnodes.size();
        _pnodes.push_back(tmp);
        MeshLocal(tmp, pethis, _pelements, EPS);
    }
}
}  // namespace DelaunayPAN3D
