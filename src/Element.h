//*********************************************************
// Title		:src/cpp/Element.h
// Author	:Tanabe Yuta
// Date		:2019/01/26
// Copyright	:(C)2019 TanabeYuta
//*********************************************************

#pragma once
#include <fenv.h>

#include <array>
#include <cmath>

#include "Surface.h"

namespace DelaunayPAN3D {
template <class N, class T>
class Surface;

template <class N, class T>
class Element {
   public:
    Element() {}
    ~Element() {
        for (auto& psurface : this->psurfaces) {
            delete psurface;
        }
    }
    Element(N* _pnode0, N* _pnode1, N* _pnode2, N* _pnode3) {
        this->IsActive = true;

        //----------Set nodes----------
        this->pnodes[0] = _pnode0;
        this->pnodes[1] = _pnode1;
        this->pnodes[2] = _pnode2;
        this->pnodes[3] = _pnode3;

        //----------Set surfaces----------
        this->psurfaces[0] =
            new Surface<N, T>(_pnode1, _pnode3, _pnode2, this, nullptr);
        this->psurfaces[1] =
            new Surface<N, T>(_pnode0, _pnode2, _pnode3, this, nullptr);
        this->psurfaces[2] =
            new Surface<N, T>(_pnode0, _pnode3, _pnode1, this, nullptr);
        this->psurfaces[3] =
            new Surface<N, T>(_pnode0, _pnode1, _pnode2, this, nullptr);

        //----------Get center and radius of external sphere----------
        N v0 = *_pnode1 - *_pnode0;
        N v1 = *_pnode2 - *_pnode0;
        N v2 = *_pnode3 - *_pnode0;

        N ABC =
            N(0.5 * ((*_pnode1).dot(*_pnode1) - (*_pnode0).dot(*_pnode0)),
              0.5 * ((*_pnode2).dot(*_pnode2) - (*_pnode0).dot(*_pnode0)),
              0.5 * ((*_pnode3).dot(*_pnode3) - (*_pnode0).dot(*_pnode0)), -1);

        T detP = v0.dot(v1.cross(v2));
        N P0 = v1.cross(v2);
        N P1 = v2.cross(v0);
        N P2 = v0.cross(v1);

        this->scenter =
            N((ABC[0] * P0[0] + ABC[1] * P1[0] + ABC[2] * P2[0]) / detP,
              (ABC[0] * P0[1] + ABC[1] * P1[1] + ABC[2] * P2[1]) / detP,
              (ABC[0] * P0[2] + ABC[1] * P1[2] + ABC[2] * P2[2]) / detP, -1);
        this->sround = (this->scenter - *_pnode0).norm();

        //----------Get center of gravity----------
        this->gcenter = (*_pnode0 + *_pnode1 + *_pnode2 + *_pnode3) / 4.0;

        //----------Get volume----------
        this->volume = ((*_pnode1 - *_pnode0).cross(*_pnode2 - *_pnode0))
                           .dot(*_pnode3 - *_pnode0);

        //----------Get aspect ratio----------
        this->aspect =
            this->volume / pow(this->sround, 3.0) / (8.0 * sqrt(3.0) / 27.0);
        if (fetestexcept(FE_DIVBYZERO)) {
            feclearexcept(FE_ALL_EXCEPT);
            this->aspect = T();
        }
    }

    bool IsActive;
    std::array<Surface<N, T>*, 4> psurfaces;
    std::array<N*, 4> pnodes;
    N scenter, gcenter;
    T sround, volume, aspect;

    Element<N, T>* GetLocateId(N* _pnode, T EPS) {
        for (auto surface : this->psurfaces) {
            if (surface->IsRayCross(this->gcenter, this->gcenter - *_pnode,
                                    EPS)) {
                return surface->pneighbor;
            }
        }
        return this;
    }
    bool IsInSphere(N* _pnode, T EPS) {
        return sround + EPS > (this->scenter - *_pnode).norm();
    }
    Surface<N, T>* GetAdjacentSurface(Element<N, T>* _pelement) {
        for (auto& psurface : this->psurfaces) {
            if (psurface->pneighbor == _pelement) {
                return psurface;
            }
        }
        return nullptr;
    }
};
}  // namespace DelaunayPAN3D
