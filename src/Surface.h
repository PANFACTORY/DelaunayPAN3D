//*********************************************************
// Title		:src/cpp/Surface.h
// Author	:Tanabe Yuta
// Date		:2019/01/26
// Copyright	:(C)2019 TanabeYuta
//*********************************************************

#pragma once
#include <array>

namespace DelaunayPAN3D {
template <class N, class T>
class Element;

template <class N, class T>
class Surface {
   public:
    Surface(){};
    ~Surface(){};
    Surface(N* _pnode0, N* _pnode1, N* _pnode2, Element<N, T>* _pparent,
            Element<N, T>* _pneighbor) {
        this->pnodes[0] = _pnode0;
        this->pnodes[1] = _pnode1;
        this->pnodes[2] = _pnode2;
        this->pparent = _pparent;
        this->pneighbor = _pneighbor;
        this->IsActive = true;
    }

    Element<N, T>* pneighbor;
    bool IsActive;
    std::array<N*, 3> pnodes;
    Element<N, T>* pparent;

    bool operator==(const Surface<N, T>& _surface) {
        if ((this->pnodes[0] == _surface.pnodes[0] &&
             this->pnodes[1] == _surface.pnodes[2] &&
             this->pnodes[2] == _surface.pnodes[1]) ||
            (this->pnodes[0] == _surface.pnodes[1] &&
             this->pnodes[1] == _surface.pnodes[0] &&
             this->pnodes[2] == _surface.pnodes[2]) ||
            (this->pnodes[0] == _surface.pnodes[2] &&
             this->pnodes[1] == _surface.pnodes[1] &&
             this->pnodes[2] == _surface.pnodes[0])) {
            return true;
        }
        return false;
    }

    bool IsRayCross(const N& _sp, const N& _dir, T EPS) const {
        N v01 = *(this->pnodes[1]) - *(this->pnodes[0]);
        N v02 = *(this->pnodes[2]) - *(this->pnodes[0]);
        N v0g = _sp - *(this->pnodes[0]);

        T det = v01.dot(v02.cross(_dir));
        if (det > EPS) {
            T u = (v0g.dot(v02.cross(_dir))) / det;
            if (-EPS < u && u < 1.0 + EPS) {
                T v = (v01.dot(v0g.cross(_dir))) / det;
                if (-EPS < v && u + v < 1.0 + EPS) {
                    T t = (v01.dot(v02.cross(v0g))) / det;
                    if (t > -EPS && t < 1.0 - EPS) {
                        return true;
                    }
                }
            }
        }
        return false;
    }
};
}  // namespace DelaunayPAN3D
