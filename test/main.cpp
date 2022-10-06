//*********************************************************
// Title		:src/cpp/main.cpp
// Author	:Tanabe Yuta
// Date		:2019/01/26
// Copyright	:(C)2019 TanabeYuta
//*********************************************************
//座標を正規化した後元に戻す作業が必要

#include <chrono>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include "../src/Delaunay.h"
#include "../src/vec3.h"

using namespace DelaunayPAN3D;

//*****************************************************************************
//	Import nodes
//*****************************************************************************
bool importnode(std::vector<Vec3<double>*>& _pnodes, std::string _fname,
                int _type) {
    std::ifstream fin(_fname);
    if (!fin) {
        std::cout << "Node File Open Error : " << _fname << "\n";
        return false;
    }

    std::string tmp;
    while (getline(fin, tmp)) {
        std::istringstream ssi(tmp);
        std::string tmpx[3];
        ssi >> tmpx[0] >> tmpx[1] >> tmpx[2];

        bool is_same_node = false;
        for (auto pnode : _pnodes) {
            if (Vec3<double>(stod(tmpx[0]), stod(tmpx[1]), stod(tmpx[2]), -1) ==
                *pnode) {
                is_same_node = true;
                break;
            }
        }
        if (!is_same_node) {
            _pnodes.push_back(new Vec3<double>(stod(tmpx[0]), stod(tmpx[1]),
                                               stod(tmpx[2]), _type));
        }
    }
    fin.close();
    return true;
}

//*****************************************************************************
//	Export to VTK
//*****************************************************************************
// void exportvtk(std::vector<Vec3<double>*> _pnodes,
//                std::vector<Element<Vec3<double>, double>*> _pelements,
//                std::string _fname) {
//     // Convert Vec3* to index
//     std::map<Vec3<double>*, int> index;
//     int idx = 0;
//     for (auto pnode : _pnodes) {
//         index[pnode] = idx++;
//     }

//     std::ofstream fout(_fname + ".vtk");

//     fout << "# vtk DataFile Version 4.1\n"
//          << _fname << "\nASCII\nDATASET UNSTRUCTURED_GRID\n";
//     fout << "POINTS\t" << _pnodes.size() << "\tfloat\n";
//     for (auto pnode : _pnodes) {
//         fout << (*pnode)[0] << "\t" << (*pnode)[1] << "\t" << (*pnode)[2]
//              << "\n";
//     }
//     fout << "CELLS\t" << _pelements.size() << "\t" << _pelements.size() * 5
//          << "\n";
//     for (auto pelement : _pelements) {
//         fout << "4\t" << index[pelement->pnodes[0]] << "\t"
//              << index[pelement->pnodes[1]] << "\t" <<
//              index[pelement->pnodes[2]]
//              << "\t" << index[pelement->pnodes[3]] << "\n";
//     }
//     fout << "CELL_TYPES\t" << _pelements.size() << "\n";
//     for (int i = 0; i < _pelements.size(); i++) {
//         fout << "10\n";
//     }

//     fout.close();
// }

//*****************************************************************************
//	main
//*****************************************************************************
int main() {
    std::string filepath = "sample/Model5";

    std::cout << "**********Delaunay trianguration**********\n";

    //----------List of nodes and elements----------
    std::vector<Vec3<double>*> pnodes;
    std::vector<Element<Vec3<double>, double>*> pelements;

    //----------Import nodes----------
    if (!importnode(pnodes, filepath + "/node.dat", 0)) {
        return -1;
    }
    bool IsCopynodeExist = importnode(pnodes, filepath + "/copynode.dat", 1);

    //----------Make mesh----------
    std::chrono::system_clock::time_point start =
        std::chrono::system_clock::now();

    MakeMesh<Vec3<double>, double>(pnodes, pelements, 5000, IsCopynodeExist,
                                   4.0, 1e-15);

    std::chrono::system_clock::time_point end =
        std::chrono::system_clock::now();

    //----------Export results----------
    std::cout << "**********The Result**********\n";
    std::cout << "time cost:\t"
              << std::chrono::duration_cast<std::chrono::milliseconds>(end -
                                                                       start)
                         .count() /
                     1000.0
              << "sec.\n";
    std::cout << "node:     \t" << pnodes.size() << "\n";
    std::cout << "element:  \t" << pelements.size() << "\n";
    std::cout << "Export to:\t" << filepath << "/mesh.vtk\n";
    std::cout << "******************************\n";
    // exportvtk(pnodes, pelements, filepath + "/mesh");

    //----------Release memory----------
    for (auto pelement : pelements) {
        delete pelement;
    }
    pelements.clear();
    for (auto pnode : pnodes) {
        delete pnode;
    }
    pnodes.clear();

    return 0;
}
