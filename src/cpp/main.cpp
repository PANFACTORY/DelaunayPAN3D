//*********************************************************
//Title		:src/cpp/main.cpp
//Author	:Tanabe Yuta
//Date		:2019/01/26
//Copyright	:(C)2019 TanabeYuta
//*********************************************************
//座標を正規化した後元に戻す作業が必要


#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <chrono>


#include "Delaunay.h"
#include "Parameter.h"


using namespace DelaunayPAN3D;


//*****************************************************************************
//	Import nodes
//*****************************************************************************
bool importnode(std::vector<Node<double>*>& _pnodes, std::string _fname, int _type) {
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
			if (Node<double>(stod(tmpx[0]), stod(tmpx[1]), stod(tmpx[2]), -1, -1) == *pnode) {
				is_same_node = true;
				break;
			}
		}
		if (!is_same_node) {
			_pnodes.push_back(new Node<double>(stod(tmpx[0]), stod(tmpx[1]), stod(tmpx[2]), _type, _pnodes.size()));
		}
	}
	fin.close();
	return true;
}


//*****************************************************************************
//	Export to VTK
//*****************************************************************************
void exportvtk(std::vector<Node<double>*> _pnodes, std::vector<Element<double>*> _pelements, std::string _fname) {
	std::ofstream fout(_fname + ".vtk");

	fout << "# vtk DataFile Version 4.1\n" << _fname << "\nASCII\nDATASET UNSTRUCTURED_GRID\n";
	fout << "POINTS\t" << _pnodes.size() << "\tfloat\n";
	for (auto pnode : _pnodes) {
		fout << pnode->x << "\t" << pnode->y << "\t" << pnode->z << "\n";
	}
	fout << "CELLS\t" << _pelements.size() << "\t" << _pelements.size() * 5 << "\n";
	for (auto pelement : _pelements) {
		fout << "4\t" << pelement->pnodes[0]->id << "\t" << pelement->pnodes[1]->id << "\t" << pelement->pnodes[2]->id << "\t" << pelement->pnodes[3]->id << "\n";
	}
	fout << "CELL_TYPES\t" << _pelements.size() << "\n";
	for (int i = 0; i < _pelements.size(); i++) {
		fout << "10\n";
	}

	fout.close();
}


//*****************************************************************************
//	main
//*****************************************************************************
int main() {
	std::string filepath = "sample/Model5";

	std::cout << "**********Delaunay trianguration**********\n";

	//----------List of nodes and elements----------
	std::vector<Node<double>*> pnodes;
	std::vector<Element<double>*> pelements;
	
	//----------Import nodes----------
	if (!importnode(pnodes, filepath + "/node.dat", 0)) {	return -1;	}
	bool IsCopynodeExist = importnode(pnodes, filepath + "/copynode.dat", 1);

	//----------Make mesh----------
	std::chrono::system_clock::time_point start = std::chrono::system_clock::now();
	
	MakeMesh<double>(pnodes, pelements, 10000, IsCopynodeExist);
	
	std::chrono::system_clock::time_point end = std::chrono::system_clock::now();
	
	//----------Export results----------
	std::cout << "**********The Result**********\n";
	std::cout << "time cost:\t" << std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count() / 1000.0 << "sec.\n";
	std::cout << "node:     \t" << pnodes.size() << "\n";
	std::cout << "element:  \t" << pelements.size() << "\n";
	std::cout << "Export to:\t" << filepath << "/mesh.vtk\n";
	std::cout << "******************************\n";
	exportvtk(pnodes, pelements, filepath + "/mesh");

	//----------Release memory----------
	for (auto pelement : pelements) {	delete pelement;	}
	pelements.clear();
	for (auto pnode : pnodes) {	delete pnode;	}
	pnodes.clear();
	
	return 0;
}