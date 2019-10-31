//*********************************************************
//Title		:Delaunay3Dver3
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
#include <time.h>


#include "Delaunay.h"
#include "Parameter.h"


using namespace DelaunayPAN3DV3;


//**********節点の取り込み**********
bool importnode(std::vector<Node<double>*> &_nlist, std::string _fname, int _type) {
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
		for (auto pnode : _nlist) {
			if (Node<double>(stod(tmpx[0]), stod(tmpx[1]), stod(tmpx[2]), -1, -1) == *pnode) {
				is_same_node = true;
				break;
			}
		}
		if (!is_same_node) {
			_nlist.push_back(new Node<double>(stod(tmpx[0]), stod(tmpx[1]), stod(tmpx[2]), _type, _nlist.size()));
		}
	}
	fin.close();
	return true;
}


//**********要素をVTKファイルに書き出す**********
void exportvtk(std::vector<Node<double>*> _nlist, std::vector<Element<double>*> _elist, std::string _fname) {
	std::ofstream fout(_fname + ".vtk");

	fout << "# vtk DataFile Version 4.1\n" << _fname << "\nASCII\nDATASET UNSTRUCTURED_GRID\n";
	fout << "POINTS\t" << _nlist.size() << "\tfloat\n";
	for (auto pnode : _nlist) {
		fout << pnode->x << "\t" << pnode->y << "\t" << pnode->z << "\n";
	}
	fout << "CELLS\t" << _elist.size() << "\t" << _elist.size() * 5 << "\n";
	for (auto pelement : _elist) {
		fout << "4\t" << pelement->pnodes[0]->id << "\t" << pelement->pnodes[1]->id << "\t" << pelement->pnodes[2]->id << "\t" << pelement->pnodes[3]->id << "\n";
	}
	fout << "CELL_TYPES\t" << _elist.size() << "\n";
	for (int i = 0; i < _elist.size(); i++) {
		fout << "10\n";
	}

	fout.close();
}


//**********メイン処理**********
int main() {
	std::cout << "**********Delaunay trianguration**********\n";

	//----------節点と要素のリスト----------
	std::vector<Node<double>*> nlist;
	std::vector<Element<double>*> elist;
	
	//----------節点の取り込み----------
	if (!importnode(nlist, MODELNAME + (std::string)"/node.dat", 0)) {	return -1;	}
	bool IsCopynodeExist = importnode(nlist, MODELNAME + (std::string)"/copynode.dat", 1);

	//----------Delaunay三角形分割を実行----------
	clock_t ts = clock();
	MakeSupertetrahedran<double>(nlist, elist);
	MakeRoughMesh<double>(nlist, elist);
	DeleteSupertetrahedran<double>(elist);
	if (IsCopynodeExist) {
		DeleteCreviceElement<double>(elist);
	}
	MakeFineMesh<double>(nlist, elist);
	clock_t te = clock();
	
	//----------分割結果の出力----------
	std::cout << "**********The Result↓**********\n";
	std::cout << "time cost:\t" << (double)(te - ts) / CLOCKS_PER_SEC << "sec.\n";
	std::cout << "node:     \t" << nlist.size() << "\n";
	std::cout << "element:  \t" << elist.size() << "\n";
	std::cout << "Export to:\t" << MODELNAME << "/mesh.vtk\n";
	exportvtk(nlist, elist, MODELNAME + (std::string)"/mesh");

	//----------メモリ開放----------
	for (auto pelement : elist) {	delete pelement;	}
	elist.clear();
	for (auto pnode : nlist) {	delete pnode;	}
	nlist.clear();
	
	return 0;
}