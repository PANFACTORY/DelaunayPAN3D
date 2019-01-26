//*********************************************************
//Title		:Surface.h
//Author	:Tanabe Yuta
//Date		:2019/01/26
//Copyright	:(C)2019 TanabeYuta
//*********************************************************


#pragma once
#include <array>

#include "Node.h"

class Element;

class Surface
{
public:
	Surface();
	~Surface();
	Surface(Node*, Node*, Node*, Element*, Element*);


	Element* pneighbor;					//隣接要素を指すポインタ
	bool IsActive;						//true:多面体の境界面
	std::array<Node*, 3> pnodes;		//節点
	Element* pparent;					//この面を持つ要素を指すポインタ
	
	
	bool operator==(const Surface &);	//接している面か判定


	bool IsRayCross(Node, Node);		//true:ベクトルが面と交差
};

