//*********************************************************
//Title		:Element.h
//Author	:Tanabe Yuta
//Date		:2019/01/26
//Copyright	:(C)2019 TanabeYuta
//*********************************************************


#pragma once
#include <array>

#include "Node.h"

class Surface;

class Element
{
public:
	Element();
	~Element();
	Element(Node*, Node*, Node*, Node*);


	bool IsActive;
	std::array<Surface*, 4> psurfaces;		//四面体の表面
	std::array<Node*, 4> pnodes;			//四面体を構成する頂点
	Node scenter;							//外接球の中心座標
	double sround;							//外接球の半径
	Node gcenter;							//四面体の重心座標
	double volume;							//四面体の体積
	double aspect;							//四面体のアスペクト比


	Element* GetLocateId(Node*);			//要素内に点があれば自身を指すポインタを返す　そうでなければ隣接要素を指すポインタを返す
	bool IsInSphere(Node*);					//true:点が外接球内
	Surface* GetAdjacentSurface(Element*);	//渡された要素と隣接する面を指すポインタを返す
};

