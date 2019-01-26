//*********************************************************
//Title		:Node.h
//Author	:Tanabe Yuta
//Date		:2019/01/26
//Copyright	:(C)2019 TanabeYuta
//*********************************************************

#pragma once
class Node
{
public:
	Node();
	~Node();
	Node(double, double, double, int);


	double x, y, z;						//座標
	int type;							//節点の種類


	Node operator+(const Node &);		//ベクトルの和を計算
	Node operator-(const Node &);		//ベクトルの差を計算
	Node operator*(const Node &);		//ベクトルの外積を計算
	double operator^(const Node &);		//ベクトルの内積を計算
	Node operator*(double);				//ベクトルのスカラー倍を計算
	Node operator/(double);				//ベクトルのスカラ―商を計算


	double Size();						//ベクトルの大きさを計算
};

