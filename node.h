#pragma once
#include "algebra.h"

class Node {
public:
	Node();
	~Node();

	int nid; // id
	double mn; // mass
	Vector2D xn; // position
	Vector2D vn; // velocity
	Vector2D pn; // momentum
	Vector2D fint; // internal force
	Vector2D fext; // external force
	Vector2D bn; // body force
};