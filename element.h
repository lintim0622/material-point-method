#pragma once

#include <memory>
#include "node.h"

class Element {
public:
	Element();
	~Element();

	std::unique_ptr<Node> n1;
	std::unique_ptr<Node> n2;
	std::unique_ptr<Node> n3;
	std::unique_ptr<Node> n4;
};