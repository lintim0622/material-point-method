#include "element.h"

Element::Element() {
	n1 = std::make_unique<Node>();
	n2 = std::make_unique<Node>();
	n3 = std::make_unique<Node>();
	n4 = std::make_unique<Node>();
}

Element::~Element() {

}