#include "node.h"

// Constructor to initialize member variables
Node::Node() : 
    nid{ 0 },
    mn{ 0.0 },
    xn{},
    vn{},
    pn{},
    fint{},
    fext{},
    bn{} {
    // Initialization list is used to set initial values for member variables
}

// Destructor to clean up resources
Node::~Node() {
    // No dynamic memory to clean up for now
}
