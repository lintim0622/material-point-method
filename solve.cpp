#include <iostream>
#include "solve.h"

Solve::Solve(mesh_list& meshs) {
	for (std::unique_ptr<Mesh>& msh : meshs) {
		_meshs.push_back(std::move(msh));
	}
	endTime = 0.0;
}

Solve::~Solve() {
	// std::cout << "~Solve()" << std::endl;
}

void Solve::particleToNode() {
	// std::cout << "particleToNode()" << std::endl;
}

void Solve::nodalSolution() {
	// std::cout << "nodalSolution()" << std::endl;
}

void Solve::nodeToParticle() {
	// std::cout << "nodeToParticle()" << std::endl;
}

void Solve::updateParticles() {
	// std::cout << "updateParticles()" << std::endl;
}