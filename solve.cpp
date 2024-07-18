#include <iostream>
#include "solve.h"

Solve::Solve(mesh_list& meshs) {
	for (std::unique_ptr<Mesh>& msh : meshs) {
		_meshs.push_back(std::move(msh));
	}
	endTime = 0.0;
}

Solve::~Solve() {

}

void Solve::particleToNode() {
	for (std::unique_ptr<Mesh>& msh : _meshs) {
		for (auto it : msh->pem) {
			const double& xo = msh->elements[it.second].n1->xn[0];
			const double& yo = msh->elements[it.second].n1->xn[1];
			static const double& l = UNITGRID, h = UNITGRID;
			// std::cout << xo << ", " << yo << std::endl;
		}
	}
}

void Solve::nodalSolution() {

}

void Solve::nodeToParticle() {

}

void Solve::updateParticles() {

}

void Solve::resetNode() {

}