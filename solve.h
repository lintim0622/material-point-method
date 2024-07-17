#pragma once
#include "mesh.h"

typedef std::vector<std::unique_ptr<Mesh>> mesh_list;

class Solve {
public:
	Solve(mesh_list& meshs);
	~Solve();

	// solution process
	void particleToNode();
	void nodalSolution();
	void nodeToParticle();
	void updateParticles();

private:
	mesh_list _meshs;
};