#pragma once
#include "mesh.h"

typedef std::vector<std::unique_ptr<Mesh>> mesh_list;

class Solve {
public:
	Solve(mesh_list& meshs);
	~Solve();

	inline void setSimulationTime(double et) { endTime = et; }

	// solution process
	void particleToNode();
	void nodalSolution();
	void nodeToParticle();
	void updateParticles();
	void resetNode();

private:
	mesh_list _meshs;
	double endTime;
};