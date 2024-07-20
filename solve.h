#pragma once
#include "mesh.h"

typedef std::vector<std::unique_ptr<Mesh>> mesh_list;

class Solve {
public:
	Solve(mesh_list& meshs);
	~Solve();

	inline void setSimulationTime(double et) { endTime = et; }

	// solution process
	void algorithm();
	

private:
	void calculateParticleInfo(std::unique_ptr<Mesh>& msh);
	void particleToNode(std::unique_ptr<Mesh>& msh);
	void nodalSolution(std::unique_ptr<Mesh>& msh);
	void nodeToParticle(std::unique_ptr<Mesh>& msh);
	void updateParticles(std::unique_ptr<Mesh>& msh);
	void resetNode(std::unique_ptr<Mesh>& msh);

	mesh_list _meshs;
	double endTime;
};