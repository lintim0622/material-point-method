#pragma once
#include <stdexcept>
#include <iomanip> // for std::setprecision

#include "mesh.h"

typedef std::vector<std::unique_ptr<Mesh>> mesh_list;

class Solve {
public:
	Solve(mesh_list& meshs);
	~Solve();

	inline void setSimulationTime(double et) { endTime = et; }

	// solution process
	void algorithm(double nowTime);
	

private:
	void calculateParticleInfo(std::unique_ptr<Mesh>& msh);
	void particleToNode(std::unique_ptr<Mesh>& msh);
	void nodalSolution(std::unique_ptr<Mesh>& msh);
	void nodeToParticle(std::unique_ptr<Mesh>& msh, double nowTime);
	void updateParticles(std::unique_ptr<Mesh>& msh);
	void resetNode(std::unique_ptr<Mesh>& msh);

	mesh_list _meshs;
	double endTime;
};

class Interpolate_Error : public std::runtime_error {
public:
	explicit Interpolate_Error(const std::string& message)
		: std::runtime_error(message) {}
};