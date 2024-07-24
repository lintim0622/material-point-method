#pragma once
#include <stdexcept>
#include <iomanip> // for std::setprecision
#include <functional> // for std::function

#include "mesh.h"

typedef std::vector<std::unique_ptr<Mesh>> mesh_list;

class Boundary;
class Solve {
public:
	Solve(mesh_list& meshs);
	~Solve();

	inline void setSimulationTime(double et) { endTime = et; }

	// solution process
	void algorithm(double nowTime, std::vector<Boundary>& bcArray,
				   const std::function<double(double)>& decayFunction);

	// frame boundary
	void frameBoundary(std::unique_ptr<Mesh>& msh,
					   std::vector<Boundary>& bcArray,
					   const std::function<double(double)>& decayFunction);

	// contact algorithm
	void contact();

	// output solution information
	void data_output(const std::string& pfile_name, const std::string& nfile_name, bool append) const;

	// reset nodal information
	void resetNode();

private:
	void calculateParticleInfo(std::unique_ptr<Mesh>& msh);
	void particleToNode(std::unique_ptr<Mesh>& msh);
	void nodalSolution(std::unique_ptr<Mesh>& msh);
	void nodeToParticle(std::unique_ptr<Mesh>& msh, double nowTime);
	void updateParticles(std::unique_ptr<Mesh>& msh);
	

	mesh_list _meshs;
	double endTime;
};

class Interpolate_Error : public std::runtime_error {
public:
	explicit Interpolate_Error(const std::string& message)
		: std::runtime_error(message) {}
};

// define Global Base Vector
static const Vector2D ex{ 1.0, 0.0 };
static const Vector2D ey{ 0.0, 1.0 };

class Boundary {
public:
	Boundary() = delete;
	Boundary(const std::string& form, const Vector2D& p1, const Vector2D& p2);
	~Boundary();

	// Setters
	void setForm(const std::string& form);
	void setMu(double mu);
	void setP1(const Vector2D& p1);
	void setP2(const Vector2D& p2);
	void setNbc(const Vector2D& nbc);

	// Getters (return by reference)
	const std::string& form() const;
	double mu() const;
	const Vector2D& p1() const;
	const Vector2D& p2() const;
	const Vector2D& nbc() const;

	// find projection point
	Vector2D transformGlobal(const Vector2D& localPosition);
	Vector2D transformLocal(const Vector2D& globalPosition);
	Vector2D newtonMethod(const Vector2D& globalPosition);

private:
	std::string _form;
	double _mu;
	Vector2D _p1;
	Vector2D _p2;
	Vector2D _nbc;
	Vector2D _exl;
	Vector2D _eyl;

	// Global to Local
	std::vector<Vector2D> _Txlx;

	// Local to Global
	std::vector<Vector2D> _Txxl;
};

// define layer function
static double layer(double r, const std::function<double(double)>& decayFunction);

static double norm(const Vector2D& vector);

void modify_normal(Node& node, Node& other_node, Vector2D& nB);