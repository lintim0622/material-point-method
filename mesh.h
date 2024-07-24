#pragma once

#include <fstream>
#include <vector>
#include <sstream>
#include <algorithm>
#include <memory>
#include <map>

#include "algebra.h"
#include "constant.h"

static const double dxi{ 2.0 / UNITGRID };
static const double deta{ 2.0 / UNITGRID };
static const double N_min{ 1e-3 };
static const double L_cut{ 1e-2 };

class Material {
public:
	Material(double rho, double K, double G);
	~Material();

	void verify_time_step(double dt);

	double rho;
	double K;
	double G;
	double E;
	double v;

	double E1;
	double E2;
};

class Particle {
public:
	Particle();
	~Particle();

	void calculateMass(const Material& material);
	void calculateMomentum();
	void calculateSpecificStress(const Material& material);
	void show() const;

	int pid; // id
	double Vol; // volume
	double mp; // mass
	Vector2D xp; // position
	Vector2D vp; // velocity
	Vector2D Pp; // momentum
	Vector2D bp; // body force
	double ep[3]; // strain
	double ssp[3]; // specific stress
	double ineps[3]; // Increment strain tensor

	// shape function (QUAD4)
	double N1, N2, N3, N4;
	Vector2D dN1, dN2, dN3, dN4;
};

class Node {
public:
	Node();
	~Node();

	int nid; // id
	double mn; // mass
	Vector2D xn; // position
	Vector2D xcn; // centroid position
	Vector2D vn; // velocity
	Vector2D an; // acceleration
	Vector2D pn; // momentum
	Vector2D fint; // internal force
	Vector2D fext; // external force
	Vector2D ftot; // external force
	Vector2D fbc; // boundary force
	Vector2D fct; // contact force
	Vector2D bn; // body force
	Vector2D normal; // normal vector
};

class Element {
public:
	Element();
	~Element();

	Element(const Element&) = delete;
	Element& operator=(const Element&) = delete;

	Element(Element&&) noexcept = default;
	Element& operator=(Element&&) noexcept = default;

	bool contains(const Particle& particle) const;

	Node* n1;
	Node* n2;
	Node* n3;
	Node* n4;
};

class Mesh {
public:
	Mesh(const std::string& particleFile, const std::string& nodeFile, const Material& material);

	Mesh(const Mesh&) = delete;
	Mesh& operator=(const Mesh&) = delete;

	Mesh(Mesh&&) noexcept = default;
	Mesh& operator=(Mesh&&) noexcept = default;

	~Mesh();

	void showParticleInitInfo() const;
	void showNodeInitInfo() const;
	void showElementInfo() const;

	void createElementParticleMap();
	void clearMapId();
	// Element* findElementForParticle(const Particle& particle);

	std::vector<Particle> particles;
	std::vector<Node> nodes;
	std::vector<Element> elements;
	const Material material;

	// Map particle index to element index
	std::map<int, int> pem; 

private:
	// read initial particle information and calculate
	void initParticleInfo(const std::string& filePath, const Material& material);

	// read initial Grid information
	void initNodeInfo(const std::string& filename);

	// initial mesh
	void createElements();
};


// shape function (QUAD4)
static double stable(double fij);
static double  stabledN(double dfij, double fij);

double shapeN1(const double xi, const double eta);
double shapeN2(const double xi, const double eta);
double shapeN3(const double xi, const double eta);
double shapeN4(const double xi, const double eta);
double shapedN1x(const double eta);
double shapedN1y(const double xi);
double shapedN2x(const double eta);
double shapedN2y(const double xi);
double shapedN3x(const double eta);
double shapedN3y(const double xi);
double shapedN4x(const double eta);
double shapedN4y(const double xi);
