#pragma once

#include <fstream>
#include <vector>
#include <sstream>
#include <algorithm>
#include <memory>
#include <map>

#include "algebra.h"
#include "constant.h"

class Particle {
public:
	Particle();
	~Particle();

	void mass(const Material& material);
	void momentum();
	void show() const;

	int pid; // id
	double Vol; // volume
	double mp; // mass
	Vector2D xp; // position
	Vector2D vp; // velocity
	Vector2D Pp; // momentum
	Vector2D bp; // body force
	double ep[3]; // strain
	double sp[3]; // stress
};

class Node {
public:
	Node();
	~Node();

	int nid; // id
	double mn; // mass
	Vector2D xn; // position
	Vector2D vn; // velocity
	Vector2D pn; // momentum
	Vector2D fint; // internal force
	Vector2D fext; // external force
	Vector2D bn; // body force
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

	std::unique_ptr<Node> n1;
	std::unique_ptr<Node> n2;
	std::unique_ptr<Node> n3;
	std::unique_ptr<Node> n4;
};

class Mesh {
public:
	Mesh(const std::string& particleFile, const std::string& nodeFile, const Material& material);
	~Mesh();

	void showParticleInitInfo() const;
	void showNodeInitInfo() const;
	void showElementInfo() const;

	void createElementParticleMap();
	Element* findElementForParticle(const Particle& particle);

	std::vector<Particle> particles;
	std::vector<Node> nodes;
	std::vector<Element> elements;
	std::map<int, int> particleElementMap; // Map particle index to element index

private:
	// read initial particle information and calculate
	void initParticleInfo(const std::string& filePath, const Material& material);

	// read initial Grid information
	void initNodeInfo(const std::string& filename);

	// initial mesh
	void createElements();
};