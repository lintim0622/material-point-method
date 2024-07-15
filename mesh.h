#pragma once

#include <fstream>
#include <vector>
#include <sstream>

#include "particle.h"
#include "element.h"

class Mesh {
public:
	Mesh(const std::string& filePath, const Material& material);
	~Mesh();

	void showInitInfo() const;

	std::vector<Particle> particles;

private:
	// read initial particle information and calculate
	void init(const std::string& filePath, const Material& material);
};