#pragma once

#include <iostream>
#include <string>
#include <cmath>
#include <cstdlib>

static const double UNITGRID = 0.25;

static const std::string PARTICLEFILE{ "PARTICLE INFO A.txt" };
static const std::string NODEFILE{ "NODAL POSITION.txt.txt" };

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
};
