#pragma once

#include <iostream>
#include <string>
#include <cmath>
#include <cstdlib>

// Grid Size (Boundary)
static const double UNITGRID = 0.25;
static const double MAXBC = 1.5;
static const double MINBC = -1.5;

// Simulation Time
static const double DT = 1e-4;
static const double ENDTIME = 2e-4;

// material property
static const double RHO = 960.0;
static const double K = 0.003E+9;
static const double G = 0.0006E+9;

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

	double E1;
	double E2;
};

