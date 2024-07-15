#pragma once

#include <string>

static const double UNITGRID = 0.25;

static const std::string FILENAME{ "PARTICLE INFO A.txt" };

class Material {
public:
	Material(double rho, double K, double G);
	~Material();

	double rho;
	double K;
	double G;
};
