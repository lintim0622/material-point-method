#pragma once

#include <iostream>
#include <string>
#include <cmath>
#include <cstdlib>

static const double PI = 3.14;

// Grid Size (Boundary)
static const double UNITGRID = 0.25;
static const double MAXBC = 1.5;
static const double MINBC = -1.5;
static const double TOL = 1e-10;

// Simulation Time
static const double DT = 1e-4;
static const double ENDTIME = 0.5;

// material property
static const double RHO = 960.0;
static const double K = 0.003E+9;
static const double G = 0.0006E+9;

static const std::string PARTICLEFILE{ "PARTICLE INFO.txt" };
static const std::string NODEFILE{ "NODAL POSITION.txt" };


