#pragma once
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