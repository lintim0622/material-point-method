#pragma once

class Solve {
public:
	Solve();
	~Solve();

	// solution process
	void particleToNode();
	void nodalSolution();
	void nodeToParticle();
	void updateParticles();

};