#include "particle.h"

// Constructor to initialize member variables
Particle::Particle() :
    pid{ 0 },
    Vol{ 0.0 },
    mp{ 0.0 },
    xp{},
    vp{},
    Pp{},
    bp{},
    ep{ 0.0, 0.0, 0.0 },
    sp{ 0.0, 0.0, 0.0 } {
    // Initialization list is used to set initial values for member variables
    // std::cout << "Particle()" << std::endl;
}

// Destructor to clean up resources
Particle::~Particle() {
    // std::cout << "~Particle()" << std::endl;
}

void Particle::mass(const Material& material) {
    mp = material.rho * Vol;
}

void Particle::momentum() {
    Pp = mp * vp;
    // std::cout << "momentum()" << std::endl;
}

void Particle::show() const {
    std::cout << "id : " << pid << std::endl;
    std::cout << "Volume : " << Vol << std::endl;
    std::cout << "mass : " << mp << std::endl;
    std::cout << "position : " << xp << std::endl;
    std::cout << "velocity : " << vp << std::endl;
    std::cout << "momentum : " << Pp << std::endl;
    std::cout << "stress : " << sp[0]  << ", " << sp[1] << ", " << sp[2] << std::endl;
}