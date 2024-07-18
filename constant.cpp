#include "constant.h"

Material::Material(double rho, double K, double G) :
    rho{ rho }, K{ K }, G{ G } 
{
    E = 9.0 * K * G / (3.0 * K + G);
    v = 0.5 * (3.0 * K - 2.0 * G) / (3.0 * K + G); // poisson
}

Material::~Material() {}

void Material::verify_time_step(double dt) {
    // P wave velocity->vel_p
    double M = K + 4.0 * G / 3.0;
    double vel_p = std::sqrt(M / rho);
    double dt_critical = UNITGRID / (10.0 * vel_p);
    std::cout << "MPM critical time step = " << dt_critical << std::endl;
    if (dt > dt_critical) {
        std::cout << "dt = " << dt << std::endl;
        std::cout << "dt_critical = " << dt_critical << std::endl;
        std::cout << "need reset MPM time step" << std::endl;
        // Equivalent to Python's sys.exit()
        std::exit(EXIT_FAILURE); 
    }
}
    