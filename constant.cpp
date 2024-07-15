#include "constant.h"

Material::Material(double rho, double K, double G) :
    rho{ rho }, K{ K }, G{ G } {}

Material::~Material() {}