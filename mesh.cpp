#include "mesh.h"

Mesh::Mesh(const std::string& filePath, const Material& material) {
    init(filePath, material);
}

Mesh::~Mesh() {}

void Mesh::init(const std::string& filePath, const Material& material) {
    std::ifstream file(filePath);
    std::string line;
    while (std::getline(file, line)) {
        // Skip lines that start with #
        if (line.empty() || line[0] == '#')
            continue;

        std::istringstream iss(line);
        Particle ip{};
        iss >> ip.pid
            >> ip.Vol
            >> ip.xp[0] >> ip.xp[1]
            >> ip.vp[0] >> ip.vp[1]
            >> ip.ep[0] >> ip.ep[1] >> ip.ep[2]
            >> ip.bp[0] >> ip.bp[1];

        // calculate particle mass and momentum
        ip.mass(material);
        ip.momentum();

        if (ip.Vol != 0.0)
            particles.push_back(ip);
    }
    // std::cout << particles.size() << std::endl;
    file.close();
}

void Mesh::showInitInfo() const {
    for (const Particle& p : particles) {
        std::cout << "Particle ID: " << p.pid << "\n"
            << "Volume: " << p.Vol << "\n"
            << "Position: (" << p.xp[0] << ", " << p.xp[1] << ")\n"
            << "Velocity: (" << p.vp[0] << ", " << p.vp[1] << ")\n"
            << "Momentum: (" << p.Pp[0] << ", " << p.Pp[1] << ")\n"
            << "Body Force: (" << p.bp[0] << ", " << p.bp[1] << ")\n"
            << "Strain: (" << p.ep[0] << ", " << p.ep[1] << ", " << p.ep[2] << ")\n"
            << std::endl;
    }
}