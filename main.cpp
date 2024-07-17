#include "mesh.h"

int main() {

    // creat material object
    Material elastic{7.85e+3, 139e+9, 75e+9};

    // creat mesh object
    Mesh msh(PARTICLEFILE, NODEFILE, elastic);
    // msh.showParticleInitInfo();
    // msh.showNodeInitInfo();
    // msh.showElementInfo();

    for (auto it : msh.particleElementMap) {
        std::cout << it.first << ", " << it.second << std::endl;
    }
    std::cout << msh.elements[4].n1->nid << std::endl;
    std::cout << msh.elements[4].n2->nid << std::endl;
    std::cout << msh.elements[4].n3->nid << std::endl;
    std::cout << msh.elements[4].n4->nid << std::endl;

    std::cout << "end" << std::endl;
    return 0;
}