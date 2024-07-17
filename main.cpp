#include "mesh.h"

int main() {

    // creat material object
    Material elastic{7.85e+3, 139e+9, 75e+9};

    // creat mesh object
    Mesh msh(PARTICLEFILE, NODEFILE, elastic);
    // msh.showParticleInitInfo();
    // msh.showNodeInitInfo();
    msh.showElementInfo();

    std::cout << "end" << std::endl;
    return 0;
}