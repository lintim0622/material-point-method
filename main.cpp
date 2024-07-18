#include "mesh.h"
#include "solve.h"

int main() {

    // creat material object
    Material elastic{7.85e+3, 70e+9, 5e+7};
    elastic.verify_time_step(1e-4);

    // creat mesh object
    Mesh msh(PARTICLEFILE, NODEFILE, elastic);
    for (auto it : msh.pem) {
        std::cout << it.first << ", " << it.second << std::endl;
    }
    std::cout << msh.elements[4].n1->nid << std::endl;
    std::cout << msh.elements[4].n2->nid << std::endl;
    std::cout << msh.elements[4].n3->nid << std::endl;
    std::cout << msh.elements[4].n4->nid << std::endl;

    mesh_list mshs;
    mshs.push_back(std::make_unique<Mesh>(std::move(msh)));

    Solve sol{ mshs };

    std::cout << "end" << std::endl;
    return 0;
}