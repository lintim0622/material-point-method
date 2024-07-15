#include "constant.h"
#include "mesh.h"

int main() {

    Material elastic{7.85e+3, 139e+9, 75e+9};

    Mesh msh(FILENAME, elastic);
    msh.showInitInfo();
    
    std::cout << "end" << std::endl;
    return 0;
}