#include "mesh.h"
#include "solve.h"

int main() {

    // creat material object
    Material elastic{ RHO, K, G };
    elastic.verify_time_step(DT);

    // creat mesh object
    Mesh msh(PARTICLEFILE, NODEFILE, elastic);
    mesh_list mshs;
    mshs.push_back(std::make_unique<Mesh>(std::move(msh)));

    // creat solver object
    Solve sol{ mshs };
    sol.setSimulationTime(ENDTIME);

    // run
    double t = 0.0;
    int step = 0;
    int print_interval = static_cast<int>(round(1 / DT) / 10);
    while (t <= ENDTIME)
    {
        // main process
        sol.algorithm(t);

        if (step % print_interval == 0) {
            std::cout << "t = " << t << std::endl;
        }
        // std::cout << "now time: " << t << std::endl;
        t += DT;
        step++;
    }

    std::cout << "end" << std::endl;
    return 0;
}