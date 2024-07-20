#include "mesh.h"
#include "solve.h"

int main() {

    // creat material object
    Material elastic{ RHO, K, G };
    elastic.verify_time_step(DT);

    // creat mesh and solver object
    Mesh msh(PARTICLEFILE, NODEFILE, elastic);
    //for (auto it : msh.pem) {
    //    std::cout << it.first << ", " << it.second << std::endl;
    //}

    mesh_list mshs;
    mshs.push_back(std::make_unique<Mesh>(std::move(msh)));

    Solve sol{ mshs };
    sol.setSimulationTime(ENDTIME);

    // run
    double t = 0.0;
    int step = 0;
    int print_interval = static_cast<int>(round(1 / DT) / 10);
    while (t <= ENDTIME) {

        std::cout << "now time: " << t << std::endl;

        // particle information
        sol.calculateParticleInfo();

        // particle to node
        sol.particleToNode();

        // nodal solution
        sol.nodalSolution();

        // node to particle
        sol.nodeToParticle();

        // update particle
        sol.updateParticles();

        // reset node
        sol.resetNode();

      /*  if (step % print_interval == 0) {
            std::cout << "t = " << t << std::endl;
        }*/
        t += DT;
        step++;
    }


    std::cout << "end" << std::endl;
    return 0;
}