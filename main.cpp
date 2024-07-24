#include <chrono>

#include "mesh.h"
#include "solve.h"

std::vector<Boundary> bcSet();

// define decay function
std::function<double(double)> decayFunction = [](double x)
{
    return (1.0 - x / UNITGRID) * (1.0 - x / UNITGRID);
};

int main() {

    // Get the start time
    auto start = std::chrono::high_resolution_clock::now();

    // creat material object
    Material elastic{ RHO, K, G };
    elastic.verify_time_step(DT);

    // creat mesh object
    Mesh mshA(PARTICLEFILEA, NODEFILE, elastic);
    Mesh mshB(PARTICLEFILEB, NODEFILE, elastic);
    mesh_list mshs;
    mshs.push_back(std::make_unique<Mesh>(std::move(mshA)));
    mshs.push_back(std::make_unique<Mesh>(std::move(mshB)));

    // creat boundary object
    std::vector<Boundary> bcArray{ bcSet() };

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
        sol.algorithm(t, bcArray, decayFunction);

        // output
        sol.data_output("particle_output.txt", "node_output.txt", true);

        // reset nodal value
        sol.resetNode();

        if (step % print_interval == 0) {
            std::cout << "t = " << t << std::endl;
        }
        // std::cout << "now time: " << t << std::endl;
        t += DT;
        step++;
    }

    // Get the end time
    auto end = std::chrono::high_resolution_clock::now();

    // Calculate the execution time (in seconds)
    std::chrono::duration<double> duration = end - start;

    // Output the execution time
    std::cout << "Execution time: " << duration.count() << " seconds\n";
    return 0;
}

std::vector<Boundary> bcSet() 
{
    std::vector<Boundary> bcArray;

    Boundary bc1{ "slip", Vector2D(-1.5, -1.5), Vector2D(1.5,  -1.5) };
    Boundary bc2{ "slip", Vector2D(1.5,  -1.5), Vector2D(1.5,   1.5) };
    Boundary bc3{ "slip", Vector2D(1.5,   1.5), Vector2D(-1.5,  1.5) };
    Boundary bc4{ "slip", Vector2D(-1.5,  1.5), Vector2D(-1.5, -1.5) };

    bc1.setNbc(Vector2D(0.0,  1.0));
    bc2.setNbc(Vector2D(-1.0, 0.0));
    bc3.setNbc(Vector2D(0.0, -1.0));
    bc4.setNbc(Vector2D(1.0,  0.0));

    bc1.setMu(0.0);
    bc2.setMu(0.0);
    bc3.setMu(0.0);
    bc4.setMu(0.0);

    bcArray.push_back(bc1);
    bcArray.push_back(bc2);
    bcArray.push_back(bc3);
    bcArray.push_back(bc4);

    return bcArray;
}