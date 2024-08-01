#include <chrono>

#include "mesh.h"
#include "solve.h"

static std::vector<Boundary> bcSet();

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

    // creat boundary object
    std::vector<Boundary> bcArray{ bcSet() };

    // creat solver object
    Solve sol{ PARTICLEFILE, NODEFILE, elastic };
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
        sol.data_output("particle_output", "node_output", true);

        // reset nodal value
        sol.resetNode();

        if (step % print_interval == 0)
            std::cout << "t = " << t << std::endl;

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

static std::vector<Boundary> bcSet()
{
    std::vector<Boundary> bcArray;

    double bcvalue = 1.5;
    Boundary bc1{ "slip", Vector2D(-bcvalue, -bcvalue), Vector2D(bcvalue,  -bcvalue) };
    Boundary bc2{ "slip", Vector2D(bcvalue,  -bcvalue), Vector2D(bcvalue,   bcvalue) };
    Boundary bc3{ "slip", Vector2D(bcvalue,   bcvalue), Vector2D(-bcvalue,  bcvalue) };
    Boundary bc4{ "slip", Vector2D(-bcvalue,  bcvalue), Vector2D(-bcvalue, -bcvalue) };

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