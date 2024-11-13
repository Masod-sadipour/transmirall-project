#include <iostream>
#include "libmesh/libmesh.h"
#include "libmesh/mesh.h"
#include "libmesh/equation_systems.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/getpot.h"
#include "Poisson.hpp"


using namespace libMesh;
using namespace BeatIt;

int main(int argc, char** argv) {
    // Initialize libMesh
    LibMeshInit init(argc, argv);

    // Check command-line arguments
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " -f config.in -m mesh.e" << std::endl;
        return 1;
    }

    // Parse command-line arguments
    std::string config_file;
    std::string mesh_file;
    for (int i = 1; i < argc; ++i) {
        if (std::string(argv[i]) == "-f") {
            config_file = argv[++i];
        } else if (std::string(argv[i]) == "-m") {
            mesh_file = argv[++i];
        }
    }

    // Load the mesh
    Mesh mesh(init.comm());
    ExodusII_IO exo_io(mesh);
    exo_io.read(mesh_file);

    // Set up Equation Systems
    EquationSystems equation_systems(mesh);
    libMesh::GetPot data("input/config.in");

    // Create and set up the Poisson solver
    Poisson poisson_solver(equation_systems, "poisson");
    poisson_solver.setup(data);

    // Assemble the system
    poisson_solver.assemble_system();

    // Solve the system
    poisson_solver.solve_system();

    // Export the results
    poisson_solver.save_exo("poisson_output.exo");

    std::cout << "Simulation completed successfully!" << std::endl;
    return 0;
}
