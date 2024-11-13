/*
 * Poisson.hpp
 *
 *  Created on: Sep 14, 2016
 *      Author: srossi
 */

#ifndef SRC_POISSONSOLVER_POISSON_HPP_
#define SRC_POISSONSOLVER_POISSON_HPP_

// Include file that defines (possibly multiple) systems of equations.
#include "libmesh/equation_systems.h"
#include "libmesh/linear_solver.h"
#include "libmesh/petsc_linear_solver.h"
#include "libmesh/getpot.h" // Ensure you are using the libMesh version of GetPot

#include <memory>
#include "/nas/longleaf/home/masods/transmirall/BoundaryConditions/BCHandler.hpp"
#include "/nas/longleaf/home/masods/transmirall/Util/SpiritFunction.hpp"

namespace libMesh
{
class Mesh;
class ExodusII_IO;
class VTKIO;
class GMVIO;
class MeshRefinement;
class TimeData;
class ExplicitSystem;
class LinearImplicitSystem;
class PacingProtocol;
class ErrorVector;
template <class T> class DenseMatrix;
template <class T> class DenseVector;
class QGauss;
class MeshBase;
}

namespace BeatIt {

class BCHandler;
class SpiritFunction;

class Poisson {
    typedef libMesh::ExodusII_IO EXOExporter;

public:
    Poisson(libMesh::EquationSystems& es, std::string system_name = "poisson");
    virtual ~Poisson();
    
    void setup(const libMesh::GetPot& data, std::string section = "poisson");
    void assemble_system();
    void solve_system();
    void save_exo(const std::string& output_filename = "poisson.exo");

    const std::unique_ptr<libMesh::NumericVector<libMesh::Number>>& get_gradient();
    const std::unique_ptr<libMesh::NumericVector<libMesh::Number>>& get_solution();
    const std::unique_ptr<libMesh::NumericVector<libMesh::Number>>& get_P0_solution();

    void write_equation_system(const std::string& es = "poisson.dat");
    void read_equation_system(const std::string& es = "poisson.dat");

    double get_solution_norm();
    void compute_elemental_solution_gradient();
    void deleteSystems();

    // Input file
    libMesh::EquationSystems& M_equationSystems;
    BCHandler M_bch;

    std::set<std::string> M_parametersExporterNames;
    std::unique_ptr<EXOExporter> M_exporter;
    std::unique_ptr<libMesh::PetscLinearSolver<libMesh::Number>> M_linearSolver;

    std::string M_outputFolder;
    SpiritFunction M_rhsFunction;
    std::string M_myName;
    std::string M_myNameGradient;
    std::string M_myNameP0;
    std::set<libMesh::subdomain_id_type> M_active_subdomains;
    double M_D;

private:
    libMesh::GetPot M_datafile;
    void apply_BC(const libMesh::Elem*& elem,
                  libMesh::DenseVector<libMesh::Number>& Fe,
                  std::unique_ptr<libMesh::FEBase>& fe_face,
                  libMesh::QGauss& qface,
                  const libMesh::MeshBase& mesh);
};

} /* namespace BeatIt */

#endif /* SRC_POISSONSOLVER_POISSON_HPP_ */
