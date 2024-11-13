/*
 * Poisson.cpp
 *
 * Updated version with fixes for API changes, missing includes, and BCType issues.
 */

#include "Poisson.hpp"
#include <string>
#include <Eigen/Core>
#include <libmesh/getpot.h>
#include "/nas/longleaf/home/masods/transmirall/BoundaryConditions/BCHandler.hpp"

// Basic include files needed for the mesh functionality.
#include "libmesh/mesh.h"
#include "libmesh/type_tensor.h"

// Include files for linear systems and FE objects
#include "libmesh/linear_implicit_system.h"
#include "libmesh/vector_value.h"
#include "libmesh/linear_solver.h"
#include "libmesh/fe.h"
#include "libmesh/quadrature_gauss.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_vector.h"
#include "libmesh/elem.h"
#include "libmesh/dof_map.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/enum_preconditioner_type.h"
#include "libmesh/enum_solver_type.h"
#include "libmesh/perf_log.h"
#include "libmesh/dirichlet_boundaries.h"

#include "/nas/longleaf/home/masods/transmirall/Util/IO/io.hpp"

namespace BeatIt {

typedef libMesh::LinearImplicitSystem PoissonSystem;

Poisson::Poisson(libMesh::EquationSystems& es, std::string system_name)
    : M_datafile(),
      M_equationSystems(es),
      M_bch(),
      M_parametersExporterNames(),
      M_exporter(),
      M_linearSolver(),
      M_outputFolder(),
      M_rhsFunction(),
      M_myName(system_name),
      M_myNameGradient(system_name + "_gradient"),
      M_myNameP0(system_name + "_P0"),
      M_active_subdomains(),
      M_D(1.0) {}

Poisson::~Poisson() {}

void Poisson::setup(const libMesh::GetPot& data, std::string section) {
    M_datafile = data;
    std::string output_folder = M_datafile(section + "/output_folder", "Output");
    M_outputFolder = "./" + output_folder + "/";
    std::cout << "* POISSON: Output folder: " << M_outputFolder << std::endl;

    PoissonSystem& system = M_equationSystems.add_system<PoissonSystem>(M_myName);

    if (!system.is_initialized()) {
        std::string blockIDs_list = data.get(section + "/blockIDs", "666");
        BeatIt::readList(blockIDs_list, M_active_subdomains);
        if (*M_active_subdomains.begin() != 666)
            system.add_variable(M_myName + "_phi", libMesh::FIRST, libMesh::LAGRANGE, &M_active_subdomains);
        else {
            M_active_subdomains.clear();
            system.add_variable(M_myName + "_phi", libMesh::FIRST);
        }
    }

    std::cout << "* POISSON: Setup Homogeneous Dirichlet BC ... " << std::flush;
    auto* dirichlet_boundaries = system.get_dof_map().get_dirichlet_boundaries();
    if (dirichlet_boundaries) {
        for (auto& d : *dirichlet_boundaries) {
            system.get_dof_map().remove_dirichlet_boundary(*d);
        }
    }

    M_bch.clear();
    M_bch.readBC(data, section);
    M_bch.showMe();

    for (auto&& bc_ptr : M_bch.M_bcs) {
        auto bc_type = bc_ptr->get_type();
        if (bc_type == BCType::Dirichlet) {
            std::set<libMesh::boundary_id_type> dirichlet_boundary_ids;
            std::vector<unsigned int> variables{0};

            auto num_flags = bc_ptr->size();
            for (int nflag = 0; nflag < num_flags; ++nflag) {
                dirichlet_boundary_ids.insert(bc_ptr->get_flag(nflag));
            }

            libMesh::DirichletBoundary dirichlet_bc(std::move(dirichlet_boundary_ids), std::move(variables), &(bc_ptr->get_function()));
            system.get_dof_map().add_dirichlet_boundary(std::move(dirichlet_bc));
        }
    }

    if (!system.is_initialized()) {
        system.init();
    } else {
        system.reinit();
    }

    M_linearSolver.reset(new libMesh::PetscLinearSolver<libMesh::Number>(M_equationSystems.comm()));
    M_linearSolver->set_solver_type(libMesh::CG);
    M_linearSolver->set_preconditioner_type(libMesh::AMG_PRECOND);
    M_linearSolver->init();
    KSPSetOptionsPrefix(M_linearSolver->ksp(), "poisson_");
    PCSetOptionsPrefix(M_linearSolver->pc(), "poisson_");
    KSPSetFromOptions(M_linearSolver->ksp());

    M_exporter.reset(new EXOExporter(M_equationSystems.get_mesh()));
}

void Poisson::apply_BC(const libMesh::Elem*& elem,
                       // libMesh::DenseMatrix<libMesh::Number>& Ke,
                       libMesh::DenseVector<libMesh::Number>& Fe,
                       std::unique_ptr<libMesh::FEBase>& fe_face,
                       libMesh::QGauss& qface,
                       const libMesh::MeshBase& mesh) {
    for (unsigned int side = 0; side < elem->n_sides(); ++side) {
        if (!elem->neighbor_ptr(side)) {
            std::vector<short int> boundary_ids;
	    mesh.get_boundary_info().boundary_ids(elem, side, boundary_ids);
	    const unsigned int boundary_id = boundary_ids.empty() ? 0 : boundary_ids[0];
            auto bc = M_bch.get_bc(boundary_id);

            if (bc) {
                fe_face->reinit(elem, side);
                const auto& JxW_face = fe_face->get_JxW();
                const auto& phi_face = fe_face->get_phi();
                const auto& qface_point = fe_face->get_xyz();
                //const auto& normals = fe_face->get_normals();
		const auto& qfat = qface;
                auto bc_type = bc->get_type();

                if (bc_type == BCType::Neumann) {
                    for (unsigned int qp = 0; qp < qface.n_points(); ++qp) {
                        const double value = bc->get_function()(0.0, qface_point , qfat , 0);
                 for (unsigned int i = 0; i < phi_face.size(); ++i) {
                            Fe(i) += JxW_face[qp] * value * phi_face[i][qp];
                        }
                    }
                }
            }
        }
    }
}

}  // namespace BeatIt
