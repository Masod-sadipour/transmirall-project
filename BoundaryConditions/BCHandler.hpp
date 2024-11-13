#ifndef BCDATA_HPP
#define BCDATA_HPP

#include <string>
#include <iostream>
#include <libmesh/getpot.h>

#ifndef BC_TYPE_HPP
#define BC_TYPE_HPP

enum class BCType {
    Dirichlet,
    Neumann,
    Penalty,
    NitscheSymmetric,
    NitscheUnsymmetric,
    NodalDirichlet
};

#endif // BC_TYPE_HPP

class BCHandler {
public:
    // Method declarations
    void readBC(const libMesh::GetPot& data, const std::string& section);
    void addDirichletBC(int flag, double value);
    void showMe() const;

private:
    std::string endocardium_bc;
    std::string epicardium_bc;
};

#endif // BCDATA_HPP
