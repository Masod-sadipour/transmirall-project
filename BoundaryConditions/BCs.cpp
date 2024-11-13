#include "BCHandler.hpp"

// Implementation of BCHandler methods
void BCHandler::readBC(const GetPot& data, const std::string& section) {
    endocardium_bc = data(section + "/bc_endocardium", "");
    epicardium_bc = data(section + "/bc_epicardium", "");

    if (!endocardium_bc.empty()) {
        std::cout << "Endocardium BC: " << endocardium_bc << std::endl;
    }

    if (!epicardium_bc.empty()) {
        std::cout << "Epicardium BC: " << epicardium_bc << std::endl;
    }
}

void BCHandler::addDirichletBC(int flag, double value) {
    std::cout << "Setting Dirichlet BC on sideset flag " << flag
              << " with value " << value << std::endl;
}

void BCHandler::showMe() const {
    std::cout << "Endocardium BC: " << endocardium_bc << std::endl;
    std::cout << "Epicardium BC: " << epicardium_bc << std::endl;
}
