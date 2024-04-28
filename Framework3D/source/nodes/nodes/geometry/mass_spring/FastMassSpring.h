#pragma once 
#include "MassSpring.h"
#include <memory>

namespace USTC_CG::node_mass_spring {
// Impliment the Liu13's paper: https://tiantianliu.cn/papers/liu13fast/liu13fast.pdf
class FastMassSpring : public MassSpring {
   public:
    FastMassSpring() = default;
    ~FastMassSpring() = default; 

    FastMassSpring(const Eigen::MatrixXd& X, const EdgeSet& E, const float stiffness);
    void step() override;
    unsigned max_iter = 100; // (HW Optional) add UI for this parameter

   protected:
    // Custom variables, like prefactorized A 
    Eigen::SparseLU<Eigen::SparseMatrix<double>> A_LU;
    Eigen::SparseMatrix<double> J;
};
}  // namespace USTC_CG::node_mass_spring