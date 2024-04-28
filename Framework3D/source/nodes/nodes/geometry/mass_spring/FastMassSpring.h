#pragma once 
#include "MassSpring.h"
#include <memory>

namespace USTC_CG::node_mass_spring {
// Impliment the Liu13's paper: https://tiantianliu.cn/papers/liu13fast/liu13fast.pdf
class FastMassSpring : public MassSpring {
   public:
    FastMassSpring() = default;
    ~FastMassSpring() = default; 

    FastMassSpring(const Eigen::MatrixXd& X, const EdgeSet& E, const float stiffness, const float h);
    void step() override;
    unsigned max_iter = 100; // (HW Optional) add UI for this parameter

   protected:
    // Custom variables, like prefactorized A 
    Eigen::SparseLU<Eigen::SparseMatrix<double>> A_LU;
    Eigen::SparseMatrix<double> J;
    Eigen::SparseMatrix<double> M_h2L;
    Eigen::MatrixXd b;  // x-K*Kt*x
    Eigen::SparseMatrix<double> K;
    Eigen::SparseMatrix<double> Kt;  // K transpose
    Eigen::MatrixXd y_acc;
    Eigen::MatrixXd d;
    Eigen::MatrixXd x_nacc;
};
}  // namespace USTC_CG::node_mass_spring