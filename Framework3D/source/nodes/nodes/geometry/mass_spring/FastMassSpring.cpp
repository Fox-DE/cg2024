#include "FastMassSpring.h"
#include <iostream>


namespace USTC_CG::node_mass_spring {
FastMassSpring::FastMassSpring(const Eigen::MatrixXd& X, const EdgeSet& E, const float stiffness, const float h): 
MassSpring(X, E){
    // construct L and J at initialization
    std::cout << "init fast mass spring" << std::endl;

    unsigned n_vertices = X.rows();
    this->stiffness = stiffness; 
    this->h = h; 

    M_h2L.resize(3 * n_vertices, 3 * n_vertices);

    // (HW Optional) precompute A and prefactorize
    // Note: one thing to take care of: A is related with stiffness, if stiffness changes, A need to be recomputed
    std::vector<Triplet<double>> tripletlist;
    for (const auto& e : E) {
        for (int j = 0; j < 3; j++)
        {          
                tripletlist.push_back(Triplet<double>(3 * e.first + j, 3 * e.first + j, h * h * stiffness));
                tripletlist.push_back(Triplet<double>(3 * e.first + j, 3 * e.second + j, -h * h * stiffness));            
                tripletlist.push_back(
                    Triplet<double>(3 * e.second + j, 3 * e.first + j, -h * h * stiffness));
                tripletlist.push_back(
                    Triplet<double>(3 * e.second + j, 3 * e.second + j, h * h * stiffness));           
        }
    }
    double mass_per_vertex = mass/n_vertices;
    for (size_t i = 0; i < 3 * n_vertices; i++) {
        tripletlist.push_back(Triplet<double>(i, i, mass_per_vertex));
    }
    M_h2L.setFromTriplets(tripletlist.begin(), tripletlist.end());
    M_h2L.makeCompressed();
    //Found J
    J.resize(3 * n_vertices, 3 * E.size());
    std::vector<Triplet<double>> tripletlistJ;
    unsigned i = 0;
    for (const auto& e : E) {
        for (int j = 0; j < 3; j++) {
                tripletlistJ.push_back(Triplet<double>(3 * e.first + j, 3 * i + j, stiffness));
                tripletlistJ.push_back(Triplet<double>(3 * e.second + j, 3 * i + j, -stiffness));
        }
        i++;
    }
    J.setFromTriplets(tripletlistJ.begin(), tripletlistJ.end());

    //Found K and b
    b.resize(3 * n_vertices, 1);
    b.setZero();
    std::vector<Triplet<double>> tripletlistK;
    int m_size = 0;
    for (int fix = 0; fix < n_vertices;fix++) {
        if (!dirichlet_bc_mask[fix])
        {
                m_size++;
        }
    }
    K.resize(3 * m_size, 3 * n_vertices);
    int j = 0;
    for (size_t i = 0; i < n_vertices; i++) {
        if (!dirichlet_bc_mask[i]) {
                for (int k = 0; k < 3; k++) {
                    tripletlistK.push_back(Triplet<double>(3 * j + k, 3 * i + k, 1.0));
                }
                j++;
        }
        else {
                for (int k = 0; k < 3; k++) {
                    b(3 * i + k, 0) = X(i,k);
                }
        }
    }
    K.setFromTriplets(tripletlistK.begin(), tripletlistK.end());
    Kt = K.transpose();
    //Found A
    Eigen::SparseMatrix<double> A;
    if (n_vertices-m_size>0)
    {
        A = K * M_h2L * Kt;
    }
    else
    {
        A = M_h2L;
    }
    A.makeCompressed();
    A_LU.compute(A);
    y_acc.resize(3 * n_vertices, 1);
    d.resize(3 * E.size(), 1);
    x_nacc.resize(3 * n_vertices, 1);
}

void FastMassSpring::step()
{
    TIC(acc_implicit)
    // (HW Optional) Necessary preparation
    // ...
    Eigen::Vector3d acceleration_ext = gravity + wind_ext_acc;
    unsigned n_vertices = X.rows();
    double mass_per_vertex = mass/n_vertices ;
    //get x_na and y
    for (size_t i = 0; i < n_vertices; i++) {
        for (int j = 0; j < 3; j++) {
                x_nacc(3 * i + j, 0) = X(i,j);
                y_acc(3 * i + j, 0) = X(i, j) + h * vel(i,j) + h * h* acceleration_ext[j];
        }
    }
    //Get d
    unsigned spring = 0;
    for (const auto& e : E) {
        auto p1 = X.row(e.first);
        auto p2 = X.row(e.second);
        auto norm = (p1 - p2).norm();
        auto l = E_rest_length[spring];
        for (int j = 0; j < 3; j++) {
                d(3 * spring + j, 0) = l / norm * (p1[j] - p2[j]);
        }
        spring++;
    }

    for (unsigned iter = 0; iter < max_iter; iter++) {
        // (HW Optional)
        // local_step and global_step alternating solving
        //std::cout << K << std::endl;
        MatrixXd part1 = h * h * J * d;
        MatrixXd part2 = mass_per_vertex * y_acc;
        MatrixXd part3 = M_h2L * b;
        MatrixXd xr0 = (h * h * J * d + mass_per_vertex * y_acc - M_h2L * b);
        MatrixXd xr = K * xr0;
        MatrixXd xa = A_LU.solve(xr);
        //std::cout << M_h2L << std::endl;
        //std::cout << b << std::endl;

        //update x_nacc
       size_t j = 0;
        for (size_t i = 0; i < n_vertices; i++) {
                if (!dirichlet_bc_mask[i]) {
                    for (int k = 0; k < 3; k++) {
                        x_nacc(3 * i + k, 0) = xa(3 * j + k, 0);
                    }
                    j++;
                }
        }

        //update d for next itr
        unsigned i = 0;
        for (const auto& e : E) {
                Vector3d p1, p2;
                for (int k = 0; k < 3; k++) {
                    p1[k] = x_nacc(3 * e.first + k, 0);
                    p2[k] = x_nacc(3 * e.second + k, 0);
                }
                double dist = (p1 - p2).norm();
                double l = E_rest_length[i];
                for (int j = 0; j < 3; j++) {
                    d(3 * i + j, 0) = l / dist * (p1[j] - p2[j]);
                }
                i++;
        }
    }

	for (int i = 0; i < n_vertices; i++) {
        if (!dirichlet_bc_mask[i]) {
                for (int j = 0; j < 3; j++) {
                    vel(i,j) = damping*(x_nacc(3 * i + j, 0) - X(i,j)) / h;
                    X(i,j) = x_nacc(3 * i + j, 0);
                }
        }
    }
    TOC(acc_implicit)
}

}  // namespace USTC_CG::node_mass_spring