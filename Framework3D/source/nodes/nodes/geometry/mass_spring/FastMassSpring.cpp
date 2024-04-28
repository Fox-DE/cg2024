#include "FastMassSpring.h"
#include <iostream>


namespace USTC_CG::node_mass_spring {
FastMassSpring::FastMassSpring(const Eigen::MatrixXd& X, const EdgeSet& E, const float stiffness): 
MassSpring(X, E){
    // construct L and J at initialization
    std::cout << "init fast mass spring" << std::endl;

    unsigned n_vertices = X.rows();
    this->stiffness = stiffness; 

    Eigen::SparseMatrix<double> A(n_vertices * 3, n_vertices * 3);
    A.setZero();

    // (HW Optional) precompute A and prefactorize
    // Note: one thing to take care of: A is related with stiffness, if stiffness changes, A need to be recomputed
    std::vector<Triplet<double>> tripletlist;
    for (const auto& e : E) {
        for (int j = 0; j < 3; j++)
        {
            if (true)//!dirichlet_bc_mask[e.first]) 
            {
                tripletlist.push_back(
                    Triplet<double>(3 * e.first + j, 3 * e.first + j, h * h * stiffness));
                tripletlist.push_back(
                Triplet<double>(3 * e.first + j, 3 * e.second + j, -h * h * stiffness));
            }    
            if (true)//!dirichlet_bc_mask[e.second]) 
            {
                tripletlist.push_back(
                    Triplet<double>(3 * e.second + j, 3 * e.first + j, -h * h * stiffness));
                tripletlist.push_back(
                    Triplet<double>(3 * e.second + j, 3 * e.second + j, h * h * stiffness));
            }
        }
    }
    double mass_per_vertex = mass / n_vertices;
    for (int i = 0; i < n_vertices ; i++)
    {
        if (false)//dirichlet_bc_mask[i])
        {
            for (int j = 0; j < 3; j++)
            {
                tripletlist.push_back(
                    Triplet<double>(3 * i + j, 3 * i + j, 1.0));
            }
        }
        else
        {
            for (int j = 0; j < 3; j++) 
            {
                tripletlist.push_back(Triplet<double>(3 * i + j, 3 * i + j, mass_per_vertex));
            }
        }
    }
    A.setFromTriplets(tripletlist.begin(), tripletlist.end());
    A.makeCompressed();
    A_LU.compute(A);
    std::cout << A << std::endl;
}

void FastMassSpring::step()
{
    // (HW Optional) Necessary preparation
    // ...
    Eigen::Vector3d acceleration_ext = gravity + wind_ext_acc;
    unsigned n_vertices = X.rows();
    double mass_per_vertex = mass / n_vertices;

    // Found J   
    J.resize(3 * n_vertices, 3 * E_rest_length.size());
    std::vector<Triplet<double>> tripletlistJ;
    unsigned ik = 0;
    for (const auto& e : E) {
        for (int j = 0; j < 3; j++) {
            tripletlistJ.push_back(Triplet<double>(3 * e.first + j, 3 * ik + j, stiffness));
            tripletlistJ.push_back(Triplet<double>(3 * e.second + j, 3 * ik + j, -stiffness));
        }
        ik++;
    }
    J.setFromTriplets(tripletlistJ.begin(), tripletlistJ.end());

    Eigen::MatrixXd x = MatrixXd::Zero(3 * n_vertices, 1);
    Eigen::MatrixXd Y = MatrixXd::Zero(3 * n_vertices, 1);
    for (int i = 0; i < n_vertices; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            x(3 * i + j, 0) = X(i, j);
            Y(3 * i + j, 0) = X(i, j) + h * vel(i, j) + h * h * acceleration_ext[j]/mass_per_vertex;
        }
    }
    //x = Y;
    Eigen::MatrixXd d = MatrixXd::Zero(3 * E_rest_length.size(), 1);
    for (unsigned iter = 0; iter < max_iter; iter++) {
        // (HW Optional)
        // local_step and global_step alternating solving

        //Get d
        unsigned i = 0;
        for (const auto& e : E)
        {
            Vector3d x1 = { x(3 * e.first, 0),
                            x(3 * e.first+1, 0),
                            x(3 * e.first+2, 0)                              
            };
            Vector3d x2 = { x(3 * e.second, 0),
                            x(3 * e.second + 1, 0),
                            x(3 * e.second + 2, 0)
            };
            Vector3d xi = x1 - x2;
            double norm_xi = xi.norm();
            double Li = E_rest_length[i];
            d(3 * i, 0) = Li * xi[0] / norm_xi;
            d(3 * i + 1, 0) = Li * xi[1] / norm_xi;
            d(3 * i + 2, 0) = Li * xi[2] / norm_xi;            
            i++;
        }
        MatrixXd b = h * h * J * d + mass_per_vertex * Y;
        //fix point
        /* for (int i = 0; i < n_vertices; i++)
        {
            if (dirichlet_bc_mask[i]) {
                for (int j = 0; j < 3; j++) {
                    b(3 * i + j, 0) = X(i,j);
                }
            }
        }*/
        //std::cout << b << std::endl;
        x = A_LU.solve(b);
    }
    //std::cout << X << std::endl;

    //std::cout << x <<std::endl;
    //set X
    /*for (int i = 0; i < n_vertices; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            X(i, j) = x(3 * i + j);
        }
    }*/

    //set vel and X
    for (int i = 0; i < n_vertices; i++) {
        if (!dirichlet_bc_mask[i]) {
            for (int j = 0; j < 3; j++) {
               vel(i, j) = (x(3 * i + j, 0) - X(i, j)) / h;
                X(i, j) = x(3 * i + j, 0);
            }
        }
    }    
    //std::cout << X << std::endl;

}

}  // namespace USTC_CG::node_mass_spring