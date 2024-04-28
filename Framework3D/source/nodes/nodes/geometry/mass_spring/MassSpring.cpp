#include "MassSpring.h"
#include <iostream>

namespace USTC_CG::node_mass_spring {
MassSpring::MassSpring(const Eigen::MatrixXd& X, const EdgeSet& E)
{
    this->X = this->init_X = X;
    this->vel = Eigen::MatrixXd::Zero(X.rows(), X.cols());
    this->E = E;

    std::cout << "number of edges: " << E.size() << std::endl;
    std::cout << "init mass spring" << std::endl;

    // Compute the rest pose edge length
    for (const auto& e : E) {
        Eigen::Vector3d x0 = X.row(e.first);
        Eigen::Vector3d x1 = X.row(e.second);
        this->E_rest_length.push_back((x0 - x1).norm());
    }

    // Initialize the mask for Dirichlet boundary condition
    dirichlet_bc_mask.resize(X.rows(), false);

    // (HW_TODO) Fix two vertices, feel free to modify this 
    unsigned n_fix = sqrt(X.rows());  // Here we assume the cloth is square
    dirichlet_bc_mask[0] = true;
    dirichlet_bc_mask[n_fix - 1] = true;
}

void MassSpring::step()
{
    Eigen::Vector3d acceleration_ext = gravity + wind_ext_acc;

    unsigned n_vertices = X.rows();

    // The reason to not use 1.0 as mass per vertex: the cloth gets heavier as we increase the resolution
    double mass_per_vertex =
        mass / n_vertices; 

    //----------------------------------------------------
    // (HW Optional) Bonus part: Sphere collision
    Eigen::MatrixXd acceleration_collision =
        getSphereCollisionForce(sphere_center.cast<double>(), sphere_radius);
    //----------------------------------------------------

    if (time_integrator == IMPLICIT_EULER) {
        // Implicit Euler
        TIC(step)

        // (HW TODO)
        
        //compute Y
        MatrixXd Y = MatrixXd::Zero(n_vertices * 3, 1);
        MatrixXd currentX = MatrixXd::Zero(n_vertices * 3, 1);
        
        for (size_t i = 0; i < n_vertices; i++) {
            for (int j = 0; j < 3; j++) {    
                currentX(3 * i + j, 0) = X(i, j);
                if (dirichlet_bc_mask[i]) {
                    Y(3 * i + j, 0) = X(i,j);
                }
                else {
                    Y(3 * i + j, 0) = X(i, j) + h * vel(i,j) + h * h *acceleration_ext[j];
                    if (enable_sphere_collision)
                    {
                        Y(3 * i + j, 0) += h * h / mass_per_vertex * acceleration_collision(i, j);
                    }
                }
            }
        }
        MatrixXd x = Y;
        //init f int
        MatrixXd f_int = MatrixXd::Zero(n_vertices, 3);
        Eigen::SparseMatrix<double> g_m;
        Eigen::MatrixXd g;
        g_m.resize(3 * n_vertices, 3 * n_vertices);
        g.resize(3 * n_vertices, 1);
        double maxerror = 0;
        int itr = 0;
        do {
            g_m.setZero();
            //Get f_int and Found Coff Matrix
            std::vector<Triplet<double>> tripletlist;                    
            unsigned i = 0;           
            const auto I = Eigen::MatrixXd::Identity(3, 3);           
            f_int = -computeGrad(stiffness);
            g_m = computeHessianSparse(stiffness);
            Eigen::SparseLU<Eigen::SparseMatrix<double>> g_LU;
            g_LU.compute(g_m);
            //Get grad g(x)            
            for (size_t i = 0; i < n_vertices; i++) {
                for (int j = 0; j < 3; j++) {
                    if (dirichlet_bc_mask[i]) {
                        g(3 * i + j, 0) = 0;
                    }
                    else {
                        g(3 * i + j, 0) =
                            mass_per_vertex * (x(3 * i + j, 0) - Y(3 * i + j, 0)) - h * h * f_int(i,j);
                    }
                }
            }
            Eigen::MatrixXd U = g_LU.solve(g);
            x -= U;
            maxerror = 0;
            for (size_t i = 0; i < 3 * n_vertices; i++) {
                if (maxerror < fabs(U(i, 0))) {
                    maxerror = fabs(U(i, 0));
                }
            }                               
            itr++;
        } while (maxerror > 10e-6 && itr <= 3);

        //set vel
        for (int i = 0; i < n_vertices; i++) {
            if (!dirichlet_bc_mask[i]) {
                for (int j = 0; j < 3; j++) {
                    vel(i, j) = damping*(x(3 * i + j, 0) - X(i,j)) / h;
                    X(i, j) = x(3 * i + j, 0);
                }
            }
        }    
        // Solve Newton's search direction with linear solver 
        
        // update X and vel 

        TOC(step)
    }
    else if (time_integrator == SEMI_IMPLICIT_EULER) {
        // Semi-implicit Euler
        // -----------------------------------------------


        // -----------------------------------------------

        Eigen::MatrixXd acceleration = -computeGrad(stiffness);   
        // (HW Optional)
        if (enable_sphere_collision) {
            acceleration += acceleration_collision;
        }
        acceleration = acceleration / mass_per_vertex;
        acceleration.rowwise() += acceleration_ext.transpose();      
        Vector3d Fix(0.0, 0.0, 0.0);
        for (int i = 0; i < X.rows(); i++)
        {
            if (dirichlet_bc_mask[i])
            {
                acceleration.row(i) = Fix;
            }
        }
        vel += h  * acceleration;
        vel *= damping;
        X += h * vel;
        



        // (HW TODO): Implement semi-implicit Euler time integration

        // Update X and vel 
        
    }
    else {
        std::cerr << "Unknown time integrator!" << std::endl;
        return;
    }
}

// There are different types of mass spring energy:
// For this homework we will adopt Prof. Huamin Wang's energy definition introduced in GAMES103
// course Lecture 2 E = 0.5 * stiffness * sum_{i=1}^{n} (||x_i - x_j|| - l)^2 There exist other
// types of energy definition, e.g., Prof. Minchen Li's energy definition
// https://www.cs.cmu.edu/~15769-f23/lec/3_Mass_Spring_Systems.pdf
double MassSpring::computeEnergy(double stiffness)
{
    double sum = 0.;
    unsigned i = 0;
    for (const auto& e : E) {
        auto diff = X.row(e.first) - X.row(e.second);
        auto l = E_rest_length[i];
        sum += 0.5 * stiffness * std::pow((diff.norm() - l), 2);
        i++;
    }
    return sum;
}

Eigen::MatrixXd MassSpring::computeGrad(double stiffness)
{
    Eigen::MatrixXd g = Eigen::MatrixXd::Zero(X.rows(), X.cols());
    unsigned i = 0;
    

    for (const auto& e : E) {
        // --------------------------------------------------
        // (HW TODO): Implement the gradient computation
        auto xi = X.row(e.first) - X.row(e.second);
        /*auto x1 = xi[0];
        auto x2 = xi[1];
        auto x3 = xi[2];*/
        double l = E_rest_length[i];
        double currentLen = xi.norm();
        auto GradE = stiffness * (currentLen - l) * xi / currentLen;
        g.row(e.first) += GradE;
        g.row(e.second) += -GradE;
        // --------------------------------------------------
        i++;
    }
    return g;
}

Eigen::SparseMatrix<double> MassSpring::computeHessianSparse(double stiffness)
{
    unsigned n_vertices = X.rows();
    Eigen::SparseMatrix<double> H(n_vertices * 3, n_vertices * 3);
    unsigned i = 0;
    auto k = stiffness;
    const auto I = Eigen::MatrixXd::Identity(3, 3);
    std::vector<Triplet<double>> tripletlist;
    Eigen::MatrixXd f_int = Eigen::MatrixXd::Zero(n_vertices, 3);
    for (const auto& e : E) {
        // --------------------------------------------------
        // (HW TODO): Implement the sparse version Hessian computation
        // Remember to consider fixed points 
        // You can also consider positive definiteness here
        MatrixXd f1(3, 1);
        for (int k = 0; k < 3; k++) {
            f1(k, 0) = (X(e.second, k) - X(e.first, k));
        }
        double norm12 = f1.norm();
        double k1 = stiffness * (1 - E_rest_length[i] / norm12);
        for (int k = 0; k < 3; k++) {
            f_int(e.first, k) += k1 * f1(k, 0);
            f_int(e.second, k) -= k1 * f1(k, 0);
        }
        auto f1t = f1.transpose();
        MatrixXd f1_d_x1 = -stiffness * E_rest_length[i] / (norm12 * norm12 * norm12) * f1 * f1t;
        f1_d_x1(0, 0) -= k1;
        f1_d_x1(1, 1) -= k1;
        f1_d_x1(2, 2) -= k1;

        // Save d_f_int
        if (!dirichlet_bc_mask[e.first]) {
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    tripletlist.push_back(
                        Triplet<double>(3 * e.first + i, 3 * e.second + j, h * h * f1_d_x1(i, j)));
                    tripletlist.push_back(
                        Triplet<double>(3 * e.first + i, 3 * e.first + j, -h * h * f1_d_x1(i, j)));
                }
            }
        }
        if (!dirichlet_bc_mask[e.second]) {
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    tripletlist.push_back(
                        Triplet<double>(3 * e.second + i, 3 * e.first + j, h * h * f1_d_x1(i, j)));
                    tripletlist.push_back(Triplet<double>(
                        3 * e.second + i, 3 * e.second + j, -h * h * f1_d_x1(i, j)));
                }
            }
        }
        // --------------------------------------------------

        i++;
    }
    auto mass_per_vertex = mass / n_vertices;
    for (size_t k = 0; k < 3 * n_vertices; k++) {
        tripletlist.push_back(Triplet<double>(k, k, mass_per_vertex));
    }

    H.setFromTriplets(tripletlist.begin(), tripletlist.end());
    H.makeCompressed();
    return H;
}


bool MassSpring::checkSPD(const Eigen::SparseMatrix<double>& A)
{
    // Eigen::SimplicialLDLT<SparseMatrix_d> ldlt(A);
    // return ldlt.info() == Eigen::Success;
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(A);
    auto eigen_values = es.eigenvalues();
    return eigen_values.minCoeff() >= 1e-10;
}

void MassSpring::reset()
{
    std::cout << "reset" << std::endl;
    this->X = this->init_X;
    this->vel.setZero();
}

// ----------------------------------------------------------------------------------
// (HW Optional) Bonus part
Eigen::MatrixXd MassSpring::getSphereCollisionForce(Eigen::Vector3d center, double radius)
{
    Eigen::MatrixXd force = Eigen::MatrixXd::Zero(X.rows(), X.cols());
    for (int i = 0; i < X.rows(); i++) {
        Vector3d x_c = { X(i, 0) - center[0], X(i, 1) - center[1], X(i, 2) - center[2] };
        //auto x_c = X.row(i) - center;
         double norm = x_c.norm();
         auto dis = collision_scale_factor * radius - norm;
        if (dis > 0)
        {
            Vector3d f_sphere = collision_penalty_k * dis * x_c / norm;
            force.row(i) = f_sphere;
        }

       // (HW Optional) Implement penalty-based force here 
    }
    return force;
}
// ----------------------------------------------------------------------------------


}  // namespace USTC_CG::node_mass_spring

