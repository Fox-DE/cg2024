#include "algorithm.h"
#include <omp.h>
using namespace Eigen;

namespace USTC_CG {

void ARAP::get_data(
        std::shared_ptr<PolyMesh> original_mesh_in,
        std::shared_ptr<PolyMesh> halfedge_mesh_in)
{
    original_mesh = original_mesh_in;
    halfedge_mesh = halfedge_mesh_in;
}

void ARAP::compute()
{
    init();
    import_u_data();
    build_A();
    for (int itr_tmp = 0; itr_tmp < itr; itr_tmp++)
    {
        MatrixXd L(2 * nF + 1, 2);
        MatrixXd edge;
        MatrixXd L_t;
        #pragma omp parallel for
        for (int i = 0; i < nF; i++)
        {
            edge = edge_s.block(2 * i, 0, 2, 2);
            Matrix2d temp;
            temp << u(mat_index(i, 1), 0) - u(mat_index(i, 0), 0), u(mat_index(i, 1), 1) - u(mat_index(i, 0), 1),
                u(mat_index(i, 2), 0) - u(mat_index(i, 0), 0), u(mat_index(i, 2), 1) - u(mat_index(i, 0), 1);
            temp = edge * temp;
            JacobiSVD<Matrix2d> svd(temp, ComputeFullU | ComputeFullV);
            auto U = svd.matrixU();
            auto V = svd.matrixV();
            if (temp.determinant() < 0)
            {
                U(0, 1) = -U(0, 1);
                U(1, 1) = -U(1, 1);
            }
            L_t = sqrt(area(i, 0)) * U * V.transpose();
            L(2 * i, 0) = L_t(0, 0);
            L(2 * i, 1) = L_t(0, 1);
            L(2 * i + 1, 0) = L_t(1, 0);
            L(2 * i + 1, 1) = L_t(1, 1);
        }
        L(2 * nF, 0) = area.sum();
        L(2 * nF, 1) = 0;
        
        u = solver.solve(mat_A.transpose() * L);
        
    }
    update_u_data();
}
void ARAP::set_itr(int i)
{
    itr = i>1?i:1;
}

void ARAP::init()
{
    nV = original_mesh->n_vertices();
    nF = original_mesh->n_faces();
    area = MatrixXd::Zero(nF, 1);
    angle = MatrixXd::Zero(nF, 1);
    edge_s = MatrixXd::Zero(2 * nF, 2);
    mat_index = MatrixXi::Zero(nF, 3);
    mat_A .resize(2 * nF + 1, nV);

}

void ARAP::build_A()
{
    std::vector<Eigen::Triplet<double>> triplet_list;
    OpenMesh::DefaultTraits::Point p1, p2, p3;
    int i = 0, idx1, idx2, idx3;
    for (const auto& face_handle : original_mesh->faces())
    {
        auto v_itr = face_handle.vertices().begin();
        p1 = original_mesh->point(v_itr);
        idx1 = v_itr->idx();
        v_itr++;
        p2 = original_mesh->point(v_itr);
        idx2 = v_itr->idx();
        v_itr++;
        p3 = original_mesh->point(v_itr);
        idx3 = v_itr->idx();
        //obtain the point
        mat_index(i, 0) = idx1;
        mat_index(i, 1) = idx2;
        mat_index(i, 2) = idx3;
        //set the idx
        auto e21 = p2 - p1;
        auto e31 = p3 - p1;
        double e21_norm = e21.norm();
        double e31_norm = e31.norm();
        area(i, 0) = e21.cross(e31).norm() / 2;
        double sqrt_area = sqrt(area(i, 0));
        auto tmp = e21.dot(e31);        
        angle(i, 0) = acos(e21.dot(e31) / (e21_norm*e31_norm));
        edge_s(2 * i, 0) = 1 / e21_norm;
        edge_s(2 * i, 1) = 0;
        edge_s(2 * i + 1, 0) = -cos(angle(i, 0)) / (sin(angle(i, 0)) * e21_norm);
        edge_s(2 * i + 1, 1) = 1 / (sin(angle(i, 0)) * e31_norm);
        triplet_list.push_back(Eigen::Triplet<double>(2 * i, idx1, -edge_s(2 * i, 0) * sqrt_area));
        triplet_list.push_back(Eigen::Triplet<double>(2 * i, idx2, edge_s(2 * i, 0) * sqrt_area));
        triplet_list.push_back(Eigen::Triplet<double>(
            2 * i + 1, idx1, (-edge_s(2 * i + 1, 0) - edge_s(2 * i + 1, 1)) * sqrt_area));
        triplet_list.push_back(
            Eigen::Triplet<double>(2 * i + 1, idx2, edge_s(2 * i + 1, 0) * sqrt_area));
        triplet_list.push_back(
            Eigen::Triplet<double>(2 * i + 1, idx3, edge_s(2 * i + 1, 1) * sqrt_area));
        i++;
    }
    triplet_list.push_back(Eigen::Triplet<double>(2 * nF, 0, area.sum()));
    mat_A.setFromTriplets(triplet_list.begin(), triplet_list.end());
    mat_A.makeCompressed();
    

    solver.compute(mat_A.transpose() * mat_A);
}

void ARAP::import_u_data()
{
    u = MatrixXd::Zero(nV, 2);
#pragma omp parallel for
    for (const auto& vertex_handle:halfedge_mesh->vertices())
    {
        auto p = halfedge_mesh->point(vertex_handle);
        u(vertex_handle.idx(), 0) = p[0];
        u(vertex_handle.idx(), 1) = p[1];
    }
    
}
void ARAP::update_u_data()
{
#pragma omp parallel for
    for (const auto& vertex_handle : halfedge_mesh->vertices())
    {
        auto& po = halfedge_mesh->point(vertex_handle);
        auto idx = vertex_handle.idx();
        po[0] = u(idx, 0);
        po[1] = u(idx, 1);
        po[2] = 0.0;

    }
}
std::shared_ptr<PolyMesh> ARAP::output_mesh()
{
    return halfedge_mesh;
}

}
