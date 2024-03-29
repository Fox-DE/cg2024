#include <Eigen/sparse>

#include "Eigen/Eigen"
#include "GCore/Components/MeshOperand.h"
#include "Nodes/node.hpp"
#include "Nodes/node_declare.hpp"
#include "Nodes/node_register.h"
#include "geom_node_base.h"
#include "utils/util_openmesh_bind.h"
/*
** @brief HW5_ARAP_Parameterization
**
** This file presents the basic framework of a "node", which processes inputs
** received from the left and outputs specific variables for downstream nodes to
** use.
**
** - In the first function, node_declare, you can set up the node's input and
** output variables.
**
** - The second function, node_exec is the execution part of the node, where we
** need to implement the node's functionality.
**
** - The third function generates the node's registration information, which
** eventually allows placing this node in the GUI interface.
**
** Your task is to fill in the required logic at the specified locations
** within this template, especially in node_exec.
*/

namespace USTC_CG::node_hao_arap {
static void node_arap_declare(NodeDeclarationBuilder& b)
{
    // Input-1: Original 3D mesh with boundary
    // Maybe you need to add another input for initialization?
    b.add_input<decl::Geometry>("Geometry");
    b.add_input<decl::Float2Buffer>("InputUV");

    // Output-1: The UV coordinate of the mesh, provided by ARAP algorithm
    b.add_output<decl::Geometry>("Output");
}

static void node_arap_exec(ExeParams params)
{
    // Get the input from params
    auto Geometry = params.get_input<GOperandBase>("Geometry");
    auto InputUV = params.get_input<pxr::VtArray<pxr::GfVec2f>>("InputUV");
    // Avoid processing the node when there is no input
    if (!Geometry.get_component<MeshComponent>()) {
        throw std::runtime_error("Need Geometry Input.");
    }
    auto halfedge_mesh = operand_to_openmesh(&Geometry);
    // halfedge_mesh means the original version
    int nv = halfedge_mesh->n_vertices(), nf = halfedge_mesh->n_faces();
    // nv is the num of vertex,nf similar
    Eigen::MatrixXd area(nf, 1), angle(nf, 1), edge1(2 * nf, 2);
    Eigen::MatrixXi t(nf, 3);
    std::vector<Eigen::Triplet<double>> triplet_list;
    OpenMesh::DefaultTraits::Point p1, p2, p3;
    int i = 0, idx1, idx2, idx3;
    int n = 200;
    for (const auto& face_handle : halfedge_mesh->faces()) {
        auto temp = face_handle.vertices().begin();
        // temp is the vertex of the face
        p1 = halfedge_mesh->point(temp);
        idx1 = (temp++)->idx();
        // represent temp's value before it moves
        p2 = halfedge_mesh->point(temp);
        idx2 = (temp++)->idx();
        p3 = halfedge_mesh->point((temp));
        idx3 = (temp)->idx();
        // i represent the face
        t(i, 0) = idx1;
        t(i, 1) = idx2;
        t(i, 2) = idx3;
        auto e21 = p2 - p1, e31 = p3 - p1;

        double e21_norm = e21.norm(), e31_norm = e31.norm();

        area(i, 0) = e21.cross(e31).norm() / 2;
        // face i 's area
        double sqrt_area = sqrt(area(i, 0));
        auto tmp = e21.dot(e31);
        angle(i, 0) = acos(e21.dot(e31) / (e21_norm * e31_norm));
        // angle represent the angle of angle 1
        edge1(2 * i, 0) = 1 / e21_norm;
        edge1(2 * i, 1) = 0;
        edge1(2 * i + 1, 0) = -cos(angle(i, 0)) / (sin(angle(i, 0)) * e21_norm);
        edge1(2 * i + 1, 1) = 1 / (sin(angle(i, 0)) * e31_norm);
        triplet_list.push_back(Eigen::Triplet<double>(2 * i, idx1, -edge1(2 * i, 0) * sqrt_area));
        triplet_list.push_back(Eigen::Triplet<double>(2 * i, idx2, edge1(2 * i, 0) * sqrt_area));
        triplet_list.push_back(Eigen::Triplet<double>(
            2 * i + 1, idx1, (-edge1(2 * i + 1, 0) - edge1(2 * i + 1, 1)) * sqrt_area));
        triplet_list.push_back(
            Eigen::Triplet<double>(2 * i + 1, idx2, edge1(2 * i + 1, 0) * sqrt_area));
        triplet_list.push_back(
            Eigen::Triplet<double>(2 * i + 1, idx3, edge1(2 * i + 1, 1) * sqrt_area));
        i++;
    }
    triplet_list.push_back(Eigen::Triplet<double>(2 * nf, 0, area.sum()));
    Eigen::SparseMatrix<double> M(2 * nf + 1, nv), MT;
    M.setFromTriplets(triplet_list.begin(), triplet_list.end());
    MT = M.transpose();
    // auto MM=MT* M;

    // debug use
    // std::cout << M<<std::endl;

    Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> solver;
    solver.compute(MT * M);
    // build Matrix

    Eigen::MatrixXd u(nv, 2);

    Eigen::Matrix2d edge;
    Eigen::Matrix2d temp;
    Eigen::Matrix2d U;
    Eigen::Matrix2d V;
    Eigen::Matrix2d L1;
    Eigen::MatrixXd L(2 * nf + 1, 2);
    L(2 * nf, 0) = area.sum();
    L(2 * nf, 1) = 0;
    // import the data
    for (i = 0; i < nv; i++) {
        u(i, 0) = InputUV[i].data()[0];
        u(i, 1) = InputUV[i].data()[1];
    }
    for (int k = 0; k < n; k++) {
        for (i = 0; i < nf; i++) {
            edge = edge1.block(2 * i, 0, 2, 2);
            temp << u(t(i, 1), 0) - u(t(i, 0), 0), u(t(i, 1), 1) - u(t(i, 0), 1),
                u(t(i, 2), 0) - u(t(i, 0), 0), u(t(i, 2), 1) - u(t(i, 0), 1);
            temp = edge * temp;
            Eigen::JacobiSVD<Eigen::Matrix2d> svd(temp, Eigen::ComputeFullU | Eigen::ComputeFullV);
            U = svd.matrixU();
            V = svd.matrixV();
            if (temp.determinant() < 0) {
                U(0, 1) = -U(0, 1);
                U(1, 1) = -U(1, 1);
            }
            L1 = sqrt(area(i, 0)) * U * V.transpose();
            // std::cout << L1 << std::endl << std::endl;
            L(2 * i, 0) = L1(0, 0);
            L(2 * i, 1) = L1(0, 1);
            L(2 * i + 1, 0) = L1(1, 0);
            L(2 * i + 1, 1) = L1(1, 1);
        }
        std::cout << L << std::endl;
        u = solver.solve(MT * L);
        auto MTL = MT * L;
        std::cout << MTL << std::endl;
        // std::cout << u << std::endl;
    }

    // update vertex
    for (const auto& vertex_handle : halfedge_mesh->vertices()) {
        auto& position = halfedge_mesh->point(vertex_handle);
        int vertex_idx = vertex_handle.idx();
        position[0] = u(vertex_idx, 0);
        position[1] = u(vertex_idx, 1);
        position[2] = 0;
    }
    // The result UV coordinates
    auto operand_base = openmesh_to_operand(halfedge_mesh.get());
    // Set the output of the nodes
    params.set_output("Output", std::move(*operand_base));
}

static void node_register()
{
    static NodeTypeInfo ntype;

    strcpy(ntype.ui_name, "ARAP Parameterization 1");
    strcpy_s(ntype.id_name, "geom_arap_hao");

    geo_node_type_base(&ntype);
    ntype.node_execute = node_arap_exec;
    ntype.declare = node_arap_declare;
    nodeRegisterType(&ntype);
}

NOD_REGISTER_NODE(node_register)
}  // namespace USTC_CG::node_arap
