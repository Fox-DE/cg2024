#pragma once
#include "GCore/Components/MeshOperand.h"
#include "Nodes/node.hpp"
#include "Nodes/node_declare.hpp"
#include "Nodes/node_register.h"
#include "geom_node_base.h"
#include "utils/util_openmesh_bind.h"
#include "utils/algorithm.h"
#include "Eigen/Sparse"
#include "Eigen/Dense"

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

namespace USTC_CG::node_arap {
static void node_arap_declare(NodeDeclarationBuilder& b)
{
    // Input-1: Original 3D mesh with boundary
    // Maybe you need to add another input for initialization?
    b.add_input<decl::Geometry>("Input");
    b.add_input<decl::Geometry>("OriginalMesh");
    /*
    ** NOTE: You can add more inputs or outputs if necessary. For example, in
    ** some cases, additional information (e.g. other mesh geometry, other 
    ** parameters) is required to perform the computation.
    **
    ** Be sure that the input/outputs do not share the same name. You can add
    ** one geometry as
    **
    **                b.add_input<decl::Geometry>("Input");
    **
    ** Or maybe you need a value buffer like:
    **
    **                b.add_input<decl::Float1Buffer>("Weights");
    */

    // Output-1: The UV coordinate of the mesh, provided by ARAP algorithm
    b.add_output<decl::Geometry>("OutputMesh");
    //b.add_output<decl::Float2Buffer>("OutputUV");
}

static void node_arap_exec(ExeParams params)
{
    // Get the input from params
    auto input = params.get_input<GOperandBase>("Input");
    auto input2 = params.get_input<GOperandBase>("OriginalMesh");
    // Avoid processing the node when there is no input
    if (!input.get_component<MeshComponent>()) {
        throw std::runtime_error("Need Geometry Input.");
    }
    if (!input2.get_component<MeshComponent>()) {
        throw std::runtime_error("Minimal Surface: Need Geometry Input.");
    }
    //throw std::runtime_error("Not implemented");

    /* ----------------------------- Preprocess -------------------------------
    ** Create a halfedge structure (using OpenMesh) for the input mesh. The
    ** half-edge data structure is a widely used data structure in geometric
    ** processing, offering convenient operations for traversing and modifying
    ** mesh elements.
    */
    auto halfedge_mesh = operand_to_openmesh(&input);
    auto original_mesh = operand_to_openmesh(&input2);

   /* ------------- [HW5_TODO] ARAP Parameterization Implementation -----------
   ** Implement ARAP mesh parameterization to minimize local distortion.
   **
   ** Steps:
   ** 1. Initial Setup: Use a HW4 parameterization result as initial setup.
   **
   ** 2. Local Phase: For each triangle, compute local orthogonal approximation
   **    (Lt) by computing SVD of Jacobian(Jt) with fixed u.
   **
   ** 3. Global Phase: With Lt fixed, update parameter coordinates(u) by solving
   **    a pre-factored global sparse linear system.
   **
   ** 4. Iteration: Repeat Steps 2 and 3 to refine parameterization.
   **
   ** Note:
   **  - Fixed points' selection is crucial for ARAP and ASAP.
   **  - Encapsulate algorithms into classes for modularity.
   */
    init_Parameterization initPara;
    initPara.init(original_mesh);

    //build global phase Coefficient Matrix
    auto nV = halfedge_mesh->n_vertices();
    Eigen::SparseMatrix<double, Eigen::RowMajor> mat_A(nV, nV);
    std::vector<Eigen::Triplet<double>> ARAP_coeff;   
    for (const auto& vertex_handle : halfedge_mesh->vertices())
    {
        if (vertex_handle.idx() != nV) {
            auto start = vertex_handle.halfedge();
            auto he = start;
            double cot_sum = 0;
            do {
                auto v_to = he.to();
                auto cot_ij_ji = initPara.get_cot(vertex_handle.idx(), v_to.idx()) +
                                 initPara.get_cot(v_to.idx(), vertex_handle.idx());
                ARAP_coeff.push_back(Eigen::Triplet<double>(vertex_handle.idx(), v_to.idx(), -cot_ij_ji));
                cot_sum += cot_ij_ji;
                he = he.prev().opp();
            } while (he != start);
            ARAP_coeff.push_back(
                Eigen::Triplet<double>(vertex_handle.idx(), vertex_handle.idx(), cot_sum));
        }
    }
    ARAP_coeff.push_back(Eigen::Triplet<double>(nV - 1, nV - 1, 1));

    mat_A.setFromTriplets(ARAP_coeff.begin(), ARAP_coeff.end());
    Eigen::SparseLU<Eigen::SparseMatrix<double, Eigen::RowMajor>> Global_solver;
    mat_A.makeCompressed();
    Global_solver.compute(mat_A.transpose()*mat_A);


    for (int itr = 0; itr < 5; itr++) {
        // Local Phase
        std::vector<Eigen::Matrix2d> Lt(original_mesh->n_faces(), Eigen::Matrix2d::Zero());
        for (const auto& face_handle : original_mesh->faces()) {
            auto t_face = face_handle.idx();
            auto e1 = face_handle.halfedge();
            auto e2 = e1.next();
            auto e3 = e2.next();
            std::vector<OpenMesh::SmartVertexHandle> vertex(4);
            vertex[0] = e1.from();
            vertex[1] = e2.from();
            vertex[2] = e3.from();
            vertex[3] = vertex[0];
            // get the related vertex and edges
            std::vector<Eigen::Vector2d> x(4, Eigen::Vector2d::Zero());
            x[0](0) = 0.0;
            x[0](1) = 0.0;
            x[1](0) = initPara.get_iso(t_face)[0];
            x[1](1) = 0.0;
            x[2](0) = initPara.get_iso(t_face)[1];
            x[2](1) = initPara.get_iso(t_face)[2];
            x[3] = x[0];
            Eigen::MatrixXd U_t(2, 2);
            Eigen::MatrixXd V_t(2, 2);
            auto u2_1 = halfedge_mesh->point(vertex[1]) - halfedge_mesh->point(vertex[0]);
            auto u3_1 = halfedge_mesh->point(vertex[2]) - halfedge_mesh->point(vertex[0]);
            U_t(0, 0) = u2_1[0];
            U_t(1, 0) = u2_1[1];
            U_t(0, 1) = u3_1[0];
            U_t(1, 1) = u3_1[1];
            V_t(0, 0) = x[1](0) - x[0](0);
            V_t(1, 0) = x[1](1) - x[0](1);
            V_t(0, 1) = x[2](0) - x[0](0);
            V_t(1, 1) = x[2](1) - x[0](1);
            auto V_inv = V_t.inverse();
            auto Lt_tmp = U_t * V_inv;
            /*if (Lt_tmp.determinant() < 0)
            {
                U_t(0, 1) = -U_t(0, 1);
                U_t(1, 1) = -U_t(1, 1);
                auto Lt_tmp = U_t * V_inv;
            }*/
            Eigen::JacobiSVD<Eigen::MatrixXd> svd(Lt_tmp, Eigen::ComputeFullV | Eigen::ComputeFullU);
            auto U = svd.matrixU();
            auto V = svd.matrixV();
            Lt[t_face] = U * V.transpose();
            //debug use
            if (Lt[t_face].determinant() < 0)
            {
                printf("<0");
            }
        }
        // Global Phase
        Eigen::SparseMatrix<double, Eigen::RowMajor> b(nV, 2);
        std::vector<Eigen::Triplet<double>> b_entries;
        b_entries.push_back(Eigen::Triplet<double>(
            nV - 1, 0, halfedge_mesh->point(halfedge_mesh->vertex_handle(nV - 1))[0]));
        b_entries.push_back(Eigen::Triplet<double>(
            nV - 1, 1, halfedge_mesh->point(halfedge_mesh->vertex_handle(nV - 1))[1]));
        for (const auto& face_handle : original_mesh->faces()) {
            auto t_face = face_handle.idx();
            auto e1 = face_handle.halfedge();
            auto e2 = e1.next();
            auto e3 = e2.next();
            std::vector<OpenMesh::SmartVertexHandle> vertex(4);
            vertex[0] = e1.from();
            vertex[1] = e2.from();
            vertex[2] = e3.from();
            vertex[3] = vertex[0];
            // get the related vertex and edges
            std::vector<Eigen::Vector2d> x(4, Eigen::Vector2d::Zero());
            x[0](0) = 0.0;
            x[0](1) = 0.0;
            x[1](0) = initPara.get_iso(t_face)[0];
            x[1](1) = 0.0;
            x[2](0) = initPara.get_iso(t_face)[1];
            x[2](1) = initPara.get_iso(t_face)[2];
            x[3] = x[0];
            for (int i = 0; i < 3; i++) {
                auto cot = initPara.get_cot(vertex[i].idx(), vertex[i + 1].idx());
                Eigen::MatrixXd x_matrix(2, 1);
                x_matrix(0, 0) = cot * (x[i](0) - x[i + 1](0));
                x_matrix(1, 0) = cot * (x[i](1) - x[i + 1](1));
                auto L_x = Lt[t_face] * x_matrix;
                if (vertex[i].idx() != nV - 1) {
                    b_entries.push_back(Eigen::Triplet<double>(vertex[i].idx(), 0, L_x(0, 0)));
                    b_entries.push_back(Eigen::Triplet<double>(vertex[i].idx(), 1, L_x(1, 0)));
                }
                
                if (vertex[i + 1].idx() != nV - 1) {
                    b_entries.push_back(Eigen::Triplet<double>(vertex[i + 1].idx(), 0, -1.0*L_x(0, 0)));
                    b_entries.push_back(Eigen::Triplet<double>(vertex[i + 1].idx(), 1, -1.0*L_x(1, 0)));
                }
                
            }
        }

        
        b.setFromTriplets(b_entries.begin(), b_entries.end());
        Eigen::MatrixXd b_dense = b;
        auto solve = Global_solver.solve(mat_A.transpose()*b_dense);
        for (const auto& vertex_handle : halfedge_mesh->vertices()) {
            
                auto index = vertex_handle.idx();
                halfedge_mesh->point(vertex_handle)[0] = solve(vertex_handle.idx(), 0);
                halfedge_mesh->point(vertex_handle)[1] = solve(vertex_handle.idx(), 1);
                halfedge_mesh->point(vertex_handle)[2] = 0.0;
                if (vertex_handle.idx() == nV - 1)
                {
                printf("debug use");
                }
        }
    }





    // The result UV coordinates 
    //pxr::VtArray<pxr::GfVec2f> uv_result;
    
    auto operand_base = openmesh_to_operand(halfedge_mesh.get());
    // Set the output of the node
    params.set_output("OutputMesh", std::move(*operand_base));
    //params.set_output("OutputUV", uv_result);
}

static void node_register()
{
    static NodeTypeInfo ntype;

    strcpy(ntype.ui_name, "ARAP Parameterization");
    strcpy_s(ntype.id_name, "geom_arap");

    geo_node_type_base(&ntype);
    ntype.node_execute = node_arap_exec;
    ntype.declare = node_arap_declare;
    nodeRegisterType(&ntype);
}

NOD_REGISTER_NODE(node_register)
}  // namespace USTC_CG::node_arap
