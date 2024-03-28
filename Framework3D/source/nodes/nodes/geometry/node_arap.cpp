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
    auto vertices_num = halfedge_mesh->n_vertices();
    Eigen::SparseMatrix<double, Eigen::RowMajor> GlobalCoeff(vertices_num-2,vertices_num-2);  
    std::vector<Eigen::Triplet<double>> entries;
    for (const auto& vertex_handle : halfedge_mesh->vertices())
    {
        if (vertex_handle.idx() != vertices_num - 1 && vertex_handle.idx() != vertices_num - 2) {
            auto start = vertex_handle.halfedge();
            auto he = start;
            do {
                auto vertex_j = he.to();
                auto cot_ij_ji = initPara.get_cot(vertex_handle.idx(), vertex_j.idx()) +
                                 initPara.get_cot(vertex_j.idx(), vertex_handle.idx());
                entries.push_back(
                    Eigen::Triplet<double>(vertex_handle.idx(), vertex_handle.idx(),cot_ij_ji ));
                if (vertex_j.idx() != vertices_num - 1 && vertex_j.idx() != vertices_num - 2)
                {
                    entries.push_back(Eigen::Triplet<double>(vertex_handle.idx(), vertex_j.idx(), -cot_ij_ji));
                }
                he = he.prev().opp();
            } while (he != start);
        }
    }
    GlobalCoeff.setFromTriplets(entries.begin(), entries.end());
    Eigen::SparseLU<Eigen::SparseMatrix<double, Eigen::RowMajor>> Global_solver;
    Global_solver.compute(GlobalCoeff);

    //Local Phase
    std::vector<Eigen::Matrix2d> Lt(original_mesh->n_faces(), Eigen::Matrix2d::Zero());
    for (const auto& face_handle : original_mesh->faces())
    {
        auto t_face = face_handle.idx();
        auto e1 = face_handle.halfedge();
        auto e2 = e1.next();
        auto e3 = e2.next();
        std::vector<OpenMesh::SmartVertexHandle> vertex(4);
        vertex[0] = e1.from();
        vertex[1] = e2.from();
        vertex[2] = e3.from();
        vertex[3] = vertex[0];
        //get the related vertex and edges
        std::vector<Eigen::Vector2d> x(4, Eigen::Vector2d::Zero());
        x[0](0) = 0.0;
        x[0](1) = 0.0;        
        x[1](0) = initPara.get_iso(t_face)[0];
        x[1](1) = 0.0;
        x[2](0) = initPara.get_iso(t_face)[1];
        x[2](1) = initPara.get_iso(t_face)[2];
        x[3] = x[0];
        Eigen::Matrix2d S_u = Eigen::Matrix2d::Zero();
        for (int i = 0; i < 3; i++)
        {
            Eigen::MatrixXd u(2, 1);
            Eigen::MatrixXd v(1, 2);
            u(0, 0) = halfedge_mesh->point(vertex[i])[0] - halfedge_mesh->point(vertex[i + 1])[0];
            u(1, 0) = halfedge_mesh->point(vertex[i])[1] - halfedge_mesh->point(vertex[i + 1])[1];
            auto i_idx = vertex[i].idx();
            auto i_1_idx = vertex[i + 1].idx();
            auto cot = initPara.get_cot(vertex[i].idx(), vertex[i + 1].idx());
            v(0, 0) =  cot*(x[i](0) - x[i + 1](0));
            v(0, 1) =  cot*(x[i](1) - x[i + 1](1));
            auto S_tmp = u * v;
            S_u = S_u + S_tmp;
        }
        Eigen::JacobiSVD<Eigen::MatrixXd> svd(S_u, Eigen::ComputeFullV | Eigen::ComputeFullU);
        auto U = svd.matrixU();
        auto V = svd.matrixV();       
        Lt[t_face] = U * V.transpose();
    }
    //Global Phase
    Eigen::SparseMatrix<double, Eigen::RowMajor> b(vertices_num - 2, 2);
    std::vector<Eigen::Triplet<double>> b_entries;
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
        for (int i = 0; i < 3; i++)
        {
            auto cot=initPara.get_cot(vertex[i].idx(), vertex[i + 1].idx());
            Eigen::MatrixXd x_matrix(2, 1);
            x_matrix(0, 0) = cot * (x[i](0) - x[i + 1](0));
            x_matrix(1, 0) = cot * (x[i](1) - x[i + 1](1));
            auto L_x = Lt[t_face] * x_matrix;
            if (vertex[i].idx() != vertices_num - 1 && vertex[i].idx() != vertices_num - 2) {
                b_entries.push_back(Eigen::Triplet<double>(vertex[i].idx(), 0, L_x(0, 0)));
                b_entries.push_back(Eigen::Triplet<double>(vertex[i].idx(), 1, L_x(1, 0)));
            }
            else if (
                vertex[i + 1].idx() == vertices_num - 1 &&
                vertex[i + 1].idx() == vertices_num - 2)
            {
                Eigen::MatrixXd u(2, 1);               
                u(0, 0) =halfedge_mesh->point(vertex[i])[0];
                u(1, 0) =halfedge_mesh->point(vertex[i])[1];
                auto cot_ij_ji = initPara.get_cot(vertex[i].idx(), vertex[i+1].idx()) +
                                 initPara.get_cot(vertex[i + 1].idx(), vertex[i].idx());
                b_entries.push_back(Eigen::Triplet<double>(vertex[i+1].idx(), 0, cot_ij_ji*u(0,0)));
                b_entries.push_back(Eigen::Triplet<double>(vertex[i+1].idx(), 1, cot_ij_ji*u(1,0)));

            }
            if (vertex[i+1].idx() != vertices_num - 1 && vertex[i+1].idx() != vertices_num - 2) {
                b_entries.push_back(Eigen::Triplet<double>(vertex[i + 1].idx(), 0, -L_x(0, 0)));
                b_entries.push_back(Eigen::Triplet<double>(vertex[i + 1].idx(), 1, -L_x(1, 0)));
            }
            else if (vertex[i].idx() == vertices_num - 1 && vertex[i].idx() == vertices_num - 2)
            {
                Eigen::MatrixXd u(2, 1);
                u(0, 0) = halfedge_mesh->point(vertex[i+1])[0];
                u(1, 0) = halfedge_mesh->point(vertex[i+1])[1];
                auto cot_ij_ji = initPara.get_cot(vertex[i].idx(), vertex[i + 1].idx()) +
                                 initPara.get_cot(vertex[i + 1].idx(), vertex[i].idx());
                b_entries.push_back(
                    Eigen::Triplet<double>(vertex[i].idx(), 0, cot_ij_ji * u(0, 0)));
                b_entries.push_back(
                    Eigen::Triplet<double>(vertex[i].idx(), 1, cot_ij_ji * u(1, 0)));
            }
            
        }
    }
    b.setFromTriplets(b_entries.begin(), b_entries.end());
    Eigen::MatrixXd b_dense = b;
    auto solve = Global_solver.solve(b_dense);
    for (const auto& vertex_handle : halfedge_mesh->vertices()) {
        if (vertex_handle.idx() != vertices_num - 1 && vertex_handle.idx() != vertices_num - 2) {
            auto index = vertex_handle.idx();
            halfedge_mesh->point(vertex_handle)[0] = solve(vertex_handle.idx(), 0);
            halfedge_mesh->point(vertex_handle)[1] = solve(vertex_handle.idx(), 1);
            halfedge_mesh->point(vertex_handle)[2] = 0.0;
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
