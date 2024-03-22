#include "GCore/Components/MeshOperand.h"
#include "Nodes/node.hpp"
#include "Nodes/node_declare.hpp"
#include "Nodes/node_register.h"
#include "geom_node_base.h"
#include "utils/util_openmesh_bind.h"
#include "Eigen/Sparse"
#include <cmath>


namespace USTC_CG::node_min_surf_cot_weights {
static void node_min_surf_cotangent_weights_declare(NodeDeclarationBuilder& b)
{
    // Input-1: Original 3D mesh with boundary
    b.add_input<decl::Geometry>("Input");
    b.add_input<decl::Geometry>("OriginalMesh");

    // Output-1: Minimal surface with fixed boundary
    b.add_output<decl::Geometry>("Output");
}

static void node_min_surf_cotangent_weights_exec(ExeParams params)
{
    // Get the input from params
    auto input = params.get_input<GOperandBase>("Input"); 
    auto input2 = params.get_input<GOperandBase>("OriginalMesh");
    // (TO BE UPDATED) Avoid processing the node when there is no input 
    if (!input.get_component<MeshComponent>()) {
        throw std::runtime_error("Minimal Surface: Need Geometry Input.");
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
    /* ---------------- [HW4_TODO] TASK 1: Minimal Surface --------------------
    ** In this task, you are required to generate a 'minimal surface' mesh with
    ** the boundary of the input mesh as its boundary.
    **
    ** Specifically, the positions of the boundary vertices of the input mesh
    ** should be fixed. By solving a global Laplace equation on the mesh,
    ** recalculate the coordinates of the vertices inside the mesh to achieve
    ** the minimal surface configuration.
    **
    ** (Recall the Poisson equation with Dirichlet Boundary Condition in HW3)
    */

    /*
    ** Algorithm Pseudocode for Minimal Surface Calculation
    ** ------------------------------------------------------------------------
    ** 1. Initialize mesh with input boundary conditions.
    **    - For each boundary vertex, fix its position.
    **    - For internal vertices, initialize with initial guess if necessary.
    **
    ** 2. Construct Laplacian matrix for the mesh.
    **    - Compute weights for each edge based on the chosen weighting scheme
    **      (e.g., uniform weights for simplicity).
    **    - Assemble the global Laplacian matrix.
    **
    ** 3. Solve the Laplace equation for interior vertices.
    **    - Apply Dirichlet boundary conditions for boundary vertices.
    **    - Solve the linear system (Laplacian * X = 0) to find new positions
    **      for internal vertices.
    **
    ** 4. Update mesh geometry with new vertex positions.
    **    - Ensure the mesh respects the minimal surface configuration.
    **
    ** Note: This pseudocode outlines the general steps for calculating a
    ** minimal surface mesh given fixed boundary conditions using the Laplace
    ** equation. The specific implementation details may vary based on the mesh
    ** representation and numerical methods used.
    */
    Eigen::VectorXi Index_(halfedge_mesh->n_vertices());
    int count_num = 0;
    for (const auto& vertex_handle : halfedge_mesh->vertices()) {
        if (!vertex_handle.is_boundary())
        {
            auto index = vertex_handle.idx();
            Index_(index) = count_num;
            count_num++;
        }
    }
    Eigen::SparseMatrix<double, Eigen::RowMajor> A(count_num,count_num);
    Eigen::MatrixXd b =Eigen::MatrixXd::Zero(count_num,3);
    for (const auto& vertex_handle : original_mesh->vertices())
    {
        if (!vertex_handle.is_boundary()) {
            auto index = vertex_handle.idx();
            double sum_w = 0.0;
            auto start = vertex_handle.halfedge();
            auto he = start;
            const auto start_point = original_mesh->point(vertex_handle);
            do {
                auto Vertex_near = he.to();
                auto last = he.opp().next();
                auto next = he.prev().opp();
                const auto vec_tmp = original_mesh->point(Vertex_near); 
                const auto last_point = original_mesh->point(last.to());
                const auto next_point = original_mesh->point(next.to());
                const auto vec1_last = start_point - last_point;
                const auto vec2_last = vec_tmp - last_point;
                const auto vec1_next = start_point - next_point;
                const auto vec2_next = vec_tmp - next_point;

                double cos_last = vec1_last.dot(vec2_last) / (vec1_last.norm() * vec2_last.norm());
                double cos_next = vec1_next.dot(vec2_next) / (vec1_next.norm() * vec2_next.norm());
                double cot_last = cos_last / sqrt(1.0 - cos_last * cos_last);
                double cot_next = cos_next / sqrt(1.0 - cos_next * cos_next);
                double w = cot_last + cot_next;
                sum_w += w;
                if (Vertex_near.is_boundary())
                {
                    
                    const auto vec_tmp_input1 = halfedge_mesh->point(Vertex_near);
                    
                    b(Index_(index), 0) += w * vec_tmp_input1[0];
                    b(Index_(index), 1) += w * vec_tmp_input1[1];
                    b(Index_(index), 2) += w * vec_tmp_input1[2];
                }
                else
                {
                    A.insert(Index_(index), Index_(Vertex_near.idx())) = -w;
                }

                he = he.prev().opp();
            } while (he != start);
            A.insert(Index_(index), Index_(index)) = sum_w;
        }
    }
    Eigen::SparseLU<Eigen::SparseMatrix<double, Eigen::RowMajor>> solver;
    solver.compute(A);
    Eigen::MatrixXd f=solver.solve(b);
    for (const auto& vertex_handle : halfedge_mesh->vertices()) {
        if (!vertex_handle.is_boundary()) {
            auto index = vertex_handle.idx();
            halfedge_mesh->point(vertex_handle)[0] = f(Index_(vertex_handle.idx()), 0);
            halfedge_mesh->point(vertex_handle)[1] = f(Index_(vertex_handle.idx()), 1);
            halfedge_mesh->point(vertex_handle)[2] = f(Index_(vertex_handle.idx()), 2);
        }
    }
    /* ----------------------------- Postprocess ------------------------------
    ** Convert the minimal surface mesh from the halfedge structure back to
    ** GOperandBase format as the node's output.
    */
    auto operand_base = openmesh_to_operand(halfedge_mesh.get());

    // Set the output of the nodes
    params.set_output("Output", std::move(*operand_base));
}

static void node_register()
{
    static NodeTypeInfo ntype;

    strcpy(ntype.ui_name, "Minimal Surface Cotangent Weights");
    strcpy_s(ntype.id_name, "geom_min_surf_cot");

    geo_node_type_base(&ntype);
    ntype.node_execute = node_min_surf_cotangent_weights_exec;
    ntype.declare = node_min_surf_cotangent_weights_declare;
    nodeRegisterType(&ntype);
}

NOD_REGISTER_NODE(node_register)
}  // namespace USTC_CG::node_min_surf_cot_weights
