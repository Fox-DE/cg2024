#pragma once
#include "Eigen/Dense"
#include "Eigen/Sparse"
#include "GCore/Components/MeshOperand.h"
#include "Nodes/node.hpp"
#include "Nodes/node_declare.hpp"
#include "Nodes/node_register.h"
#include "geom_node_base.h"
#include "utils/algorithm.h"
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

namespace USTC_CG::node_arap {
static void node_arap_declare(NodeDeclarationBuilder& b)
{
    // Input-1: Original 3D mesh with boundary
    // Maybe you need to add another input for initialization?
    b.add_input<decl::Geometry>("Input");
    b.add_input<decl::Geometry>("OriginalMesh");
    b.add_input<decl::Int>("Iteration Times").min(1).max(200).default_val(100);
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
    // b.add_output<decl::Float2Buffer>("OutputUV");
}

static void node_arap_exec(ExeParams params)
{
    // Get the input from params
    auto input = params.get_input<GOperandBase>("Input");
    auto input2 = params.get_input<GOperandBase>("OriginalMesh");
    auto itr = params.get_input<int>("Iteration Times");
    // Avoid processing the node when there is no input
    if (!input.get_component<MeshComponent>()) {
        throw std::runtime_error("Need Geometry Input.");
    }
    if (!input2.get_component<MeshComponent>()) {
        throw std::runtime_error("Minimal Surface: Need Geometry Input.");
    }
    // throw std::runtime_error("Not implemented");

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
    ARAP arap;
    arap.get_data(original_mesh, halfedge_mesh);
    arap.set_itr(itr);
    arap.compute();
    halfedge_mesh = arap.output_mesh();
    

    // The result UV coordinates
    // pxr::VtArray<pxr::GfVec2f> uv_result;

    auto operand_base = openmesh_to_operand(halfedge_mesh.get());
    // Set the output of the node
    params.set_output("OutputMesh", std::move(*operand_base));
    // params.set_output("OutputUV", uv_result);
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
