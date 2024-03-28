#include "algorithm.h"


double cos2cot(double input)
{
    return input / sqrt(1.0 - input * input);
}

namespace USTC_CG {
void init_Parameterization::init(std::shared_ptr<PolyMesh> input)
{
    original_mesh = input;
    Isometric_parameterization();
    

}

void init_Parameterization::Isometric_parameterization()
{
    auto face_num = original_mesh->n_faces();
    auto vertex_num = original_mesh->n_vertices();
    Iso_index = Eigen::MatrixXd::Zero(face_num, 3);
    cot_index = Eigen::MatrixXd::Zero(vertex_num,vertex_num);
    for (const auto& face_handle : original_mesh->faces())
    {
        auto e1 = face_handle.halfedge();
        auto e2 = e1.next();
        auto e3 = e2.next();
        auto vertex1 = e1.from();
        auto vertex2 = e2.from();
        auto vertex3 = e3.from();
        auto vec1 = original_mesh->point(vertex2) - original_mesh->point(vertex1);
        auto vec2 = original_mesh->point(vertex3) - original_mesh->point(vertex2);
        auto vec3 = original_mesh->point(vertex1) - original_mesh->point(vertex3);
        auto face_idx = face_handle.idx();
        Iso_index(face_idx, 0) = vec1.norm();
        double cos1 = -vec1.dot(vec3) / (vec1.norm() * vec3.norm());
        double cos2 = -vec1.dot(vec2) / (vec1.norm() * vec2.norm());
        double cos3 = -vec3.dot(vec2) / (vec3.norm() * vec2.norm());
        Iso_index(face_idx, 1) = vec3.norm() * cos1;
        Iso_index(face_idx, 2) = vec3.norm() * sqrt(1.0 - cos1 * cos1);
        cot_index(vertex1.idx(), vertex2.idx()) = cos2cot(cos3);
        cot_index(vertex2.idx(), vertex3.idx()) = cos2cot(cos1);
        cot_index(vertex3.idx(), vertex1.idx()) = cos2cot(cos2);
        
    }
}

double init_Parameterization::get_cot(int i, int j)
{
    return cot_index(i, j);
}

std::vector<double> init_Parameterization::get_iso(int i)
{
    std::vector<double> output(3, 0.0);
    output[0] = Iso_index(i, 0);
    output[1] = Iso_index(i, 1);
    output[2] = Iso_index(i, 2);
    return output;
}

}

