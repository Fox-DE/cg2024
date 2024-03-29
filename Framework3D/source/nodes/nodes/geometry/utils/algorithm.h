#include <cmath>
#include "GCore/Components/MeshOperand.h"
#include "util_openmesh_bind.h"
#include "Eigen/Dense"
#include "Eigen/Sparse"
#include <vector>

using namespace Eigen;
namespace USTC_CG {
	class ARAP
	{
	public:
        void get_data(
                std::shared_ptr<PolyMesh> original_mesh_in,
                std::shared_ptr<PolyMesh> halfedge_mesh_in);
        void compute();
        void set_itr(int i);
        std::shared_ptr<PolyMesh> output_mesh();
        
        //import initial coordinate       

	private: 
         void init();
         void build_A();
         void import_u_data();
         void update_u_data();
         int itr = 1;
         std::shared_ptr<PolyMesh> original_mesh;
         std::shared_ptr<PolyMesh> halfedge_mesh;
         int nV;
         int nF;
         MatrixXd area;
         MatrixXd angle;
         //edge_s means edge stored
         MatrixXd edge_s;
         MatrixXi mat_index;
         SparseMatrix<double> mat_A;
         Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> solver;
         MatrixXd u;
	};

}