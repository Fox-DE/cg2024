#include <cmath>
#include "GCore/Components/MeshOperand.h"
#include "util_openmesh_bind.h"
#include "Eigen/Dense"
#include <vector>

double cos2cot(double input);

namespace USTC_CG {
	class init_Parameterization
	{
	public:
         void init(std::shared_ptr<PolyMesh> input);
         double get_cot(int i, int j);
         void Isometric_parameterization();
         std::vector<double> get_iso(int i);
	private:
         std::shared_ptr<PolyMesh> original_mesh;
         Eigen::MatrixXd cot_index;
         Eigen::MatrixXd Iso_index;
	};

}