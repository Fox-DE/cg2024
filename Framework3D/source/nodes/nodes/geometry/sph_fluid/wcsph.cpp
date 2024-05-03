#include "wcsph.h"
#include <iostream>
using namespace Eigen;

namespace USTC_CG::node_sph_fluid {

WCSPH::WCSPH(const MatrixXd& X, const Vector3d& box_min, const Vector3d& box_max)
    : SPHBase(X, box_min, box_max)
{
}

void WCSPH::compute_density()
{
	// -------------------------------------------------------------
	// (HW TODO) Implement the density computation
    // You can also compute pressure in this function 
	// -------------------------------------------------------------
    //override this function to compute density and pressure at the same time
    for (auto& p : ps_.particles()) {
            
        auto density = ps_.mass() * W_zero(ps_.h());       
               
        for (auto& q : p->neighbors()) {
                density += ps_.mass() * W(p->x() - q->x(), ps_.h());
               
        }
        p->density_ = density;
        double pi = stiffness_ * (pow(density / ps_.density0(), exponent_) - 1.0);
        p->pressure_ = std::max(pi, 0.0);
    }


}

void WCSPH::step()
{
    TIC(step)
    // -------------------------------------------------------------
    // (HW TODO) Follow the instruction in documents and PPT,
    // implement the pipeline of fluid simulation
    // -------------------------------------------------------------

    // Search neighbors, compute density, advect, solve pressure acceleration, etc.
    ps_.assign_particles_to_cells();
    ps_.searchNeighbors();
    compute_density();
    this->compute_non_pressure_acceleration();
    for (auto& p : ps_.particles()) {
        p->vel_ += dt() * p->acceleration();
    }
    this->compute_pressure_gradient_acceleration();
    for (auto& p : ps_.particles()) {
        p->vel_ += dt() * p->acceleration();
        p->x() += p->vel() * dt();
        check_collision(p);
    }
    advect();



    TOC(step)
}
}  // namespace USTC_CG::node_sph_fluid