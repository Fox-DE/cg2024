#pragma once
#include"ImageWarping.h"
#include"Eigen/Dense"
#include<vector>

using namespace Eigen;

class RBFalgo : public ImageWarping
{
   public:
    ImVec2 Warping(ImVec2 p) override;
    double RadialBasis(double d);
    double dist(ImVec2 p, ImVec2 q);
    void Init_A_b();
    void Init_Alpha();
    void Init_r();

   private:
    double u=1;
    Matrix2d A;
    Vector2d b;
    double r = 1;
    Matrix<double, Dynamic, 2> Alpha;
};
