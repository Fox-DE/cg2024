#include"RBF.h"
#include<cmath>

double RBFalgo::RadialBasis(double d)
{
    return pow(sqrt(d * d + r * r), u);
}

double RBFalgo::dist(ImVec2 p, ImVec2 q)
{
    return sqrt((p.x - q.x) * (p.x - q.x) + (p.y - q.y) * (p.y - q.y));
}

void RBFalgo::Init_A_b()
{
    // Init A and b
    const int n = start_points_.size();
    A(0, 0) = 1;
    A(0, 1) = 0;
    A(1, 0) = 0;
    A(1, 1) = 1;
    b[0] = 0;
    b[1] = 0;
    if (n == 0)
    {
        A(0, 0) = 1;
        A(0, 1) = 0;
        A(1, 0) = 0;
        A(1, 1) = 1;
        b[0] = 0;
        b[1] = 0;
    }
    else if (n == 1)
    {
        A(0, 0) = 1;
        A(0, 1) = 0;
        A(1, 0) = 0;
        A(1, 1) = 1;
        b[0] = end_points_[0].x - start_points_[0].x;
        b[1] = end_points_[0].y - start_points_[0].y;
    }
    else if (n == 2)
    {
        A(0, 1) = 0;
        A(1, 0) = 0;
        A(0, 0) = (start_points_[1].x - start_points_[0].x) == 0
                      ? 1
                      : (end_points_[1].x - end_points_[0].x) /
                            double(start_points_[1].x - start_points_[0].x);
        A(1, 1) = (start_points_[1].y - start_points_[0].y) == 0
                      ? 1
                      : (end_points_[1].y - end_points_[0].y) /
                            double(start_points_[1].y - start_points_[0].y);
        b[0] = end_points_[0].x - A(0, 0) * start_points_[0].x;
        b[1] = end_points_[1].y - A(1, 1) * start_points_[1].y;
    }
    else
    {
        Matrix<double, Dynamic, 3> X;
        Matrix<double, Dynamic, 2> Y;
        X.resize(n, 3);
        Y.resize(n, 2);
        Matrix<double, 3, 2> parameters;
        parameters.fill(0.0);
        for (int i = 0; i < n; i++)
        {
            X(i, 0) = 1;
            X(i, 1) = start_points_[i].x;
            X(i, 2) = start_points_[i].y;
            Y(i, 0) = end_points_[i].x;
            Y(i, 1) = end_points_[i].y;
        }
        parameters = (X.transpose() * X).inverse() * X.transpose() * Y;
        b[0] = parameters(0, 0);
        b[1] = parameters(0, 1);
        A(0, 0) = parameters(1, 0);
        A(1, 0) = parameters(1, 1);
        A(0, 1) = parameters(2, 0);
        A(1, 1) = parameters(2, 1);
    }
    
}

void RBFalgo::Init_Alpha()
{
    Matrix<double, Dynamic, Dynamic> Coff;
    Matrix<double, Dynamic, 2> Constant;
    size_t n = start_points_.size();
    Coff.resize(n, n);
    Constant.resize(n, 2);
    Alpha.resize(n, 2);
    if (n == 1)
        Alpha.fill(0.0);
    else
    {
        // Construct Coff and Constant preparing for the equation
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                Coff(i, j) =
                    RadialBasis(dist(start_points_[i], start_points_[j]));
            }
            Vector2d temp =
                Vector2d(end_points_[i].x, end_points_[i].y) - b -
                A * Vector2d(start_points_[i].x, start_points_[i].y);
            Constant(i, 0) = temp[0];
            Constant(i, 1) = temp[1];
        }

        Alpha = Coff.colPivHouseholderQr().solve(Constant);
    }
}

ImVec2 RBFalgo::Warping(ImVec2 p)
{    
    int n = start_points_.size();
    if (n== 0) 
        return p;
    Init_A_b();
    Init_Alpha();


    ImVec2 q;
    Vector2d p_input(double(p.x), double(p.y));
    Vector2d q_output(0, 0);
    q_output = A * p_input + b;
    for (int i = 0; i < n; i++)
    {
        q_output[0] += Alpha(i, 0) * RadialBasis(dist(p, start_points_[i]));
        q_output[1] += Alpha(i, 1) * RadialBasis(dist(p, start_points_[i]));
    }
    q.x = q_output[0];
    q.y = q_output[1];
    return q;
}