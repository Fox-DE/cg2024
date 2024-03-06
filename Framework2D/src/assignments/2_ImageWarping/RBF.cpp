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
    int n = start_points_.size();
    A(0, 0) = 1;
    A(0, 1) = 0;
    A(1, 0) = 0;
    A(1, 1) = 1;
    b[0] = 0;
    b[1] = 0;
    
    if (n == 1)
    {
        
        b[0] = end_points_[0].x - start_points_[0].x;
        b[1] = end_points_[0].y - start_points_[0].y;
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

void RBFalgo::Init_r()
{ 
    int n = start_points_.size();
    if (n >= 2)
    {
        double temp = dist(start_points_[0], start_points_[1]);
        double tmp;
        for (int i = 0; i < n; i++)
        {
            for (int j = i+1; j < n; j++)
            {
                tmp = dist(start_points_[i], start_points_[j]);
                if (tmp < temp) temp = tmp;
            }
        }
        r = temp;
    }
}

ImVec2 RBFalgo::Warping(ImVec2 p)
{    
    int n = start_points_.size();
    if (n== 0) 
        return p; 
    Init_r();
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