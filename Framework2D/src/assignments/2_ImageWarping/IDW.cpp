#include"IDW.h"

ImVec2 IDWalgo::Warping(ImVec2 p)
{
    double x = (double)p.x;
    double y = (double)p.y;
    std::vector<double> pq_tmp(2);
    std::vector<int> pq1(2);  // the final output data
    int n = start_points_.size();
    sigma.clear();
    w.clear();
    double tmp = NULL;
    double xi, yi;
    for (int i = 0; i < n; i++)
    {
        xi = start_points_[i].x * 1.0;
        yi = start_points_[i].y * 1.0;  //(xi,yi) belongs to the P points
        tmp = 1 / ((x - xi) * (x - xi) + (y - yi) * (y - yi));
        sigma.push_back(tmp);
        tmp = NULL;
    }
    tmp = 0;
    for (int i = 0; i < n; i++)
    {
        tmp += sigma[i];
    }
    for (int i = 0; i < n; i++)
    {
        w.push_back(sigma[i] / tmp);
    }
    pq_tmp[0] = 0;
    pq_tmp[1] = 0;
    std::vector<double> D;
    for (int i = 0; i < n; i++)
    {
        D = MinimumEnergy(i);
        pq_tmp[0] += w[i] * (end_points_[i].x + D[0] * (x - start_points_[i].x) +
                    D[1] * (y - start_points_[i].y));
        pq_tmp[1] += w[i] * (end_points_[i].y + D[2] * (x - start_points_[i].x) +
                    D[3] * (y - start_points_[i].y));
    }
    pq1[0] = (int)pq_tmp[0];
    pq1[1] = (int)pq_tmp[1];
    ImVec2 q;
    q.x = pq1[0];
    q.y = pq1[1];
    return q;

}

std::vector<double> IDWalgo::MinimumEnergy(int i)
{
    //minimum the energy to get the matrix D
    std::vector<double> linearD(4);
    int n = start_points_.size();
    // std::cout << n << endl;
    std::vector<double> Wij(n, 0);
    std::vector<double> t1(n, 0);
    std::vector<double> t2(n, 0);
    std::vector<double> s1(n, 0);
    std::vector<double> s2(n, 0);
    for (size_t j = 0; j < n; j++)
    {
        if (j != i)
        {
            double tmp =
                ((start_points_[j].x - start_points_[i].x) *
                     (start_points_[j].x - start_points_[i].x) +
                 (start_points_[j].y - start_points_[i].y) *
                     (start_points_[j].y - start_points_[i].y));
            Wij[j] = 1 / tmp;
            t1[j] = start_points_[j].x - start_points_[i].x;
            t2[j] = start_points_[j].y - start_points_[i].y;
            s1[j] = end_points_[i].x - end_points_[j].x;
            s2[j] = end_points_[i].y - end_points_[j].y;
        }
    }
    double A = 0, B = 0, C = 0, M = 0, N = 0, P = 0, Q = 0;
    for (size_t j = 0; j < n; j++)
    {
        if (j != i)
        {
            A += Wij[j] * t1[j] * t1[j];
            B += Wij[j] * t2[j] * t2[j];
            C += Wij[j] * t2[j] * t1[j];
            M += Wij[j] * s1[j] * t1[j];
            N += Wij[j] * s1[j] * t2[j];
            P += Wij[j] * s2[j] * t1[j];
            Q += Wij[j] * s2[j] * t2[j];
        }
    }
    // std::cout << A << " " << B << " " << C << " " << " " << M << " " << N <<
    // " " << P << " " << Q << " " << endl;
    if (fabs(A * B - C * C) > 1e-10)
    {
        linearD[0] = (C * N - B * M) / (A * B - C * C);
        linearD[1] = (C * M - A * N) / (A * B - C * C);
        linearD[2] = (C * Q - B * P) / (A * B - C * C);
        linearD[3] = (C * P - A * Q) / (A * B - C * C);
    }
    else  // if the denominator is too small, set D to Eye
    {
        linearD[0] = 1;
        linearD[1] = 0;
        linearD[2] = 0;
        linearD[3] = 1;
    }
    return linearD;
}