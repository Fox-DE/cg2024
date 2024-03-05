#include"ImageWarping.h"
#include<cassert>

void ImageWarping::Init()
{
    start_points_.clear();
    end_points_.clear();
}

void ImageWarping::GetData(
    std::vector<ImVec2> start_points,
    std::vector<ImVec2> end_points)
{
    Init();
    assert(start_points.size() == end_points.size());
    for (int i = 0; i < start_points.size(); i++)
    {
        start_points_.push_back(start_points[i]);
        end_points_.push_back(end_points[i]);
    }
}


