#pragma once
#include<vector>
#include<imgui.h>

class ImageWarping
{
public:
    void GetData(std::vector<ImVec2>start_points,std::vector<ImVec2>end_points);
    void Init();
    virtual ImVec2 Warping(ImVec2 p)=0;

public:
    std::vector<ImVec2> start_points_, end_points_;
};