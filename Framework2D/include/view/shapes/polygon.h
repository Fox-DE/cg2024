#pragma once

#include "shape.h"
#include<vector>

namespace USTC_CG
{
class Polygon : public Shape
{
   public:
    Polygon() = default;

    Polygon(ImVec2 point);
    //i implies the status
    //1 means add the first vertex
    //2 means add new vertex
    //3 means add the last vertex

    void SetStatus(bool status) override;

    virtual ~Polygon() = default;

    void draw(const Config& config) const override;

    void update(float x, float y) override;

    void AddVertex(float x,float y) override;

   private:
    float start_point_x_, start_point_y_, end_point_x_, end_point_y_;
    std::vector<ImVec2> Vertex;
    bool Closure_Status=false;
};
}  // namespace USTC_CG
