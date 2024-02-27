#pragma once

#include "shape.h"
#include<vector>

namespace USTC_CG
{
class Freehand : public Shape
{
   public:
    Freehand() = default;

    Freehand(ImVec2 point);

    virtual ~Freehand() = default;

    void draw(const Config& config) const override;

    void update(float x, float y) override;

   private:
    float start_point_x_, start_point_y_, end_point_x_, end_point_y_;
    std::vector<ImVec2> Trace;
};
}  // namespace USTC_CG
