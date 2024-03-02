#include <imgui.h>

#include <cmath>

#include "view/shapes/freehand.h"

namespace USTC_CG
{

Freehand::Freehand(ImVec2 point)
{

    Trace.clear();
    Trace.push_back(point);
}

void Freehand::draw(const USTC_CG::Shape::Config& config) const
{
    ImDrawList* draw_list = ImGui::GetWindowDrawList();


    for (int i = 0; i < Trace.size() - 1; i++)
    {
        draw_list->AddLine(
            ImVec2(
                config.bias[0] + Trace[i].x, config.bias[1] + Trace[i].y),
            ImVec2(
                config.bias[0] + Trace[i + 1].x, config.bias[1] + Trace[i+1].y),
            IM_COL32(
                config.line_color[0],
                config.line_color[1],
                config.line_color[2],
                config.line_color[3]),
            config.line_thickness);
    }
    //draw lines between points, which makes the figure seem continuous
}

void Freehand::update(float x, float y)
{
    ImVec2 point;
    point.x = x;
    point.y = y;
    Trace.push_back(point);
}

}  // namespace USTC_CG
