#include <imgui.h>

#include <cmath>

#include "view/shapes/polygon.h"

namespace USTC_CG
{

Polygon::Polygon(ImVec2 point)
    : start_point_x_(point.x),
      start_point_y_(point.y),
      end_point_x_(point.x),
      end_point_y_(point.y)
{
    Vertex.push_back(point);
}

void Polygon::draw(const USTC_CG::Shape::Config& config) const
{
    ImDrawList* draw_list = ImGui::GetWindowDrawList();


    for (int i = 0; i < Vertex.size() - 1; i++)
    {
        draw_list->AddLine(
            ImVec2(
                config.bias[0] + Vertex[i].x, config.bias[1] + Vertex[i].y),
            ImVec2(
                config.bias[0] + Vertex[i + 1].x, config.bias[1] + Vertex[i+1].y),
            IM_COL32(
                config.line_color[0],
                config.line_color[1],
                config.line_color[2],
                config.line_color[3]),
            config.line_thickness);
    }
    draw_list->AddLine(
        ImVec2(
            config.bias[0] + Vertex[Vertex.size() - 1].x,
            config.bias[1] + Vertex[Vertex.size() - 1].y),
        ImVec2(
            config.bias[0] + end_point_x_, config.bias[1] + end_point_y_),
        IM_COL32(
            config.line_color[0],
            config.line_color[1],
            config.line_color[2],
            config.line_color[3]),
        config.line_thickness);
    if (Closure_Status)
    {
        draw_list->AddLine(
            ImVec2(config.bias[0] + Vertex[0].x, config.bias[1] + Vertex[0].y),
            ImVec2(
                config.bias[0] + Vertex[Vertex.size() - 1].x,
                config.bias[1] + Vertex[Vertex.size() - 1].y),
            IM_COL32(
                config.line_color[0],
                config.line_color[1],
                config.line_color[2],
                config.line_color[3]),
            config.line_thickness);
    }
    //draw the line between the first vertex and the last one, make the polygon closed
}

void Polygon::update(float x, float y)
{
    end_point_x_ = x;
    end_point_y_ = y;
}

void Polygon::SetStatus(bool status)
{
    if (status)
    {
        Closure_Status = true;
    }
    else
    {
        Closure_Status = false;
    }
}

void Polygon::AddVertex(float x,float y)
{
    ImVec2 point;
    point.x = x;
    point.y = y;
    Vertex.push_back(point);
}

}  // namespace USTC_CG
