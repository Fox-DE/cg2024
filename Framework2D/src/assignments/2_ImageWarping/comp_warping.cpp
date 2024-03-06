#include "comp_warping.h"

#include <cmath>

#include"ImageWarping.h"

#include"IDW.h"

#include"annoy/annoylib.h"
#include"annoy/kissrandom.h"
#include"annoy/mman.h"

namespace USTC_CG
{
using uchar = unsigned char;

CompWarping::CompWarping(const std::string& label, const std::string& filename)
    : ImageEditor(label, filename)
{
    if (data_)
        back_up_ = std::make_shared<Image>(*data_);
}

void CompWarping::draw()
{
    // Draw the image
    ImageEditor::draw();
    // Draw the canvas
    if (flag_enable_selecting_points_)
        select_points();
}

void CompWarping::invert()
{
    for (int i = 0; i < data_->width(); ++i)
    {
        for (int j = 0; j < data_->height(); ++j)
        {
            const auto color = data_->get_pixel(i, j);
            data_->set_pixel(
                i,
                j,
                { static_cast<uchar>(255 - color[0]),
                  static_cast<uchar>(255 - color[1]),
                  static_cast<uchar>(255 - color[2]) });
        }
    }
    // After change the image, we should reload the image data to the renderer
    update();
}
void CompWarping::mirror(bool is_horizontal, bool is_vertical)
{
    Image image_tmp(*data_);
    int width = data_->width();
    int height = data_->height();

    if (is_horizontal)
    {
        if (is_vertical)
        {
            for (int i = 0; i < width; ++i)
            {
                for (int j = 0; j < height; ++j)
                {
                    data_->set_pixel(
                        i,
                        j,
                        image_tmp.get_pixel(width - 1 - i, height - 1 - j));
                }
            }
        }
        else
        {
            for (int i = 0; i < width; ++i)
            {
                for (int j = 0; j < height; ++j)
                {
                    data_->set_pixel(
                        i, j, image_tmp.get_pixel(width - 1 - i, j));
                }
            }
        }
    }
    else
    {
        if (is_vertical)
        {
            for (int i = 0; i < width; ++i)
            {
                for (int j = 0; j < height; ++j)
                {
                    data_->set_pixel(
                        i, j, image_tmp.get_pixel(i, height - 1 - j));
                }
            }
        }
    }

    // After change the image, we should reload the image data to the renderer
    update();
}
void CompWarping::gray_scale()
{
    for (int i = 0; i < data_->width(); ++i)
    {
        for (int j = 0; j < data_->height(); ++j)
        {
            const auto color = data_->get_pixel(i, j);
            uchar gray_value = (color[0] + color[1] + color[2]) / 3;
            data_->set_pixel(i, j, { gray_value, gray_value, gray_value });
        }
    }
    // After change the image, we should reload the image data to the renderer
    update();
}

void CompWarping::IDW()
{
   
    IDWalgo IDW_;
    IDW_.GetData(start_points_, end_points_);

    Image Image_tmp(*data_);
    int width = data_->width();
    int height = data_->height();
    
    ImVec2 q_tmp,p_tmp;
    //initialize the graph
    for (int y = 0; y < data_->height(); ++y)
    {
        for (int x = 0; x < data_->width(); ++x)
        {
            data_->set_pixel(x, y, { 255, 255, 255 });
        }
    }

     for (size_t i = 0; i < width; i++)
    {
        for (int j = 0; j < height; j++)
        {
            p_tmp.x = i;
            p_tmp.y = j;
            q_tmp = IDW_.Warping(p_tmp);


            if (q_tmp.x < width && q_tmp.y < height &&
                q_tmp.x >= 0 && q_tmp.y >= 0)
                data_->set_pixel(q_tmp.x,q_tmp.y, Image_tmp.get_pixel(i, j));
        }
    }
    update();
}

void CompWarping::RBF()
{
    RBFalgo RBF_;
    RBF_.GetData(start_points_, end_points_);

    Image Image_tmp(*data_);
    int width = data_->width();
    int height = data_->height();

    ImVec2 q_tmp, p_tmp;
    // initialize the graph
    for (int y = 0; y < data_->height(); ++y)
    {
        for (int x = 0; x < data_->width(); ++x)
        {
            data_->set_pixel(x, y, { 255, 255, 255 });
        }
    }

    for (size_t i = 0; i < width; i++)
    {
        for (int j = 0; j < height; j++)
        {
            p_tmp.x = i;
            p_tmp.y = j;
            q_tmp = RBF_.Warping(p_tmp);

            if (q_tmp.x < width && q_tmp.y < height && q_tmp.x >= 0 &&
                q_tmp.y >= 0)
                data_->set_pixel(q_tmp.x, q_tmp.y, Image_tmp.get_pixel(i, j));
        }
    }
    update();
}
void CompWarping::warping()
{
    // HW2_TODO: You should implement your own warping function that interpolate
    // the selected points.
    // You can design a class for such warping operations, utilizing the
    // encapsulation, inheritance, and polymorphism features of C++. More files
    // like "*.h", "*.cpp" can be added to this directory or anywhere you like.

    // Create a new image to store the result
    Image warped_image(*data_);
    // Initialize the color of result image
    for (int y = 0; y < data_->height(); ++y)
    {
        for (int x = 0; x < data_->width(); ++x)
        {
            warped_image.set_pixel(x, y, { 0, 0, 0 });
        }
    }

    // Example: (simplified) "fish-eye" warping
    // For each (x, y) from the input image, the "fish-eye" warping transfer it
    // to (x', y') in the new image:
    // Note: For this transformation ("fish-eye" warping), one can also
    // calculate the inverse (x', y') -> (x, y) to fill in the "gaps".
    for (int y = 0; y < data_->height(); ++y)
    {
        for (int x = 0; x < data_->width(); ++x)
        {
            // Apply warping function to (x, y), and we can get (x', y')
            auto [new_x, new_y] =
                fisheye_warping(x, y, data_->width(), data_->height());
            // Copy the color from the original image to the result image
            if (new_x >= 0 && new_x < data_->width() && new_y >= 0 &&
                new_y < data_->height())
            {
                std::vector<unsigned char> pixel = data_->get_pixel(x, y);
                warped_image.set_pixel(new_x, new_y, pixel);
            }
        }
    }

    *data_ = std::move(warped_image);
    update();
}
void CompWarping::restore()
{
    *data_ = *back_up_;
    update();
}
void CompWarping::enable_selecting(bool flag)
{
    flag_enable_selecting_points_ = flag;
}
void CompWarping::select_points()
{
    /// Invisible button over the canvas to capture mouse interactions.
    ImGui::SetCursorScreenPos(position_);
    ImGui::InvisibleButton(
        label_.c_str(),
        ImVec2(
            static_cast<float>(image_width_),
            static_cast<float>(image_height_)),
        ImGuiButtonFlags_MouseButtonLeft);
    // Record the current status of the invisible button
    bool is_hovered_ = ImGui::IsItemHovered();
    // Selections
    ImGuiIO& io = ImGui::GetIO();
    if (is_hovered_ && ImGui::IsMouseClicked(ImGuiMouseButton_Left))
    {
        draw_status_ = true;
        start_ = end_ =
            ImVec2(io.MousePos.x - position_.x, io.MousePos.y - position_.y);
    }
    if (draw_status_)
    {
        end_ = ImVec2(io.MousePos.x - position_.x, io.MousePos.y - position_.y);
        if (!ImGui::IsMouseDown(ImGuiMouseButton_Left))
        {
            start_points_.push_back(start_);
            end_points_.push_back(end_);
            draw_status_ = false;
        }
    }
    // Visualization
    auto draw_list = ImGui::GetWindowDrawList();
    for (size_t i = 0; i < start_points_.size(); ++i)
    {
        ImVec2 s(
            start_points_[i].x + position_.x, start_points_[i].y + position_.y);
        ImVec2 e(
            end_points_[i].x + position_.x, end_points_[i].y + position_.y);
        draw_list->AddLine(s, e, IM_COL32(255, 0, 0, 255), 2.0f);
        draw_list->AddCircleFilled(s, 4.0f, IM_COL32(0, 0, 255, 255));
        draw_list->AddCircleFilled(e, 4.0f, IM_COL32(0, 255, 0, 255));
    }
    if (draw_status_)
    {
        ImVec2 s(start_.x + position_.x, start_.y + position_.y);
        ImVec2 e(end_.x + position_.x, end_.y + position_.y);
        draw_list->AddLine(s, e, IM_COL32(255, 0, 0, 255), 2.0f);
        draw_list->AddCircleFilled(s, 4.0f, IM_COL32(0, 0, 255, 255));
    }
}
void CompWarping::init_selections()
{
    start_points_.clear();
    end_points_.clear();
}

std::pair<int, int>
CompWarping::fisheye_warping(int x, int y, int width, int height)
{
    float center_x = width / 2.0f;
    float center_y = height / 2.0f;
    float dx = x - center_x;
    float dy = y - center_y;
    float distance = std::sqrt(dx * dx + dy * dy);

    // Simple non-linear transformation r -> r' = f(r)
    float new_distance = std::sqrt(distance) * 10;

    if (distance == 0)
    {
        return { static_cast<int>(center_x), static_cast<int>(center_y) };
    }
    // (x', y')
    float ratio = new_distance / distance;
    int new_x = static_cast<int>(center_x + dx * ratio);
    int new_y = static_cast<int>(center_y + dy * ratio);

    return { new_x, new_y };
}

void CompWarping::FixImage()
{
    int search_num = 16;
    int f = 2;
    Annoy::AnnoyIndex<
        int,
        float,
        Annoy::Euclidean,
        Annoy::Kiss32Random,
        Annoy::AnnoyIndexSingleThreadedBuildPolicy>
        index(f);
    int width = data_->width();
    int height = data_->height();
    for (int j = 0; j < height; j++)
    {
        for (int i = 0; i < width; i++)
        {
            std::vector<float> vec(f);
            vec[0] = i;//ºá×ø±ê
            vec[1] = j;//×Ý×ø±ê
            index.add_item(j * width + i, vec.data());
        }
    }
    index.build(10);
    std::vector<int> closet_items; 
    std::vector<float> distances;
    for (int j = 0; j < height; j++)
    {
        for (int i = 0; i < width; i++)
        {
            const auto color = data_->get_pixel(i, j);
            if (color[0] == 255 && color[1] == 255 && color[2] == 255)
            {
                closet_items.clear();
                distances.clear();
                index.get_nns_by_item(
                    j * width + i, search_num, -1, &closet_items,&distances);
                
                const auto color_tmp = data_->get_pixel(
                    closet_items[0] - int(closet_items[0] / width)*width,
                    int(closet_items[0] / width));
                int value_r = (int)color_tmp[0];
                int value_g = (int)color_tmp[1];
                int value_b = (int)color_tmp[2];
                
                for (int k = 1; k < closet_items.size(); k++)
                {
                    const auto color_tmp = data_->get_pixel(
                        closet_items[k] - int(closet_items[k] / width) * width,
                        int(closet_items[k] / width));
                    value_r = value_r + (int)color_tmp[0];
                    value_g = value_g + (int)color_tmp[1];
                    value_b = value_b + (int)color_tmp[2];
                }
                value_r = (int)value_r / search_num;
                value_g = (int)value_g / search_num;
                value_b = (int)value_b / search_num;
                data_->set_pixel(
                    i, j, { (uchar)value_r, (uchar)value_g, (uchar)value_b });
                
            }
        }
    }
    update();
}

}  // namespace USTC_CG