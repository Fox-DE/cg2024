#include "comp_target_image.h"

#include <cmath>

namespace USTC_CG
{
using uchar = unsigned char;

CompTargetImage::CompTargetImage(
    const std::string& label,
    const std::string& filename)
    : ImageEditor(label, filename)
{
    if (data_)
        back_up_ = std::make_shared<Image>(*data_);
}

void CompTargetImage::draw()
{
    // Draw the image
    ImageEditor::draw();
    // Invisible button for interactions
    ImGui::SetCursorScreenPos(position_);
    ImGui::InvisibleButton(
        label_.c_str(),
        ImVec2(
            static_cast<float>(image_width_),
            static_cast<float>(image_height_)),
        ImGuiButtonFlags_MouseButtonLeft);
    bool is_hovered_ = ImGui::IsItemHovered();
    // When the mouse is clicked or moving, we would adapt clone function to
    // copy the selected region to the target.
    ImGuiIO& io = ImGui::GetIO();
    if (is_hovered_ && ImGui::IsMouseClicked(ImGuiMouseButton_Left))
    {
        edit_status_ = true;
        mouse_position_ =
            ImVec2(io.MousePos.x - position_.x, io.MousePos.y - position_.y);
        clone();
    }
    if (edit_status_)
    {
        mouse_position_ =
            ImVec2(io.MousePos.x - position_.x, io.MousePos.y - position_.y);
        if (flag_realtime_updating)
            clone();
        if (!ImGui::IsMouseDown(ImGuiMouseButton_Left))
        {
            edit_status_ = false;
        }
    }
}

void CompTargetImage::set_source(std::shared_ptr<CompSourceImage> source)
{
    source_image_ = source;
}

void CompTargetImage::set_realtime(bool flag)
{
    flag_realtime_updating = flag;
}

void CompTargetImage::restore()
{
    *data_ = *back_up_;
    update();
}

void CompTargetImage::set_paste()
{
    clone_type_ = kPaste;
}

void CompTargetImage::set_seamless()
{
    clone_type_ = kSeamless;
}

void CompTargetImage::clone()
{
    // The implementation of different types of cloning
    // HW3_TODO: In this function, you should at least implement the "seamless"
    // cloning labeled by `clone_type_ ==kSeamless`.
    //
    // The realtime updating (update when the mouse is moving) is only available
    // when the checkboard is selected. It is required to improve the efficiency
    // of your seamless cloning to achieve realtime editing. (Use decomposition
    // of sparse matrix before solve the linear system)
    std::shared_ptr<Image> mask = source_image_->get_region();
    int width = mask->width();
    int height = mask->height();
    switch (clone_type_)
    {
        case USTC_CG::CompTargetImage::kDefault: break;
        case USTC_CG::CompTargetImage::kPaste:
        {
            restore();

            for (int i = 0; i < mask->width(); ++i)
            {
                for (int j = 0; j < mask->height(); ++j)
                {
                    int tar_x =
                        static_cast<int>(mouse_position_.x) + i -
                        static_cast<int>(source_image_->get_position().x);
                    int tar_y =
                        static_cast<int>(mouse_position_.y) + j -
                        static_cast<int>(source_image_->get_position().y);
                    if (0 <= tar_x && tar_x < image_width_ && 0 <= tar_y &&
                        tar_y < image_height_ && mask->get_pixel(i, j)[0] > 0)
                    {
                        data_->set_pixel(
                            tar_x,
                            tar_y,
                            source_image_->get_data()->get_pixel(i, j));
                    }
                }
            }
            break;
        }
        case USTC_CG::CompTargetImage::kSeamless:
        {
            // You should delete this block and implement your own seamless
            // cloning. For each pixel in the selected region, calculate the
            // final RGB color by solving Poisson Equations.
            restore();

            /*for (int i = 0; i < mask->width(); ++i)
            {
                for (int j = 0; j < mask->height(); ++j)
                {
                    int tar_x =
                        static_cast<int>(mouse_position_.x) + i -
                        static_cast<int>(source_image_->get_position().x);
                    int tar_y =
                        static_cast<int>(mouse_position_.y) + j -
                        static_cast<int>(source_image_->get_position().y);
                    if (0 <= tar_x && tar_x < image_width_ && 0 <= tar_y &&
                        tar_y < image_height_ && mask->get_pixel(i, j)[0] > 0)
                    {
                        data_->set_pixel(
                            tar_x,
                            tar_y,
                            source_image_->get_data()->get_pixel(i, j));
                    }
                }
            }*/
            int point_count = source_image_->get_point_count();
            Eigen::MatrixXd b(point_count,3);
            //Eigen::SparseMatrix<double, Eigen::RowMajor> b(point_count, 3);
            for (int i = 0; i < mask->width(); ++i)
            {
                for (int j = 0; j < mask->height(); ++j)
                {
                    if (mask->get_pixel(i, j)[0] != 0)//this pixel is interior or boundary
                    {
                        int tar_x =
                            static_cast<int>(mouse_position_.x) + i -
                            static_cast<int>(source_image_->get_position().x);
                        int tar_y =
                            static_cast<int>(mouse_position_.y) + j -
                            static_cast<int>(source_image_->get_position().y);
                        //tar_x,y为对应的target_image中的点
                        int index = source_image_->get_index(j, i);
                        b(index, 0) = 0;
                        b(index, 1) = 0;
                        b(index, 2) = 0;
                        int row_neighbor[4] = { 1, -1, 0, 0 };
                        int col_neighbor[4] = { 0, 0, -1, 1 };
                        auto color_tmp = source_image_->get_data()->get_pixel(
                            i,j);
                        //color_tmp为源图像该点处的rgb值
                        for (int k = 0; k < 4; k++)
                        {
                            int tar_tmp_x = tar_x+row_neighbor[k];
                            int tar_tmp_y = tar_y+col_neighbor[k];
                            int row_tmp = j + col_neighbor[k];
                            int col_tmp = i + row_neighbor[k];
                            int index_tmp =
                                source_image_->get_index(row_tmp, col_tmp);
                            if ((row_tmp < 0) || (col_tmp < 0) ||
                                (row_tmp > height - 1) || (col_tmp > width - 1))
                            {
                                continue;
                            }//越界则跳过

                            auto color = source_image_->get_data()->get_pixel(
                                col_tmp, row_tmp);
                            //color即为该邻点对应的rgb值

                            b(index, 0) += (int)color_tmp[0] - (int)color[0];
                            b(index, 1) += (int)color_tmp[1] - (int)color[1];
                            b(index, 2) += (int)color_tmp[2] - (int)color[2];

                            if (mask->get_pixel(col_tmp, row_tmp)[0] == 100)
                                //即该点为边界项
                            {
                                auto color =
                                    data_->get_pixel(tar_tmp_x, tar_tmp_y);
                                b(index, 0) += (int)color[0];
                                b(index, 1) += (int)color[1];
                                b(index, 2) += (int)color[2];

                            }
                            
                                /*auto color = data_->get_pixel(tar_tmp_x,tar_tmp_y);
                                b(index, 0) += (int)color[0];
                                b(index, 1) += (int)color[1];
                                b(index, 2) += (int)color[2];*/
                            
                        }

                    }
                }
            }
            Eigen::MatrixXd f = solver.solve(b);

            for (int i = 0; i < mask->width(); ++i)
            {
                for (int j = 0; j < mask->height(); ++j)
                {
                    int tar_x =
                        static_cast<int>(mouse_position_.x) + i -
                        static_cast<int>(source_image_->get_position().x);
                    int tar_y =
                        static_cast<int>(mouse_position_.y) + j -
                        static_cast<int>(source_image_->get_position().y);
                    if (0 <= tar_x && tar_x < image_width_ && 0 <= tar_y &&
                        tar_y < image_height_ && mask->get_pixel(i, j)[0] != 0)
                    {
                        int index_tmp = source_image_->get_index(j, i);
                        double color_check[3];
                        color_check[0] = f(index_tmp, 0);
                        color_check[1] = f(index_tmp, 1);
                        color_check[2] = f(index_tmp, 2);
                        for (int i_ = 0; i_ < 3; i_++)
                        {
                            if (color_check[i_] < 0)
                                color_check[i_] = 0;
                            if (color_check[i_] > 255)
                                color_check[i_] = 255;
                        }
                        unsigned char color[3];
                        color[0] = color_check[0];
                            //;  >255?f(index_tmp,0):255;
                        color[1] = color_check[1];  
                        //> 255?  : 255;
                        color[2] = color_check[2];  
                        // > 255? f(index_tmp, 2): 255;
                        
                        data_->set_pixel(
                            tar_x, tar_y, { color[0], color[1], color[2] });
                    }
                }
            }




            break;
        }
        default: break;
    }

    update();
}

void CompTargetImage::preDecomposition()
{
    std::shared_ptr<Image> mask = source_image_->get_region();
    int width = mask->width();
    int height = mask->height();
    //width and height are the original size of the source_image
    //if mask->get_pixel(i, j)[0]==255, (i,j) is in the point;
    int point_count = source_image_->get_point_count();
    //total point number to be Calculated
    Eigen::SparseMatrix<double, Eigen::RowMajor> A(point_count, point_count);
    //RowMajor is more efficient for this operation
    for (int i = 0; i < mask->height(); i++)
    {
        for (int j = 0; j < mask->width();j++)
        {
            if (mask->get_pixel(j, i)[0] != 0)//this means interior or boundary
                //不等于外部，边界或者内部，边界亦为需要求解的点
            {
                A.insert(
                    source_image_->get_index(i,j),
                    source_image_->get_index(i,j)) = 4;
                if (i == 0 || i == height - 1)
                {
                    A.coeffRef(
                        source_image_->get_index(i,j),
                        source_image_->get_index(i,j))--;
                }
                if (j == 0 || j == width - 1)
                {
                    A.coeffRef(
                        source_image_->get_index(i,j),
                        source_image_->get_index(i,j))--;
                }
                if (mask->get_pixel(j, i)[0] == 100)
                {
                    A.coeffRef(
                        source_image_->get_index(i, j),
                        source_image_->get_index(i, j))--;
                }
                

                // operate the border points
                int row_neighbor[4] = { 1, -1, 0, 0 };
                int col_neighbor[4] = { 0, 0, -1, 1 };
                for (int k = 0; k < 4; k++)
                {
                    int row_tmp = i + row_neighbor[k];
                    int col_tmp = j + col_neighbor[k];
                    int index_tmp = source_image_->get_index(row_tmp, col_tmp);
                    //get_index先行后列
                    //get_pixel先列后行
                    if ((row_tmp < 0) || (col_tmp < 0) ||
                        (row_tmp > height - 1) || (col_tmp > width - 1))
                    {
                        continue;
                    }
                    if (mask->get_pixel(col_tmp, row_tmp)[0] == 255)
                        //255即为内部
                    {
                        A.insert(
                            source_image_->get_index(i,j),
                            source_image_->get_index( row_tmp,col_tmp)) = -1;
                        
                    }
                }
                /* A.insert(
                   source_image_->get_index(i, j),
                    source_image_->get_index(i, j)) = count;*/
            }
        }
    }
    Eigen::MatrixXd A_dense = Eigen::MatrixXd(A);
    solver.compute(A);

}

}  // namespace USTC_CG