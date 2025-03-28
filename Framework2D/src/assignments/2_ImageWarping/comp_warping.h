#pragma once

#include "view/comp_image.h"
#include"ImageWarping.h"
#include"IDW.h"
#include"RBF.h"
#include"annoy/annoylib.h"
#include<vector>

namespace USTC_CG
{
// Image component for warping and other functions
class CompWarping : public ImageEditor
{
   public:
    explicit CompWarping(const std::string& label, const std::string& filename);
    virtual ~CompWarping() noexcept = default;

    void draw() override;

    // Simple edit functions
    void invert();
    void mirror(bool is_horizontal, bool is_vertical);
    void gray_scale();
    void warping();
    void restore();
    void IDW();
    void RBF();
    void FixImage();

    // Point selecting interaction
    void enable_selecting(bool flag);
    void select_points();
    void init_selections();

    void init_FixedList();
   private:
    // Store the original image data
    std::shared_ptr<Image> back_up_;
    // The selected point couples for image warping
    std::vector<ImVec2> start_points_, end_points_;

    ImVec2 start_, end_;
    bool flag_enable_selecting_points_ = false;
    bool draw_status_ = false;

   private:
    // A simple "fish-eye" warping function
    std::pair<int, int> fisheye_warping(int x, int y, int width, int height);

    private:
    std::vector<std::vector<bool>> ToBeFixed;
     bool FixedListStatus = false;
};

}  // namespace USTC_CG