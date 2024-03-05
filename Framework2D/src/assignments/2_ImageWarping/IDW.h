#pragma once
#include"ImageWarping.h"

class IDWalgo :public ImageWarping
{
   public:
        ImVec2 Warping(ImVec2 p) override;
    std::vector<double> MinimumEnergy(int i);
   private:
    std::vector<double> sigma;
    std::vector<double> w;
};