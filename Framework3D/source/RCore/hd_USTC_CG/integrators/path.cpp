#include "path.h"

#include <random>

#include "surfaceInteraction.h"
USTC_CG_NAMESPACE_OPEN_SCOPE
using namespace pxr;

VtValue PathIntegrator::Li(const GfRay& ray, std::default_random_engine& random)
{
    std::uniform_real_distribution<float> uniform_dist(
        0.0f, 1.0f - std::numeric_limits<float>::epsilon());
    std::function<float()> uniform_float = std::bind(uniform_dist, random);

    auto color = EstimateOutGoingRadiance(ray, uniform_float, 0);

    return VtValue(GfVec3f(color[0], color[1], color[2]));
}

GfVec3f PathIntegrator::EstimateOutGoingRadiance(
    const GfRay& ray,
    const std::function<float()>& uniform_float,
    int recursion_depth)
{
    
    if (recursion_depth >= 50) {
        return {};
    }
     float RR = uniform_float();
    if (recursion_depth > 1 && RR > 0.6)
    {
        return {};
    }
    //Russian Roullete
    SurfaceInteraction si;
    if (!Intersect(ray, si)) {
        if (recursion_depth == 0) {
            return IntersectDomeLight(ray);
        }

        return GfVec3f{ 0, 0, 0 };
    }

    // This can be customized : Do we want to see the lights? (Other than dome lights?)
    if (recursion_depth == 0) {
    }

    // Flip the normal if opposite
    if (GfDot(si.shadingNormal, ray.GetDirection()) > 0) {
        si.flipNormal();
        si.PrepareTransforms();
    }
    
    GfVec3f color{ 0 };
    GfVec3f directLight = EstimateDirectLight(si, uniform_float);

    // HW7_TODO: Estimate global lighting here.
    GfVec3f globalLight = {0,0,0};
    
    GfRay Ray_next;
    auto normal = si.shadingNormal;

    //pick a random dir
    auto basis = constructONB(normal);
    GfVec3f dir_next;
    float theta = 2.0f * M_PI * uniform_float();
    float eta = uniform_float();
    float sqrteta = sqrtf(eta);
    dir_next[0] = cosf(theta) * sqrteta;
    dir_next[1] = sinf(theta) * sqrteta;
    dir_next[2] = sqrtf(1.0f - eta);
    float pdf = dir_next[2] / M_PI;
    //transform to world coordinate
    auto worldSampledDir = basis * dir_next;
    worldSampledDir = worldSampledDir.GetNormalized();
    
    //BRDF
    auto eval = si.Eval(-worldSampledDir);
    Ray_next.SetPointAndDirection(si.position + 0.000001 * normal, worldSampledDir);
    auto cosVal = GfDot(worldSampledDir, normal);
    
    globalLight = cosVal*EstimateOutGoingRadiance(Ray_next, uniform_float, recursion_depth + 1)/pdf;
    globalLight = GfCompMult(globalLight, eval);

    color = directLight + globalLight;

    return color;
}

USTC_CG_NAMESPACE_CLOSE_SCOPE
