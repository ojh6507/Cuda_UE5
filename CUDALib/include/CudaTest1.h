#include "cuda_runtime.h"
#include "SphCu.cuh"

void updateParticlesGPU(
    Particle* particles, double3* particleTransforms,
    const size_t particleCount, SPHSettings& settings, float deltaTime);