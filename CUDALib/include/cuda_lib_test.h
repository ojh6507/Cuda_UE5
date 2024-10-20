#pragma once
 
#include <string>
#include <vector>
#include "cuda_runtime.h"
#include "vector_types.h"
#include "vector_functions.h"
#include "device_launch_parameters.h"
#include "lookup_list.h"

cudaError_t calculateInside(const double3 _bound, double gridsize, const double3 base, const double3* particlePositions,
    int numParticles, int* inside,double Particle_Size, std::string* error_message);

cudaError_t calculateIntersections(const double3 _bound, double gridsize, double3* intersections, const double3 base, const int* inside,
    double* verticesX, double* verticesY, double* verticesZ, std::string* error_message);

cudaError_t  calculateTriangles(const double3 _bound, double gridsize, int* tri_index_array, const int* inside, const int* mapping,
    const int* offset, const int* offset_edge, std::string* error_message);

cudaError_t findNeighbors(const double3* ParticlePositions, int numParticles, int MAX_NEIGHBORS,
    int* neighbors, int* numNeighbors, float* neighborDistances, std::string* error_message);

cudaError_t SPH(const double3 bound ,double3* particlePositions, float* densities, float* pressures, double3* forces, int numParticles, double3* viscosities,
    const int* neighbors, const int* numNeighbors, const float* neighborDistances, float mass,
    double3* velocities, double3* tenssions, int MAX_NEIGHBORS, std::string* error_message);