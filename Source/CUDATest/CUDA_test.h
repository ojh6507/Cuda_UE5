// Fill out your copyright notice in the Description page of Project Settings.

#pragma once
//
//#include "CoreMinimal.h"
//#include "GameFramework/Actor.h"
//#include "ProceduralMeshComponent.h"
//#include "../../CUDALib/include/cuda_lib_test.h"
//#include "cuda_runtime.h"
//#include <algorithm>
//#include <list>
//#include "CUDA_test.generated.h"
//
//UCLASS()
//class CUDATEST_API ACUDA_test : public AActor
//{
//    GENERATED_BODY()
//
//public:
//    // Sets default values for this actor's properties
//    ACUDA_test();
//    
//protected:
//    // Called when the game starts or when spawned
//    virtual void BeginPlay() override;
//    virtual void BeginDestroy() override;
//public:
//    // Called every frame
//    virtual void Tick(double DeltaTime) override;
//    UPROPERTY(EditAnywhere)
//    FVector Frame_Length;
//    UPROPERTY(EditAnywhere)
//    double Grid_Length;
//    UPROPERTY(EditAnywhere)
//    double Particle_Size;
//
//    UPROPERTY(EditAnywhere)
//    double mass;
//    UPROPERTY(EditAnywhere)
//    int MAX_NEIGHBORS= 50;
//    
//    
//    /////////////////////////////////////////////////////////////
//    UFUNCTION(BlueprintCallable, Category = "CUDATest")
//    bool CUDA_MarchingCube() 
//    {
//     //   vertices.Empty();
//     //   Triangles.Empty();
//     //   mapping = new int[sum_e] {};
//     //   ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//     //   // cal particle position
//     //   //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//     //   double3 bound = make_double3(Frame_Length.X, Frame_Length.Y, Frame_Length.Z);
//     //   double3 base = make_double3(-Frame_Length.X, -Frame_Length.Y, -Frame_Length.Z);
//     //   std::string error_message;
//
//     //   cudaError_t cuda_status = findNeighbors(host_particlePositions, particleCount, MAX_NEIGHBORS,
//     //       neighbors.get(), numNeighbors.get(), neighborDistances.get(), &error_message);
//     //   if (cuda_status != cudaSuccess) {
//     //       UE_LOG(LogTemp, Warning, TEXT("findNeighbors failed!\n"));
//     //       UE_LOG(LogTemp, Warning, TEXT("%s"), *FString(error_message.c_str()));
//     //       return false;
//     //   }
//
//     //   cuda_status = SPH(bound ,host_particlePositions, densities, pressures, forces, particleCount, viscosities,
//     //       neighbors.get(), numNeighbors.get(), neighborDistances.get() , mass, velocities, tenssions, MAX_NEIGHBORS, &error_message);
//     //   if (cuda_status != cudaSuccess) {
//     //       UE_LOG(LogTemp, Warning, TEXT("SPH failed!\n"));
//     //       UE_LOG(LogTemp, Warning, TEXT("%s"), *FString(error_message.c_str()));
//     //       return false;
//     //   }
//
//     //   ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//     //   // marching cube
//     //   //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//     //   cuda_status = calculateInside(bound, Grid_Length, base, host_particlePositions, particleCount, inside.get(), Particle_Size, &error_message);
//     //   if (cuda_status != cudaSuccess) {
//     //       UE_LOG(LogTemp, Warning, TEXT("calculateInside failed!\n"));
//     //       UE_LOG(LogTemp, Warning, TEXT("%s"), *FString(error_message.c_str()));
//     //       return false;
//     //   }
//     //   
//     //   temp_verticesX = new double[sum_e];
//     //   temp_verticesY = new double[sum_e];
//     //   temp_verticesZ = new double[sum_e];
//     //   
//     //   cuda_status = calculateIntersections(bound, Grid_Length, intersections.get(), base, inside.get(), temp_verticesX, temp_verticesY, temp_verticesZ, &error_message);
//     //   if (cuda_status != cudaSuccess) {
//     //       UE_LOG(LogTemp, Warning, TEXT("calculateIntersections failed!\n"));
//     //       UE_LOG(LogTemp, Warning, TEXT("%s"), *FString(error_message.c_str()));
//     //       return false;
//     //   }
//     //   std::vector<double> distVerticesX(temp_verticesX, temp_verticesX + sum_e);
//     //   std::vector<double> distVerticesY(temp_verticesY, temp_verticesY + sum_e);
//     //   std::vector<double> distVerticesZ(temp_verticesZ, temp_verticesZ + sum_e);
//
// 
//     //   distVerticesX.erase(std::remove(distVerticesX.begin(), distVerticesX.end(), -9999), distVerticesX.end());
//     //   distVerticesY.erase(std::remove(distVerticesY.begin(), distVerticesY.end(), -9999), distVerticesY.end());
//     //   distVerticesZ.erase(std::remove(distVerticesZ.begin(), distVerticesZ.end(), -9999), distVerticesZ.end());
//     // 
//     //   for (int i = 0; i < distVerticesX.size(); ++i) {
//     //       mapping[i] = vertices.Num();
//     //       vertices.Add(FVector(distVerticesX[i], distVerticesY[i], distVerticesZ[i]));
//     //   }
//
//
//     //   cuda_status = calculateTriangles(bound, Grid_Length, tri_index_array, inside.get(), mapping, offset, offset_edge, &error_message);
//     //   if (cuda_status != cudaSuccess) {
//     //       UE_LOG(LogTemp, Warning, TEXT("calculateTriangles failed!\n"));
//     //       UE_LOG(LogTemp, Warning, TEXT("%s"), *FString(error_message.c_str()));
//     //       return false;
//     //   }
//
//     //   std::list<int> dest(tri_index_array, tri_index_array + max_size);
//     //   dest.remove(-9999);
//     //   for (const int idx : dest) {
//     //       Triangles.Add(idx);
//     //   }
//     ///*    UE_LOG(LogTemp, Warning, TEXT("====================================="));
//     //   if (!printDebug) {
//     //       for (int i{}; i < vertices.Num(); ++i) {
//     //           UE_LOG(LogTemp, Warning, TEXT("%d : %f, %f, %f"),i, vertices[i].X, vertices[i].Y, vertices[i].Z);
//     //       }
//     //     
//     //   }
//     //   UE_LOG(LogTemp, Warning, TEXT("====================================="));*/
//     //   delete[] mapping;
//     //  
//        return true;
//    }
//    
//private:
//    void CreateProceduralMesh();
//private:
//    
//    TArray<FVector> ParticlePositions;
//    TArray<FVector> vertices;
//    TArray<int> Triangles;
//
//    double* temp_verticesX;
//    double* temp_verticesY;
//    double* temp_verticesZ;
//
//    int* tri_index_array;
//    int* mapping;
//
//    std::unique_ptr<int[]> neighbors;
//    std::unique_ptr<int[]> numNeighbors;
//    std::unique_ptr<double[]> neighborDistances;
//    int particleCount;
//
//
//    bool printDebug{};
//   std::unique_ptr<int[]> inside;
//   std::unique_ptr<double3[]>  intersections;
//    UProceduralMeshComponent* ProceduralMesh;
//    double3* host_particlePositions;
//
//    bool init{};
//    int dx{};
//    int dy{};
//    int dz{};
//
//    int3 total_edge;
//    int offset[8]{};
//    int offset_edge[12]{};
//
//    double* densities;
//    double* pressures;
//    int max_size{};
//    double3* forces;
//    double3* viscosities;
//    double3* velocities;
//    double3* tenssions;
//    double time{};
//
//    int sum_v{};
//    int sum_e{};
//};
