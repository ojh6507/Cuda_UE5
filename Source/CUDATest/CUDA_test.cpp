// Fill out your copyright notice in the Description page of Project Settings.


#include "CUDA_test.h"

//
//// Sets default values
//ACUDA_test::ACUDA_test()
//{
// 	// Set this actor to call Tick() every frame.  You can turn this off to improve performance if you don't need it.
//	PrimaryActorTick.bCanEverTick = true;
//    ProceduralMesh = CreateDefaultSubobject<UProceduralMeshComponent>(TEXT("CustomMesh"));
//    RootComponent = ProceduralMesh;
//}
//
//// Called when the game starts or when spawned
//void ACUDA_test::BeginPlay()
//{
//	Super::BeginPlay();
//    sum_v = (1 + (int)(Frame_Length.X * 2 / Grid_Length)) *
//        (1 + (int)(Frame_Length.Y * 2 / Grid_Length)) * (1 + (int)(Frame_Length.Z * 2 / Grid_Length));
//    sum_e = (1 + (int)(Frame_Length.X * 2 / Grid_Length)) * (1 + (int)(Frame_Length.Y * 2 / Grid_Length)) * (1 + (int)(Frame_Length.Z * 2 / Grid_Length)) * 3;
//    
//    total_edge = make_int3((int)(Frame_Length.X * 2 / Grid_Length), (int)(Frame_Length.Y * 2 / Grid_Length), (int)(Frame_Length.Z * 2 / Grid_Length));
//    max_size = total_edge.x * total_edge.y * total_edge.z * 16;
//    
//    for (int i = 0; i < 12; i += 1)
//        for (int j = 0; j < 12; j += 1)
//            for (int k = 0; k < 12; k += 1)
//                ParticlePositions.Add(FVector(-0.50 + i * 0.03, -0.2 + j * 0.03, -0.05 + k * 0.03));
//
//    particleCount = ParticlePositions.Num();
//    neighbors = std::make_unique<int[]>(MAX_NEIGHBORS * particleCount);
//    numNeighbors = std::make_unique<int[]>(particleCount);
//    neighborDistances = std::make_unique<double[]>(MAX_NEIGHBORS * particleCount);
//    inside = std::make_unique<int[]>(sum_v);
//    intersections = std::make_unique<double3[]>(sum_e);
//    densities = new double[particleCount];
//    pressures = new double[particleCount];
//    forces = new double3[particleCount];
//    viscosities = new double3[particleCount];
//    velocities = new double3[particleCount];
//    tenssions = new double3[particleCount];  
//    host_particlePositions = new double3[particleCount];
//    
//
//   temp_verticesX = new double[sum_e];
//   temp_verticesY = new double[sum_e];
//   temp_verticesZ = new double[sum_e];
//
//    mapping = new int[sum_e];
//    tri_index_array = new int[max_size];
//
//    dx = (total_edge.z + 1) * (total_edge.y + 1);
//    dy = total_edge.z + 1;
//    dz = 1;
//
//
//    offset[0] = 0;			  offset[1] = dx;
//    offset[2] = dx + dz;	  offset[3] = dz;
//    offset[4] = dy;			  offset[5] = dy + dx;
//    offset[6] = dx + dy + dz; offset[7] = dy + dz;
//
//    offset_edge[0] = 0;
//    offset_edge[1] = dx * 3 + 2;
//    offset_edge[2] = dz * 3;
//    offset_edge[3] = 2;
//    offset_edge[4] = dy * 3;
//    offset_edge[5] = (dy + dx) * 3 + 2;
//    offset_edge[6] = (dy + dz) * 3;
//    offset_edge[7] = dy * 3 + 2;
//    offset_edge[8] = 1;
//    offset_edge[9] = dx * 3 + 1;
//    offset_edge[10] = (dx + dz) * 3 + 1;
//    offset_edge[11] = dz * 3 + 1;
//
//    for (int i = 0; i < particleCount; ++i) {
//        host_particlePositions[i] = make_double3(ParticlePositions[i].X, ParticlePositions[i].Y, ParticlePositions[i].Z);
//    }
//
//    for (int i{}; i < particleCount; ++i)
//        velocities[i] = make_double3(0, 0, 0);
//    init = true;
// 
//}
//
//
//// Called every frame
//void ACUDA_test::Tick(double DeltaTime)
//{
//    Super::Tick(DeltaTime);
//    if (init) {
//        bool result = CUDA_MarchingCube();
//        if (!result) {
//            UE_LOG(LogTemp, Warning, TEXT("ERROR"));
//        }
//        else {
//            CreateProceduralMesh();
//           
//        }
//    }
//}
//
//void ACUDA_test::CreateProceduralMesh()
//{
//    /*ProceduralMesh->CreateMeshSection_LinearColor(0, vertices, Triangles, TArray<FVector>(), TArray<FVector2D>(), TArray<FLinearColor>(), TArray<FProcMeshTangent>(), true);
//    ProceduralMesh->UpdateBounds();*/
//}
//
//
//void ACUDA_test::BeginDestroy()
//{
//    Super::BeginDestroy();
//    
//    delete[] host_particlePositions;
//    
//    delete[] densities;
//    delete[] pressures;
//
//    delete[] forces;
//    delete[] viscosities;
//    delete[] velocities;
//    delete[] tenssions;
//    
//
// 
//    delete[] temp_verticesX;
//    delete[] temp_verticesY;
//    delete[] temp_verticesZ;
//
//    delete[] tri_index_array;
//}