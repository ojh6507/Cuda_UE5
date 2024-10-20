// Fill out your copyright notice in the Description page of Project Settings.

#pragma once

#include "CoreMinimal.h"
#include "GameFramework/Actor.h"
#include "Components/BoxComponent.h"
#include "Components/InstancedStaticMeshComponent.h"
#include "ProceduralMeshComponent.h"
#include "Particle.generated.h"

struct particle_1 {
	FVector Velocity{};
	FVector Position{};
	float density{};
	float properties{};
	bool bound{};
	float dst{};
};


struct FEntry
{
	int32 Index;
	uint32 CellKey;

	FEntry(int32 InIndex = 0, uint32 InCellKey = 0)
		: Index(InIndex), CellKey(InCellKey) {}
};

UCLASS()
class CUDATEST_API AParticle : public AActor
{
	GENERATED_BODY()
	
public:	
	// Sets default values for this actor's properties
	AParticle();

protected:
	// Called when the game starts or when spawned
	virtual void BeginPlay() override;
	virtual void Tick(float DeltaTime) override;
public:	
	UProceduralMeshComponent* ProcMesh;
	UPROPERTY(EditAnyWhere)
	float particleSpacing;
	UPROPERTY(EditAnyWhere)
	int numParticles = 10;
	
	UPROPERTY(EditAnyWhere)
	FRotator Rotation = FRotator(0, 0, 0);
	UPROPERTY(EditAnyWhere)
	FVector Scale = FVector(0.2f, 0.2f, 0.2f);
	UPROPERTY(EditAnyWhere)
	float Gravity = 10;
	UPROPERTY(EditAnyWhere)
	float smoothingRadius = 200;
	UPROPERTY(EditAnyWhere)
	float CollisionDapping;
	UPROPERTY(EditAnyWhere)
	float targetDensity = 10;
	UPROPERTY(EditAnyWhere)
	float pressureMultiplier = 10;
	float PARTICLE_RADIUS;
	UPROPERTY(EditAnyWhere)
	UMaterial* MyMaterial;
	UPROPERTY(EditAnyWhere)
	UBoxComponent* BoundaryBox;
	UPROPERTY(VisibleAnywhere)
	UInstancedStaticMeshComponent* Spheres;

private:
	void CreateParticles(int seed);
	void ResolveCollision(FVector& Position, FVector& velocity, bool& bound);
	float SmoothingKernel(float radius, float dst);
	float SmoothingKernelDerivative(float dst, float rad);

	void UpdateDensities();
	float CalculateDensity(FVector samplePoint);

	float CalculateProperty(FVector samplePoint);
	float ConvertDensityToPressure(float density);
	float  CalculateSharedPressure(float densityA, float densityB);
	FVector CalculatePropertyGradient(FVector samplePoint);
	FVector CalculatePressureForce(int particleIndex);
	float ExampleFunction(FVector Position);

	void SimulationStep(float deltaTime);
	void UpdateSpatialLookup(const TArray< particle_1>& Points);
	void PositionToCellCoord(const FVector& Position, int32& OutCellX, int32& OutCellY, int32& OutCellZ) const;
	uint32 HashCell(int32 CellX, int32 CellY, int32 CellZ) const;
	uint32 GetKeyFromHash(uint32 Hash) const;
	void ForeachPointWithinRadius(const FVector& SamplePoint);
	FVector CalculateViscosityForce(int particleIndex);
	float ViscositySmoothingKernel(float distance);
private:
	void MarchingCubeInit(FVector bound, float length);
	void Calculate();
	FVector checkClosetParticle(FVector v);
	FVector countInter(const FVector& a, const FVector& b);
private:
	TArray<FVector> vertices;
	TArray<int> triangleIdx;
	TArray<FVector> normals;
	FVector bound;
	float l;
	int offset_edge[12];
	int m_offset[8];
	TArray<FVector> intersections;
	TArray<FVector> inside;
	TArray<int> mapping;
	FVector base;
	FVector total_edge;
	int dx;
	int dy;
	int dz;

	int sum_e;
	int sum_v;

private:
	TArray <FIntVector>CellOffsets;
	
	TArray< particle_1> particles;
	TArray< particle_1> predict_particles;
	float particleSize;
	const float mass = 1;
	float min_z;
	TArray < FVector> GravityVector;
	// 공간 해시를 위한 배열
	TArray<FEntry> SpatialLookup;
	TArray<int32> StartIndices;
	bool init{};
};
