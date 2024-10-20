// Fill out your copyright notice in the Description page of Project Settings.

#pragma once

#include "CoreMinimal.h"
#include "GameFramework/Actor.h"
#include "ProceduralMeshComponent.h"
#include "Components/InstancedStaticMeshComponent.h"
#include "CudaTest1.h"
#include "MyActor.generated.h"

UCLASS()
class CUDATEST_API AMyActor : public AActor
{
	GENERATED_BODY()
	
public:	
	// Sets default values for this actor's properties
	AMyActor();
public:
	UPROPERTY(EditAnyWhere)
	FRotator Rotation = FRotator(0, 0, 0);
	UPROPERTY(EditAnyWhere)
	FVector Scale = FVector(0.2f, 0.2f, 0.2f);
	UPROPERTY(VisibleAnywhere)
	UInstancedStaticMeshComponent* Spheres;

protected:
	virtual void BeginPlay() override;
	virtual void Tick(float DeltaTime) override;
private:	
	void initParticles();
	void generateParticleSphere();
	FVector Transdouble3ToFvector(double3 a);
	void updateParticles(double deltaTime);
	void updateParticlesSphere();
	
private:
	bool init{};
	int particleCount{1000};
	int particleCubeWidth = 12;
	Particle* particles;
	UProceduralMeshComponent* ProceduralMesh;
	SPHSettings sphSettings;
	TArray <FTransform> sphereModelTransform;
	double3* particlesPosition;
};
