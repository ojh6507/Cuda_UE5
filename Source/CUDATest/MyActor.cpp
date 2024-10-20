// Fill out your copyright notice in the Description page of Project Settings.


#include "MyActor.h"
#include <fstream>
// Sets default values
AMyActor::AMyActor()
{
 	// Set this actor to call Tick() every frame.  You can turn this off to improve performance if you don't need it.
	PrimaryActorTick.bCanEverTick = true;
	ProceduralMesh = CreateDefaultSubobject<UProceduralMeshComponent>(TEXT("CustomMesh"));
	RootComponent = ProceduralMesh;

	Spheres = CreateDefaultSubobject<UInstancedStaticMeshComponent>(TEXT("Spheres"));
	Spheres->SetupAttachment(RootComponent);

	static ConstructorHelpers::FObjectFinder<UStaticMesh> SphereMesh(TEXT("/Game/StarterContent/Shapes/Shape_Sphere.Shape_Sphere"));
	if (SphereMesh.Succeeded()) {
		Spheres->SetStaticMesh(SphereMesh.Object);
	}

}

// Called when the game starts or when spawned
void AMyActor::BeginPlay()
{
	Super::BeginPlay();
	//SPHSettings(double mass, double restDensity, double gasConst, double viscosity,double h, double g, double tension)
	sphSettings = SPHSettings(0.02f, 1000, 100, 1.04f, 0.8f, -0.8f, 0.5f);
	particleCount = particleCubeWidth * particleCubeWidth * particleCubeWidth;
	initParticles();
	generateParticleSphere();
	//for (int i{}; i < particleCount; ++i) {
	//	FVector pos = Transdouble3ToFvector(particles[i].Position);
	//		//UE_LOG(LogTemp, Warning, TEXT("index: %d ] %f, %f, %f"), i, pos.X, pos.Y, pos.Z);
	//}
	init = true;
}

// Called every frame
void AMyActor::Tick(float DeltaTime)
{
	Super::Tick(DeltaTime);
	if (!init) return;

	updateParticles(DeltaTime / 2);
	updateParticlesSphere();
}

void AMyActor::initParticles()
{
	std::srand(1024);
	double particleSeperation = sphSettings.h + 0.01f;
	
	particles = (Particle*)malloc(sizeof(Particle) * particleCount);
	if (!particles) return;
	particlesPosition = (double3*)malloc(sizeof(double3) * particleCount);
	memset(particlesPosition, 0, sizeof(double3) * particleCount);

	sphereModelTransform.Empty();
	sphereModelTransform.Reserve(particleCubeWidth * particleCubeWidth * particleCubeWidth);
	for (int i{}; i< particleCubeWidth; ++i) {
		for (int j{}; j < particleCubeWidth; ++j) {
			for (int k{}; k < particleCubeWidth; ++k) {
				double ranX
					= (double(rand()) / double((RAND_MAX)) * 0.5f - 1)
					* sphSettings.h / 10;
				double ranY
					= (double(rand()) / double((RAND_MAX)) * 0.5f - 1)
					* sphSettings.h / 10;
				double ranZ
					= (double(rand()) / double((RAND_MAX)) * 0.5f - 1)
					* sphSettings.h / 10;
				double3 nParticlePos = make_double3(
					i * particleSeperation + ranX - 1.5f,
					j * particleSeperation + ranY + sphSettings.h + 4.1f,
					k * particleSeperation + ranZ - 1.5f);


				size_t particleIndex = i + (j + particleCubeWidth * k) * particleCubeWidth;
				Particle* tparticle = &particles[particleIndex];
				tparticle->Position = nParticlePos;
				tparticle->velocity = make_double3(0,0,0);
				
				// Transform 추가 sphere
				FVector ParticlePosition = Transdouble3ToFvector(tparticle->Position);
				FVector SphereScale(sphSettings.sphereScale, sphSettings.sphereScale, sphSettings.sphereScale);
				FTransform ParticleTransform;
				ParticleTransform.SetLocation(ParticlePosition * 20);
				ParticleTransform.SetScale3D(SphereScale/10);
				sphereModelTransform.Add(ParticleTransform);
			}
		}
	}

	//std::ifstream in{ "C:/Users/ojh65/OneDrive/Documents/JELLY.obj" };
	//char ch;
	//float v1;
	//float v2;
	//float v3;

	//if (!in)
	//	UE_LOG(LogTemp, Warning, TEXT("Error!"));
	//TArray<double3> positions;
	//while (in) {
	//	in >> ch >> v1 >> v2 >> v3;
	//	if (ch == 'v') {
	//		positions.Add(make_double3(v1 * 0.01, v2 * 0.01, v3 * 0.01));
	//		positions.Add(make_double3(v1 * 0.01, v2 * 0.01, v3 * 0.01));
	//		positions.Add(make_double3(v1 * 0.01, v2 * 0.01, v3 * 0.01));
	//	}
	//}
	//particleCubeWidth = positions.Num();
	//particleCount = positions.Num();;
	

	//for (int i = 0; i < particleCubeWidth; i += 1)
	//{
	//	particles[i].Position = positions[i];
	//	particles[i].velocity = make_double3(0, 0, 0);

	//	// Transform 추가 sphere
	//	FVector ParticlePosition = Transdouble3ToFvector(particles[i].Position);
	//	FVector SphereScale(sphSettings.sphereScale, sphSettings.sphereScale, sphSettings.sphereScale);
	//	FTransform ParticleTransform;
	//	ParticleTransform.SetLocation(ParticlePosition * 20);
	//	ParticleTransform.SetScale3D(SphereScale/8);
	//	sphereModelTransform.Add(ParticleTransform);
	//}
}

void AMyActor::generateParticleSphere()
{
	for (int i{}; i < particleCount; ++i) {
		Spheres->AddInstance(sphereModelTransform[i]);
		//auto t = sphereModelTransform[i].GetLocation();
		//UE_LOG(LogTemp, Warning, TEXT("%f, %f, %f"), t.X, t.Y, t.Z);
	}
}

FVector AMyActor::Transdouble3ToFvector(double3 a)
{
	return FVector(a.x, a.z, a.y);
}

void AMyActor::updateParticles(double deltaTime)
{
	updateParticlesGPU(particles, particlesPosition, particleCount, sphSettings, deltaTime);
	/*for (int i{}; i < particleCount; ++i) {
		FVector pos = Transdouble3ToFvector(particles[i].Position);
		if (pos.X > 10 || pos.Y > 10 || pos.Z > 10 || pos.X < -10 || pos.Y < -10 || pos.Z < 0) {
			UE_LOG(LogTemp, Warning, TEXT("index: %d ] %f, %f, %f"), i, pos.X, pos.Y, pos.Z);
			continue;
		}
	}*/
}

void AMyActor::updateParticlesSphere()
{
	FTransform CurrentTransform;
	for (int i{}; i < particleCount; ++i) {
		Spheres->GetInstanceTransform(i, CurrentTransform);
		FVector pos = Transdouble3ToFvector(particlesPosition[i]);
		if (pos.X > 10 || pos.Y > 10 || pos.Z > 10 || pos.X < -10 || pos.Y < -10 || pos.Z < 0) {
			//UE_LOG(LogTemp, Warning, TEXT("index: %d ] %f, %f, %f"), i, pos.X, pos.Y, pos.Z);
			continue;
		}
		CurrentTransform.SetLocation(pos * 20);
		Spheres->UpdateInstanceTransform(i,CurrentTransform,true, true);
	}

	//UE_LOG(LogTemp, Warning, TEXT("START START START START START START START START START START START START "));
	//for (int i{}; i < particleCount; ++i) {
	//	FVector pos = Transdouble3ToFvector(particles[i].force);
	//	float d = particles[i].density;
	//	UE_LOG(LogTemp, Warning, TEXT("index: %d ] %f, %f, %f"), i, pos.X/d, pos.Y / d, pos.Z / d);
	//}
	//UE_LOG(LogTemp, Warning, TEXT("END END END END END END END END END END END END END END END END END END "));
	//


}
