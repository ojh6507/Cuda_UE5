#include "Particle.h"
#include <fstream>
#include "lookup_list.h"
#include "Algo/MinElement.h"

AParticle::AParticle()
{
 	// Set this actor to call Tick() every frame.  You can turn this off to improve performance if you don't need it.
	PrimaryActorTick.bCanEverTick = true;
	BoundaryBox = CreateDefaultSubobject<UBoxComponent>(TEXT("BoundaryBox"));
	RootComponent = BoundaryBox;

	BoundaryBox->InitBoxExtent(FVector(100, 100, 100));

	ProcMesh = CreateDefaultSubobject<UProceduralMeshComponent>(TEXT("CustomMesh"));
	ProcMesh->SetupAttachment(RootComponent);
	Spheres = CreateDefaultSubobject<UInstancedStaticMeshComponent>(TEXT("Spheres"));
	Spheres->SetupAttachment(RootComponent);

	static ConstructorHelpers::FObjectFinder<UStaticMesh> SphereMesh(TEXT("/Game/StarterContent/Shapes/Shape_Sphere.Shape_Sphere"));
	if (SphereMesh.Succeeded()) {
		Spheres->SetStaticMesh(SphereMesh.Object);
	}
	else {
		UE_LOG(LogTemp, Warning, TEXT("Error Can't find"));
	}
	
}
void AParticle::BeginPlay()
{
	Super::BeginPlay();
	
	particleSize = Scale.X;
	//CellOffsets = {
	//	FIntVector(-1, -1, -1), FIntVector(0, -1, -1), FIntVector(1, -1, -1),
	//	FIntVector(-1, 0, -1), FIntVector(0, 0, -1), FIntVector(1, 0, -1),
	//	FIntVector(-1, 1, -1), FIntVector(0, 1, -1), FIntVector(1, 1, -1),

	//	FIntVector(-1, -1, 0), FIntVector(0, -1, 0), FIntVector(1, -1, 0),
	//	FIntVector(-1, 0, 0), /* FIntVector(0, 0, 0), */ FIntVector(1, 0, 0),
	//	FIntVector(-1, 1, 0), FIntVector(0, 1, 0), FIntVector(1, 1, 0),

	//	FIntVector(-1, -1, 1), FIntVector(0, -1, 1), FIntVector(1, -1, 1),
	//	FIntVector(-1, 0, 1), FIntVector(0, 0, 1), FIntVector(1, 0, 1),
	//	FIntVector(-1, 1, 1), FIntVector(0, 1, 1), FIntVector(1, 1, 1)
	//};


	bound = FVector(0.8,0.8,0.8);
	l = 0.08;
	total_edge = FVector((int)(bound.X * 2 / l), (int)(bound.Y * 2 / l), (int)(bound.Z * 2 / l));
	dx = (total_edge.Z + 1) * (total_edge.Y + 1);
	dy = total_edge.Z + 1;
	dz = 1;
	sum_v = (1 + (int)(bound.X * 2 / l)) * (1 + (int)(bound.Y * 2 / l)) * (1 + (int)(bound.Z * 2 / l));
	sum_e = (1 + (int)(bound.X * 2 / l)) * (1 + (int)(bound.Y * 2 / l)) * (1 + (int)(bound.Z * 2 / l)) * 3;

	base = -bound;


	m_offset[0] = 0;			m_offset[1] = dx;
	m_offset[2] = dx + dz;	m_offset[3] = dz;
	m_offset[4] = dy;			m_offset[5] = dy + dx;
	m_offset[6] = dx + dy + dz; m_offset[7] = dy + dz;

	offset_edge[0] = 0;
	offset_edge[1] = dx * 3 + 2;
	offset_edge[2] = dz * 3;
	offset_edge[3] = 2;
	offset_edge[4] = dy * 3;
	offset_edge[5] = (dy + dx) * 3 + 2;
	offset_edge[6] = (dy + dz) * 3;
	offset_edge[7] = dy * 3 + 2;
	offset_edge[8] = 1;
	offset_edge[9] = dx * 3 + 1;
	offset_edge[10] = (dx + dz) * 3 + 1;
	offset_edge[11] = dz * 3 + 1;

	for (int x = -1; x <= 1; ++x) {
		for (int y = -1; y <= 1; ++y) {
			for (int z = -1; z <= 1; ++z) {
				if (x == 0 && y == 0 && z == 0) continue;
				CellOffsets.Add(FIntVector(x * dx, y * dy, z * dz));
			}
		}
	}

	
	CreateParticles(10);
	init = true;
	
}
void AParticle::Tick(float DeltaTime)
{
	Super::Tick(DeltaTime);
	if (init) {


		Calculate();
		if (vertices.Num() && triangleIdx.Num()) {
			normals.SetNumZeroed(vertices.Num());
			for (int32 i = 0; i < triangleIdx.Num(); i += 3) {
				// 삼각형을 구성하는 세 정점
				FVector vertex1 = vertices[triangleIdx[i]];
				FVector vertex2 = vertices[triangleIdx[i + 1]];
				FVector vertex3 = vertices[triangleIdx[i + 2]];

				// 두 벡터를 구하고 외적을 사용하여 면의 법선을 계산
				FVector edge1 = vertex2 - vertex1;
				FVector edge2 = vertex3 - vertex1;
				FVector triangleNormal = FVector::CrossProduct(edge1, edge2).GetSafeNormal();

				// 각 정점에 대한 법선을 누적
				normals[triangleIdx[i]] += triangleNormal;
				normals[triangleIdx[i + 1]] += triangleNormal;
				normals[triangleIdx[i + 2]] += triangleNormal;
			}

			// 각 정점의 법선을 정규화
			for (FVector& normal : normals) {
				normal.Normalize();
			}
			TArray<FVector2D>uv{};
			TArray<FColor> VertexColors{};
			TArray<FProcMeshTangent> Tangents{};

			ProcMesh->CreateMeshSection(0, vertices, triangleIdx, normals, uv, VertexColors, Tangents, true);
			ProcMesh->SetMaterial(0, MyMaterial);
		}
		SimulationStep(DeltaTime);
		
		for (int i{}; i < numParticles; ++i) {
			//FTransform CurrentTransform;

			//Spheres->GetInstanceTransform(i, CurrentTransform, true);
			/*CurrentTransform.SetLocation(particles[i].Position + GetActorLocation());
			Spheres->UpdateInstanceTransform(i, CurrentTransform, true, true);*/
			particles[i].Position += particles[i].Velocity * DeltaTime;
			ResolveCollision(particles[i].Position, particles[i].Velocity, particles[i].bound );
		}
		//Spheres->MarkRenderStateDirty();
	}
}
void AParticle::CreateParticles(int seed) 
{
	FRandomStream RandomStream;
	RandomStream.Initialize(seed);
	int particlesPerRow = (int)sqrt(numParticles);
	int particlesPerCol = (numParticles - 1) / particlesPerRow + 1;
	float spacing = particleSize * 2 + particleSpacing;

	//std::ifstream in{ "C:/Users/ojh65/OneDrive/Documents/JELLY.obj" };
	//char ch;
	//float v1;
	//float v2;
	//float v3;

	//if(!in) 
	//	UE_LOG(LogTemp, Warning, TEXT("Error!"));
	//while (in) {
	//	in >> ch >> v1 >> v2 >> v3;
	//	if (ch == 'v') {
	//		particles.Add({});
	//		particles[particles.Num() -1].Position = FVector(v1 * 0.0045, v3 * 0.0045, v2 * 0.0035);
	//		
	//	}
	//}

    for (int i = 0; i < 12; i += 1)
        for (int j = 0; j < 12; j += 1)
            for (int k = 0; k < 12; k += 1){
				particles.Add({});
                particles[i].Position = (FVector(-0.50 + i * 0.03, -0.2 + j * 0.03, -0.05 + k * 0.03));
			}
	if (particles.Num() > 0) {
		const  particle_1* lowestZParticle = Algo::MinElement(particles, [](const  particle_1& A, const  particle_1& B) {
			return A.Position.Z < B.Position.Z;
			});
		min_z = lowestZParticle->Position.Z;
	}
	numParticles = particles.Num();
	predict_particles.Init({}, numParticles);
	SpatialLookup.Init({}, numParticles);
	StartIndices.Init({}, numParticles);
	GravityVector.Init({}, numParticles);

	FTransform SphereTransform;
	for (int i{}; i < particles.Num(); ++i) {
		//float x = (RandomStream.FRand() - 0.5) * BoundaryBox->GetScaledBoxExtent().X;
		//float z = (RandomStream.FRand() - 0.5) * BoundaryBox->GetScaledBoxExtent().Z;
		
		/*float x = (i % particlesPerRow - particlesPerRow / 2.f + 0.5f) * spacing;
		float y = (i % particlesPerRow - particlesPerRow / 2.f + 0.5f) * spacing;
		float z = (i / particlesPerRow - particlesPerCol / 2.f + 0.5f) * spacing;*/
		//particles[i].Position = FVector(x, 0, z);
		particles[i].properties = ExampleFunction(particles[i].Position);

		//SphereTransform.SetLocation(particles[i].Position);
		//SphereTransform.SetRotation(Rotation.Quaternion());
		//SphereTransform.SetScale3D(Scale);
		//Spheres->AddInstance(SphereTransform);
	}

}
void AParticle::ResolveCollision(FVector& Position, FVector& velocity, bool& collide)
{
	FVector halfBoxSize = bound * 2;

	if (abs(Position.X) > halfBoxSize.X -1e-6) {
		Position.X = bound.X * (Position.X / abs(Position.X)) * 2 - Position.X;
		velocity.X *= -1 * CollisionDapping;
		collide = true;
	}
	if (abs(Position.Y) > halfBoxSize.Y - 1e-6) {
		Position.Y = bound.Y * (Position.Y / abs(Position.Y)) * 2 - Position.Y;
		velocity.Y *= -1 * CollisionDapping;
		collide = true;
	}
	/*if (Position.Z < min_z - 1e-6) {
		Position.Z = min_z * FMath::Sign(Position.Z);
		velocity.Z *= -1 * CollisionDapping;
	}*/
	if (abs(Position.Z) > halfBoxSize.Z - 1e-6) {
		Position.Z= bound.Z * (Position.Z / abs(Position.Z)) * 2 - Position.Z;
		velocity.Z *= -1 * CollisionDapping;
	}
}
float AParticle::SmoothingKernel(float radius, float dst) 
{
	if (dst >= radius) return 0;

	float volume = PI * FMath::Pow(radius, 4) / 6;
	return (radius - dst) * (radius - dst) / volume;
}
float AParticle::SmoothingKernelDerivative(float dst, float rad)
{
	if (dst >= rad) return 0;

	float scale = 12 / (PI * FMath::Pow(rad, 4));
	return (dst - rad) * scale;
}
void AParticle::UpdateDensities()
{
	ParallelFor(numParticles,
		[&](int32 i) {
			particles[i].density = CalculateDensity(particles[i].Position);
		}
	);
}

float AParticle::CalculateDensity(FVector samplePoint)
{
	float density = 0;
	
	for (const auto& p : particles) {
		float dst = (p.Position - samplePoint).Length();
		float influence = SmoothingKernel(smoothingRadius, dst);
		density += mass * influence;
	}
	return density;
}
float AParticle::CalculateProperty(FVector samplePoint)
{
	float property{};
	for (int i{}; i < particles.Num(); ++i) {
		float dst = (particles[i].Position - samplePoint).Length();
		float influence = SmoothingKernel(dst, smoothingRadius);
		float density = particles[i].density;
		if(density * influence)
			property += particles[i].properties * mass / density * influence;
	}
	return property;
}
float AParticle::ExampleFunction(FVector position)
{
	return FMath::Cos(position.Y - 3 + FMath::Sin(position.X));
}

FVector AParticle::CalculatePressureForce(int particleIndex)
{
	FVector pressureForce = FVector::Zero();
	for (int i{}; i < particles.Num(); ++i) {
		if (particleIndex == i) continue;

		FVector offset = particles[i].Position - particles[particleIndex].Position;
		float dst = offset.Length();
		FVector dir = dst == 0? FMath::VRand() : offset / dst;
		
		float slope = SmoothingKernelDerivative(dst, smoothingRadius);
		float density = particles[i].density;
		float sharedPressure = CalculateSharedPressure(density, particles[particleIndex].density);
		pressureForce += sharedPressure * dir * slope * mass / density;
	}
	return pressureForce;
}
float  AParticle::CalculateSharedPressure(float densityA, float densityB) {
	float pressureA = ConvertDensityToPressure(densityA);
	float pressureB = ConvertDensityToPressure(densityB);
	return (pressureA + pressureB) / 2;
}
void AParticle::SimulationStep(float deltaTime)
{
	ParallelFor(numParticles,
		[&](int32 i)
		{
				predict_particles[i].bound = particles[i].bound;
			if (!predict_particles[i].bound) {
				particles[i].Velocity += FVector::DownVector * Gravity * deltaTime;
				predict_particles[i].Velocity = particles[i].Velocity;
				predict_particles[i].Position = particles[i].Position + particles[i].Velocity * 1 / 120.f;
			}
			else {
				GravityVector[i] += FVector::DownVector * Gravity * deltaTime;
				predict_particles[i].Position = particles[i].Position + GravityVector[i] * 1 / 120.f;
			}
		});
	UpdateSpatialLookup(predict_particles);
	ParallelFor(numParticles,
		[&](int32 i)
		{
			particles[i].density = CalculateDensity(predict_particles[i].Position);
		}
	);
	ParallelFor(numParticles,
		[&](int32 i) {
			if (!particles[i].bound) {
				FVector pressureForce = CalculatePressureForce(i);
				FVector pressureAcceleration = pressureForce / particles[i].density;
				particles[i].Velocity = pressureAcceleration * deltaTime;
			}
		});
}
FVector AParticle::CalculateViscosityForce(int particleIndex) 
{
	/*FVector viscosityForce = FVector::ZeroVector;
	FVector position = particles[particleIndex].Position;
	for ( int otherIndex : GetNeighbours(position)) {
		float dst = (position - particles[otherIndex].Position).Length();
		float influence = ViscositySmoothingKernel(dst);
	}*/

	return FVector();
}
float AParticle::ViscositySmoothingKernel(float distance)
{

	return 0.f;
}
void AParticle::MarchingCubeInit(FVector _bound, float length)
{

	
}
void AParticle::Calculate()
{
	FVector nonVector = FVector(-9999, -9999, -9999);
	ProcMesh->ClearAllMeshSections();
	vertices.Empty();
	triangleIdx.Empty();
	mapping.Init(0, sum_e);
	inside.Init(nonVector, sum_v);

	intersections.Init(nonVector, sum_e);

	// count inside
	int p = 0;
	for (int i = 0; i <= total_edge[0]; ++i) {
		for (int j = 0; j <= total_edge[1]; ++j) {
			for (int k = 0; k <= total_edge[2]; ++k, ++p) {
				FVector vetex(base[0] + i * l, base[1] + j * l, base[2] + k * l);
				inside[p] = (checkClosetParticle(vetex));
				
			}
		}
	}

	// cout intersection
	p = 0; 
	 particle_1* tmp = nullptr;
	for (int i = 0; i <= total_edge[0]; ++i) {
		for (int j = 0; j <= total_edge[1]; ++j) {
			for (int k = 0; k <= total_edge[2]; ++k, p += 3) {
				FVector v0(base[0] + i * l, base[1] + j * l, base[2] + k * l);
				FVector v1(base[0] + i * l + l, base[1] + j * l, base[2] + k * l);
				FVector v2(base[0] + i * l, base[1] + j * l + l, base[2] + k * l);
				FVector v3(base[0] + i * l, base[1] + j * l, base[2] + k * l + l);

				if ((p / 3 + dx) < inside.Num()) {
					if (inside[p / 3] != nonVector || inside[p / 3 + dx] != nonVector) {
						intersections[p] = countInter(v0, v1);
					}
				}
				if ((p / 3 + dy) < inside.Num()) {
					if (inside[p / 3] != nonVector || inside[p / 3 + dy] != nonVector) {
						intersections[p + 1] = countInter(v0, v2);
					}
				}
				if ((p / 3 + dz) < inside.Num()) {
					if (inside[p / 3] != nonVector || inside[p / 3 + dz] != nonVector) {
						intersections[p + 2] = countInter(v0, v3);
					}
				}
			}
		}
	}


	for (int i{}; i < intersections.Num() - 1; ++i) {
		if (intersections[i] != nonVector) {
			mapping[i] = vertices.Num();	
			vertices.Add(intersections[i]);
			
		}
	}
		

	// lookup
	p = 0;
	// printf("rang: %d, %d, %d\n", total_edge[0], total_edge[1], total_edge[2]);
	for (int i = 0; i < total_edge[0]; ++i) {
		for (int j = 0; j < total_edge[1]; ++j) {
			for (int k = 0; k < total_edge[2]; ++k, ++p) {
				// count status
				int status = 0;
				for (int sta_i = 0, tw = 1; sta_i < 8; ++sta_i, tw = tw << 1) {
					status |= (inside[p + m_offset[sta_i]]!= nonVector) ? tw : 0;
				}
				if (status == 0 || status == 255) continue;

				// count triangle
				for (int tri_p = 0; tri_p < 16 && Host_TRI_TABLE[status][tri_p] >= 0; tri_p+=3) {
					int edge1 = Host_TRI_TABLE[status][tri_p];
					int edge2 = Host_TRI_TABLE[status][tri_p + 1];
					int edge3 = Host_TRI_TABLE[status][tri_p + 2];

					triangleIdx.Add(mapping[offset_edge[edge1] + p * 3]);
					triangleIdx.Add(mapping[offset_edge[edge2] + p * 3]);
					triangleIdx.Add(mapping[offset_edge[edge3] + p * 3]);
				}

			}
			p += dz;
		}
		p += dy;
	}
}

FVector AParticle::checkClosetParticle(FVector v)
{
	double min_dis2 = l * l + 0.1;
	FVector result= FVector(-9999,-9999,-9999);
	for (auto  particle_1 : particles) {
		double dis2 = ( particle_1.Position - v).SizeSquared();
		
		if (dis2  < min_dis2) {
			min_dis2 = dis2;
			result =  particle_1.Position;
		}
	}
	if (min_dis2 < smoothingRadius * smoothingRadius) return result;
	return FVector(-9999,-9999,-9999);
}

FVector AParticle::countInter(const FVector& a, const FVector& b)
{
	float alpha = 0.2;
	//return FMath::Lerp(a, b, alpha);;
	return (a+b)*0.1;
}
FVector AParticle::CalculatePropertyGradient(FVector samplePoint)
{
	FVector propertyGradient = FVector::Zero();
	for (int i{}; i < particles.Num(); ++i) {
		if (!particles[i].bound) {
			float dst = (particles[i].Position - samplePoint).Length();

			FVector dir = FVector::ZeroVector;
			if (dst != 0)
				dir = (particles[i].Position - samplePoint) / dst;

			float slope = SmoothingKernelDerivative(dst, smoothingRadius);
			float density = CalculateDensity(particles[i].Position);
			propertyGradient += -particles[i].properties * dir * slope * mass / density;
		}
	}
	return propertyGradient;

}
float AParticle::ConvertDensityToPressure(float density)
{
	float densityError = density - targetDensity;
	float pressure = densityError * pressureMultiplier;
	return pressure;

}
void AParticle::UpdateSpatialLookup(const TArray< particle_1>& Points)
{
	particles = Points;

	ParallelFor(numParticles,
		[&](int32 i)
		{
			int32 CellX, CellY, CellZ;
			PositionToCellCoord(particles[i].Position, CellX, CellY, CellZ);
			uint32 CellKey = GetKeyFromHash(HashCell(CellX, CellY, CellZ));
			SpatialLookup[i] = FEntry(i, CellKey);
			StartIndices[i] = INT_MAX;
		}
	);

	SpatialLookup.Sort([](const FEntry& A, const FEntry& B)
		{
			return A.CellKey < B.CellKey;
		});

	ParallelFor(numParticles,
		[&](int32 i)
		{
			uint32 Key = SpatialLookup[i].CellKey;
			uint32 KeyPrev = i == 0 ? UINT32_MAX : SpatialLookup[i - 1].CellKey;
			if (Key != KeyPrev) {
				StartIndices[Key] = i;
			}
		}
	);
}

uint32 AParticle::HashCell(int32 CellX, int32 CellY, int32 CellZ) const
{
	uint32 a = static_cast<uint32>(CellX) * 15823;
	uint32 b = static_cast<uint32>(CellY) * 9737333;
	uint32 c = static_cast<uint32>(CellZ) * 76231;

	return a + b + c;
}
uint32 AParticle::GetKeyFromHash(uint32 Hash) const
{
	return Hash% (uint32)SpatialLookup.Num();
}
void AParticle::ForeachPointWithinRadius(const FVector& SamplePoint)
{
	int32 CellX, CellY, CellZ;
	PositionToCellCoord(SamplePoint, CellX, CellY, CellZ);
	float SqrRadius = smoothingRadius * smoothingRadius;

	for (const FIntVector& Offset : CellOffsets) {  // CellOffsets는 3차원 벡터로 변경
		uint32 Key = GetKeyFromHash(HashCell(CellX + Offset.X, CellY + Offset.Y, CellZ + Offset.Z));
		int32 CellStartIndex = StartIndices.IsValidIndex(Key) ? StartIndices[Key] : -1;

		if (CellStartIndex != -1) {
			for (int32 i = CellStartIndex; i < SpatialLookup.Num(); ++i) {
				if (SpatialLookup[i].CellKey != Key) {
					break;
				}

				int32 ParticleIndex = SpatialLookup[i].Index;
				float SqrDst = FVector::DistSquared(particles[ParticleIndex].Position, SamplePoint);
				if (SqrDst <= SqrRadius) {
					
				}
			}
		}
	}
}
void AParticle::PositionToCellCoord(const FVector& Position, int32& OutCellX, int32& OutCellY, int32& OutCellZ) const
{
	OutCellX = static_cast<int32>(Position.X / smoothingRadius);
	OutCellY = static_cast<int32>(Position.Y / smoothingRadius);
	OutCellZ = static_cast<int32>(Position.Z / smoothingRadius);
}




