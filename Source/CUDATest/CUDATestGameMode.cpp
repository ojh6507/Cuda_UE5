// Copyright Epic Games, Inc. All Rights Reserved.

#include "CUDATestGameMode.h"
#include "CUDATestCharacter.h"
#include "UObject/ConstructorHelpers.h"

ACUDATestGameMode::ACUDATestGameMode()
{
	// set default pawn class to our Blueprinted character
	static ConstructorHelpers::FClassFinder<APawn> PlayerPawnBPClass(TEXT("/Game/ThirdPerson/Blueprints/BP_ThirdPersonCharacter"));
	if (PlayerPawnBPClass.Class != NULL)
	{
		DefaultPawnClass = PlayerPawnBPClass.Class;
	}
}
