// Copyright Epic Games, Inc. All Rights Reserved.

using System.IO;
using UnrealBuildTool;
using UnrealBuildTool.Rules;

public class CUDATest : ModuleRules
{
    private string poject_root_path
    {
        get { return Path.Combine(ModuleDirectory, "../.."); }
    }
    public CUDATest(ReadOnlyTargetRules Target) : base(Target)
	{
		PCHUsage = PCHUsageMode.UseExplicitOrSharedPCHs;

		PublicDependencyModuleNames.AddRange(new string[] { "Core", "CoreUObject", "Engine", "InputCore", "EnhancedInput", "ProceduralMeshComponent" });
        string custom_cuda_lib_include = "CUDALib/include";
        string custom_cuda_lib_lib = "CUDALib/lib";

        PublicIncludePaths.Add(Path.Combine(poject_root_path, custom_cuda_lib_include));
        PublicAdditionalLibraries.Add(Path.Combine(poject_root_path, custom_cuda_lib_lib, "CudaTest1.lib"));

        string cuda_path = "C:/Program Files/NVIDIA GPU Computing Toolkit/CUDA/v12.3";
        string cuda_include = "include";
        string cuda_lib = "lib/x64";

        PublicIncludePaths.Add(Path.Combine(cuda_path, cuda_include));
        //PublicAdditionalLibraries.Add(Path.Combine(cuda_path, cuda_lib, "cudart.lib"));
        PublicAdditionalLibraries.Add(Path.Combine(cuda_path, cuda_lib, "cudart_static.lib"));
    }
}
