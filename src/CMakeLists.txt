if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
	cmake_minimum_required(VERSION 3.20.0)
	project(CADTest LANGUAGES CXX)
endif()

set(CADTEST_FILES
		
	Geometry.cpp Geometry.h
	GeometryUtil.cpp GeometryUtil.h
	Transform.cpp Transform.h
	Face.cpp Face.h
	Solid.cpp Solid.h
	
	Mesh.cpp Mesh.h
	MeshFactory.cpp MeshFactory.h
	MeshKernelLinkedTriangles.cpp MeshKernelLinkedTriangles.h
	MeshMeasurements.cpp MeshMeasurements.h
	MeshUtil.cpp MeshUtil.h
	
	BezierUtil.cpp BezierUtil.h
	BezierSurface.cpp BezierSurface.h
	
	MeshTessellate.cpp MeshTessellate.h
	MeshBoolean.cpp MeshBoolean.h
	NurbsBoolean.cpp NurbsBoolean.h
	NurbsIntersection.cpp NurbsIntersection.h
	NurbsTrimmedSurface.cpp NurbsTrimmedSurface.h
	STLFile.cpp STLFile.h
	OBJFile.cpp OBJFile.h
	
	NurbsBasis.cpp NurbsBasis.h
	NurbsKnots.cpp NurbsKnots.h
	NurbsCurve.cpp NurbsCurve.h
	NurbsSurface.cpp NurbsSurface.h
	NurbsExtrude.cpp NurbsExtrude.h
	NurbsRuled.cpp NurbsRuled.h
	NurbsFillet.cpp NurbsFillet.h
	NurbsRevolve.cpp NurbsRevolve.h
	NurbsSweep.cpp NurbsSweep.h
	NurbsFitting.cpp NurbsFitting.h
	NurbsCurveJoin.cpp NurbsCurveJoin.h
	NurbsUtil.cpp NurbsUtil.h
	NurbsOffset.cpp NurbsOffset.h
	NurbsFactory.cpp NurbsFactory.h
	NurbsSolid.cpp NurbsSolid.h
	StepFile.cpp StepFile.h
	
	Renderer.cpp Renderer.h
	RendererCamera.cpp RendererCamera.h
	RendererLight.cpp RendererLight.h
	Image.h
	ImageGenerator.cpp ImageGenerator.h
	ImageIoBmp.cpp ImageIoBmp.h
	ImageUtil.cpp ImageUtil.h
	Color.h	
	BoundingBox.cpp BoundingBox.h
	LinearAlgebra.cpp LinearAlgebra.h
)

add_library(libCADTest STATIC ${CADTEST_FILES})
target_include_directories(libCADTest PUBLIC .)
target_compile_features(libCADTest PUBLIC cxx_std_11)